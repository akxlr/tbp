/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2011, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


#include <iostream>
#include <dai/matlab/matlab.h>
#include "mex.h"
#include <dai/jtree.h>


using namespace std;
using namespace dai;


/* Convert cell vector of Matlab sets to vector<VarSet> */
vector<VarSet> mx2VarSets(const mxArray *vs, const FactorGraph &fg, long verbose, vector<Permute> &perms) {
    vector<VarSet> varsets;

    int n1 = mxGetM(vs);
    int n2 = mxGetN(vs);
    if( n2 != 1 && n1 != 1 )
        mexErrMsgTxt("varsets should be a Nx1 or 1xN cell matrix.");
    size_t nr_vs = n1;
    if( n1 == 1 )
        nr_vs = n2;

    // interpret vs, linear cell array of varsets
    varsets.reserve( nr_vs );
    perms.clear();
    perms.reserve( nr_vs );
    for( size_t cellind = 0; cellind < nr_vs; cellind++ ) {
        if( verbose >= 3 )
            cerr << "reading varset " << cellind << ": " << endl;
        mxArray *cell = mxGetCell(vs, cellind);
        if( verbose >= 3 )
            cerr << "  got cell " << endl;
        size_t nr_mem = mxGetN(cell);
        if( verbose >= 3 )
            cerr << "  number members: " << nr_mem << endl;
        double *members = mxGetPr(cell);
        if( verbose >= 3 )
            cerr << "  got them! " << endl;

        // add variables
        VarSet vsvars;
        if( verbose >= 3 )
            cerr << "  vars: ";
        vector<long> labels(nr_mem,0);
        vector<size_t> dims(nr_mem,0);
        for( size_t mi = 0; mi < nr_mem; mi++ ) {
            labels[mi] = (long)members[mi];
            dims[mi] = fg.var(labels[mi]).states();
            if( verbose >= 3 )
                cerr << labels[mi] << " ";
            vsvars.insert( fg.var(labels[mi]) );
        }
        if( verbose >= 3 )
            cerr << endl;
        DAI_ASSERT( nr_mem == vsvars.size() );
        varsets.push_back(vsvars);

        // calculate permutation matrix
        vector<size_t> perm(nr_mem,0);
        VarSet::iterator j = vsvars.begin();
        for( size_t mi = 0; mi < nr_mem; mi++,j++ ) {
            long gezocht = j->label();
            vector<long>::iterator piet = find(labels.begin(),labels.end(),gezocht);
            perm[mi] = piet - labels.begin();
        }
        if( verbose >= 3 ) {
            cerr << endl << "  perm: ";
            for( vector<size_t>::iterator r=perm.begin(); r!=perm.end(); r++ )
                cerr << *r << " ";
            cerr << endl;
        }
        // create Permute object
        vector<size_t> di(nr_mem,0);
        size_t prod = 1;
        for( size_t k = 0; k < nr_mem; k++ ) {
            di[k] = dims[k];
            prod *= dims[k];
        }
        Permute permindex( di, perm );
        perms.push_back( permindex );
    }

    if( verbose >= 3 ) {
        for(vector<VarSet>::const_iterator I=varsets.begin(); I!=varsets.end(); I++ )
            cerr << *I << endl;
    }

    return( varsets );
}


/* Input Arguments */

#define PSI_IN          prhs[0]
#define VARSETS_IN      prhs[1]
#define OPTS_IN         prhs[2]
#define NR_IN           3
#define NR_IN_OPT       0


/* Output Arguments */

#define LOGZ_OUT        plhs[0]
#define Q_OUT           plhs[1]
#define QV_OUT          plhs[2]
#define QF_OUT          plhs[3]
#define QMAP_OUT        plhs[4]
#define MARGS_OUT       plhs[5]
#define NR_OUT          3
#define NR_OUT_OPT      3


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ) {
    // Check for proper number of arguments
    if( ((nrhs < NR_IN) || (nrhs > NR_IN + NR_IN_OPT)) || ((nlhs < NR_OUT) || (nlhs > NR_OUT + NR_OUT_OPT)) ) {
        mexErrMsgTxt("Usage: [logZ,q,qv,qf,qmap,margs] = dai_jtree(psi,varsets,opts)\n\n"
        "\n"
        "INPUT:  psi        = linear cell array containing the factors\n"
        "                     (psi{i} should be a structure with a Member field\n"
        "                     and a P field).\n"
        "        varsets    = linear cell array containing varsets for which marginals\n"
        "                     are requested.\n"
        "        opts       = string of options.\n"
        "\n"
        "OUTPUT: logZ       = logarithm of the partition sum.\n"
        "        q          = linear cell array containing all calculated marginals.\n"
        "        qv         = linear cell array containing all variable marginals.\n"
        "        qf         = linear cell array containing all factor marginals.\n"
        "        qmap       = linear array containing the MAP state.\n"
        "        margs      = linear cell array containing all requested marginals.\n");
    }

    // Get psi and construct factorgraph
    vector<Factor> factors = mx2Factors(PSI_IN, 0);
    FactorGraph fg(factors);

    // Get varsets
    vector<Permute> perms;
    vector<VarSet> varsets = mx2VarSets(VARSETS_IN,fg,0,perms);

    // Get options string
    char *opts;
    size_t buflen = mxGetN( OPTS_IN ) + 1;
    opts = (char *)mxCalloc( buflen, sizeof(char) );
    mxGetString( OPTS_IN, opts, buflen );
    // Convert to options object props
    stringstream ss;
    ss << opts;
    PropertySet props;
    ss >> props;

    // Construct InfAlg object, init and run
    JTree jt = JTree( fg, props );
    jt.init();
    jt.run();

    // Save logZ
	double logZ = NAN;
    logZ = jt.logZ();

    // Hand over results to MATLAB
    LOGZ_OUT = mxCreateDoubleMatrix(1,1,mxREAL);
    *(mxGetPr(LOGZ_OUT)) = logZ;

    Q_OUT = Factors2mx(jt.beliefs());

    if( nlhs >= 3 ) {
        vector<Factor> qv;
        qv.reserve( fg.nrVars() );
        for( size_t i = 0; i < fg.nrVars(); i++ )
            qv.push_back( jt.belief( fg.var(i) ) );
        QV_OUT = Factors2mx( qv );
    }

    if( nlhs >= 4 ) {
        vector<Factor> qf;
        qf.reserve( fg.nrFactors() );
        for( size_t I = 0; I < fg.nrFactors(); I++ )
            qf.push_back( jt.belief( fg.factor(I).vars() ) );
        QF_OUT = Factors2mx( qf );
    }

    if( nlhs >= 5 ) {
        std::vector<size_t> map_state;
        bool supported = true;
        try {
            map_state = jt.findMaximum();
        } catch( Exception &e ) {
            if( e.getCode() == Exception::NOT_IMPLEMENTED )
                supported = false;
            else
                throw;
        }
        if( supported ) {
            QMAP_OUT = mxCreateNumericMatrix(map_state.size(), 1, mxUINT32_CLASS, mxREAL);
            uint32_T* qmap_p = reinterpret_cast<uint32_T *>(mxGetPr(QMAP_OUT));
            for (size_t n = 0; n < map_state.size(); ++n)
                qmap_p[n] = map_state[n];
        } else {
            mexErrMsgTxt("Calculating a MAP state is not supported by this inference algorithm.");
        }
    }

    if( nlhs >= 6 ) {
        vector<Factor> margs;
        margs.reserve( varsets.size() );
        for( size_t s = 0; s < varsets.size(); s++ ) {
            Factor marg;
            jt.init();
            jt.run();
            marg = jt.calcMarginal( varsets[s] );

            // permute entries of marg
            Factor margperm = marg;
            for( size_t li = 0; li < marg.nrStates(); li++ )
                margperm.set( li, marg[perms[s].convertLinearIndex(li)] );
            margs.push_back( margperm );
        }
        MARGS_OUT = Factors2mx( margs );
    }

    return;
}
