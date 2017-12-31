/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2012, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


#include <dai/dai_config.h>
#ifdef DAI_WITH_GLC


#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <dai/glc.h>
#include <dai/util.h>
#include <dai/alldai.h>


namespace dai {


using namespace std;


void GLC::setProperties( const PropertySet &opts ) {
    DAI_ASSERT( opts.hasKey("tol") );
    DAI_ASSERT( opts.hasKey("maxiter") );
    DAI_ASSERT( opts.hasKey("cavity") );
    DAI_ASSERT( opts.hasKey("updates") );
    DAI_ASSERT( opts.hasKey("rgntype") );
    DAI_ASSERT( opts.hasKey("cavainame") );
    DAI_ASSERT( opts.hasKey("cavaiopts") );
    DAI_ASSERT( opts.hasKey("inainame") );
    DAI_ASSERT( opts.hasKey("inaiopts") );

    if( opts.hasKey("maxtime") )
        props.maxtime = opts.getStringAs<Real>("maxtime");
    else
        props.maxtime = INFINITY;

    props.tol = opts.getStringAs<Real>("tol");
    props.maxiter = opts.getStringAs<size_t>("maxiter");
    props.cavity = opts.getStringAs<Properties::CavityType>("cavity");
    props.updates = opts.getStringAs<Properties::UpdateType>("updates");
    props.rgntype = opts.getStringAs<Properties::RegionType>("rgntype");
    if(opts.hasKey("neighbors"))
        props.neighbors = opts.getStringAs<CobwebGraph::NeighborType>("neighbors");
    else
        props.neighbors = CobwebGraph::NeighborType::CLOSEST;
    if( opts.hasKey("verbose") )
        props.verbose = opts.getStringAs<size_t>("verbose");
    else
        props.verbose = 0;
    if( opts.hasKey("loopdepth") )
        props.loopdepth = opts.getStringAs<size_t>("loopdepth");
    else
        props.loopdepth = 4;
    if( opts.hasKey("cavainame") )
        props.cavainame = opts.getStringAs<string>("cavainame");
    if( opts.hasKey("cavaiopts") )
        props.cavaiopts = opts.getStringAs<PropertySet>("cavaiopts");
    if( opts.hasKey("inainame") )
        props.inainame = opts.getStringAs<string>("inainame");
    if( opts.hasKey("inaiopts") )
        props.inaiopts = opts.getStringAs<PropertySet>("inaiopts");
    if( opts.hasKey("reinit") )
        props.reinit = opts.getStringAs<bool>("reinit");
}


PropertySet GLC::getProperties() const {
    PropertySet opts;
    opts.set( "tol", props.tol );
    opts.set( "maxiter", props.maxiter );
    opts.set( "verbose", props.verbose );
    opts.set( "loopdepth", props.loopdepth );
    opts.set( "cavity", props.cavity );
    opts.set( "neighbors", props.neighbors );
    opts.set( "updates", props.updates );
    opts.set( "rgntype", props.rgntype );
    opts.set( "cavainame", props.cavainame );
    opts.set( "cavaiopts", props.cavaiopts );
    opts.set( "inainame", props.inainame );
    opts.set( "inaiopts", props.inaiopts );
    opts.set( "reinit", props.reinit );
    opts.set( "maxtime", props.maxtime );

    return opts;
}


string GLC::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "tol=" << props.tol << ",";
    s << "maxiter=" << props.maxiter << ",";
    s << "verbose=" << props.verbose << ",";
    s << "loopdepth=" << props.loopdepth << ",";
    s << "cavity=" << props.cavity << ",";
    s << "updates=" << props.updates << ",";
    s << "neighbors=" << props.neighbors << ",";
    s << "rgntype=" << props.rgntype << ",";
    s << "cavainame=" << props.cavainame << ",";
    s << "cavaiopts=" << props.cavaiopts << ",";
    s << "inainame=" << props.inainame << ",";
    s << "inaiopts=" << props.inaiopts << ",";
    s << "reinit=" << props.reinit << ",";
    s << "maxtime=" << props.maxtime << "]";
    return s.str();
}


GLC::GLC( const FactorGraph& fg, const PropertySet &opts ) : DAIAlgCG(fg), _CWs(),_outOffset(), _cavitydists(), _beliefs(), _factorBeliefs(), _maxdiff(0.0), _iters(0), props() {
    setProperties( opts );
    setRgn( calcRegions(), props.neighbors, false );
    if( props.verbose >= 3 ){
        cerr << "#regions: " << nrCWs() << endl;
        for( size_t i = 0; i < nrCWs(); i++ )
            cerr << INRs(i).elements() << endl;
    }
    if( props.verbose >= 3 )
        cerr << "initializing cavity dists" << endl;
    _cavitydists.resize( nrCWs() );
    _maxdiff = InitCavityDists( props.cavainame, props.cavaiopts );
    // build the vector of Cobweb Regions
    setCWs( props.inainame, props.inaiopts );
    if( props.verbose >= 3 )
        cerr << "Regions are built" << endl;
    // create beliefs
    _beliefs.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); i++ )
        _beliefs.push_back( Factor(var(i)) );
}


vector<SmallSet<size_t> > GLC::calcRegions() {
    // contains the variable indices for each region
    vector<SmallSet<size_t> > regions;
    // single variables
    if( props.rgntype == Properties::RegionType::SINGLE ) {
        regions.resize( nrVars() );
        for( size_t i = 0; i < nrVars(); i++ )
            regions[i] = SmallSet<size_t>(i);
    // partitioning factors
    } else if( props.rgntype == Properties::RegionType::FACTOR ) {
        SmallSet<size_t> remainingVars;
        for( size_t i = 0; i < nrVars(); i++ )
            remainingVars.insert( i );
        vector<SmallSet<size_t> > facVars( nrFactors() );
        for( size_t I = 0; I < nrFactors(); I++ )
            facVars[I] = findVars( factor(I).vars() );
        bool canadd = true;
        while( canadd ) {
            canadd = false;
            for( size_t I = 0; I < nrFactors(); I++ ) {
                if( facVars[I] << remainingVars ) {
                    regions.push_back( facVars[I] );
                    remainingVars /= facVars[I];
                    canadd = true;
                    break;
                }
            }
        }
        // add the variables that are not covered as single variable regions
        bforeach(size_t i, remainingVars)
            regions.push_back( SmallSet<size_t>(i) );
    // overlapping factors (non-partitioning)
    } else if( props.rgntype == Properties::RegionType::OVFACTOR ) {
        SmallSet<size_t> remainingVars;
        for( size_t I = 0; I < nrFactors(); I++ )
            regions.push_back( findVars( factor(I).vars() ) );
    // partitioning: adds loops over variables that are not covered recursively
    } else if( props.rgntype == Properties::RegionType::LOOP ) {
        set<SmallSet<size_t> > scl; // to contain clusters
        SmallSet<size_t> remaining; // variables not in a cluster yet
        // build the set of all var inds
        for( size_t v = 0; v < nrVars(); v++ )
            remaining.insert( v );
        // while there for a remaining var we can find a loop
        bool changing = true;
        while( changing ) {
            changing = false;
            for( size_t i0 = 0; i0 < nrVars(); i0++ ) {
                if( !remaining.contains(i0) )
                    continue;
                size_t currentsize = scl.size();
                if( props.loopdepth > 1 )
                    findLoopClusters( remaining, scl, SmallSet<size_t>(i0),i0, props.loopdepth - 1, deltai(i0) & remaining );
                if(scl.size() > currentsize)
                    changing = true;
            }
        }
        copy( scl.begin(), scl.end(), back_inserter(regions) );
        // add the vars that are not covered as single variable regions
        bforeach( size_t i, remaining )
            regions.push_back( SmallSet<size_t>(i) );
    // adds overlapping loops of maximum size loopdepth recursively
    } else if( props.rgntype == Properties::RegionType::OVLOOP ) {
        for( size_t I = 0; I < nrFactors(); I++ )
            regions.push_back( findVars( factor(I).vars() ) );

        SmallSet<size_t> remaining; // variables not in a cluster yet
        // build the set of all var inds
        for( size_t v = 0; v < nrVars(); v++ )
            remaining.insert( v );

        set<SmallSet<size_t> > scl; // to contain clusters
        // while there for a remaining var we can find a loop

        for( size_t i0 = 0; i0 < nrVars(); i0++ )
            findOVLoopClusters( remaining, scl, SmallSet<size_t>(i0),i0, props.loopdepth - 1, deltai(i0) );
        copy( scl.begin(), scl.end(), back_inserter(regions) );
    // non-partitioning: overlapping (variable+markov blanket)s
    } else if( props.rgntype == Properties::RegionType::OVDELTA ) {
        for( size_t I = 0; I < nrFactors(); I++ )
            regions.push_back( findVars(factor(I).vars()) );
        SmallSet<size_t> remaining; // variables not in a cluster yet
        // build the set of all var inds
        for( size_t v = 0; v < nrVars(); v++ )
            remaining.insert( v );
        set<SmallSet<size_t> > scl; // to contain clusters
        // while there for a remaining var we can find a loop
        for( size_t i0 = 0; i0 < nrVars(); i0++ )
            findOVLoopClusters( remaining, scl, SmallSet<size_t>(i0),i0, props.loopdepth - 1, deltai(i0) );
        copy( scl.begin(), scl.end(), back_inserter(regions) );
    } else if( props.rgntype == Properties::RegionType::DELTA ) {
        SmallSet<size_t> remaining; // variables not in a cluster yet
        // build the set of all var inds
        for( size_t v = 0; v < nrVars(); v++ )
            remaining.insert( v );
        while( remaining.size() > 0 ) {
            size_t cv = remaining.elements()[ rand() % remaining.size() ];
            regions.push_back( Deltai(cv) & remaining );
            remaining /= Deltai(cv);
        }
    }
    return regions;
}


void GLC::findOVLoopClusters( SmallSet<size_t>& remaining, set<SmallSet<size_t> > &allcl, SmallSet<size_t> newcl, const size_t& root, size_t length, SmallSet<size_t> vars ) {
    for( SmallSet<size_t>::const_iterator in = vars.begin(); in != vars.end(); in++ ) {
        SmallSet<size_t> ind = deltai( *in );
        if( (newcl.size()) >= 2 && ind.contains( root ) ) {
            allcl.insert( newcl | *in );
            remaining /= (newcl | *in);
        } else if( length > 1 )
            findOVLoopClusters( remaining, allcl, newcl | *in, root, length - 1, ind / newcl );
    }
}


void GLC::findLoopClusters( SmallSet<size_t>& remaining, set<SmallSet<size_t> > &allcl, SmallSet<size_t> newcl, const size_t& root, size_t length, SmallSet<size_t> vars ) {
    for( SmallSet<size_t>::const_iterator in = vars.begin(); in != vars.end(); in++ ) {
        SmallSet<size_t> ind = deltai( *in ) & remaining;
        if( (newcl.size()) >= 2 && ind.contains( root ) ) {
            allcl.insert( newcl | *in );
            remaining /= (newcl | *in);
            break;
        } else if( length > 1 )
            findLoopClusters( remaining, allcl, newcl | *in, root, length - 1, ind / newcl );
    }
}


void GLC::CalcBelief( size_t i, bool isFinal ) {
    if( !isFinal )
        _beliefs[i] = CW( var2CW(i)[0] ).belief( var(i) );
    else{ // take geometric average
        _beliefs[i] = Factor(var(i), (Real)1);
        for( size_t rind = 0; rind < var2CW(i).size(); rind++ )
            _beliefs[i] *= CW( var2CW(i)[rind] ).belief( var(i) );
        _beliefs[i] ^= (Real)1 / var2CW(i).size();
    }
}


void GLC::CalcFactorBelief( size_t I ) {
    VarSet ns = factor(I).vars();
    Factor tmp( ns, 1.0 );
    size_t counter = 0;
    for( size_t R = 0; R < nrCWs(); R++ )
        if( ns << inds2vars( INRs(R).elements() ) ) {
            tmp *= CW(R).belief(ns);
            counter++;
        }
    if( counter > 0 )
        _factorBeliefs[ns] = tmp ^ (1.0 / counter);
}


Factor GLC::belief( const VarSet &ns ) const {
    if( ns.size() == 0 )
        return Factor();
    else if( ns.size() == 1 )
        return beliefV( findVar( *(ns.begin()) ) );
    else {
        map<VarSet, Factor>::const_iterator it = _factorBeliefs.find(ns);
        if( it != _factorBeliefs.end() )
            return it->second;
        for( map<VarSet, Factor>::const_iterator it = _factorBeliefs.begin(); it != _factorBeliefs.end(); it++ )
            if( ns << it->second.vars() )
                return it->second.marginal( ns );
    }
    DAI_THROW(BELIEF_NOT_AVAILABLE);
    return Factor();
}


vector<Factor> GLC::beliefs() const {
    vector<Factor> result = _beliefs;
    for( map<VarSet, Factor>::const_iterator it = _factorBeliefs.begin(); it != _factorBeliefs.end(); it++ )
        result.push_back( it->second );
    return result;
}


Real GLC::CalcCavityDist( size_t R, const std::string &name, const PropertySet& opts ) {
    vector<Factor> pairbeliefs;
    Real maxdiff = 0;

    if( props.cavity != Properties::CavityType::UNIFORM ) {
        InfAlg *cav = newInfAlg( name, *this, opts );
        cav->makeRegionCavity( Rfs(R).elements() );

        if( props.cavity == Properties::CavityType::FULL )
            pairbeliefs.push_back( calcMarginal( *cav, inds2vars(EXRs(R).elements()), props.reinit ) );
        if( props.cavity == Properties::CavityType::PAIR )
            pairbeliefs = calcPairBeliefs( *cav, inds2vars(EXRs(R).elements()), props.reinit, false );
        else if( props.cavity == Properties::CavityType::PAIR2 )
            pairbeliefs = calcPairBeliefs( *cav, inds2vars(EXRs(R).elements()), props.reinit, true );
        maxdiff = cav->maxDiff();
        delete cav;
    }
    if( props.verbose >= 3 )
        cerr << "R:" << R << "cavity of size " << pairbeliefs.size() << endl;
    _cavitydists[R] = pairbeliefs;
    return maxdiff;
}


Real GLC::InitCavityDists( const std::string &name, const PropertySet &opts ) {
    double tic = toc();

    if( props.verbose >= 2 ) {
        cerr << this->name() << "::InitCavityDists:  ";
        if( props.cavity == Properties::CavityType::UNIFORM )
            cerr << "Using uniform initial cavity distributions" << endl;
        else if( props.cavity == Properties::CavityType::FULL )
            cerr << "Using full " << name << opts << "...";
        else if( props.cavity == Properties::CavityType::PAIR )
            cerr << "Using pairwise " << name << opts << "...";
        else if( props.cavity == Properties::CavityType::PAIR2 )
            cerr << "Using pairwise(new) " << name << opts << "...";
    }

    Real maxdiff = 0.0;
    for( size_t R = 0; R < nrCWs(); R++ ) {
        Real md = CalcCavityDist(R, name, opts);
        if( md > maxdiff )
            maxdiff = md;
    }
    if( props.verbose >= 2 )
        cerr << this->name() << "::InitCavityDists used " << toc() - tic << " seconds." << endl;
    return maxdiff;
}


void GLC::NewPancake( size_t R, size_t _R2 ) {
    Factor newmsg;
    // calculate the new msg
    newmsg = CW(M(R, _R2).his).marginal( M(R, _R2).msg.vars(), M(R, _R2).fc) / ((CW(R).marginal( M(R, _R2).msg.vars(), M(R, _R2).fc))*M(R, _R2).msg.inverse());
    newmsg.normalize();
    // update the cobweb-graph with this new msg
    M(R, _R2).msg = newmsg;
    // update the corresponding factor in the cobweb region so to be used in future marginalizations (-_R2 - 1) is the index (see cobweb class documents)
    CW(R).updateFactor( -_R2 - 1 ,newmsg);
}


void GLC::OVNewPancake( size_t R ) {
    Factor newmsg, allmsgs, marg;
    for( size_t _R2 = 0; _R2 < M(R).size(); _R2++ ) { // for all _R2 that send messages to R
        // calculate the message
        newmsg = (CW(M(R, _R2).his).marginal( M(R, _R2).msg.vars(), M(R, _R2).fc) / ((CW(R).marginal( M(R, _R2).msg.vars(), M(R, _R2).fc))))*M(R, _R2).msg;
        newmsg.normalize();
        // allmsgs takes care of the double-counts the way using the region-graph does in the paper. It is more compact but not as efficient
        allmsgs = newmsg;
        for( size_t sub = 0; sub < M(R,_R2).subregions.size(); sub++ ) { //for all descendants of rho of (\ominus r_{R, R2})
            Real c = (Real)cn(R)[M(R, _R2).subregions[sub]].first; // counting number for rho
            size_t nrParents = cn(R)[M(R, _R2).subregions[sub]].second.size(); // number of rho's parents
            marg = newmsg.marginal(M(R, _R2).subregions[sub], false); // marginal of the message from R2 to R over rho
            allmsgs *= (marg^(c / nrParents)); // discount the double counting of rho from the message
            for( size_t pind = 0; pind < nrParents; pind++ ) // this operation corresponds to the upward pass
                if( cn(R)[M(R, _R2).subregions[sub]].second[pind] != _R2 )
                    M(R, cn(R)[M(R, _R2).subregions[sub]].second[pind]).newmsg *= marg^(-c / (nrParents * (nrParents - (Real)1)));
        }
        // set the message after taking double-countings into account
        M(R, _R2).msg = allmsgs.normalized();
    }

    for( size_t _R2 = 0; _R2 < M(R).size(); _R2++ ) {
        M(R, _R2).newmsg *= M(R, _R2).msg;
        M(R, _R2).newmsg.normalize();
        // update the corresponding factor in the cobweb region: to be used in future marginalizations
        CW(R).updateFactor( -_R2 - 1 , M(R, _R2).msg, false);
        M(R, _R2).msg = M(R, _R2).newmsg;
        M(R, _R2).newmsg.fill( (Real)1 );
    }
}


Real GLC::run() {
    if( props.verbose >= 2 )
        cerr << "Starting " << identify() << "...";
    if( props.verbose >= 2 )
        cerr << endl;

    double tic = toc();

    if( props.verbose >= 3 )
        cerr << "Initializing older beliefs" << endl;
    vector<vector<Factor> > oldBeliefsM;
    vector<Factor> oldBeliefsP;
    if( isPartition ) {
        oldBeliefsM.reserve( nrCWs() );
        for( size_t R = 0; R < nrCWs(); R++ ) {
            oldBeliefsM.push_back( vector<Factor>() );
            oldBeliefsM[R].reserve( M(R).size() );
            for( size_t m = 0; m < M(R).size(); m++ )
                oldBeliefsM[R].push_back( M(R,m).msg );
        }
    } else {
        oldBeliefsP.reserve( nrVars() );
        for( size_t i = 0; i < nrVars(); i++ ) {
            CalcBelief( i, true );
            oldBeliefsP.push_back( _beliefs[i] );
        }
    }
    if( props.verbose >= 3 )
        cerr << "Setting the update sequence" <<endl;

    size_t nredges = 0;
    vector<Edge> update_seq;
    update_seq.reserve( nredges );
    for( size_t R = 0; R < nrCWs(); ++R )
        for( size_t R2 = 0; R2 < M(R).size(); R2++ ) {
            update_seq.push_back( Edge( R, R2 ) );
            nredges++;
        }

    vector<size_t> update_seq_node;
    update_seq_node.reserve( nrCWs() );
    for( size_t R = 0; R < nrCWs(); ++R )
        update_seq_node.push_back( R );
    if( props.verbose >= 3 )
        cerr << "Going into the main loop" << endl;

    Real maxDiff = INFINITY;
    for( _iters = 0; _iters < props.maxiter && maxDiff > props.tol && (toc() - tic) < props.maxtime; _iters++ ) {
        if( props.updates == Properties::UpdateType::SEQRND ) {
            random_shuffle( update_seq.begin(), update_seq.end() );
            random_shuffle( update_seq_node.begin(), update_seq_node.end() );
        }
        if( isPartition ) {
            for( size_t t = 0; t < nredges; t++ ) {
                size_t R = update_seq[t].first;
                size_t _R2 = update_seq[t].second;
                NewPancake( R, _R2);
            }
        } else {
            for( size_t t = 0; t < nrCWs(); t++ )
                OVNewPancake( update_seq_node[t] );
        }

        if( !isPartition )
            for( size_t i = 0; i < nrVars(); i++ )
                CalcBelief( i , true );

        maxDiff = -INFINITY;
        if( isPartition ) {
            for( size_t R = 0; R < nrCWs(); R++ ) {
                for( size_t m = 0; m < M(R).size(); m++ ) {
                    maxDiff = std::max( maxDiff, dist( M(R, m).msg, oldBeliefsM[R][m], DISTLINF ) );
                    oldBeliefsM[R][m] = M(R,m).msg;
                }
            }
        } else {
            for( size_t i = 0; i < nrVars(); i++ ) {
                maxDiff = std::max( maxDiff, dist( _beliefs[i], oldBeliefsP[i], DISTLINF ) );
                oldBeliefsP[i] = _beliefs[i];
            }
        }

        if( props.verbose >= 2 )
            cerr << this->name() << "::run:  maxdiff " << maxDiff << " after " << _iters+1 << " passes" << endl;
    }

    if( isPartition )
        for( size_t i = 0; i < nrVars(); i++ )
            CalcBelief( i, true );

    if( maxDiff > _maxdiff )
        _maxdiff = maxDiff;

    if( props.verbose >= 2 ) {
        if( maxDiff > props.tol ) {
            if( props.verbose == 2 )
                cerr << endl;
                cerr << this->name() << "::run:  WARNING: not converged within " << props.maxiter << " passes (" << toc() - tic << " seconds)...final maxdiff:" << maxDiff << endl;
        } else {
            if( props.verbose >= 2 )
                cerr << this->name() << "::run:  ";
                cerr << "converged in " << _iters << " passes (" << toc() - tic << " seconds)." << endl;
        }
    }

    for( size_t I = 0; I < nrFactors(); I++ )
        CalcFactorBelief( I );

    return maxDiff;
}


void GLC::setCWs( const std::string &name, const PropertySet &opts ) {
    // needs to generate a vector of relevant factors, and a dictionary for mapping between local and global factors
    _CWs.clear();
    _CWs.reserve( _INRs.size() );
    _outOffset.clear();
    _outOffset.reserve( _INRs.size() );
    // main factors
    if( props.verbose >= 3 )
        cerr << "Setting CWs..." << endl;

    for( size_t R = 0; R < nrCWs(); R++ ) { // for each region
        map<int, size_t> g2l;
        vector<Factor> Rfactors; // contains all factors
        size_t potsize = Rfs(R).size() + M(R).size() + _cavitydists[R].size();
        Rfactors.reserve( potsize );
        size_t _I = 0;
        // adding main factors
        for( SmallSet<size_t>::const_iterator it = Rfs(R).begin(); it != Rfs(R).end(); it++ ) {
            Rfactors.push_back( factor(*it) );
            g2l[(int)(*it)] = _I;
            _I++;
        }

        // adding factors for incoming messages
        size_t m = 0;
        size_t offset = Rfs(R).size();
        for( vector<Connection>::const_iterator it = M(R).begin(); it != M(R).end(); it++ ) {
            Rfactors.push_back( (*it).msg );
            g2l[(int)-m-1] = m + offset;
            m++;
        }
        // cavities
        offset = Rfs(R).size() + M(R).size();
        for( size_t c = 0; c < _cavitydists[R].size(); c++ ) {
            Rfactors.push_back( _cavitydists[R][c] );
            g2l[(int)-M(R).size() - c - 1] = offset + c;
        }
        // outgoing msgs
        offset = Rfs(R).size() + M(R).size() + _cavitydists[R].size();
//      _outOffset.push_back(-M(R).size() - _cavitydists[R].size() - 1);
        size_t counter = 0;
        if( props.verbose >= 3 )
            cerr << "going to the loop!" << endl;
        for( size_t c = 0; c < outM(R).size(); c++ ) {
            bool isAvailable = false;
            for( size_t I = 0; I < Rfactors.size(); I++ )
                if( outM(R,c) << Rfactors[I].vars() ) {
                    isAvailable = true;
                    break;
                }
            if( !isAvailable ) {
                Rfactors.push_back( Factor( outM(R)[c], (Real)1 ) );
                g2l[(int)-M(R).size() - _cavitydists[R].size() - 1 - counter] = offset + counter;
                counter++;
            }
        }

        _CWs.push_back( Cobweb(g2l) );
        InfAlg *alg = newInfAlg( name, FactorGraph(Rfactors), opts );
        DAIAlgFG *dalg = static_cast<DAIAlgFG*> (alg);
        _CWs[R].setInfAlg( dalg );
    }
}


void GLC::initCWs(){
    for( size_t R = 0; R < nrCWs(); R++ ) {
        _CWs[R].initialize();
        for( size_t m = 0; m < M(R).size(); m++ ) {
            M(R)[m].msg.fill( (Real)1 );
            M(R)[m].newmsg.fill( (Real)1 );
        }
    }
}


} // end of namespace dai


#endif
