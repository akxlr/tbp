/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2012, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


/// \file
/// \brief Defines classes GLC and Cobweb, which implement the "Generalized Loop Correction method"
/// \todo Fix the init of GLC


#ifndef __defined_libdai_glc_h
#define __defined_libdai_glc_h


#include <dai/dai_config.h>
#ifdef DAI_WITH_GLC


#include <algorithm>
#include <set>
#include <string>
#include <dai/util.h>
#include <dai/daialg.h>
#include <dai/enum.h>
#include <dai/cobwebgraph.h>
#include <dai/properties.h>
#include <dai/exceptions.h>


namespace dai {


/// Represents a region for the GLC algorithm.
/** This region contains all the proper factors from the factor graph. It also contains a factor for incoming messages and the
 *  cavity distribution over this region. The outgoing messages also receive a factor mainly so that we can easily calculate the marginal
 *  of the region over the outgoing message variables.
 *
 *  \author Siamak Ravanbakhsh
 */
class Cobweb {
    public:
        /// Pointer to the approximate inference algorithm and factor graph over region R (provided in parameters of GLC)
        DAIAlgFG *cav;

        /// Global to local factor index conversion.
        /** Factor(I) in original graph corresponds to cav->factor(_g2l[I]) in this region
         *  the j'th incoming message to this region M(R,j) corresponds to cav->factor(_g2l[-j-1])
         *  the j'th cavity factor (for pairwise case),_cavitydists[R][j] corresponds to cav->factor(_g2l[-m-j-1])
         *      when m = M(R).size()
         *      note that for full cavity the _cavitydists[R] has a single factor and for uniform cavity it is empty
         *  the j'th outgoing message from this region outM(R,j) corresponds to cav->factor(_g2l[-c-m-j-1])
         *      when m  = M(R).size() and c = _cavitydists[R].size()
         */
        std::map<int, size_t> _g2l;

        /// Default constructor
        Cobweb(){}

        /// Construct from global-to-local factor index conversion structure
        Cobweb( const std::map<int, size_t>& g2l ): cav(0), _g2l(g2l) {}

        /// Default destructor
        virtual ~Cobweb(){ delete cav; }

        /// Sets the factor graph and inference algorithm for this region
        void setInfAlg( DAIAlgFG* alg ) {
            cav = alg;
        }

        /// Returns the factor of this region globaly indexed by \a I (mainly used to access messages to and from this region)
        Factor factor( size_t I ) {
            return cav->factor( _g2l[I] );
        }

        /// Initializes the inference algorithm over this region
        void initialize() {
            cav->init();
        }

        /// Calculates the marginal over \a ns after temporarily removing the factors indexed by \a rmgind
        Factor marginal( const VarSet& ns, const std::vector<size_t>& rmgind ) {
            // Set the factors in rmgind to one
            VarSet vs;
            std::vector<size_t> rmlinds; // local index
            rmlinds.reserve( rmgind.size() );
            for( std::vector<size_t>::const_iterator it = rmgind.begin(); it != rmgind.end(); it++ ) {
                vs |= cav->factor( _g2l[*it] ).vars();
                rmlinds.push_back( _g2l[*it] );
            }
            cav->makeRegionCavity( rmlinds, true );
            // Initialize if necessary
            cav->init( vs );
            cav->run();
            Factor result;
            result = cav->belief( ns );
            cav->restoreFactors();
            return result;
        }

        /// Calculates the marginal over the variables in outgoing message indexed by \a outmsgind after temporarily removing the factors indexed by \a rmgind
        /** This function is used to calculate outgoing messages from this region
         */
        Factor marginal( const int& outmsgind, const std::vector<size_t>& rmgind ) {
            // set the factors in rmgind to one
            VarSet vs;
            std::vector<size_t> rmlinds; // local index
            rmlinds.reserve( rmgind.size() );
            for( std::vector<size_t>::const_iterator it = rmgind.begin(); it != rmgind.end(); it++ ) {
                vs |= cav->factor( _g2l[*it] ).vars();
                rmlinds.push_back( _g2l[*it] );
            }
            cav->makeRegionCavity( rmlinds, true );
            // Initialize if necessary
            cav->init( vs );
            cav->run();
            Factor result;
            result = cav->beliefF( _g2l[outmsgind] );
            cav->restoreFactors();
            return result;
        }

        /// Sets (or multiplies, if \a multiply == \c true) the factor indexed by \a gind by factor \a msg
        void updateFactor( int gind, const Factor& msg, bool multiply=false ) {
            if( !multiply )
                cav->setFactor( _g2l[gind], msg, false );
            else
                cav->setFactor( _g2l[gind], (msg*(cav->factor( _g2l[gind] ))).normalized(), false );
        }

        /// Runs inference and returns belief of variable \a v
        Factor belief( const Var v ) {
            cav->run();
            return cav->belief( v );
        }

        /// Runs inference and returns belief of variables \a vs
        Factor belief( const VarSet& vs ) {
            cav->run();
            return cav->belief( vs );
        }
};


/// Approximate inference algorithm "Generalized Loop Correction" [\ref RYG12]
/** Inherits from a CobwebGraph which in turn inherits from a FactorGraph.
 *  Although a CobwebGraph represents the graph with cobweb clusters (for design reasons?), 
 *  the actual vector of cobwebs is contained in GLC class. All the other functionalities 
 *  needed to represent the cobwebgraph is in its own class. The cobweb provides the function 
 *  of marginalization using the provided inference method. This can be beneficial when using 
 *  uniform cavities (or a subset of pairwise cavities, but this is not implemented yet). 
 *  However, in general the performance does not improve using Junction Tree instead of exact marginalization...
 *
 *  \author Siamak Ravanbakhsh
 */
class GLC : public DAIAlgCG {
    private:
        /// Stores for each variable the approximate cavity distribution
        std::vector<Cobweb> _CWs;

        /// Gives the starting global index for the factors of outgoing messages of each region
        std::vector<int> _outOffset; 

        /// Cavity distributions
        /** For each region it contains:
         *  - Nothing if using UniCAV
         *  - One Factor if using FullCAV
         *  - All pairwise caivity factors if using PairCAV or Pair2CAV
         */
        std::vector<std::vector<Factor> > _cavitydists;

        /// Single variable beliefs
        std::vector<Factor> _beliefs;

        /// Beliefs over factors
        std::map<VarSet, Factor> _factorBeliefs;

        /// Maximum difference encountered so far
        Real _maxdiff;

        /// Number of iterations needed
        size_t _iters;

    public:
        /// Parameters for GLC
        struct Properties {
            /// Enumeration of possible ways to initialize the cavities
            /** The following initialization methods are defined:
             *  - FULL calculates the marginal using calcMarginal()
             *  - PAIR calculates only second order interactions using calcPairBeliefs() with \a accurate == \c false
             *  - PAIR2 calculates only second order interactions using calcPairBeliefs() with \a accurate == \c true
             *  - UNIFORM uses a uniform distribution
             */
            DAI_ENUM(CavityType,FULL,PAIR,PAIR2,UNIFORM);

            /// Enumeration of different possible types of region (\ominus r in the paper)
            /** The following initialization methods are defined:
             *  - SINGLE using single variable
             *  - FACTOR partitioning regions of that form corresponding to GLC (in the paper)
             *  - LOOP partitioning regions of that form corresponding to GLC (in the paper)
             *  - DELTA partitioning regions of that form corresponding to GLC (in the paper)
             *  - OVFACTOR, OVLOOP, OVDELTA the corresponding form of region in the overlapping form corresponding to GLC+ (in the paper)
             *  For example, using OVFACTOR each factor becomes a region (\ominus r), while using FACTOR the regions contain only
             *  non-overlapping factors and if some variables can't be covered by any remaining factors single variable regions
             *  are constructed for them.
             *  With this convention OVDELTA here corresponds to DELTA using CVM (HAK).
             */
            DAI_ENUM(RegionType,SINGLE,FACTOR,OVFACTOR,LOOP,OVLOOP,DELTA,OVDELTA)

            /// Enumeration of different update schedules
            /** The following update schedules are defined:
             *  - SEQFIX sequential fixed schedule
             *  - SEQRND sequential random schedule
             */
            DAI_ENUM(UpdateType,SEQFIX,SEQRND);

            /// Verbosity (amount of output sent to stderr)
            size_t verbose;

            // Maximum length of the loops to use in construction of loopy regions
            size_t loopdepth;

            /// Maximum number of iterations
            size_t maxiter;

            /// The algorithm will quit and report the current result at reaching this runtime
            Real maxtime;

            /// Tolerance for convergence test
            Real tol;

            /// Complete or partial reinitialization of the messages in calculating the cavity distribution?
            bool reinit;

            /// Whether or not to include auxiliary factors in a region. 
            /** The auxilary factors are the ones that only depend on variables on the boundary of a region (\oplus r - \ominus r).
             *  This sometimes produces orders of magnitude improvement (e.g., in grids), but it is not discussed in the paper.
             */
            bool auxfactors;

            /// How to initialize the cavities
            CavityType cavity;

            // Type of region
            RegionType rgntype;

            /// What update schedule to use
            UpdateType updates;

            /** Indicates what happens if a subset of variables in the boundary of a region (\ominus r_p) is shared by
             *  some neighbors such that one (\ominus r_{p,q1}) is a subset of another (\ominus r_{p,q2}).
             *  - ALL     all such neighbors are included in the the updates.
             *  - TOP     only (\ominus r_{p,q2}) is included unless (\ominus r_{p,q2} = \ominus r_{p,q1}) in which case both are included
             *  - CLOSEST (default): similar to TOP but in case of a tie the the region r_q with largest r_q \cap r_p is considered
             *  \note Not important in perfomance!
             */
            CobwebGraph::NeighborType neighbors;

            /// Name of the algorithm used to initialize the cavity distributions
            std::string cavainame;

            /// Parameters for the algorithm used to initialize the cavity distributions
            PropertySet cavaiopts;

            // Name of the algorithm used for marginalization over a region (EXACT is fine)
            std::string inainame;

            // Parameters for algorithm used for marginalization over a region
            PropertySet inaiopts;
        } props;


    public:
        /// Default constructor
        GLC() : DAIAlgCG(), _CWs(),_outOffset(), _cavitydists(), _beliefs(), _factorBeliefs(), _maxdiff(), _iters(), props() {}

        /// Construct from FactorGraph \a fg and PropertySet \a opts
        /** \param fg Factor graph.
         *  \param opts Parameters @see Properties
         */
        GLC( const FactorGraph &fg, const PropertySet &opts );

    /// \name General InfAlg interface
    //@{
        virtual GLC* clone() const { return new GLC(*this); }
        virtual GLC* construct( const FactorGraph &fg, const PropertySet &opts ) const { return new GLC( fg, opts ); }
        virtual std::string name() const { return "GLC"; }
        virtual Factor belief( const Var &v ) const { return beliefV( findVar( v ) ); }
        virtual Factor belief( const VarSet &vs ) const;
        virtual Factor beliefV( size_t i ) const { return _beliefs[i]; }
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const { DAI_THROW(NOT_IMPLEMENTED); return 0.0; }
        virtual void init(){ initCWs(); }
        virtual void init( const VarSet &/*ns*/ ) { init(); }
        virtual Real run();
        virtual Real maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
        virtual void setMaxIter( size_t maxiter ) { props.maxiter = maxiter; }
        virtual void setProperties( const PropertySet &opts );
        virtual PropertySet getProperties() const;
        virtual std::string printProperties() const;
    //@}

    /// \name Additional interface specific for GLC
    //@{
        /// Returns constant reference to outer region \a alpha
        const Cobweb& CW( size_t alpha ) const {
            DAI_DEBASSERT( alpha < nrCWs() );
            return _CWs[alpha];
        }

        /// Returns reference to outer region \a alpha
        Cobweb& CW( size_t alpha ) {
            DAI_DEBASSERT( alpha < nrCWs() );
            return _CWs[alpha];
        }

        /// Approximates the cavity distribution of variable \a i, using the inference algorithm \a name with parameters \a opts
        Real CalcCavityDist( size_t i, const std::string &name, const PropertySet &opts );

        /// Approximates all cavity distributions using inference algorithm \a name with parameters \a opts
        Real InitCavityDists( const std::string &name, const PropertySet &opts );

        /// Updates the belief of the Markov blanket of region \a R using its intersection with \a _R2'th neighbor (partitioning case).
        void NewPancake( size_t R, size_t _R2 );

        /// Updates the belief of the Markov blanket of region \a R using all its overlapping neighbors (non-partitioning case).
        /** \note The implementation uses a more compact form than the message passing over neighborhood region-graph explained in the paper.
         */
        void OVNewPancake( size_t R );

        /// Calculates the belief of variable \a i: if \a isFinal == \c true, it returns the geometric average given by all regions that include \a i
        void CalcBelief( size_t i, bool isFinal = false );

        /// Calculates the belief of factor \a I
        void CalcFactorBelief( size_t I );

        /// Designates the variables in each region of the graph based on the rgntype
        std::vector<SmallSet<size_t> > calcRegions();

        /// Initializes the cobweb regions and all the messages between regions.
        void initCWs();

        /// Initializes all the cobweb regions.
        /** Defines their factors using the cavity distribution (as factor(s))
         *  and all the anticipated messages to be exchanged with other regions (treating a message as a factor),
         *  and also it builds a mapping from cobweb factors to factors of the original factor graph and messages.
         */
        void setCWs( const std::string &name, const PropertySet &opts );

        /// Helper function for partitioning case to build the regions by finding loops of maximum given size.
        /** If a variable is contained in a loop it is not considered for inclusion in any other loop.
         */
        void findLoopClusters( SmallSet<size_t>& remaining, std::set<SmallSet<size_t> > &allcl, SmallSet<size_t> newcl, const size_t& root, size_t length, SmallSet<size_t> vars );

        /// Helper function for general case to build the regions by finding loops of maximum given size.
        /** All loops of max-size are found and variables that are not covered are contained in smaller loops or in factor-sized or single-variable regions.
         */
        void findOVLoopClusters( SmallSet<size_t>& remaining, std::set<SmallSet<size_t> > &allcl, SmallSet<size_t> newcl, const size_t& root, size_t length, SmallSet<size_t> vars );
    //@}
};


} // end of namespace dai


#endif


#endif
