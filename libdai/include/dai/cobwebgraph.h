/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2012, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


/// \file
/// \brief Defines class CobwebGraph, which implements a type of region graph used by GLC.


#ifndef __defined_libdai_cobwebgraph_h
#define __defined_libdai_cobwebgraph_h


#include <iostream>
#include <dai/factorgraph.h>
#include <dai/weightedgraph.h>
#include <dai/smallset.h>
#include <dai/enum.h>
#include <algorithm>
#include <map>
#include <set>


namespace dai {


/// A CobwebGraph is a special type of region graph used by the GLC algorithm
/** \author Siamak Ravanbakhsh
 *  \todo Implement unit test for Cobwebgraph
 */
class CobwebGraph : public FactorGraph {
    public:
        /// The information in connection between two regions
        struct Connection {
            /// Index of the first region (p)
            size_t my;
            /// Index of the second region (q)
            size_t his;
            /// Index of this connection in the connections of the first region
            size_t iter;
            /// Index of the mirror of this connection in the connections of the second region. (reference of this message in _outM)
            size_t dual;
            /// The message sent from region q (his) to p (my)
            Factor msg;
            /// "Index" of variables in the message
            std::vector<size_t> varinds;
            /// Used as a temporary factor only
            /** \todo Remove CobwebGraph::Connection::newmsg
             */
            Factor newmsg;
            /// Index of factors in common
            std::vector<size_t> fc;
            /// Regions rho that are descendents of msg from q to p in p's msg-region graph (not used in partitioning case)
            std::vector<VarSet> subregions;
        };

        /** Indicates what happens if a subset of variables in the boundary of a region (\ominus r_p) is shared by
         *  some neighbors such that one (\ominus r_{p,q1}) is a subset of another (\ominus r_{p,q2}).
         *  - ALL     all such neighbors are included in the the updates.
         *  - TOP     only (\ominus r_{p,q2}) is included unless (\ominus r_{p,q2} = \ominus r_{p,q1}) in which case both are included
         *  - CLOSEST (default): similar to TOP but in case of a tie the the region r_q with largest r_q \cap r_p is considered
         *  \note Not important in perfomance!
         */
        DAI_ENUM(NeighborType,ALL,TOP,CLOSEST);

    protected:
        /// Vector of variable indices internal to each region (r)
        std::vector<SmallSet<size_t> >        _INRs;
        /// Vector of variable indices on the boundary of each region (\ominus r)
        std::vector<SmallSet<size_t> >        _EXRs;
        /// Index of factors in each region
        std::vector<SmallSet<size_t> >        _Rfs;
        /// Index of factors internal to each region, i.e., all its variables are internal to the region
        std::vector<SmallSet<size_t> >        _Rifs;
        /// Index of factors that bridge each region, i.e., not all its variables are internal to the region
        std::vector<SmallSet<size_t> >        _Rxfs;
        /// The vector of domain of messages leaving each region (\ominus r_{p,q})
        std::vector<std::vector<VarSet> >     _outM;
        /// Vector of all connections to each region
        std::vector<std::vector<Connection> > _M;

        /// For each i the index of (cobweb) regions that contain variable i
        std::vector<std::vector<size_t> > _var2CW;

        /** For each region r_p a mapping from all variables rho in its msg-region graph to a pair:
         *  first: counting number
         *  second: the index of top-regions that contain rho
         */
        std::vector<std::map<VarSet, std::pair<int,std::vector<size_t> > > > _cn;

        /// Whether a given set of region vars makes a partitioning or not
        bool isPartition;

    public:
    /// \name Constructors and destructors
    //@{
        /// Default constructor
        CobwebGraph() : FactorGraph(), _INRs(), _EXRs(),_Rfs(),_Rifs(),_Rxfs(), _M(), _var2CW(), _cn(), isPartition(true) {}

        /// Constructs a cobweb graph from a factor graph
        CobwebGraph( const FactorGraph& fg ): FactorGraph(), _INRs(), _EXRs(),_Rfs(), _M(), _var2CW(), _cn(), isPartition(true) {
            // Copy factor graph structure
            FactorGraph::operator=( fg );
        }

        /// Clone \c *this (virtual copy constructor)
        virtual CobwebGraph* clone() const { return new CobwebGraph(*this); }
    //@}

    /// \name Accessors and mutators
    //@{
        /// Returns the number of regions
        size_t nrCWs() const { return _INRs.size(); }

        /// Returns constant reference to the _cn for region \a R
        const std::map<VarSet, std::pair<int,std::vector<size_t> > >& cn( size_t R ) const {
            DAI_DEBASSERT( R < nrCWs() );
            return _cn[R];
        }

        /// Returns a reference to the _cn for region \a R
        std::map<VarSet, std::pair<int,std::vector<size_t> > >& cn( size_t R ) {
            DAI_DEBASSERT( R < nrCWs() );
            return _cn[R];
        }

        /// Returns constant reference the vector of domain of all outgoing messages from region \a R
        const std::vector<VarSet>& outM( size_t R ) const {
            DAI_DEBASSERT( R < _outM.size() );
            return _outM[R];
        }

        /// Returns reference the vector of domain of all outgoing messages from region \a R
        std::vector<VarSet>& outM( size_t R ) {
            DAI_DEBASSERT( R < _outM.size() );
            return _outM[R];
        }

        /// Returns constant reference to the variables in the outgoing message from region \a R to its \a j'th neighbor 
        /** \a j corresponds to dual in the connection construct
         */
        const VarSet& outM( size_t R, size_t j ) const {
            DAI_DEBASSERT( R < _outM.size() );
            DAI_DEBASSERT( j < _outM[R].size() );
            return _outM[R][j];
        }

        /// Returns a reference to the variables in the outgoing message from region \a R to its \a j'th neighbor 
        /** \a j corresponds to dual in the connection construct
         */
        VarSet& outM( size_t R, size_t j ) {
            DAI_DEBASSERT( R < _outM.size() );
            DAI_DEBASSERT( j < _outM[R].size() );
            return _outM[R][j];
        }

        /// Returns constant reference to the index of factors of region \a R
        const SmallSet<size_t>& Rfs( size_t R ) const {
            DAI_DEBASSERT( R < _Rfs.size() );
            return _Rfs[R];
        }

        /// Returns reference to the index of factors of region \a R
        SmallSet<size_t>& Rfs( size_t R ) {
            DAI_DEBASSERT( R < _Rfs.size() );
            return _Rfs[R];
        }

        /// Returns constant reference to the index of variables on the boundary of region \a R (\ominus r)
        const SmallSet<size_t>& EXRs( size_t R ) const {
            DAI_DEBASSERT( R < _EXRs.size() );
            return _EXRs[R];
        }

        /// Returns reference to the index of variables on the boundary of region \a R (\ominus r)
        SmallSet<size_t>& EXRs( size_t R ) {
            DAI_DEBASSERT( R < _EXRs.size() );
            return _EXRs[R];
        }

        /// Returns constant reference to the index of variables inside region \a R (r)
        const SmallSet<size_t>& INRs( size_t R ) const {
            DAI_DEBASSERT( R < _INRs.size() );
            return _INRs[R];
        }

        /// Returns reference to the index of variables inside region \a R (r)
        SmallSet<size_t>& INRs( size_t R ) {
            DAI_DEBASSERT( R < _INRs.size() );
            return _INRs[R];
        }

        /// Returns constant reference to the connection structure from region \a R to its \a i'th neighbour
        const Connection& M( size_t R, size_t i  ) const {
            DAI_DEBASSERT(R < _M.size());
            DAI_DEBASSERT(i < _M[R].size());
            return _M[R][i]; }

        /// Returns reference to the connection structure from region \a R to its \a i'th neighbour
        Connection& M( size_t R, size_t i  ) {
            DAI_DEBASSERT(R < _M.size());
            DAI_DEBASSERT(i < _M[R].size());
            return _M[R][i]; }

        /// Returns constant reference to the vector of connection structure from region \a R to all its neighbours
        const std::vector<Connection>& M( size_t R ) const {
            DAI_DEBASSERT(R < _M.size());
            return _M[R]; 
        }

        /// Returns vector of all connections to region \a R
        std::vector<Connection>& M( size_t R ) { return _M[R]; }

        /// Returns the vector of region indices that contain \a i as internal variable
        const std::vector<size_t>& var2CW( size_t i ) const { return _var2CW[i]; }
    //@}

    /// \name Operations
    //@{
        /// Sets up all the regions and messages
        void setRgn( std::vector<SmallSet<size_t> > regions, NeighborType neighbors, bool debugging = false );
    //@}

    /// \name Input/output
    //@{
        /// Reads a cobweb graph from a file
        /** \note Not implemented yet
         */
        virtual void ReadFromFile( const char* /*filename*/ ) {
            DAI_THROW(NOT_IMPLEMENTED);
        }

        /// Writes a cobweb graph to a file
        /** \note Not implemented yet
         */
        virtual void WriteToFile( const char* /*filename*/, size_t /*precision*/=15 ) const {
            DAI_THROW(NOT_IMPLEMENTED);
        }

        /// Writes a cobweb graph to an output stream
        friend std::ostream& operator<< ( std::ostream& /*os*/, const CobwebGraph& /*rg*/ ) {
            DAI_THROW(NOT_IMPLEMENTED);
        }

        /// Formats a cobweb graph as a string
        std::string toString() const {
            std::stringstream ss;
            ss << *this;
            return ss.str();
        }

        /// Writes a cobweb graph to a GraphViz .dot file
        /** \note Not implemented yet
         */
        virtual void printDot( std::ostream& /*os*/ ) const {
            DAI_THROW(NOT_IMPLEMENTED);
        }
    //@}


    protected:
        /// The function to check for partitioning
        bool checkPartition( const std::vector<SmallSet<size_t> >& regions ) const;

        /// Helper function that sets the regions containing each variable using the values in _INRs
        void setVar2CW();

        /// Setup the _INRs, _EXRs and all the factor indices (e.g. Rifs)
        void setExtnFact();

        /// Helper function that setups the msgs (_M, _outM) using the values in _INRs and _EXRs
        void setMSGs( NeighborType neighbors );

        /// Sets _cn
        void setCountingNumbers( bool debugging = false );

        /// For the given set of variables for each region, removes the regions that are non-maximal
        void eraseNonMaximal( std::vector<SmallSet<size_t> >& regions );
};


} // end of namespace dai


#endif
