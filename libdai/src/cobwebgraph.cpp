/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2012, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


#include <dai/cobwebgraph.h>
#include <algorithm>
#include <map>
#include <set>


namespace dai {


using namespace std;


void CobwebGraph::setRgn( vector<SmallSet<size_t> > partition, NeighborType neighbors, bool debugging ) {
    if( checkPartition( partition ) )
        isPartition = true;
    else
        isPartition = false;

    eraseNonMaximal( partition );
    if( debugging )
        cerr << "isPartition:" << isPartition << endl;
    _INRs = partition;
    if( debugging )
        cerr << "setting external factors and regional factors" << endl;
    setExtnFact();
    if( debugging )
        cerr << "setting msgs" << endl;
    setMSGs( neighbors );
//    if(!isPartition)
    if( debugging )
        cerr << "setting the counting numbers" << endl;
    if( !isPartition )
        setCountingNumbers( debugging );
    if( debugging )
        cerr << "setting var2dr" << endl;
    setVar2CW();
}


void CobwebGraph::eraseNonMaximal( vector<SmallSet<size_t> >& partition ) {
    vector<SmallSet<size_t> >::iterator del = partition.begin();
    while( del != partition.end() ){
        bool isNonMax = false;
        for( vector<SmallSet<size_t> >::iterator it = partition.begin(); it != partition.end(); it++ )
            if( (*it) >> (*del) && it != del ) {
                isNonMax = true;
                break;
            }
        if( isNonMax )
            del = partition.erase( del );
        else
            del++;
    }
}


void CobwebGraph::setCountingNumbers( bool debugging ) {
    // should handle the case when one region in the msg-region graph is subset of another or even to top regions are the same
    _cn.reserve( nrCWs() );
    for( size_t R = 0; R < _INRs.size(); R++ ) {
        vector<VarSet> topR;
        // finding the intersection of messages
        SmallSet<VarSet> betas;
        for( size_t m1 = 0; m1 < M(R).size(); m1++ ) {
            topR.push_back( M(R,m1).msg.vars() );
            for( size_t m2 = 0; m2 < M(R).size(); m2++ ) {
                if( m1 == m2 )
                    continue;
                VarSet b = M(R,m1).msg.vars() & M(R,m2).msg.vars();
                if( b.size() == 0 )
                    continue;
                betas.insert( b );
            }
        }
        if( debugging )
            cerr << "topR: " << topR << endl;

        // finding the intersections of intersections
        bool somechange = true;
        while( somechange ) {
            SmallSet<VarSet> newbetas;
            somechange = false;
            bforeach( const VarSet &b1, betas ) {
                bforeach( const VarSet &b2, betas ) {
                    if( b1 == b2 )
                        continue;
                    VarSet b3 = b1 & b2;
                    if( betas.contains(b3) || b3.size() == 0 )
                        continue;
                    newbetas |= b3;
                    somechange = true;
                }
            }
            betas |= newbetas;
        }
        if( debugging )
            cerr << "betas: " << betas << endl;

        // set the counting number
        _cn.push_back( map<VarSet, pair<int, vector<size_t> > >() );
        // adding sub-regions of every message
        for( size_t i = 0; i < topR.size(); i++ ) {
            M(R,i).subregions.clear();
            bforeach( const VarSet& b, betas )
                if( b << topR[i] )
                    M(R,i).subregions.push_back( b );
        }
        SmallSet<VarSet> subVisited;
        SmallSet<VarSet> topRSet( topR.begin(), topR.end(), topR.size() );
        // looks to see if all parents of a sub-region got their counting number
        while( !betas.empty() ) {
            bforeach( const VarSet &beta, betas ) {
                bool allparentsset = true;
                bforeach( const VarSet &beta2, betas )
                    if( beta2 >> beta && beta2 != beta ) {
                        allparentsset = false;
                        break;
                    }
                if( allparentsset ) {
                    // the first in the pair is cn and the second the index of top regions containing it
                    _cn[R][beta] = make_pair( 1, vector<size_t>() );
                    for( size_t TR = 0; TR < topR.size(); TR++ ) {
                        if( topR[TR] >> beta ) {
                            _cn[R][beta].first--;
                            _cn[R][beta].second.push_back( TR );
                        }
                    }
                    bforeach( const VarSet& possibleparent, subVisited )
                        if( possibleparent >> beta )
                            _cn[R][beta].first -= _cn[R][possibleparent].first;

                    if( debugging )
                        cerr << "cn[" << R << "][" << beta << "] <- " << _cn[R][beta] << endl;
                    subVisited.insert( beta );
                    betas /= beta;
                    break; // since betas has changed we need to enter the loop again
                }
            }
        }
//        if( debugging )
//            cerr << "cn[" << R << "] " << _cn[R] << endl;
    }
}


bool CobwebGraph::checkPartition( const vector<SmallSet<size_t> >& regions ) const {
    for( size_t i = 0; i < regions.size(); i++ )
        for( size_t j = 0; j < regions.size(); j++ ) {
            if( j == i )
                continue;
            if( regions[i].intersects( regions[j] ) )
                return false;
        }
    return true;
}


void CobwebGraph::setVar2CW() {
    _var2CW.clear();
    _var2CW.resize( nrVars() );
    for( size_t R = 0; R < _INRs.size(); R++ ) {
        bforeach( size_t i, _INRs[R] )
            _var2CW[i].push_back( R );
    }
}


// assumes availablity of externals and internals
void CobwebGraph::setMSGs( NeighborType neighbors ) {
    // initialize
    _outM.clear();
    _outM.reserve( nrCWs() );
    for( size_t R = 0; R < nrCWs(); R++ )
        _outM.push_back( vector<VarSet>() );
    _M.clear();
    _M.reserve( _INRs.size() );
    for( size_t R = 0; R < _INRs.size(); R++ ) {
        _M.push_back( vector<Connection>() );
        for( size_t R2 = 0; R2 < _INRs.size(); R2++ ) {
            if( R2 == R )
                continue;
            SmallSet<size_t> tmpSet = _EXRs[R] & (_INRs[R2]); // if R2 has boundary of R1 inside it
            if( tmpSet.size() == 0 )
                continue;

            Connection tmp;
            tmp.fc = (_Rfs[R] & _Rfs[R2]).elements();
            tmp.my = R;
            tmp.iter = _M[R].size();
            tmp.his = R2;
            tmp.dual = _outM[R2].size();
            VarSet tmpVarSet;
            bforeach( const size_t &i, tmpSet )
                tmpVarSet.insert( var(i) );
            tmp.msg = Factor( tmpVarSet, (Real)1.0 ); //*
            _outM[R2].push_back( tmpVarSet ); //*
            // tmp.varinds = tmpSet.elements();
            _M[R].push_back( tmp );
        }
        if( neighbors == NeighborType::CLOSEST || neighbors == NeighborType::TOP ) {
            bool foundone = true;
            while( foundone ) {
                for( size_t Rfirst = 0; Rfirst < _M[R].size(); Rfirst++ ) {
                    foundone = false;
                    for( size_t Rsec = 0; Rsec < _M[R].size(); Rsec++ ) {
                        if( Rfirst == Rsec )
                            continue;
                        if( _M[R][Rfirst].msg.vars() >> _M[R][Rsec].msg.vars() ) { // one message is subseteq of the other
                            // only if one exclusively contains the other remove it
                            if( ((neighbors == NeighborType::TOP || neighbors == NeighborType::CLOSEST) && _M[R][Rfirst].msg.vars().size() > _M[R][Rsec].msg.vars().size()) ) {
                                _M[R].erase( _M[R].begin() + Rsec );
                                foundone = true;
                                break;
                            } else if( neighbors == NeighborType::CLOSEST && // if both have the same size then in case of CLOSEST delete the second one if it has smaller inner region overlap
                                    (INRs(_M[R][Rfirst].his) & INRs(R)).size() > (INRs(_M[R][Rsec].his) & INRs(R)).size() ) {
                                _M[R].erase( _M[R].begin() + Rsec );
                                foundone = true;
                                break;
                            }
                        }
                    }
                    if( foundone )
                        break;
                }
            }
        }
    }
}


void CobwebGraph::setExtnFact() {
    // finds and sets all the externals using _INRs
    SmallSet<size_t> allpots;
    _EXRs.clear();
    _Rfs.clear();
    _EXRs.reserve( _INRs.size() );
    _Rfs.reserve( _INRs.size() );
    _Rifs.reserve( _INRs.size() );
    _Rxfs.reserve( _INRs.size() );
    for( size_t R = 0; R < _INRs.size(); R++ ) {
        _EXRs.push_back( SmallSet<size_t>() );
        _Rfs.push_back( SmallSet<size_t>() );
        _Rifs.push_back( SmallSet<size_t>() );
        _Rxfs.push_back( SmallSet<size_t>() );
        bforeach( const size_t &i, _INRs[R] ) {
            _EXRs[R] |= deltai(i) / _INRs[R];
            bforeach( const Neighbor &I, nbV(i) ) {
                _Rfs[R] |= I.node;
                if( factor(I).vars() << inds2vars(INRs(R).elements()) )
                    _Rifs[R] |= I.node;
                else
                    _Rxfs[R] |= I.node;
            }
        }
    }
}


} // end of namespace dai
