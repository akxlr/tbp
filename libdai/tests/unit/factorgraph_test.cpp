/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2011, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


#include <dai/bipgraph.h>
#include <dai/factorgraph.h>
#include <vector>
#include <strstream>


using namespace dai;


const double tol = 1e-8;


#define BOOST_TEST_MODULE FactorGraphTest


#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


BOOST_AUTO_TEST_CASE( ConstructorsTest ) {
    FactorGraph G;
    BOOST_CHECK_EQUAL( G.vars(), std::vector<Var>() );
    BOOST_CHECK_EQUAL( G.factors(), std::vector<Factor>() );

    std::vector<Factor> facs;
    facs.push_back( Factor( VarSet( Var(0, 2), Var(1, 2) ) ) );
    facs.push_back( Factor( VarSet( Var(0, 2), Var(2, 2) ) ) );
    facs.push_back( Factor( VarSet( Var(1, 2), Var(2, 2) ) ) );
    facs.push_back( Factor( VarSet( Var(1, 2) ) ) );
    std::vector<Var> vars;
    vars.push_back( Var( 0, 2 ) );
    vars.push_back( Var( 1, 2 ) );
    vars.push_back( Var( 2, 2 ) );

    FactorGraph G1( facs );
    BOOST_CHECK_EQUAL( G1.vars(), vars );
    BOOST_CHECK_EQUAL( G1.factors(), facs );

    FactorGraph G2( facs.begin(), facs.end(), vars.begin(), vars.end(), facs.size(), vars.size() ); 
    BOOST_CHECK_EQUAL( G2.vars(), vars );
    BOOST_CHECK_EQUAL( G2.factors(), facs );

    FactorGraph *G3 = G2.clone();
    BOOST_CHECK_EQUAL( G3->vars(), vars );
    BOOST_CHECK_EQUAL( G3->factors(), facs );
    delete G3;

    FactorGraph G4 = G2;
    BOOST_CHECK_EQUAL( G4.vars(), vars );
    BOOST_CHECK_EQUAL( G4.factors(), facs );

    FactorGraph G5( G2 );
    BOOST_CHECK_EQUAL( G5.vars(), vars );
    BOOST_CHECK_EQUAL( G5.factors(), facs );
}


BOOST_AUTO_TEST_CASE( AccMutTest ) {
    std::vector<Factor> facs;
    facs.push_back( Factor( VarSet( Var(0, 2), Var(1, 2) ) ) );
    facs.push_back( Factor( VarSet( Var(0, 2), Var(2, 2) ) ) );
    facs.push_back( Factor( VarSet( Var(1, 2), Var(2, 2) ) ) );
    facs.push_back( Factor( VarSet( Var(1, 2) ) ) );
    std::vector<Var> vars;
    vars.push_back( Var( 0, 2 ) );
    vars.push_back( Var( 1, 2 ) );
    vars.push_back( Var( 2, 2 ) );

    FactorGraph G( facs );
    BOOST_CHECK_EQUAL( G.var(0), Var(0, 2) );
    BOOST_CHECK_EQUAL( G.var(1), Var(1, 2) );
    BOOST_CHECK_EQUAL( G.var(2), Var(2, 2) );
    BOOST_CHECK_EQUAL( G.vars(), vars );
    BOOST_CHECK_EQUAL( G.factor(0), facs[0] );
    BOOST_CHECK_EQUAL( G.factor(1), facs[1] );
    BOOST_CHECK_EQUAL( G.factor(2), facs[2] );
    BOOST_CHECK_EQUAL( G.factor(3), facs[3] );
    BOOST_CHECK_EQUAL( G.factors(), facs );
    BOOST_CHECK_EQUAL( G.nbV(0).size(), 2 );
    BOOST_CHECK_EQUAL( G.nbV(0,0), 0 );
    BOOST_CHECK_EQUAL( G.nbV(0,1), 1 );
    BOOST_CHECK_EQUAL( G.nbV(1).size(), 3 );
    BOOST_CHECK_EQUAL( G.nbV(1,0), 0 );
    BOOST_CHECK_EQUAL( G.nbV(1,1), 2 );
    BOOST_CHECK_EQUAL( G.nbV(1,2), 3 );
    BOOST_CHECK_EQUAL( G.nbV(0).size(), 2 );
    BOOST_CHECK_EQUAL( G.nbV(2,0), 1 );
    BOOST_CHECK_EQUAL( G.nbV(2,1), 2 );
    BOOST_CHECK_EQUAL( G.nbF(0).size(), 2 );
    BOOST_CHECK_EQUAL( G.nbF(0,0), 0 );
    BOOST_CHECK_EQUAL( G.nbF(0,1), 1 );
    BOOST_CHECK_EQUAL( G.nbF(1).size(), 2 );
    BOOST_CHECK_EQUAL( G.nbF(1,0), 0 );
    BOOST_CHECK_EQUAL( G.nbF(1,1), 2 );
    BOOST_CHECK_EQUAL( G.nbF(2).size(), 2 );
    BOOST_CHECK_EQUAL( G.nbF(2,0), 1 );
    BOOST_CHECK_EQUAL( G.nbF(2,1), 2 );
    BOOST_CHECK_EQUAL( G.nbF(3).size(), 1 );
    BOOST_CHECK_EQUAL( G.nbF(3,0), 1 );
}


BOOST_AUTO_TEST_CASE( QueriesTest ) {
    Var v0( 0, 2 );
    Var v1( 1, 2 );
    Var v2( 2, 2 );
    VarSet v01( v0, v1 );
    VarSet v02( v0, v2 );
    VarSet v12( v1, v2 );
    VarSet v012 = v01 | v2;

    FactorGraph G0;
    BOOST_CHECK_EQUAL( G0.nrVars(), 0 );
    BOOST_CHECK_EQUAL( G0.nrFactors(), 0 );
    BOOST_CHECK_EQUAL( G0.nrEdges(), 0 );
    BOOST_CHECK_THROW( G0.findVar( v0 ), Exception );
    BOOST_CHECK_THROW( G0.findVars( v01 ), Exception );
    BOOST_CHECK_THROW( G0.findFactor( v01 ), Exception );
#ifdef DAI_DBEUG
    BOOST_CHECK_THROW( G0.delta( 0 ), Exception );
    BOOST_CHECK_THROW( G0.Delta( 0 ), Exception );
    BOOST_CHECK_THROW( G0.delta( v0 ), Exception );
    BOOST_CHECK_THROW( G0.Delta( v0 ), Exception );
#endif
    BOOST_CHECK( G0.isConnected() );
    BOOST_CHECK( G0.isTree() );
    BOOST_CHECK( G0.isBinary() );
    BOOST_CHECK( G0.isPairwise() );
    BOOST_CHECK( G0.MarkovGraph() == GraphAL() );
    BOOST_CHECK( G0.bipGraph() == BipartiteGraph() );
    BOOST_CHECK_EQUAL( G0.maximalFactorDomains().size(), 1 );
    BOOST_CHECK_CLOSE( G0.logScore( std::vector<size_t>() ), (Real)0.0, tol );

    std::vector<Factor> facs;
    facs.push_back( Factor( v01 ) );
    facs.push_back( Factor( v12 ) );
    facs.push_back( Factor( v1 ) );
    std::vector<Var> vars;
    vars.push_back( v0 );
    vars.push_back( v1 );
    vars.push_back( v2 );
    GraphAL H(3);
    H.addEdge( 0, 1 );
    H.addEdge( 1, 2 );
    BipartiteGraph K(3, 3);
    K.addEdge( 0, 0 );
    K.addEdge( 1, 0 );
    K.addEdge( 1, 1 );
    K.addEdge( 2, 1 );
    K.addEdge( 1, 2 );

    FactorGraph G1( facs );
    BOOST_CHECK_EQUAL( G1.nrVars(), 3 );
    BOOST_CHECK_EQUAL( G1.nrFactors(), 3 );
    BOOST_CHECK_EQUAL( G1.nrEdges(), 5 );
    BOOST_CHECK_EQUAL( G1.findVar( v0 ), 0 );
    BOOST_CHECK_EQUAL( G1.findVar( v1 ), 1 );
    BOOST_CHECK_EQUAL( G1.findVar( v2 ), 2 );
    BOOST_CHECK_EQUAL( G1.findVars( v01 ), SmallSet<size_t>( 0, 1 ) );
    BOOST_CHECK_EQUAL( G1.findVars( v02 ), SmallSet<size_t>( 0, 2 ) );
    BOOST_CHECK_EQUAL( G1.findVars( v12 ), SmallSet<size_t>( 1, 2 ) );
    BOOST_CHECK_EQUAL( G1.findFactor( v01 ), 0 );
    BOOST_CHECK_EQUAL( G1.findFactor( v12 ), 1 );
    BOOST_CHECK_EQUAL( G1.findFactor( v1 ), 2 );
    BOOST_CHECK_THROW( G1.findFactor( v02 ), Exception );
    BOOST_CHECK_EQUAL( G1.delta( 0 ), v1 );
    BOOST_CHECK_EQUAL( G1.delta( 1 ), v02 );
    BOOST_CHECK_EQUAL( G1.delta( 2 ), v1 );
    BOOST_CHECK_EQUAL( G1.Delta( 0 ), v01 );
    BOOST_CHECK_EQUAL( G1.Delta( 1 ), v012 );
    BOOST_CHECK_EQUAL( G1.Delta( 2 ), v12 );
    BOOST_CHECK_EQUAL( G1.delta( v0 ), v1 );
    BOOST_CHECK_EQUAL( G1.delta( v1 ), v02 );
    BOOST_CHECK_EQUAL( G1.delta( v2 ), v1 );
    BOOST_CHECK_EQUAL( G1.delta( v01 ), v2 );
    BOOST_CHECK_EQUAL( G1.delta( v02 ), v1 );
    BOOST_CHECK_EQUAL( G1.delta( v12 ), v0 );
    BOOST_CHECK_EQUAL( G1.delta( v012 ), VarSet() );
    BOOST_CHECK_EQUAL( G1.Delta( v0 ), v01 );
    BOOST_CHECK_EQUAL( G1.Delta( v1 ), v012 );
    BOOST_CHECK_EQUAL( G1.Delta( v2 ), v12 );
    BOOST_CHECK_EQUAL( G1.Delta( v01 ), v012 );
    BOOST_CHECK_EQUAL( G1.Delta( v02 ), v012 );
    BOOST_CHECK_EQUAL( G1.Delta( v12 ), v012 );
    BOOST_CHECK_EQUAL( G1.Delta( v012 ), v012 );
    BOOST_CHECK( G1.isConnected() );
    BOOST_CHECK( G1.isTree() );
    BOOST_CHECK( G1.isBinary() );
    BOOST_CHECK( G1.isPairwise() );
    BOOST_CHECK( G1.MarkovGraph() == H );
    BOOST_CHECK( G1.bipGraph() == K );
    BOOST_CHECK(  G1.isMaximal( 0 ) );
    BOOST_CHECK(  G1.isMaximal( 1 ) );
    BOOST_CHECK( !G1.isMaximal( 2 ) );
    BOOST_CHECK_EQUAL( G1.maximalFactor( 0 ), 0 );
    BOOST_CHECK_EQUAL( G1.maximalFactor( 1 ), 1 );
    BOOST_CHECK_EQUAL( G1.maximalFactor( 2 ), 0 );
    BOOST_CHECK_EQUAL( G1.maximalFactorDomains().size(), 2 );
    BOOST_CHECK_EQUAL( G1.maximalFactorDomains()[0], v01 );
    BOOST_CHECK_EQUAL( G1.maximalFactorDomains()[1], v12 );
    BOOST_CHECK_CLOSE( G1.logScore( std::vector<size_t>(3,0) ), -dai::log((Real)32.0), tol ); 

    facs.push_back( Factor( v02 ) );
    H.addEdge( 0, 2 );
    K.addNode2();
    K.addEdge( 0, 3 );
    K.addEdge( 2, 3 );
    FactorGraph G2( facs );
    BOOST_CHECK_EQUAL( G2.nrVars(), 3 );
    BOOST_CHECK_EQUAL( G2.nrFactors(), 4 );
    BOOST_CHECK_EQUAL( G2.nrEdges(), 7 );
    BOOST_CHECK_EQUAL( G2.findVar( v0 ), 0 );
    BOOST_CHECK_EQUAL( G2.findVar( v1 ), 1 );
    BOOST_CHECK_EQUAL( G2.findVar( v2 ), 2 );
    BOOST_CHECK_EQUAL( G2.findVars( v01 ), SmallSet<size_t>( 0, 1 ) );
    BOOST_CHECK_EQUAL( G2.findVars( v02 ), SmallSet<size_t>( 0, 2 ) );
    BOOST_CHECK_EQUAL( G2.findVars( v12 ), SmallSet<size_t>( 1, 2 ) );
    BOOST_CHECK_EQUAL( G2.findFactor( v01 ), 0 );
    BOOST_CHECK_EQUAL( G2.findFactor( v12 ), 1 );
    BOOST_CHECK_EQUAL( G2.findFactor( v1 ), 2 );
    BOOST_CHECK_EQUAL( G2.findFactor( v02 ), 3 );
    BOOST_CHECK_EQUAL( G2.delta( 0 ), v12 );
    BOOST_CHECK_EQUAL( G2.delta( 1 ), v02 );
    BOOST_CHECK_EQUAL( G2.delta( 2 ), v01 );
    BOOST_CHECK_EQUAL( G2.Delta( 0 ), v012 );
    BOOST_CHECK_EQUAL( G2.Delta( 1 ), v012 );
    BOOST_CHECK_EQUAL( G2.Delta( 2 ), v012 );
    BOOST_CHECK( G2.isConnected() );
    BOOST_CHECK( !G2.isTree() );
    BOOST_CHECK( G2.isBinary() );
    BOOST_CHECK( G2.isPairwise() );
    BOOST_CHECK( G2.MarkovGraph() == H );
    BOOST_CHECK( G2.bipGraph() == K );
    BOOST_CHECK(  G2.isMaximal( 0 ) );
    BOOST_CHECK(  G2.isMaximal( 1 ) );
    BOOST_CHECK( !G2.isMaximal( 2 ) );
    BOOST_CHECK(  G2.isMaximal( 3 ) );
    BOOST_CHECK_EQUAL( G2.maximalFactor( 0 ), 0 );
    BOOST_CHECK_EQUAL( G2.maximalFactor( 1 ), 1 );
    BOOST_CHECK_EQUAL( G2.maximalFactor( 2 ), 0 );
    BOOST_CHECK_EQUAL( G2.maximalFactor( 3 ), 3 );
    BOOST_CHECK_EQUAL( G2.maximalFactorDomains().size(), 3 );
    BOOST_CHECK_EQUAL( G2.maximalFactorDomains()[0], v01 );
    BOOST_CHECK_EQUAL( G2.maximalFactorDomains()[1], v12 );
    BOOST_CHECK_EQUAL( G2.maximalFactorDomains()[2], v02 );
    BOOST_CHECK_CLOSE( G2.logScore( std::vector<size_t>(3,0) ), -dai::log((Real)128.0), tol );

    Var v3( 3, 3 );
    VarSet v03( v0, v3 );
    VarSet v13( v1, v3 );
    VarSet v23( v2, v3 );
    VarSet v013 = v01 | v3;
    VarSet v023 = v02 | v3;
    VarSet v123 = v12 | v3;
    VarSet v0123 = v012 | v3;
    vars.push_back( v3 );
    facs.push_back( Factor( v3 ) );
    H.addNode();
    K.addNode1();
    K.addNode2();
    K.addEdge( 3, 4 );
    FactorGraph G3( facs );
    BOOST_CHECK_EQUAL( G3.nrVars(), 4 );
    BOOST_CHECK_EQUAL( G3.nrFactors(), 5 );
    BOOST_CHECK_EQUAL( G3.nrEdges(), 8 );
    BOOST_CHECK_EQUAL( G3.findVar( v0 ), 0 );
    BOOST_CHECK_EQUAL( G3.findVar( v1 ), 1 );
    BOOST_CHECK_EQUAL( G3.findVar( v2 ), 2 );
    BOOST_CHECK_EQUAL( G3.findVar( v3 ), 3 );
    BOOST_CHECK_EQUAL( G3.findVars( v01 ), SmallSet<size_t>( 0, 1 ) );
    BOOST_CHECK_EQUAL( G3.findVars( v02 ), SmallSet<size_t>( 0, 2 ) );
    BOOST_CHECK_EQUAL( G3.findVars( v12 ), SmallSet<size_t>( 1, 2 ) );
    BOOST_CHECK_EQUAL( G3.findFactor( v01 ), 0 );
    BOOST_CHECK_EQUAL( G3.findFactor( v12 ), 1 );
    BOOST_CHECK_EQUAL( G3.findFactor( v1 ), 2 );
    BOOST_CHECK_EQUAL( G3.findFactor( v02 ), 3 );
    BOOST_CHECK_EQUAL( G3.findFactor( v3 ), 4 );
    BOOST_CHECK_THROW( G3.findFactor( v23 ), Exception );
    BOOST_CHECK_EQUAL( G3.delta( 0 ), v12 );
    BOOST_CHECK_EQUAL( G3.delta( 1 ), v02 );
    BOOST_CHECK_EQUAL( G3.delta( 2 ), v01 );
    BOOST_CHECK_EQUAL( G3.delta( 3 ), VarSet() );
    BOOST_CHECK_EQUAL( G3.Delta( 0 ), v012 );
    BOOST_CHECK_EQUAL( G3.Delta( 1 ), v012 );
    BOOST_CHECK_EQUAL( G3.Delta( 2 ), v012 );
    BOOST_CHECK_EQUAL( G3.Delta( 3 ), v3 );
    BOOST_CHECK( !G3.isConnected() );
    BOOST_CHECK( !G3.isTree() );
    BOOST_CHECK( !G3.isBinary() );
    BOOST_CHECK( G3.isPairwise() );
    BOOST_CHECK( G3.MarkovGraph() == H );
    BOOST_CHECK( G3.bipGraph() == K );
    BOOST_CHECK(  G3.isMaximal( 0 ) );
    BOOST_CHECK(  G3.isMaximal( 1 ) );
    BOOST_CHECK( !G3.isMaximal( 2 ) );
    BOOST_CHECK(  G3.isMaximal( 3 ) );
    BOOST_CHECK(  G3.isMaximal( 4 ) );
    BOOST_CHECK_EQUAL( G3.maximalFactor( 0 ), 0 );
    BOOST_CHECK_EQUAL( G3.maximalFactor( 1 ), 1 );
    BOOST_CHECK_EQUAL( G3.maximalFactor( 2 ), 0 );
    BOOST_CHECK_EQUAL( G3.maximalFactor( 3 ), 3 );
    BOOST_CHECK_EQUAL( G3.maximalFactor( 4 ), 4 );
    BOOST_CHECK_EQUAL( G3.maximalFactorDomains().size(), 4 );
    BOOST_CHECK_EQUAL( G3.maximalFactorDomains()[0], v01 );
    BOOST_CHECK_EQUAL( G3.maximalFactorDomains()[1], v12 );
    BOOST_CHECK_EQUAL( G3.maximalFactorDomains()[2], v02 );
    BOOST_CHECK_EQUAL( G3.maximalFactorDomains()[3], v3 );
    BOOST_CHECK_CLOSE( G3.logScore( std::vector<size_t>(4,0) ), -dai::log((Real)384.0), tol );

    facs.push_back( Factor( v123 ) );
    H.addEdge( 1, 3 );
    H.addEdge( 2, 3 );
    K.addNode2();
    K.addEdge( 1, 5 );
    K.addEdge( 2, 5 );
    K.addEdge( 3, 5 );
    FactorGraph G4( facs );
    BOOST_CHECK_EQUAL( G4.nrVars(), 4 );
    BOOST_CHECK_EQUAL( G4.nrFactors(), 6 );
    BOOST_CHECK_EQUAL( G4.nrEdges(), 11 );
    BOOST_CHECK_EQUAL( G4.findVar( v0 ), 0 );
    BOOST_CHECK_EQUAL( G4.findVar( v1 ), 1 );
    BOOST_CHECK_EQUAL( G4.findVar( v2 ), 2 );
    BOOST_CHECK_EQUAL( G4.findVar( v3 ), 3 );
    BOOST_CHECK_EQUAL( G4.findVars( v01 ), SmallSet<size_t>( 0, 1 ) );
    BOOST_CHECK_EQUAL( G4.findVars( v02 ), SmallSet<size_t>( 0, 2 ) );
    BOOST_CHECK_EQUAL( G4.findVars( v12 ), SmallSet<size_t>( 1, 2 ) );
    BOOST_CHECK_EQUAL( G4.findFactor( v01 ), 0 );
    BOOST_CHECK_EQUAL( G4.findFactor( v12 ), 1 );
    BOOST_CHECK_EQUAL( G4.findFactor( v1 ), 2 );
    BOOST_CHECK_EQUAL( G4.findFactor( v02 ), 3 );
    BOOST_CHECK_EQUAL( G4.findFactor( v3 ), 4 );
    BOOST_CHECK_EQUAL( G4.findFactor( v123 ), 5 );
    BOOST_CHECK_THROW( G4.findFactor( v23 ), Exception );
    BOOST_CHECK_EQUAL( G4.delta( 0 ), v12 );
    BOOST_CHECK_EQUAL( G4.delta( 1 ), v023 );
    BOOST_CHECK_EQUAL( G4.delta( 2 ), v013 );
    BOOST_CHECK_EQUAL( G4.delta( 3 ), v12 );
    BOOST_CHECK_EQUAL( G4.Delta( 0 ), v012 );
    BOOST_CHECK_EQUAL( G4.Delta( 1 ), v0123 );
    BOOST_CHECK_EQUAL( G4.Delta( 2 ), v0123 );
    BOOST_CHECK_EQUAL( G4.Delta( 3 ), v123 );
    BOOST_CHECK( G4.isConnected() );
    BOOST_CHECK( !G4.isTree() );
    BOOST_CHECK( !G4.isBinary() );
    BOOST_CHECK( !G4.isPairwise() );
    BOOST_CHECK( G4.MarkovGraph() == H );
    BOOST_CHECK( G4.bipGraph() == K );
    BOOST_CHECK(  G4.isMaximal( 0 ) );
    BOOST_CHECK( !G4.isMaximal( 1 ) );
    BOOST_CHECK( !G4.isMaximal( 2 ) );
    BOOST_CHECK(  G4.isMaximal( 3 ) );
    BOOST_CHECK( !G4.isMaximal( 4 ) );
    BOOST_CHECK(  G4.isMaximal( 5 ) );
    BOOST_CHECK_EQUAL( G4.maximalFactor( 0 ), 0 );
    BOOST_CHECK_EQUAL( G4.maximalFactor( 1 ), 5 );
    BOOST_CHECK_EQUAL( G4.maximalFactor( 2 ), 0 );
    BOOST_CHECK_EQUAL( G4.maximalFactor( 3 ), 3 );
    BOOST_CHECK_EQUAL( G4.maximalFactor( 4 ), 5 );
    BOOST_CHECK_EQUAL( G4.maximalFactor( 5 ), 5 );
    BOOST_CHECK_EQUAL( G4.maximalFactorDomains().size(), 3 );
    BOOST_CHECK_EQUAL( G4.maximalFactorDomains()[0], v01 );
    BOOST_CHECK_EQUAL( G4.maximalFactorDomains()[1], v02 );
    BOOST_CHECK_EQUAL( G4.maximalFactorDomains()[2], v123 );
    BOOST_CHECK_CLOSE( G4.logScore( std::vector<size_t>(4,0) ), -dai::log((Real)4608.0), tol );
}


BOOST_AUTO_TEST_CASE( BackupRestoreTest ) {
    Var v0( 0, 2 );
    Var v1( 1, 2 );
    Var v2( 2, 2 );
    VarSet v01( v0, v1 );
    VarSet v02( v0, v2 );
    VarSet v12( v1, v2 );
    VarSet v012 = v01 | v2;

    std::vector<Factor> facs;
    facs.push_back( Factor( v01 ) );
    facs.push_back( Factor( v12 ) );
    facs.push_back( Factor( v1 ) );
    std::vector<Var> vars;
    vars.push_back( v0 );
    vars.push_back( v1 );
    vars.push_back( v2 );

    FactorGraph G( facs );
    FactorGraph Gorg( G );

    BOOST_CHECK_THROW( G.setFactor( 0, Factor( v0 ), false ), Exception );
    G.setFactor( 0, Factor( v01, 2.0 ), false );
    BOOST_CHECK_THROW( G.restoreFactor( 0 ), Exception );
    G.setFactor( 0, Factor( v01, 3.0 ), true );
    G.restoreFactor( 0 );
    BOOST_CHECK_EQUAL( G.factor( 0 )[0], 2.0 );
    G.setFactor( 0, Gorg.factor( 0 ), false );
    G.backupFactor( 0 );
    BOOST_CHECK_EQUAL( G.factor(0), Gorg.factor(0) );
    G.setFactor( 0, Factor( v01, 2.0 ), false );
    BOOST_CHECK_EQUAL( G.factor( 0 )[0], 2.0 );
    G.restoreFactor( 0 );
    BOOST_CHECK_EQUAL( G.factor(0), Gorg.factor(0) );

    std::map<size_t, Factor> fs;
    fs[0] = Factor( v01, 3.0 );
    fs[2] = Factor( v1, 2.0 );
    G.setFactors( fs, false );
    BOOST_CHECK_EQUAL( G.factor(0), fs[0] );
    BOOST_CHECK_EQUAL( G.factor(2), fs[2] );
    BOOST_CHECK_EQUAL( G.factor(1), Gorg.factor(1) );
    G.restoreFactors();
    BOOST_CHECK_EQUAL( G.factor(0), fs[0] );
    BOOST_CHECK_EQUAL( G.factor(2), fs[2] );
    BOOST_CHECK_EQUAL( G.factor(1), Gorg.factor(1) );
    G = Gorg;
    BOOST_CHECK_EQUAL( G.factor(0), Gorg.factor(0) );
    BOOST_CHECK_EQUAL( G.factor(1), Gorg.factor(1) );
    BOOST_CHECK_EQUAL( G.factor(2), Gorg.factor(2) );
    G.setFactors( fs, true );
    BOOST_CHECK_EQUAL( G.factor(0), fs[0] );
    BOOST_CHECK_EQUAL( G.factor(2), fs[2] );
    BOOST_CHECK_EQUAL( G.factor(1), Gorg.factor(1) );
    G.restoreFactors();
    BOOST_CHECK_EQUAL( G.factor(0), Gorg.factor(0) );
    BOOST_CHECK_EQUAL( G.factor(1), Gorg.factor(1) );
    BOOST_CHECK_EQUAL( G.factor(2), Gorg.factor(2) );
    std::set<size_t> fsind;
    fsind.insert( 0 );
    fsind.insert( 2 );
    G.backupFactors( fsind );
    G.setFactors( fs, false );
    BOOST_CHECK_EQUAL( G.factor(0), fs[0] );
    BOOST_CHECK_EQUAL( G.factor(2), fs[2] );
    BOOST_CHECK_EQUAL( G.factor(1), Gorg.factor(1) );
    G.restoreFactors();
    BOOST_CHECK_EQUAL( G.factor(0), Gorg.factor(0) );
    BOOST_CHECK_EQUAL( G.factor(1), Gorg.factor(1) );
    BOOST_CHECK_EQUAL( G.factor(2), Gorg.factor(2) );

    G.backupFactors( v2 );
    G.setFactor( 1, Factor(v12, 5.0) );
    BOOST_CHECK_EQUAL( G.factor(0), Gorg.factor(0) );
    BOOST_CHECK_EQUAL( G.factor(1), Factor(v12, 5.0) );
    BOOST_CHECK_EQUAL( G.factor(2), Gorg.factor(2) );
    G.restoreFactors( v2 );
    BOOST_CHECK_EQUAL( G.factor(0), Gorg.factor(0) );
    BOOST_CHECK_EQUAL( G.factor(1), Gorg.factor(1) );
    BOOST_CHECK_EQUAL( G.factor(2), Gorg.factor(2) );

    G.backupFactors( v1 );
    fs[1] = Factor( v12, 5.0 );
    G.setFactors( fs, false );
    BOOST_CHECK_EQUAL( G.factor(0), fs[0] );
    BOOST_CHECK_EQUAL( G.factor(1), fs[1] );
    BOOST_CHECK_EQUAL( G.factor(2), fs[2] );
    G.restoreFactors();
    BOOST_CHECK_EQUAL( G.factor(0), Gorg.factor(0) );
    BOOST_CHECK_EQUAL( G.factor(1), Gorg.factor(1) );
    BOOST_CHECK_EQUAL( G.factor(2), Gorg.factor(2) );
    G.setFactors( fs, true );
    BOOST_CHECK_EQUAL( G.factor(0), fs[0] );
    BOOST_CHECK_EQUAL( G.factor(1), fs[1] );
    BOOST_CHECK_EQUAL( G.factor(2), fs[2] );
    G.restoreFactors( v1 );
    BOOST_CHECK_EQUAL( G.factor(0), Gorg.factor(0) );
    BOOST_CHECK_EQUAL( G.factor(1), Gorg.factor(1) );
    BOOST_CHECK_EQUAL( G.factor(2), Gorg.factor(2) );
}


BOOST_AUTO_TEST_CASE( TransformationsTest ) {
    Var v0( 0, 2 );
    Var v1( 1, 2 );
    Var v2( 2, 2 );
    VarSet v01( v0, v1 );
    VarSet v02( v0, v2 );
    VarSet v12( v1, v2 );
    VarSet v012 = v01 | v2;

    std::vector<Factor> facs;
    facs.push_back( Factor( v01 ).randomize() );
    facs.push_back( Factor( v12 ).randomize() );
    facs.push_back( Factor( v1 ).randomize() );
    std::vector<Var> vars;
    vars.push_back( v0 );
    vars.push_back( v1 );
    vars.push_back( v2 );

    FactorGraph G( facs );

    FactorGraph Gsmall = G.maximalFactors();
    BOOST_CHECK_EQUAL( Gsmall.nrVars(), 3 );
    BOOST_CHECK_EQUAL( Gsmall.nrFactors(), 2 );
    BOOST_CHECK_EQUAL( Gsmall.factor( 0 ), G.factor( 0 ) * G.factor( 2 ) );
    BOOST_CHECK_EQUAL( Gsmall.factor( 1 ), G.factor( 1 ) );

    size_t i = 0;
    for( size_t x = 0; x < 2; x++ ) {
        FactorGraph Gcl = G.clamped( i, x );
        BOOST_CHECK_EQUAL( Gcl.nrVars(), 3 );
        BOOST_CHECK_EQUAL( Gcl.nrFactors(), 3 );
        BOOST_CHECK_EQUAL( Gcl.factor(0), createFactorDelta(vars[i], x) );
        BOOST_CHECK_EQUAL( Gcl.factor(1), G.factor(0).slice(vars[i], x) * G.factor(2) );
        BOOST_CHECK_EQUAL( Gcl.factor(2), G.factor(1) );
    }
    i = 1;
    for( size_t x = 0; x < 2; x++ ) {
        FactorGraph Gcl = G.clamped( i, x );
        BOOST_CHECK_EQUAL( Gcl.nrVars(), 3 );
        BOOST_CHECK_EQUAL( Gcl.nrFactors(), 3 );
        BOOST_CHECK_EQUAL( Gcl.factor(0), createFactorDelta(vars[i], x) * G.factor(2).slice(vars[i],x) );
        BOOST_CHECK_EQUAL( Gcl.factor(1), G.factor(0).slice(vars[i], x) );
        BOOST_CHECK_EQUAL( Gcl.factor(2), G.factor(1).slice(vars[i], x) );
    }
    i = 2;
    for( size_t x = 0; x < 2; x++ ) {
        FactorGraph Gcl = G.clamped( i, x );
        BOOST_CHECK_EQUAL( Gcl.nrVars(), 3 );
        BOOST_CHECK_EQUAL( Gcl.nrFactors(), 3 );
        BOOST_CHECK_EQUAL( Gcl.factor(0), createFactorDelta(vars[i], x) );
        BOOST_CHECK_EQUAL( Gcl.factor(1), G.factor(0) );
        BOOST_CHECK_EQUAL( Gcl.factor(2), G.factor(1).slice(vars[i], x) * G.factor(2) );
    }
}


BOOST_AUTO_TEST_CASE( OperationsTest ) {
    Var v0( 0, 2 );
    Var v1( 1, 2 );
    Var v2( 2, 2 );
    VarSet v01( v0, v1 );
    VarSet v02( v0, v2 );
    VarSet v12( v1, v2 );
    VarSet v012 = v01 | v2;

    std::vector<Factor> facs;
    facs.push_back( Factor( v01 ).randomize() );
    facs.push_back( Factor( v12 ).randomize() );
    facs.push_back( Factor( v1 ).randomize() );
    std::vector<Var> vars;
    vars.push_back( v0 );
    vars.push_back( v1 );
    vars.push_back( v2 );

    FactorGraph G( facs );

    // clamp
    FactorGraph Gcl = G;
    for( size_t i = 0; i < 3; i++ )
        for( size_t x = 0; x < 2; x++ ) {
            Gcl.clamp( i, x, true );
            Factor delta = createFactorDelta( vars[i], x );
            BOOST_CHECK_EQUAL( Gcl.nrVars(), 3 );
            BOOST_CHECK_EQUAL( Gcl.nrFactors(), 3 );
            for( size_t j = 0; j < 3; j++ )
                if( G.factor(j).vars().contains( vars[i] ) )
                    BOOST_CHECK_EQUAL( Gcl.factor(j), G.factor(j) * delta );
                else
                    BOOST_CHECK_EQUAL( Gcl.factor(j), G.factor(j) );

            Gcl.restoreFactors();
            for( size_t j = 0; j < 3; j++ )
                BOOST_CHECK_EQUAL( Gcl.factor(j), G.factor(j) );
        }

    // clampVar
    for( size_t i = 0; i < 3; i++ )
        for( size_t x = 0; x < 2; x++ ) {
            Gcl.clampVar( i, std::vector<size_t>(1, x), true );
            Factor delta = createFactorDelta( vars[i], x );
            BOOST_CHECK_EQUAL( Gcl.nrVars(), 3 );
            BOOST_CHECK_EQUAL( Gcl.nrFactors(), 3 );
            for( size_t j = 0; j < 3; j++ )
                if( G.factor(j).vars().contains( vars[i] ) )
                    BOOST_CHECK_EQUAL( Gcl.factor(j), G.factor(j) * delta );
                else
                    BOOST_CHECK_EQUAL( Gcl.factor(j), G.factor(j) );

            Gcl.restoreFactors();
            for( size_t j = 0; j < 3; j++ )
                BOOST_CHECK_EQUAL( Gcl.factor(j), G.factor(j) );
        }
    for( size_t i = 0; i < 3; i++ )
        for( size_t x = 0; x < 2; x++ ) {
            Gcl.clampVar( i, std::vector<size_t>(), true );
            BOOST_CHECK_EQUAL( Gcl.nrVars(), 3 );
            BOOST_CHECK_EQUAL( Gcl.nrFactors(), 3 );
            for( size_t j = 0; j < 3; j++ )
                if( G.factor(j).vars().contains( vars[i] ) )
                    BOOST_CHECK_EQUAL( Gcl.factor(j), G.factor(j) * 0.0 );
                else
                    BOOST_CHECK_EQUAL( Gcl.factor(j), G.factor(j) );

            Gcl.restoreFactors();
            for( size_t j = 0; j < 3; j++ )
                BOOST_CHECK_EQUAL( Gcl.factor(j), G.factor(j) );
        }
    std::vector<size_t> both;
    both.push_back( 0 );
    both.push_back( 1 );
    for( size_t i = 0; i < 3; i++ )
        for( size_t x = 0; x < 2; x++ ) {
            Gcl.clampVar( i, both, true );
            BOOST_CHECK_EQUAL( Gcl.nrVars(), 3 );
            BOOST_CHECK_EQUAL( Gcl.nrFactors(), 3 );
            for( size_t j = 0; j < 3; j++ )
                BOOST_CHECK_EQUAL( Gcl.factor(j), G.factor(j) );
            Gcl.restoreFactors();
        }

    // clampFactor
    for( size_t x = 0; x < 4; x++ ) {
        Gcl.clampFactor( 0, std::vector<size_t>(1,x), true );
        Factor mask( v01, 0.0 );
        mask.set( x, 1.0 );
        BOOST_CHECK_EQUAL( Gcl.factor(0), G.factor(0) * mask );
        BOOST_CHECK_EQUAL( Gcl.factor(1), G.factor(1) );
        BOOST_CHECK_EQUAL( Gcl.factor(2), G.factor(2) );
        Gcl.restoreFactor( 0 );
    }
    for( size_t x = 0; x < 4; x++ ) {
        Gcl.clampFactor( 1, std::vector<size_t>(1,x), true );
        Factor mask( v12, 0.0 );
        mask.set( x, 1.0 );
        BOOST_CHECK_EQUAL( Gcl.factor(0), G.factor(0) );
        BOOST_CHECK_EQUAL( Gcl.factor(1), G.factor(1) * mask );
        BOOST_CHECK_EQUAL( Gcl.factor(2), G.factor(2) );
        Gcl.restoreFactor( 1 );
    }
    for( size_t x = 0; x < 2; x++ ) {
        Gcl.clampFactor( 2, std::vector<size_t>(1,x), true );
        Factor mask( v1, 0.0 );
        mask.set( x, 1.0 );
        BOOST_CHECK_EQUAL( Gcl.factor(0), G.factor(0) );
        BOOST_CHECK_EQUAL( Gcl.factor(1), G.factor(1) );
        BOOST_CHECK_EQUAL( Gcl.factor(2), G.factor(2) * mask );
        Gcl.restoreFactors();
    }

    // makeCavity
    FactorGraph Gcav( G );
    Gcav.makeCavity( 0, true );
    BOOST_CHECK_EQUAL( Gcav.factor(0), Factor( v01, 1.0 ) );
    BOOST_CHECK_EQUAL( Gcav.factor(1), G.factor(1) );
    BOOST_CHECK_EQUAL( Gcav.factor(2), G.factor(2) );
    Gcav.restoreFactors();
    Gcav.makeCavity( 1, true );
    BOOST_CHECK_EQUAL( Gcav.factor(0), Factor( v01, 1.0 ) );
    BOOST_CHECK_EQUAL( Gcav.factor(1), Factor( v12, 1.0 ) );
    BOOST_CHECK_EQUAL( Gcav.factor(2), Factor( v1, 1.0 ) );
    Gcav.restoreFactors();
    Gcav.makeCavity( 2, true );
    BOOST_CHECK_EQUAL( Gcav.factor(0), G.factor(0) );
    BOOST_CHECK_EQUAL( Gcav.factor(1), Factor( v12, 1.0 ) );
    BOOST_CHECK_EQUAL( Gcav.factor(2), G.factor(2) );
    Gcav.restoreFactors();
}


BOOST_AUTO_TEST_CASE( IOTest ) {
    Var v0( 0, 2 );
    Var v1( 1, 2 );
    Var v2( 2, 2 );
    VarSet v01( v0, v1 );
    VarSet v02( v0, v2 );
    VarSet v12( v1, v2 );
    VarSet v012 = v01 | v2;

    std::vector<Factor> facs;
    facs.push_back( Factor( v01 ).randomize() );
    facs.push_back( Factor( v12 ).randomize() );
    facs.push_back( Factor( v1 ).randomize() );
    std::vector<Var> vars;
    vars.push_back( v0 );
    vars.push_back( v1 );
    vars.push_back( v2 );

    FactorGraph G( facs );

    G.WriteToFile( "factorgraph_test.fg" );
    FactorGraph G2;
    G2.ReadFromFile( "factorgraph_test.fg" );

    BOOST_CHECK( G.vars() == G2.vars() );
    BOOST_CHECK( G.bipGraph() == G2.bipGraph() );
    BOOST_CHECK_EQUAL( G.nrFactors(), G2.nrFactors() );
    for( size_t I = 0; I < G.nrFactors(); I++ ) {
        BOOST_CHECK( G.factor(I).vars() == G2.factor(I).vars() );
        for( size_t s = 0; s < G.factor(I).nrStates(); s++ )
            BOOST_CHECK_CLOSE( G.factor(I)[s], G2.factor(I)[s], tol );
    }

    std::stringstream ss;
    std::string s;
    G.printDot( ss );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "graph FactorGraph {" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "node[shape=circle,width=0.4,fixedsize=true];" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "\tv0;" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "\tv1;" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "\tv2;" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "node[shape=box,width=0.3,height=0.3,fixedsize=true];" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "\tf0;" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "\tf1;" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "\tf2;" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "\tv0 -- f0;" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "\tv1 -- f0;" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "\tv1 -- f1;" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "\tv1 -- f2;" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "\tv2 -- f1;" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "}" );

    G.setFactor( 0, Factor( G.factor(0).vars(), 1.0 ) );
    G.setFactor( 1, Factor( G.factor(1).vars(), 2.0 ) );
    G.setFactor( 2, Factor( G.factor(2).vars(), 3.0 ) );
    ss << G;
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "3" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "2" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "0 1 " );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "2 2 " );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "4" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "0          1" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "1          1" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "2          1" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "3          1" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "2" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "1 2 " );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "2 2 " );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "4" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "0          2" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "1          2" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "2          2" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "3          2" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "1" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "1 " );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "2 " );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "2" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "0          3" );
    std::getline( ss, s ); BOOST_CHECK_EQUAL( s, "1          3" );

    BOOST_CHECK_EQUAL( G.toString(), "3\n\n2\n0 1 \n2 2 \n4\n0          1\n1          1\n2          1\n3          1\n\n2\n1 2 \n2 2 \n4\n0          2\n1          2\n2          2\n3          2\n\n1\n1 \n2 \n2\n0          3\n1          3\n" );

    ss << G;
    FactorGraph G3;
    ss >> G3;
    BOOST_CHECK( G.vars() == G3.vars() );
    BOOST_CHECK( G.bipGraph() == G3.bipGraph() );
    BOOST_CHECK_EQUAL( G.nrFactors(), G3.nrFactors() );
    for( size_t I = 0; I < G.nrFactors(); I++ ) {
        BOOST_CHECK( G.factor(I).vars() == G3.factor(I).vars() );
        for( size_t s = 0; s < G.factor(I).nrStates(); s++ )
            BOOST_CHECK_CLOSE( G.factor(I)[s], G3.factor(I)[s], tol );
    }

    FactorGraph G4;
    G4.fromString( G.toString() );
    BOOST_CHECK( G.vars() == G4.vars() );
    BOOST_CHECK( G.bipGraph() == G4.bipGraph() );
    BOOST_CHECK_EQUAL( G.nrFactors(), G4.nrFactors() );
    for( size_t I = 0; I < G.nrFactors(); I++ ) {
        BOOST_CHECK( G.factor(I).vars() == G4.factor(I).vars() );
        for( size_t s = 0; s < G.factor(I).nrStates(); s++ )
            BOOST_CHECK_CLOSE( G.factor(I)[s], G4.factor(I)[s], tol );
    }
}
