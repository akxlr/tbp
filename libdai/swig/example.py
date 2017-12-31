# This file is part of libDAI - http:#www.libdai.org/
#
# Copyright (c) 2006-2011, The libDAI authors. All rights reserved.
#
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.


# This example program illustrates how to read a factrograph from
# a file and run Belief Propagation, Max-Product and JunctionTree on it.
# This version uses the SWIG python wrapper of libDAI


import dai
import sys

a = dai.IntVector()

if len(sys.argv) != 2 and len(sys.argv) != 3:
    print 'Usage:', sys.argv[0], "<filename.fg> [maxstates]"
    print 'Reads factor graph <filename.fg> and runs'
    print 'Belief Propagation, Max-Product and JunctionTree on it.'
    print 'JunctionTree is only run if a junction tree is found with'
    print 'total number of states less than <maxstates> (where 0 means unlimited).'
    sys.exit(1)
else:
    # Report inference algorithms built into libDAI
#   print 'Builtin inference algorithms:', dai.builtinInfAlgNames()
#   TODO THIS CRASHES

    # Read FactorGraph from the file specified by the first command line argument
    fg = dai.FactorGraph()
    fg.ReadFromFile(sys.argv[1])
    maxstates = 1000000
    if len(sys.argv) == 3:
        maxstates = int(sys.argv[2])

    # Set some constants
    maxiter = 10000
    tol = 1e-9
    verb = 1

    # Store the constants in a PropertySet object
    opts = dai.PropertySet()
    opts["maxiter"] = str(maxiter)   # Maximum number of iterations
    opts["tol"] = str(tol)           # Tolerance for convergence
    opts["verbose"] = str(verb)      # Verbosity (amount of output generated)

    # Bound treewidth for junctiontree
    do_jt = True
    try:
        dai.boundTreewidth( fg, dai.eliminationCost_MinFill, maxstates )
        # TODO fix memory leak: no destructor for BigInt
    except: # TODO cannot catch libDAI exceptions yet
#       if( e.getCode() == Exception::OUT_OF_MEMORY ) TODO
        do_jt = False
        print "Skipping junction tree (need more than", maxstates, "states)."

    if do_jt:
        # Construct a JTree (junction tree) object from the FactorGraph fg
        # using the parameters specified by opts and an additional property
        # that specifies the type of updates the JTree algorithm should perform
        jtopts = opts
        jtopts["updates"] = "HUGIN"
        jt = dai.JTree( fg, jtopts )
        # Initialize junction tree algorithm
        jt.init()
        # Run junction tree algorithm
        jt.run()

        # Construct another JTree (junction tree) object that is used to calculate
        # the joint configuration of variables that has maximum probability (MAP state)
        jtmapopts = opts
        jtmapopts["updates"] = "HUGIN"
        jtmapopts["inference"] = "MAXPROD"
        jtmap = dai.JTree( fg, jtmapopts )
        # Initialize junction tree algorithm
        jtmap.init()
        # Run junction tree algorithm
        jtmap.run()
        # Calculate joint state of all variables that has maximum probability
        jtmapstate = jtmap.findMaximum()

    # Construct a BP (belief propagation) object from the FactorGraph fg
    # using the parameters specified by opts and two additional properties,
    # specifying the type of updates the BP algorithm should perform and
    # whether they should be done in the real or in the logdomain
    bpopts = opts
    bpopts["updates"] = "SEQRND"
    bpopts["logdomain"] = "0"
    bp = dai.BP( fg, bpopts )
    # Initialize belief propagation algorithm
    bp.init()
    # Run belief propagation algorithm
    bp.run()

    # Construct a BP (belief propagation) object from the FactorGraph fg
    # using the parameters specified by opts and two additional properties,
    # specifying the type of updates the BP algorithm should perform and
    # whether they should be done in the real or in the logdomain
    #
    # Note that inference is set to MAXPROD, which means that the object
    # will perform the max-product algorithm instead of the sum-product algorithm
    mpopts = opts
    mpopts["updates"] = "SEQRND"
    mpopts["logdomain"] = "0"
    mpopts["inference"] = "MAXPROD"
    mpopts["damping"] = "0.1"
    mp = dai.BP( fg, mpopts )
    # Initialize max-product algorithm
    mp.init()
    # Run max-product algorithm
    mp.run()
    # Calculate joint state of all variables that has maximum probability
    # based on the max-product result
    mpstate = mp.findMaximum()

    # Construct a decimation algorithm object from the FactorGraph fg
    # using the parameters specified by opts and three additional properties,
    # specifying that the decimation algorithm should use the max-product
    # algorithm and should completely reinitalize its state at every step
    decmapopts = opts
    decmapopts["reinit"] = "1"
    decmapopts["ianame"] = "BP"
    decmapopts["iaopts"] = "[damping=0.1,inference=MAXPROD,logdomain=0,maxiter=1000,tol=1e-9,updates=SEQRND,verbose=1]"
    decmap = dai.DecMAP( fg, decmapopts )
    decmap.init()
    decmap.run()
    decmapstate = decmap.findMaximum()

    if do_jt:
        # Report variable marginals for fg, calculated by the junction tree algorithm
        print 'Exact variable marginals:'
        for i in range(fg.nrVars()):                 # iterate over all variables in fg
            print jt.belief(dai.VarSet(fg.var(i)))   # display the "belief" of jt for that variable

    # Report variable marginals for fg, calculated by the belief propagation algorithm
    print 'Approximate (loopy belief propagation) variable marginals:'
    for i in range(fg.nrVars()):                     # iterate over all variables in fg
        print bp.belief(dai.VarSet(fg.var(i)))       # display the belief of bp for that variable

    if do_jt:
        # Report factor marginals for fg, calculated by the junction tree algorithm
        print 'Exact factor marginals:'
        for I in range(fg.nrFactors()):              # iterate over all factors in fg
            print jt.belief(fg.factor(I).vars())     # display the "belief" of jt for the variables in that factor

    # Report factor marginals for fg, calculated by the belief propagation algorithm
    print 'Approximate (loopy belief propagation) factor marginals:'
    for I in range(fg.nrFactors()):                  # iterate over all factors in fg
        print bp.belief(fg.factor(I).vars())         # display the belief of bp for the variables in that factor

    if do_jt:
        # Report log partition sum (normalizing constant) of fg, calculated by the junction tree algorithm
        print 'Exact log partition sum:', jt.logZ()

    # Report log partition sum of fg, approximated by the belief propagation algorithm
    print 'Approximate (loopy belief propagation) log partition sum:', bp.logZ()

    if do_jt:
        # Report exact MAP variable marginals
        print 'Exact MAP variable marginals:'
        for i in range(fg.nrVars()):
            print jtmap.belief(dai.VarSet(fg.var(i)))

    # Report max-product variable marginals
    print 'Approximate (max-product) MAP variable marginals:'
    for i in range(fg.nrVars()):
        print mp.belief(dai.VarSet(fg.var(i)))

    if do_jt:
        # Report exact MAP factor marginals
        print 'Exact MAP factor marginals:'
        for I in range(fg.nrFactors()):
            print jtmap.belief(fg.factor(I).vars()), '==', jtmap.belief(fg.factor(I).vars())

    # Report max-product factor marginals
    print 'Approximate (max-product) MAP factor marginals:'
    for I in range(fg.nrFactors()):
        print mp.belief(fg.factor(I).vars()), '==', mp.belief(fg.factor(I).vars())

    if do_jt:
        # Report exact MAP joint state
        print 'Exact MAP state (log score =', fg.logScore( jtmapstate ), '):'
        for i in range(len(jtmapstate)):
            print fg.var(i), ':', jtmapstate[i]

    # Report max-product MAP joint state
    print 'Approximate (max-product) MAP state (log score =', fg.logScore( mpstate ), '):'
    for i in range(len(mpstate)):
        print fg.var(i), ':', mpstate[i]

    # Report DecMAP joint state
    print 'Approximate DecMAP state (log score =', fg.logScore( decmapstate ), '):'
    for i in range(len(decmapstate)):
        print fg.var(i), ':', decmapstate[i]
