"""
Tests to run after installing TBP.

Author: Andrew Wrigley, National University of Singapore and Australian National University
"""

import os
import tbp

BASE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../')

def test_small_ising():
    """
    Check that for a small test Ising model, approximate marginals are close to exact marginals for large sample size.
    """
    g = tbp.load_fg_graph(os.path.join(BASE_DIR, 'tests/ising_8x8.fg'))
    g.apply_evidence(os.path.join(BASE_DIR, 'tests/ising_8x8.evid'))
    print("Computing true marginals with elimination algorithm")
    true_mar = g.exact_marg_elim()
    print("Decomposing graph")
    dg = g.decompose(r=2)
    print("Running TBP")
    mar = dg.tbp_marg(k=100000)
    # Error should be around 0.002-0.003
    assert tbp.l1_error(mar, true_mar) <= 0.01

