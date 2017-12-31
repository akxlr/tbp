"""
Tests to run after installing TBP.

Author: Andrew Wrigley, National University of Singapore and Australian National University
"""

import tbp

ERROR_THRESHOLD = 1E-12

def run_tests():

    true_mar = tbp.load_mar('tests/sample_true.MAR')
    g = tbp.load_fg_graph('tests/sample.fg')
    g.apply_evidence('tests/sample.evid')
    dg = g.decompose(r=2)
    mar = dg.tbp_marg(k=10000)
    assert tbp.l1_error(mar, true_mar) <= ERROR_THRESHOLD, tbp.l1_error(mar, true_mar)

if __name__ == '__main__':
    run_tests()

