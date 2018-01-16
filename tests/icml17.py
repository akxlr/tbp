"""
Code to run the experiments in Wrigley, Lee, and Ye (2017). See https://github.com/akxlr/tbp.

Author: Andrew Wrigley, National University of Singapore and Australian National University
"""

import glob
import os
from typing import List, Tuple
import tbp

BASE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../')

def load_fig1_tests() -> List[Tuple]:
    """
    Load the Ising model tests used to produce Wrigley, Lee and Ye (2017), Figure 1.
    """

    unary_str = 1
    pairwise_str = 2
    N = 10
    tests = []

    test_name = 'ising_{}_mixed'.format(N)
    print("Generating test case {}".format(test_name))
    g, dg = tbp.ising_g(N, -unary_str, unary_str, -pairwise_str, pairwise_str)
    true_marg = g.exact_marg_elim()
    tests.append((test_name, g, dg, true_marg))

    test_name = 'ising_{}_attractive'.format(N)
    print("Generating test case {}".format(test_name))
    g, dg = tbp.ising_g(N, -unary_str, unary_str, 0, pairwise_str)
    true_marg = g.exact_marg_elim()
    tests.append((test_name, g, dg, true_marg))

    return tests

def load_uai_tests(globs, r, force=False) -> List[Tuple]:
    """
    Load all UAI tests matching any glob patterns in globs.
    :param r: Number of rank-1 terms to use for tensor decomposition
    :param force: Force all input graphs to be decomposed again, otherwise files with an existing decomposition are
    skipped.
    :return: list of tests of the form (test_name, g, dg, true_marg) where g is the non-decomposed graph and dg the
    corresponding decomposed graph previously saved with decompose_all().
    """

    # Path that contains .uai files
    prob_folder = os.path.join(BASE_DIR, 'tests/uai/MAR_prob')
    # Path that contains .uai.MAR files
    sol_folder = os.path.join(BASE_DIR, 'tests/uai/MAR_sol')

    filenames = []
    for pattern in globs:
        filenames.extend(glob.glob(os.path.join(prob_folder, pattern)))

    tests = []
    for prob_filename in filenames:

        test_name = os.path.basename(prob_filename)
        print("Generating test case {}".format(test_name))

        # Load graph
        try:
            g = tbp.load_uai_graph(prob_filename)
        except tbp.BadGraphException:
            print("Graph {} could not be loaded, skipping".format(prob_filename))
            continue

        # Apply evidence
        evid_filename = '{}.evid'.format(prob_filename)
        g.apply_evidence(evid_filename)

        # Decompose graph
        decomp_filename = '{}.r={}.dfg'.format(prob_filename.replace('.uai', ''), r)
        if os.path.isfile(decomp_filename) and not force:
            print("Using existing decomposition {}".format(decomp_filename))
            dg = tbp.load_decomposed_graph(decomp_filename)
        else:
            dg = g.decompose(r)
            dg.save_to_file(decomp_filename)
            print("Saved {} ({} components)".format(decomp_filename, r))

        # Load true marginals
        sol_filename = os.path.join(sol_folder, '{}.MAR'.format(os.path.basename(prob_filename)))
        true_marg = tbp.load_mar(sol_filename)
        tests.append((test_name, g, dg, true_marg))

    return tests


def main():
    tests_ising = load_fig1_tests()
    tests_uai_2 = load_uai_tests(['linkage_*.uai', 'Promedus_*.uai'], r=2)
    tests_uai_4 = load_uai_tests(['linkage_*.uai', 'Promedus_*.uai'], r=4)

    tbp.run_tests(tests_ising, [10, 100, 1000, 10000, 100000], binary_err=True)
    tbp.run_tests(tests_uai_2, [10, 100, 1000, 10000])
    tbp.run_tests(tests_uai_4, [10, 100, 1000, 10000])


if __name__ == '__main__':
    main()

