"""
Code to run the experiments in Wrigley, Lee, and Ye (2017). See https://github.com/akxlr/tbp.

Author: Andrew Wrigley, National University of Singapore and Australian National University
"""

import glob
import json
import os
from typing import List, Tuple
import tbp
import time
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt

BASE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../')


def run_tests(tests, marg_params=None, binary_err=False):
    """
    :param tests: List of tests, each of which is a dict with keys {name, g, dg, ks, true_marg, plot_series}.
    :return: None; adds a new key `results` to each test, whose value is a dict {k -> {infer_time, pred_marg, err}}.
    """

    for i, test in enumerate(tests):

        print("Running test {} of {}: {}".format(i+1, len(tests), test['name']))

        results = {}
        for k in test['ks']:
            results[k] = {}
            t0 = time.time()
            results[k]['pred_marg'] = test['dg'].tbp_marg(k=k,
                **(marg_params if marg_params else {}))
            results[k]['infer_time'] = time.time() - t0
            results[k]['marg_err'] = tbp.l1_error(results[k]['pred_marg'], test['true_marg'], binary=binary_err)
            print("  (k={}) Mean L1 error: {}".format(k, results[k]['marg_err']))

        test['results'] = results


def save_plot(results, name):
    """
    :param results: Dict prepared for plotting in the form {series_label -> {k -> err}}
    :return:
    """

    if not os.path.exists('plots'):
        os.makedirs('plots')

    # # Add padding so we can see the markers
    # xticks, xticklabels = plt.xticks()
    # xmin = (3 * xticks[0] - xticks[1]) / 2.
    # xmax = (3 * xticks[-1] - xticks[-2]) / 2.
    # plt.xlim(xmin, xmax)
    # plt.xticks(xticks)

    # Draw plot
    plt.title(name)
    plt.xlabel("Sample size K")
    plt.xscale('log')
    plt.ylabel("Average classification error")
    # markers = ['-.co', '-ks', ':g^', '-ks', '-.ks', ':g', '-m', '-.y', '--r^']

    for series_label in results.keys():
        ks = sorted(results[series_label].keys())
        plt.plot(ks, [results[series_label][k] for k in ks], label=series_label)

    plt.legend(loc='best', ncol=3)
    file_id = time.time()
    plot_filename = "plots/%s-%s.png" % (name, file_id)
    plt.savefig(plot_filename)
    plt.gcf().clear()
    print("%s saved\n" % plot_filename)
    with open("plots/sourcedata_%s.json" % file_id, 'w') as f:
        print(json.dumps(results), file=f)


def plot_tests(tests, name):
    """
    Plot classification error versus sample size k.
    :param tests: As per run_tests, after invoking run_tests to fill results.
    """

    # Collect into groups according to plot_series key
    series = defaultdict(list)
    for test in tests:
        if 'plot_series' not in test:
            test['plot_series'] = '(no label)'
        series[test['plot_series']].append(test)

    results = {}
    for series_label in series.keys():
        # Average together the error of all tests with the same plot_series
        results[series_label] = defaultdict(list)
        for test in series[series_label]:
            for k in test['ks']:
                results[series_label][k].append(test['results'][k]['marg_err'])
        for k in results[series_label].keys():
            results[series_label][k] = np.mean(results[series_label][k])

    save_plot(results, name)


def plot_from_files(series, out_name):
    """
    Combine a set of sourcedata_XXX.json files into a single plot.
    :param series: List of series, each a tuple of the form (filename, series_name, new_series_name).
    """
    results = {}
    for (filename, series_name, new_series_name) in series:
        if not filename.startswith('plots/'):
            filename = 'plots/' + filename
        with open(filename, 'r') as f:
            results[new_series_name] = {int(k): err for k, err in json.loads(f.read())[series_name].items()}
    save_plot(results, out_name)


def load_uai_tests(globs, r, ks, force=False) -> List[Tuple]:
    """
    Load all UAI tests matching any glob pattern in globs.
    :param r: Number of rank-1 terms to use for tensor decomposition
    :param force: Force all input graphs to be decomposed again, otherwise files with an existing decomposition are
    skipped.
    :return: list of tests, each of which is a dict as per run_tests.
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
        tests.append({
            'name': test_name,
            'g': g,
            'dg': dg,
            'ks': ks,
            'true_marg': true_marg,
        })

    return tests


def load_ising_tests(ks, sample_size=100) -> List[Tuple]:
    """
    Load the Ising model tests used to produce Wrigley, Lee and Ye (2017), Figure 1.
    """

    unary_str = 1
    pairwise_str = 2
    N = 10
    tests = []

    for i in range(sample_size):
        test_name = 'ising_{}_mixed_{}'.format(N, i)
        print("Generating test case {}".format(test_name))
        g, dg = tbp.ising_g(N, -unary_str, unary_str, -pairwise_str, pairwise_str)
        true_marg = g.exact_marg_elim()
        tests.append({
            'name': test_name,
            'g': g,
            'dg': dg,
            'ks': ks,
            'plot_series': 'mixed',
            'true_marg': true_marg,
        })

    for i in range(sample_size):
        test_name = 'ising_{}_attractive_{}'.format(N, i)
        print("Generating test case {}".format(test_name))
        g, dg = tbp.ising_g(N, -unary_str, unary_str, 0, pairwise_str)
        true_marg = g.exact_marg_elim()
        tests.append({
            'name': test_name,
            'g': g,
            'dg': dg,
            'ks': ks,
            'plot_series': 'attractive',
            'true_marg': true_marg,
        })

    return tests


def load_random_tests(ks, sample_size=100) -> List[Tuple]:
    """
    Load the random MRF tests used to produce Wrigley, Lee and Ye (2017), Figure 2.
    """

    unary_str = 1
    pairwise_str = 2
    edge_prob = 0.5
    N = 15
    tests = []

    # First two tests - varying k
    for i in range(sample_size):
        test_name = 'random_k_{}_mixed_{}'.format(N, i)
        print("Generating test case {}".format(test_name))
        g, dg = tbp.random_g(N, edge_prob, -unary_str, unary_str, -pairwise_str, pairwise_str)
        true_marg = g.exact_marg_elim()
        tests.append({
            'name': test_name,
            'g': g,
            'dg': dg,
            'ks': ks,
            'plot_series': 'mixed',
            'true_marg': true_marg,
        })

    for i in range(sample_size):
        test_name = 'random_k_{}_attractive_{}'.format(N, i)
        print("Generating test case {}".format(test_name))
        g, dg = tbp.random_g(N, edge_prob, -unary_str, unary_str, 0, pairwise_str)
        true_marg = g.exact_marg_elim()
        tests.append({
            'name': test_name,
            'g': g,
            'dg': dg,
            'ks': ks,
            'plot_series': 'attractive',
            'true_marg': true_marg,
        })

    # TODO Second two tests - varying N


    return tests

def run_ising():
    tests_ising = load_ising_tests(ks=[10, 100, 1000, 10000, 100000], sample_size=100)
    # tests_ising = load_ising_tests(ks=[10, 100, 1000], sample_size=2)
    run_tests(tests_ising, binary_err=True)
    plot_tests(tests_ising, 'icml17-ising')

def run_random():
    tests_random = load_random_tests(ks=[10, 100, 1000, 10000, 100000], sample_size=100)
    # tests_random = load_random_tests(ks=[10, 100, 1000], ns=[10, 15], sample_size=2)
    run_tests(tests_random, binary_err=True)
    plot_tests(tests_random, 'icml17-random')

def run_uai():
    tests_uai_2 = load_uai_tests(['linkage_*.uai', 'Promedus_*.uai'], r=2, ks=[10, 100, 1000, 10000])
    tests_uai_4 = load_uai_tests(['linkage_*.uai', 'Promedus_*.uai'], r=4, ks=[10, 100, 1000, 10000])
    run_tests(tests_uai_2)
    run_tests(tests_uai_4)

def run_all():
    run_ising()
    run_uai()

if __name__ == '__main__':
    run_ising()
    # run_random()
    # plot_from_files([
    #     ('sourcedata_1529829179.731084.json', 'mixed', 'ising-marg-minfill'),
    # ], 'mixed')


