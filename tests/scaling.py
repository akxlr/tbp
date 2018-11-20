from typing import List, Tuple
import tbp
import json
import copy
import time
import pickle
import os
import glob
import pandas as pd
import numpy as np


def transpose_elim_order(M, N, vars_list):
    origin = list(range(len(vars_list)))
    return np.array(origin).reshape((M,N)).T.reshape(-1).tolist()


def sample_ising_models(M, N, interaction, sample_size=100) -> List[Tuple]:
    unary_str = 1
    pairwise_str = 2
    models = []

    try:
        for i in range(sample_size):
            model_name = 'ising_{}x{}_{}_{}'.format(M, N, interaction, i)
            print("Generating test case {}".format(model_name))
            if interaction == 'attractive':
                params = (M, N) + (-unary_str, unary_str, 0, pairwise_str)
            elif interaction == 'mixed':
                params = (M, N) + (-unary_str, unary_str, -pairwise_str, pairwise_str)
            g, dg = tbp.asymmetric_ising_g(*params)
            elim_order = None if M > N else transpose_elim_order(M, N, g.get_vars_list()) # better elimination ordering
            true_marg = g.exact_marg_elim()
            models.append({
                'name': model_name,
                'sample_id': i,
                'M': M,
                'N': N,
                'interaction': interaction, 
                'g': g,
                'dg': dg,
                'true_marg': true_marg,
            })
    except MemoryError:
        print("MEMORY ERROR")
        models = []

    return models


def save_models(models, name):
    os.makedirs ('models', exist_ok=True)
    with open("models/{}.p".format(name), 'wb') as f:
        pickle.dump(models, f)


def run_tests(tests, k, marg_params=None, binary_err=False):
    for i, test in enumerate(tests):
        try:
            print("Running test {} of {}: {}".format(i+1, len(tests), test['name']))
            result = {}
            t0 = time.time()
            result['pred_marg'] = test['dg'].tbp_marg(k=k, **(marg_params if marg_params else {}))
            result['infer_time'] = time.time() - t0
            result['marg_err'] = tbp.l1_error(result['pred_marg'], test['true_marg'], binary=binary_err)
            print("  (k={}) Mean L1 error: {}".format(k, result['marg_err']))
        except MemoryError:
            print("MEMORY ERROR")
            result['pred_marg'] = -1
            result['infer_time'] = -1
            result['marg_err'] = -1
        test['result'] = result


def create_tests(models, params):
    test_params = {
        'n_run': params['n_run'],
        'k': params['k'],
        'processes': params['processes']
    }
    tests = copy.copy(models)
    for test in tests:
        test.update(test_params)
    return tests


def save_results(tests, name):
    os.makedirs ('results', exist_ok=True)
    results = []
    for test in tests:
        results.append({
            'name': test['name'],
            'sample_id': test['sample_id'],
            'M': test['M'],
            'N': test['N'],
            'interaction': test['interaction'], 
            'k': test['k'],
            'n_run': test['n_run'],
            'processes': test['processes'],
            'infer_time': test['result']['infer_time'],
            'marg_err': test['result']['marg_err'],
            })
    with open("results/{}.json".format(name), 'w') as f:
        print(json.dumps(results), file=f)


def load_results_as_dataframe(name):
    results = []
    for filepath in glob.glob("results/{}.json".format(name)):
        with open(filepath, "rb") as f:
            results += json.load(f)
        
    return pd.DataFrame.from_records(results)


def group_dataframe(df):
    return df.drop(columns=['sample_id']).groupby(['interaction', 'M', 'N', 'k', 'n_run', 'processes']).agg([np.mean, np.std])


def save_to_xlsx(dfs, name):
    writer = pd.ExcelWriter('{}'.format(name))
    for sheet, df in dfs:
        df.to_excel(writer, sheet)
    writer.save()


def compile_results():
    df1 = group_dataframe(load_results_as_dataframe('scaling-ising-icml17-*'))
    df2 = group_dataframe(load_results_as_dataframe('scaling-ising-runs-*'))
    save_to_xlsx([('icml17', df1), ('runs', df2)], 'results/compiled.xlsx')


def run_experiments():
    # DIMENSIONS = [(10, 10), (20, 10), (40, 10), (80, 10), (160, 10), (320, 10)]
    DIMENSIONS = [(10, 5)]
    # DIMENSIONS = [(15, 10)]
    # KS = [10, 100, 1000, 10000, 100000]
    KS = [10, 100, 1000]
    # K_RUNS = [(1000, 100), (10000, 10), (100000, 1)]
    K_RUNS = [(100, 100), (1000, 10), (10000, 1)]

    dimensions = DIMENSIONS + [(N, M) for M, N in DIMENSIONS]
    interactions = ['attractive', 'mixed']
    # sample_size = 100
    sample_size = 1

    for M, N in dimensions:
        for interaction in interactions:
            t0 = time.time()
            models = sample_ising_models(M, N, interaction, sample_size=sample_size)
            print(time.time() - t0)

            save_models(models, 'ising_{}x{}_{}_{}'.format(M, N, interaction, sample_size))

            # ICML with different dimensions
            for k in KS:
                params = {'k': k, 'n_run': 1, 'processes': 1}
                tests = create_tests(models, params)
                run_tests(tests, params['k'], binary_err=True, marg_params={ 'n_run': params['n_run'], 'processes': params['processes'] })
                save_results(tests, "scaling-ising-icml17-{}x{}-{}-{}".format(M, N, interaction, params['k']))

            # Compare k vs runs
            k_runs = K_RUNS
            for k, n_run in k_runs:
                params = {'k': k, 'n_run': n_run, 'processes': 1}
                tests = create_tests(models, params)
                run_tests(tests, k, binary_err=True, marg_params={ 'n_run': params['n_run'], 'processes': params['processes'] })
                save_results(tests, "scaling-ising-runs-{}x{}-{}-k={}-run={}".format(M, N, interaction, params['k'], params['n_run']))


if __name__ == '__main__':
    run_experiments()
    compile_results()