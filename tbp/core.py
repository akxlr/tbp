"""
Python component of Tensor Belief Propagation. See https://github.com/akxlr/tbp.

Author: Andrew Wrigley, National University of Singapore and Australian National University
"""

import functools
import os
import sys
import subprocess
from operator import mul
from typing import List
import numpy as np
from numpy.linalg import LinAlgError

BASE_DIR = os.path.dirname(os.path.realpath(__file__))
TBP_BINARY = os.path.join(BASE_DIR, 'dfgmarg')
HOMEPAGE = "https://github.com/akxlr/tbp"

# How many rank-1 terms to decompose each initial potential function into
DEFAULT_COMPONENTS = 4
# Default sample size K for sampling step of TBP
DEFAULT_SAMPLE_SIZE = 10000
# Limits for truncating values when writing .dfg files (one order of magnitude inside C++ DBL_MIN, DBL_MAX)
DBL_MIN = 2.22507e-307
DBL_MAX = 1.79769e+307

verbosity = 1

# By default, tensor decomposition uses tensorly library. Alternatively, the MATLAB Tensor Toolbox by Kolda et al can
# be used (http://www.sandia.gov/tgkolda/TensorToolbox/). To use this:
#  - Install MATLAB normally
#  - Install MATLAB Python API:
#    https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html
#  - Install Tensor Toolbox from the above link
#  - Pass method='matlab' in the appropriate places below.
try:
    import matlab
    import matlab.engine
    eng = None
    MATLAB_PATH = os.path.join(BASE_DIR, 'matlab')
    with_matlab = True
except ImportError:
    with_matlab = False


def import_tensorly():
    """
    Importing TensorLy takes a few seconds and prints things (which is annoying when using `import core`
    for other things), so we only import when needed.
    """
    globals()['tensorly'] = __import__('tensorly')
    globals()['tensorly.decomposition'] = __import__('tensorly.decomposition')
    tensorly.set_backend('numpy')


def status(msg, level):
    if verbosity >= level:
        print(msg)


def safedouble(number):
    """
    Return number as a string, making sure it's within C++ bounds [DBL_MIN, DBL_MAX]. This is necessary to avoid
    underflow when reading .dfg files in tbp.cpp.
    """
    if number < DBL_MIN:
        number = DBL_MIN
    if number > DBL_MAX:
        number = DBL_MAX
    return str(number)


class SuppressOutput:
    """
    Context manager to suppress all output.
    """
    def __enter__(self):
        self._stdout = sys.stdout
        self._stderr = sys.stderr
        sys.stdout = self
        sys.stderr = self
        return self

    def write(self, string):
        pass

    def flush(self):
        pass

    def __exit__(self, *args):
        sys.stdout = self._stdout
        sys.stderr = self._stderr

class BadGraphException(Exception):
    pass

class FactorBase:

    vars = None

    @property
    def n_vars(self):
        return len(self.vars)

    @property
    def cardinalities(self):
        raise NotImplementedError()


class Factor(FactorBase):
    """
    Non-decomposed potential function represented as a multidimensional array.
    """
    def __init__(self, vars: List[int], table: np.ndarray):
        assert len(vars) == table.ndim, "Cannot create Factor: {} vars but table shape is {}".format(len(vars), table.shape)
        self.vars = vars
        self.table = table

    def __str__(self):
        return "Factor({},\n{})".format(self.vars, self.table)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (
            sorted(self.vars) == sorted(other.vars) and
            np.allclose(
                # Transpose axes so they match, in case each table is the same but variable order is different
                np.transpose(self.table, np.argsort(self.vars)), np.transpose(other.table, np.argsort(other.vars))
            )
        )

    def __mul__(self, other):
        """
        Multiply two Factors together using np.einsum.
        """
        # Resulting factor contains self.vars union other.vars sorted in numerical order
        combined_vars = sorted(set(self.vars).union(set(other.vars)))

        # Map variables to letters for np.einsum; works for up to 52 variables
        letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
        assert len(combined_vars) <= len(letters)
        letter_map = dict(zip(combined_vars, letters[:len(combined_vars)]))

        self_letters = ''.join([letter_map[var] for var in self.vars])
        other_letters = ''.join([letter_map[var] for var in other.vars])
        combined_letters = ''.join([letter_map[var] for var in combined_vars])

        return Factor(combined_vars, np.einsum('{},{}->{}'.format(self_letters, other_letters, combined_letters),
                                               self.table, other.table))

    @property
    def cardinalities(self):
        return self.table.shape

    def decompose(self, r, method='tensorly') -> 'DecomposedFactor':
        """
        Decompose this factor into r rank-1 components using non-negative CP decomposition. In certain trivial cases,
        returns a decomposition with one component rather than r.
        :param method: Library to use to compute CP decomposition, either 'tensorly' or 'matlab' (see note about MATLAB
        support at the top of this file).
        """

        assert method in ['tensorly', 'matlab']
        if method == 'matlab':
            assert with_matlab, 'MATLAB Python API not detected'

        # Variables with cardinality 1 cause problems with tensorly non_negative_decomp, so we only decompose the
        # nontrivial dimensions of the table (and set the matrices for trivial dimensions to [[1, ..., 1]]).
        table_nt = np.squeeze(self.table)

        if table_nt.ndim == 1:
            # If the nontrivial table has just one variable, the resulting decomposition has just one term (not r terms)
            # and we don't call the decomposition algorithm at all
            n_terms = 1
            weights = np.ones(n_terms)
            # Reshape is to change np.array([1,2,3]) to np.array([[1],[2],[3]])
            matrices_nt = [table_nt.reshape(table_nt.shape[0], 1)]
        else:
            n_terms = r

            if method == 'tensorly':

                if 'tensorly' not in sys.modules:
                    if verbosity == 3:
                        import_tensorly()
                    else:
                        with SuppressOutput():
                            import_tensorly()

                weights = np.ones(n_terms)
                matrices_nt = None
                table_nt_orig = table_nt
                i = 0
                while not matrices_nt:
                    try:
                        matrices_nt = tensorly.decomposition.non_negative_parafac(table_nt, n_terms)
                    except LinAlgError:
                        i += 1
                        status("Warning: singular matrix encountered, perturbing tensor and attempting decomposition "
                              "again ({})".format(i), 3)
                        table_nt = perturb(table_nt_orig)

            elif method == 'matlab':

                # Initialise MATLAB
                global eng
                if not eng:
                    eng = matlab.engine.start_matlab()
                    eng.addpath(MATLAB_PATH)

                w, T = eng.cp_kolda(matlab.double(table_nt.flatten(order='F').tolist()), matlab.int64(table_nt.shape),
                                    r, nargout=2)
                matrices_nt = [np.array(t) for t in T]
                weights = np.array(w).flatten()

        # Add the trivial dimensions back in
        matrices = []
        j = 0
        for i in range(self.n_vars):
            if self.cardinalities[i] == 1:
                matrices.append(np.ones((1, n_terms)))
            else:
                matrices.append(matrices_nt[j])
                j += 1
            assert matrices[i].shape == (self.cardinalities[i], n_terms)
        assert j == len(matrices_nt)

        df = DecomposedFactor(self.vars, weights, matrices)

        # Make sure the decomposition is close to the original factor
        # diff = df.expand().table - self.table
        # print("Reconstruction error: norm=%s, largest=%s" % (
        #     np.linalg.norm(diff),
        #     np.max(diff),
        # ))

        return df


class Graph:
    def __init__(self, factors: List[FactorBase]):
        self.factors = factors
        # Check variable cardinalities inferred by table shapes are consistent
        _ = self.get_cardinality_dict()

    def __str__(self):
        res = ''

        # Header line
        res += str(self.n_factors) + '\n'

        # Block for each factor
        for f in self.factors:
            res += '\n'
            res += str(f)

        return res

    @property
    def n_factors(self):
        return len(self.factors)

    def get_cardinality_list(self):
        return [x[1] for x in sorted(self.get_cardinality_dict().items(), key=lambda x: x[0])]

    def get_cardinality_dict(self):
        cardinalities = {}
        for f in self.factors:
            for var, cardinality in zip(f.vars, f.cardinalities):
                if var in cardinalities:
                    msg = "Variable cardinalities implied by table shapes are inconsistent"
                    assert cardinalities[var] == cardinality, msg
                else:
                    cardinalities[var] = cardinality
        return cardinalities

    def get_vars_list(self):
        v = set()
        for f in self.factors:
            v.update(f.vars)
        return sorted(v)

    def apply_evidence(self, evid_filename):
        """
        Apply evidence to graph by taking slices of factors containing clamped variables and adding new delta factors.
        Evidence format (.evid files): http://www.hlt.utdallas.edu/~vgogate/uai14-competition/evidformat.html
        :param evid_filename: .evid file
        """

        # Read evidence file to get evidence in the form of (var, val) pairs
        with open(evid_filename, 'r') as f:
            lines = [x.strip() for x in f.readlines() if x.strip()]
            # There is sometimes a '1' as the first line which we ignore
            assert len(lines) in [1, 2]
            if len(lines) == 1:
                line = lines[0]
            else:
                line = lines[1]
            fields = [int(x) for x in line.split()]
            evid_list = fields[1:]
            assert len(evid_list) == 2 * fields[0]
            evidence = []
            i = 0
            while i + 1 < len(evid_list):
                evidence.append((evid_list[i], evid_list[i+1]))
                i += 2

        # Replace factors containing clamped variables with slices
        new_factors = []
        for f in self.factors:
            new_f = f
            add_to_new = True
            for (var, val) in evidence:
                if var in new_f.vars:
                    axis = new_f.vars.index(var)
                    slices = [slice(None)] * len(new_f.vars)
                    slices[axis] = val
                    new_vars = new_f.vars[:axis] + new_f.vars[axis+1:]
                    if new_vars:
                        new_f = Factor(new_vars, new_f.table[slices])
                    else:
                        # We have clamped all variables in this factor, so it can be removed
                        add_to_new = False
                        break

            # Only include clamped factor if it's not empty
            if add_to_new:
                new_factors.append(new_f)

        # Add delta factors
        cardinalities = self.get_cardinality_dict()
        for (var, val) in evidence:
            table = np.zeros(cardinalities[var])
            table[val] = 1
            new_factors.append(Factor([var], table))

        self.factors = new_factors

    def elim_var(self, var):
        """
        Eliminate a single variable.
        """

        # Multiply all factors containing var
        prod = None
        for f in self.factors:
            if var in f.vars:
                if prod is None:
                    prod = f
                else:
                    prod *= f

        # Marginalise out var
        axis = prod.vars.index(var)
        new_f = Factor([v for v in prod.vars if v != var], np.sum(prod.table, axis=axis))

        self.factors = [f for f in self.factors if var not in f.vars] + [new_f]


    def exact_marg_elim(self):
        """
        Run a naive elimination algorithm for each variable in turn (inefficient).
        :return: Marginals as list of lists.
        """

        all_marg = []
        vars_list = self.get_vars_list()
        for var in vars_list:
            temp_g = Graph(list(self.factors))
            # Sequential elimination order, skipping var (for Ising models, this is row-wise starting at top left)
            elim_order = list(range(len(vars_list)))
            elim_order.remove(var)

            for var_to_elim in elim_order:
                temp_g.elim_var(var_to_elim)

            # All factors are now just functions of var, multiply them together and normalise
            marg_unnormalised = functools.reduce(mul, temp_g.factors).table
            all_marg.append((marg_unnormalised / np.sum(marg_unnormalised)).tolist())

        return all_marg

    def decompose(self, r=DEFAULT_COMPONENTS, method='tensorly') -> 'DecomposedGraph':
        """
        Decompose all factors in this graph into sums of r rank-1 tensors.
        """
        return DecomposedGraph([f.decompose(r, method) for f in self.factors])

    def get_connected_components(self):
        """
        Split into a list of graphs, each of which is fully connected.
        """

        h = list(self.factors)
        components = []
        while len(h) > 0:
            connected_h = [h[0]]
            any_added = True
            while any_added == True:
                any_added = False
                for fact in h:
                    if fact not in connected_h:
                        # Try to add fact to connected_h
                        added = False
                        for cfact in connected_h:
                            for var in fact.vars:
                                if var in cfact.vars:
                                    connected_h.append(fact)
                                    added = True
                                    any_added = True
                                    break
                            if added:
                                break

            components.append(Graph(connected_h))
            for fact in connected_h:
                h.remove(fact)

        return components


    def show_stats(self):
        """
        Print some stats for this graph.
        """
        cardinalities = self.get_cardinality_list()
        components = self.get_connected_components()
        print("Graph info:")
        print("  n_factors: {}".format(self.n_factors))
        print("  n_vars: {}".format(len(cardinalities)))
        print("  cardinalities: min {} max {}".format(np.min(cardinalities), np.max(cardinalities)))
        print("  connected_components: {}".format(len(components)))
        print("  component n_factors: {}".format(', '.join([str(x.n_factors) for x in components])))
        return components


class DecomposedGraph(Graph):
    """
    Graph with all factors represented in decomposed form.
    """
    def save_to_file(self, filename):
        """
        File format as follows:

        n_factors

        for each factor:
            n_terms
            <weights>
            n_variables
            <variable indices>
            <variable cardinalities>
            for each matrix:
                n_nonzero
                1 0.5
                3 0.1
                4 0.1
                ...
                n_nonzero
                1 0.5
                3 0.1
                4 0.3
                ...
        """
        with open(filename, 'w') as f:
            f.write(str(self))

    def tbp_marg(self, k=DEFAULT_SAMPLE_SIZE):
        """
        Call libdai/utils/dfgmarg binary via pipe to get approximate marginals.
        :param k: TBP sample size
        by sampling, and then marginalise once.
        :return: Marginals as list of lists.
        """
        p = subprocess.run([TBP_BINARY, str(k)], stdout=subprocess.PIPE,
                           input=str(self), encoding='ascii')
        assert p.returncode == 0

        # Marginals are returned in .MAR format, convert this to list of lists
        marginals = []
        tokens = iter(p.stdout.split())
        n_vars = int(next(tokens))
        for i in range(n_vars):
            marginals.append([])
            cardinality = int(next(tokens))
            for j in range(cardinality):
                marginals[-1].append(float(next(tokens)))

        return marginals


class DecomposedFactor(FactorBase):
    """
    Decomposed potential function.
    """

    def __init__(self, vars: List[int], weights: np.ndarray, matrices: List[np.ndarray]):
        """
        :param vars: List of variables
        :param weights: Numpy array containing weights
        :param matrices: List of matrices (2d np.array), where each matrix corresponds to a single variable. Rows in
        each matrix correspond to variable configurations, columns to term number.
        """
        assert len(matrices) == len(vars)
        assert all([m.ndim == 2 for m in matrices])
        assert all([m.shape[1] == len(weights) for m in matrices])

        self.vars = vars
        self.weights = weights
        self.matrices = matrices

    def __str__(self):
        """
        Return this DecomposedFactor as a string per .dfg file format.
        """
        res = '\n'.join([
            str(self.n_terms),
            ' '.join(safedouble(x) for x in self.weights),
            str(self.n_vars),
            ' '.join(str(x) for x in self.vars),
            ' '.join(str(x) for x in self.cardinalities),
        ]) + '\n'
        for m in self.matrices:
            # 'F' = column-major order
            table_flat = m.flatten(order='F')
            nonzero = [(i, safedouble(table_flat[i])) for i in range(len(table_flat)) if table_flat[i]]
            res += str(len(nonzero)) + '\n'
            for x in nonzero:
                res += '{} {}\n'.format(*x)
        return res

    @property
    def n_vars(self):
        return len(self.vars)

    @property
    def cardinalities(self):
        return [m.shape[0] for m in self.matrices]

    @property
    def n_terms(self):
        return len(self.weights)

    def expand(self) -> 'Factor':
        """
        Reconstruct decomposition into a dense Factor represented in tabular form.
        """
        dims = tuple([m.shape[0] for m in self.matrices])
        table = np.zeros(dims)
        r = self.matrices[0].shape[1]
        for i in range(r):
            vs = [m[:, i] for m in self.matrices]
            table += self.weights[i] * functools.reduce(np.multiply.outer, vs)
        return Factor(self.vars, table)


def perturb(arr: np.ndarray, max_change=1e-12) -> np.ndarray:
    """
    Add a very small random positive number to each entry in arr.
    """
    return arr + (np.random.rand(*arr.shape) * max_change)


# def read_tokens(f):
#     for line in f:
#         if line[0] == '#':
#             continue
#         for token in line.split():
#             yield token

def read_next_line(f):
    x = None
    while not x or x[0] == '#':
        x = f.readline().strip()
    return x


def load_fg_graph(filename) -> Graph:
    """
    Load .fg graph (libDAI's file format for factor graphs) into memory.
    Model format (.fg files): https://staff.fnwi.uva.nl/j.m.mooij/libDAI/doc/fileformats.html
    :param filename: .fg file
    :return: Graph
    """

    with open(filename, 'r') as f:
        n_factors = int(read_next_line(f))
        g = []
        for i in range(n_factors):
            n_vars = int(read_next_line(f))
            vars = [int(x) for x in read_next_line(f).split()]
            cardinalities = [int(x) for x in read_next_line(f).split()]
            assert len(vars) == n_vars
            assert len(cardinalities) == n_vars
            n_nonzero = int(read_next_line(f))
            table_flat = np.zeros(functools.reduce(mul, cardinalities))
            for entry in range(n_nonzero):
                var_index, val = read_next_line(f).split()
                var_index = int(var_index)
                val = float(val)
                table_flat[var_index] = val
            g.append(Factor(vars, table_flat.reshape(cardinalities, order='F')))

    return Graph(g)


def load_uai_graph(filename) -> Graph:
    """
    Load .uai graph into memory.
    Model format (.uai files): http://www.hlt.utdallas.edu/~vgogate/uai14-competition/modelformat.html
    :param filename: .uai file
    :return: Graph
    """

    # TODO rewrite this to treat newlines the same as other whitespace - see read_tokens

    with open(filename, 'r') as f:

        # Preamble
        lines = [x.strip() for x in f.readlines() if x.strip()]
        assert lines[0] == 'MARKOV'
        n_vars = int(lines[1])
        cardinalities = [int(x) for x in lines[2].split()]
        n_cliques = int(lines[3])
        clique_varsets = []
        for line in lines[4:4+n_cliques]:
            fields = [int(x) for x in line.split()]
            clique_vars = fields[1:]
            assert len(clique_vars) == fields[0]
            clique_varsets.append(clique_vars)
        assert len(clique_varsets) == n_cliques

        # Function tables
        # For each table, we need to read a single list (table_flat) which we then reshape. The size
        # of the list is given as the first entry (sometimes on its own line and sometimes not), and table_flat
        # itself is not always on a single line. We need to read lines until we have read all of the table entries
        # (line breaks are arbitrary and mean nothing - one line can contain any number of table entries).

        factors = []
        j = 0 # Current clique
        try:
            table_size = int(lines[4+n_cliques])
        except ValueError:
            # TODO table_size is on the same line as table_flat
            raise BadGraphException

        table_flat = []
        i = 4+n_cliques+1 # Current file line
        while True:

            table_flat += [float(x) for x in lines[i].split()]
            if len(table_flat) == table_size:

                # We have read the entire table - reshape and add to factors list
                shape = [cardinalities[var] for var in clique_varsets[j]]
                # Number of entries in table should equal product of variable cardinalities
                assert table_size == functools.reduce(mul, shape, 1)
                # Row-major order
                table = np.reshape(table_flat, shape)
                factors.append(Factor(clique_varsets[j], table))

                if i+1 < len(lines):
                    j += 1
                    table_flat = []
                    # The next line contains the size of the next table
                    table_size = int(lines[i+1])
                    i += 2
                else:
                    break
            else:
                i += 1

        g = Graph(factors)
        assert cardinalities == g.get_cardinality_list()
        return g


def load_mar(filename):
    """
    Load .MAR file containing true marginals.
    :param filename: .MAR file
    :return: Marginals as a list of lists [[<x1 marginals>], [<x2 marginals], ...]
    """

    # .MAR file contains a header line "MAR", then a single line of the form:
    #  <n_vars> <cardinality_1>, [<x1 marginals>], <cardinality_2>, [<x2 marginals>], ...

    with open(filename, 'r') as f:
        lines = [x.strip() for x in f.readlines() if x.strip()]
    assert len(lines) == 2
    assert lines[0] == 'MAR'
    fields = lines[1].strip().split()
    n_vars = int(fields[0])

    res = []
    for i, field in enumerate(fields[1:]):
        if i == 0 or len(res[-1]) == cardinality:
            cardinality = int(field)
            res.append([])
        else:
            res[-1].append(float(field))

    assert len(res) == n_vars
    return res



def load_decomposed_graph(filename) -> DecomposedGraph:
    """
    Load a decomposed graph from a .dfg file.
    """
    with open(filename, 'r') as f:
        n_factors = int(read_next_line(f))
        dg = []
        for i in range(n_factors):
            n_terms = int(read_next_line(f))
            weights = np.array([float(x) for x in read_next_line(f).split()])
            n_vars = int(read_next_line(f))
            vars = [int(x) for x in read_next_line(f).split()]
            cardinalities = [int(x) for x in read_next_line(f).split()]
            assert len(vars) == n_vars
            assert len(cardinalities) == n_vars
            factor_matrices = []
            for var, cardinality in zip(vars, cardinalities):
                n_nonzero = int(read_next_line(f))
                table_flat = np.zeros(cardinality * n_terms)
                for entry in range(n_nonzero):
                    var_index, val = read_next_line(f).split()
                    var_index = int(var_index)
                    val = float(val)
                    table_flat[var_index] = val
                factor_matrices.append(table_flat.reshape((cardinality, n_terms), order='F'))
            dg.append(DecomposedFactor(vars, weights, factor_matrices))

    return DecomposedGraph(dg)


def factorise_ising(t):
    # Obtained by solving {a^2 + b^2 = t, 2ab = 1/t} for a, b

    flip = False
    if t < 1:
        t = 1/t
        flip = True

    x = np.sqrt((t**3 - np.sqrt(t**2 * (t**4 - 1))) / t**2)
    a = x / np.sqrt(2)
    b = (np.sqrt(2) * x * t**3 + np.sqrt(2) * np.sqrt(t**2 * (t**4 - 1)) * x) / (2 * t)

    if flip:
        G = np.array([[a, b],
                      [b, a]])
        H = np.array([[b, a],
                      [a, b]])
    else:
        G = np.array([[a, b],
                      [b, a]])
        H = np.array([[a, b],
                      [b, a]])

    return [G, H]

def ising_g(N, unary_min, unary_max, pairwise_min, pairwise_max):
    """
    NxN Ising model with random potentials uniformly chosen within the given limits.
    :return: Tuple (g, dg) containing the ising model as a Graph and DecomposedGraph instance respectively. The
    tensor decomposition used is described in Wrigley, Lee and Ye (2017), supplementary material.
    """
    g = []
    dg = []
    for var in range(N*N):
        # Unary potential
        theta = np.random.uniform(unary_min, unary_max)
        g.append(Factor([var], np.exp([theta, -theta])))
        dg.append(DecomposedFactor([var], np.array([1.0]), [np.exp([[theta], [-theta]])]))

        # Pairwise potentials
        edges = []
        if var % N != 0:
            # Not first col
            edges.append([var-1, var])
        if var >= N:
            # Not first row
            edges.append([var-N, var])

        for vars in edges:

            # Non-decomposed
            theta = np.random.uniform(pairwise_min, pairwise_max)
            f = Factor(vars, np.exp([[theta, -theta], [-theta, theta]]))

            # Decomposed
            G, H = factorise_ising(np.exp(theta))
            df = DecomposedFactor(vars, np.array([1.0, 1.0]), [G, H])

            # Check decomposed potential is equal to original
            assert df.expand() == f

            g.append(f)
            dg.append(df)

    return Graph(g), DecomposedGraph(dg)


def random_g(N, edge_prob, unary_min, unary_max, pairwise_min, pairwise_max):
    """
    Connected N-variable MRF with Ising potentials uniformly chosen within the given limits. Each possible edge is
    present with independent probability edge_prob.
    :return: Tuple (g, dg) containing the MRF as a Graph and DecomposedGraph instance respectively. The tensor
    decomposition used is described in Wrigley, Lee and Ye (2017), supplementary material.
    """
    g = []
    dg = []
    while not g or len(Graph(g).get_connected_components()) > 1:
        for i in range(N):

            # Unary potential
            theta = np.random.uniform(unary_min, unary_max)
            g.append(Factor([i], np.exp([theta, -theta])))
            dg.append(DecomposedFactor([i], np.array([1.0]), [np.exp([[theta], [-theta]])]))

            for j in range(i + 1, N):

                if np.random.rand() < edge_prob:

                    # Pairwise potential
                    theta = np.random.uniform(pairwise_min, pairwise_max)
                    f = Factor([i, j], np.exp([[theta, -theta], [-theta, theta]]))
                    G, H = factorise_ising(np.exp(theta))
                    df = DecomposedFactor([i, j], np.array([1.0, 1.0]), [G, H])

                    # Check decomposed potential is equal to original
                    assert df.expand() == f

                    g.append(f)
                    dg.append(df)

    return Graph(g), DecomposedGraph(dg)


def l1_error(marg_1, marg_2, binary=False):
    """
    Mean L1 marginal error.
    :param marg_1: List of lists containing marginals
    :param marg_2: List of lists containing marginals
    :param binary: Set to true if computing errors for binary variables to return the average error in P(X=0), ignoring
    P(X=1) (no difference to result, just slightly more efficient).
    :return: Mean error over all marginals.
    """

    # Results sets must contain the same number of variables
    assert len(marg_1) == len(marg_2)

    err = 0
    n = 0
    for i in range(len(marg_1)):
        # Variable i must contain the same number of marginals in both result sets
        assert len(marg_1[i]) == len(marg_2[i])
        if binary:
            assert len(marg_1[i]) == 2
            err += np.abs(marg_1[i][0] - marg_2[i][0])
            n += 1
        else:
            err += np.sum([np.abs(marg_1[i][j] - marg_2[i][j]) for j in range(len(marg_1[i]))])
            n += len(marg_1[i])
    return err / n


def format_mar(marg):
    """
    :param marg: List of lists containing marginals
    :return: string in .MAR format (http://www.hlt.utdallas.edu/~vgogate/uai14-competition/resformat.html)
    """
    fields = [len(marg)]
    for m in marg:
        fields.append(len(m))
        fields += m
    return ' '.join(str(x) for x in fields)




