# Tensor Belief Propagation 


[Tensor Belief Propagation](http://proceedings.mlr.press/v70/wrigley17a/wrigley17a.pdf) (TBP) is an algorithm for
approximate inference in discrete graphical models [1]. At a high level, it takes a factor graph in [.uai](#other-file-formats) or [.fg](#other-file-formats) format and outputs approximate marginals for
each variable.

This code should currently be considered experimental.

[1] [Wrigley, Andrew, Wee Sun Lee, and Nan Ye. "Tensor Belief Propagation." International Conference on Machine Learning. 2017.](http://proceedings.mlr.press/v70/wrigley17a/wrigley17a.pdf)

## Requirements

 * Linux or OSX
 * Python 3 (tested with 3.6)
 
## Installation

Install with the Python package manager `pip`:
```bash
$ pip install tbp
TODO Example output
```

## Usage:
### Command line
```bash
$ tbp in.any -o out.any
```

### Python library

```python
import tbp

# Load a factor graph in .uai format
g = tbp.load_uai_graph('tests/uai/MAR_prob/linkage_11.uai')

# Apply evidence (fixed variable assignments)
g.apply_evidence('tests/uai/MAR_prob/linkage_11.evid')

# Decompose each factor into a weighted sum of 4 rank-1 tensors
dg = g.decompose(r=4)

# Run TBP to find marginals with sample size of 10000
mar = dg.tbp_marg(K=10000)

```



## Computing marginals for decomposed graphs

A "decomposed graph" is a factor graph with all potential functions represented as sums of rank-1 tensors (rather than 
multidimensional tables). We use the [.dfg file format](#dfg-file-format) to represent decomposed graphs.

If you have a graph in this form, the command line utility `utils/dfgmarg` can be used to compute approximate
marginals:

```bash
$ utils/dfgmarg tests/ising_8x8.dfg 1000
64 2 0.607513 0.392487 2 ... 0.633576 0.366424
```
This runs TBP on the graph `tests/ising_8x8.dfg` with a sample size of `K = 1000`, and outputs marginal estimates for
 each variable in [.MAR format](#other-file-formats). `K` controls the accuracy and
running time of the algorithm - a higher `K` gives more accurate marginals at the cost of increased running time.

## Computing marginals for arbitrary factor graphs


To use TBP to estimate marginals for an arbitrary factor graph, the graph first needs to be decomposed. This requires
installation of the included Python scripts.

### Installation

First, make sure you have the Python packaging tool [Pipenv](https://github.com/pypa/pipenv) installed on your
system (i.e. `pip install pipenv`). Then, install the Python dependencies:

```bash
$ pipenv install
Installing dependencies from Pipfile.lock (1c8ebb)…
...
```


### Usage

The python utility `marg.py` allows you to run TBP on arbitrary factor graphs. It uses the
[TensorLy library](https://tensorly.github.io/stable/index.html) to perform a
non-negative [CP decomposition](https://en.wikipedia.org/wiki/Tensor_rank_decomposition) for each
potential function in the input graph, and then invokes TBP to find marginals.
Supported input file formats for factor graphs are [libDAI's .fg format](#other-file-formats) and
[UAI format](#other-file-formats). Some example usages are as follows.


Find marginals for the factor graph `ising_8x8.fg`:
```bash
$ pipenv run python python/tbp/marg.py tests/ising_8x8.fg
64 2 0.594961 0.405039 2 ... 0.608573 0.391427
```

Find marginals for the factor graph `Promedus_11.uai` including evidence, with sample size 1000, and save the output
to `out.MAR`:
```bash
$ pipenv run python python/tbp/marg.py tests/uai/MAR_prob/Promedus_11.uai \
-e tests/uai/MAR_prob/Promedus_11.uai.evid -k 1000 -o out.MAR --verbosity 2
Reading graph tests/uai/MAR_prob/Promedus_11.uai (UAI format)...
Applying evidence file tests/uai/MAR_prob/Promedus_11.uai.evid...
Decomposing input graph (r=4 terms per factor)...
Running TBP with sample size K=1000...
Successfully saved marginals to out.MAR.
```

Perform the decomposition with 3 rank-1 components and save the resulting decomposed graph:
```bash
$ pipenv run python python/tbp/marg.py tests/ising_8x8.fg -d -r 3 \
-o tests/ising_8x8.dfg --verbosity 2
Reading graph tests/ising_8x8.fg (libDAI format)...
Decomposing input graph (r=3 terms per factor)...
Successfully saved decomposed graph to tests/ising_8x8.dfg.
```

For documentation on the full range of options, see `pipenv run python python/tbp/marg.py --help`.

### Using MATLAB for the decomposition
As an alternative to TensorLy, the [MATLAB Tensor Toolbox](http://www.sandia.gov/~tgkolda/TensorToolbox) can be used
to perform the initial tensor decompositions (this was what was used in [1]). To use this instead of Tensorly:

 * Install MATLAB normally
 * Install the
 [MATLAB Python API](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html)
 * Install the [MATLAB Tensor Toolbox](http://www.sandia.gov/~tgkolda/TensorToolbox)
 
You can now replace `method='tensorly'` with `method='matlab'` in [core.py](python/tbp/core.py) wherever you wish.

## ICML experiments
The results from [1] can be reproduced with:

```
$ pipenv run python python/tbp/icml17.py
```
Note that these tests take considerable time to finish.

## File formats

Descriptions of file formats used in this project are as follows.

### .dfg (decomposed factor graph)
We created the `.dfg` file format based on
[libDAI's .fg file format](https://staff.fnwi.uva.nl/j.m.mooij/libDAI/doc/fileformats.html)
to represent decomposed factor graphs. A decomposed factor graph is a
factor graph with all factors represented as sums of rank-1 tensors rather than multidimensional tables.

The first line of a `.dfg` file contains the number of factors in the graph, followed by a blank line. Then, factors
are described in turn by blocks separated by a single blank line. Each factor block is structured as follows:

     1. n_terms
     2. <weights>
     3. n_variables
     4. <variable indices>
     5. <variable cardinalities>
     6. n_nonzero_1
     7. 1 0.5
     8. 3 0.1
     9. 4 0.1
    10. ...
    11. n_nonzero_2
    12. 1 0.5
    13. 3 0.1
    14. 4 0.3
    15. ...

In the header section of the factor block (lines 1-5), `n_terms` is the number of terms in the decomposition and
`<weights>`, `<variable indices>` and `<variable cardinalities>` are self-explanatory space-separated lists of length `n_terms`,
`n_variables` and `n_variables` respectively. 

The remainder of the factor block (line 6 onwards) describes
a series of `n_variables` 2D matrices that together describe the `n_terms` rank-1 tensors.
Each matrix corresponds to a single variable and has shape `(cardinality, n_terms)`, where `cardinality` is
the cardinality of the variable and `n_terms` is the number of rank-1 terms in the decomposition (constant
for all variables). Each matrix begins with the
number of nonzero values in the matrix, followed by a series of `index value` pairs describing the nonzero
entries of the matrix in column-major order. See
[libDAI's documentation](https://staff.fnwi.uva.nl/j.m.mooij/libDAI/doc/fileformats.html) for examples of how to
reshape these lists back into matrices.

 The __i__th rank-1 tensor is constructed by taking the outer product of the __i__th columns of
all matrices. The complete factor is then reconstructed by adding up these rank-1 tensors and weighting
according to `<weights>`.

### Other file formats
Other file formats used in this project are:

 * `.fg` (libDAI factor graph): https://staff.fnwi.uva.nl/j.m.mooij/libDAI/doc/fileformats.html
 * `.uai` (UAI factor graph): http://www.hlt.utdallas.edu/~vgogate/uai14-competition/modelformat.html
 * `.MAR` (marginals): http://www.hlt.utdallas.edu/~vgogate/uai14-competition/resformat.html
 * `.evid` (evidence): http://www.hlt.utdallas.edu/~vgogate/uai14-competition/evidformat.html

## To do
 * Add experiments from Figure 2 (random pairwise MRFs)
 * Document all file formats
 * Rewrite code that loads .uai files to handle all problems (currently breaks on some)
 * Deal with Z <= 0 warning from C++ code
 * Add Python tests to Travis

## Feedback
For bug reports, suggestions and comments please email [andrew@wrigley.io](mailto:andrew@wrigley.io) or use the issue tracker.
If you find TBP useful we would love to know.

## License

TODO (also link to libDAI and Eigen licences)

## Acknowledgments

* [libDAI](https://staff.fnwi.uva.nl/j.m.mooij/libDAI/) (modified version included in [libdai](libdai) folder, used extensively for core junction tree algorithm)
* [Eigen](http://eigen.tuxfamily.org/) (version 3.3.4 included in [libdai/vendor/include](libdai/vendor/include) folder)
* [TensorLy](https://github.com/tensorly/tensorly)
* [MATLAB Tensor Toolbox](http://www.sandia.gov/~tgkolda/TensorToolbox)





