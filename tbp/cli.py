#!/usr/bin/env python3

"""
Command line utilities for Tensor Belief Propagation. See https://github.com/akxlr/tbp.

Author: Andrew Wrigley, National University of Singapore and Australian National University
"""

import os
import sys
from . import core
from .core import status
import argparse


class HelpfulParser(argparse.ArgumentParser):
    """
    Display full help message if arguments are not supplied correctly.
    """
    def error(self, message):
        self.print_help(sys.stderr)
        sys.stderr.write('\nerror: {}\n'.format(message))
        sys.exit(2)


def _marg(in_file, r=None, K=None, evid_file=None, out_file=None, exact=False, verbosity=None):

    if r is not None:
        r = int(r)

    if K is not None:
        K = int(K)

    if verbosity:
        core.verbosity = int(verbosity)

    in_extn = os.path.splitext(in_file)[1]
    out_extn = os.path.splitext(out_file)[1] if out_file else None
    if out_extn == '.mar':
        out_file = os.path.splitext(out_file)[0] + '.MAR'
        out_extn = '.MAR'


    if in_extn not in ['.uai', '.fg', '.dfg']:
        raise ValueError("Unknown graph file format: {} (supported: .fg, .uai, .dfg)".format(in_extn))
    if out_extn and out_extn not in ['.dfg', '.MAR']:
        raise ValueError("Unknown output format: {} (supported: .dfg, .MAR)".format(out_extn))

    if in_extn == '.dfg':

        assert not evid_file, "Cannot apply evidence to .dfg - supply the original graph instead"
        assert not exact, "Cannot compute exact marginals for .dfg - supply the original graph instead"
        status("Reading decomposed graph {}...".format(in_file), 2)
        dg = core.load_decomposed_graph(in_file)

    else:

        if in_extn == '.fg':
            status("Reading graph {} (libDAI format)...".format(in_file), 2)
            g = core.load_fg_graph(in_file)
        elif in_extn == '.uai':
            status("Reading graph {} (UAI format)...".format(in_file), 2)
            g = core.load_uai_graph(in_file)

        if evid_file:
            status("Applying evidence file {}...".format(evid_file), 2)
            g.apply_evidence(evid_file)

        if not exact:
            status("Decomposing input graph (r={} terms per factor)...".format(r), 2)
            dg = g.decompose(r)


    if out_extn == '.dfg':
        dg.save_to_file(out_file)
        status("Successfully saved decomposed graph to {}.".format(out_file), 2)
        return

    if exact:
        status("Computing exact marginals...".format(r), 2)
        mar = g.exact_marg_elim()

    else:
        status("Running TBP with sample size K={}...".format(K), 2)
        mar = dg.tbp_marg(K)

    if out_extn == '.MAR':
        with open(out_file, 'w') as f:
            f.write(core.format_mar(mar))
            status("Successfully saved marginals to {}.".format(out_file), 2)

    else:
        print(core.format_mar(mar))


def main():

    parser = HelpfulParser(
        description="Compute approximate marginals for a factor graph using Tensor Belief Propagation. For more "
                    "information and a description of the various file formats, see the project homepage at "
                    "{}.".format(core.HOMEPAGE),
    )

    parser.add_argument(
        'in_file',
        help="File describing a factor graph (*.uai, *.fg) or decomposed factor graph (*.dfg)"
    )

    parser.add_argument(
        '-v',
        '--verbosity',
        required=False,
        help="How much progress info to show (1=none, 2=some, 3=all) (default: 1)",
        default=1
    )

    parser.add_argument(
        '-r',
        '--rank',
        required=False,
        help="Number of rank-1 terms to decompose each initial potential function into (default: {})".format(core.DEFAULT_COMPONENTS),
        default=core.DEFAULT_COMPONENTS,
    )

    parser.add_argument(
        '-k',
        '--sample-size',
        required=False,
        help="Sample size K to use for sampling step of TBP. Higher K gives more accurate marginals at the cost of "
             "increased running time. (default: {})".format(core.DEFAULT_SAMPLE_SIZE),
        default=core.DEFAULT_SAMPLE_SIZE,
    )

    parser.add_argument(
        '-e',
        '--evidence',
        required=False,
        help="File containing variable assignments to apply before running inference (*.evid)."
    )

    parser.add_argument(
        '-o',
        '--output',
        help="File to save output to. The action taken depends on the file extension - supply *.MAR to save marginals,"
             "or *.dfg to save the decomposed graph and not compute marginals. If not supplied, marginals are printed "
             "to stdout in MAR format."
    )
    parser.add_argument(
        '-x',
        '--exact',
        action='store_true',
        help="Instead of computing approximate marginals with TBP, compute exact marginals by running the elimination "
             "algorithm on each variable in turn (very inefficient, use only for very small graphs)."
    )

    args = vars(parser.parse_args())
    _marg(in_file=args['in_file'], r=args['rank'], K=args['sample_size'], evid_file=args['evidence'],
          out_file=args['output'], exact=args['exact'], verbosity=args['verbosity'])






