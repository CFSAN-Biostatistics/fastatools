#!/usr/bin/env python

"""Console script for fastatools."""

from __future__ import print_function
from __future__ import absolute_import

import argparse
import logging
import sys

from fastatools import fastatools
from fastatools.__init__ import __version__


def parse_arguments(system_args):
    """Parse command line arguments.

    Parameters
    ----------
    system_args : list
        List of command line arguments, usually sys.argv[1:].

    Returns
    -------
    Namespace
        Command line arguments are stored as attributes of a Namespace.
    """
    def non_negative_int(value):
        try:
            ivalue = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError("Must be a number >= 0")
        if ivalue < 0:
            raise argparse.ArgumentTypeError("Must be >= 0")
        return ivalue

    def positive_int(value):
        try:
            ivalue = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError("Must be a number greater than 0")
        if ivalue <= 0:
            raise argparse.ArgumentTypeError("Must be greater than 0")
        return ivalue

    description = """Fasta Tools."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--version', action='version', version='%(prog)s version ' + __version__)
    subparsers = parser.add_subparsers(dest="subparser_name", help=None, metavar="subcommand")
    subparsers.required = True

    formatter_class = argparse.ArgumentDefaultsHelpFormatter

    description = "Print the lengths of sequences."
    subparser = subparsers.add_parser("length", formatter_class=formatter_class, description=description, help=description)
    subparser.add_argument(dest="fasta_paths", type=str, metavar="FILE", help="Fasta file.", nargs='+')
    subparser.set_defaults(func=length_command)

    help_str = "Determine if two fasta files are equivalent."
    description = "Determine if two fasta files are equivalent, ignoring the line lengths, uppercase / lowercase characters, and sequence order."
    subparser = subparsers.add_parser("equiv", formatter_class=formatter_class, description=description, help=help_str)
    subparser.add_argument(dest="fasta_path1", type=str, metavar="FILE", help="Fasta file 1.")
    subparser.add_argument(dest="fasta_path2", type=str, metavar="FILE", help="Fasta file 2.")
    subparser.add_argument("--order", action="store_true", default=False, dest="enforce_order", help="Require same sequence order when there are multiple sequences.")
    subparser.add_argument("--ignore_defline", action="store_true", default=False, dest="ignore_defline", help="Ignore the sequence description lines.")
    subparser.set_defaults(func=equivalent_command)

    description = "Fix inconsistent line lengths and inconsistent lowercase."
    subparser = subparsers.add_parser("rewrite", formatter_class=formatter_class, description=description, help=description)
    subparser.add_argument(dest="fasta_path", type=str, metavar="FILE", help="Fasta file.")
    subparser.add_argument("--upper", action="store_true", default=False, dest="force_upper", help="Convert sequences to uppercase.")
    subparser.set_defaults(func=rewrite_command)

    description = "Generate a reverse complement of a fasta file."
    subparser = subparsers.add_parser("reverse", formatter_class=formatter_class, description=description, help=description)
    subparser.add_argument(dest="fasta_path", type=str, metavar="FILE", help="Fasta file.")
    subparser.set_defaults(func=reverse_command)

    description = "Extract a slice from a fasta file delimited by primers."
    subparser = subparsers.add_parser("between", formatter_class=formatter_class, description=description, help=description)
    subparser.add_argument(type=str, dest="fwd_primer", metavar="FWD", help="Forward primer.")
    subparser.add_argument(type=str, dest="rev_primer", metavar="REV", help="Reverse primer.")
    subparser.add_argument(type=str, dest="fasta_path", metavar="FILE", help="Fasta file.")
    subparser.add_argument("--no_rev_comp", action="store_true", default=False, dest="no_reverse_complement", help="Do not automatically reverse complement the reverse primer.")
    subparser.set_defaults(func=between_command)

    description = "Extract a range of positions from a fasta file."
    subparser = subparsers.add_parser("range", formatter_class=formatter_class, description=description, help=description)
    subparser.add_argument(type=str, dest="contig", metavar="id", help="Fasta contig id.")
    subparser.add_argument(type=positive_int, dest="start", help="Starting position (1-based).")
    subparser.add_argument(type=positive_int, dest="end", help="Ending position (1-based).")
    subparser.add_argument(type=str, dest="fasta_path", metavar="FILE", help="Fasta file.")
    subparser.set_defaults(func=range_command)

    args = parser.parse_args(system_args)
    return args


def length_command(args):
    """Print the lengths of sequences.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv
    """
    return fastatools.length(args.fasta_paths)


def equivalent_command(args):
    """Determine if two sequences are equivalent, ignoring the line lengths and uppercase / lowercase characters.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv
    """
    return fastatools.equivalent(args.fasta_path1, args.fasta_path2, ignore_defline=args.ignore_defline, enforce_order=args.enforce_order)


def rewrite_command(args):
    """Fix inconsistent line lengths and inconsistent lowercase.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv
    """
    return fastatools.rewrite(args.fasta_path, force_upper=args.force_upper)


def reverse_command(args):
    """Generate a reverse complement of a fasta file.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv
    """
    return fastatools.reverse(args.fasta_path)


def between_command(args):
    """Extract a slice from a fasta file delimited by primers.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv
    """
    return fastatools.between(args.fasta_path, args.fwd_primer, args.rev_primer, no_reverse_complement=args.no_reverse_complement)


def range_command(args):
    """Extract a range of positions from a fasta file.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv
    """
    return fastatools.range_command(args.fasta_path, args.contig, args.start, args.end)


def run_command_from_args(args):
    """Run a subcommand with previously parsed arguments in an argparse namespace.

    This function is intended to be used for unit testing.

    Parameters
    ----------
    args : Namespace
        Command line arguments are stored as attributes of a Namespace.
        The args are obtained by calling parse_argument_list().

    Returns
    -------
    Returns 0 on success if it completes with no exceptions.
    """
    return args.func(args)  # this executes the function previously associated with the subparser with set_defaults


def run_from_line(line):
    """Run a command with a command line.

    This function is intended to be used for unit testing.

    Parameters
    ----------
    line : str
        Command line.

    Returns
    -------
    Returns 0 on success if it completes with no exceptions.
    """
    argv = line.split()
    args = parse_arguments(argv)
    return args.func(args)  # this executes the function previously associated with the subparser with set_defaults


def main():
    """This is the main function which is turned into an executable
    console script by the setuptools entry_points.  See setup.py.

    To run this function as a script, first install the package:
        $ python setup.py develop
        or
        $ pip install --user fastatools

    Parameters
    ----------
    This function must not take any parameters

    Returns
    -------
    The return value is passed to sys.exit().
    """
    enable_log_timestamps = False
    if enable_log_timestamps:
        logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S", level=logging.INFO)
    else:
        logging.basicConfig(format="%(message)s", level=logging.INFO)
    args = parse_arguments(sys.argv[1:])
    return args.func(args)  # this executes the function previously associated with the subparser with set_defaults


# This snippet lets you run the cli without installing the entrypoint.
if __name__ == "__main__":
    sys.exit(main())
