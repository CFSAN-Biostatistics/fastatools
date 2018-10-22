#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import argparse
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

__version__ = '1.0.0'

def parse_arguments(system_args):
    """
    Parse command line arguments.

    Parameters
    ----------
    system_args : list
        List of command line arguments, usually sys.argv[1:].

    Returns
    -------
    Namespace
        Command line arguments are stored as attributes of a Namespace.
    """
    description = """Fasta Tools."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--version',          action='version', version='%(prog)s version ' + __version__)
    subparsers = parser.add_subparsers(dest="subparser_name", help=None, metavar="subcommand       ")
    subparsers.required = True

    formatter_class = argparse.ArgumentDefaultsHelpFormatter

    description = "Print the lengths of sequences."
    subparser = subparsers.add_parser("length", formatter_class=formatter_class, description=description, help=description)
    subparser.add_argument(dest="fasta_path", type=str, metavar="FILE", help="Fasta file.", nargs='+')
    subparser.set_defaults(func=length)

    description = "Generate a reverse complement of a fasta file."
    subparser = subparsers.add_parser("reverse", formatter_class=formatter_class, description=description, help=description)
    subparser.add_argument(dest="fasta_path", type=str, metavar="FILE", help="Fasta file.")
    subparser.set_defaults(func=reverse)

    description = "Extract a slice from a fasta file delimited by primers."
    subparser = subparsers.add_parser("slice", formatter_class=formatter_class, description=description, help=description)
    subparser.add_argument(type=str, dest="fwd_primer", metavar="FWD", help="Forward primer.")
    subparser.add_argument(type=str, dest="rev_primer", metavar="REV", help="Reverse primer.")
    subparser.add_argument(type=str, dest="fasta_path", metavar="FILE", help="Fasta file.")
    subparser.add_argument("--no_rev_comp", action="store_true", default=False, dest="no_reverse_complement", help="Do not automatically reverse complement the reverse primer.")
    subparser.set_defaults(func=slice)

    description = "Extract a range of positions from a fasta file."
    subparser = subparsers.add_parser("range", formatter_class=formatter_class, description=description, help=description)
    subparser.add_argument(type=str, dest="contig", metavar="id", help="Fasta contig id.")
    subparser.add_argument(type=int, dest="start", help="Starting position.")
    subparser.add_argument(type=int, dest="end", help="Ending position.")
    subparser.add_argument(type=str, dest="fasta_path", metavar="FILE", help="Fasta file.")
    subparser.set_defaults(func=range)

    args = parser.parse_args(system_args)
    return args


def length(args):
    """Read a fasta file and print the lengths of all the sequences.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv
    """
    seqrecords = []
    for path in args.fasta_path:
        for seqrecord in SeqIO.parse(path, "fasta"):
            print(path, len(seqrecord.seq), seqrecord.id);


def reverse(args):
    """Read a fasta file and write its reverse complement to stdout.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv
    """
    seqrecords = []
    for seqrecord in SeqIO.parse(args.fasta_path, "fasta"):
        seqrecord.seq = seqrecord.seq.reverse_complement()
        seqrecord.description = seqrecord.description + " reverse complement"
        seqrecords.append(seqrecord)
    SeqIO.write(seqrecords, sys.stdout, "fasta")


def slice(args):
    """Extract a slice from a fasta file.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv
    """
    fwd_primer = args.fwd_primer.upper()
    rev_primer = args.rev_primer.upper()

    if not args.no_reverse_complement:
        rev_primer = Seq(rev_primer, generic_dna).reverse_complement()

    fwd_found = False
    rev_found = False
    for seqrecord in SeqIO.parse(args.fasta_path, "fasta"):
        seq = seqrecord.seq.upper()
        findex = seq.find(fwd_primer)
        if findex == -1:
            continue
        fwd_found = True
        rindex = seq.find(rev_primer)
        if rindex == -1:
            continue
        rev_found = True
        print(seqrecord.seq[findex : rindex + len(rev_primer)])
        return

    if not fwd_found:
        print("Forward primer not found.", file=sys.stderr)
    if not rev_found:
        print("Reverse primer not found.", file=sys.stderr)


def range(args):
    """Extract a range of positions from a fasta file.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv
    """
    contig = args.contig
    start_pos = args.start
    end_pos = args.end

    contig_found = False
    for seqrecord in SeqIO.parse(args.fasta_path, "fasta"):
        if seqrecord.id == contig:
            seq = seqrecord.seq
            print(seqrecord.seq[start_pos : 1 + end_pos])
            return

    print("Contig %s not found." % contig, file=sys.stderr)


def main():
    args = parse_arguments(sys.argv[1:])
    args.func(args)


if __name__ == '__main__':
    main()

