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

    help = "Determine if two sequences are equivalent."
    description = "Determine if two sequences are equivalent, ignoring the line lengths and uppercase / lowercase characters."
    subparser = subparsers.add_parser("equiv", formatter_class=formatter_class, description=description, help=help)
    subparser.add_argument(dest="fasta_path1", type=str, metavar="FILE", help="Fasta file 1.")
    subparser.add_argument(dest="fasta_path2", type=str, metavar="FILE", help="Fasta file 2.")
    subparser.add_argument("--ignore_defline", action="store_true", default=False, dest="ignore_defline", help="Ignore the sequence description lines.")
    subparser.set_defaults(func=equivalent)

    description = "Fix inconsistent line lengths and inconsistent lowercase."
    subparser = subparsers.add_parser("rewrite", formatter_class=formatter_class, description=description, help=description)
    subparser.add_argument(dest="fasta_path", type=str, metavar="FILE", help="Fasta file.")
    subparser.add_argument("--upper", action="store_true", default=False, dest="force_upper", help="Convert sequences to uppercase.")
    subparser.set_defaults(func=rewrite)

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
    subparser.set_defaults(func=range_command)

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


def equivalent(args):
    """Determine if two fasta files are equivalent.  This comparison
    ignores sequence line lengths and uppercase / lowercase.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv
    """
    seqrecords1 = [seqrecord for seqrecord in SeqIO.parse(args.fasta_path1, "fasta")]
    seqrecords2 = [seqrecord for seqrecord in SeqIO.parse(args.fasta_path2, "fasta")]
    if len(seqrecords1) != len(seqrecords2):
        print("Not equivalent -- the number of sequences is different (%d and %d)." % (len(seqrecords1), len(seqrecords2)))
        return 1

    equiv = True

    if not args.ignore_defline:
        for i, seq1, seq2 in zip(range(1, 1+len(seqrecords1)), seqrecords1, seqrecords2):
            if len(seq1.description) != len(seq2.description):
                print('Not equivalent -- sequence %d has different descriptions ("%s" and "%s").' % (i, seq1.description, seq2.description))
                equiv = False
        if not equiv:
            return 1

    for i, seq1, seq2 in zip(range(1, 1+len(seqrecords1)), seqrecords1, seqrecords2):
        if len(seq1.seq) != len(seq2.seq):
            print("Not equivalent -- sequence %d has different lengths (%d and %d)." % (i, len(seq1.seq), len(seq2.seq)))
            equiv = False
    if not equiv:
        return 1

    for i, seq1, seq2 in zip(range(1, 1+len(seqrecords1)), seqrecords1, seqrecords2):
        if seq1.seq.upper() != seq2.seq.upper():
            print("Not equivalent -- sequence %d has different contents." % (i))
            equiv = False
    if not equiv:
        return 1

    print("Equivalent")


def rewrite(args):
    """Fix inconsistent line lengths and uppercase lowercase.  Write to stdout.

    Parameters
    ----------
    args : Namespace
        Command line arguments stored as attributes of a Namespace, usually
        parsed from sys.argv
    """
    seqrecords = []
    for seqrecord in SeqIO.parse(args.fasta_path, "fasta"):
        if args.force_upper:
            seqrecord.seq = seqrecord.seq.upper()
        seqrecords.append(seqrecord)
    SeqIO.write(seqrecords, sys.stdout, "fasta")


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


def range_command(args):
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
    return args.func(args)


if __name__ == '__main__':
    main()

