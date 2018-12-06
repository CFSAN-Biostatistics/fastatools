# -*- coding: utf-8 -*-

"""This module is part of fastatools.
"""

from __future__ import print_function
from __future__ import absolute_import


from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import logging
import sys


def length(fasta_paths):
    """Read fasta files and print the lengths of all the sequences.

    Parameters
    ----------
    fasta_paths : list of str
        List of fasta file paths to process
    """
    for path in fasta_paths:
        for seqrecord in SeqIO.parse(path, "fasta"):
            print(path, len(seqrecord.seq), seqrecord.id)


def equivalent(fasta_path1, fasta_path2, ignore_defline=False):
    """Determine if two fasta files are equivalent.  This comparison
    ignores sequence line lengths and uppercase / lowercase.

    Parameters
    ----------
    fasta_path1 : str
        Fasta file path.
    fasta_path2 : str
        Fasta file path.
    ignore_defline : bool, optional
        When true, ignore the sequence description lines.  Defaults to False.
    """
    seqrecords1 = [seqrecord for seqrecord in SeqIO.parse(fasta_path1, "fasta")]
    seqrecords2 = [seqrecord for seqrecord in SeqIO.parse(fasta_path2, "fasta")]
    if len(seqrecords1) != len(seqrecords2):
        print("Not equivalent -- the number of sequences is different (%d and %d)." % (len(seqrecords1), len(seqrecords2)))
        return 1

    equiv = True

    if not ignore_defline:
        for i, seq1, seq2 in zip(range(1, 1+len(seqrecords1)), seqrecords1, seqrecords2):
            if seq1.description != seq2.description:
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


def rewrite(fasta_path, force_upper=True):
    """Fix inconsistent line lengths and uppercase lowercase.  Write to stdout.

    Parameters
    ----------
    fasta_path : str
        Fasta file path.
    force_upper : bool, optional
        Convert sequences to uppercase. Defaults to True.
    """
    seqrecords = []
    for seqrecord in SeqIO.parse(fasta_path, "fasta"):
        if force_upper:
            seqrecord.seq = seqrecord.seq.upper()
        seqrecords.append(seqrecord)
    SeqIO.write(seqrecords, sys.stdout, "fasta")


def reverse(fasta_path):
    """Read a fasta file and write its reverse complement to stdout.

    Parameters
    ----------
    fasta_path : str
        Fasta file path.
    """
    seqrecords = []
    for seqrecord in SeqIO.parse(fasta_path, "fasta"):
        seqrecord.seq = seqrecord.seq.reverse_complement()
        seqrecord.description = seqrecord.description + " reverse complement"
        seqrecords.append(seqrecord)
    SeqIO.write(seqrecords, sys.stdout, "fasta")


def between(fasta_path, fwd_primer, rev_primer, no_reverse_complement=False):
    """Extract a slice from a fasta file.

    Parameters
    ----------
    fasta_path : str
        Fasta file path.
    fwd_primer : str
        Forward primer.
    rev_primer : str
        Reverse primer.
    no_reverse_complement : bool, optional
        Do not automatically reverse complement the reverse primer.  Defaults to False.
    """
    fwd_primer = fwd_primer.upper()
    rev_primer = rev_primer.upper()

    if not no_reverse_complement:
        rev_primer = Seq(rev_primer, generic_dna).reverse_complement()

    fwd_found = False
    rev_found = False
    for seqrecord in SeqIO.parse(fasta_path, "fasta"):
        seq = seqrecord.seq.upper()
        findex = seq.find(fwd_primer)
        if findex == -1:
            continue
        fwd_found = True
        rindex = seq.find(rev_primer)
        if rindex == -1:
            continue
        rev_found = True
        print(seqrecord.seq[findex: rindex + len(rev_primer)])
        return

    if not fwd_found:
        print("Forward primer not found.", file=sys.stderr)
    if not rev_found:
        print("Reverse primer not found.", file=sys.stderr)


def range_command(fasta_path, contig, start_pos, end_pos):
    """Extract a range of positions from a fasta file.

    Parameters
    ----------
    fasta_path : str
        Fasta file path.
    contig : str
        Fasta contig id.
    start_pos : int
        Starting position (1-based).
    end_pos : int
        Ending position (1-based).
    """
    for seqrecord in SeqIO.parse(fasta_path, "fasta"):
        if seqrecord.id == contig:
            print(seqrecord.seq[start_pos - 1: end_pos])
            return

    print("Contig %s not found." % contig, file=sys.stderr)
