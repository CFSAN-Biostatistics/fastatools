# -*- coding: utf-8 -*-

"""This module is part of fastatools.
"""

from __future__ import print_function
from __future__ import absolute_import


from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import collections
import logging
import os
import sys


FastaSequence = collections.namedtuple("FastaSequence", ["length", "id", "seq", "original_order"])


def verify_file_exists(path):
    """Exit with error message if specified file does not exist.

    Parameters
    ----------
    path : str
        File path.
    """
    if not os.path.isfile(path):
        logging.error("Error: the file %s does not exist." % path)
        sys.exit(1)


def length(fasta_paths):
    """Read fasta files and print the lengths of all the sequences.

    Parameters
    ----------
    fasta_paths : list of str
        List of fasta file paths to process
    """
    for path in fasta_paths:
        if not os.path.isfile(path):
            logging.info("Error: the file %s does not exist." % path)
            continue

        for seqrecord in SeqIO.parse(path, "fasta"):
            print(path, len(seqrecord.seq), seqrecord.id)


def equivalent(fasta_path1, fasta_path2, ignore_defline=False, enforce_order=False):
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
    enforce_order : bool, optional
        When true, require sequences to be in the same order.  Defaults to False.
    """
    verify_file_exists(fasta_path1)
    verify_file_exists(fasta_path2)

    seqrecords1 = [seqrecord for seqrecord in SeqIO.parse(fasta_path1, "fasta")]
    seqrecords2 = [seqrecord for seqrecord in SeqIO.parse(fasta_path2, "fasta")]

    if len(seqrecords1) != len(seqrecords2):
        print("Not equivalent -- the number of sequences is different (%d and %d)." % (len(seqrecords1), len(seqrecords2)))
        return 1

    FastaSequence = collections.namedtuple("FastaSequence", ["length", "description", "seq", "original_order"])
    seqrecords1 = [FastaSequence(len(seq.seq), seq.description, seq.seq, i + 1) for i, seq in enumerate(seqrecords1)]
    seqrecords2 = [FastaSequence(len(seq.seq), seq.description, seq.seq, i + 1) for i, seq in enumerate(seqrecords2)]

    if not enforce_order:
        # sort by sequence length, then description
        seqrecords1 = sorted(seqrecords1)
        seqrecords2 = sorted(seqrecords2)

    equiv = True

    if not ignore_defline:
        for seq1, seq2 in zip(seqrecords1, seqrecords2):
            if seq1.description != seq2.description:
                print('Not equivalent -- sequence %d has different descriptions ("%s" and "%s").' % (seq1.original_order, seq1.description, seq2.description))
                equiv = False
        if not equiv:
            return 1

    for seq1, seq2 in zip(seqrecords1, seqrecords2):
        if len(seq1.seq) != len(seq2.seq):
            print("Not equivalent -- sequence %d has different lengths (%d and %d)." % (seq1.original_order, len(seq1.seq), len(seq2.seq)))
            equiv = False
    if not equiv:
        return 1

    for seq1, seq2 in zip(seqrecords1, seqrecords2):
        if seq1.seq.upper() != seq2.seq.upper():
            print("Not equivalent -- sequence %d has different contents." % (seq1.original_order))
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
    verify_file_exists(fasta_path)

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
    verify_file_exists(fasta_path)

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
    verify_file_exists(fasta_path)

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
        logging.info("Forward primer not found.")
    if not rev_found:
        logging.info("Reverse primer not found.")


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
    verify_file_exists(fasta_path)

    for seqrecord in SeqIO.parse(fasta_path, "fasta"):
        if seqrecord.id == contig:
            print(seqrecord.seq[start_pos - 1: end_pos])
            return

    logging.info("Contig %s not found." % contig)
