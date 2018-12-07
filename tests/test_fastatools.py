# -*- coding: utf-8 -*-

"""
test_fastatools
----------------------------------

Tests for `fastatools` module.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import logging
import pytest

from fastatools import fastatools


def write_fasta(seq_strings, tmpdir, file_name, defline_prefix=""):
    """Write sequences to a temporary fasta file.

    Parameters
    ----------
    seq_strings : list of str
        Sequences of nucleotide characters to write to the fasta file
    tmpdir : py.path.local
        Pytest tmpdir fixture
    file_name : str
        File name of the fasta file to write
    defline_prefix : str, optional
        Prefix prepended to id, name, description.  Defaults to "".

    Returns
    -------
    str
        Path of fasta file
    """
    records = []
    for i, seq_string in enumerate(seq_strings):
        seq_num_str = defline_prefix + str(1 + i)
        ident = "Id" + seq_num_str
        name = "Name" + seq_num_str
        description = "Description" + seq_num_str
        seq = Seq(seq_string)
        record = SeqRecord(seq, ident, name, description)
        records.append(record)

    file_path = str(tmpdir.join(file_name))
    SeqIO.write(records, file_path, "fasta")
    return file_path


def test_verify_file_exists_when_exists(tmpdir, caplog):
    """Verify no error message and no system exit upon existing file."""
    existing_file = write_fasta(['A'], tmpdir, "existing")
    fastatools.verify_file_exists(existing_file)
    assert(caplog.text == "")


def test_verify_file_exists_when_missing(tmpdir, caplog):
    """Verify error message and system exit upon non-existing file."""
    missing_file = str(tmpdir.join("missing"))
    with pytest.raises(SystemExit):
        fastatools.verify_file_exists(missing_file)
        assert("does not exist" in caplog.text)


def test_length(tmpdir, capsys):
    """Verify correct chrom length calculations."""
    seq_strings = ['A' * 1, 'T' * 10, 'C' * 100, 'G' * 1000]
    path1 = write_fasta(seq_strings, tmpdir, "A1.T10.C100.G1000")
    seq_strings = ['A' * 2, 'T' * 20, 'C' * 200, 'G' * 2000]
    path2 = write_fasta(seq_strings, tmpdir, "A2.T20.C200.G2000")
    fastatools.length([path1, path2])
    captured = capsys.readouterr()
    assert("A1.T10.C100.G1000 1 Id1" in captured.out)
    assert("A1.T10.C100.G1000 10 Id2" in captured.out)
    assert("A1.T10.C100.G1000 100 Id3" in captured.out)
    assert("A1.T10.C100.G1000 1000 Id4" in captured.out)
    assert("A2.T20.C200.G2000 2 Id1" in captured.out)
    assert("A2.T20.C200.G2000 20 Id2" in captured.out)
    assert("A2.T20.C200.G2000 200 Id3" in captured.out)
    assert("A2.T20.C200.G2000 2000 Id4" in captured.out)


def test_equivalent_different_num_sequences(tmpdir, capsys):
    """Verify the correct message when different number of sequences."""
    seq_strings = ['A' * 1, 'T' * 10, 'C' * 100]
    path1 = write_fasta(seq_strings, tmpdir, "A1.T10.C100")
    seq_strings = ['A' * 1, 'T' * 10, 'C' * 100, 'G' * 1000]
    path2 = write_fasta(seq_strings, tmpdir, "A2.T20.C200.G2000")
    fastatools.equivalent(path1, path2)
    captured = capsys.readouterr()
    assert(captured.out == "Not equivalent -- the number of sequences is different (3 and 4).\n")


def test_equivalent_different_deflines(tmpdir, capsys):
    """Verify the correct message when different sequence descriptions."""
    seq_strings = ['A' * 1]
    path1 = write_fasta(seq_strings, tmpdir, "file1", '1')
    seq_strings = ['A' * 1]
    path2 = write_fasta(seq_strings, tmpdir, "file2", '2')
    fastatools.equivalent(path1, path2, ignore_defline=False)
    captured = capsys.readouterr()
    assert(captured.out == 'Not equivalent -- sequence 1 has different descriptions ("Id11 Description11" and "Id21 Description21").\n')


def test_equivalent_different_deflines_ignore(tmpdir, capsys):
    """Verify the correct message when different sequence descriptions, but deflines are ignored."""
    seq_strings = ['A' * 1]
    path1 = write_fasta(seq_strings, tmpdir, "file1", '1')
    seq_strings = ['A' * 1]
    path2 = write_fasta(seq_strings, tmpdir, "file2", '2')
    fastatools.equivalent(path1, path2, ignore_defline=True)
    captured = capsys.readouterr()
    assert(captured.out == "Equivalent\n")


def test_equivalent_different_lengths(tmpdir, capsys):
    """Verify the correct message when different sequence lengths."""
    seq_strings = ['A' * 1]
    path1 = write_fasta(seq_strings, tmpdir, "file1")
    seq_strings = ['A' * 2]
    path2 = write_fasta(seq_strings, tmpdir, "file2")
    fastatools.equivalent(path1, path2, ignore_defline=False)
    captured = capsys.readouterr()
    assert(captured.out == 'Not equivalent -- sequence 1 has different lengths (1 and 2).\n')


def test_equivalent_different_contents(tmpdir, capsys):
    """Verify the correct message when different sequence contents."""
    seq_strings = ['A' * 10]
    path1 = write_fasta(seq_strings, tmpdir, "file1")
    seq_strings = ['T' * 10]
    path2 = write_fasta(seq_strings, tmpdir, "file2")
    fastatools.equivalent(path1, path2, ignore_defline=False)
    captured = capsys.readouterr()
    assert(captured.out == 'Not equivalent -- sequence 1 has different contents.\n')


def test_equivalent_enforce_order(tmpdir, capsys):
    """Verify the sequences must be in the same order"""
    seq_strings = ['A' * 10, 'T' * 100]
    path1 = write_fasta(seq_strings, tmpdir, "file1")
    seq_strings = ['T' * 100, 'A' * 10]
    path2 = write_fasta(seq_strings, tmpdir, "file2")
    fastatools.equivalent(path1, path2, ignore_defline=True, enforce_order=True)
    captured = capsys.readouterr()
    assert(captured.out == 'Not equivalent -- sequence 1 has different lengths (10 and 100).\n'
                           'Not equivalent -- sequence 2 has different lengths (100 and 10).\n')


def test_equivalent_ignore_order(tmpdir, capsys):
    """Verify the sequences can be in different order"""
    seq_strings = ['A' * 10, 'T' * 100]
    path1 = write_fasta(seq_strings, tmpdir, "file1")
    seq_strings = ['T' * 100, 'A' * 10]
    path2 = write_fasta(seq_strings, tmpdir, "file2")
    fastatools.equivalent(path1, path2, ignore_defline=True, enforce_order=False)
    captured = capsys.readouterr()
    assert(captured.out == "Equivalent\n")


def test_rewrite_line_length(tmpdir, capsys):
    """Verify line length is fixed."""
    seq_strings = ['A']
    path1 = write_fasta(seq_strings, tmpdir, "file1")
    with open(path1, 'w') as f:  # overwrite
        f.write(">1\n")
        f.write("a" * 20 + "\n")
        f.write("a" * 20 + "\n")
        f.write("a" * 20 + "\n")
        f.write("t" * 20 + "\n")
    fastatools.rewrite(path1, force_upper=False)
    captured = capsys.readouterr()
    lines = captured.out.split('\n')
    assert(lines[1] == "a" * 60)
    assert(lines[2] == "t" * 20)


def test_rewrite_upper(tmpdir, capsys):
    """Verify line length is fixed."""
    seq_strings = ['A']
    path1 = write_fasta(seq_strings, tmpdir, "file1")
    with open(path1, 'w') as f:  # overwrite
        f.write(">1\n")
        f.write("a" * 20 + "\n")
        f.write("a" * 20 + "\n")
        f.write("a" * 20 + "\n")
        f.write("t" * 20 + "\n")
    fastatools.rewrite(path1, force_upper=True)
    captured = capsys.readouterr()
    lines = captured.out.split('\n')
    assert(lines[1] == "A" * 60)
    assert(lines[2] == "T" * 20)


def test_reverse(tmpdir, capsys):
    """Verify the correct reverse complement."""
    seq_strings = ["ATCG"]
    path1 = write_fasta(seq_strings, tmpdir, "file1")
    fastatools.reverse(path1)
    captured = capsys.readouterr()
    lines = captured.out.split('\n')
    assert(lines[1] == "CGAT")


def test_between(tmpdir, capsys):
    """Verify extraction between primers."""
    seq_strings = ['A' * 20 + "AATTCCGGA" + "GATACA" + "AATTCCGGT" + 'A' * 20]
    path1 = write_fasta(seq_strings, tmpdir, "file1")
    fastatools.between(path1, "AATTCCGGA", "AATTCCGGT", no_reverse_complement=True)
    captured = capsys.readouterr()
    assert(captured.out == "AATTCCGGA" + "GATACA" + "AATTCCGGT\n")


def test_between_not_found(tmpdir, caplog):
    """Verify message when missing primers."""
    caplog.set_level(logging.INFO)
    seq_strings = ['A' * 20 + "AATTCCGGA" + "GATACA" + "AATTCCGGT" + 'A' * 20]
    path1 = write_fasta(seq_strings, tmpdir, "file1")
    fastatools.between(path1, "TTTTTT", "GGGGGG", no_reverse_complement=True)
    assert(caplog.text == "Forward primer not found.\nReverse primer not found.\n")


def test_range(tmpdir, capsys):
    """Verify extracting a range of positions."""
    seq_strings = ["ATTCCGGA"]
    path1 = write_fasta(seq_strings, tmpdir, "file1")
    fastatools.range_command(path1, "Id1", 2, 7)
    captured = capsys.readouterr()
    assert(captured.out == "TTCCGG\n")
