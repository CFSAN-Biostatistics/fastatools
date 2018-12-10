========
Usage
========

.. highlight:: bash

To print the lengths of sequences::

    fastatools length FILE [FILE ...]

To determine if two fasta files are equivalent, ignoring the line lengths, uppercase / lowercase characters, and sequence order::

    fastatools equiv [--order] [--ignore_defline] FILE1 FILE2

To fix inconsistent line lengths and inconsistent lowercase::

    fastatools rewrite --upper FILE

To generate a reverse complement of a fasta file::

    fastatools reverse FILE

To extract a slice from a fasta file delimited by primers::

    fastatools between --no_rev_comp FWD REV FILE

To extract a range of positions from a fasta file::

    fastatools range id start end FILE

