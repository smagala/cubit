#!/usr/bin/env python
#file sequence.py

# Copyright (C) 2008 James Smagala

# Contact:
#   James Smagala
#   smagala@gmail.com

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

"""Utility functions for handling sequence data"""

import re


CODONS = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L',
          'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'AUU': 'I', 'AUC': 'I',
          'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
          'GUG': 'V', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
          'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACU': 'T',
          'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A',
          'GCA': 'A', 'GCG': 'A', 'UAU': 'Y', 'UAC': 'Y', 'UAA': '*',
          'UAG': '*', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
          'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D',
          'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'UGU': 'C', 'UGC': 'C',
          'UGA': '*', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R',
          'CGG': 'R', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
          'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', '---': '-'}

DNA = "ACGT"
RNA = "ACGU"
GAPS = '-'

IUPAC_DNA = {'R': ('A', 'G'), 'Y': ('C', 'T'), 'S': ('C', 'G'),
             'W': ('A', 'T'), 'K': ('G', 'T'), 'M': ('A', 'C'),
             'B': ('C', 'G', 'T'), 'D': ('A', 'G', 'T'),
             'H': ('A', 'C', 'T'), 'V': ('A', 'C', 'G'),
             'N': ('A', 'C', 'G', 'T')}

DNA_COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'R': 'Y', 'Y': 'R',
            'W': 'W', 'S': 'S', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
            'D': 'H', 'H': 'D', '-': '-'}

IUPAC_RNA = {'R': ('A', 'G'), 'Y': ('C', 'U'), 'S': ('C', 'G'),
             'W': ('A', 'U'), 'K': ('G', 'U'), 'M': ('A', 'C'),
             'B': ('C', 'G', 'U'), 'D': ('A', 'G', 'U'),
             'H': ('A', 'C', 'U'), 'V': ('A', 'C', 'G'),
             'N': ('A', 'C', 'G', 'U')}

RNA_COMP = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A', 'R': 'Y', 'Y': 'R',
            'W': 'W', 'S': 'S', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
            'D': 'H', 'H': 'D', '-': '-'}

PROTEIN = "ACDEFGHIKLMNPQRSTVWY*"

PROTEIN_DEGENERICIES = {'X': ('X')}

NEWICK_CHARS_RE = re.compile(r"['\"\[\]\(\),:;]")


class SequenceError(Exception):
    """Base Exception for all sequence objects"""
    pass


class DNASequenceError(Exception):
    """Exception for DNASequence Objects"""
    pass


class RNASequenceError(Exception):
    """Exception for RNASequence Objects"""
    pass


class Sequence(object):
    """Base class for biological sequence data"""
    def __init__(self, name='', seq=''):
        """Initalize a sequence object"""
        super(Sequence, self).__init__()
        # sequences have an additional name property
        self.name = name
        self._seq = seq
        self.bioedit_re = re.compile("(.*?)\s+\d+\sbases\Z")

    def upper(self):
        """Make the sequence uppercase"""
        self._seq = self._seq.upper()

    def lower(self):
        """Make the sequence lowercase"""
        self._seq = self._seq.lower()

    def to_fasta(self):
        """Get a fasta string representation of a Sequence"""
        return ">%s\n%s" % (self.name, self._seq)

    def to_string(self):
        """Get a string representation of just the sequence"""
        return self._seq

    def __str__(self):
        """The string representation of a Sequence is just the sequence"""
        return self._seq

    def __len__(self):
        return len(self._seq)

    def __repr__(self):
        return self.to_fasta()

    def __iter__(self):
        return iter(self._seq)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self._seq[index.start, index.stop, index.step]
        else:
            return self._seq[index]

    def clean_sequence(self):
        """Ensure that the sequence is valid"""
        self._seq = self._seq.upper()
        self._seq = self.bad_char_re.sub('-', self._seq)

    def clean_bioedit_comment(self):
        """Get rid of bioedit junk in comments"""
        match = self.bioedit_re.match(self.name)
        if match:
            self.name = match.group(1)

    def is_gapped(self):
        """Check that the sequence contains gap characters"""
        for gap in self.gaps:
            if gap in self._seq:
                return True
        return False

    def is_ungapped(self):
        """Check that the sequence contains no gap characters"""
        for gap in self.gaps:
            if gap in self._seq:
                return False
        return True

    def is_valid(self):
        """The sequence is valid if it contains no unrecognized chars"""
        for c in self._seq:
            if not c in self.alphabet:
                return False
        return True

    def is_ambiguous(self):
        """The sequence is ambiguous if it has one or more ambiguous chars"""
        for c in self._seq:
            if c in self.degenericies:
                return True
        return False


class ProteinSequence(Sequence):
    """Represents a Protein sequence"""

    def __init__(self, name='', seq=''):
        """Initalize a ProteinSequence object"""
        super(ProteinSequence, self).__init__(name, seq)
        self.gaps = GAPS
        self.alphabet = PROTEIN
        self.degenericies = PROTEIN_DEGENERICIES
        self.valid_chars = self.gaps + self.alphabet + \
                           ''.join(self.degenericies.keys())
        self.bad_char_re = re.compile(r"[^%s]" % self.valid_chars)

        self.clean_sequence()


class RNASequence(Sequence):
    """Represents a RNA Sequence"""

    def __init__(self, name='', seq=''):
        """Initalize a RNASequence object"""
        super(RNASequence, self).__init__(name, seq)
        self.gaps = GAPS
        self.alphabet = RNA
        self.degenericies = IUPAC_RNA
        self.valid_chars = self.gaps + self.alphabet + \
                           ''.join(self.degenericies.keys())
        self.codons = CODONS
        self.bad_char_re = re.compile(r"[^%s]" % self.valid_chars)
        self.t2u_re = re.compile(r"[tT]")

        self.clean_sequence()

    def clean_sequence(self):
        """Convert T->U before doing any further cleanup"""
        self._seq = self.t2u_re.sub('U', self._seq)
        super(RNASequence, self).clean_sequence()

    def disambiguate(self):
        """Return a list of non-ambiguous RNA sequences.

        Uses a new, faster algorithm to handle highly-degenerate seqs.
        """
        try:
            for seq in disambiguate(self._seq, self.alphabet + self.gaps, \
                                    self.degenericies):
                yield RNASequence(self.name, seq)
        except SequenceError, msg:
            raise RNASequenceError(msg)

    def to_dna(self):
        """Convert this sequence to DNA"""
        return DNASequence(self.name, self._seq)

    def translate(self):
        """Translate RNA -> protein"""
        start = 0
        stop = len(self)
        step = 3
        protein = []
        for pos in range(start, stop, step):
            codon = self._seq[pos:pos+3]
            if codon in self.codons:
                protein.append(self.codons[codon])
            elif len(codon) < 3:
                continue
            else:
                codons = [c for c in disambiguate(codon, \
                        self.alphabet + self.gaps, \
                        self.degenericies)]
                aa_set = set([self.codons.get(codon, 'X') for codon in codons])
                if len(aa_set) == 1:
                    protein.append(aa_set.pop())
                else:
                    protein.append('X')
        return ProteinSequence(self.name, ''.join(protein))

    def translate_disambiguate(self):
        """Translate ambiguous RNA, returning a list of AA's by position"""
        start = 0
        stop = len(self)
        step = 3
        protein = []
        for pos in range(start, stop, step):
            degen_codon = self._seq[pos:pos + 3]
            codons = [c for c in disambiguate(degen_codon, \
                                              self.alphabet + self.gaps, \
                                              self.degenericies)]
            aa_set = set([self.codons.get(codon, 'X') for codon in codons])
            protein.append(aa_set)
        return protein


class DNASequence(Sequence):
    """Represents a DNA Sequence"""

    def __init__(self, name='', seq=''):
        """Initalize a DNASequence object"""
        super(DNASequence, self).__init__(name, seq)
        self.gaps = GAPS
        self.alphabet = DNA
        self.degenericies = IUPAC_DNA
        self.valid_chars = self.gaps + self.alphabet + \
                           ''.join(self.degenericies.keys())
        self.bad_char_re = re.compile(r"[^%s]" % self.valid_chars)
        self.u2t_re = re.compile(r"[uU]")

        self.clean_sequence()

    def clean_sequence(self):
        """Convert U->T before doing any further cleanup"""
        self._seq = self.u2t_re.sub('T', self._seq)
        super(DNASequence, self).clean_sequence()

    def disambiguate(self):
        """Return a list of non-ambiguous DNA sequences.

        Uses a new, faster algorithm to handle highly-degenerate seqs.
        """
        try:
            for seq in disambiguate(self._seq, self.alphabet + self.gaps, \
                                    self.degenericies):
                yield DNASequence(self.name, seq)
        except SequenceError, msg:
            raise DNASequenceError(msg)

    def to_rna(self):
        """Convert the sequence to RNA"""
        return RNASequence(self.name, self._seq)

    def translate(self):
        """Translate the sequence to get amino acids"""
        return self.to_rna().translate()

    def translate_disambiguate(self):
        """Get a list of all possible amino acids at every position"""
        return self.to_rna().translate_disambiguate()


def disambiguate(seq, alphabet, degenericies):
    """Generate all possible non-ambiguous sequences"""
    # build an array containing all choices for all sequence positions
    degen_seq = []
    pos_count = []
    current_indicies = [0] * len(seq)    
    for char in seq:
        if char in alphabet:
            char_set = (char, )
        elif char in degenericies:
            char_set = degenericies[char]
        else:
            raise SequenceError("Unknown Character: '%s'" % char)
        degen_seq.append(char_set)
        pos_count.append(len(char_set))

    # loop until all possible combinations have been generated
    while True:
        # generate a sequence
        unambiguous_seq = ''.join([char_set[index] for char_set, index in \
                                   zip(degen_seq, current_indicies)])
        yield unambiguous_seq
        for i in range(len(pos_count)):
            # increment the current position, then take the mod of the
            # length of the list of choices - if the number rolls over,
            # continue on to the next position, if not, generate a new
            # sequence
            carryover, current_indicies[i] = divmod(current_indicies[i] + 1,
                                                    pos_count[i])
            if not carryover:
                break
        else:
            # break the while loop when all possible sequences have been
            # generated
            break


def diff2seqs(reference, seq):
    """Compare two sequences, and record the differences"""

    def _find_start(s):
        for n, c in enumerate(s):
            if not c.issubset(set('-X')):
                break
        return n

    def _find_stop(s):
        return len(s) - 1 - _find_start(s[::-1])

    # find the start and end of actual data - this way we aren't marking
    # missing data as a deletion
    start = max(_find_start(reference), _find_start(seq))
    stop = min(_find_stop(reference), _find_stop(seq))

    changes = {}
    for pos, pair in enumerate(zip(reference, seq)):
        if not start <= pos <= stop:
            continue
        if not pair[0] == pair[1]:
            changes[pos] = pair[1]
    return changes


def line_len(line, length=79):
    """Format long lines to a given length"""
    lines = []
    while line:
        lines.append(line[:length])
        line = line[length:]
    return '\n'.join(lines)


class FastaFormatError(Exception):
    """Error raised by parse_fasta for poorly formatted records"""
    pass


def format_fasta(comment, seq, length=79):
    """Format comment, sequence pairs into fasta records."""
    return ">%s\n%s" % (comment, line_len(seq.upper()))


def parse_fasta(infile, bioedit_cleanup=False):
    """Convert FASTA files into a comment,sequence pair.

    Requires a file-like object as input, in FASTA form:
    >comment
    sequence

    Yields (comment,sequence) tuples.

    """
    # get the first comment line
    comment = infile.readline()
    if comment.startswith('>'):
        comment = comment[1:].strip()
    else:
        raise FastaFormatError("First character of FASTA file must be '>'")

    seq = []
    bioedit_re = re.compile('\s+\d+\sbases\Z')

    def compile_seq(seq):
        """Build a complete sequence from a list of strings"""
        # ensure that the sequence data exists
        seq = ''.join(seq)
        if not seq:
            raise FastaFormatError("FASTA record contains no sequence data")
        return seq

    for line in infile:
        # check for the start of a new record
        if line.startswith('>'):
            if bioedit_cleanup:
                comment = bioedit_re.split(comment)[0]
            seq = compile_seq(seq)
            yield comment, seq
            seq = []
            comment = line[1:].strip()
            continue
        # add the line to the sequence data, removing white space
        seq.append(''.join(line.split()))

    if bioedit_cleanup:
        comment = bioedit_re.split(comment)[0]
    seq = compile_seq(seq)

    yield comment, seq


def read_fasta_file(path, bioedit_cleanup=False):
    """Open, read and close a fasta file, returning a list of (comment, seq)"""
    try:
        infile = open(path, 'r')
        records = [record for record in parse_fasta(infile, bioedit_cleanup)]
    finally:
        try:
            # ensure that the file is always closed
            infile.close()
        except:
            pass
    return records


def clean_for_tree(comment):
    """clean up any bad newick characters"""
    comment = '_'.join(comment.split('.'))
    comment = NEWICK_CHARS_RE.sub('', comment)
    return '_'.join(comment.split())
