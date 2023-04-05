#!/usr/bin/env python
# file fix_frame.py

# Copyright (C) 2005 James Smagala

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

"""Concatenate ORFs for MP, NS genes

Revision History
    20050916 James Smagala: Written
    20090128 JS: updated to work with B viruses
"""

from optparse import OptionParser
import sys
from os.path import splitext
import re
from core.sequence import parse_fasta, format_fasta, FastaFormatError, \
        DNASequence


def splice_pb1(seq):
    """Carefully splice PB1 and PB1-F2 if it exists"""
    # start with the full length PB1 gene
    result = [seq[0:2274], ]

    def translate_fragment(seq):
        """Translate a DNA string into a protein string"""
        s = DNASequence('', seq)
        return s.translate()._seq

    # check if the first potential stop codon was a stop
    # if it was, don't show any more amino acids after it
    stop1 = translate_fragment(seq[2271:2274])
    stop2 = translate_fragment(seq[2274:])
    if stop1 == '*' and stop2 != '*':
        result.append('---')
    else:
        result.append(seq[2274:])

    # check for the presence of a start codon for the F2 protein
    # if it is missing, refuse to guess
    start2 = translate_fragment(seq[94:97])
    #if start2 == 'M':
    #    f2_aa = translate_fragment(seq[94:]).split('-')[0]
    #    try:
    #        index = f2_aa.index('*') + 1
    #    except ValueError:
    #        index = len(f2_aa)
    #    result.append(seq[94:94+index*3])

    if start2 == 'M':
        result.append(seq[94:94+102*3])
    else:
        result.append(seq[94:94+102*3])

    return ''.join(result)


if __name__ == "__main__":

    def parse_command_line():
        """Get and check command line options and arguments."""
        usage = "usage: %prog FILE"
        version = "%prog 0.1"
        description = ("Concatenate ORFs for any of flu A/B, genes PB1/MP/NS.",
                       "Genes must be truncated to contain only the ORF.",
                       "A PB1 genes must be 2277 nt (include *both* stop codons).",
                       "A MP genes must be 982 nt.",
                       "A NS genes must be 838 nt.",
                       "B MP genes must be 1076 nt.",
                       "B NS genes must be 1024 nt.", )
        description = '\n'.join(description)
        parser = OptionParser(usage=usage, version=version, \
                              description=description)
        args = parser.parse_args()[1]

        if len(args) < 1:
            parser.error("Missing input file")
        infile = ''.join(args[0].split('"'))
        return infile

    def main():
        """Run as a script"""
        infile = parse_command_line()
        root = splitext(infile)[0]
        outfile = root + "_orf.fas"

        try:
            infile = open(infile, 'r')
            outfile = open(outfile, 'w')
        except IOError, msg:
            print msg
            sys.exit(1)

        try:
            for comment, seq in parse_fasta(infile, bioedit_cleanup=True):
                if len(seq) == 2277:
                    seq = splice_pb1(seq)
                elif len(seq) == 982:
                    seq = ''.join([seq[0:759], seq[0:26], seq[714:982], ])
                elif len(seq) == 838:
                    seq = ''.join([seq[0:693], seq[0:30], seq[502:838], ])
                elif len(seq) == 1076:
                    seq = ''.join([seq[0:747], seq[746:1076], ])
                elif len(seq) == 1024:
                    seq = ''.join([seq[0:846], seq[0:33], seq[688:1024], ])
                elif len(seq) == 1027:
                    seq = ''.join([seq[0:849], seq[0:36], seq[691:1027], ])
                else:
                    print "Sequence '%s' does not appear to be PB1, MP or NS.  " \
                          "(length: %s, expected: 2277, 982, 838, 1076, 1024, 1027)" \
                          % (comment, len(seq))
                    sys.exit(1)
                print >> outfile, format_fasta(comment, seq)
        except FastaFormatError, msg:
            print msg
            sys.exit(1)

        outfile.close()

    main()
