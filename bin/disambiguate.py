#!/usr/bin/env python
# file disambiguate.py

# Copyright (C) 2006 James Smagala

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

"""Create multiple non-ambiguous sequences for each ambiguous one

Revision History
  20060724 James Smagala: file creation

"""

import sys
from optparse import OptionParser
from core.sequence import parse_fasta, FastaFormatError, DNASequence
import os.path
import re


def disambiguate(infile):
    """Disambiguate sequences in a file, returning multiple near-copies"""
    # get the base file name for later use
    basefile = os.path.splitext(os.path.basename(infile))[0]

    infile = open(infile, 'r')
    outfile = open(basefile + "_unamb.fas", 'w')
    for comment, seq in parse_fasta(infile, bioedit_cleanup=True):
        amb_seq = DNASequence(comment, seq)
        amb_seq.upper()
        # don't modify unambiguous sequences
        if not amb_seq.is_ambiguous():
            print >> outfile, amb_seq.to_fasta()
            continue
        # for ambiguous sequences, append a count for easier tracking
        for label, unamb_seq in enumerate(amb_seq.disambiguate()):
            unamb_seq.name = '%s_%d' % (unamb_seq.name, label + 1)
            print >> outfile, unamb_seq.to_fasta()

    infile.close()
    outfile.close()


if __name__ == "__main__":

    def parse_command_line():
        """Get and check command line options and arguments."""
        usage = "usage: %prog INFILE"
        version = "%prog 0.1"
        description = "Disambiguate all sequences"
        parser = OptionParser(usage=usage, version=version, \
                              description=description)
        args = parser.parse_args()[1]

        try:
            infile = ''.join(args[0].split('"'))
        except IndexError:
            parser.error("Missing input")
            sys.exit(1)
        return infile

    def main():
        """Keep the module namespace clean"""
        infile = parse_command_line()
        try:
            disambiguate(infile)
        except IOError, msg:
            print msg
        except FastaFormatError, msg:
            print msg

    main()
