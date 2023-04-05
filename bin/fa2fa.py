#!/usr/bin/env python
# file fa2fa.py

# Copyright (C) 2005 Regents of the University of Colorado

# Contact:
#   James Smagala
#   smagala@colorado.edu

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

"""Read many .fa files to build one large .fa for use as a local BLAST db

Can also exclude seqs.

Revision History
  050426 James Smagala: Rework of dnd2fa
  050627 JS: Complete overhaul

"""

from optparse import OptionParser
import sys
from core.sequence import parse_fasta, FastaFormatError, line_len


if __name__ == "__main__":

    def parse_command_line():
        """Get and check command line options and arguments."""
        usage = "usage: %prog [OPTIONS] FILE1 FILE2 FILE3..."
        version = "%prog 0.2"
        description = "Compile multiple fasta files into a larger file " \
                      "and/or exclude sequences from a fasta file."
        parser = OptionParser(usage=usage, version=version, \
                              description=description)

        parser.add_option("-e", "--exclude", dest="exfile", metavar="FILE",
                          type="string",
                          help="path to a fasta FILE of sequences to be "
                          "excluded from the output")
        parser.set_defaults(exfile="")
        parser.add_option("-o", "--output", dest="outfile", metavar="FILE",
                          type="string",
                          help="path to the output FILE to be created or "
                               "overwritten [default: %default]")
        parser.set_defaults(outfile="out.fas")

        options, args = parser.parse_args()

        # Get rid of any quotes, expand any wildcards
        infiles = [''.join(infile.split('"')) for infile in args]
        if not infiles:
            parser.error("Missing input file")
        outfile = ''.join(options.outfile.split('"'))
        exfile = ''.join(options.exfile.split('"'))

        return infiles, exfile, outfile

    def main():
        """The program to run if this is run as a script"""
        infiles, exfile, outfile = parse_command_line()

        # get the identifiers for the exceptions
        ex_db = {}
        if exfile:
            try:
                exfile = open(exfile, 'r')
            except IOError, msg:
                print msg
                sys.exit(1)
            for comment, seq in parse_fasta(exfile):
                ex_db[comment] = None
            exfile.close()

        # open the output file for writing
        try:
            outfile = open(outfile, 'w')
        except IOError, msg:
            print msg
            sys.exit(1)

        # loop through all input files
        for next_file in infiles:
            try:
                infile = open(next_file, 'r')
            except IOError, msg:
                print msg
                continue

            try:
                for comment, seq in parse_fasta(infile):
                    if not comment in ex_db:
                        print >> outfile, ">%s\n%s" % (comment, line_len(seq))
                infile.close()
            except FastaFormatError, msg:
                print "%s: %s" % (next_file, msg)
                infile.close()

        # clean up open files
        outfile.close()

    main()
