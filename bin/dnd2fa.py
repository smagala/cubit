#!/usr/bin/env python
# file dnd2fa.py

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

"""Parse .dnd files, get accession numbers, and retrieve fasta records

Assumes that all sequences of interest are in a single large fasta file

Revision History
  050418 James Smagala: Rework of parse_fa
  050818 JS: Cleaned up, use optparse, fasta parser copied from confind.py
  051108 JS: Fix to work with longer names

"""

from core.sequence import parse_fasta, FastaFormatError, line_len
from optparse import OptionParser
import sys
import re


if __name__ == "__main__":

    def parse_command_line():
        """Get and check command line options and arguments."""
        usage = "usage: %prog [OPTIONS] SOURCE_FASTA_FILE DND_FILE"
        version = "%prog 0.3"
        description = "Create a FASTA file containing all sequences with an " \
            "identifier found in both SOURCE_FASTA_FILE and DND_FILE."
        parser = OptionParser(usage=usage, version=version, \
                              description=description)

        parser.add_option("-o", "--output", dest="outfile",
                          help="output file location [default: %default]")
        parser.set_defaults(outfile="out.fas")

        options, args = parser.parse_args()

        try:
            dbfile = ''.join(args[0].split('"'))
            infile = ''.join(args[1].split('"'))
        except IndexError:
            parser.error("Missing input file")
        outfile = ''.join(options.outfile.split('"'))

        return infile, outfile, dbfile

    def main():
        """Executed when run as a script"""
        infile, outfile, dbfile = parse_command_line()

        try:
            infile = open(infile, 'r')
            dbfile = open(dbfile, 'r')
            outfile = open(outfile, 'w')
        except IOError, msg:
            print msg

        # parse the fa flat db into a dict
        try:
            seq_db = dict([record for record in parse_fasta(dbfile)])
        except FastaFormatError, msg:
            print msg
            sys.exit(1)

        tree = ''.join([line.strip() for line in infile])
        exp = re.compile(r',[^)(:,;]*?:|\([^)(:,;]*?:')
        labels = [label[1:-1] for label in exp.findall(tree)]

        for label in labels:
            try:
                print >> outfile, ">%s\n%s" % (label, line_len(seq_db[label]))
            except KeyError:
                print "Label '%s' is not in database" % label

        # clean up open files
        dbfile.close()
        infile.close()
        outfile.close()

    main()
