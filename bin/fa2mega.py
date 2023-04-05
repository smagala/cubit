#!/usr/bin/env python
# file fa2mega.py

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

"""Convert fasta formatted files to mega formatted files.

History
    20080211 James Smagala: File created.

"""

import os.path
from optparse import OptionParser
from core.sequence import parse_fasta, FastaFormatError


def to_mega(title, seqs):
    """Convert sequences to a MEGA formatted string"""
    lines = ["#Mega", "Title: %s" % title, '']
    for seq in seqs:
        lines.append("#%s\n%s" % seq)

    return '\n'.join(lines)


if __name__ == "__main__":

    def parse_cmd_line():
        """Get command line options and arguments"""
        usage = "usage: %prog [OPTIONS] INFILE.fas"
        version = "%prog 0.1"
        description = "Convert FASTA formatted file to MEGA formatted file"
        parser = OptionParser(usage=usage, version=version,
                              description=description)

        #opts, args = parser.parse_args()
        args = parser.parse_args()[1]

        if len(args) < 1:
            parser.error("Missing input file")

        seqfiles = (''.join(f.split('"')) for f in args)

        return seqfiles

    def main():
        """Handles filenames if run as a script from the command line"""
        seqfiles = parse_cmd_line()
        for seqfile in seqfiles:
            infile = open(seqfile, 'r')
            try:
                seqs = [seq for seq in parse_fasta(infile)]
            except FastaFormatError, msg:
                print msg
            finally:
                infile.close()
            megafile = os.path.splitext(seqfile)[0] + ".meg"
            outfile = open(megafile, 'w')
            # use the seqfile name as the title in the mega file
            print >> outfile, to_mega(os.path.basename(seqfile), seqs)
            outfile.close()

    main()
