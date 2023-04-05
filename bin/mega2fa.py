#!/usr/bin/env python
# file mega2fa.py

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

"""Clean FASTA record names to correct format for Becky

Revision History
    20050912 James Smagala: Written
    20080328 JS: Cleaned up

"""

from optparse import OptionParser
import os.path
import sys
import re


if __name__ == "__main__":

    def parse_command_line():
        """Get and check command line options and arguments."""
        usage = "usage: %prog FILE"
        version = "%prog 0.1"
        description = "Convert MEGA files to FASTA."
        parser = OptionParser(usage=usage,version=version,description=description)
        options, args = parser.parse_args()

        try:
            infile = ''.join(args[0].split('"'))
        except IndexError:
            parser.error("Missing input file")
        return infile

    def main():
        """Avoid cluttering the global namespace."""
        infile = parse_command_line()
        root,ext = os.path.splitext(infile)
        outfile = root + ".fas"

        try:
            infile = open(infile,'r')
            outfile = open(outfile,'w')
        except IOError,msg:
            print msg
            sys.exit(1)

        fa = re.compile(r'#')
        start = re.compile(r'#mega',re.I)
        for line in infile:
            if start.match(line):
                continue
            if line.startswith('#'):
                print >>outfile, fa.sub('>',line),
                break

        for line in infile:
            print >>outfile, fa.sub('>',line),

    main()
