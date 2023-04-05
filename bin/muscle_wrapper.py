#!/usr/bin/env python
# file muscle_wrapper.py

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

"""Clean up i/o from BioEdit

Revision History
    20080606 James Smagala: File created
    20130319 JAS: Fixed assumed directory to check current path instead

"""

from optparse import OptionParser
import os
import sys
from core.sequence import parse_fasta


if __name__ == "__main__":

    def parse_cmd_line():
        """Get and check command line options and arguments."""
        usage = "usage: %prog [OPTIONS] FILE"
        version = "%prog 0.1"
        description = "Align seqs using muscle, avoiding BioEdit stupidity"
        parser = OptionParser(usage=usage,version=version,description=description)
        parser.add_option("-i", "--first-iteration-only", dest="iteration",
                          action="store_true",
                          help="first alignment step only")
        parser.set_defaults(align=False)
        parser.add_option("-f", "--fast", dest="fast",
                          action="store_true",
                          help="optimize for speed")
        parser.set_defaults(fast=False)
        options, args = parser.parse_args()

        try:
            infile = ''.join(args[0].split('"'))
        except IndexError:
            parser.error("Missing input file")
        
        optstr = '-maxmb 1500'
        if options.fast:
            optstr = '-diags ' + optstr
        if options.iteration:
            optstr = '-maxiters 1 ' + optstr

        return infile, optstr

    def main():
        # directory issue reported by Rubin Donis in CITGO - don't
        # assume a directory, and calculate it differently depending on
        # if we are running as a script or executable
        cubit_dir = os.path.dirname(sys.executable if hasattr(sys, 'frozen')
                else sys.argv[0])
        muscle_dir = os.path.join(os.path.split(cubit_dir)[0], 'muscle')
        muscle_path = os.path.join(muscle_dir, 'muscle.exe')

        infile, optstr = parse_cmd_line()

        # Attempt to open the file for reading
        try:
            infile = open(infile, 'r')
            tmpfile = open('in_clean.fas', 'w')
            outfile = open('out.fas', 'w')
        except IOError,msg:
            print msg
            sys.exit(1)

        for comment, seq in parse_fasta(infile):
            comment = comment.strip()
            comment = ' '.join(comment.split()[:-2])
            print >>tmpfile, ">%s\n%s" % (comment, seq)

        infile.close()
        tmpfile.close()

        ifile = 'in_clean.fas'
        ofile = 'out.fas'

        tmpfile = os.popen('%s %s -in %s -out %s' % 
                (muscle_path, optstr, ifile, ofile))
        for line in tmpfile:
            print line,
        tmpfile.close()
      
        os.remove('in_clean.fas')

    main()
