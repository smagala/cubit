#!/usr/bin/env python
# file label_tree.py

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

"""Parse .dnd/.nwk files and label internal nodes.

Takes an unlabled Newick tree as input, and outputs a tree with labeled nodes.

This program is not robust - it is just a quick hack to get labels on a tree
to simplify further viewing and manipulation.

Revision History
  050420 James Smagala: Rework of dnd2fa
  050817 JS: Cleaned up, added OptionParser

"""

from optparse import OptionParser
from sys import exit

def parseCommandLine():
  """Get and check command line options and arguments."""
  usage = "usage: %prog [OPTIONS] FILE"
  version = "%prog 0.2"
  description = "Add unique node labels to a .dnd/.nwk FILE."
  parser = OptionParser(usage=usage,version=version,description=description)

  parser.add_option("-o", "--output", dest="outfile",
                    help="output file location [default: %default]")
  parser.set_defaults(outfile="out.dnd")

  options, args = parser.parse_args()

  try:
    infile = ''.join(args[0].split('"'))
  except IndexError:
    parser.error("Missing input file")
  outfile = ''.join(options.outfile.split('"'))

  return infile,outfile


def main():
  # read the input and possibly output files from the command line
  infile,outfile = parseCommandLine()

  try:
    infile = open(infile,'r')
    outfile = open(outfile,'w')
  except IOError, msg:
    print msg
    exit (1)

  # Read in a tree
  tree = ''.join((''.join([line for line in infile])).split())
  infile.close()

  # Split on )
  split_tree = tree.split(')')

  # Add ')' back in, add node number after that
  new_tree = []
  count = 0
  for node in split_tree:
    new_tree.append(node)
    new_tree.append(')')
    new_tree.append(repr(count))
    count = count + 1
  tree = ''.join(new_tree)

  print >>outfile, tree
  outfile.close()

if __name__ == '__main__':
  main()
