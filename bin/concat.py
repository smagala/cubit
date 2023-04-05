#!/usr/bin/env python
# file concat.py

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

"""Concatenate flu genomes for easy analysis

Takes a flu_genome folder, reads the inventory.csv and the concat_files folder,
then generates a concat.fas file with the shortened name and all the genes
concatenated.

Revision History
    20060412 James Smagala: rewrite of flu_genome

"""

from core.sequence import parse_fasta, FastaFormatError
from optparse import OptionParser
from sys import exit
import os

def lineLen(line,length=79):
  """Format long lines to a given length"""
  lines = []
  while line:
    lines.append(line[:length])
    line = line[length:]
  return '\n'.join(lines)

def parseCommandLine():
  """Get and check command line options and arguments."""
  usage = "usage: %prog DATA_DIR"
  version = "%prog 0.1"
  description = "Read a flu_genome data directory, concatenate the\
 concat_files data according to the inventory.csv"
  parser = OptionParser(usage=usage,version=version,description=description)

  options, args = parser.parse_args()

  try:
    data_dir = ''.join(args[0].split('"'))
  except IndexError:
    print "No DATA_DIR argument on command line"
    exit(1)

  return data_dir

def concat(data_dir):
  os.chdir(os.path.join(data_dir,"results"))
  infile = 'inventory.csv'

  infile = open(infile,'r')
  outfile = open("concat.fas",'w')

  os.chdir("concat_files")

  # read in all the fasta records available in the directory
  records = []
  fasta_files = [f for f in os.listdir(os.getcwd()) if f.endswith(".fas")]
  for fafile in fasta_files:
    ff = open(fafile)
    records.extend([(comment,seq) for comment,seq in parse_fasta(ff)])
    ff.close()
  records = dict(records)

  seq = ''
  for line in infile:
    if line.startswith("full name"):
      continue
    if seq:
      print >>outfile, ">%s\n%s" % (strain,lineLen(seq))
      seq = ''
    line = line.strip().split(',')
    for c in range(10):
      # column 0 is the name of the virus, column 1 the short name
      if c == 0:
        continue
      if c == 1:
        strain = line[c]
        continue
      # columns 2-9 are gene segments 1-8 in order (see column:gene map above)
      acc = line[c].strip().split()[0]
      if not acc:
        print "Not all genes present for %s" % strain
        exit(1)
      seq = seq + records[acc]

  print >>outfile, ">%s\n%s" % (strain,lineLen(seq))

  # clean up all open files 
  infile.close()
  outfile.close()


if __name__ == "__main__":
  # we should only check command line args if we were actually called from
  # the command line
  data_dir = parseCommandLine()

  # handle all errors here
  try:
    concat(data_dir)
  except IOError,msg:
    print msg
    exit(1)
  except OSError,msg:
    print msg
    exit(1)
  except FastaFormatError,msg:
    print msg
    exit(1)
