#!/usr/bin/env python
# file mark_resistant.py

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

"""Mark records that have antiviral resistance

Revision History
    20060725 James Smagala: file creation

"""

from sys import exit
from optparse import OptionParser
from core.sequence import parse_fasta, FastaFormatError, DNASequence
import os.path
import re

def parse_command_line():
    """Get and check command line options and arguments."""
    usage = "usage: %prog GENE INFILE [GENE must be M, N1 or N2]"
    version = "%prog 0.1"
    description = "Generates an alignment with copies of all sequences in INFILE " \
                  "with only the codons in CODON_LIST"
    parser = OptionParser(usage=usage,version=version,description=description)

    options, args = parser.parse_args()

    gene = args[0]
    if not gene in ['M','N1','N2']:
        parser.error('GENE must be M, N1 or N2')
        exit(1)

    try:
        infile = ''.join(args[1].split('"'))
    except IndexError:
        parser.error("Missing input")
        exit(1)

    return infile, gene

def mark_res(infile, gene):
    outfile = os.path.splitext(os.path.basename(infile))[0] + "_marked.fas"

    # mutations known to cause antiviral resistance
    genes = {'M':{26:'F',
                  27:'TA',
                  30:'TV',
                  31:'DN',
                  34:'E'},
             'N1':{119:'V',
                   152:'K',
                   275:'Y',
                   293:'K'},
             'N2':{119:'V',
                   152:'K',
                   274:'Y',
                   292:'K'}}

    # select the mutation list for the appropriate gene
    mutations = genes[gene]

    # read in the sequence data and clean it up
    seqs = []
    # check if bioedit added its tag to the end of the comment line
    bioedit = re.compile(".*\s+\d+\sbases\Z")
    infile = open(infile,'r')
    for comment, seq in parse_fasta(infile):
        seq = ''.join(seq.upper().split())
        if bioedit.search(comment):
            comment = ''.join(comment.split()[:-2])
        seqs.append((comment,seq))
    infile.close()

    # mark the records and write them out
    outfile = open(outfile,'w')
    for comment, seq in seqs:
        res = False
        unk = False
        for codon in mutations.keys():
            # get the start of the nt sequence
            start = (codon-1)*3
            try:
                # get the codon and translate it
                aa = str(DNASequence('', seq[start:start+3]).translate())
            except IndexError:
                print "%s is missing data at codon %s" % (comment,codon)
                continue
            # if there is a resistant mutation, mark it
            if aa in mutations[codon]:
                res = True
                comment = comment + "_%s%s" % (codon,aa)
            # if there is an unknown position, mark it and look at other positions
            if aa == 'X':
                unk = True
        if res:
            comment = comment + "_res"
        elif unk:
            comment = comment + "_unk"
        else:
            comment = comment + "_sen"
        print >>outfile, ">%s\n%s" % (comment,seq)

    outfile.close()

if __name__ == "__main__":
    infile, gene = parse_command_line()

    try:
        mark_res(infile,gene)
    except IOError,msg:
        print msg
    except FastaFormatError,msg:
        print msg
