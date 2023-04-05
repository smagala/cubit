#!/usr/bin/env python
# file clean_fa_comments.py

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
  20060410 JS: Heavy modification to allow use with ncbi as well as
               lanl, better handeling of weird records and a general
               trend away from regex based detection of important parts
               of comment lines
  20080314 JS: Heavy modification to conform to a sane coding style,
               and to leave full strain names instead of generating older-
               style identifiers.
  20090320 JS: Modification to fix NCBI changing the strain name style.
               Commented out LANL code now that it is private.

"""

# standard library imports
import re
import sys
from os.path import splitext
from optparse import OptionParser

# local imports
from core.sequence import parse_fasta, FastaFormatError, line_len, \
                          clean_for_tree, format_fasta


class CommentCleaner(object):
    """Identify and clean comments from various databases"""

    def __init__(self, keep_gi, keep_gb, keep_id, keep_strain, keep_serotype,
            keep_segment):
        # output formatting options
        self.keep_gi = keep_gi
        self.keep_gb = keep_gb
        self.keep_id = keep_id
        self.keep_strain = keep_strain
        self.keep_serotype = keep_serotype
        self.keep_segment = keep_segment

        # ncbi regexs
        self.ncbi_gi_re = re.compile('gi\|(\d+)')
        self.ncbi_gb_re = re.compile('(?:gb|emb|dbj|ref)\|([A-Z0-9_.]+)')
        self.ncbi_strain_re = re.compile('(?:virus|strain) (?:H\d+(?:N\d+)? )?\(?([AB]/.*/\d{2,4})', re.I)
        self.ncbi_serotype_re = re.compile('\(((?:H\d+)?(?:N\d+)?)\)', re.I)
        self.ncbi_type_re = re.compile('Influenza ([ABC]) virus', re.I)
        self.ncbi_segment_re = re.compile(' /\w+/(\d)[\(\s]', re.I)
        self.ncbi_protein_re = re.compile('/([ABFHMNPS12-]{2,6})/')

    def clean_ncbi(self, comment):
        """Clean all comments that come from NCBI"""
        clean_comment = []
        if self.keep_gi:
            try:
                clean_comment.append(self.ncbi_gi_re.search(comment).group(1))
            except AttributeError:
                pass
        if self.keep_gb or self.keep_id:
            try:
                clean_comment.append(self.ncbi_gb_re.search(comment).group(1))
            except AttributeError:
                pass
        if self.keep_strain:
            try:
                clean_comment.append(
                        self.ncbi_strain_re.search(comment).group(1))
            except AttributeError:
                pass
        if self.keep_serotype:
            try:
                clean_comment.append(
                    self.ncbi_serotype_re.search(comment).group(1).upper())
            except AttributeError:
                try:
                    clean_comment.append(
                        self.ncbi_type_re.search(comment).group(1).upper())
                except AttributeError:
                    pass
        if self.keep_segment:
            try:
                clean_comment.append('seg%s' %
                        self.ncbi_segment_re.search(comment).group(1))
            except AttributeError:
                try:
                    clean_comment.append(
                            self.ncbi_protein_re.search(
                            comment).group(1).upper())
                except AttributeError:
                    pass
        return '_'.join(clean_comment)

    def clean_gisaid(self, comment):
        """Clean comment lines from GISAID/EPIFLUDB"""
        parts = comment.split(' | ')
        clean_comment = []
        if self.keep_gb:
            clean_comment.append(parts[4])
        if self.keep_id:
            clean_comment.append(parts[0])
        if self.keep_strain:
            clean_comment.append(parts[2])
        if self.keep_serotype:
            clean_comment.append(parts[5])
        if self.keep_segment:
            clean_comment.append(parts[1])
        return '_'.join(clean_comment)

    def __call__(self, comment):
        dirty_comment = comment
        if self.ncbi_gb_re.search(dirty_comment):
            dirty_comment = self.clean_ncbi(dirty_comment)
        elif len(dirty_comment.split(' | ')) == 6:
            dirty_comment = self.clean_gisaid(dirty_comment)
        clean_comment = clean_for_tree(dirty_comment)
        return clean_comment


def clean_fasta(*args):
    """Clean up the fasta comments from NCBI"""
    # Use a sane output file name
    infile, args = args[0], args[1:]
    root = splitext(infile)[0]
    outfile = root + '_mod.fas'

    # open files for reading, writing
    infile = open(infile, 'r')
    outfile = open(outfile, 'w')

    # loop through each record in the file, clean up the comment, and
    # print the corrected version to the outfile in fasta form
    comment_cleaner = CommentCleaner(*args)
    for comment, seq in parse_fasta(infile, bioedit_cleanup=True):
        comment = comment_cleaner(comment)
        print >> outfile, format_fasta(comment, seq)

    # clean up open files and exit
    infile.close()
    outfile.close()


if __name__ == "__main__":

    def parse_command_line():
        """Get and check command line options and arguments."""

        def _parser_setup():
            """Define all of the parser options"""
            usage = 'usage: %prog [OPTIONS] FILE'
            version = '%prog 0.2'
            description = 'Clean FASTA flu names.  Handles NCBI data'
            parser = OptionParser(usage=usage, version=version, \
                                  description=description)
            parser.add_option('-g', '--keep-gi-number', dest='keep_gi',
                    action='store_true',
                    help='keep gi number on comment line (NCBI)')
            parser.set_defaults(keep_id=False)

            parser.add_option('-b', '--keep-gb-number', dest='keep_gb',
                    action='store_true',
                    help='keep gb number on comment line')

            parser.add_option('-i', '--keep-id-number', dest='keep_id',
                    action='store_true',
                    help='keep accession number on comment line')
            parser.set_defaults(keep_id=False)

            parser.add_option('-s', '--keep-strain', dest='keep_strain',
                    action='store_true',
                    help='keep strain name on comment line')
            parser.set_defaults(keep_strain=False)

            parser.add_option('-t', '--keep-serotype', dest='keep_serotype',
                    action='store_true',
                    help='keep serotype on comment line')
            parser.set_defaults(keep_serotype=False)

            parser.add_option('-n', '--keep-segment', dest='keep_segment',
                    action='store_true',
                    help='keep segment number on comment line')
            parser.set_defaults(keep_segment=False)

            return parser

        parser = _parser_setup()
        opts, args = parser.parse_args()
        try:
            infile = ''.join(args[0].split('"'))
        except IndexError:
            parser.error('Missing input file')
        return (infile, opts.keep_gi, opts.keep_gb, opts.keep_id,
                opts.keep_strain, opts.keep_serotype, opts.keep_segment)

    def main():
        """Avoid cluttering the module namespace"""
        args = parse_command_line()

        try:
            clean_fasta(*args)
        except IOError, msg:
            print msg
            sys.exit(1)
        except FastaFormatError, msg:
            print msg,
            sys.exit(1)

    main()
