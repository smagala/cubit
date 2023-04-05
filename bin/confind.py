#!/usr/bin/env python
# file confind.py

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

"""Find conserved regions in multiple sequence alignments.

ConFind (CONserved region FINDer) is a tool to detect regions of conservation
in multiple sequence alignments.  Conservation is defined in terms of:
  - the maximum entropy allowed per position
  - the number of exceptions to the maximum entropy allowed
  - the minimum region length
  - the minimum number of sequences that must contain a non-ambiguous character
    for a position to be considered.

Revision History
    20050802 James Smagala: File created - cleanup of 'conserved'
    20050808 JS: parseFasta tested and cleaned up
    20050810 JS: ConservedRegionFinder tested and cleaned up
    20050811 JS: copyright, GPL and contact info added

"""

import os.path
import sys
from math import log # note that log is log base e, or ln
from optparse import OptionParser
from core.sequence import parse_fasta, FastaFormatError, line_len


class ConservedRegionFinderError(Exception):
    """Custom error class for the conserved region finder"""
    pass


class ConservedRegionFinder(object):
    """Locate conserved regions in a multiple sequence alignment"""
    def __init__(self, records, min_region_len, max_entropy, exceptions,
                 num_seqs):
        self.min_region_len = min_region_len
        self.max_entropy = max_entropy
        self.exceptions = exceptions
        self.num_seqs = num_seqs
        self.records = [(comment, seq.upper()) for comment, seq in records]
        self.entropies = self._get_entropies()
        self.conserved_regions = self._get_conserved_regions()

    def _get_entropies(self):
        """Calculate entropy at each position (in Shannon bits)

        Does not include non-ACGT characters in the character count.
        seqs is a list of aligned sequences.  min_seqs is the minimum
        number of sequences needed to calculate an entropy for a given
        position.  If the entropy cannot be calculated, it will be set
        to None.  This parameter should always be greater than 0.
        Returns a list of entropies.

        """

        def get_counts(*chars):
            """Count characters by position

            Takes an list of characters found at a given position of
            a number of sequences, and construct a count dictionary.

            """
            # construct a new count dictionary for this position
            count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

            # loop over the characters, counting only ACGT characters.
            # This ignores gaps and degenericies
            for char in chars:
                # Check for unrecognized characters
                try:
                    count[char] += 1
                except KeyError:
                    if char == None:
                        # This is an invalid alignment
                        raise ConservedRegionFinderError(
                            "Sequences are not all the same length")
            return count

        def calc_entropy(count):
            """Calculate the shannon entropy"""
            total_chars = sum(count)

            def get_contrib(single_char):
                """Calculate the contribution of a single character"""
                # Calculates x*log2(x), where x is the fractional
                # contribution of a single character to the sum of all
                # characters counts in count

                # Fail on divide by zero - this is intentional.  You
                # shouldn't be calculating informational entropy if you
                # have no information.
                fractional_abundance = float(single_char) / total_chars

                # 0 * log2(0) is undefined, but approaches 0 from the right
                # set it to zero if fractional_abundance is 0
                try:
                    return fractional_abundance * \
                           log(fractional_abundance) / log(2)
                except ValueError:
                    return 0

            # Sum over the contributions of all characters in count.
            # abs value instead of multiplying by -1 prevents values
            # like -0.0
            return abs(sum([get_contrib(pos) for pos in count]))


        # For each position, pass a tuple with the character from that
        # position in each sequence to get_counts.
        counts = map(get_counts, *[record[1] for record in self.records])

        # Short circut and/or trick sets the value of the new list to
        # either the calculated shannon entropy, or to None, depending
        # on whether there are at least num_seqs characters in a given
        # column.
        return [(sum(count.values()) >= self.num_seqs and
                [calc_entropy(count.values())] or [None])[0]
                for count in counts]

    def _get_conserved_regions(self):
        """Find all conserved regions that meet the parameters"""
        length = len(self.entropies)

        # Find all possible start sites If this is the first position,
        # or if the position before this one could have been the end
        # point for another conserved region, and this position meets
        # all the requirements, add the position to the list.
        starts = [pos for pos in range(length-self.min_region_len+1)
                  if (pos == 0
                  or self.entropies[pos-1] > self.max_entropy
                  or self.entropies[pos-1] is None)
                  and self.entropies[pos] is not None
                  and self.entropies[pos] < self.max_entropy]

        # Greedy algorithm to get the longest possible regions with no
        # more that the specified maxium number of exceptions to the
        # entropy.
        csvd = []
        skip = None
        for start in starts:
            # if we hit a None, it is the end of the region.
            if skip and start < skip:
                continue
            exceptions = 0
            curr = start
            # Get as many positions as possible
            while curr < length - 1:
                curr += 1
                if self.entropies[curr] > self.max_entropy:
                    exceptions += 1
                # Got as many positions as we could
                if self.entropies[curr] is None:
                    csvd.append((start, curr-1))
                    skip = curr
                    break
                if exceptions > self.exceptions:
                    csvd.append((start, curr-1))
                    break
                # Got as many positions as we could, up to the end of
                # the sequence.  Don't try any more start sites, as they
                # will be redundant.
                if curr == length - 1:
                    csvd.append((start, curr))
                    skip = curr

        # Get only the regions that are long enough
        csvd = [region for region in csvd
                if region[1]-region[0]+1 >= self.min_region_len]
        return csvd

    def get_summary_string(self):
        """Create a string describing the conserved regions"""
        summary = ["Conserved regions\n",
                   "Parameters:",
                   "  Minimum segment length: %s" % self.min_region_len,
                   "  Maximum entropy per position: %s" % self.max_entropy,
                   "  Exceptions allowed: %s" % self.exceptions,
                   "\n%s conserved regions found" %
                   len(self.conserved_regions)]

        # For each position, write the entropy score
        for region in self.conserved_regions:
            summary.append("\n\nPosition %s to %s (%s bases in length)\n" %
                           (region[0] + 1, region[1] + 1,
                            region[1] - region[0] + 1))
            for i in range(region[0], region[1]+1):
                summary.append("Position %5d\t:%.4f" %
                               (i + 1, self.entropies[i]))
        return '\n'.join(summary)

    def get_fasta_strings(self):
        """Create a string with sequence information for each region"""
        for region in self.conserved_regions:
            fasta = []
            for record in self.records:
                fasta.append(">%s: %s to %s\n%s" %
                             (record[0], region[0] + 1, region[1] + 1,
                             line_len(record[1][region[0]:region[1] + 1])))
            yield (region[0] + 1, region[1] + 1, '\n'.join(fasta))


if __name__ == "__main__":

    def parse_command_line():
        """Get and check command line options and arguments."""
        usage = "usage: %prog [OPTIONS] FILE"
        version = "%prog 1.0.01"
        description = "Find conserved regions of sequences in FILE.  FILE " \
                      "should contain aligned sequences in FASTA format."
        parser = OptionParser(usage=usage, version=version,
                              description=description)
        parser.add_option("-l", "--length", dest="min_region_len",
                          metavar="LENGTH", type="int",
                          help="minimum LENGTH of conserved regions "
                               "[default: %default]")
        parser.set_defaults(min_region_len=45)
        parser.add_option("-e", "--entropy", dest="max_entropy",
                          metavar="ENTROPY", type="float",
                          help="maximum bits of Shannon ENTROPY allowed per "
                               "position [default: %default]")
        parser.set_defaults(max_entropy=0.1)
        parser.add_option("-x", "--exceptions", dest="exceptions", type="int",
                          help="number of EXCEPTIONS to the max entropy value "
                               "allowed per conserved region "
                               "[default: %default]")
        parser.set_defaults(exceptions=2)
        parser.add_option("-n", "--min-seqs-exist", dest="num_seqs",
                          type="int",
                          help="minimum NUMBER OF SEQUENCES that must "
                               "contain unambiguous characters to be "
                               "considered for inclusion in a conserved "
                               "region [default: %default]")
        parser.set_defaults(num_seqs=10)
        parser.add_option("-o", "--outfile-prefix", dest="outfile_prefix",
                          help="explicily specify an OUTPUT FILE NAME "
                               "PREFIX - this is a work-around for BioEdit")
        parser.add_option("-f", "--fasta", dest="fasta", action="store_true",
                          help="output FASTA files of each conserved region")
        parser.set_defaults(fasta=False)
        parser.add_option("-s", "--extra-summary", dest="summary",
                          action="store_true",
                          help="output an EXTRA SUMMARY file - "
                               "this is a work-around for BioEdit")
        parser.set_defaults(summary=False)

        options, args = parser.parse_args()
        if options.num_seqs < 1:
            parser.error("The minimum number of sequences (-n) must be >0")
        # Make sure we've got a single input file
        try:
            infile = ''.join(args[0].split('"'))
        except IndexError:
            parser.error("Missing input file")
        # Set outfile_prefix to something sane if the user didn't pass it in.
        if not options.outfile_prefix:
            options.outfile_prefix = os.path.splitext(
                                            os.path.basename(infile))[0]
        return options, infile

    def main():
        """Main program executed if script is run from cmd line"""
        options, file_name = parse_command_line()
        # Attempt to open the file for reading
        try:
            infile = open(file_name, 'r')
        except IOError, msg:
            print msg
            sys.exit(1)
        # If given bad input, this will raise an exception
        try:
            crf = ConservedRegionFinder(parse_fasta(infile),
                                        options.min_region_len,
                                        options.max_entropy,
                                        options.exceptions,
                                        options.num_seqs)
        except FastaFormatError, msg:
            print msg
            sys.exit(1)
        except ConservedRegionFinderError, msg:
            print msg
            sys.exit(1)
        infile.close()

        # Output the summary string to appropriate files in the current
        # directory
        summary = crf.get_summary_string()
        file_name = options.outfile_prefix + "_csvd.txt"
        while True:
            try:
                outfile = open(file_name, 'w')
            except IOError, msg:
                print msg
                sys.exit(1)
            print >> outfile, summary
            outfile.close()
            # Check to see if a constant name file should be output -
            # for use with BioEdit, allowing the summary to pop up in a
            # window
            if not options.summary:
                break
            file_name = "out_csvd.txt"
            options.summary = False

        # If requested, output fasta files of each region found
        if options.fasta:
            for start, end, fasta in crf.get_fasta_strings():
                file_name = options.outfile_prefix + "_csvd_%s_%s.fas" % \
                                                      (start, end)
                try:
                    outfile = open(file_name, 'w')
                except IOError, msg:
                    print msg
                    sys.exit(1)
                print >> outfile, fasta
                outfile.close()

    main()
