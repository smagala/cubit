#!/usr/bin/env python
#file aadiff.py

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

"""Construct an amino acid difference table with aligned coding nt seqs.

The first sequence in the file will be treated as a reference consensus
sequence, and all other sequences will be compared against it.

Degenerate codons will be disambiguated to detect all changes.

History:
    20080206 James Smagala: File created
    20081003 JS: Added wx GUI
    20090629 JS: Modified to allow transposed table

"""

import wx
import os
import os.path
from core.sequence import read_fasta_file, DNASequence, diff2seqs
from gui.aadiff_gui import GUIApp, GUIMainFrame
from gui.common import pick_save_file, pick_infile, error_dialog, save_file


def aadiff(reference, seqs):
    """Construct an amino acid diff table from aligned coding nucleotides"""
    # translate the reference sequence
    reference_aa = reference.translate_disambiguate()

    # translate each remaining sequence and compare it to the reference
    # note - translate_disambiguate disambiguates codons, then
    # translates those.  This ensures that there can never be more than
    # 64 things to translate.  This makes it *much* faster than trying
    # to disambiguate the entire sequence, then translate all possible
    # resulting sequneces, which can grow exponentially
    seqs_aa = []
    changesets = []
    for seq in seqs:
        seq_aa = seq.translate_disambiguate()
        seqs_aa.append(seq_aa)
        changeset = diff2seqs(reference_aa, seq_aa)
        changesets.append(changeset)

    return reference_aa, seqs_aa, changesets


class AADiffMainFrame(GUIMainFrame):
    """The wx main frame for aadiff - implements the event handlers"""

    def __init__(self, *args, **kwargs):
        # call this first so that the underlying C++ object exists
        # call the parent __init__ to handle additional setup
        super(AADiffMainFrame, self).__init__(*args, **kwargs)

        # get a reference to the app object
        self.app = wx.GetApp()

        # store the base file name to derive other file names from
        self.base_file_name = ''

        # bind events to handlers
        self.Bind(wx.EVT_CLOSE, self.on_exit)

        # set an icon if desired
        icon = wx.Icon('icons/aadiff.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

    def _fix_gaps(self):
        """correct for deletions 'del' that are noted as missing data '-'"""
        # add back in '-' chars for only positions that already have a
        # reported change
        self.positions = []
        for changeset in self.changesets:
            self.positions.extend(changeset.keys())
        self.positions = dict([(pos, None) for pos in self.positions]).keys()
        self.positions.sort()
        for seq_aa, changes in zip(self.seqs_aa, self.changesets):
            for pos in self.positions:
                if pos in changes:
                    if '-' in changes[pos]:
                        changes[pos].remove('-')
                        changes[pos].add('del')
                else:
                    aa_set = seq_aa[pos]
                    if '-' in aa_set:
                        changes[pos] = ['-']

    def _calc_aa_strings(self):
        """Pre-calculate the reference amino acid strings"""
        # we only need strings for positions we know we will use
        self.reference_str = []
        for pos, aa_set in enumerate(self.reference_aa):
            if pos in self.positions:
                aa_set = list(aa_set)
                aa_set.sort()
                self.reference_str.append('/'.join(aa_set))
            else:
                self.reference_str.append('')

    def _gen_change_report(self):
        """Show all positions that differ across all sequences"""
        lines = ["Isolate,Number of Differences,Positions Different",
                 "%s,," % self.reference.name]
        for seq, seq_aa, changeset in zip(
                self.seqs, self.seqs_aa, self.changesets):
            positions = changeset.keys()
            positions.sort()
            changes = []
            for pos in positions:
                ref = self.reference_str[pos]
                aa_set = list(seq_aa[pos])
                aa_set.sort()
                aa_str = '/'.join(aa_set)
                # do not include trailing gaps a changes in the report
                if aa_str == '-':
                    continue
                if self.report_aa_option.IsChecked():
                    changes.append('%s%s%s' % (ref, str(pos + 1), aa_str))
                else:
                    changes.append(str(pos + 1))
            if changes:
                csv_positions = '"%s"' % ', '.join(changes)
            else:
                csv_positions = ''
            lines.append('%s,%d,%s' % (seq.name, len(changes), csv_positions))
        return '\n'.join(lines)

    def _gen_diff_table(self):
        """List positions that differ for each virus"""
        # build a list of lists so it can be reversed later if needed
        lines = []
        header = ['', self.reference.name, ]
        header.extend([seq.name for seq in self.seqs])
        lines.append(header)
        for pos in self.positions:
            line = []
            line.append(str(pos + 1))
            line.append(self.reference_str[pos])
            for changeset in self.changesets:
                if not pos in changeset:
                    line.append('')
                else:
                    changes = changeset[pos]
                    changes = list(changes)
                    changes.sort()
                    changes = '"%s"' % '/'.join(changes)
                    line.append(changes)
            lines.append(line)
        if self.transpose_table_option.IsChecked():
            lines = map(None, *lines)
            #lines = [[pos if pos is not None else '' for pos in line] 
            #        for line in lines]
        return '\n'.join([','.join(line) for line in lines])

    def _set_button_state(self, enabled):
        """Enable the analysis buttons once a file is selected"""
        self.table_button.Enable(bool(enabled))
        self.report_button.Enable(bool(enabled))

    def on_select_alignment(self, event):
        """Read in a fasta file containing an alignment and diff it"""
        path = pick_infile(parent=self, msg='Pick a nt alignment:')
        updated = False
        if path:
            # try to read in the alignment and diff it
            # only the results are actually stored
            try:
                updated = False
                self.file_box.SetValue('')
                self._set_button_state(enabled=updated)
                seqs = [DNASequence(comment, seq) for comment, seq in
                             read_fasta_file(path)]
                self.reference = seqs[0]
                self.seqs = seqs[1:]
                self.reference_aa, self.seqs_aa, self.changesets = \
                        aadiff(self.reference, self.seqs)
                self._fix_gaps()
                self._calc_aa_strings()
                updated = True
            except Exception, msg:
                error_dialog(self, msg=msg)
            if updated:
                self.base_file_name = os.path.splitext(path)[0]
                self.file_box.SetValue(path)
                self._set_button_state(enabled=updated)
        event.Skip()

    def on_generate_table(self, event):
        """Generate a difference table"""
        table = self._gen_diff_table()
        path = pick_save_file(self, msg='Save difference table:',
                default_file=self.base_file_name + '_aadiff_table.csv')
        if path:
            save_file(path, table)
        event.Skip()

    def on_generate_report(self, event):
        """Generate a change report"""
        report = self._gen_change_report()
        path = pick_save_file(self, msg='Save change report:',
                default_file=self.base_file_name + '_aadiff_report.csv')
        if path:
            save_file(path, report)
        event.Skip()

    def on_exit(self, event):
        """Explicitly destroy the main frame, ending the program"""
        self.Destroy()
        self.app.Exit()
        event.Skip()


class AADiffApp(GUIApp):
    """The wx GUI app class for aadiff"""
    def OnInit(self):
        """Required setup for the wx app"""
        wx.InitAllImageHandlers()

        # set a top level frame and show it
        main_frame = AADiffMainFrame(None, -1, '')
        self.SetTopWindow(main_frame)
        main_frame.Show()

        # must return true or the app will exit
        return True


if __name__ == "__main__":

    def init_gui():
        """GUI setup and initalization"""
        app = AADiffApp(redirect=False)
        app.MainLoop()

    init_gui()
