#!/usr/bin/env python
# file count_gs.py

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

"""Count glycosylation sites on HA

Revision History
    20060727 James Smagala: file creation
    20080314 JS: Major refactor to attempt to remove cogent entirely.  Also
                 reformatted to match typical python code style.
"""

import wx
import os
import os.path
from core.sequence import read_fasta_file, DNASequence, PROTEIN
from gui.count_gs_gui import GUIApp, GUIMainFrame
from gui.common import pick_save_file, pick_infile, error_dialog, save_file


class CountGSMainFrame(GUIMainFrame):
    """wx main frame for count_gs, implements the event handlers"""
    def __init__(self, *args, **kwargs):
        # call this first so that the underlying C++ object exists
        # call the parent __init__ to handle additional setup
        super(CountGSMainFrame, self).__init__(*args, **kwargs)

        # get a reference to the app object
        self.app = wx.GetApp()

        # set some other variables
        self.base_file_name = ''
        self.fasta = ''
        self.table = ''
        self.pos2 = set(PROTEIN) - set('P')
        self.pos3 = set('ST')

        # bind events to handlers
        self.Bind(wx.EVT_CLOSE, self.on_exit)

        # set an icon if desired
        icon = wx.Icon('icons/count_gs.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

    def _set_button_state(self, enabled):
        """Enable the analysis buttons once a file is selected"""
        self.table_button.Enable(bool(enabled))
        self.fasta_button.Enable(bool(enabled))

    def count_gs(self, seqs):
        """Count glycosylation sites by translating DNA"""
        # make a list of glycosylation site positions
        sites = []
        for seq in seqs:
            sites.append([])
            aa_seq = seq.translate_disambiguate()
            for start in range(0, len(aa_seq) - 2):
                if not 'N' in aa_seq[start]:
                    continue
                if not aa_seq[start+1].intersection(self.pos2):
                    continue
                if not aa_seq[start+2].intersection(self.pos3):
                    continue
                sites[-1].append(start)
            seq.name = seq.name + '_%dgs' % len(sites[-1])

        # compile a list of all sites in all sequences
        all_sites = {}
        for per_seq_sites in sites:
            for site in per_seq_sites:
                all_sites[site] = None

        # sort the list
        all_sites = all_sites.keys()
        all_sites.sort()

        # make strings of the fasta records and the results table
        table = 'Strain,# glycosylation sites,%s' % \
                ','.join([str(pos + 1) for pos in all_sites])
        for i in range(len(seqs)):
            tmp = []
            for site in all_sites:
                if site in sites[i]:
                    tmp.append('X')
                else:
                    tmp.append('')
            tmp = ','.join(tmp)
            table += "\n%s,%s,%s" % (seqs[i].name, len(sites[i]), tmp)
        self.fasta = '\n'.join([seq.to_fasta() for seq in seqs])
        self.table = table

    def on_select_alignment(self, event):
        """Read in a fasta file containing an alignment and count gs sites"""
        path = pick_infile(parent=self, msg='Pick a nt alignment:')
        if path:
            # try to read in the alignment and diff it
            # only the results are actually stored
            try:
                progress = wx.ProgressDialog(title='Progress',
                    message='',
                    maximum=2,
                    parent=self,
                    style=wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT|wx.PD_REMAINING_TIME)
                progress.Update(0, 'Reading FASTA records')
                seqs = [DNASequence(comment, seq) for comment, seq in
                        read_fasta_file(path, bioedit_cleanup=True)]
                progress.Update(1, 'Counting glycosylation sites')
                self.count_gs(seqs)
                progress.Update(2, 'Finished')
            except Exception, msg:
                error_dialog(self, msg=msg)
            finally:
                progress.Destroy()
            self.base_file_name = os.path.splitext(path)[0]
            self.file_box.SetValue(path)
            self._set_button_state(enabled=path)
        event.Skip()

    def on_generate_table(self, event):
        """Generate a glycosylation table"""
        path = pick_save_file(self, msg='Save glycosylation table:',
                              default_file=self.base_file_name + \
                                           '_gs_summary.csv')
        if path:
            save_file(path, self.table)
        event.Skip()

    def on_generate_fasta(self, event):
        """Generate a glycosylation table"""
        path = pick_save_file(self, msg='Save labeled FASTA:',
                              default_file=self.base_file_name +\
                                           '_gs.fas')
        if path:
            save_file(path, self.fasta)
        event.Skip()

    def on_exit(self, event):
        """Explicitly destroy the main frame, ending the program"""
        self.Destroy()
        self.app.Exit()
        event.Skip()


class CountGSApp(GUIApp):
    """The wx GUI app class for CountGS"""
    def OnInit(self):
        """Required setup for the wx app"""
        wx.InitAllImageHandlers()

        # set a top level frame and show it
        main_frame = CountGSMainFrame(None, -1, '')
        self.SetTopWindow(main_frame)
        main_frame.Show()

        # must return true or the app will exit
        return True


if __name__ == "__main__":

    def init_gui():
        """GUI setup and initalization"""
        app = CountGSApp(redirect=False)
        app.MainLoop()

    init_gui()
