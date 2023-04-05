#!/usr/bin/env python
# file fasta_tree_order.py

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

"""Parse .dnd files, get accession numbers, and retrieve fasta records

Assumes that all sequences of interest are in a single large fasta file

Revision History
    20050418 James Smagala: Rework of parse_fa
    20050818 JS: Cleaned up, use optparse, fasta parser copied from confind.py
    20051108 JS: Fix to work with longer names
    20081008 JS: Added a GUI

"""

import re
import os.path

import wx

from core.sequence import line_len, read_fasta_file
from gui.fasta_tree_order_gui import GUIApp, GUIMainFrame
from gui.common import pick_save_file, pick_infile, error_dialog, save_file


class FastaTreeOrderMainFrame(GUIMainFrame):
    """wx main frame for fasta_tree_order, implements event handlers"""
    def __init__(self, *args, **kwargs):
        # call this first so that the underlying C++ object exists
        # call the parent __init__ to handle additional setup
        super(FastaTreeOrderMainFrame, self).__init__(*args, **kwargs)

        # get a reference to the app object
        self.app = wx.GetApp()

        # set some other variables
        self.base_file_name = ''
        self.seqs = ''
        self.labels = ''
        self.tree_re = re.compile(r',[^)(:,;]*?:|\([^)(:,;]*?:')

        # bind events to handlers
        self.Bind(wx.EVT_CLOSE, self.on_exit)

        # set an icon if desired
        icon = wx.Icon('icons/fasta_tree_order.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

    def _check_button_state(self):
        """Enable the analysis buttons once a file is selected"""
        if self.seqs and self.labels:
            self.new_fasta_button.Enable(True)
        else:
            self.new_fasta_button.Enable(False)

    def on_select_fasta(self, event):
        """Read in a fasta file"""
        path = pick_infile(parent=self, msg='Pick a FASTA file:')
        if path:
            try:
                self.seqs = dict(read_fasta_file(path, bioedit_cleanup=True))
                self.base_file_name = os.path.splitext(path)[0]
                self.fasta_box.SetValue(path)
                self._check_button_state()
            except Exception, msg:
                error_dialog(self, msg=msg)
        event.Skip()

    def on_select_tree(self, event):
        """Read in a fasta file"""
        path = pick_infile(parent=self, msg='Pick a tree file:')
        if path:
            try:
                infile = open(path, 'r')
                tree = ''.join([line.strip() for line in infile])
                infile.close()
                self.labels = [label[1:-1] for label in
                               self.tree_re.findall(tree)]
                self.tree_box.SetValue(path)
                self._check_button_state()
            except Exception, msg:
                error_dialog(self, msg=msg)
        event.Skip()

    def on_generate_fasta(self, event):
        """Generate a new fasta file from the old fasta and tree labels"""
        result = []
        for label in self.labels:
            try:
                result.append(">%s\n%s" % (label, line_len(self.seqs[label])))
            except KeyError:
                #print "Label '%s' is not in database" % label
                pass
        result = '\n'.join(result)
        path = pick_save_file(self, msg='Save new FASTA file:',
                              default_file=self.base_file_name +\
                                           '_dnd.fas')
        if path:
            save_file(path, result)
        event.Skip()

    def on_exit(self, event):
        """Explicitly destroy the main frame, ending the program"""
        self.Destroy()
        self.app.Exit()
        event.Skip()


class FastaTreeOrderApp(GUIApp):
    """The wx GUI app class for aadiff"""
    def OnInit(self):
        """Required setup for the wx app"""
        wx.InitAllImageHandlers()

        # set a top level frame and show it
        main_frame = FastaTreeOrderMainFrame(None, -1, '')
        self.SetTopWindow(main_frame)
        main_frame.Show()

        # must return true or the app will exit
        return True


if __name__ == "__main__":

    def init_gui():
        """GUI setup and initalization"""
        app = FastaTreeOrderApp(redirect=False)
        app.MainLoop()

    init_gui()
