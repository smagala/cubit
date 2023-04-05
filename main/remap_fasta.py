#!/usr/bin/env python
#file remap_fasta.py

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

"""Remap the comment lines of fasta files using an xls key

History:
    20081028 James Smagala: File created

"""

import wx
import os.path
from core.sequence import read_fasta_file, format_fasta
from gui.remap_fasta_gui import GUIApp, GUIMainFrame
from gui.common import pick_save_file, pick_infile, error_dialog, save_file
from core.xlsreader import parse_xls


class RemapFastaMainFrame(GUIMainFrame):
    """The wx main frame for remap_fasta - implements the event handlers"""
    def __init__(self, *args, **kwargs):
        # call this first so that the underlying C++ object exists
        # call the parent __init__ to handle additional setup
        super(RemapFastaMainFrame, self).__init__(*args, **kwargs)

        # get a reference to the app object
        self.app = wx.GetApp()

        # bind events to handlers
        self.Bind(wx.EVT_CLOSE, self.on_exit)

        # set some other variables
        self.base_file_name = ''
        self.seqs = []
        self.map_dict = []

        # set an icon if desired
        icon = wx.Icon('icons/remap_fasta.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)

    def _set_button_state(self):
        """Enable the analysis buttons once a file is selected"""
        self.run_button.Enable(bool(self.seqs and self.map_dict))

    def _remap_fasta(self):
        """Generate a new fasta file, using comments found in the xls key"""
        new_seqs = []
        max_row = max([cell[0] for cell in self.map_dict])
        max_col = max([cell[1] for cell in self.map_dict])

        for comment, seq in self.seqs:
            new_comment = ''
            value = comment.split()[0]
            value = comment.split('_')[0]
            for row in range(max_row + 1):
                # try to find the identifier in the xls map
                # if we found it, use it to construct a new comment
                if (row, 0) in self.map_dict:
                    row_key = self.map_dict[(row, 0)]
                    try:
                        int_row_key = int(row_key)
                    except ValueError:
                        int_row_key = None
                    if value == str(row_key) or value == int_row_key:
                        bits = []
                        for col in range(1, max_col + 1):
                            if (row, col) in self.map_dict:
                                value = self.map_dict[(row, col)]
                                try:
                                    bits.append(str(int(value)))
                                except ValueError:
                                    bits.append(str(value))
                        new_comment = '_'.join(bits)
                else:
                    continue
            if new_comment:
                comment = new_comment
            new_seqs.append((comment, seq))

        return '\n'.join([format_fasta(comment, seq)
                          for comment, seq in new_seqs])

    def on_get_fasta(self, event):
        """Read in a fasta file"""
        path = pick_infile(parent=self, msg='Pick a FASTA File:')

        if path:
            try:
                self.seqs = read_fasta_file(path)
                self.base_file_name = os.path.splitext(path)[0]
                self.fasta_box.SetValue(path)
                self._set_button_state()
            except Exception, msg:
                error_dialog(self, msg=msg)
        event.Skip()

    def on_get_xls(self, event):
        """Read in a .xls map file"""
        path = pick_infile(parent=self, msg='Pick an xls Map File:')

        if path:
            try:
                # this is a dict with entries like {(x, y): value}
                self.map_dict = parse_xls(path)[0][1]
                self.xls_box.SetValue(path)
                self._set_button_state()
            except Exception, msg:
                error_dialog(self, msg=msg)
        event.Skip()

    def on_remap(self, event):
        """Generate a new fasta file with remapped comment lines"""
        path = pick_save_file(self, msg='Save Remapped FASTA File:',
                              default_file=self.base_file_name + \
                                           '_remapped.fas')
        if path:
            save_file(path, self._remap_fasta())
        event.Skip()

    def on_exit(self, event):
        """Explicitly destroy the main frame, ending the program"""
        self.Destroy()
        self.app.Exit()
        event.Skip()


class RemapFastaApp(GUIApp):
    """The wx GUI app class for remap_fasta"""
    def OnInit(self):
        """Required setup for the wx app"""
        wx.InitAllImageHandlers()

        # set a top level frame and show it
        main_frame = RemapFastaMainFrame(None, -1, '')
        self.SetTopWindow(main_frame)
        main_frame.Show()

        # must return true or the app will exit
        return True


if __name__ == "__main__":

    def init_gui():
        """GUI setup and initalization"""
        app = RemapFastaApp(redirect=False)
        app.MainLoop()

    init_gui()
