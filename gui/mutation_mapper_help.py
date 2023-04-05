'''
#!/usr/bin/env python
#file mutation_mapper_help.py

# Copyright (C) 2008 James Smagala/Melih Gunay

# Contact:
#   James Smagala, Melih Gunay
#   smagala@gmail.com, gmelih@gmail.com

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
    2008.02.06 James Smagala: File created
    2008.10.03 JS: Added wx GUI
    2009.06.29 JS: Modified to allow transposed table
    2011.08.26 Melih Gunay: Modified to add tree order and mutation mapping functionality
    2011.09.15 Melih Gunay: Modified to add phyloxml support for output tree
    2011.12.25 Melih Gunay: Newick output format is changed to support 
                            mutations on branches as oposed to at subtree nodes
"""

'''
import wx.html
import os, sys

class GUIHelp():

    def __init__(self, parent, id, title, logger=None):
        hlpfile = 'docs'+os.sep+'mutation_mapper.hlp'
        f = get_file(hlpfile)
        if f is None:
            logger.debug("Help logger file not found in path...")
            return
        dialog = wx.Dialog(parent, id, title, size=(800, 600))
        content = ""
        for line in f:
            line = line.replace("\\\\",  os.sep)
            content += str(line + "\n")
        f.close()
        html = wx.html.HtmlWindow(dialog)
        html.SetPage(content)
        hbox = wx.BoxSizer(wx.VERTICAL)
        hbox.Add(html, 1, wx.EXPAND | wx.ALL, 20)
        dialog.SetSizer(hbox)
        dialog.ShowModal()
        dialog.Destroy()


def get_file(fname):
    '''Try to agresively locate file, in unix and windows'''
    path = os.getcwd() + os.sep + fname
    exists = os.access(path, os.R_OK)
    if not exists:
        '''Ugly -- check if it is bundled with windows executable..'''
        path = os.path.dirname(sys.executable) + os.sep + fname
    exists = os.access(path, os.R_OK)
    if not exists:
        return None
    try:
        return open(path, 'r')
    except IOError:
        return None        
