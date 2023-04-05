#!/usr/bin/env python
#file common.py

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

"""Common GUI utilities

History:
    20081003 James Smagala: File created

"""

import wx
import os
import os.path

import logging

def pick_save_file(parent=None, msg='Save a file:', default_file='', wild_card=""):
    """Select a location to save a file to"""
    dlg = wx.FileDialog(parent, msg, defaultFile=default_file, wildcard=wild_card,
                        style=wx.SAVE|wx.CHANGE_DIR|wx.OVERWRITE_PROMPT)
    new_file = ''
    if dlg.ShowModal() == wx.ID_OK:
        # this dialog was not canceled
        new_file = dlg.GetPath()
    dlg.Destroy()
    return new_file


def pick_infile(parent=None, msg='Select a file:', wild_card=""):
    """Select a file"""
    dlg = wx.FileDialog(parent, msg, wildcard=wild_card, style=wx.OPEN|wx.CHANGE_DIR)
    new_file = ''
    if dlg.ShowModal() == wx.ID_OK:
        # this dialog was not canceled
        new_file = dlg.GetPath()
    dlg.Destroy()
    return new_file


def error_dialog(parent=None, msg='An unknown error occured'):
    """Display an error dialog"""
    msg = unicode(msg)
    dlg = wx.MessageDialog(parent, msg, caption='Error',
                           style=wx.OK|wx.ICON_ERROR)
    dlg.ShowModal()
    dlg.Destroy()


def save_file(path, contents):
    """Save the contents to a file located at path"""
    try:
        outfile = open(path, 'w')
        print >> outfile, contents
        outfile.close()
        try:
            os.startfile(path)
        except AttributeError:
            pass
        except Exception, errorcode:
            if errorcode[0] == 1155: # this is when OS does not know how to open the file, often not be associated with a type
                pass
            else:
                error_dialog(None, "Failed to open file '%s'" %os.path.basename(path))
    except IOError:
        error_dialog(None, "Failed to write to file '%s'" %\
                           os.path.basename(path))
    finally:
        try:
            outfile.close()
        except NameError:
            pass
