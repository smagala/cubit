#!/usr/bin/env python
# file swap_plates.py

"""A GUI to manage the names of files from two mixed-up sequencing plates.

History:
    20080813 James Smagala: Project created, initial gui design and
                            functionality

"""

import os
import os.path
import shutil

import wx
from gui.swap_plates_gui import SwapPlateGUIMainFrame


class SwapPlateError(Exception):
    """Custom error class for swap_plates specific errors"""
    pass


def pick_dir(parent=None, current_dir=''):
    """Select a directory"""
    dlg = wx.DirDialog(parent, "Select a directory:", defaultPath=current_dir)
    dir = current_dir
    if dlg.ShowModal() == wx.ID_OK:
        # this dialog was not canceled
        dir = dlg.GetPath()
    dlg.Destroy()
    return dir


def error_dialog(parent=None, message='An unknown error occured'):
    """Display an error dialog"""
    message = unicode(message)
    dlg = wx.MessageDialog(parent, message, caption='Error',
                           style=wx.OK|wx.ICON_ERROR)
    dlg.ShowModal()
    #dlg.Destroy()


class SwapPlateMainFrame(SwapPlateGUIMainFrame):
    """Inherit from wxGlade-generated frame and add functionality"""
    def __init__(self, *args, **kwargs):
        super(SwapPlateMainFrame, self).__init__(*args, **kwargs)

        # setup an icon for the window
        icon = wx.Icon('icons/swap_plates.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon) 

        # get a reference to the app object
        self.app = wx.GetApp()

        # bind events to handlers
        self.Bind(wx.EVT_CLOSE, self.OnExit)

    def OnGetDir1(self, event):
        """Get the first dir using a directory picker"""
        current_dir = self.dir1.GetValue()
        if not current_dir:
            current_dir = self.dir2.GetValue()
            if current_dir:
                head, tail = os.path.split(current_dir)
                current_dir = head
        new_dir = pick_dir(self, current_dir)
        self.dir1.SetValue(new_dir)
        event.Skip()

    def OnGetDir2(self, event):
        """Get the second dir using a directory picker"""
        current_dir = self.dir2.GetValue()
        if not current_dir:
            current_dir = self.dir1.GetValue()
            if current_dir:
                head, tail = os.path.split(current_dir)
                current_dir = head
        new_dir = pick_dir(self, current_dir)
        self.dir2.SetValue(new_dir)
        event.Skip()

    def OnRunSwap(self, event):
        """Actually swap files"""
        dir1 = self.dir1.GetValue()
        dir2 = self.dir2.GetValue()
        try:
            # check for problems with user-specified directories
            if not os.path.isdir(dir1):
                raise SwapPlateError('Pick a valid directory for plate 1')
            if not os.path.isdir(dir2):
                raise SwapPlateError('Pick a valid directory for plate 2')
            if dir1 == dir2:
                raise SwapPlateError("Pick two different directories")
            # check for problems with directories to be created
            dir1_fix = dir2 + '_swapped'
            dir2_fix = dir1 + '_swapped'
            if os.path.exists(dir1_fix):
                raise SwapPlateError("'%s' exists" % dir1_fix)
            if os.path.exists(dir2_fix):
                raise SwapPlateError("'%s' exists" % dir2_fix)
            # create the target directories
            os.mkdir(dir1_fix)
            os.mkdir(dir2_fix)
            # get a list of files to be swapped
            plate1 = [(f.split('.')[-3], f) for f in os.listdir(dir1)
                      if os.path.isfile(os.path.join(dir1, f))
                      and f.endswith('.ab1')]
            plate2 = [(f.split('.')[-3], f) for f in os.listdir(dir2)
                      if os.path.isfile(os.path.join(dir2, f))
                      and f.endswith('.ab1')]
            plate1.sort()
            plate2.sort()
            # check if the plate layouts match
            if not all([well1 == well2 for (well1, f1), (well2, f2) in
                        zip(plate1, plate2)]):
                raise SwapPlateError("The plate layouts don't match")
            # actually perform the swap
            progress = wx.ProgressDialog(title='Progress',
                message='Copying files', maximum=len(plate1), parent=self,
                style=wx.PD_AUTO_HIDE|wx.PD_CAN_ABORT|wx.PD_REMAINING_TIME)
            for n, ((well1, f1), (well2, f2)) in enumerate(zip(plate1, plate2)):
                # note - this is a change from the original version
                # before, the files swapped
                # now, the files stay the same but the name swaps
                source1 = os.path.join(dir1, f1)
                dest1 = os.path.join(dir1_fix, f2)
                source2 = os.path.join(dir2, f2)
                dest2 = os.path.join(dir2_fix, f1)
                cont, skip = progress.Update(n, "Swapping %s" % well1)
                if not cont:
                    break
                shutil.copyfile(source1, dest1)
                shutil.copyfile(source2, dest2)
        except Exception, msg:
            error_dialog(self, msg)
        finally:
            try:
                progress.Destroy()
            except:
                pass
            event.Skip()

    def OnExit(self, event):
        """Explicitly destroy this frame"""
        self.Destroy()
        self.app.Exit()


class SwapPlateApp(wx.App):
    """The wx gui application manager for the swap plate app"""
    def OnInit(self):
        """Required setup for the applicaiton"""
        wx.InitAllImageHandlers()

        # set a top level frame and show it
        main_frame = SwapPlateMainFrame(None, -1, "")
        self.SetTopWindow(main_frame)
        main_frame.Show()

        # must retrun True or the applicaiton will exit
        return True


if __name__ == "__main__":

    def init_gui():
        """Run the swap_plate gui interface"""
        app = SwapPlateApp(redirect=False)
        app.MainLoop()

    init_gui()
