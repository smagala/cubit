#!/usr/bin/env python
# -*- coding: utf-8 -*-
# generated by wxGlade 0.5 on Wed Feb  4 09:38:53 2009 from /home/smagala/tasks/swap_plates/gui/swap_plates_gui.wxg

import wx

class SwapPlateGUIMainFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: SwapPlateGUIMainFrame.__init__
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.panel_1 = wx.Panel(self, -1)
        
        # Menu Bar
        self.main_frame_menubar = wx.MenuBar()
        self.SetMenuBar(self.main_frame_menubar)
        self.file_menu = wx.Menu()
        self.file_menu.AppendSeparator()
        self.exit_menu_option = wx.MenuItem(self.file_menu, wx.NewId(), "Exit", "Quit the Swap Plate Application", wx.ITEM_NORMAL)
        self.file_menu.AppendItem(self.exit_menu_option)
        self.main_frame_menubar.Append(self.file_menu, "File")
        # Menu Bar end
        self.main_frame_statusbar = self.CreateStatusBar(1, 0)
        self.label_1 = wx.StaticText(self.panel_1, -1, "Dir for plate 1:", style=wx.ALIGN_RIGHT)
        self.dir1 = wx.TextCtrl(self.panel_1, -1, "", style=wx.TE_READONLY)
        self.dir1_button = wx.Button(self.panel_1, -1, "Select")
        self.label_2 = wx.StaticText(self.panel_1, -1, "Dir for plate 2:", style=wx.ALIGN_RIGHT)
        self.dir2 = wx.TextCtrl(self.panel_1, -1, "", style=wx.TE_READONLY)
        self.dir2_button = wx.Button(self.panel_1, -1, "Select")
        self.run_button = wx.Button(self.panel_1, -1, "Swap Plates")

        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_MENU, self.OnExit, self.exit_menu_option)
        self.Bind(wx.EVT_BUTTON, self.OnGetDir1, self.dir1_button)
        self.Bind(wx.EVT_BUTTON, self.OnGetDir2, self.dir2_button)
        self.Bind(wx.EVT_BUTTON, self.OnRunSwap, self.run_button)
        # end wxGlade

    def __set_properties(self):
        # begin wxGlade: SwapPlateGUIMainFrame.__set_properties
        self.SetTitle("Swap Plates")
        self.main_frame_statusbar.SetStatusWidths([-1])
        # statusbar fields
        main_frame_statusbar_fields = [""]
        for i in range(len(main_frame_statusbar_fields)):
            self.main_frame_statusbar.SetStatusText(main_frame_statusbar_fields[i], i)
        self.dir1.SetMinSize((240, 27))
        self.dir2.SetMinSize((240, 27))
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: SwapPlateGUIMainFrame.__do_layout
        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_2 = wx.FlexGridSizer(2, 1, 5, 0)
        sizer_2 = wx.BoxSizer(wx.VERTICAL)
        grid_sizer_1 = wx.FlexGridSizer(3, 3, 5, 5)
        grid_sizer_1.Add(self.label_1, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.dir1, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.dir1_button, 0, 0, 0)
        grid_sizer_1.Add(self.label_2, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.dir2, 0, wx.ALIGN_CENTER_VERTICAL, 0)
        grid_sizer_1.Add(self.dir2_button, 0, 0, 0)
        sizer_2.Add(grid_sizer_1, 1, wx.ALL, 5)
        grid_sizer_2.Add(sizer_2, 1, wx.EXPAND, 0)
        grid_sizer_2.Add(self.run_button, 0, wx.ALIGN_CENTER_HORIZONTAL, 0)
        self.panel_1.SetSizer(grid_sizer_2)
        sizer_1.Add(self.panel_1, 1, wx.EXPAND, 5)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        self.Layout()
        # end wxGlade

    def OnExit(self, event): # wxGlade: SwapPlateGUIMainFrame.<event_handler>
        print "Event handler `OnExit' not implemented!"
        event.Skip()

    def OnGetDir1(self, event): # wxGlade: SwapPlateGUIMainFrame.<event_handler>
        print "Event handler `OnGetDir1' not implemented!"
        event.Skip()

    def OnGetDir2(self, event): # wxGlade: SwapPlateGUIMainFrame.<event_handler>
        print "Event handler `OnGetDir2' not implemented!"
        event.Skip()

    def OnRunSwap(self, event): # wxGlade: SwapPlateGUIMainFrame.<event_handler>
        print "Event handler `OnRunSwap' not implemented!"
        event.Skip()

# end of class SwapPlateGUIMainFrame


class SwapPlateApp(wx.App):
    def OnInit(self):
        wx.InitAllImageHandlers()
        main_frame = SwapPlateGUIMainFrame(None, -1, "")
        self.SetTopWindow(main_frame)
        main_frame.Show()
        return 1

# end of class SwapPlateApp

if __name__ == "__main__":
    app = SwapPlateApp(0)
    app.MainLoop()
