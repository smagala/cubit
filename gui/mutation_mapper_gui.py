#!/usr/bin/env python
# -*- coding: utf-8 -*-
# generated by wxGlade 0.6.3 on Thu Sep 10 12:38:48 2009

import wx

from newick.tree import *

from gui.common import pick_save_file, pick_infile, error_dialog, save_file


class GUIMainFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: GUIMainFrame.__init__
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)
        self.panel_1 = wx.Panel(self, -1)
        
        # Menu Bar
        self.main_frame_menubar = wx.MenuBar()
        self.file_menu = wx.Menu()
        self.file_menu.AppendSeparator()
        self.exit_menu_option = wx.MenuItem(self.file_menu, wx.NewId(), "Exit", "Quit the application", wx.ITEM_NORMAL)
        self.file_menu.AppendItem(self.exit_menu_option)
        self.main_frame_menubar.Append(self.file_menu, "File")

        wxglade_tmp_menu = wx.Menu()
        self.transpose_table_option = wx.MenuItem(wxglade_tmp_menu, wx.NewId(), "Transpose table", "Work around Excel's column limit for large alignments", wx.ITEM_CHECK)
        wxglade_tmp_menu.AppendItem(self.transpose_table_option)
        self.report_aa_option = wx.MenuItem(wxglade_tmp_menu, wx.NewId(), "Report AA changes", "Show amino acids as well as positions in report", wx.ITEM_CHECK)
        wxglade_tmp_menu.AppendItem(self.report_aa_option)
        self.main_frame_menubar.Append(wxglade_tmp_menu, "Options")
        
        self.help_menu = wx.Menu()
        self.help_menu_tutorial = wx.MenuItem(self.help_menu, wx.NewId(), "Tutorial", "Information on program", wx.ITEM_NORMAL)
        self.help_menu.AppendItem(self.help_menu_tutorial)
        self.main_frame_menubar.Append(self.help_menu, "Help")

        self.SetMenuBar(self.main_frame_menubar)
        # Menu Bar end


        self.main_frame_statusbar = self.CreateStatusBar(1, 0)

        
        self.file_box = wx.TextCtrl(self.panel_1, -1, "", style=wx.TE_READONLY)


        self.reference_option = wx.TreeCtrl(self.panel_1, 1, wx.DefaultPosition, (350,200), 
                                         wx.TR_HIDE_ROOT|wx.TR_HAS_BUTTONS)
        
        #self.reference_option = wx.ComboBox(self.panel_1, pos=(150, 90), size=(95, -1), 
        #                                         choices=[], 
        #style=wx.CB_DROPDOWN|wx.CB_READONLY|wx.CB_SORT)
        #self.Bind(wx.EVT_COMBOBOX, self.OnSelect)

        self.referenceSelect_button = wx.Button(self.panel_1, -1, "Select References")
        
        self.nwkfile_box = wx.TextCtrl(self.panel_1, -1, "", style=wx.TE_READONLY)
        
        # Buttons
        self.nwkSelect_button = wx.Button(self.panel_1, -1, "Select")
        self.fasSelect_button = wx.Button(self.panel_1, -1, "Select")
        self.loadReferences_button = wx.Button(self.panel_1, -1, "Load References")
        
        self.consensus_button = wx.Button(self.panel_1, -1, "Generate Reference Consensus")
        self.table_button = wx.Button(self.panel_1, -1, "Generate AA Diff Table")
        self.report_button = wx.Button(self.panel_1, -1, "Generate AA Change Report")
        self.mutation_button = wx.Button(self.panel_1, -1, "Generate Mutation Tree")
        self.orderFasta_button = wx.Button(self.panel_1, -1, "Download Ordered Fasta")
        
        self.__set_properties()
        self.__do_layout()

        self.Bind(wx.EVT_MENU, self.on_exit, self.exit_menu_option)
        self.Bind(wx.EVT_MENU, self.on_tutorial, self.help_menu_tutorial)
        self.Bind(wx.EVT_BUTTON, self.on_select_alignment, self.fasSelect_button)
        self.Bind(wx.EVT_BUTTON, self.on_select_newick, self.nwkSelect_button)        
        self.Bind(wx.EVT_BUTTON, self.on_generate_consensus, self.consensus_button)
        self.Bind(wx.EVT_BUTTON, self.on_generate_table, self.table_button)
        self.Bind(wx.EVT_BUTTON, self.on_generate_report, self.report_button)
        self.Bind(wx.EVT_BUTTON, self.on_generate_mutation, self.mutation_button)
        self.Bind(wx.EVT_BUTTON, self.on_select_reference, self.referenceSelect_button)
        self.Bind(wx.EVT_BUTTON, self.on_load_references, self.loadReferences_button)
        self.Bind(wx.EVT_BUTTON, self.on_order_fasta, self.orderFasta_button)

    def __set_properties(self):
        self.SetTitle("Mutation Mapping Tool")
        self.main_frame_statusbar.SetStatusWidths([-1])
        # statusbar fields
        main_frame_statusbar_fields = [""]
        for i in range(len(main_frame_statusbar_fields)):
            self.main_frame_statusbar.SetStatusText(main_frame_statusbar_fields[i], i)
        self.file_box.SetMinSize((250, 27))
        self.nwkfile_box.SetMinSize((250, 27))        
        self.table_button.Enable(False)
        self.consensus_button.Enable(False)
        self.report_button.Enable(False)
        self.mutation_button.Enable(False)
        self.orderFasta_button.Enable(False)

    def __do_layout(self):
        afntSizer = wx.StaticBoxSizer(wx.StaticBox(self.panel_1, wx.NewId(), 
                                                   'PICK'), wx.VERTICAL)
        afntGridSizer = wx.FlexGridSizer(3, 3, 5, 5)
        afntGridSizer.AddGrowableCol(1)
        
        afntGridSizer.Add(wx.StaticText(self.panel_1, -1, "nt Alignment File:"), 0, wx.RIGHT|wx.ALIGN_CENTER_VERTICAL, 5)
        afntGridSizer.Add(self.file_box, 0, wx.EXPAND, 0)
        afntGridSizer.Add(self.fasSelect_button, 0, wx.LEFT|wx.ALIGN_CENTER, 5)

        afntGridSizer.Add(wx.StaticText(self.panel_1, -1, "Newick Tree:"), 0, wx.RIGHT|wx.ALIGN_CENTER_VERTICAL, 5)
        afntGridSizer.Add(self.nwkfile_box, 0, wx.EXPAND, 0)
        afntGridSizer.Add(self.nwkSelect_button, 0, wx.LEFT|wx.ALIGN_CENTER, 5)

        afntGridSizer.Add(self.loadReferences_button, 0, wx.ALL|wx.ALIGN_CENTER, 5)

        sizer_6 = wx.BoxSizer(wx.VERTICAL)
        sizer_6.Add(self.reference_option, 0, wx.EXPAND, 0)
        sizer_5 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_5.Add(wx.StaticText(self.panel_1, -1, "Pick Reference:"), 0, wx.RIGHT|wx.ALIGN_CENTER_VERTICAL, 5)
        sizer_5.Add(sizer_6, 1, wx.ALIGN_CENTER_VERTICAL, 0)
        sizer_5.Add(self.referenceSelect_button, 0, wx.LEFT|wx.ALIGN_CENTER_VERTICAL, 5)
        
        afntSizer.Add(afntGridSizer, 0, wx.ALL|wx.ALIGN_CENTER|wx.EXPAND, 5)

        sizer_2 = wx.BoxSizer(wx.VERTICAL)
        sizer_2.Add(afntSizer, 0, wx.ALL|wx.EXPAND, 5)
        sizer_2.Add(sizer_5, 0, wx.ALL|wx.EXPAND, 5)
        
        grid_sizer_1 = wx.GridSizer(1, 5, 5, 5)
        grid_sizer_1.Add(self.consensus_button, 1, wx.ALL|wx.EXPAND, 5)
        grid_sizer_1.Add(self.table_button, 1, wx.ALL|wx.EXPAND, 5)
        grid_sizer_1.Add(self.report_button, 1, wx.ALL|wx.EXPAND, 5)
        grid_sizer_1.Add(self.mutation_button, 1, wx.ALL|wx.EXPAND, 5)
        grid_sizer_1.Add(self.orderFasta_button, 1, wx.ALL|wx.EXPAND, 5)
        sizer_2.Add(grid_sizer_1, 0, wx.ALL|wx.ALIGN_CENTER, 10)
        self.panel_1.SetSizer(sizer_2)

        sizer_1 = wx.BoxSizer(wx.VERTICAL)
        sizer_1.Add(self.panel_1, 1, wx.ALL|wx.ALL|wx.EXPAND, 10)
        self.SetSizer(sizer_1)
        sizer_1.Fit(self)
        self.Layout()

    def on_exit(self, event):
        print "Event handler `on_exit' not implemented!"
        event.Skip()

    def on_tutorial(self, event):
        print "Event handler `on_exit' not implemented!"
        event.Skip()

    def on_select_alignment(self, event): # wxGlade: GUIMainFrame.<event_handler>
        print "Event handler `on_select_alignment' not implemented!"
        event.Skip()

    def on_generate_table(self, event): # wxGlade: GUIMainFrame.<event_handler>
        print "Event handler `on_generate_table' not implemented!"
        event.Skip()

    def on_generate_consensus(self, event): # wxGlade: GUIMainFrame.<event_handler>
        print "Event handler `on_generate_consensus' not implemented!"
        event.Skip()

    def on_generate_report(self, event): # wxGlade: GUIMainFrame.<event_handler>
        print "Event handler `on_generate_report' not implemented!"
        event.Skip()

    def on_generate_mutation(self, event): # wxGlade: GUIMainFrame.<event_handler>
        print "Event handler `on_generate_mutation' not implemented!"
        event.Skip()

    def on_select_newick(self, event): # wxGlade: GUIMainFrame.<event_handler>
        print "Event handler `on_select_newick' not implemented!"
        event.Skip()

    def on_order_fasta(self, event):
        print "Event handler 'on_order_fasta' not implemented!"
        event.Skip()
		
    def on_load_references(self, event):
        print "Event handler for load references was not implemented..."
        event.Skip()
    
    def on_select_reference(self, event):
        try:
            referenceIndex = self.reference_option.GetSelection()
            referenceTree = self.reference_option.GetItemData(referenceIndex).GetData()
            self._referencesSelected(referenceTree.get_leaves_identifiers())
            self._set_button_state(A=True, B=False, C=False, D=False, E=True, F=True, G=True, H=True, I=True)
        except Exception, msg:
            error_dialog(self, msg=msg)
        event.Skip()
 
    def _set_button_state(self, A=True, B=False, C=False, D=False, E=False, F=False, G=False, H=False, I=False):
        ''' Make sure the order matches...'''
        self.fasSelect_button.Enable(A)
        self.nwkSelect_button.Enable(B)
        self.loadReferences_button.Enable(C)
        self.referenceSelect_button.Enable(D)
        self.table_button.Enable(E)
        self.report_button.Enable(F)
        self.mutation_button.Enable(G)
        self.orderFasta_button.Enable(H)
        self.consensus_button.Enable(I)
        
    def _referencesSelected(self, referenceList):
        print "Rerefence selected..."

    def _addReferences(self, tree):
        self.clear_tree()
        parent = self.reference_option.AddRoot('Root')
        self._constructTree(parent, tree)

    def _constructTree(self, parent, tree):
        pytreeitem = wx.TreeItemData()
        pytreeitem.SetData(tree)
        if not isinstance(tree, Leaf):
            np = self.reference_option.AppendItem(parent, '+', -1, -1, pytreeitem)
            for t,_,_ in tree.get_edges():
                self._constructTree(np, t)
        else:
            self.reference_option.AppendItem(parent, tree.identifier, -1, -1, pytreeitem)

    def clear_tree(self):
        self.reference_option.DeleteAllItems()

# end of class GUIMainFrame


class GUIApp(wx.App):
    def OnInit(self):
        wx.InitAllImageHandlers()
        main_frame = GUIMainFrame(None, -1, "")
        self.SetTopWindow(main_frame)
        main_frame.Show()
        return 1

# end of class GUIApp

if __name__ == "__main__":
    app = GUIApp(0)
    app.MainLoop()
