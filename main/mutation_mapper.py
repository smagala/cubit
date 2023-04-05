#!/usr/bin/env python
#file mutation_mapper.py

# Copyright (C) 2012 James Smagala/Melih Gunay

# Contact:
#   James Smagala: smagala@gmail.com
#   Melih Gunay:   gmelih@gmail.com

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

import wx
import os
import os.path
from core.sequence import read_fasta_file, DNASequence, diff2seqs, disambiguate, DNA, IUPAC_DNA
from gui.mutation_mapper_gui import GUIApp, GUIMainFrame
from gui.mutation_mapper_help import GUIHelp
from gui.common import pick_save_file, pick_infile, error_dialog, save_file
import re
import logging

import sys, traceback

from newick import parse_tree
from newick.tree import *

prsMut = re.compile('\d+')

logger = logging.getLogger('mutation-mapper')
if not os.path.exists('log'):
    os.mkdir('log')
hdlr = logging.FileHandler('log/mutation-mapper.log', mode='w')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)
logger.info("GUI setup and initalization")

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

def output_phyloxml(T, length=None):
    tree_str = "<clade>\n"
    if isinstance(T, Leaf):
        tree_str +=" <name>"+T.identifier+"</name>\n"
    if length:
        tree_str+="  <branch_length>"+str(length)+"</branch_length>\n"
    if T.label is not None:
        tree_str+=" <taxonomy>\n"
        tree_str+="   <scientific_name>"+T.label+"</scientific_name>\n"
        tree_str+=" </taxonomy>\n"
    if isinstance(T, Tree):
        for (n,b,l) in T.get_edges():
            tree_str += output_phyloxml(n, l)
    tree_str +='</clade>\n'
    return tree_str

def output_newick(T):
    if isinstance(T, Leaf):
        return "'"+T.identifier+"'"
    tree_str = '('
    sep = ''
    for (n,b,l) in T.get_edges():
        tree_str += sep+str(output_newick(n))
        if b:
            tree_str += str(b) + ' '
        if l:
            tree_str += ':' + str(l)
            if n.label is not None:
                tree_str += "[&Mutations="+n.label+"]"
        sep = ', '
    return tree_str+')'

def output_nexus(T):
    tree_str = "#Nexus\n\n"
    tree_str += "Begin Trees;\n\n"
    tree_str += "tree Annotated_tree ="
    tree_str += output_newick(T)
    tree_str += ";\nend"
    return tree_str

class AADiffMainFrame(GUIMainFrame):
    """The wx main frame for aadiff - implements the event handlers"""

    def __init__(self, *args, **kwargs):
        # call this first so that the underlying C++ object exists
        # call the parent __init__ to handle additional setup
        super(AADiffMainFrame, self).__init__(*args, **kwargs)

        # get a reference to the app object
        self.app = wx.GetApp()

        # bind events to handlers
        self.Bind(wx.EVT_CLOSE, self.on_exit)

        # set an icon if desired
        icon = wx.Icon('icons/mutation_mapper.ico', wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon)
        self.initialize()

    def initialize(self):
        # store the base file name to derive other file names from
        self.base_file_name = ''
        # store the nwk file name to derive other file names from
        self.nwk_file_name = ''
        self.positions = None
        self.reference_str = None
        self.reference = None
        self.seqs = None
        self.nwktree = None
        self.seqs_aa = None
        self.changesets = None
        self.mutations = None   
        self._set_button_state()
        self.nwkfile_box.SetValue('')
        self.file_box.SetValue('')
        self.clear_tree()

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

    def _gen_consensus_table(self):
        """generate the consensus reference sequence"""
        return self.reference.to_string()

    def _gen_mutation_tree(self):
        #nwkstr = str(self.nwktree)+";"
        nwkstr = output_newick(self.nwktree)+";"
        logger.debug( "nwktree : "+nwkstr)
        nexusstr = output_nexus(self.nwktree)+";"
        pxml ="<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        pxml += "<phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd\" xmlns=\"http://www.phyloxml.org\">\n"
        pxml +="<phylogeny rooted=\"false\">\n"
        pxml += output_phyloxml(self.nwktree)
        pxml += "</phylogeny>\n </phyloxml>"
        logger.debug( "phyloxml : "+pxml)
        return nwkstr, nexusstr, pxml
 
    def _determineConsensus(self, seqs):
        if len(seqs) < 1:
            error_dialog(self, "No sequences were provided...")
            return
        seqLength = len(seqs[0])
        cseq = [{} for i in range(seqLength)]
        for i in range(seqLength):
            for each in seqs:
                if not each[i] in cseq[i]:
                    cseq[i][each[i]]=0
                cseq[i][each[i]] += 1
        referenceConsensus = ["" for i in range(seqLength)]
        i = 0
        for each in cseq:
            maxCount = 0
            dna = None
            for key, value in each.items():
                if value > maxCount:
                    maxCount = value
                    dna = key
                elif value == maxCount:
                    dna=dna+"/"+key
            referenceConsensus[i] = dna
            i += 1
            #logger.debug( str(i)+" : "+dna+"-"+str(maxCount)
        return referenceConsensus

    def _referencesSelected(self, referenceNameList):
        self.mutations = {}
        sequencesForSelectedOutgroup = []
        for each in referenceNameList:
            for seq in self.seqs:
                if seq.name == each:
                    #logger.debug( "Sequence matched to the reference..."
                    sequencesForSelectedOutgroup.append(seq)
                    break

        logger.debug("begin outgroup sequences.. >>>>>")
        logger.debug(sequencesForSelectedOutgroup)
        logger.debug("<<< end outgroup sequences.")

        outGroupReferenceConsensus = self._determineConsensus(sequencesForSelectedOutgroup)
        outgroupseq=''
        i = -1
        for each in outGroupReferenceConsensus:
            i += 1
            if "/" in each:
                outgroupseq += self.__decideOnConsensusWhenConflict(each, 
                                                                    self.consensusSequence[i])
            else:
                outgroupseq += each
            
        #outgroupseq = ''.join(outGroupReferenceConsensus)
        #logger.debug( "outgroup seq : "+outgroupseq

        self.reference = DNASequence('OutgroupConsensusSequence', outgroupseq)
        self._set_button_state(A=True, B=True, C=False, D=False, E=False, F=False, G=False, H=False, I=True)
        self._referenceSelected()
        self._addLabels(self.nwktree)

    def __decideOnConsensusWhenConflict(self, current, csDNA):
        #decide on the consensus
        #logger.debug( "Current seq="+current+" , consensus="+csDNA
        for sdna in current.split("/"):
            if sdna == csDNA:
                logger.debug( "matched..."+sdna+"="+csDNA)
                return csDNA
            elif sdna in IUPAC_DNA:
                if csDNA in IUPAC_DNA[sdna]:
                    logger.debug( csDNA+" is in "+ str(IUPAC_DNA[sdna]))
                    return csDNA
            else:
                logger.warn( "Did not match...current="+current+", "+sdna+"="+csDNA)

        #open up all sequences 
        cseq = current.split("/")
        cseq.extend(csDNA.split("/"))
        logger.debug( "All sequences.."+str(cseq))
        winningDNA = self.__getWinningDNAs(cseq)
        logger.debug( "winning dna="+winningDNA)
        if winningDNA is not None and len(winningDNA) == 1:
            return winningDNA 
        # still could not decide, use all sequences
        nseqs = []
        for each in cseq:
            if each in IUPAC_DNA:
                nseqs.extend(IUPAC_DNA[each])
            else:
                nseqs.append(each)
        logger.debug("nseqs="+str(nseqs))
        winningDNA = self.__getWinningDNAs(nseqs)
        logger.debug("winning nseqs="+winningDNA)
        if winningDNA is None:
            logger.error("something is very wrong... check the implementation...")
        if len(winningDNA) == 1:
            return winningDNA
        wseq = winningDNA.split("/")
        wseqS = set(wseq)
        logger.debug("wseqS :"+str(wseqS))
        for key, value in IUPAC_DNA.items():
            if len(wseqS.symmetric_difference(set(value)))==0:
                return key
        raise Exception("Should not arrive here..")

    def __buildOccuranceMap(self, seq):
        '''Given a list of sequences builds occurance map'''
        occuranceMap = {}
        for each in seq:
            if each not in occuranceMap:
                occuranceMap[each]=0
            occuranceMap[each] += 1
        return occuranceMap

    def __getWinningDNAs(self, seq):
        '''Given a list of sequences returns the winning seq
        may include / '''
        occuranceMap = self.__buildOccuranceMap(seq)
        i = 0
        maxCount = 0
        dna = None
        for key, value in occuranceMap.items():
            if value > maxCount:
                maxCount = value
                dna = key
            elif value == maxCount:
                dna=dna+"/"+key
                
        return dna                     

    def _referenceSelected(self):
        self.reference_aa, self.seqs_aa, self.changesets = aadiff(self.reference, self.seqs)
        self._fix_gaps()
        self._calc_aa_strings()
        self.missingSeq = {}
        i = 0
        for seq in self.seqs:
            locations = []
            for k,v in self.changesets[i].iteritems():
                mtns = list(v)
                if '-' in mtns:
                    locations.append(k+1)
            i += 1
            locations.sort()
            self.missingSeq[seq.name]=locations
        i = 0
        for seq in self.seqs:
            self.mutations[seq.name] = self.changesets[i]
            i += 1
            
        self._addMutations(self.nwktree)

    def _addMutations(self, nwktree):
        #Add mutations to the leaves first
        for l in nwktree.get_leaves():
            if l.identifier not in self.mutations:
                raise LookupError(l.identifier, " could not be found!\n Perhaps mismatch between newick and fasta...")
            changes = self.mutations[l.identifier]
            for k,v in changes.iteritems():
                mtns = list(v)
                mtns.sort()
                if '-' in mtns:
                    mtns.remove('-')
                mtns = '%s' % '/'.join(mtns)
                if len(mtns)>0:
                    l.mutation.add(self.reference_str[k]+str(k+1)+mtns)
          
    def _markParentMutation(self, current, leaf, m):
        '''Tree and mutation (T, m) '''
        T = current.parent
        if T is None:
            return current
        for l in T.get_leaves():
            if l == leaf:
                continue
            missingLocations = self.missingSeq[l.identifier]
            mloc = prsMut.search(m)
            skip = (mloc is not None and int(mloc.group()) in missingLocations)
            if not (skip or m in l.mutation):
                return current
        return self._markParentMutation(T, leaf, m)

    def _removeMutationFromLeaves(self, T, m):
        for l in T.get_leaves():
            if m in l.mutation:
                l.mutation.remove(m)

    def _addLabels(self, nwktree):
        #Propage mutations to upper nodes/tree
        for l in nwktree.get_leaves():
            mutset = l.mutation.copy()
            while True:
                try:
                    m = mutset.pop()
                    T = self._markParentMutation(l, l, m)
                    if T.label is None:
                        T.label = str(m)
                    else:
                        T.label = T.label+"+"+str(m)

                    logger.info(str(m) +" in "+str(T))
                    logger.info(T.label)
                    self._removeMutationFromLeaves(T, m)
                    logger.info("mutation "+str(m)+" is removed from "+str(T))
                    #logger.debug( l.mutation
                except KeyError, msg:
                    #logger.error(msg)
                    break
        #logger.debug( nwktree

    def _load_seq(self, label):
        for seq in self.seqs:
            if seq.name == label:
                return seq
        return None
	
    def _order_fasta(self):
        fasFile = []
        ordSeqs = []
        if self.nwktree is None:
            return
        for l in self.nwktree.get_leaves():
            seq = self._load_seq(l.identifier)
            if seq:
                fasFile.append(seq.to_fasta())
                ordSeqs.append(seq)
            else:
                logger.warn("Label '%s' is not in original fasta file" % l.identifier)
                pass
        fasFile = '\n'.join(fasFile)
        return fasFile, ordSeqs

    def on_select_alignment(self, event):
        self.initialize()
        """Read in a fasta file containing an alignment and diff it"""
        path = pick_infile(parent=self, msg='Pick a nt alignment:', wild_card="Fasta (*.fas)|*.fas")
        updated = False
        if path:
            # try to read in the alignment and diff it
            # only the results are actually stored
            try:
                self.seqs = [DNASequence(comment, seq) for comment, seq in
                             read_fasta_file(path, bioedit_cleanup=True)]
                self.consensusSequence = self._determineConsensus(self.seqs)
                updated = True
            except Exception, msg:
                error_dialog(self, msg=msg)
            if updated:
                self.base_file_name = os.path.splitext(path)[0]
                self.file_box.SetValue(path)
                self._set_button_state(A=True, B=True, C=False, D=False, E=False, F=False, G=False, H=False)
        logger.debug("Alignment was selected : "+path)
        event.Skip()
       
    def on_load_references(self, event):
        self.seqs = self._order_fasta()[1]
        self._addReferences(self.nwktree)
        refButtonState = self.reference != None
        self._set_button_state(A=True, B=False, C=False, D=True, E=False, F=False, G=False, H=True, I=refButtonState)
        event.Skip()

    def on_generate_table(self, event):
        """Generate a difference table"""
        table = self._gen_diff_table()
        path = pick_save_file(self, msg='Save difference table:',
                default_file=self.base_file_name + '_aadiff_table.csv')
        if path:
            save_file(path, table)
        event.Skip()

    def on_generate_consensus(self, event):
        """Generate a consensus table"""
        table = self._gen_consensus_table()
        path = pick_save_file(self, msg='Save consensus table:',
                default_file=self.base_file_name + '_reference_consensus.csv')
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

    def on_generate_mutation(self, event):
        """Generate a mutation report"""
        nwk_tr, nexus_tr, xml_tr = self._gen_mutation_tree()
        path = pick_save_file(self, msg='Save mutation tree:',
                default_file=self.nwk_file_name + '_aadiff_mutationtree.trees', 
                              wild_card="Nexus (*.trees)|*.trees|Phyloxml (*.xml)|*.xml|Newick (*.nwk)|*.nwk")
        if path:
            if path.endswith(".xml"):
                save_file(path, xml_tr)
            if path.endswith(".nwk"):
                save_file(path, nwk_tr)
            if path.endswith(".trees"):
                save_file(path, nexus_tr)    
        event.Skip()

    def on_select_newick(self, event):
        """Read in a mutation file containing newick tree"""
        path = pick_infile(parent=self, msg='Pick the newick tree:', wild_card="Newick (*.nwk)|*.nwk")
        updated = False
        if path:
            # try to read the newick tree
            try:
                infile = open(path, 'r')
                content = infile.read()
                nwktree = parse_tree(content)
                add_parent_links(nwktree)
                self.nwktree = nwktree
                updated = True
            except Exception, msg:
                #traceback.logger.debug(_exc(file=sys.stdout))
                #traceback.print_exc(file=sys.stdout)
                logger.error("Could not parse newick file", exc_info=1)
                error_dialog(self, msg=msg)
            finally:
                try:
                    infile.close
                except:
                    pass
            if updated:
                self.nwk_file_name = os.path.splitext(path)[0]
                self.nwkfile_box.SetValue(path)
                self._set_button_state(A=True, B=True, C=True, D=False, E=False, F=False, G=False, H=True)
        event.Skip()
		
    def on_order_fasta(self, event):
        """Generate a new fasta file from the old fasta and tree labels"""
        result = self._order_fasta()[0]
        path = pick_save_file(self, msg='Save new FASTA file:',
                              default_file=self.base_file_name +'_dnd.fas', wild_card="Fasta (*.fas)|*.fas")
        if path:
            save_file(path, result)
        event.Skip()
	
    def on_exit(self, event):
        """Explicitly destroy the main frame, ending the program"""
        self.Destroy()
        self.app.Exit()
        logger.debug("Exiting the flumuation tool ...")
        event.Skip()

    def on_tutorial(self, event):
        GUIHelp(self, -1, "User Manual", logger)

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
        app = AADiffApp(redirect=False)
        app.MainLoop()

    init_gui()
