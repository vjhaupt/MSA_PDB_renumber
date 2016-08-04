#! /usr/bin/python
# -*- coding: utf-8 -*-
#@author: V. Joachim Haupt
#@contact: joachimh@biotec.tu-dresden.de
#@organization: BIOTEC/TU Dresden




from __future__ import division
from kendrew.toolchain.files import read, fileprefix

from Bio import AlignIO
import getopt
import sys, os
from argparse import ArgumentParser

from kendrew.toolchain.mp import parallel_fn
from kendrew.alignLib.Align import Needle
from kendrew.alignLib.Evaluate import NeedleSingle
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def usage():
    """Prints the help message for the script."""
    print """USAGE:

    NAME
       msaPDBrenumber.py

    SYNOPSIS
        msaPDBrenumber.py -a msa.fasta
                          -p 1abc.pdb

    OPTIONS
       -a           multiple sequence alignment in fasta format
       -p           pdb file to renumber
       -v           be verbose
       --renumligs  Activate renumbering of ligand chain positions (to 0)

    """


def main(verbose=False, msafile=None, pdbfile=None, renumligs=False):

    onelettercode = {'ASP':'D', 'GLU':'E', 'ASN':'N', 'GLN':'Q',  'ARG':'R', 'LYS':'K', 'PRO':'P', 'GLY':'G', 'CYS':'C', 'THR':'T', 'SER':'S', 'MET':'M', 'TRP':'W', 'PHE':'F', 'TYR':'Y', 'HIS':'H', 'ALA':'A', 'VAL':'V', 'LEU':'L', 'ILE':'I'}

    wd = os.getcwd()
    msa = AlignIO.read(read(msafile), 'fasta')
    querypdbid = fileprefix(pdbfile)
    querypdb = PDBParser(PERMISSIVE=1).get_structure(querypdbid, pdbfile)
    #for each chain, get it's sequence
    qseqdicts = dict()
    chains = dict( (chain.get_id(), chain) for chain in querypdb.get_list()[0].get_list() )
    for chainid in chains:
        qseqdicts[chainid] = dict()
        residuelist = chains[chainid].get_list()
        for i, residue in enumerate(residuelist):
            if residue.get_resname() not in onelettercode:
                oletter = 'X'
            else:
                oletter = onelettercode[residue.get_resname()]
            qseqdicts[chainid][int(residue.get_id()[1])] = oletter
    #get seqs as strings
    for chainid in qseqdicts:
        print "\nCurrent chain is '%s'." % chainid
        qseqdict = qseqdicts[chainid]
        qseq = "".join( [qseqdict[resnr] for resnr in sorted(qseqdict.keys())] )
        qseqkeys = [resnr for resnr in sorted(qseqdict.keys())]
        qseq = SeqRecord(Seq(qseq,IUPAC.protein),id=querypdbid)
        # determine the sequence closest to the one in the PDB file
        i=0
        maxid = 0 # maximum sequence identity
        maxal = maxmsaseq = None
        # align qseq to each seq in the msa and find the one with highest identity
        for alignment in msa:
            msaseq = SeqRecord(alignment.seq, id=alignment.id)
            ids = (qseq.id, msaseq.id)
            ali = Needle(pair=ids,records=(qseq, msaseq), name=ids, confopts={'gapopen':10, 'gapextend':0, 'brief':0}) # wd='/tmp/needle'
            ali.align()
            seqid = NeedleSingle(outfile=ali.outfile).getResults()['Longest_Identity']
            if seqid > maxid:
                maxid = seqid
                maxal = AlignIO.read(ali.outfile, 'emboss')
                maxmsaseq = msaseq
        print "Highest seq id is %.2f for:\n%s\n================" % (maxid, str(maxal))
        qseq = str(maxal[0].seq)
        tseq = str(maxal[1].seq)
        msaseq = str(maxmsaseq.seq)

        # number letters in the msa sequence
        # get the positions without gaps, this gives the actual positions in tseq
        msaseqindices = [  i for i in range(1,len(msaseq)+1) if msaseq[i-1] != '-' ]


        #   get positions with gaps
        # qseqgappos = [ i for i in range(len(qseq)) if tseq[i] == '-' ]
        # qseqgappos.reverse()
        # qseqkeys = range(len(qseq))
        # for i in qseqgappos:
        #     del qseqkeys[i]


        # trim qseq if some pos there are mapped to gaps in tseq
        # expand the list of qseqkeys to align with MSA
        tmp = list(qseqkeys)
        qseqkeys = []
        tmp.reverse()
        for i in range(len(qseq)):
            if qseq[i] != '-':
                key = tmp.pop()
            else:
                key = None
            if tseq[i] != '-':
                qseqkeys.append(key)

        qseq = [ qseq[i] for i in range(len(qseq)) if tseq[i] != '-' ]
        # and trim gaps from tseq
        tseq = [ tseq[i] for i in range(len(tseq)) if tseq[i] != '-' ]
        # now, tseq is guaranteed gap-free

        #replace each position in qseqindices by None if there is a gap in qseq
        newqseqindices = list(msaseqindices)
        for i in range(len(qseq)):
            if qseq[i] == '-':
                newqseqindices[i] = None
        #print newqseqindices
        # now newqseqindices is the list of indices for the renumbering in the order of the query sequence

        #sanity check
        assert(len(qseqkeys) == len(qseq) == len(tseq) == len(newqseqindices))
        for i in range(len(qseqkeys)):
            if qseqkeys[i] is None:
                assert(newqseqindices[i] is None)

        #remove empty positions and create a dict of the remapping
        newqseqindices = [index for index in newqseqindices if index is not None]
        qseqkeys = [index for index in qseqkeys if index is not None]

        newqseqindices = dict(zip(qseqkeys,newqseqindices))

        #renumber residues
        chain = chains[chainid]
        for residue in chain.get_list():
            resid = list(residue.get_id())
            index = resid[1]
            if index in newqseqindices:
                resid[1] = newqseqindices[index]
            elif renumligs:
                 resid[1] = 0
            residue.id = tuple(resid)

    # write renumbered PDB file
    io = PDBIO()
    io.set_structure(querypdb)
    io.save(os.path.join(wd, '%s_renum.pdb' % querypdbid))


# Parse command line options
if __name__ == '__main__':
    parser = ArgumentParser(description="MSA PDB Renumbering")
    parser.add_argument('-v', help="Activate verbose mode", dest="verbose", action="store_true", default=False)
    parser.add_argument('-a', help="Specify MSA file", dest="msafile", required=True)
    parser.add_argument('-p', help="Specifiy PDB file", dest="pdbfile", required=True)
    parser.add_argument('--renumligs', help="Renumber ligand chain positions (to 0)", dest="renumligs", action="store_true", default=False)
    args = parser.parse_args()
    
    main(verbose=args.verbose, msafile=args.msafile, pdbfile=args.pdbfile, renumligs=args.renumligs)
#end
