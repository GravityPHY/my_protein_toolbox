import os
import csv
import argparse
import itertools

from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import PDBIO, PDBParser, MMCIFParser, PPBuilder, CaPPBuilder
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1

import warnings

warnings.filterwarnings('ignore')


def get_structure_cif(cif_file):
    parser = MMCIFParser()
    structure = parser.get_structure('cif', cif_file)
    return structure


def get_structure_pdb(pdb_file):
    parser = PDBParser()
    structure = parser.get_structure('PDB', pdb_file)
    return structure


def pdb_to_seq(pdb_file):
    """
    Give a path of .pdb with single chain, parse structure to get the sequence
    :param pdb_file: path of pdb file, should end with .pdb
    :return: amino acid sequence in Seq object
    """
    parser = PDBParser()
    PDB_structure = parser.get_structure('PDB', pdb_file)
    for model in PDB_structure:
        for chain_name in model:
            PDB_chain = PDB_structure[0][chain_name.id]
    PDB_chain_seq = PPBuilder().build_peptides(PDB_chain, aa_only=False)[0].get_sequence()
    return PDB_chain_seq


def get_Seq(structure):
    """
    Retrieve the AA sequence from a PDB structure as a Seq
    :param structure:
    :return: Seq object defined in Biopython
    """
    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
    seq = [_aainfo(r)[1] for r in structure.get_residues() if is_aa(r)]
    return Seq(''.join(seq))


def get_sequence_bfactor(structure, sum=False):
    """
    Retrieve the AA sequence from a structure and the b-factor of the residue (CA atom)
    :param structure:
    :param sum: the option of summing all b-factor values in a residue. True enable the summation, False disabled.
    :return: list of (amino acid, bfactor) tuple. For example [('A',0.1),('G',0.8),('H',0)]
    """
    _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
    seq_bfactor = []
    for r in structure.get_residues():
        bfactor = 0
        counter = 0
        if sum:
            for atom in r.get_atoms():
                bfactor += atom.get_bfactor()
            # if atom.get_id() == 'CA':
            #    bfactor = atom.get_bfactor()
            # seq_bfactor = [(_aainfo(r)[1], next(r.get_atoms()).get_bfactor()) for r in structure.get_residues() if is_aa(r)]
            seq_bfactor.append((_aainfo(r)[1], bfactor))
        else:
            seq_bfactor = [(_aainfo(r)[1], next(r.get_atoms()).get_bfactor()) for r in structure.get_residues()
                           if is_aa(r)]
    return seq_bfactor


def sum_atom_bfactor_residue(structure, pdb_file):
    """
    Only apply to AbEMap pdb at the moment, will save the structure in .pdb to pdb_file
    :param structure:
    :return:None
    """
    for r in structure.get_residues():
        bfactor = 0
        for atom in r.get_atoms():
            bfactor += atom.get_bfactor()
        for atom in r.get_atoms():
            atom.set_bfactor(bfactor)
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)


def align_sequence_bfactor(aligned_seq_list, seq_bfactor_list):
    """
    update the seq_bfactor list with gap residue (represent in '0')
    :param aligned_seq_list:
    :param seq_bfactor_list: [('res', bfactor value)]
    :return:
    """
    index = 0
    for token in aligned_seq_list:
        if token == '-':
            seq_bfactor_list.insert(index, ('0', 0.0))
            index += 1
        else:
            index += 1


def align_seq_bfactor(structure_pred, structure_true, sum):
    """
    Align the unbound sequence with bound sequence, and return the aligned (res, bfactor value) list
    :param structure_pred: structure with predict epitope labels
    :param structure_true: structure with true epitope labels
    :param sum: option to sum the atom bfactor or not. True - sum, False get the residue bfactor
    :return: both the prediction and true (res, bfactor value) lists
    """
    seq_pred = get_Seq(structure_pred)
    bfactor_pred = get_sequence_bfactor(structure_pred, sum=sum)
    seq_true = get_Seq(structure_true)
    bfactor_true = get_sequence_bfactor(structure_true)
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(seq_pred, seq_true)
    # print(alignments[0])
    align_sequence_bfactor(alignments[0][0, :], bfactor_pred)
    align_sequence_bfactor(alignments[0][1, :], bfactor_true)
    return bfactor_pred, bfactor_true


def combine_bfactor(bfactor_true, bfactor_pred1, bfactor_pred2):
    new = []
    for i, j, k in zip(bfactor_true, bfactor_pred1, bfactor_pred2):
        new.append((i[0], max(j[1], k[1])))
    return new