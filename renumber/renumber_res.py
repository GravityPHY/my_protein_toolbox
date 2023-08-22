"""

This version only works for single chain renumbering
"""
import sys
import Bio.PDB.StructureBuilder
from Bio.PDB import PDBIO, PDBParser, PPBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio import SeqIO, Align

import warnings
warnings.filterwarnings('ignore')

amino_acid_map = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

def in_range(number, range_list):
    is_in_range = False
    range_index = -1
    for index, r in enumerate(range_list):
        if number in r:
            is_in_range = True
            return is_in_range, index
    return is_in_range, range_index


def renumber(fasta_seq, PDB_structure, new_PDB_name):
    # check number of chain in the structure
    chain_id = [chain.id for chain in structure.get_chains()]
    assert len(chain_id) == 1
    chain_name = chain_id[0]
    PDB_chain = PDB_structure[0][chain_name]
    PDB_chain_seq = PPBuilder().build_peptides(PDB_chain)[0].get_sequence()
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(PDB_chain_seq, fasta_seq)
    query_target_diff = [alignments[0].aligned[1][:][0][0] - alignments[0].aligned[0][:][0][0]]
    range_list = [range(inner_list[0], inner_list[-1] + 1) for inner_list in alignments[0].aligned[0]]
    PDB_chain.parent.detach_child(PDB_chain.id)
    new_chain = Chain(chain_name)
    for idx, res in enumerate(PDB_chain):
        is_in_range, range_index = in_range(idx, range_list)
        res = res.copy()
        h, num, ins = res.id
        if is_in_range:
            res.idx = (h, query_target_diff[range_index] + idx, ins)
            new_chain.add(res)
        else:
            res.id = (h, num, ins)
            new_chain.add(res)

    io = PDBIO()
    io.set_structure(PDB_structure)
    io.save(f"{new_PDB_name}.pdb")