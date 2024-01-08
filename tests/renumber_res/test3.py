"""

This version only works for single chain renumbering
"""
import sys
import Bio.PDB.StructureBuilder
from Bio.PDB import PDBIO, PDBParser, PPBuilder, CaPPBuilder
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


def add_res_name(res_name, res_short):
    """
    In case any abnormal amino acid in the PDB file
    :param res_name: 3 capital letters of the residue
    :param res_short:  shor name of the residue
    :return:
    """
    if len(res_name) != 3:
        raise Exception("Invalid Key")
    if res_name in amino_acid_map.keys():
        raise Exception("Key already exist")
    amino_acid_map[res_name] = res_short


def in_range(number, range_list):
    is_in_range = False
    range_index = -1
    for index, r in enumerate(range_list):
        if number in r:
            is_in_range = True
            return is_in_range, index
    return is_in_range, range_index


def renumber(fasta_seq, PDB_structure, chain_name=None, new_PDB_name=None):
    # check number of chain in the structure
    if not chain_name:
        chain_id = [chain.id for chain in structure.get_chains()]
        assert len(chain_id) == 1
        chain_name = chain_id[0]
    else:
        chain_name = chain_name

    PDB_chain = PDB_structure[0][chain_name]
    PDB_chain_seq = PPBuilder().build_peptides(PDB_chain,aa_only=False)[0].get_sequence()
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
            res.id = (h, query_target_diff[range_index] + idx + 1, ins)
            new_chain.add(res)
        else:
            res.id = (h, num, ins)
            new_chain.add(res)
    PDB_structure[0].add(new_chain)
    assert new_PDB_name is not None
    io = PDBIO()
    io.set_structure(PDB_structure)
    io.save(f"{new_PDB_name}.pdb")


if __name__ == '__main__':
    fasta = next(
        SeqIO.parse('/projectnb2/docking/imhaoyu/my_protein_toolbox/tests/renumber_res/case_3/7POW.fa', 'fasta'))
    fasta_seq = fasta.seq
    parser = PDBParser()
    structure = parser.get_structure('7POW',
                                     '/projectnb2/docking/imhaoyu/my_protein_toolbox/tests/renumber_res/case_3/7POW.pdb')
    renumber(fasta_seq, structure, chain_name='A', new_PDB_name='/projectnb2/docking/imhaoyu/my_protein_toolbox/tests/renumber_res/case_3/7POW_A')
   # with open('/projectnb2/docking/imhaoyu/my_protein_toolbox/tests/renumber_res/case_2/7MDE.pdb') as handle:
    #    pdb_seq = next(SeqIO.parse(handle, "pdb-atom"))
    #renumber(fasta_seq, structure, chain_name='A', new_PDB_name='7MDE_AB')
