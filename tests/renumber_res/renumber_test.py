import sys
import Bio.PDB.StructureBuilder
from Bio.PDB import PDBIO, PDBParser, PPBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio import SeqIO
from Bio import Align

import warnings

warnings.filterwarnings('ignore')


def three2one(amino_acid):
    amino_acid_map = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
    }
    return amino_acid_map.get(amino_acid.upper(), "")


fasta_seq = next(SeqIO.parse('/projectnb2/docking/imhaoyu/my_protein_toolbox/renumber_res/case_0/4DKL.fa', 'fasta'))
print(fasta_seq.seq)

# with open('/projectnb2/docking/imhaoyu/my_protein_toolbox/renumber_res/case_1/4YO0.pdb') as handle:
#    pdb_seq = next(SeqIO.parse(handle, "pdb-atom"))

parser = PDBParser()
structure = parser.get_structure('4YO0', '/projectnb2/docking/imhaoyu/my_protein_toolbox/renumber_res/case_0/4DKL.pdb')



# chain.parent.detach_child(chain.id)
# new_chain = Chain("A")

# for idx, res in enumerate(chain):
#    res = res.copy()
#    h, num, ins = res.id
#    if 0 <= idx <= 442:
#        res.id = (h, 13 + idx, ins)
#        new_chain.add(res)
# structure[0].add(new_chain)
# io = PDBIO()
# io.set_structure(structure)
# io.save("4DKL.pdb")


# aligner = Align.PairwiseAligner()
# alignments = aligner.align(fasta_seq, chain_sequence)
# print(alignments[0].aligned)

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
    with open('/projectnb2/docking/imhaoyu/my_protein_toolbox/tests/renumber_res/case_2/7MDE.pdb') as handle:
        pdb_seq = next(SeqIO.parse(handle, "pdb-atom"))
    PDB_chain_seq = pdb_seq.seq #PPBuilder().build_peptides(PDB_chain)[0].get_sequence()
    PDB_chain_seq=(PDB_chain_seq.replace('X',''))

    aligner = Align.PairwiseAligner()
    alignments = aligner.align(PDB_chain_seq, fasta_seq)
    print(alignments[0])
    query_target_diff = [alignments[0].aligned[1][:][0][0] - alignments[0].aligned[0][:][0][0]]
    range_list = [range(inner_list[0], inner_list[-1] +1) for inner_list in alignments[0].aligned[0]]
    PDB_chain.parent.detach_child(PDB_chain.id)
    new_chain = Chain(chain_name)
    for idx, res in enumerate(PDB_chain):
        is_in_range, range_index = in_range(idx, range_list)
        res = res.copy()
        h, num, ins = res.id
        if is_in_range:
            res.id = (h, query_target_diff[range_index] + idx+1, ins)
            new_chain.add(res)
        else:
            res.id = (h, num, ins)
            new_chain.add(res)
    PDB_structure[0].add(new_chain)
    assert new_PDB_name is not None
    io = PDBIO()
    io.set_structure(PDB_structure)
    io.save(f"{new_PDB_name}.pdb")


if __name__=='__main__':
    fasta = next(SeqIO.parse('/projectnb2/docking/imhaoyu/my_protein_toolbox/tests/renumber_res/case_2/7MDE.fa', 'fasta'))
    fasta_seq = fasta.seq
    parser = PDBParser()
    structure = parser.get_structure('7MDE',
                                     '/projectnb2/docking/imhaoyu/my_protein_toolbox/tests/renumber_res/case_2/7MDE.pdb')
    renumber(fasta_seq,structure,chain_name='A',new_PDB_name='7MDE_AB')
    with open('/projectnb2/docking/imhaoyu/my_protein_toolbox/tests/renumber_res/case_2/7MDE.pdb') as handle:
        pdb_seq = next(SeqIO.parse(handle, "pdb-atom"))
    renumber(fasta_seq, structure, chain_name='A', new_PDB_name='7MDE_AB')
