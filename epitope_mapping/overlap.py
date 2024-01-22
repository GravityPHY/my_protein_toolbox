import sys

sys.path.append("/projectnb2/docking/imhaoyu/my_protein_toolbox/epitope_mapping")
import Bio.PDB.StructureBuilder
from Bio.PDB import PDBIO, PDBParser, PPBuilder, CaPPBuilder, Superimposer
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio import SeqIO, Align

import seq_align, metrics
import warnings

warnings.filterwarnings('ignore')


def cutoff_bfactor(structure, path, rank, deterministic, sum):
    """

    :param structure:
    :param path: saving path
    :param rank: Top ranking of residue likelihood score (e.g. result from AbEMap),
    None if using deterministic
    :param deterministic: cutoff value of a probability score in [0,1]
    None if using rank
    :return:
    """
    bfactor_combined = seq_align.get_sequence_bfactor(structure, sum)
    assert (rank is None and isinstance(deterministic, float)) or (
            isinstance(rank, int) or isinstance(rank, float) and deterministic is None)
    if isinstance(rank, int) or isinstance(rank, float):
        cutoff = metrics.find_top_cutoff(bfactor_combined, int(rank))
    elif isinstance(deterministic, float):
        cutoff = deterministic
    for model in structure:
        for chain in model:
            for resdiue in chain:
                for atom in resdiue:
                    if atom.get_bfactor() >= cutoff:
                        atom.set_bfactor(1.0)
                    else:
                        atom.set_bfactor(0.0)
    io = PDBIO()
    io.set_structure(structure)
    io.save(path)


def overlap_bfactor(structure1, structure2, path):
    """

    :param structure1:
    :param structure2:
    :param path: new bfactor saving path
    :return:
    """
    superimposer = Superimposer()
    structure1_residues = []
    ref_atoms = []
    sample_atoms = []

    for model1 in structure1:
        for chain1 in model1:
            for residue1 in chain1:
                try:
                    ref_atoms.append(residue1['CA'])
                except:
                    pass

    for model2 in structure2:
        for chain2 in model2:
            for residue2 in chain2:
                try:
                    sample_atoms.append(residue2['CA'])
                except:
                    pass
    print(len(ref_atoms), len(sample_atoms))

    # superimposer.set_atoms(ref_atoms, sample_atoms)
    # superimposer.apply(structure2[0].get_atoms())
    for model1, model2 in zip(structure1, structure2):
        for chain1, chain2 in zip(model1, model2):
            for residue1, residue2 in zip(chain1, chain2):
                for atom1, atom2 in zip(residue1, residue2):
                    if int(atom1.get_bfactor()) == int(atom2.get_bfactor() == 1):
                        atom1.set_bfactor(atom1.get_bfactor())
                    else:
                        atom1.set_bfactor(0.0)
    io = PDBIO()
    io.set_structure(structure1)
    io.save(path)
