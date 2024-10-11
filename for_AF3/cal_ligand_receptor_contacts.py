import os
import pandas as pd
from pathlib import Path

import pymol
from pymol import cmd

from Bio.PDB import MMCIFParser, PDBIO, PDBParser


def analyze_contacts(structure_files,
                     receptor_chains="R",
                     ligand_chains="L",
                     d_cutoff=5.0):
    """

    :param structure_files:
    :param receptor_chains: string of names of receptor chains
    :param ligand_chains: string of names of ligand chains
    :param d_cutoff:
    :return:
    """
    contact_count = {}
    total_structures = len(structure_files)
    for structure in structure_files:
        cmd.load(structure,"complex")
        structure_name=Path(structure).stem
        selection_name = f"ligand_contacts_{structure_name}"
        cmd.select(selection_name,f"(br. (chain {'+'.join(ligand_chains)})) within {d_cutoff} of chain {'+'.join(receptor_chains)}")
        contacts = cmd.get_model(selection_name)
        for atom in contacts.atom:
            residue_id = (atom.chain, atom.resn, atom.resi)  # Tuple of residue name and ID
            if residue_id not in contact_count:
                contact_count[residue_id] = 0
            contact_count[residue_id] += 1

        cmd.delete("all")  # Clear the structure from PyMOL
    contact_frequency = {residue: count / total_structures for residue, count in contact_count.items()}
    return contact_frequency


def clean_bfactors(pdb_file, output_file=None, default_bfactor=0.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.set_bfactor(default_bfactor)
    return structure


def update_bfactor(pdb_file,contact_frequency,ligand_chains="L"):
    #parser = PDBParser(QUIET=True)
    structure = clean_bfactors(pdb_file)
    for model in structure:
        for chain in model:
            if chain.id in ligand_chains:
                for residue in chain:
                    res_id = (chain.id, residue.get_resname(), str(residue.get_id()[1]))
                    if res_id in contact_frequency:
                        print("yes")
                        frequency = contact_frequency[res_id]
                        for atom in residue:
                            atom.set_bfactor(frequency)
    io = PDBIO()
    io.set_structure(structure)
    updated_pdb_file=os.path.splitext(pdb_file)[0]+"_bfactor.pdb"
    io.save(updated_pdb_file)

pdb_dir="/projectnb2/docking/imhaoyu/24_epitope_mapping/database/ABAG/AlphaFold3/multiseed/pdb/7SU1"
structure_files=os.listdir(pdb_dir)
structure_paths=[os.path.join(pdb_dir,i) for i in structure_files]
contact_frequency = analyze_contacts(structure_paths, receptor_chains='AB', ligand_chains='C')
print(contact_frequency)
update_bfactor("/projectnb2/docking/imhaoyu/24_epitope_mapping/database/ABAG/AlphaFold3/multiseed/mapping_test/fold_7su1_1_model_0.pdb",
               contact_frequency,
               ligand_chains="C")