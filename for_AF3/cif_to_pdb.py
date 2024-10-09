import os
import json
import glob
import argparse

from pathlib import Path
from Bio import PDB
from Bio.PDB import MMCIFParser, PDBIO


def convert_cif_to_pdb(directory):
    parser = MMCIFParser(QUIET=True)
    io = PDBIO()
    for filename in os.listdir(directory):
        if filename.endswith(".cif"):
            cif_file=os.path.join(directory,filename)
            structure_id = os.path.splitext(filename)[0]
            structure = parser.get_structure(structure_id, cif_file)
            pdb_file = os.path.splitext(cif_file)[0]+".pdb"
            io.set_structure(structure)
            io.save(pdb_file)

parser = argparse.ArgumentParser(description="Script that convert all .cif in a given directory to .pdb")

parser.add_argument("--path", required=True, help="Path to a directory")

args = parser.parse_args()
dir_path=args.path
convert_cif_to_pdb(dir_path)