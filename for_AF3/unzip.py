import os
import json
import glob
import argparse

from pathlib import Path
from Bio import PDB
from Bio.PDB import MMCIFParser, PDBIO


def unzip(zip_dir, save_dir):
    paths=glob.glob(os.path.join(zip_dir,"*.zip"))
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    for path_to_zip in paths:
        os.system(f"unzip {path_to_zip} -d {save_dir}")


parser = argparse.ArgumentParser(description="Script unzip all .zip files to another folder")
parser.add_argument("--zip_dir", required=True, help="Path to a directory contains .zip")
parser.add_argument("--save_dir", required=True, help="Path to a saving directory")
args = parser.parse_args()
zip_dir=args.zip_dir
save_dir=args.save_dir
unzip(zip_dir,save_dir)

