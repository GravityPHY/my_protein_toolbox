"""
Prepare input amino acid sequence by providing
 - PDB id or
 - structure in PDB format (TODO)
"""
import os
import csv
import requests as re
import pandas as pd
from Bio import PDB
from Bio.PDB import PPBuilder


def fasta_from_rcsb(pdb_id, output_dir):
    rcsb_fasta_url = f'https://www.rcsb.org/fasta/entry/{pdb_id}/display'
    try:
        response = re.get(rcsb_fasta_url)
        response.raise_for_status()
        file_dir = os.path.join(output_dir, pdb_id)
        file_name = f'{pdb_id}.fa'
        file_path = f"{file_dir}/{file_name}"
        with open(file_path, 'wb') as file:
            file.write(response.content)
        print(f"Fasta file has been save to {file_path}")
    except re.exceptions.RequestException as e:
        print(f"{e}")
