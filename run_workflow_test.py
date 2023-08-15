import os
import sys
import subprocess
sys.path.insert(0,'/projectnb2/docking/imhaoyu/my_protein_toolbox/multiseed/inputs')
from get_fasta import fasta_from_rcsb
pdb_id='4DKL'
fasta_from_rcsb(pdb_id, '/projectnb2/docking/imhaoyu/my_protein_toolbox/example')
fasta_dir=f'/projectnb2/docking/imhaoyu/my_protein_toolbox/example/{pdb_id}'
command=f'python /projectnb2/docking/imhaoyu/my_protein_toolbox/multiseed/inputs/msa_mmseq.py -fas {fasta_dir} -out {fasta_dir}'.split()
subprocess.call(command)
