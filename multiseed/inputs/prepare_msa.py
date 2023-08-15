import os
import argparse
import subprocess
#python msa_mmseq.py -fa [directory where your fasta is]#

parser=argparse.ArgumentParser(description="Give the directory path where your fasta is")
parser.add_argument('-fa','--fasta_dir')
args=parser.parse_args()
fasta_dir=args.fasta_dir

#fasta_dir='/projectnb2/docking/imhaoyu/AFM-tests/example/input_files/4DKL/'
command=f'python msa_mmseq.py -fas {fasta_dir} -out {fasta_dir}'.split()
subprocess.call(command)
