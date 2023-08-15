#!/bin/sh
# Join stout and sterr to a single file
#$ -j y
# Output log files
#$ -o :logs/
#$ -e :logs/
#$ -P docking
#$ -N afm_no_temp
# Execute the job from the current working directory.
#$ -cwd
# Do not email.
#$ -m n
#$ -V
#$ -l h_rt=16:00:00
#$ -l gpus=1
#$ -l gpu_c=3.5
#$ -t 1-1

module load miniconda
module load cuda/11.1
conda activate alpha_mul_v3
cd /projectnb2/docking/imhaoyu/my_protein_toolbox/multiseed/config

echo $PDB_chain
#PDB_chain=3OMI_A
python3 af_multimer_multiseed.py -output_dir /projectnb2/docking/imhaoyu/my_protein_toolbox/example/{PDB_chain}/${PDB_chain}.fa --pre_align_dir /projectnb2/docking/imhaoyu/my_protein_toolbox/example/${PDB_chain}/${PDB_chain} --num_prediction 10 --model_preset multimer_v2 --random_seed 10 --enable_single_seed=False
