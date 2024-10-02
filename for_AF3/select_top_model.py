import os
import json
import glob
import argparse
from pathlib import Path
from Bio import PDB

parser = argparse.ArgumentParser(description="Script to select the Top ranking model from AF3 server prediction")

parser.add_argument("--path_to_zip", required=True, help="Path to .zip")
parser.add_argument("--path_to_unzip", required=True, help="Path to folder to save the unzip files")
parser.add_argument("--pdb_save_path", required=True, help="Path to save the Top ranking model in PDB format")


args = parser.parse_args()

path_to_zip = args.path_to_zip #"/projectnb2/docking/imhaoyu/my_protein_toolbox/for_AF3/8WSQ_antibody.zip"
path_to_unzip = args.path_to_unzip #"/projectnb2/docking/imhaoyu/my_protein_toolbox/for_AF3/8WSQ_test"
os.system(f"unzip {path_to_zip} -d {path_to_unzip}")
job_request_json_path = glob.glob(f"{path_to_unzip}/*_job_request.json")[0]
# find *_job_request.json, take prefix of *
prefix = Path(job_request_json_path).stem.strip("_job_request")
print(prefix)
# read every prefix_summary_confidences_*.json, find the one with the highest ranking_score, ptm
path_to_summary_confidences = glob.glob(f"{path_to_unzip}/{prefix}_summary_confidences_*.json")
models_summary = {}
for path in path_to_summary_confidences:
    model_num = Path(path).stem.split("_")[-1]
    model_name = f"{prefix}_model_{model_num}.cif"
    content = json.load(open(path))
    models_summary[model_name] = [content['ranking_score'],
                                  content['ptm'],
                                  content['iptm']]
sorted_models_summary = sorted(models_summary.items(), key=lambda item: item[1], reverse=True)

# select top ranking model and save as pdb
top_model_path = os.path.join(path_to_unzip, sorted_models_summary[0][0])

cif_parser = PDB.MMCIFParser()
pdb_io = PDB.PDBIO()

structure = cif_parser.get_structure("structure", top_model_path)
pdb_io.set_structure(structure)
pdb_save_path = args.pdb_save_path
pdb_io.save(pdb_save_path)
# print(top_model_path)
#
# print(content)
# "ranking_score": 0.9
# "ptm": 0.85,

# unzip
# .zip path
# unzip folder name and save dir
# parse confidence information in json
# find model with highest confidence and save it
# save path
# save name
