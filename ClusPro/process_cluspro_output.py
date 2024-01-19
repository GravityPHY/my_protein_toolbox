import os
import glob
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('--base_path', help='the directory that save the ClusPro outputs model.000.00*.pdb')
parser.add_argument('--save_path', help='the directory that save the ClusPro models after combining rec and lig in one .pdb file')
args = parser.parse_args()
if not args.save_path:
    args.save_path=os.path.join(args.base_path,'cluspro_models')
model_list=glob.glob(os.path.join(args.base_path,'model.00*.0**.pdb'))

for model in model_list:
    cluspro_pdb_path=model
    save_pdb_name=os.path.basename(cluspro_pdb_path)
    save_pdb_path=os.path.join(args.save_path,save_pdb_name)
    with open(cluspro_pdb_path, "r") as file, open(save_pdb_path, 'w+') as out:
        lines = file.readlines()
        total_lines = len(lines)
        switch_flag = False
        out.write(f"HEADER {model}\n")
        for i, line in enumerate(lines):
            if line.startswith("HEADER") and 'lig' in line:
                switch_flag = True
                continue

            if line.startswith("ATOM"):
                if not switch_flag:
                    s = list(line)
                    s[21] = 'A'
                    out.write("".join(s))
                    # row_data = line.split()
                    # row_data[4] = 'A'
                    # print(row_data)
                    # out.write("\t".join(row_data) + '\n')
                elif switch_flag:
                    s = list(line)
                    s[21] = 'B'
                    out.write("".join(s))

        out.write("END")
#args.save_path

