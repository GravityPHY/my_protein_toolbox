import os
import argparse
import subprocess
pdb_id=['6Z70']#['6Z70','7A0W','7DRJ','7M93','7SEL','7U55','7VR5','8A6E']
for pdb_id in pdb_id:
    directory = f'/projectnb2/docking/imhaoyu/ClusPro_TM/showcases/{pdb_id}/ClusProPDB'
    if not os.path.exists(directory):
        os.makedirs(directory)
    if not os.path.exists(f'/projectnb2/docking/imhaoyu/ClusPro_TM/showcases/{pdb_id}/clusterPDB'):
        os.makedirs(f'/projectnb2/docking/imhaoyu/ClusPro_TM/showcases/{pdb_id}/clusterPDB')
    current_path=directory
    model_list = os.listdir(current_path)

    for model in model_list:
        input_file_path = os.path.join(current_path, model)
        output_file_path = os.path.join(f'/projectnb2/docking/imhaoyu/ClusPro_TM/showcases/{pdb_id}/clusterPDB', model)
        #os.makedirs(output_dir, exist_ok=True)
        #output_file_path = os.path.join(output_dir, model)
        print(output_file_path)
        with open(input_file_path, "r") as file, open(output_file_path, 'w+') as out:
            lines = file.readlines()
            total_lines = len(lines)
            switch_flag = False
            #PDBID = name[-7::]
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