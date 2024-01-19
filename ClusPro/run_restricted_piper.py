import os
import csv
import shutil
import subprocess


def read_csv_worksheet(csv_path):
    rec_lig = []
    with open(csv_path, 'r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)

        for row in csv_reader:
            if len(row) == 2:
                rec_lig.append(row[0][0:4])
    return rec_lig


csv_path = '/projectnb2/docking/imhaoyu/my_protein_toolbox/run_example/collect/rec_lig_name.csv'
work_list = read_csv_worksheet(csv_path)
workspace_base_path = '/projectnb2/docking/imhaoyu/ClusPro_transmembrane/piper-workspace'
for item in work_list[2:]:
    if item == '7TAK':
        work_path = os.path.join(workspace_base_path, item)
        shutil.copy('/projectnb2/docking/imhaoyu/my_protein_toolbox/ClusPro/cluster.sh', work_path)
        os.chdir(work_path)
        print(os.getcwd())
        os.system(
            '/projectnb2/docking/imhaoyu/piper-test/build/piper -vv -c1.0 -p /usr4/spclpgm/imhaoyu/mol-prms/atom/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec -f /usr4/spclpgm/imhaoyu/mol-prms/coeffs/coeffs.0.0.6.cs1.0.serverdemo -r /projectnb2/docking/imhaoyu/ClusPro_TM/rot_mat/rot70k.0.0.6.jm.max_10deg_z-axis_tilt.mol2 -k4 rec.pdb lig.pdb')
        os.system('chmod 744 cluster.sh')
        os.system('./cluster.sh')
        os.chdir('/projectnb2/docking/imhaoyu/my_protein_toolbox/ClusPro')
