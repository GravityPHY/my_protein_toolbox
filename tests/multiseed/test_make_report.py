import os
import sys
import subprocess
sys.path.insert(0,'/projectnb2/docking/imhaoyu/my_protein_toolbox/multiseed')
import analysis.alphafold_report as report
pdb_list=['1E12','1EHK','1M0L']
version=['v1','v2','v3']
for pdb in pdb_list:
    for v in version:

        reporter=report.Multimer(f'/projectnb2/docking/imhaoyu/transmebrane/dataset1/trimer/{v}/{pdb}',
                         '/projectnb2/docking/imhaoyu/my_protein_toolbox/tests/multiseed/reports')
        reporter.write_csv(name=f'{v}-{pdb}.csv')