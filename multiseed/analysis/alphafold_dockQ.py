import os
import csv
import argparse
import subprocess


class DockQ:
    def __init__(self):
        self.dockq_path = '/projectnb2/docking/imhaoyu/DockQ/DockQ.py'

    def set_DockQ_path(self, path):
        self.dockq_path = path

    def get_DockQ_path(self):
        return self.dockq_path

    def set_native(self, native_path):
        self.native_path = native_path

    def get_native(self):
        return self.native_path

    def run_DockQ(self,model_path):
        output = subprocess.run([self.dockq_path, '-short', model_path, self.native_path],
                                capture_output=True, text=True)
        return output

    def write_DockQ(self,save_path, model_paths):
        if isinstance(model_paths,str):
            model_paths=list(model_paths)
        assert isinstance(model_paths,list)

        fieldnames = ['DockQ', 'Fnat', 'iRMS', 'LRMS', 'Model', 'Native']
        writer = csv.DictWriter(save_path, fieldnames=fieldnames)
        for path in model_paths:
            DockQ_output=self.run_DockQ(path)
            cml_output = DockQ_output.stdout.split()
            try:
                record_dict = {'DockQ': cml_output[1], 'Fnat': cml_output[3], 'iRMS': cml_output[5], 'LRMS': cml_output[7],
                           'Model': cml_output[10], 'Native': cml_output[11]}
                writer.writerow(record_dict)
            except IndexError:
                print(path)


