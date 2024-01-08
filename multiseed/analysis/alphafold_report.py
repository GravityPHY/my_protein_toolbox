
"""Need some tests"""
import os
import csv
import json
import glob
import subprocess
import argparse
from abc import ABC, abstractmethod

class Report(ABC):
    """
    Template for AlphaFold output report
    """

    def get_txt_list(self, read_dir):
        pass

    @abstractmethod
    def txt_to_dict(self, txt_path):
        pass

    @abstractmethod
    def get_PDBname(self, txtname):
        pass

    def write_csv(self, name, sort:bool):
        pass


class Multimer(Report):
    def __init__(self, read_dir, write_dir):
        super(Multimer, self).__init__()
        self.read_dir = read_dir
        self.write_dir = write_dir

    def get_txt_list(self, read_dir):
        """
        The confidence score of each sampling predictions is in a .txt file
        Get a list of .txt of all predictions in the directory
        """
        txt_list = glob.glob(read_dir + '/*.txt')
        return txt_list

    def txt_to_dict(self, txt_path):
        """
        There is a dictionary in the .txt saving the confidence score,
        this function load this dictionary
        :param txt_path:
        :return: dict
        """
        file = open(txt_path)
        d = json.load(file)
        return d

    def get_PDBname(self, txtname):
        """
        Each confidence score match a structure saved in .pdb, save under the same directory
        The name of .pdb and .txt is only different at beginning and ending
        :param txtname: str
        :return: pdbname: str
        """
        words = ['confidence', 'txt']
        replace = ['unrelaxed', 'pdb']

        for word, rep in zip(words, replace):
            txtname = txtname.replace(word, rep)
        return txtname

    def write_csv(self, name, sort=True, *args, **kwargs):
        """

        :param name: the name of the .csv file
        :param sort: boolean, if True write the file of confidence score in decreasing order
        :return: None
        """
        txt_list = self.get_txt_list(self.read_dir)
        dict0 = json.load(open(txt_list[0]))
        dict0['model'] = os.path.basename(txt_list[0])
        os.makedirs(self.write_dir, exist_ok=True)
        with open(os.path.join(self.write_dir, name), 'w') as f:
            w = csv.DictWriter(f, dict0.keys())
            w.writeheader()
            for txt in txt_list:
                _dict = json.load(open(txt))
                model_add = os.path.join(self.read_dir, self.get_PDBname(os.path.basename(txt)))
                _dict['model'] = model_add
                w.writerow(_dict)
        if sort:
            with open(os.path.join(self.write_dir, name), 'r') as input_file:
                csv_reader = csv.DictReader(input_file)
                data = list(csv_reader)
            sorted_data = sorted(data, key=lambda row: row['ranking_confidence'], reverse=True)

            with open(os.path.join(self.write_dir, name), 'w') as output_file:
                csv_writer = csv.DictWriter(output_file, dict0.keys())
                csv_writer.writeheader()
                for row in sorted_data:
                    csv_writer.writerow(row)







