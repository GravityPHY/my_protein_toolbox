import os
import json
import pathlib
import numpy as np
import pandas as pd
from pathlib import Path

from Bio.PDb import PDBIO
from Bio.PDB import MMCIFParser
from Bio.PDB import NeighborSearch, Selection


class AF3model:
    def __init__(self, output_path):
        self.model_path = output_path
        self.basedir = os.path.dirname(output_path)
        self.basename = os.path.basename(output_path)
        self.prefix = self.get_prefix(output_path)
        self.model_num = self.get_model_num(output_path)

    def get_file_name(self, path):
        file_name = Path(path).stem
        return file_name

    def get_prefix(self, path):
        file_name = self.get_file_name(path)
        prefix = "_".join(file_name.split("_")[0:-2])
        return prefix

    def get_suffix(self, path):
        suffix = pathlib.Path(path).suffix
        return suffix

    def get_model_num(self, path):
        filename = self.get_file_name(path)
        model_num_str = filename.split("_")[-1]
        return int(model_num_str)

    def get_score(self, choice):
        """

        :param choice: choose from ["ranking_score", "iptm", "ptm"]
        :return:
        """
        assert choice in ["ranking_score", "iptm", "ptm"]
        d=self.read_json(self.get_json_path("confidence"))
        return d[choice]

    def read_json(self, path):
        assert self.get_suffix(path) == '.json'
        with open(path) as file:
            d = json.load(file)
        return d

    def get_json_path(self, choice):
        """

        :param choice: choose from ["full_data","confidence","job_request"]
        :return:
        """
        assert choice in ["full_data", "confidence", "job_request"]
        if choice == 'full_data':
            midfix = "_full_data_"
        elif choice == "confidence":
            midfix = "_summary_confidences_"
        elif choice == "job_request":
            midfix = "_job_request_"

        model_dir = self.basedir
        prefix = self.prefix
        model_num = self.model_num
        suffix = ".json"
        filename = prefix + midfix + str(model_num) + suffix
        full_path = os.path.join(model_dir, filename)
        return full_path

    def select_interface_atoms(self, cif_path, part_1, part_2):
        parser = MMCIFParser()
        structure = parser.get_structure('complex', cif_path)

        chain_1 = [structure[0][name] for name in part_1]
        chain_2 = [structure[0][name] for name in part_2]

        atoms_1 = Selection.unfold_entities(chain_1, 'A')  # 'A' is for atoms
        atoms_2 = Selection.unfold_entities(chain_2, 'A')

        cutoff_distance = 5.0
        neighbor_search_2 = NeighborSearch(atoms_2)
        interface_atoms_2 = [atom for atom in atoms_1 if neighbor_search_2.search(atom.coord, cutoff_distance)]

        neighbor_search_1 = NeighborSearch(atoms_1)
        interface_atoms_1 = [atom for atom in atoms_2 if neighbor_search_1.search(atom.coord, cutoff_distance)]

        interface_atoms_ = interface_atoms_1 + interface_atoms_2
        return interface_atoms_

    def get_interface_bfactors(self, cif_path, part_1, part_2):
        interface_atoms=self.select_interface_atoms(cif_path,part_1,part_2)
        b_factors = [atom.get_bfactor() for atom in interface_atoms]
        average_bfactor = np.mean(b_factors)
        return average_bfactor

    def get_antigen_interface_atoms(self,cif_path,antigen_ch, antibody_ch):
        parser = MMCIFParser()
        structure = parser.get_structure('complex', cif_path)

        antigen = [structure[0][name] for name in antigen_ch]
        antibody = [structure[0][name] for name in antibody_ch]

        atoms_antigen = Selection.unfold_entities(antigen, 'A')  # 'A' is for atoms
        atoms_antibody = Selection.unfold_entities(antibody, 'A')

        cutoff_distance = 5.0
        neighbor_search = NeighborSearch(atoms_antibody)
        antigen_interface_atoms_ = [(atom.get_serial_number(),atom) for  atom in atoms_antigen if neighbor_search.search(atom.coord, cutoff_distance)]
        return antigen_interface_atoms_

    def cif_reader(self):
        raise NotImplemented




if __name__ == "__main__":
    model = AF3model("/projectnb2/docking/imhaoyu/my_protein_toolbox/for_AF3/AF3_unzip/7SGM/fold_7sgm_1_model_0.cif")
    print(model.prefix)
    print(model.get_score("ranking_score"))
    print(model.select_interface_atoms(model.model_path,"BC","A"))
    print(model.get_antigen_interface_atoms(model.model_path,"A","BC"))
