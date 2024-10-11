import os
import json
import pathlib
from pathlib import Path


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
        prefix = "_".join(file_name.split("_")[0:3])
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

    def cif_reader(self):
        raise NotImplemented




if __name__ == "__main__":
    model = AF3model("/projectnb2/docking/imhaoyu/my_protein_toolbox/for_AF3/AF3_unzip/7SGM/fold_7sgm_1_model_0.cif")
    print(model.prefix)
    print(model.get_score("ranking_score"))
