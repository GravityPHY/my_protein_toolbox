import os
import csv
import glob
import json
import shutil


def get_txt_address(dir_save_prediction_files, pdb_id):
    """

    :param path_save_prediction_files:
    :param pdb_id:
    :return: list of addresses
    """
    assert os.path.isdir(dir_save_prediction_files) == True
    pdb_dir = os.path.join(dir_save_prediction_files, pdb_id)
    assert os.path.isdir(pdb_dir) == True
    txt_files = glob.glob(os.path.join(pdb_dir, '*.txt'))
    return txt_files


def get_pdb_address(dir_save_prediction_files, pdb_id):
    """

    :param dir_save_prediction_files:
    :param pdb_id:
    :return: list of addresses
    """
    assert os.path.isdir(dir_save_prediction_files) == True
    pdb_dir = os.path.join(dir_save_prediction_files, pdb_id)
    assert os.path.isdir(pdb_dir) == True
    pdb_files = glob.glob(os.path.join(pdb_dir, '*.pdb'))
    return pdb_files


def get_match_txt_pdb(dir_save_prediction_files, pdb_id, model_num):
    """

    :param dir_save_prediction_files:
    :param pdb_id:
    :param model_num: can be an integer in [1,2,3,4,5]
    :return:
    """
    assert os.path.isdir(dir_save_prediction_files) == True
    pdb_dir = os.path.join(dir_save_prediction_files, pdb_id)

    assert os.path.isdir(pdb_dir) == True
    if isinstance(model_num, int):
        model_num = str(model_num)

    txt_file = glob.glob(os.path.join(pdb_dir, f'*model_{model_num}*.txt'))
    pdb_file = glob.glob(os.path.join(pdb_dir, f'*model_{model_num}*.pdb'))
    return txt_file[0], pdb_file[0]


def txt_to_json(txt_path):
    assert os.path.exists(txt_path)
    with open(txt_path) as f:
        data = f.read()
    js = json.loads(data)
    return js


def get_highest_confidence_pdb(dir_save_prediction_files, pdb_id, total_model=5, key="ranking_confidence"):
    """

    :param dir_save_prediction_files:
    :param pdb_id:
    :return:
    """
    max_confidence = -float('inf')
    max_conf_pdb = None
    assert isinstance(total_model, int)
    for i in range(1, total_model + 1):
        txt_path, pdb_path = get_match_txt_pdb(dir_save_prediction_files, pdb_id, i)
        confidence_dict = txt_to_json(txt_path)
        score = confidence_dict[key]
        if score > max_confidence:
            max_confidence = score
            max_conf_pdb = pdb_path
    return max_conf_pdb, max_confidence


def copy_pdb(pdb_path, destination_dir, destination_filename):
    """

    :param pdb_path:
    :param destination_dir:
    :param destination_filename:
    :return:
    """
    assert os.path.isdir(destination_dir)
    dst_path = os.path.join(destination_dir, destination_filename)
    os.makedirs(os.path.dirname(dst_path), exist_ok=True)
    shutil.copyfile(pdb_path, dst_path)
    print(f'Finish copying {pdb_path} to {dst_path}')


def prepare_piper_input(pdb_id,
                        AFM_models_dir,
                        piper_workspace_path,
                        piper_input_type='rec'):
    max_conf_pdb, max_confidence = get_highest_confidence_pdb(AFM_models_dir,
                                                              pdb_id, )
    copy_pdb(max_conf_pdb, piper_workspace_path, f'{pdb_id[0:4]}/{piper_input_type}.pdb')


def read_csv_worksheet(csv_path):
    rec_lig = []
    with open(csv_path, 'r') as file:
        csv_reader = csv.reader(file)
        next(csv_reader)

        for row in csv_reader:
            if len(row) == 2:
                rec_lig.append(('rec', row[0]))
                rec_lig.append(('lig', row[1]))
    return rec_lig


if __name__ == '__main__':
    csv_path = '/projectnb2/docking/imhaoyu/my_protein_toolbox/run_example/collect/rec_lig_name.csv'
    work_list=read_csv_worksheet(csv_path)
    for item in work_list:
        piper_input_type, pdb_id= item
        prepare_piper_input(pdb_id,
                            '/projectnb2/docking/imhaoyu/transmebrane/_data_ALLDIMERTMproteins/AF2monomer/output_models',
                            '/projectnb2/docking/imhaoyu/ClusPro_transmembrane/piper-workspace',
                            piper_input_type=piper_input_type)
    # prepare_piper_input('1KPL_1',
    #                    '/projectnb2/docking/imhaoyu/transmebrane/_data_ALLDIMERTMproteins/AF2monomer/output_models',
    #                    '/projectnb2/docking/imhaoyu/ClusPro_transmembrane/piper-workspace')
