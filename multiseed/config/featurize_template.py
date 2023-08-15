import numpy as np
from Bio import Align
from Bio.Seq import Seq
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.PDBParser import PDBParser
from Bio.Align import substitution_matrices
from Bio.Data.SCOPData import protein_letters_3to1
import sys
sys.path.insert(0,'/projectnb2/sc3dm/imhaoyu/alphafold/')
from alphafold.common import protein, residue_constants
from alphafold.data import pipeline, templates

def read_pdb_structure(pdb_path):
    parser = PDBParser(PERMISSIVE=0)
    cif_object = parser.get_structure('', pdb_path)
    model = next(cif_object.get_models())
    return(model)


def seq_mapping(query_seq, template_seq):  # can be replaced by any sequence alignment tool
    match_list = []
    query_len = len(query_seq)
    template_len = len(template_seq)
    for shift in range(1 - template_len, query_len):
        query_start = max(0, shift)
        template_start = max(0, -shift)
        match_len = min(template_len - template_start, query_len - query_start)

        match = {'count': 0, 'mapping': {}}
        for i in range(match_len):
            query_i = query_start + i
            template_i = template_start + i
            if query_seq[query_i] == template_seq[template_i]:
                match['count'] += 1
                match['mapping'][query_i] = template_i
            elif template_seq[template_i] == 'X':  # unresolved or super non-standard residue
                match['count'] += 1
                # Maybe we should differentiate between unresolved residue "." and
                # super non-standard residue "X" (perhaps we can still extract some
                # template atomic coordinates from it, for example, backbone atoms)
            else:
                break

        if match['count'] == match_len:
            match_list.append(match)
    if match_list:
        max_count = max(match_list, key=lambda x: x['count'])['count']
        max_match_list = [x for x in match_list if x['count'] == max_count]
        if len(max_match_list) == 1:
            return max_match_list[0]['mapping']
        else:
            # TODO: raise exception & logging
            print(f'Error: Ambiguous sequence mapping of max size {max_count}!')
            print(f'Query: {query_seq}')
            print(f'Template: {template_seq}')
            print(max_match_list)
            exit(1)
    else:
        # TODO: raise exception & logging
        print(f'Error: no match between two sequences!')
        print(f'Query: {query_seq}')
        print(f'Template: {template_seq}')
        exit(1)


def alignment_seq_mapping(query_seq, template_seq):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignments = aligner.align(Seq(query_seq), Seq(template_seq))
    if len(alignments) != 1:
        # TODO: raise exception & logging
        print(f'Warning: More than one sequence alignment!')
        print(f'Query: {query_seq}')
        print(f'Template: {template_seq}')
        print(len(alignments))
        for alignment in alignments:
            print(alignment)
        print("The first one will be used for further processing")
    alignment = alignments[0]
    mapping = {}
    for qrange, trange in zip(*alignment.aligned):
        for qi, ti in zip(range(*qrange), range(*trange)):
            mapping[qi] = ti
    return mapping

def get_positions_masks_seq_for_chain(model, chain_id, atom_types='all'):
    chain = model[chain_id]

    positions = []
    masks = []
    seq = ''
    last_resi = None
    for residue in chain:
        if is_aa(residue):
            resi = residue.id[1]
            if last_resi and resi > last_resi + 1:
                gap_len = resi - last_resi - 1
                seq += 'X' * gap_len
                positions.append(np.zeros((gap_len, residue_constants.atom_type_num, 3), dtype=np.float64))
                masks.append(np.zeros((gap_len, residue_constants.atom_type_num), dtype=np.int64))
            resn = residue.get_resname()
            code = protein_letters_3to1[resn]
            seq += (code if len(code) == 1 else 'X')

            pos = np.zeros((1, residue_constants.atom_type_num, 3), dtype=np.float64)
            mask = np.zeros((1, residue_constants.atom_type_num), dtype=np.int64)
            for atom in residue:
                atom_name = atom.id  # get_name()
                if atom_types != 'all' and (atom_name not in atom_types):
                    continue
                coords = atom.get_coord()
                if atom_name in residue_constants.atom_order:
                    pos[0][residue_constants.atom_order[atom_name]] = coords
                    mask[0][residue_constants.atom_order[atom_name]] = 1
                elif atom_name.upper() == 'SE' and resn == 'MSE':  # upper?
                    # Put the coordinates of the selenium atom in the sulphur column.
                    pos[0][residue_constants.atom_order['SD']] = coords
                    mask[0][residue_constants.atom_order['SD']] = 1
            positions.append(pos)
            masks.append(mask)
            last_resi = resi
    positions = np.concatenate(positions) #, dtype=np.float64)
    masks = np.concatenate(masks) #, dtype=np.int64)
    templates._check_residue_distances(positions, masks, max_ca_ca_distance=150)
    return positions, masks, seq

def get_positions_masks_seqs(model, chain_ids):
    all_positions = []
    all_masks = []
    seqs = []
    for chain_id in chain_ids:
        positions, masks, seq = get_positions_masks_seq_for_chain(model, chain_id)
        all_positions.append(positions)
        all_masks.append(masks)
        seqs.append(seq)
    all_positions = np.concatenate(all_positions) #, dtype=np.float64)
    all_masks = np.concatenate(all_masks) #, dtype=np.int64)
    templates._check_residue_distances(all_positions, all_masks, max_ca_ca_distance=150)
    return all_positions, all_masks, seqs


def featurize_custom_template(
        seq_data_list,
        custom_template):
    template_model = read_pdb_structure(custom_template['pdb_file'])
    query_to_template = custom_template.get('query_to_template', None)
    # if query-to-template chain mapping non specified,
    # assume same chain names.
    if query_to_template == None:
        query_to_template = {}
        for seq_data in seq_data_list:
            for chain in seq_data['chains']:
                query_to_template[chain] = chain
    mapping_method = custom_template.get("mapping_method", None)
    if mapping_method == None:
        print("Mapping method not specified for template. Will use a default one.")
        mapping_method = "custom"
    atom_types = custom_template.get("atom_types", 'all')
    alaninize = custom_template.get("alaninize", False)
    #query_chain_ids = custom_template['query_to_template'].keys()
    #template_chain_ids = custom_template['query_to_template'].values()
    #all_positions, all_masks, template_seqs = get_positions_masks_seqs(template_model, template_chain_ids)

    per_chain_template_features = {}
    for seq_data in seq_data_list:
        for query_chain_id in seq_data['chains']:
            query_seq = seq_data['sequence']


            query_seq_len = len(query_seq)
            output_seq = ['-'] * query_seq_len
            output_confidence_scores = np.full(query_seq_len, -1)
            output_positions = np.zeros((query_seq_len, residue_constants.atom_type_num, 3), dtype=np.float64)
            output_masks = np.zeros((query_seq_len, residue_constants.atom_type_num), dtype=np.int64)

            template_chain_id = query_to_template.get(query_chain_id, None)
            if template_chain_id != None:
                positions, masks, template_seq = get_positions_masks_seq_for_chain(
                                                        template_model,
                                                        template_chain_id,
                                                        atom_types=atom_types)

                if mapping_method == 'custom':
                    mapping = seq_mapping(query_seq, template_seq)
                elif mapping_method == 'alignment':
                    mapping = alignment_seq_mapping(query_seq, template_seq)
                else:
                    print(f"Unknown mapping method specified: {mapping_method}")
                    exit()
                for q, t in mapping.items():
                    output_positions[q] = positions[t]
                    output_masks[q] = masks[t]
                    output_seq[q] = template_seq[t] if not alaninize else 'A'
                    output_confidence_scores[q] = 5
                output_seq = ''.join(output_seq)
                templates_aatype = residue_constants.sequence_to_onehot(output_seq, residue_constants.HHBLITS_AA_TO_ID)
                template_features = {'template_all_atom_positions': output_positions,
                                     'template_all_atom_masks': output_masks,
                                     'template_sequence': output_seq.encode(),
                                     'template_aatype': np.array(templates_aatype),
                                     'template_confidence_scores': output_confidence_scores,
                                     'template_domain_names': f'none'.encode(),
                                     'template_release_date': f'none'.encode()}
                for feature in template_features:
                    template_features[feature] = np.expand_dims(template_features[feature], axis=0)
                    per_chain_template_features[query_chain_id] = template_features
    return per_chain_template_features

