import collections
import contextlib
import copy
import dataclasses
import json
import os
import tempfile
from typing import Mapping, MutableMapping, Sequence, Optional, Any, Dict, Union, List
import numpy as np
import featurize_template as ft
import pathlib
import pickle
import random
import shutil
import sys
import time
from absl import logging
from absl import app
from absl import flags
from alphafold.common import protein
from alphafold.common import residue_constants
from alphafold.data import feature_processing
from alphafold.data import msa_pairing
from alphafold.data import parsers
from alphafold.data import pipeline
from alphafold.data import templates
from alphafold.data.tools import hhsearch
from alphafold.data.tools import hmmsearch
from alphafold.data import pipeline_multimer
#from alphafold.model import config
from alphafold.model import data
from alphafold.model import model
from alphafold.data import parsers
from alphafold.data import msa_identifiers
import codecs
import config as config

flags.DEFINE_enum('model_preset', 'monomer',
                  ['monomer', 'monomer_casp14', 'monomer_ptm', 'multimer_v1', 'multimer_v2', 'multimer_v3', 'multimer'],
                  'Choose preset model configuration - the monomer model, '
                  'the monomer model with extra ensembling, monomer model with '
                  'pTM head, or multimer model')
flags.DEFINE_integer('random_seed', None, 'The random seed for the data '
                     'pipeline. By default, this is randomly generated. Note '
                     'that even if this is set, Alphafold may still not be '
                     'deterministic, because processes like GPU inference are '
                     'nondeterministic.')
flags.DEFINE_string('output_dir', None, 'Path to a directory that will '
                    'store the results.')
flags.DEFINE_string('data_dir', None, 'Path to directory of supporting data.')

flags.DEFINE_string('fasta_dir', None, 'Path to fasta including all chain sequences')

flags.DEFINE_string('pre_align_dir', None, 'Path to prealignment including hhr, a3m')

#flags.DEFINE_string('id', None, 'id of design test')
flags.DEFINE_integer('num_prediction', 5, 'Number of runs per model parameter')

flags.DEFINE_integer('num_recycle', 3, 'Number of recycles')

flags.DEFINE_float('score_threshold', None, 'Threshold to exit the run')

flags.DEFINE_bool('use_mock_temp', False, 'Choose to use mock templates or not')

flags.DEFINE_bool('use_float16', True, 'Choose to use float16 or not')

flags.DEFINE_bool('use_pdb_temp', False, 'Choose to use pdb templates or not')

flags.DEFINE_string('pdb_template_path', None, 'Path to pdb template')

flags.DEFINE_bool('enable_single_seed', True, 'Choose single seed mode')

flags.DEFINE_bool('allow_msa_duplicate', True, 'Allow duplication in msa sequences')

flags.DEFINE_integer('single_seed', 100, 'Single seed for all predictions')

FLAGS = flags.FLAGS


def _check_flag(flag_name: str,
                other_flag_name: str,
                should_be_set: bool):
    if should_be_set != bool(FLAGS[flag_name].value):
        verb = 'be' if should_be_set else 'not be'
        raise ValueError(f'{flag_name} must {verb} set when running with '
                     f'"--{other_flag_name}={FLAGS[other_flag_name].value}".')

@dataclasses.dataclass(frozen=True)
class _FastaChain:
  sequence: str
  description: str


def _make_chain_id_map(*,
                       sequences: Sequence[str],
                       descriptions: Sequence[str],
                       ) -> Mapping[str, _FastaChain]:
  """Makes a mapping from PDB-format chain ID to sequence and description."""
  if len(sequences) != len(descriptions):
    raise ValueError('sequences and descriptions must have equal length. '
                     f'Got {len(sequences)} != {len(descriptions)}.')
  if len(sequences) > protein.PDB_MAX_CHAINS:
    raise ValueError('Cannot process more chains than the PDB format supports. '
                     f'Got {len(sequences)} chains.')
  chain_id_map = {}
  for chain_id, sequence, description in zip(
      protein.PDB_CHAIN_IDS, sequences, descriptions):
    chain_id_map[chain_id] = _FastaChain(
        sequence=sequence, description=description)
  return chain_id_map

@contextlib.contextmanager
def temp_fasta_file(fasta_str: str):
  with tempfile.NamedTemporaryFile('w', suffix='.fasta') as fasta_file:
    fasta_file.write(fasta_str)
    fasta_file.seek(0)
    yield fasta_file.name


def convert_monomer_features(
    monomer_features: pipeline.FeatureDict,
    chain_id: str) -> pipeline.FeatureDict:
  """Reshapes and modifies monomer features for multimer models."""
  converted = {}
  converted['auth_chain_id'] = np.asarray(chain_id, dtype=np.object_)
  unnecessary_leading_dim_feats = {
      'sequence', 'domain_name', 'num_alignments', 'seq_length'}
  for feature_name, feature in monomer_features.items():
    if feature_name in unnecessary_leading_dim_feats:
      # asarray ensures it's a np.ndarray.
      feature = np.asarray(feature[0], dtype=feature.dtype)
    elif feature_name == 'aatype':
      # The multimer model performs the one-hot operation itself.
      feature = np.argmax(feature, axis=-1).astype(np.int32)
    elif feature_name == 'template_aatype':
      feature = np.argmax(feature, axis=-1).astype(np.int32)
      new_order_list = residue_constants.MAP_HHBLITS_AATYPE_TO_OUR_AATYPE
      feature = np.take(new_order_list, feature.astype(np.int32), axis=0)
    elif feature_name == 'template_all_atom_masks':
      feature_name = 'template_all_atom_mask'
    converted[feature_name] = feature
  return converted

def int_id_to_str_id(num: int) -> str:
  """Encodes a number as a string, using reverse spreadsheet style naming.

  Args:
    num: A positive integer.

  Returns:
    A string that encodes the positive integer using reverse spreadsheet style,
    naming e.g. 1 = A, 2 = B, ..., 27 = AA, 28 = BA, 29 = CA, ... This is the
    usual way to encode chain IDs in mmCIF files.
  """
  if num <= 0:
    raise ValueError(f'Only positive integers allowed, got {num}.')

  num = num - 1  # 1-based indexing.
  output = []
  while num >= 0:
    output.append(chr(num % 26 + ord('A')))
    num = num // 26 - 1
  return ''.join(output)

def add_assembly_features(
    all_chain_features: MutableMapping[str, pipeline.FeatureDict],
    ) -> MutableMapping[str, pipeline.FeatureDict]:
  """Add features to distinguish between chains.

  Args:
    all_chain_features: A dictionary which maps chain_id to a dictionary of
      features for each chain.

  Returns:
    all_chain_features: A dictionary which maps strings of the form
      `<seq_id>_<sym_id>` to the corresponding chain features. E.g. two
      chains from a homodimer would have keys A_1 and A_2. Two chains from a
      heterodimer would have keys A_1 and B_1.
  """
  # Group the chains by sequence
  seq_to_entity_id = {}
  grouped_chains = collections.defaultdict(list)
  for chain_id, chain_features in all_chain_features.items():
    seq = str(chain_features['sequence'])
    if seq not in seq_to_entity_id:
      seq_to_entity_id[seq] = len(seq_to_entity_id) + 1
    grouped_chains[seq_to_entity_id[seq]].append((chain_features, protein.PDB_CHAIN_IDS.index(chain_id)))

  new_all_chain_features = {}
  chain_id = 1
  for entity_id, group_chain_features in grouped_chains.items():
    for sym_id, (chain_features, chain_id) in enumerate(group_chain_features, start=1):
      new_all_chain_features[
          f'{int_id_to_str_id(entity_id)}_{sym_id}'] = chain_features
      seq_length = chain_features['seq_length']
      chain_features['asym_id'] = chain_id * np.ones(seq_length)
      chain_features['sym_id'] = sym_id * np.ones(seq_length)
      chain_features['entity_id'] = entity_id * np.ones(seq_length)
      chain_id += 1
  return new_all_chain_features

def pad_msa(np_example, min_num_seq):
  np_example = dict(np_example)
  num_seq = np_example['msa'].shape[0]
  if num_seq < min_num_seq:
    for feat in ('msa', 'deletion_matrix', 'bert_mask', 'msa_mask'):
      np_example[feat] = np.pad(
          np_example[feat], ((0, min_num_seq - num_seq), (0, 0)))
    np_example['cluster_bias_mask'] = np.pad(
        np_example['cluster_bias_mask'], ((0, min_num_seq - num_seq),))
  return np_example

def run_msa_tool( msa_out_path: str,
                 msa_format: str,
                 max_sto_sequences: Optional[int] = None
                 ) -> Mapping[str, Any]:
    logging.warning('Reading MSA from file %s', msa_out_path)
    if msa_format == 'sto' and max_sto_sequences is not None:
      precomputed_msa = parsers.truncate_stockholm_msa(
          msa_out_path, max_sto_sequences)
      result = {'sto': precomputed_msa}
    else:
      with open(msa_out_path, 'r') as f:
        result = {msa_format: f.read()}
    return result

def mk_mock_template(
    query_sequence: Union[List[str], str], num_temp: int = 1
) -> Dict[str, Any]:
    ln = (
        len(query_sequence)
        if isinstance(query_sequence, str)
        else sum(len(s) for s in query_sequence)
    )
    output_templates_sequence = "A" * ln
    output_confidence_scores = np.full(ln, 1.0)

    templates_all_atom_positions = np.zeros(
        (ln, templates.residue_constants.atom_type_num, 3)
    )
    templates_all_atom_masks = np.zeros((ln, templates.residue_constants.atom_type_num))
    templates_aatype = templates.residue_constants.sequence_to_onehot(
        output_templates_sequence, templates.residue_constants.HHBLITS_AA_TO_ID
    )
    template_features = {
        "template_all_atom_positions": np.tile(
            templates_all_atom_positions[None], [num_temp, 1, 1, 1]
        ),
        "template_all_atom_masks": np.tile(
            templates_all_atom_masks[None], [num_temp, 1, 1]
        ),
        "template_sequence": [f"none".encode()] * num_temp,
        "template_aatype": np.tile(np.array(templates_aatype)[None], [num_temp, 1, 1]),
        "template_confidence_scores": np.tile(
            output_confidence_scores[None], [num_temp, 1]
        ),
        "template_domain_names": [f"none".encode()] * num_temp,
        "template_release_date": [f"none".encode()] * num_temp,
        "template_sum_probs": np.zeros([num_temp], dtype=np.float32),
    }
    return template_features

def make_msa_features_colab(msas: Sequence[parsers.Msa]) -> pipeline.FeatureDict:
  """Constructs a feature dict of MSA features."""
  if not msas:
    raise ValueError('At least one MSA must be provided.')

  int_msa = []
  deletion_matrix = []
  species_ids = []
  for msa_index, msa in enumerate(msas):
    if not msa:
      raise ValueError(f'MSA {msa_index} must contain at least one sequence.')
    for sequence_index, sequence in enumerate(msa.sequences):
      int_msa.append(
          [residue_constants.HHBLITS_AA_TO_ID[res] for res in sequence])
      deletion_matrix.append(msa.deletion_matrix[sequence_index])
      identifiers = msa_identifiers.get_identifiers(
          msa.descriptions[sequence_index])
      species_ids.append(identifiers.species_id.encode('utf-8'))

  num_res = len(msas[0].sequences[0])
  num_alignments = len(int_msa)
  features = {}
  features['deletion_matrix_int'] = np.array(deletion_matrix, dtype=np.int32)
  features['msa'] = np.array(int_msa, dtype=np.int32)
  features['num_alignments'] = np.array(
      [num_alignments] * num_res, dtype=np.int32)
  features['msa_species_identifiers'] = np.array(species_ids, dtype=np.object_)
  return features

def pair_and_merge_colab(
    all_chain_features: MutableMapping[str, pipeline.FeatureDict]
    ) -> pipeline.FeatureDict:
  """Runs processing on features to augment, pair and merge.

  Args:
    all_chain_features: A MutableMap of dictionaries of features for each chain.

  Returns:
    A dictionary of features.
  """

  feature_processing.process_unmerged_features(all_chain_features)

  np_chains_list = list(all_chain_features.values())

  pair_msa_sequences = not feature_processing._is_homomer_or_monomer(np_chains_list)

  if pair_msa_sequences:
    np_chains_list = msa_pairing.create_paired_features(
        chains=np_chains_list)
  np_chains_list = feature_processing.crop_chains(
      np_chains_list,
      msa_crop_size=feature_processing.MSA_CROP_SIZE,
      pair_msa_sequences=pair_msa_sequences,
      max_templates=feature_processing.MAX_TEMPLATES)
  np_example = msa_pairing.merge_chain_features(
      np_chains_list=np_chains_list, pair_msa_sequences=pair_msa_sequences,
      max_templates=feature_processing.MAX_TEMPLATES)
  np_example = feature_processing.process_final(np_example)
  return np_example

class DataPipeline:
  """Runs the alignment tools and assembles the input features."""

  def __init__(self):
      pass
  
  def process_monomer(self,
              input_fasta_path: str,
              msa_output_dir: str) -> pipeline.FeatureDict:
      with open(input_fasta_path) as f:
          input_fasta_str = f.read()
          input_seqs, input_descs = parsers.parse_fasta(input_fasta_str)
      a3m_file = os.path.join(msa_output_dir, 'mmseqs/uniref.a3m')
      print(a3m_file)
      #a3m_file = os.path.join(msa_output_dir,input_descs[0],'mmseqs/uniref.a3m')
      #a3m_file = '/home/thu/Downloads/alphafold/msa/' + input_descs[0] + '.a3m'
      #hhr_file = '/home/thu/Downloads/alphafold/msa/' + input_descs[0] + '.hhr'
      with open(a3m_file, "r") as fp:
          msa = parsers.parse_a3m(fp.read())
          data = {"msa": msa.sequences, "deletion_matrix": msa.deletion_matrix}
          sequence = msa.sequences[0]
          num_res = len(sequence)
          msas, deletion_matrices = zip(*[
          (data["msa"], data["deletion_matrix"])])
          chain_features = {
              **pipeline.make_sequence_features(sequence=sequence, description="none",
                                              num_res=num_res),
              **pipeline.make_msa_features((msa,))
          }
      #template_searcher = hhsearch.HHSearch(
      #  binary_path='/home/thu/miniconda3/bin/hhsearch',
      #  databases='/pool-data/data/thu/pdb70/pdb70')
      #template_featurizer = templates.HhsearchHitFeaturizer(
      #  mmcif_dir='/pool-data/data/thu/mmcif_files',
      #  max_template_date='2021-10-10',
      #  max_hits=4,
      #  kalign_binary_path=None,
      #  release_dates_path=None,
      #  obsolete_pdbs_path=None)
      #with open(hhr_file) as f:
      #    hhr = f.read()
      #pdb_template_hits = template_searcher.get_template_hits(
      #   output_string=hhr, input_sequence=sequence)
      #templates_result = template_featurizer.get_templates(
      #   query_sequence=sequence,
      #   hits=pdb_template_hits)
      #chain_features.update(templates_result.features)

      return chain_features


  def _process_single_chain(
      self,
      chain_id: str,
      sequence: str,
      description: str,
      msa_output_dir: str,
      use_mock_temp: bool,
      use_pdb_temp: bool,
      template_pdb_path: str,
      allow_msa_dup: bool) -> pipeline.FeatureDict:
      
      a3m_file = os.path.join(msa_output_dir,'mmseqs/uniref.a3m')
      #a3m_file = os.path.join(msa_output_dir,description,'mmseqs/uniref.a3m')
      # hhr_file = '/home/thu/Downloads/alphafold/msa/' + description + '.hhr'
      with open(a3m_file, "r") as fp:
          msa = parsers.parse_a3m(fp.read())
          data = {"msa": msa.sequences, "deletion_matrix": msa.deletion_matrix}
          sequence = msa.sequences[0]
          num_res = len(sequence)
          #print(sequence, num_res)
          msas, deletion_matrices = zip(*[
          (data["msa"], data["deletion_matrix"])])
          if allow_msa_dup:
              chain_features = {
                **pipeline.make_sequence_features(sequence=sequence, description="none",
                                              num_res=num_res),
                **make_msa_features_colab((msa,))
            }
          else:
              chain_features = {
                **pipeline.make_sequence_features(sequence=sequence, description="none",
                                              num_res=num_res),
                **pipeline.make_msa_features((msa,))
            }
      if use_mock_temp:
          chain_features.update(mk_mock_template(sequence, 4))
      if use_pdb_temp:
          ####################### template from pdb###############################
          template_model = ft.read_pdb_structure(template_pdb_path)
          chains = [i.get_id() for i in template_model.get_chains()]
          query_seq_len = len(sequence)
          output_seq = ['-'] * query_seq_len
          output_confidence_scores = np.full(query_seq_len, -1)
          output_positions = np.zeros((query_seq_len, residue_constants.atom_type_num, 3), dtype=np.float64)
          output_masks = np.zeros((query_seq_len, residue_constants.atom_type_num), dtype=np.int64)
          alaninize = False
          if(chain_id in chains):
              positions, masks, template_seq = ft.get_positions_masks_seq_for_chain(
                                                      template_model,
                                                      chain_id,
                                                      atom_types='all')
              mapping = ft.alignment_seq_mapping(sequence, template_seq)
              for q, t in mapping.items():
                  output_positions[q] = positions[t]
                  output_masks[q] = masks[t]
                  output_seq[q] = template_seq[t] if not alaninize else 'A'
          output_seq = ''.join(output_seq)
          templates_aatype = residue_constants.sequence_to_onehot(output_seq, residue_constants.HHBLITS_AA_TO_ID)
          template_features = {'template_all_atom_positions': output_positions,
                       'template_all_atom_masks': output_masks,
                       'template_sequence': output_seq.encode(),
                       'template_aatype': np.array(templates_aatype),
                       'template_domain_names': f'none'.encode(),
                       'template_release_date': f'none'.encode()}
          for feature in template_features:
              template_features[feature] = np.expand_dims(template_features[feature], axis=0)

          chain_features.update(template_features)

          ########################################################################
      # template_searcher = hhsearch.HHSearch(
      #   binary_path='/home/thu/miniconda3/bin/hhsearch',
      #   databases='/pool-data/data/thu/pdb70/pdb70')
      # template_featurizer = templates.HhsearchHitFeaturizer(
      #   mmcif_dir='/pool-data/data/thu/mmcif_files',
      #   max_template_date='2021-10-10',
      #   max_hits=4,
      #   kalign_binary_path=None,
      #   release_dates_path=None,
      #   obsolete_pdbs_path=None)
      # with open(hhr_file) as f:
      #     hhr = f.read()
      # pdb_template_hits = template_searcher.get_template_hits(
      #    output_string=hhr, input_sequence=sequence)
      # templates_result = template_featurizer.get_templates(
      #    query_sequence=sequence,
      #    hits=pdb_template_hits)
      # chain_features.update(templates_result.features)
      all_seq_features = pipeline.make_msa_features([msa])
      valid_feats = msa_pairing.MSA_FEATURES + (
        'msa_uniprot_accession_identifiers',
        'msa_species_identifiers',
      )
      feats = {f'{k}_all_seq': v for k, v in all_seq_features.items()
             if k in valid_feats}

      chain_features.update(feats)

      return chain_features
  
  def process(self,
              input_fasta_path: str,
              msa_output_dir: str,
              use_mock_temp: bool = False,
              use_pdb_temp: bool = False,
              template_pdb_path: str = None,
              allow_msa_dup: bool = False) -> pipeline.FeatureDict:
    """Runs alignment tools on the input sequences and creates features."""
    with open(input_fasta_path) as f:
      input_fasta_str = f.read()
    input_seqs, input_descs = parsers.parse_fasta(input_fasta_str)

    chain_id_map = _make_chain_id_map(sequences=input_seqs,
                                      descriptions=input_descs)
    chain_id_map_path = os.path.join(msa_output_dir, 'chain_id_map.json')
    with open(chain_id_map_path, 'w') as f:
      chain_id_map_dict = {chain_id: dataclasses.asdict(fasta_chain)
                           for chain_id, fasta_chain in chain_id_map.items()}
      json.dump(chain_id_map_dict, f, indent=4, sort_keys=True)

    all_chain_features = {}
    sequence_features = {}
    is_homomer_or_monomer = len(set(input_seqs)) == 1
    for chain_id, fasta_chain in chain_id_map.items():
      if fasta_chain.sequence in sequence_features:
        all_chain_features[chain_id] = copy.deepcopy(
            sequence_features[fasta_chain.sequence])
        continue
      chain_features = self._process_single_chain(
          chain_id=chain_id,
          sequence=fasta_chain.sequence,
          description=fasta_chain.description,
          msa_output_dir=msa_output_dir,
          use_mock_temp=use_mock_temp,
          use_pdb_temp=use_pdb_temp,
          template_pdb_path=template_pdb_path,
          allow_msa_dup=allow_msa_dup)

      chain_features = convert_monomer_features(chain_features,
                                                chain_id=chain_id)
      all_chain_features[chain_id] = chain_features
      sequence_features[fasta_chain.sequence] = chain_features

    all_chain_features = add_assembly_features(all_chain_features)
    
    if allow_msa_dup:
        np_example = pair_and_merge_colab(
            all_chain_features=all_chain_features,
        )
    else:
        np_example = feature_processing.pair_and_merge(
            all_chain_features=all_chain_features,
        )

    # Pad MSA to avoid zero-sized extra_msa.
    np_example = pad_msa(np_example, 512)
    return np_example

def predict_structure(
    output_dir_base,
    model_runners,
    random_seed,
    features,
    num_predictions,
    model_preset,
    enable_single_seed,
    single_seed
    ):
    timings = {}
    output_dir = os.path.join(output_dir_base, "")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    unrelaxed_pdbs = {}
    relaxed_pdbs = {}
    ranking_confidences = {}
    num_models = len(model_runners)
    for model_index, (model_name, model_runner) in enumerate(model_runners.items()):
        logging.info('Running model %s', model_name)
        t_0 = time.time()
        #model_random_seed = model_index + random_seed * num_models
        if enable_single_seed:
            model_random_seed = single_seed
        else:
            model_random_seed = int(random_seed[model_index % num_predictions])
        timings[f'process_features_{model_name}'] = time.time() - t_0

        t_0 = time.time()
        processed_feature_dict = model_runner.process_features(
         features, random_seed=model_random_seed)
        prediction_result = model_runner.predict(processed_feature_dict,
                                             random_seed=model_random_seed)
        t_diff = time.time() - t_0
        timings[f'predict_and_compile_{model_name}'] = t_diff

        plddt = prediction_result['plddt']
        ranking_confidences[model_name] = prediction_result['ranking_confidence']

        cf = {}
        confidence_path = os.path.join(output_dir, f'confidence_{model_name}.txt')
        if 'multimer' in model_preset:
            score_list = ['max_predicted_aligned_error', 'ptm', 'iptm', 'ranking_confidence']
            score_benchmark = prediction_result['ranking_confidence'] * 100
        elif model_preset == 'monomer_ptm':
            score_list = ['max_predicted_aligned_error', 'ptm', 'ranking_confidence']
            score_benchmark = prediction_result['ranking_confidence']
        else:
            score_list = ['ranking_confidence']
            score_benchmark = prediction_result['ranking_confidence']
        for k in score_list:
            cf[k] = prediction_result[k].tolist()
        json.dump(cf, codecs.open(confidence_path, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)

        plddt_b_factors = np.repeat(
            plddt[:, None], residue_constants.atom_type_num, axis=-1)
        unrelaxed_protein = protein.from_prediction(
            features=processed_feature_dict,
            result=prediction_result,
            b_factors=plddt_b_factors,
            remove_leading_feature_dimension=not model_runner.multimer_mode)

        unrelaxed_pdbs[model_name] = protein.to_pdb(unrelaxed_protein)
        unrelaxed_pdb_path = os.path.join(output_dir, f'unrelaxed_{model_name}.pdb')
        with open(unrelaxed_pdb_path, 'w') as f:
            f.write(unrelaxed_pdbs[model_name])
        if FLAGS.score_threshold != None and score_benchmark >= FLAGS.score_threshold:
            quit()

def main(argv):
    dp = DataPipeline()
    if FLAGS.use_pdb_temp:
        assert FLAGS.use_mock_temp==False
        assert FLAGS.pdb_template_path!=None
    if FLAGS.use_mock_temp:
        assert FLAGS.use_pdb_temp==False
        assert FLAGS.pdb_template_path==None
    if FLAGS.model_preset == 'monomer' or FLAGS.model_preset == 'monomer_ptm':
        multimer_feature = dp.process_monomer(FLAGS.fasta_dir,FLAGS.pre_align_dir)
    else:
        multimer_feature = dp.process(FLAGS.fasta_dir,FLAGS.pre_align_dir, use_mock_temp=FLAGS.use_mock_temp, use_pdb_temp=FLAGS.use_pdb_temp, template_pdb_path=FLAGS.pdb_template_path, allow_msa_dup=FLAGS.allow_msa_duplicate)
    for k,v in multimer_feature.items():
        print(k, v.shape)
    #np.savez('feat.npz',aatype=multimer_feature['aatype'], 
    #        residue_index=multimer_feature['residue_index'],
    #        seq_length=multimer_feature['seq_length'],
    #        msa=multimer_feature['msa'],
    #        num_alignments=multimer_feature['num_alignments'],
    #        asym_id=multimer_feature['asym_id'],
    #        sym_id=multimer_feature['sym_id'],
    #        entity_id=multimer_feature['entity_id'],
    #        deletion_matrix=multimer_feature['deletion_matrix'],
    #        deletion_mean=multimer_feature['deletion_mean'],
    #        all_atom_mask=multimer_feature['all_atom_mask'],
    #        all_atom_positions=multimer_feature['all_atom_positions'],
    #        assembly_num_chains=multimer_feature['assembly_num_chains'],
    #        entity_mask=multimer_feature['entity_mask'],
    #        cluster_bias_mask=multimer_feature['cluster_bias_mask'],
    #        bert_mask=multimer_feature['bert_mask'],
    #        seq_mask=multimer_feature['seq_mask'],
    #        msa_mask=multimer_feature['msa_mask'],
    #        template_all_atom_positions=multimer_feature['template_all_atom_positions'],
    #        template_all_atom_mask=multimer_feature['template_all_atom_mask'],
    #        template_aatype=multimer_feature['template_aatype'],
    #        num_templates=multimer_feature['num_templates'])
    #quit()
    if len(argv) > 1:
        raise app.UsageError('Too many command-line arguments.')
    run_multimer_system = 'multimer' in FLAGS.model_preset
    if FLAGS.model_preset == 'monomer_casp14':
        num_ensemble = 8
    else:
        num_ensemble = 1

    model_runners = {}
    model_names = config.MODEL_PRESETS[FLAGS.model_preset]
    for model_name in model_names:
        model_config = config.model_config(model_name)
        if run_multimer_system:
            if FLAGS.use_mock_temp or FLAGS.use_pdb_temp:
                model_config.model.embeddings_and_evoformer.template.enabled = True
            model_config.model.global_config.bfloat16 = FLAGS.use_float16
            model_config.model.num_ensemble_eval = num_ensemble
            model_config.model.num_recycle = FLAGS.num_recycle
        else:
            model_config.data.eval.num_ensemble = num_ensemble
            model_config.data.common.num_recycle = FLAGS.num_recycle
        #model_config.model.embeddings_and_evoformer.masked_msa.replace_fraction = 0.0

        model_params = data.get_model_haiku_params(model_name=model_name, data_dir=FLAGS.data_dir)
        #print(model_params['alphafold/alphafold_iteration/evoformer/preprocess_1d'])
        model_runner = model.RunModel(model_config, model_params)
        for i in range(FLAGS.num_prediction):
           model_runners[f'{model_name}_pred_{i}'] = model_runner
    #random_seed = FLAGS.random_seed
    #if random_seed is None:
    #    random_seed = random.randrange(sys.maxsize // len(model_names))
    random_seed = np.linspace(0, sys.maxsize-512, num=FLAGS.num_prediction)
    predict_structure(FLAGS.output_dir, model_runners, random_seed, multimer_feature, FLAGS.num_prediction, FLAGS.model_preset, FLAGS.enable_single_seed, FLAGS.single_seed)

if __name__ == '__main__':
    flags.mark_flags_as_required([
      'output_dir',
      'data_dir',
      'fasta_dir',
      'pre_align_dir'
    ])

    app.run(main)



