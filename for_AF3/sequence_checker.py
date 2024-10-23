####
# prepare input sequence to test AlphaFold
# For each data point from RCSB, the sequence in fasta file is usually
# different than the sequence in protein structure
# This script aims for making the fasta file consistent with sequence from structure
# to make comparison between the structure predicted by AlphaFold and native crystal structure easier
# The basic idea is
#  1. get pdb_sequence for each chain from .pdb use pymol
#  2. get fasta_sequence for each chain from .fasta (from RCSB)
#  3. for each corresponding chain,
#       a. If there is additional subsequence in fasta_sequence compare to pdb_sequence at the beginning,
#           trim the beginning extra subsequence in fasta_sequence
#       b. If there is additional subsequence in fasta_sequence compare to pdb_sequence at the end,
#           trim the end extra subsequence in fasta_sequence
# 4. save the output chain, sequence in fasta format for fasta_sequence
# the reason that we modified the fasta_sequence is it rarely contains gap, while gap(s) could exist in
# pdb_sequence
####

import os
import pymol
import pandas as pd
from pymol import cmd

from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def get_pdb_sequence(pdb_file_path):
    """

    :param pdb_file_path:
    :return: dict of (chain, sequence) pair
    """
    pdb_sequence = {}
    cmd.reinitialize()
    cmd.load(pdb_file_path)
    chains = cmd.get_chains()

    for chain in chains:
        sequence = cmd.get_fastastr(f"chain {chain}")
        lines = sequence.splitlines()
        seq = ""
        for line in lines:
            if line.startswith(">"):
                chain_name = line[-1]
            else:
                seq += line.strip("\n")
        pdb_sequence[chain_name] = seq
    return pdb_sequence


def get_fasta_from_df(df, pdb_id):
    select_sequence = df[df['PDB_id'] == pdb_id]
    fasta_sequence = {}
    for row in select_sequence.iterrows():
        fasta_sequence[row[1]['native']] = row[1]['sequence']
    return fasta_sequence


def trim_fasta_sequence(fasta_sequence, pdb_sequence):
    for chain, seq in fasta_sequence.items():
        reference_seq=pdb_sequence[chain]
        start,end = trim_function(reference_seq,seq)
        fasta_sequence[chain]=seq[start:end]
    return fasta_sequence


def get_fasta_sequence():
    raise NotImplemented


def trim_function(seq1, seq2):
    alignments = pairwise2.align.globalms(seq1, seq2,
                                          match=2, mismatch=-1, open=-10, extend=-0.5)
    aligned_seq1, aligned_seq2, score, start, end = alignments[0]
    first_alignment_pos = next(i for i, (a, b) in enumerate(zip(aligned_seq1, aligned_seq2)) if a != '-' and b != '-')
    last_align_pos = max(i for i,(a,b) in enumerate(zip(aligned_seq1, aligned_seq2)) if a!='-' and b!='-')
    return first_alignment_pos,last_align_pos


#seq1 = "ENLWVTVYYGVPVWKDAETTLFCASDAEKHNVWATHACVPTDPNPQEIHLENVTEEFNMWKNNMVEQMHEDIISLWDQSLKPCVKLTPLCVTLNCTNVTNNITDDMRGELKNCSFNATTELRNKRVKRYSLFYRLDIVQIDSNRTKSHYRLINCNTSAITQACPKVSFEPIPIHYCAPAGFAILKCKDKKFNGTGPCPSVSTVQCTHGIKPVVSTQLLLNGSLAEEEVIIRSENITNNAKNILVQLNTPVQINCTRPNNNTVKSIRIGPGQAFYYTGDIIGDIRQAHCNVSKATWNETLGKVVKQLRKHFGNNTIIRFAQSSGGDLEVTTHSFNCGGEFFYCNTSGLFNSTWISNNDSITLPCRIKQIINMWQRIGQAMYAPPIQGVIRCVSNITGLILTRDGGSTNSTTETFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTRCK"
#seq2 = "MGILPSPGMPALLSLVSLLSVLLMGCVAETGAENLWVTVYYGVPVWKDAETTLFCASDAKAYETEKHNVWATHACVPTDPNPQEIHLENVTEEFNMWKNNMVEQMHEDIISLWDQSLKPCVKLTPLCVTLNCTNVTNNITDDMRGELKNCSFNATTELRNKRVKRYSLFYRLDIVQIDSNRTKSHYRLINCNTSAITQACPKVSFEPIPIHYCAPAGFAILKCKDKKFNGTGPCPSVSTVQCTHGIKPVVSTQLLLNGSLAEEEVIIRSENITNNAKNILVQLNTPVQINCTRPNNNTVKSIRIGPGQAFYYTGDIIGDIRQAHCNVSKATWNETLGKVVKQLRKHFGNNTIIRFAQSSGGDLEVTTHSFNCGGEFFYCNTSGLFNSTWISNTSVQGSNSTGSNDSITLPCRIKQIINMWQRIGQAMYAPPIQGVIRCVSNITGLILTRDGGSTNSTTETFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTRCKRRVVGRRRRRR"


pdb_sequence=get_pdb_sequence("/projectnb2/docking/imhaoyu/24_epitope_mapping/database/ABAG/ABAG-complex-structure-benchmark-dataset/ABAG-Docking_benchmark_dataset/87cases/7SJO/7SJO_l_u.pdb")
print(pdb_sequence)
#print(trim_function(seq2, seq1))
#start,end=trim_function(seq2, seq1)
#print(seq2[start:end])
# print(get_pdb_sequence("/projectnb2/docking/imhaoyu/24_epitope_mapping/database/ABAG/complexes_labeled/bound/7T77.pdb"))
#df = pd.read_csv(
#    "/projectnb2/docking/imhaoyu/24_epitope_mapping/database/ABAG/AlphaFold3/benchmark/AfterSept302021.csv", sep="\t")

#df = pd.read_csv('/projectnb2/docking/imhaoyu/24_epitope_mapping/database/ABAG/AlphaFold3/benchmark/AfterSept302021_unbound.csv')
fasta={"H":"EVQLVQSGAEVKKPGASVKVSCKASGYKFTDSEMHWVRQAPGQGLEWIGGVDPETEGAAYNQKFKGRATITRDTSTSTAYLELSSLRSEDTAVYYCTRGYDYDYALDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHT",
       "L":"DIQMTQSPSSLSASVGDRVTITCRASSSVEFIHWYQQKPGKAPKPLISATSNLASGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQWSSAPWTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"} #get_fasta_from_df(df, "7T77")
fasta={"A":"MGQEDPNSLRHKYNFIADVVEKIAPAVVHIELFRKLPFSKREVPVASGSGFIVSEDGLIVTNAHVVTNKHRVKVELKNGATYEAKIKDVDEKADIALIKIDHQGKLPVLLLGRSSELRPGEFVVAIGSPFSLQNTVTTGIVSTTQRGGKELGLRNSDMDYIQTDAIINYGNAGGPLVNLDGEVIGINTLKVTAGISFAIPSDKIKKFLTESHDRQAKGKAITKKKYIGIRMMSLTSSKAKELKDRHRDFPDVISGAYIIEVIPDTPAEAGGLKENDVIISINGQSVVSANDVSDVIKRESTLNMVVRRGNEDIMITVIPEEIDPLEHHHHHH"}
print(trim_fasta_sequence(fasta,pdb_sequence))