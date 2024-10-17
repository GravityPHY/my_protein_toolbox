from Bio import pairwise2
from Bio.pairwise2 import format_alignment


query_sequence = "HHHHHMYLSEQHH"
sequence_dict = {"A":"MYLSEQTKG",
                 "B":"MYLSEQ",
                 "C":"AMYLSEQ"}


def find_most_similar_sequence(query_seq, sequences):
    """

    :param query_seq:
    :param sequences:
    :return:
    """

    best_match=None
    best_score=float("-inf")
    for chain_id, seq in sequences.items():
        alignments = pairwise2.align.globalxx(query_seq, seq)
        score = alignments[0][2] if alignments else 0
        if score>best_score:
            best_score=score
            best_match=(chain_id, seq,score)
    return best_match
print(find_most_similar_sequence(query_sequence, sequence_dict))
