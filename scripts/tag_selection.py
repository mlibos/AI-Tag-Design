from Bio import *
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio import SeqIO
from Bio import Align
from Bio import Blast
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA
import time
from pathlib import Path
import csv
start = time.time()


def load_tags(files):
    """
    This function takes in a list of fasta files with viable filtered tag
    designs and loads them into a dictionary of labeled sequences.
    """
    sequences = {}
    for file in files:
        file_name = file[23:-14]
        for seq_record in SeqIO.parse(file, 'fasta'):
            sequences[f'{file_name}_{seq_record.id[13:]}'] = seq_record.seq
    return sequences

def diversity_selection(sequences, candidates=100):
    """
    This function takes in a dictionary of labeled sequences and filters
    down to the most diverse candidate sequences.
    """
    aligner = Align.PairwiseAligner(mismatch_score = 0, match_score = 1)
    aligner.open_gap_score = -1000000000000

    diversity_candidates = [seq_id for seq_id in sequences.keys()]
    diversity_scores = {}
    for candidate in diversity_candidates:
        diversity_scores[candidate] = 0
    for i1, c1 in enumerate(diversity_candidates[:-1]):
        for i2, c2 in enumerate(diversity_candidates[i1+1:]):
            hd = len(sequences[c1]) - aligner.score(sequences[c2], sequences[c1])
            diversity_scores[c1] += hd
            diversity_scores[c2] += hd
    for candidate in diversity_scores:
        diversity_scores[candidate] = round(diversity_scores[candidate] / (len(sequences)-1), 2)
    top_candidates = sorted(diversity_scores, key=lambda score: diversity_scores[score], reverse=True)
    top_seqs = []
    esm_candidates = 0
    mpnn_candidates = 0
    for candidate in top_candidates:
        if mpnn_candidates + esm_candidates >= candidates:
            break
        elif candidate[0] == 'h' and esm_candidates < 50:
            top_seqs.append([candidate, sequences[candidate], diversity_scores[candidate]])
            esm_candidates += 1
        elif candidate[0] == 'h' and esm_candidates >= 50:
            pass
        else:
            top_seqs.append([candidate, sequences[candidate], diversity_scores[candidate]])
            mpnn_candidates += 1
    return top_seqs

def convert_to_fasta(sequences):
    """
    Given a dictionary of sequences this function creates a fasta file
    of each sequence id:sequence.
    """
    fasta_file_name = f'top_tag_candidates.fasta'
    new_fasta = []
    for sequence in sequences:
        new_fasta.append('>%s\n%s' % (sequence, sequences[sequence]))
    with open(fasta_file_name, 'w') as f:
        f.write('\n'.join(new_fasta))

def convert_to_csv(seqs_list):
    """
    Given a list of sequences in this format:
    [ [id1, seq1, diversity_score1], [id2, seq2, diversity_score2], ... ]
    converts this list into a csv file in the same format.
    Headers = seq_id, seq, diversity_score
    """
    with open('top_tag_candidates.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['seq_id', 'seq', 'diversity_score'])
        writer.writerows(seqs_list)

def main(files):
    """
    This function takes in a list of fasta files with viable filtered tags
    designs and loads them into a dictionary of labeled sequences. It then
    selects top sequences based on diversity scores and ensures it selects the top
    50 candidates from each sequence bank. It then converts this data to a csv file.
    """
    for index in range(len(files)):
        files[index] = 'filtered_designs/' + files[index]
    sequences = load_tags(files)
    top_seqs = diversity_selection(sequences)
    convert_to_csv(top_seqs)
    return top_seqs

if __name__ == '__main__':
    tag_files = ['final_high_temp_only_tag_designs.fasta', 'final_only_tag_designs.fasta']
    main(tag_files)
    end = time.time()
    print(end - start)