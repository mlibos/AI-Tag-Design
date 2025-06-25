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
from pprint import pprint
from sympy.physics.units import percent
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import colors
from matplotlib.ticker import PercentFormatter


def build_record(fasta_file):
    """
    This function builds a record from a fasta file of sequences. The
    record is a dictionary object where keys are sequence IDs and values
    are sequences
    """
    record = {}
    for index, seq_record in enumerate(SeqIO.parse(fasta_file, 'fasta')):
        record[f'{seq_record.id[0:-2]}{index}'] = seq_record.seq
    return record

def extract_values(fasta_file):
    """
    This function extracts the ligand_confidence score and the sequence recovery score
    and places them into their own distribution sets.
    """
    confidence_scores = []
    recovery_scores = []
    for index, seq_record in enumerate(SeqIO.parse(fasta_file, 'fasta')):
        if index == 0:
            continue
        description = seq_record.description.split(', ')
        confidence_scores.append(float(description[-2][18:]))
        recovery_scores.append(float(description[-1][8:]))
    return confidence_scores, recovery_scores

def plot_scores(confidence_scores, recovery_scores, tag_name):
    """
    This function plots the confidence scores and the recovery scores
    in a histogram.
    """
    n_bins = 10
    fig, axs = plt.subplots(1, 2)
    axs[0].hist(confidence_scores, bins=n_bins, color='skyblue')
    axs[0].set_xlabel('Confidence Scores')
    axs[0].set_ylabel('Count')
    axs[1].hist(recovery_scores, bins=n_bins, color='skyblue')
    axs[1].set_xlabel('Recovery Scores')
    plt.title(f'Tag Confidence and Recovery Scores: Tag {tag_name}', loc='right')
    plt.savefig(f'figures/{tag_name}_scores.png', dpi=300)




def extract_original(fasta_file='6i2g.fasta'):
    original_record = {}
    chains = 'AB'
    for index, chain in enumerate(SeqIO.parse(fasta_file, 'fasta')):
        original_record[f'6i2g_{chains[index]}'] = chain.seq
    return original_record

def nanobody_statistics(sequences):
    """
    Generate statistics on the nanobodies generated for each unique tag.
    """
    # get avg hamming distance
    aligner = Align.PairwiseAligner(mismatch_score = 0, match_score = 1)
    aligner.open_gap_score = -1000000000000
    ids = list(sequences.keys())
    total_hamming_distance = 0
    for index, nb in enumerate(ids):
        if index + 1 < len(ids):
            for nb2 in ids[index+1:]:
                hd = len(sequences[nb]) - aligner.score(sequences[nb], sequences[nb2])
                total_hamming_distance += hd
    total_pairs = round(len(ids)*(len(ids)-1)/2, 2)
    avg_hamming_distance_self = round(total_hamming_distance / total_pairs, 2)

    # get avg hamming distance from original sequence
    total_hamming_distance = 0
    original_seq = extract_original()
    for index, nb in enumerate(ids):
        hd = len(original_seq['6i2g_A']) - aligner.score(original_seq['6i2g_A'], sequences[nb][:124])
        total_hamming_distance += hd
    avg_hamming_distance_original = round(total_hamming_distance / len(ids), 2)

    # get average residue per pos
    residue_probabilities = [(i, {}) for i in range(124)]
    for index, nb in enumerate(ids):
        for pos, residue in enumerate(sequences[nb][:124]):
            if residue in residue_probabilities[pos][1]:
                residue_probabilities[pos][1][residue] += 1
            else:
                residue_probabilities[pos][1][residue] = 1
    deletions = []
    for pos in range(len(residue_probabilities)):
        if len(residue_probabilities[pos][1]) <= 1:
            deletions.append(pos)
        for residue in residue_probabilities[pos][1]:
            residue_probabilities[pos][1][residue] = round(residue_probabilities[pos][1][residue] / len(ids), 4)
    # filter out unchanged residues
    residue_probabilities = [residue_probabilities[pos] for pos in range(len(residue_probabilities)) if pos not in deletions]
    return residue_probabilities, avg_hamming_distance_self, avg_hamming_distance_original



def main(files):
    """
    Main run file to produce statistics for each fasta file input.
    """
    for fasta_file in files:
        fn = fasta_file[:-3]
        seqs = build_record(fasta_file)
        residues, self_ham, original_ham = nanobody_statistics(seqs)
        print(f'{fn} residue probabilities: {residues}')
        print(f'{fn} average hamming score amongst self: {self_ham}')
        print(f'{fn} average hamming score vs original: {original_ham}')

def scores_main(files):
    """
    Creates plots for each file given.
    """
    for fasta_file in files:
        tag_name = fasta_file[:-3]
        confidence_scores, recovery_scores =extract_values(fasta_file)
        print(f'{tag_name} average confidence score: {round(sum(confidence_scores)/len(confidence_scores), 4)}')
        print(f'{tag_name} average recovery score: {round(sum(recovery_scores)/len(recovery_scores), 4)}')
        plot_scores(confidence_scores, recovery_scores, tag_name)



if __name__ == '__main__':
    start = time.time()
    files = ['nanobodies_1181.fa', 'nanobodies_3352.fa', 'nanobodies_6813.fa', 'nanobodies_8829.fa', 'nanobodies_9901.fa']
    main(files)
    scores_main(files)
    end = time.time()
    print(end - start)