# This python file is intended to generate statistical data on the
# different outputs for generating tag design. The three given input
# models currently are: generating tags with fixed residues in the A chain,
# generating tags with an entirely fixed A chain, and generating tags
# with only the original Alfa Tag as given input (no A chain).

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
import pprint
from sympy.physics.units import percent

def build_record(fasta_file, tag_name, tag_length=13):
    """
    # build our records dictionary with the original 6i2g sequence
    # and each of our initially generated new alfa tags
    """
    records = {}
    index = 0
    for seq_record in SeqIO.parse(fasta_file,'fasta'):
        if index == 0:
            records[seq_record.id] = [seq_record.seq[-tag_length:]]
        else:
            seq_record.id = f"{tag_name}_{index}"
            records[seq_record.id] = [seq_record.seq[-tag_length:]]
        index += 1
    return records

def positional_statistics(sequences, sequence_length=13):
    """
    Given sequences of residues, this function generates a list of dictionaries
    where each index's dictionary contains the residue statistics for that index
    in the sequences.
    """
    positions = {}
    for index in range(sequence_length):
        index_dict = {}
        for seq_record in sequences:
            value = sequences[seq_record][0][index]
            if value in index_dict:
                index_dict[value] += 1
            else:
                index_dict[value] = 1
        for residue in index_dict:
            index_dict[residue] = str(round(100*index_dict[residue] / len(sequences), 3)) + '%'
        positions[index] = index_dict
    return positions


def structural_data(sequences, sequence_length=13):
    """
    This function calculates the structural data of each sequence.
    Currently, this includes:
    - GRAVY Score
    - Charge at Physiological pH (7.4)
    - Presence of Unwanted Residues (K, W, M, C)
    """
    for seq_record in sequences:
        sequence = sequences[seq_record][0]
        if 'C' in sequence or 'M' in sequence or 'W' in sequence or 'K' in sequence:
            unwanted_residues = {'Unwanted Residues:': True}
        else:
            unwanted_residues = {'Unwanted Residues:': False}
        protein = PA(sequence)
        charge_at_ph = {'Charge at pH 7.4:': round(protein.charge_at_pH(7.4), 3)}
        gravy_score = {'GRAVY Score:': round(protein.gravy(), 3)}
        sequences[seq_record].append(unwanted_residues)
        sequences[seq_record].append(charge_at_ph)
        sequences[seq_record].append(gravy_score)
    return sequences

def collect_structural_data(sequences, data_structure, sequence_length=13):
    """
    Here we collect the average GRAVY score, charge at physiological pH, and
    frequency of unwanted residues and store the information into our data structure.
    """
    avg_gravy_score = 0
    avg_charge_at_pH = 0
    percent_unwanted_residues = 0
    for seq_record in sequences:
        avg_gravy_score += sequences[seq_record][3]['GRAVY Score:']
        avg_charge_at_pH += sequences[seq_record][2]['Charge at pH 7.4:']
        if sequences[seq_record][1]['Unwanted Residues:']:
            percent_unwanted_residues += 1
    stats =    {'Average GRAVY Score:': round(avg_gravy_score/len(sequences), 3),
                'Average Charge at PH 7.4:': round(avg_charge_at_pH/len(sequences), 3),
                'Frequency of unwanted residues': round(percent_unwanted_residues/len(sequences), 3),}
    data_structure.insert(0, stats)
    return data_structure


if __name__ == '__main__':
    data_streams = ['fixed_chain_designs.fa', 'fixed_residues_designs.fa', 'only_tag_designs.fa', 'high_temp_only_tag_designs.fa']
    fixed_chain_sequences = structural_data(build_record(f'tag_design_outputs/{data_streams[0]}', 'fixed_chain_tags'))
    fixed_chain_data = collect_structural_data(fixed_chain_sequences, [positional_statistics(fixed_chain_sequences)])
    print('Fixed Chain Design Data:', fixed_chain_data)
    fixed_residues_sequences = structural_data(build_record(f'tag_design_outputs/{data_streams[1]}', 'fixed_residues_tags'))
    fixed_residues_data = collect_structural_data(fixed_residues_sequences, [positional_statistics(fixed_residues_sequences)])
    print('Fixed Residues Design Data:',fixed_residues_data)
    only_tag_sequences = structural_data(build_record(f'tag_design_outputs/{data_streams[2]}', 'only_tag_tags'))
    only_tag_data = collect_structural_data(only_tag_sequences, [positional_statistics(only_tag_sequences)])
    print('Only Tag Design Data:',only_tag_data)
    high_temp_only_tag_sequences = structural_data(build_record(f'tag_design_outputs/{data_streams[3]}', 'high_temp_only_tag_tags'))
    high_temp_only_tag_data = collect_structural_data(high_temp_only_tag_sequences, [positional_statistics(high_temp_only_tag_sequences)])
    print('High Temp Only Tag Design Data:', high_temp_only_tag_data)