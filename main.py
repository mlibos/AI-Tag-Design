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
start = time.time()

def build_record(fasta_file, tag_length=13):
    """
    # build our records dictionary with the original 6i2g sequence
    # and each of our initially generated new alfa tags
    """
    records = {}
    index = 0
    for seq_record in SeqIO.parse(fasta_file, 'fasta'):
        if index == 0:
            records[seq_record.id] = seq_record.seq[-tag_length:]
        else:
            seq_record.id = f"New_Alfa_tag_{index}"
            records[seq_record.id] = seq_record.seq[-tag_length:]
        index += 1
    return records

def filter_similar_sequences(sequences, initial_protein, passing_score=4, protein_length=13):
    """
    perform alignment comparison between each sequence and each subsequent
    sequence to filter out similar sequences. We want a minimum of 5 residue
    difference between individual sequences. This assumes proteins are all of4
    length 13 unless otherwise specified.
    """
    # simple aligner used since all proteins are same length
    aligner = Align.PairwiseAligner(mismatch_score = 0, match_score = 1)
    aligner.open_gap_score = -1000000000000
    # protein must have >= 5 residue changes from all other considered proteins
    proteins = sequences.keys()
    similarity_filtered_proteins = [initial_protein]
    for i1, p1 in enumerate(proteins):
        # check each protein in our list to any protein added to our
        # filtered list
        add_protein = True
        for i2, p2 in enumerate(similarity_filtered_proteins):
            # check score between p1, p2, convert to a residue change score
            score = protein_length - aligner.score(sequences[p1], sequences[p2])
            # print(similarity_filtered_proteins)
            # print(f"{p1}\t{p2}\t{score}")
            if score < passing_score:
                add_protein = False
                break
        if add_protein:
            similarity_filtered_proteins.append(p1)
    return similarity_filtered_proteins

def filter_protein_features(sequences):
    """
    This function takes an input list of sequences and filters out proteins based
    on the following criteria:
    -Net charge should be between +1 and -1 inclusive.
    -No cysteines, methionines, lysines or tryptophans
    -Should have a negative GRAVY score (hydrophilic)
    """
    # initially filter out the sequences with lysines (K), cysteines (C), methionines (M), Tryptophans (W)
    # initially all proteins are valid
    valid_proteins = list(sequences.keys())
    for sequence in sequences:
        if 'C' in sequences[sequence] or 'M' in sequences[sequence] or 'W' in sequences[sequence] or 'K' in sequences[sequence]:
            valid_proteins.remove(sequence)
    # now filter out anything containing a net charge below -1.0 or above 1.0 at neutral pH
    electrically_charged_proteins = []
    for valid_protein in valid_proteins:
        seq = sequences[valid_protein]
        protein = PA(seq)
        charge_at_ph_7 = protein.charge_at_pH(7.4)
        if charge_at_ph_7 < -1.0 or charge_at_ph_7 > 1.0:
            electrically_charged_proteins.append(valid_protein)
        # print(f"{valid_protein}: {charge_at_ph_7}")
        # print(f'{valid_protein}: {protein.isoelectric_point()}')
    valid_proteins = [protein for protein in valid_proteins if protein not in electrically_charged_proteins]

    hydrophobic_proteins = []
    # now we calculate GRAVY scores and filter out any protein with a hydrophobic score
    for valid_protein in valid_proteins:
        seq = sequences[valid_protein]
        protein = PA(seq)
        gravy_score = protein.gravy()
        if gravy_score >= 0.0:
            hydrophobic_proteins.append(valid_protein)
        # print(f"{valid_protein} gravy score: {gravy_score}")
    valid_proteins = [protein for protein in valid_proteins if protein not in hydrophobic_proteins]
    return valid_proteins


def blast_filter_proteins(sequences, email='munir.libos@proteininnovation.org', blast_output_file='blast_results.xml'):
    """
    This function takes an input dict of sequences and filters out proteins based
    on if they return a blast HIT or not.
    """
    def convert_sequences_to_fasta():
        """
        Given an input dictionary of sequences (in the format of protein name: sequence)
        this function writes these sequences to a fasta file format (for blast implementation)
        """
        fasta_file_name = 'filtered_sequences.fasta'
        new_fasta = []
        for sequence in sequences:
            if sequence != '6i2g,':
                new_fasta.append('>%s\n%s' % (sequence, sequences[sequence]))
        with open(fasta_file_name, 'w') as f:
            f.write('\n'.join(new_fasta))
        return fasta_file_name

    # set blast email for blasting over internet
    Blast.email = email
    found_blast_proteins = []

    blast_path = Path(blast_output_file)
    if blast_path.is_file():
        # output file exists no need to blast
        print('Blast File exists! Using local Blast results')
    else:
        print('Blast File does not exist! Using internet Blast search results....')
        fasta_file = convert_sequences_to_fasta()
        fasta_string = open(fasta_file).read()
        result_handle = Blast.qblast(
            "blastp",
            "pdb",
            fasta_string,
        )
        with open('blast_results.xml', 'wb') as out_stream:
            out_stream.write(result_handle.read())
        result_handle.close()

    # open blast file and read results
    result_stream = open(blast_output_file, 'rb')
    blast_records = Blast.parse(result_stream)
    for blast_record in blast_records:
        if len(blast_record) > 0:
            found_blast_proteins.append(blast_record.id)
    return found_blast_proteins


def generate_full_complexes(sequences, original_structure, tag_name, tag_length=13, sort=False, sort_number=0):
    """
    This function takes an input dictionary of sequences containing the 13 residue
    length tags and creates a fasta file containing each of the new tags imposed onto
    the new nbAlfa complex.

    Optionally, this function can take a sort_number and filter down the sequences
    to obtain the nth most diverse sequences
    """
    # save all tag structures in a different fasta file
    fasta_file_name = f'final_{tag_name}.fasta'
    new_fasta = []
    new_proteins = {}
    for sequence in sequences:
        if sequence != '6i2g,':
            # manually adding back in proline ends
            tag = "P" + sequences[sequence] + "P"
            new_proteins[sequence] = tag
            new_fasta.append('>%s\n%s' % (sequence, tag))
    with open(fasta_file_name, 'w') as f:
        f.write('\n'.join(new_fasta))


    # trim away 6i2g sequence
    sequences = {key: sequences[key] for key in sequences if key != '6i2g,'}
    proteins = list(sequences.keys())
    avg_hamming_distances = {key: 0 for key in sequences}
    # sort down sequences if sort
    if sort:
        # simple aligner used since all proteins are same length
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.open_gap_score = -1000000000000
        for i1, p1 in enumerate(proteins):
            for i2, p2 in enumerate(proteins[i1:]):
                if p1 == p2:
                    continue
                hamming_distance = tag_length - aligner.score(sequences[p1], sequences[p2])
                avg_hamming_distances[p1] += hamming_distance
                avg_hamming_distances[p2] += hamming_distance
        for p in avg_hamming_distances:
            avg_hamming_distances[p] /= (len(proteins)-1)
            avg_hamming_distances[p] = round(avg_hamming_distances[p], 2)
        most_diverse_proteins = proteins[:]
        most_diverse_proteins.sort(key=lambda x: avg_hamming_distances[x], reverse=True)
        most_diverse_proteins = most_diverse_proteins[:sort_number]
        sequences = {key: sequences[key] for key in sequences if key in most_diverse_proteins}

    # generate a new fasta file containing our new nbAlfa complexes
    fasta_file_name = f'final_{tag_name}_complexes.fasta'
    new_fasta = []
    new_proteins = {}
    for sequence in sequences:
        if original_structure == '6I2G.fa':
            original_sequence_without_tag = next(SeqIO.parse(original_structure, 'fasta')).seq
            # manually adding back in proline ends
            complex = original_sequence_without_tag + "/" + "P" + sequences[sequence] + "P"
            new_proteins[sequence] = complex
            new_fasta.append('>%s\n%s' % (sequence, complex))
        else:
            full_sequences = {}
            index = 0
            for seq_record in SeqIO.parse(f'tag_design_outputs/{original_structure}', 'fasta'):
                if index != 0:
                    seq_record.id = f"New_Alfa_tag_{index}"
                    full_sequences[seq_record.id] = seq_record.seq[:-tag_length]
                index += 1
            complex = full_sequences[sequence] + "P" + sequences[sequence] + "P"
            new_proteins[sequence] = complex
            new_fasta.append('>%s\n%s' % (sequence, complex))
    with open(fasta_file_name, 'w') as f:
        f.write('\n'.join(new_fasta))
    return fasta_file_name, new_proteins

def main(input_file, tag_name, original_structure, tag_length=13):
    """
    Main pipeline function to process fasta files of de novo
    tag designs. This function should take a fasta file name
    as input, process it into a record of sequences, filter the
    sequences based on similarity, structure, BLAST results, etc
    and then add prolines to tag and write the final results into
    two fasta files: one for all tags found, and one for the diversity
    filtered final six tags.
    """
    sequences = build_record(f'tag_design_outputs/{input_file}')
    initial_proteins = len(sequences)-1
    print(f'Initial Amount of {tag_name} Tags: {initial_proteins}')
    # filter by rules
    structurally_filtered_sequences = {key: sequences[key] for key in filter_protein_features(sequences)}
    # filter by similarity
    similarity_filtered_sequences = {key: structurally_filtered_sequences[key] for key in filter_similar_sequences(structurally_filtered_sequences, '6i2g,')}

    # filter by BLAST results (CURRENTLY DEPRECATED)
    # blast_filtered_sequences = {key: similarity_filtered_sequences[key] for key in blast_filter_proteins(similarity_filtered_sequences)}

    # final pass through and fasta file creation
    fasta, final_proteins = generate_full_complexes(similarity_filtered_sequences, original_structure, tag_name, tag_length, sort=True, sort_number=5)
    print(f'Final Amount of Tags: {len(similarity_filtered_sequences)-1}')
    return final_proteins

def collect_matches(files):
    """
    This function searches through each fasta file of sequences and pools
    the sequences into a set. We then check the union of each pairwise set,
    and fourwise set to see if there are any common sequences.
    """
    sets = []
    for fasta_file in files:
        sequences = set()
        for seq in SeqIO.parse(fasta_file, 'fasta'):
            if len(seq.seq) > 13:
                sequences.add(seq.seq[-13:])
            else:
                sequences.add(seq.seq)
        sets.append(sequences)
    matches = set()
    # get pairwise matches
    for i1, s1 in enumerate(sets):
        for i2, s2 in enumerate(sets):
            if i1 >= i2:
                continue
            ms = sets[i1].intersection(sets[i2])
            matches.update(ms)
    matches = list(matches)
    # now filter these down
    for seq in matches[:]:
        if 'W' in seq or 'C' in seq or 'K' in seq or 'M' in seq:
            matches.remove(seq)
            continue
        protein = PA(seq)
        charge_at_ph_7 = protein.charge_at_pH(7.4)
        if charge_at_ph_7 < -1.0 or charge_at_ph_7 > 1.0:
            matches.remove(seq)
            continue
        gravy_score = protein.gravy()
        if gravy_score >= 0.0:
            matches.remove(seq)
            continue
    print(f'Matches: {len(matches)}')
    print(matches)

    # save all tag structures in a different fasta file
    fasta_file_name = f'common_tags.fasta'
    new_fasta = []
    for index, sequence in enumerate(matches):
        # manually adding back in proline ends
        tag = "P" + sequence + "P"
        new_fasta.append('>%s\n%s' % (f'common_tag_{index}', tag))
    with open(fasta_file_name, 'w') as f:
        f.write('\n'.join(new_fasta))
    return matches


if __name__ == "__main__":
    files = ['high_temp_only_tag_designs.fa', 'fixed_residues_designs.fa',
             'only_tag_designs.fa']
    for index, file in enumerate(files):
        files[index] = 'tag_design_outputs/' + file
    collect_matches(files)
    # fixed residues output
    main('fixed_residues_designs.fa', 'fixed_residues_designs', 'fixed_residues_designs.fa')
    # fixed chain output
    main('fixed_chain_designs.fa', 'fixed_chain_designs', '6I2G.fa')
    # only tag output
    main('only_tag_designs.fa', 'only_tag_designs', '6I2G.fa')
    # high temperature only tag output
    main('high_temp_only_tag_designs.fa', 'high_temp_only_tag_designs', '6I2G.fa')
