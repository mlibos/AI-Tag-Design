#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=32g
#SBATCH --gres=gpu:rtx2080:1
#SBATCH -c 2
#SBATCH --output=example_2.out

# ensure you have created the boltz_design env
source activate boltz_design

# this is an example for running bd1 on a sample tag
# increased helix loss, 25 samples, increased binder properties, no alphafold or trajectory and included an iptm cutoff value

# here are all the optional tags you can use while running bd1
#  boltzdesign.py [-h] --target_name TARGET_NAME [--target_type {protein,rna,dna,small_molecule,metal}] [--input_type {pdb,custom}] [--pdb_path PDB_PATH]
#                     [--pdb_target_ids PDB_TARGET_IDS] [--target_mols TARGET_MOLS] [--custom_target_input CUSTOM_TARGET_INPUT]
#                     [--custom_target_ids CUSTOM_TARGET_IDS] [--binder_id BINDER_ID] [--use_msa USE_MSA] [--msa_max_seqs MSA_MAX_SEQS] [--suffix SUFFIX]
#                     [--modifications MODIFICATIONS] [--modifications_wt MODIFICATIONS_WT] [--modifications_positions MODIFICATIONS_POSITIONS]
#                     [--modification_target MODIFICATION_TARGET] [--constraint_target CONSTRAINT_TARGET] [--contact_residues CONTACT_RESIDUES]
#                     [--length_min LENGTH_MIN] [--length_max LENGTH_MAX] [--optimizer_type {SGD,AdamW}] [--pre_iteration PRE_ITERATION]
#                     [--soft_iteration SOFT_ITERATION] [--temp_iteration TEMP_ITERATION] [--hard_iteration HARD_ITERATION]
#                     [--semi_greedy_steps SEMI_GREEDY_STEPS] [--recycling_steps RECYCLING_STEPS] [--use_default_config USE_DEFAULT_CONFIG]
#                     [--mask_ligand MASK_LIGAND] [--optimize_contact_per_binder_pos OPTIMIZE_CONTACT_PER_BINDER_POS] [--distogram_only DISTOGRAM_ONLY]
#                     [--design_algorithm {3stages,3stages_extra}] [--learning_rate LEARNING_RATE] [--learning_rate_pre LEARNING_RATE_PRE] [--e_soft E_SOFT]
#                     [--e_soft_1 E_SOFT_1] [--e_soft_2 E_SOFT_2] [--inter_chain_cutoff INTER_CHAIN_CUTOFF] [--intra_chain_cutoff INTRA_CHAIN_CUTOFF]
#                     [--num_inter_contacts NUM_INTER_CONTACTS] [--num_intra_contacts NUM_INTRA_CONTACTS] [--con_loss CON_LOSS] [--i_con_loss I_CON_LOSS]
#                     [--plddt_loss PLDDT_LOSS] [--pae_loss PAE_LOSS] [--i_pae_loss I_PAE_LOSS] [--rg_loss RG_LOSS] [--helix_loss_max HELIX_LOSS_MAX]
#                     [--helix_loss_min HELIX_LOSS_MIN] [--num_designs NUM_DESIGNS] [--cutoff CUTOFF] [--i_ptm_cutoff I_PTM_CUTOFF]
#                     [--complex_plddt_cutoff COMPLEX_PLDDT_CUTOFF] [--gpu_id GPU_ID] [--design_samples DESIGN_SAMPLES] [--work_dir WORK_DIR]
#                     [--high_iptm HIGH_IPTM] [--boltz_checkpoint BOLTZ_CHECKPOINT] [--ccd_path CCD_PATH] [--alphafold_dir ALPHAFOLD_DIR]
#                     [--af3_docker_name AF3_DOCKER_NAME] [--af3_database_settings AF3_DATABASE_SETTINGS] [--af3_hmmer_path AF3_HMMER_PATH]
#                     [--run_boltz_design RUN_BOLTZ_DESIGN] [--run_ligandmpnn RUN_LIGANDMPNN] [--run_alphafold RUN_ALPHAFOLD] [--run_rosetta RUN_ROSETTA]
#                     [--redo_boltz_predict REDO_BOLTZ_PREDICT] [--show_animation SHOW_ANIMATION] [--save_trajectory SAVE_TRAJECTORY]


python ../boltzdesign.py \
        --target_name tag_8829 \ 
        --pdb_path ./inputs/tags_for_design/tag_8829.pdb \
        --target_type protein \
        --pdb_target_ids A \
        --design_samples 25 \
        --run_alphafold False \
        --length_min 50 \
        --length_max 80 \
        --num_intra_contacts 6 \
        --num_inter_contacts 4 \
        --helix_loss_min "-0.8" \
        --helix_loss_max "-0.4" \
        --save_trajectory False \
        --i_ptm_cutoff "0.85" 
