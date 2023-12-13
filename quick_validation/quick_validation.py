#!/usr/bin/env python
"""Quick validation of Protein structure Model before submitting to PDB"""
from Bio.Align import PairwiseAligner
import subprocess
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')


def get_fasta_sequence(fasta_file):
    """Get the sequence from the fasta file"""

    # Open the file
    with open(fasta_file, "r") as f:
        # Read the file
        lines = f.readlines()

    # Get the sequence
    sequence = lines[1].strip()

    return sequence

def get_pdb_sequence(model_file):
    """Get model sequence from PDB file"""

    pymol_get_seq = f"cmd.load('{model_file}');print('XXXXXXXX') ;print(cmd.get_fastastr('all'))"
    
    result = subprocess.run(["pymol", "-c", "-q", "-d", pymol_get_seq], capture_output=True)
    result = result.stdout.decode("utf-8").split("\n")
    index = result.index("XXXXXXXX")

    sequences = {}
    for i in range(index+1, len(result)):
        if result[i].startswith(">"):
            seq_name = result[i]
            sequences[seq_name] = ""
        else:
            sequences[seq_name] += result[i]

    return sequences



def check_sequence(reference_sequence, model_sequences):
    """Align Reference Sequence with Model Sequence and check for any mismatches"""

    # Create a PairwiseAligner object
    aligner = PairwiseAligner()
    errors = []

    # Perform alignment
    for model_name, model_sequence in model_sequences.items():

        alignments = aligner.align(reference_sequence, model_sequence)

        # Get the best alignment (first one in the list)
        best_alignment = alignments[0]


        ref_alignment = best_alignment[0]
        model_alignment = best_alignment[1]

        # Check for mismatches, ignoring gaps
        mismatches = [i for i in range(len(ref_alignment)) 
                      if ref_alignment[i] != model_alignment[i]  
                      and model_alignment[i] != '-']

        if mismatches:
            errors.append((model_name, mismatches, best_alignment))


    if errors:
        logging.error("Mismatches found in the following models:")
        for error in errors:
            model_name, mismatches, best_alignment = error
            logging.info(f"Mismatches in {model_name} at positions: {mismatches}\n{best_alignment}")

    else:
        logging.info("No mismatches found")

def close_contacts(model_file, cutoff=2.2):
    """Check if there are any close contacts in the model"""
    pymol_command = f"cmd.load('{model_file}'); print(cmd.distance('tmp', 'all and not inorganic', 'all and not inorganic', mode=2, cutoff={cutoff}))"
    result = subprocess.run(["pymol", "-c", "-q", "-d", pymol_command], check=True, capture_output=True, text=True)
    result = float(result.stdout.splitlines()[-1])
    
    if result > 0:
        logging.error(f"Polar contacts closer or equal to {cutoff} A found in the model.")
    else:
        logging.info(f"No polar contacts closer or equal to {cutoff} A found in the model.")

def molprobity(model_file, cif_files=None):
    """Run MolProbity validation"""
    if cif_files:
        cif_str = " ".join(str(cif) for cif in cif_files)
        ref_stats_command = f"phenix.model_statistics {model_file} {cif_str}".split()
    else:
        ref_stats_command = f"phenix.model_statistics {model_file}".split()
    
    try:
        ref_stats = subprocess.run(ref_stats_command, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running phenix.model_statistics. Did you include the cif files? Please check manually. Aborting.")
        exit(1)
    
    for index, line in enumerate(ref_stats.stdout.splitlines()): # Get all parameters from output
        if line.strip().startswith("Bond"):
            line_list = line.strip().split()
            logging.info(f"Bond RMSD: {line_list[2]}")
        if line.strip().startswith("Angle"):
            line_list = line.strip().split()
            logging.info(f"Angle RMSD: {line_list[2]}")
        if line.strip().startswith("Ramachandran Plot"):
            outl_index = index + 1
            allowed_index = index + 2
            favored_index = index + 3

            outl_list = ref_stats.stdout.splitlines()[outl_index].strip().split()
            logging.info(f"Ramachandran Outliers: {float(outl_list[2])}")

            allowed_list = ref_stats.stdout.splitlines()[allowed_index].strip().split()
            logging.info(f"Ramachandran Allowed: {float(allowed_list[2])}")

            favored_list = ref_stats.stdout.splitlines()[favored_index].strip().split()
            logging.info(f"Ramachandran Favored: {float(favored_list[2])}")

        if line.strip().startswith("Rotamer Outliers"):
            line_list = line.strip().split()
            logging.info(f"Rotamer Outliers: {float(line_list[3])}")
        if line.strip().startswith("All-atom Clashscore"):
            line_list = line.strip().split()
            logging.info(f"All-atom Clashscore: {float(line_list[-1])}")
        

def main(reference_fasta, model, cif_files=None):

    reference_sequence = get_fasta_sequence(reference_fasta)
    model_sequence = get_pdb_sequence(model)
    check_sequence(reference_sequence, model_sequence)
    close_contacts(model, cutoff=2.2)
    molprobity(model, cif_files=args.cif)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", help="Reference sequence in fasta format")
    parser.add_argument("model", help="Model file in PDB format")
    parser.add_argument("-c", "--cif", help="CIF files from refinement", required=False, nargs="+")
    args = parser.parse_args()

    main(args.reference, args.model, args.cif)