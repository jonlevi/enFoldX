import numpy as np
import pandas as pd
import os
import argparse
import mdtraj as md

##################################################################################

parser = argparse.ArgumentParser(
    description="Compute pairwise RMSD metrics of AF3-predicted ensemble.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

parser.add_argument(
    "-d",
    dest="data_dir",
    required=True,
    help="Path to directory containing AF3-predicted ensemble.",
)

parser.add_argument(
    "-o",
    dest="output",
    required=True,
    help="Path to directory to write output CSV.",
)

args = parser.parse_args()

##################################################################################


def find_cif_files(directory_path):
    """
    Find all .cif files in the given directory and its subdirectories.
    """
    cif_files = []
    basename = os.path.basename(directory_path)
    best_model = ""
    for root, _, files in os.walk(directory_path):
        for file in files:
            if file.endswith(".cif"):
                filepath = os.path.join(root, file)
                cif_files.append(filepath)
                if file.endswith("_model.cif"):
                    best_model = filepath

    if best_model in cif_files:
        cif_files.remove(best_model)
    return best_model, cif_files


def compute_rmsd(cif_path_1, cif_path_2):
    """Calculate pairwise RMSD metrics using Ca atoms:
        - overall RMSD after overall alignment (1 metric)
        - per-chain RMSD after overall alignment (4 metrics)
        - per-chain RMSD after per-chain alignment (16 metrics)

    Args:
        - cif_path_1 (str) : path to first structure (reference)
        - cif_path_2 (str) : path to second structure (target)
    Returns:
        - results_df (DataFrame) : all features for two structures
    """

    # Seed and sample numbers from file paths
    seed_sample1 = os.path.basename(os.path.dirname(cif_path_1))
    seed_sample2 = os.path.basename(os.path.dirname(cif_path_2))
    tcr_epitope = os.path.basename(os.path.dirname(os.path.dirname(cif_path_1)))

    # Load structures in MDtraj
    reference = md.load(cif_path_1)
    target = md.load(cif_path_2)

    # Overall RMSD with overall alignment
    ca_indices = reference.topology.select("name CA")
    overall_rmsd = md.rmsd(target, reference, 0, atom_indices=ca_indices)

    # Initialize results dictionary with RMSD value
    cif_results = {
        "tcr_epitope": tcr_epitope,
        "seed_sample": f"{seed_sample1}_{seed_sample2}",
        "overall_rmsd": float(overall_rmsd[0]),
    }

    chains = [c.index for c in reference.topology.chains]

    for chain_index1 in chains:

        # Select Ca atoms of current chain
        atom_indices = reference.topology.select(
            f"(chainid == {chain_index1}) and (name == CA)"
        )

        # Calculate chain RMSD with overall alignment
        rmsd_value = md.rmsd(
            target,
            reference,
            0,
            atom_indices=atom_indices,
            ref_atom_indices=atom_indices,
        )[0]
        cif_results[f"chain_{chain_index1}_rmsd"] = rmsd_value

        # Calculate chain RMSD with per-chain alignment
        for chain_index2 in chains:
            superpose_selection = reference.topology.select(f"chainid {chain_index2}")

            # Align on each chain
            target.superpose(
                reference,
                atom_indices=superpose_selection,
                ref_atom_indices=superpose_selection,
            )

            # Calulate RMSD in current chain
            rmsd_value = md.rmsd(
                target,
                reference,
                0,
                atom_indices=atom_indices,
                ref_atom_indices=atom_indices,
                superpose=False,
            )[0]
            cif_results[f"chain_{chain_index1}_{chain_index2}_aligned_rmsd"] = (
                rmsd_value
            )

    results_df = pd.DataFrame([cif_results])
    return results_df


##################################################################################

if __name__ == "__main__":
    os.makedirs(f"{args.output}", exist_ok=True)
    af3_output_path = args.data_dir

    data = []

    af3_input_name = os.path.basename(af3_output_path)
    if os.path.isfile(f"{args.output}/{af3_input_name}.csv"):
        print(f"Error: ensemble directory {af3_input_name} has already been checked.")
    else:
        ref_cif, cifs = find_cif_files(af3_output_path)
        num_samples = len(cifs)

        print(f"Computing RMSD metrics...")
        for i in range(num_samples):
            for j in range(num_samples)[i + 1 :]:
                cifs_results = compute_rmsd(cifs[i], cifs[j])
                data.append(cifs_results)

        keys = cifs_results.keys()
        ensemble_results = pd.concat(data, ignore_index=True)
        output_csv = f"{args.output}/{af3_input_name}.csv"
        ensemble_results.to_csv(output_csv, index=False)
        print(f"Completed: CSV with RMSD metrics written to output directory.")
