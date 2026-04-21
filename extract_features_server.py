import json
import pandas as pd
import tqdm
import os
import argparse
import numpy as np
import zipfile
from pathlib import Path


def unzip_file(zip_path):
    zip_path = Path(zip_path)

    output_dir = zip_path.parent

    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(output_dir)

    return output_dir


def parse_af3_results(summary_json_path, confidences_json_path):

    chains = ["A", "B", "C", "D"]
    chain_mapper = {"A": "A", "B": "B", "C": "P", "D": "M"}

    results = {}
    with open(summary_json_path, "r") as jfile1:
        summary_confidence = json.load(jfile1)
        results["af3_ranking"] = summary_confidence["ranking_score"]
        results["ptm"] = summary_confidence["ptm"]
        results["iptm"] = summary_confidence["iptm"]
        chain_pair_iptm = summary_confidence["chain_pair_iptm"]
        chain_ptm = summary_confidence["chain_ptm"]
        for i, chain1 in enumerate(chains):
            results[f"chain_ptm_{chain_mapper.get(chain1)}"] = chain_ptm[i]
            for j, chain2 in enumerate(chains):
                if chain1 < chain2:
                    results[
                        f"chain_pair_iptm_{chain_mapper.get(chain1)}_{chain_mapper.get(chain2)}"
                    ] = chain_pair_iptm[i][j]

    with open(confidences_json_path, "r") as jfile2:
        confidence = json.load(jfile2)
        plddt = np.array(confidence["atom_plddts"])
        results["avg_plddt"] = np.mean(plddt)
        pae = np.array(confidence["pae"])
        residue_chain_ids = np.array(confidence["token_chain_ids"])
        atom_chain_ids = np.array(confidence["atom_chain_ids"])
        results["avg_pae"] = np.mean(pae)
        contact_probs = np.array(confidence["contact_probs"])
        for i, chain1 in enumerate(chains):
            for j, chain2 in enumerate(chains):
                if i < j:
                    residues_1 = [
                        int(idx) for idx in np.where(residue_chain_ids == chain1)[0]
                    ]
                    residues_2 = [
                        int(idx) for idx in np.where(residue_chain_ids == chain2)[0]
                    ]
                    sub_pae1 = pae[
                        residues_1[0] : residues_1[-1] + 1,
                        residues_2[0] : residues_2[-1] + 1,
                    ]
                    sub_pae2 = pae[
                        residues_2[0] : residues_2[-1] + 1,
                        residues_1[0] : residues_1[-1] + 1,
                    ]
                    results[
                        f"avg_pae_interaction_{chain_mapper.get(chain1)}_{chain_mapper.get(chain2)}"
                    ] = 0.5 * (np.mean(sub_pae1) + np.mean(sub_pae2))
                    results[
                        f"min_pae_interaction_{chain_mapper.get(chain1)}_{chain_mapper.get(chain2)}"
                    ] = min(np.min(sub_pae1), np.min(sub_pae2))
                    results[
                        f"max_pae_interaction_{chain_mapper.get(chain1)}_{chain_mapper.get(chain2)}"
                    ] = max(np.max(sub_pae1), np.max(sub_pae2))
                    sub_contact_probs = contact_probs[
                        residues_1[0] : residues_1[-1] + 1,
                        residues_2[0] : residues_2[-1] + 1,
                    ]
                    results[
                        f"avg_contact_probs_{chain_mapper.get(chain1)}_{chain_mapper.get(chain2)}"
                    ] = np.mean(sub_contact_probs)
                    results[
                        f"min_contact_probs_{chain_mapper.get(chain1)}_{chain_mapper.get(chain2)}"
                    ] = np.min(sub_contact_probs)
                    results[
                        f"max_contact_probs_{chain_mapper.get(chain1)}_{chain_mapper.get(chain2)}"
                    ] = np.max(sub_contact_probs)
                elif i == j:
                    residues = [
                        int(idx) for idx in np.where(residue_chain_ids == chain1)[0]
                    ]
                    sub_pae = pae[
                        residues[0] : residues[-1] + 1, residues[0] : residues[-1] + 1
                    ]
                    results[f"avg_pae_{chain_mapper.get(chain1)}"] = np.mean(sub_pae)
                    results[f"min_pae_{chain_mapper.get(chain1)}"] = np.min(sub_pae)
                    results[f"max_pae_{chain_mapper.get(chain1)}"] = np.max(sub_pae)

                    atoms = [int(idx) for idx in np.where(atom_chain_ids == chain1)[0]]
                    sub_plddt = plddt[atoms[0] : atoms[-1] + 1]
                    results[f"avg_plddt_{chain_mapper.get(chain1)}"] = np.mean(
                        sub_plddt
                    )
                    results[f"min_plddt_{chain_mapper.get(chain1)}"] = np.min(sub_plddt)
                    results[f"max_plddt_{chain_mapper.get(chain1)}"] = np.max(sub_plddt)
    return results


def main(args):

    try:
        seq_df = pd.read_csv(args.sequences_file)
    except FileNotFoundError as e:
        print(f"Error: The file {args.sequences_file} was not found.")
    except Exception as e:
        print(f"Error reading CSV file: {e}")

    for name, col in [
        ("alpha", args.alpha_col),
        ("beta", args.beta_col),
        ("mhc", args.mhc_col),
        ("peptide", args.peptide_col),
    ]:
        assert (
            col in seq_df.columns
        ), f"Missing {name} sequence column which was to {col}. Please check it is set properly and try again."
        assert (
            seq_df[col].notna().all() and (seq_df[col] != "").all()
        ), f"{col} contains missing values"

    assert os.path.exists(args.zip_file), f"{rgs.zip_file} not found"

    if not os.path.exists(args.output_dir):
        print(f"{args.output_dir} does not exist... Creating new directory")
        os.makedirs(args.output_dir)

    samples = [
        0,
        1,
        2,
        3,
        4,
    ]

    # Write the files
    rows = []
    af3_dir = unzip_file(args.zip_file)
    for subdir in os.listdir(af3_dir):
        if not os.path.isdir(os.path.join(af3_dir, subdir)):
            continue

        new_row = {}

        input_json = os.path.join(af3_dir, subdir, f"fold_{subdir}_job_request.json")
        with open(input_json, "r") as ij:
            input_data = json.load(ij)[0]

        new_row["name"] = input_data["name"]
        new_row["af3_seed"] = int(input_data["modelSeeds"][0])

        chain_seqs = {}

        for i, chain in enumerate(["A", "B", "P", "M"]):

            chain_seqs[chain] = input_data["sequences"][i]["proteinChain"]["sequence"]

        original_df_row = seq_df[
            (seq_df[args.alpha_col] == chain_seqs["A"])
            & (seq_df[args.beta_col] == chain_seqs["B"])
            & (seq_df[args.mhc_col] == chain_seqs["M"])
            & (seq_df[args.peptide_col] == chain_seqs["P"])
        ]
        assert (
            original_df_row.shape[0] == 1
        ), f"Seq df must have exactly one row per input TCR-pMHC, but we found  {original_df_row.shape[0]} for {subdir}"

        for sample in samples:

            new_row["original_index"] = original_df_row.index.values[0]
            new_row["af3_sample"] = sample
            summary_path = os.path.join(
                af3_dir, subdir, f"fold_{subdir}_summary_confidences_{sample}.json"
            )
            confidences_path = os.path.join(
                af3_dir, subdir, f"fold_{subdir}_full_data_{sample}.json"
            )

            sub_results = parse_af3_results(summary_path, confidences_path)
            new_row.update(sub_results)
            rows.append(new_row)

    all_results_df = pd.DataFrame(rows)
    avg_results_df = (
        all_results_df.groupby(by=["original_index", "name"], dropna=False)
        .mean()
        .drop(columns=["af3_seed", "af3_sample", "af3_ranking"])
        .reset_index()
    )
    best_idx = all_results_df.groupby(by=["original_index"], dropna=False)[
        "af3_ranking"
    ].idxmax()
    best_results_df = (
        all_results_df.loc[best_idx]
        .drop(columns=["af3_seed", "af3_sample", "af3_ranking"])
        .reset_index()
        .drop(columns="index")
    )

    ensemble = all_results_df.groupby(by=["original_index", "name"], dropna=False).agg(
        ["mean", "std"]
    )
    ensemble.columns = [f"{col}_{stat}" for col, stat in ensemble.columns]
    columns_to_drop = [
        col
        for col in ensemble.columns
        if col.startswith(("af3_seed", "af3_sample", "af3_ranking"))
    ]
    ensemble = ensemble.drop(columns=columns_to_drop).reset_index()

    all_results_df.to_csv(os.path.join(args.output_dir, "all_structures_features.csv"))
    avg_results_df.to_csv(os.path.join(args.output_dir, "avg_features.csv"))
    best_results_df.to_csv(os.path.join(args.output_dir, "best_features.csv"))
    ensemble.to_csv(os.path.join(args.output_dir, "ensemble_features.csv"))

    print(f"Output features written to {args.output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-s",
        "--sequences-file",
        type=str,
        required=True,
        help="Path to input file with TCR-pMHC sequences",
    )

    parser.add_argument(
        "-a",
        "--alpha-col",
        type=str,
        required=False,
        default="TRA_aa",
        help="Name of column containing TCRa chain sequence in <--sequences-file>",
    )

    parser.add_argument(
        "-b",
        "--beta-col",
        type=str,
        required=False,
        default="TRB_aa",
        help="Name of column containing TCRb chain sequence in <--sequences-file>",
    )

    parser.add_argument(
        "-m",
        "--mhc-col",
        type=str,
        required=False,
        default="M_aa",
        help="Name of column containing MHC chain sequence in <--sequences-file>",
    )

    parser.add_argument(
        "-p",
        "--peptide-col",
        type=str,
        required=False,
        default="peptide",
        help="Name of column containing peptide sequence in <--sequences-file>",
    )

    parser.add_argument(
        "-z",
        "--zip-file",
        type=str,
        required=True,
        help="Zip with AlphaFold3 Server Outputs. This directory should contain one subdirectory per job downloaded",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="Directory to place output CSV files with features",
    )

    args = parser.parse_args()

    main(args)
