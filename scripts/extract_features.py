import json
import pandas as pd
import tqdm
import os
import argparse
import numpy as np
import re


def find_cdr_index(full_chain_seq, cdr3_seq):
    assert (
        cdr3_seq in full_chain_seq
    ), f"CDR3 Sequence: {cdr3_seq} must be in full chain sequence"

    return full_chain_seq.find(cdr3_seq)


def parse_af3_results(
    summary_json_path,
    confidences_json_path,
    a_cdr3_start,
    a_cdr3_end,
    b_cdr3_start,
    b_cdr3_end,
):

    chains = ["A", "B", "M", "P"]

    results = {}
    with open(summary_json_path, "r") as jfile1:
        summary_confidence = json.load(jfile1)
        results["ptm"] = summary_confidence["ptm"]
        results["iptm"] = summary_confidence["iptm"]
        chain_pair_iptm = summary_confidence["chain_pair_iptm"]
        chain_ptm = summary_confidence["chain_ptm"]
        for i, chain1 in enumerate(chains):
            results[f"chain_ptm_{chain1}"] = chain_ptm[i]
            for j, chain2 in enumerate(chains):
                if chain1 < chain2:
                    results[f"chain_pair_iptm_{chain1}_{chain2}"] = chain_pair_iptm[i][
                        j
                    ]

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
                    results[f"avg_pae_interaction_{chain1}_{chain2}"] = 0.5 * (
                        np.mean(sub_pae1) + np.mean(sub_pae2)
                    )
                    results[f"min_pae_interaction_{chain1}_{chain2}"] = min(
                        np.min(sub_pae1), np.min(sub_pae2)
                    )
                    results[f"max_pae_interaction_{chain1}_{chain2}"] = max(
                        np.max(sub_pae1), np.max(sub_pae2)
                    )
                    sub_contact_probs = contact_probs[
                        residues_1[0] : residues_1[-1] + 1,
                        residues_2[0] : residues_2[-1] + 1,
                    ]
                    results[f"avg_contact_probs_{chain1}_{chain2}"] = np.mean(
                        sub_contact_probs
                    )
                    results[f"min_contact_probs_{chain1}_{chain2}"] = np.min(
                        sub_contact_probs
                    )
                    results[f"max_contact_probs_{chain1}_{chain2}"] = np.max(
                        sub_contact_probs
                    )
                elif i == j:
                    residues = [
                        int(idx) for idx in np.where(residue_chain_ids == chain1)[0]
                    ]
                    sub_pae = pae[
                        residues[0] : residues[-1] + 1, residues[0] : residues[-1] + 1
                    ]
                    results[f"avg_pae_{chain1}"] = np.mean(sub_pae)
                    results[f"min_pae_{chain1}"] = np.min(sub_pae)
                    results[f"max_pae_{chain1}"] = np.max(sub_pae)

                    atoms = [int(idx) for idx in np.where(atom_chain_ids == chain1)[0]]
                    sub_plddt = plddt[atoms[0] : atoms[-1] + 1]
                    results[f"avg_plddt_{chain1}"] = np.mean(sub_plddt)
                    results[f"min_plddt_{chain1}"] = np.min(sub_plddt)
                    results[f"max_plddt_{chain1}"] = np.max(sub_plddt)

        # CDR3 Metrics
        residues_alpha = [int(idx) for idx in np.where(residue_chain_ids == "A")[0]]
        residues_beta = [int(idx) for idx in np.where(residue_chain_ids == "B")[0]]
        residues_alpha_cdr3 = residues_alpha[a_cdr3_start:a_cdr3_end]
        residues_beta_cdr3 = residues_beta[b_cdr3_start:b_cdr3_end]
        for k, chain3 in enumerate(["M", "P"]):

            residues_3 = [int(idx) for idx in np.where(residue_chain_ids == chain3)[0]]

            sub_pae_a_chain = pae[
                residues_alpha_cdr3[0] : residues_alpha_cdr3[-1] + 1,
                residues_3[0] : residues_3[-1] + 1,
            ]
            sub_pae_chain_a = pae[
                residues_3[0] : residues_3[-1] + 1,
                residues_alpha_cdr3[0] : residues_alpha_cdr3[-1] + 1,
            ]
            pae_submatrix_a = np.concatenate(
                (sub_pae_a_chain, sub_pae_chain_a), axis=None
            )

            results[f"avg_pae_interaction_cdr3a_{chain3}"] = np.mean(pae_submatrix_a)
            results[f"min_pae_interaction_cdr3a_{chain3}"] = np.min(pae_submatrix_a)
            results[f"max_pae_interaction_cdr3a_{chain3}"] = np.max(pae_submatrix_a)
            results[f"std_pae_interaction_cdr3a_{chain3}"] = np.std(pae_submatrix_a)

            sub_pae_b_chain = pae[
                residues_beta_cdr3[0] : residues_beta_cdr3[-1] + 1,
                residues_3[0] : residues_3[-1] + 1,
            ]
            sub_pae_chain_b = pae[
                residues_3[0] : residues_3[-1] + 1,
                residues_beta_cdr3[0] : residues_beta_cdr3[-1] + 1,
            ]
            pae_submatrix_b = np.concatenate(
                (sub_pae_b_chain, sub_pae_chain_b), axis=None
            )
            results[f"avg_pae_interaction_cdr3b_{chain3}"] = np.mean(pae_submatrix_b)
            results[f"min_pae_interaction_cdr3b_{chain3}"] = np.min(pae_submatrix_b)
            results[f"max_pae_interaction_cdr3b_{chain3}"] = np.max(pae_submatrix_b)
            results[f"std_pae_interaction_cdr3b_{chain3}"] = np.std(pae_submatrix_b)

            sub_contact_probs_a_chain = contact_probs[
                residues_alpha_cdr3[0] : residues_alpha_cdr3[-1] + 1,
                residues_3[0] : residues_3[-1] + 1,
            ]
            results[f"avg_contact_probs_cdr3a_{chain3}"] = np.mean(
                sub_contact_probs_a_chain
            )
            results[f"max_contact_probs_cdr3a_{chain3}"] = np.max(
                sub_contact_probs_a_chain
            )

            sub_contact_probs_b_chain = contact_probs[
                residues_beta_cdr3[0] : residues_beta_cdr3[-1] + 1,
                residues_3[0] : residues_3[-1] + 1,
            ]
            results[f"avg_contact_probs_cdr3b_{chain3}"] = np.mean(
                sub_contact_probs_b_chain
            )
            results[f"max_contact_probs_cdr3b_{chain3}"] = np.max(
                sub_contact_probs_b_chain
            )
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
        ("cdr3a", args.cdr3alpha_col),
        ("cdr3b", args.cdr3beta_col),
        ("complex_id", args.id_col),
    ]:
        assert (
            col in seq_df.columns
        ), f"Missing {name} sequence column which was said to be set at {col}. Please add this column or change the column key using the appropriate flag."
        assert (
            seq_df[col].notna().all() and (seq_df[col] != "").all()
        ), f"{col} contains missing values"

    assert os.path.exists(args.af_output_dir), f"{args.af_output_dir} not found"

    try:
        model_seeds = [int(s) for s in args.af3_seeds.split(",")]
    except Exception as e:
        print(
            "Could not process the input for af3-seeds correctly. Please set it to a comma-separated-string of integers like 1,2,3"
        )
        raise e

    samples = [
        0,
        1,
        2,
        3,
        4,
    ]

    # Write the files
    rows = []
    for i, (idx, row) in tqdm.tqdm(enumerate(seq_df.iterrows())):

        name = row[args.id_col]
        name = name.lower()
        name = re.sub(r"[^a-z0-9-]+", "_", name)
        name = name.strip("_")
        path = os.path.join(args.af_output_dir, name)

        ranking = f"{path}/ranking_scores.csv"
        ranking_df = pd.read_csv(ranking, header=0)

        TRA = row[args.alpha_col]
        cdr3a = row[args.cdr3alpha_col]
        TRB = row[args.beta_col]
        cdr3b = row[args.cdr3beta_col]

        a_cdr3_start = find_cdr_index(TRA, cdr3a)
        a_cdr3_end = a_cdr3_start + len(cdr3a)

        b_cdr3_start = find_cdr_index(TRB, cdr3b)
        b_cdr3_end = b_cdr3_start + len(cdr3b)

        for seed in model_seeds:
            for sample in samples:
                new_row = {}
                new_row[args.id_col] = row[args.id_col]
                new_row["af3_output_directory"] = path
                new_row["af3_seed"] = seed
                new_row["af3_sample"] = sample
                summary_path = (
                    f"{path}/seed-{seed}_sample-{sample}/summary_confidences.json"
                )
                confidences_path = (
                    f"{path}/seed-{seed}_sample-{sample}/confidences.json"
                )
                sub_results = parse_af3_results(
                    summary_path,
                    confidences_path,
                    a_cdr3_start,
                    a_cdr3_end,
                    b_cdr3_start,
                    b_cdr3_end,
                )
                new_row.update(sub_results)
                ranking = ranking_df[
                    (ranking_df["seed"] == seed) & (ranking_df["sample"] == sample)
                ]["ranking_score"].values[0]
                new_row["af3_ranking"] = ranking
                rows.append(new_row)

    if not os.path.exists(args.output_dir):
        print(f"{args.output_dir} does not exist... Creating new directory")
        os.makedirs(args.output_dir)

    print(f"Writing results to {args.output_dir}...")

    group_cols = [args.id_col, "af3_output_directory"]

    all_results_df = pd.DataFrame(rows)
    avg_results_df = (
        all_results_df.groupby(by=group_cols, dropna=False)
        .mean()
        .drop(columns=["af3_seed", "af3_sample", "af3_ranking"])
        .reset_index()
    )
    best_idx = all_results_df.groupby(by=group_cols, dropna=False)[
        "af3_ranking"
    ].idxmax()
    best_results_df = (
        all_results_df.loc[best_idx]
        .drop(columns=["af3_seed", "af3_sample", "af3_ranking"])
        .reset_index()
        .drop(columns="index")
    )

    ensemble = all_results_df.groupby(by=group_cols, dropna=False).agg(["mean", "std"])
    ensemble.columns = [f"{col}_{stat}" for col, stat in ensemble.columns]
    columns_to_drop = [
        col
        for col in ensemble.columns
        if col.startswith(("af3_seed", "af3_sample", "af3_ranking"))
    ]
    ensemble = ensemble.drop(columns=columns_to_drop).reset_index()
    ensemble = ensemble.fillna(
        0
    )  # set std columns with no variance to zero with ddof=1

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
        "--af-output-dir",
        type=str,
        required=True,
        help="Directory with AlphaFold3 outputs. This directory should contain one subdirectory per row in the sequences-file",
    )

    parser.add_argument(
        "--af3-seeds",
        type=str,
        required=False,
        help="comma-separated string of seeds that were input to AF3",
        default="1,13,17,21,42,133,177,213,315,1001",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="Directory to place output CSV files with features",
    )

    parser.add_argument(
        "-i",
        "--id-col",
        type=str,
        required=False,
        default="complex_id",
        help="Name of column containing unique complex ID in <--sequences-file>",
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
        "-cdr3a",
        "--cdr3alpha-col",
        type=str,
        required=False,
        default="TRA_CDR3",
        help="Name of column containing CDR3a chain sequence in <--sequences-file>",
    )

    parser.add_argument(
        "-cdr3b",
        "--cdr3beta-col",
        type=str,
        required=False,
        default="TRB_CDR3",
        help="Name of column containing CDR3b chain sequence in <--sequences-file>",
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

    args = parser.parse_args()

    main(args)
