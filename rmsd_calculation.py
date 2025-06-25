import pandas as pd
import numpy as np
import os
import sys
from Bio.PDB import MMCIFParser, Superimposer

import tqdm


def three_to_one(resname):
    """Convert 3-letter amino acid code to 1-letter."""
    AA_3TO1 = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLU": "E",
        "GLN": "Q",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
    }
    return AA_3TO1[resname]


def get_sequence_and_residues(chain):
    seq = []
    residues = []
    for res in chain:
        if "CA" in res:
            try:
                aa = three_to_one(res.get_resname())
                seq.append(aa)
                residues.append(res)
            except KeyError:
                continue  # skip non-standard residues
    return seq, residues


def aligned_chain_rmsd(
    reference_structure_path,
    target_structure_path,
    reference_align_chain_id,
    target_align_chain_id,
    reference_rmsd_chain_id,
    target_rmsd_chain_id,
):
    """
    A pymol-style alignment where sequence alignment is first computed in specific chains.
    Then the CA of the sequence-aligned residues between the two PDBs are structurally aligned.
    RMSD is then computed over the specified chains.

    Args:
        reference_structure_path: path to reference .cif file.
        target_structure_path: path to target .cif file.
        reference_align_chain_id: chain ID in reference structure used for alignment.
        target_align_chain_id: chain ID in align structure used for alignment.
        reference_rmsd_chain_id: chain ID in reference used for RMSD calculation.
        target_rmsd_chain_id: chain ID in target structure used for RMSD calculation.
    Returns:
        RMSD (float)
    """
    parser = MMCIFParser(QUIET=False)
    ref_structure = parser.get_structure("ref", reference_structure_path)
    target_structure = parser.get_structure("align", target_structure_path)

    # Alignment chains (used for alignment and superposition)
    ref_align_chain = ref_structure[0][reference_align_chain_id]
    target_align_chain = target_structure[0][target_align_chain_id]

    # RMSD chains (used for final RMSD calculation)
    ref_rmsd_chain = ref_structure[0][reference_rmsd_chain_id]
    target_rmsd_chain = target_structure[0][target_rmsd_chain_id]

    # Sequence and residue extraction from alignment chains
    ref_seq, ref_residues = get_sequence_and_residues(ref_align_chain)
    target_seq, target_residues = get_sequence_and_residues(target_align_chain)

    ref_seq_str = "".join(ref_seq)
    target_seq_str = "".join(target_seq)

    atoms_ref = []
    atoms_target = []

    for ref_res, target_res in zip(ref_residues, target_residues):
        if "CA" in ref_res and "CA" in target_res:
            atoms_ref.append(ref_res["CA"])
            atoms_target.append(target_res["CA"])

    # Superimpose target_structure onto ref_structure
    sup = Superimposer()
    sup.set_atoms(atoms_ref, atoms_target)
    sup.apply(target_structure.get_atoms())

    ref_coords = [res["CA"].get_coord() for res in ref_rmsd_chain if "CA" in res]
    target_coords = [res["CA"].get_coord() for res in target_rmsd_chain if "CA" in res]

    if len(ref_coords) != len(target_coords):
        raise ValueError(
            "Number of C-alpha atoms in RMSD chains must match after filtering."
        )

    ref_coords = np.array(ref_coords)
    target_coords = np.array(target_coords)

    diff = ref_coords - target_coords
    rmsd = np.sqrt((diff**2).sum() / len(ref_coords))
    return round(rmsd, 2)


def compute_rmsd_between_stuctures(structure1, structure2):

    return {
        "rmsd_M_P": aligned_chain_rmsd(structure1, structure2, "M", "M", "P", "P"),
        "rmsd_P_P": aligned_chain_rmsd(structure1, structure2, "P", "P", "P", "P"),
        "rmsd_A_P": aligned_chain_rmsd(structure1, structure2, "A", "A", "P", "P"),
        "rmsd_B_P": aligned_chain_rmsd(structure1, structure2, "B", "B", "P", "P"),
        "rmsd_M_A": aligned_chain_rmsd(structure1, structure2, "M", "M", "A", "A"),
        "rmsd_M_B": aligned_chain_rmsd(structure1, structure2, "M", "M", "B", "B"),
    }


def main():

    s1 = "path/to/some/model1.cif"
    s2 = "path/to/some/model2.cif"

    rmsd_values = compute_rmsd_between_stuctures(f1, f2)
    print(rmsd_values)


if __name__ == "__main__":
    main()
