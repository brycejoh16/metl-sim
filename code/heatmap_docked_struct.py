

import pandas as pd
from Bio import PDB
import re

def parse_pdb(pdb_file):
    """Parse the PDB and return residues with chain and position."""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue, standard=True):
                    residues.append((chain.id, residue.id[1], residue))
    return residues

def load_mutational_data(csv_file):
    """Load mutational effect scores from the CSV."""
    df = pd.read_csv(csv_file)
    df['position'] = df['variant'].apply(lambda x: int(re.findall(r'\d+', x)[0]))
    avg_effects = (df.groupby('position')['interface_delta_X'].mean()+15).to_dict()
    return avg_effects

def assign_b_factors(residues, avg_effects):
    """Assign the average mutational effect to the B-factor field of each residue."""
    for chain_id, pos, residue in residues:
        effect = avg_effects.get(pos, 0.0)  # Default to 0.0 if no data
        for atom in residue:
            atom.bfactor = effect  # Store effect in the B-factor field

def save_new_pdb(structure, output_file):
    """Save the modified PDB with new B-factors."""
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_file)

def main(pdb_file, csv_file, output_pdb):
    # Parse the PDB and load mutational data
    residues = parse_pdb(pdb_file)
    avg_effects = load_mutational_data(csv_file)

    # Assign B-factors based on the mutational effect
    assign_b_factors(residues, avg_effects)

    # Save the new PDB with modified B-factors
    structure = residues[0][2].get_parent().get_parent()  # Get the top structure object
    save_new_pdb(structure, output_pdb)

    print(f"Modified PDB saved as: {output_pdb}")


if __name__ == '__main__':

    # Example usage
    pdb_filename = "pdb_files/prepared_pdb_files/SadA_rosetta_2024_3_6_p.pdb"  # Replace with your PDB path
    csv_filename = "output/htcondor_runs/condor_energize_2024-09-23_00-16-21_sadA_singles_andres_docking/processed_run/energies_df.csv"  # Replace with your CSV path
    output_pdb = "output/htcondor_runs/condor_energize_2024-09-23_00-16-21_sadA_singles_andres_docking/processed_run/docked_with_mutants.pdb"  # Output path for modified PDB

    main(pdb_filename, csv_filename, output_pdb)
