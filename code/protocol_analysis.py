from os.path import join

import pandas as pd
from Bio.PDB import PDBParser
from Bio.SeqIO.PdbIO import AtomIterator
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde,spearmanr
from matplotlib.patches import Rectangle

CHARS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
C2I_MAPPING = {c: i for i, c in enumerate(CHARS)}



plt.rcParams['font.family'] = 'serif'  # Change 'serif' to your desired font family
plt.rcParams['font.serif'] = ['Times New Roman']  # Example: specify a font
plt.rcParams['font.size'] = 18  # Adjust the font size as needed

def get_seq_from_pdb(pdb_fn):
    # do it the hacky way that avoids parsing the PDB file header which may not be present
    structure = PDBParser().get_structure(None, pdb_fn)
    records = [record for record in AtomIterator(None, structure)]

    if len(records) > 1:
        # all our FASTAs should have single records only
        raise ValueError("pdb file has more than one record: {}".format(pdb_fn))

    # the base sequence
    seq = str(records[0].seq)

    return seq




def parity_plot(x, y, x_label='x', y_label='y', save_path=None, title=None):
    """
    Create a density-colored parity plot from two pandas Series.
    Ensures tight x/y limits with a small margin so points are visible.
    Displays Spearman correlation coefficient on the plot.
    """
    # Drop missing values
    valid = ~(x.isna() | y.isna())
    print("Omitting: ", sum(~valid), " values from parity plot")
    x = x[valid]
    y = y[valid]

    # --- Compute Spearman correlation ---
    spearman_rho, _ = spearmanr(x, y)

    # --- Compute KDE density for coloring ---
    xy = np.vstack([x, y])
    kde = gaussian_kde(xy)
    density = kde(xy)
    hue_values = (density - density.min()) / (density.max() - density.min())

    # --- Set up JointGrid manually for better control ---
    g = sns.JointGrid(x=x, y=y, space=0, height=6)
    g.plot_joint(sns.scatterplot, c=hue_values, cmap='viridis', s=20, alpha=0.6, edgecolor='none')
    g.plot_marginals(sns.histplot, kde=True, color='gray')

    # --- Add parity line ---
    min_val = min(x.min(), y.min())
    max_val = max(x.max(), y.max())
    g.ax_joint.plot([min_val, max_val], [min_val, max_val], 'k--', lw=1)

    # --- Compute limits with small wiggle (~2%) ---
    buffer = 0.02 * (max_val - min_val)
    g.ax_joint.set_xlim(min_val - buffer, max_val + buffer)
    g.ax_joint.set_ylim(min_val - buffer, max_val + buffer)

    # --- Grid, labels, and title ---
    g.ax_joint.grid(True, color='lightgrey', linestyle='-', linewidth=0.5, alpha=0.7)
    g.ax_joint.set_xlabel(x_label, fontsize=14)
    g.ax_joint.set_ylabel(y_label, fontsize=14)
    g.ax_joint.tick_params(axis='both', labelsize=12)
    if title:
        g.ax_joint.set_title(title, fontsize=14)

    # --- Add Spearman correlation text block ---
    text_x = min_val + 0.05 * (max_val - min_val)
    text_y = max_val - 0.1 * (max_val - min_val)
    g.ax_joint.text(
        text_x, text_y,
        f"Spearman ρ = {spearman_rho:.3f}",
        fontsize=12, color='black',
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="gray", alpha=0.7)
    )

    # --- Tight layout and save ---
    plt.tight_layout()
    if save_path:
        g.fig.savefig(save_path, dpi=300, bbox_inches='tight')

    return g.fig





def custom_sort(variant):

    i , mt = int(variant[1:-1]),variant[-1]
    return i*20+CHARS.index(mt)
def plot_weights_heatmap(df, save_fn, wildtype_seq, positions_per_row=100):
    '''
    Final version:
    - Clean layout
    - No colorbar in main figure
    - Save colorbar separately as its own PNG
    '''
    df = df[['interface_delta_X']]

    wt_variants = [f'{wt}{t+1}{wt}' for t, wt in enumerate(wildtype_seq)]

    # Get the wildtype row
    wt_row = df.loc[['_wt']].copy()

    # Duplicate it for each variant
    wt_rows = pd.concat([wt_row] * len(wt_variants), axis=0)
    wt_rows.index = wt_variants  # assign the new variant names

    # Append to the dataframe
    df = pd.concat([df, wt_rows])

    df = df[df.index!='_wt']



    df_sorted = df.sort_index(key=lambda x:[custom_sort(i) for i in x])





    model_weights = df_sorted['interface_delta_X'].to_numpy()

    num_amino_acids = len(CHARS)
    total_weights = len(model_weights)
    seq_length = total_weights // num_amino_acids

    if len(wildtype_seq) != seq_length:
        raise ValueError("Wildtype sequence length doesn't match model weight dimensions.")

    weights_matrix = model_weights.reshape(seq_length, num_amino_acids).T

    num_chunks = (seq_length + positions_per_row - 1) // positions_per_row
    fig_width = positions_per_row * 0.3
    fig_height = num_chunks * 5

    fig, axes = plt.subplots(nrows=num_chunks, ncols=1,
                             figsize=(fig_width, fig_height),
                             sharex=False)

    if num_chunks == 1:
        axes = [axes]

    vmin = np.min(weights_matrix)
    vmax = np.max(weights_matrix)
    cmap = "vlag"

    for chunk_idx, ax in enumerate(axes):
        start = chunk_idx * positions_per_row
        end = min((chunk_idx + 1) * positions_per_row, seq_length)
        chunk_weights = weights_matrix[:, start:end]
        width = end - start

        sns.heatmap(chunk_weights,
                    cmap=cmap,
                    cbar=False,
                    yticklabels=CHARS,
                    xticklabels=False,
                    ax=ax,
                    linewidths=0.5,
                    linecolor='white',
                    vmin=vmin,
                    vmax=vmax)

        # Wildtype boxes
        for pos in range(start, end):
            aa_index = C2I_MAPPING[wildtype_seq[pos]]
            rect = Rectangle((pos - start, aa_index), 1, 1,
                             fill=False, edgecolor='black', linewidth=1.2)
            ax.add_patch(rect)

        # X-tick labels
        tick_positions = np.arange(width)
        tick_labels = [
            f'{pos + start + 1}' if (pos + start) % 5 == 0 else ''
            for pos in tick_positions
        ]
        ax.set_xticks(tick_positions + 0.5)
        ax.set_xticklabels(tick_labels, rotation=45, ha='center')
        ax.set_ylabel("Amino Acid")
        ax.set_xlim(0, positions_per_row)

    fig.tight_layout()
    fig.savefig(save_fn, dpi=300)
    plt.close(fig)

    # Create and save separate colorbar figure
    fig_cb, ax_cb = plt.subplots(figsize=(1.5, 6))
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = fig_cb.colorbar(sm, cax=ax_cb, orientation='vertical')
    cbar.set_label("Weight")
    fig_cb.tight_layout()
    fig_cb.savefig(f"{''.join(save_fn.split('.')[:-1])}_legend.png", dpi=300)
    plt.close(fig_cb)


    return df_sorted



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

def load_mutational_data(df):
    """Load mutational effect scores from the CSV."""
    df['position'] = df.index.to_series().apply(lambda x: int(x[1:-1]))
    avg_effects = (df.groupby('position')['interface_delta_X'].mean()-df['interface_delta_X'].mean()).to_dict()
    return avg_effects

def assign_b_factors(residues, avg_effects):
    """Assign the average mutational effect to the B-factor field of each residue."""
    for chain_id, pos, residue in residues:
        effect = avg_effects.get(pos, 0.0)  # Default to 0.0 if no data
        for atom in residue:
            atom.bfactor = effect # Store effect in the B-factor field
    return residues


from Bio import PDB


def save_new_pdb(structure, output_file):
    """Save the modified PDB with new B-factors and a PyMOL remark."""
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_file)

    # Insert remark at the top of the file
    with open(output_file, "r+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write("REMARK  spectrum b, blue_white_red, minimum=-5, maximum=5\n" + content)


def save_new_pdb_image(out_fn, structure):
    from pymol import cmd

    # Start PyMOL without GUI
    cmd.reinitialize()

    # Load your structure
    cmd.load(structure)

    # Apply your spectrum command
    cmd.spectrum("b", "blue_white_red", minimum=-5, maximum=5)

    # Set a nice representation
    cmd.show("cartoon")
    cmd.bg_color("white")

    # Zoom and save the image
    cmd.zoom()
    cmd.ray(1200, 900)
    cmd.png(out_fn, width=1200, height=900, dpi=300, ray=1)

    # Optional: quit PyMOL (only needed in batch scripts)
    # cmd.quit()


def docked_on_struct_heatmap(pdb_file, df, output_pdb):
    # Parse the PDB and load mutational data
    residues = parse_pdb(pdb_file)
    avg_effects = load_mutational_data(df)

    # Assign B-factors based on the mutational effect
    residues=assign_b_factors(residues, avg_effects)

    # Save the new PDB with modified B-factors
    structure = residues[0][2].get_parent().get_parent()  # Get the top structure object
    save_new_pdb(structure, output_pdb)
    print(f"Modified PDB saved as: {output_pdb}")




def make_single_variant_replicate_analysis(out_dir):
    df=pd.read_csv(join(out_dir,'energies_df.csv'))
    pdb_file =join('pdb_files','prepared_pdb_files',df['pdb_fn'][0])
    wt_seq= get_seq_from_pdb(pdb_file)

    assert len(df)==(len(wt_seq)*19*2)+2, 'this should be 2 replicates of single mutants, (19*N*2) plus ' \
                                          'two wild type variants _wt.'

    # Group by variant
    groups = df.groupby("variant")

    # Split each group into two
    df1 = []
    df2 = []

    for _, group in groups:
        first, second = group.iloc[0], group.iloc[1]
        df1.append(first)
        df2.append(second)

    # Convert lists back to dataframes
    df1 = pd.DataFrame(df1).set_index("variant",drop=True)
    df2 = pd.DataFrame(df2).set_index("variant",drop=True)
    assert len(df1) == len(df2) == len(df) / 2

    #ordering should be the same but just in case
    df1['interface_delta_X_2'] = df2['interface_delta_X']



    parity_plot(df1['interface_delta_X'],df1['interface_delta_X_2'],x_label='Replicate 1',y_label='Replicate 2',
                # title=f"interface_delta_X across replicates for protocol\n"
                #       f"{'_'.join(out_dir.split('/')[-2].split('_')[4:])}",
                save_path=join(out_dir,'replicate1_vs_replicate2.png'))

    heatmap_df1= plot_weights_heatmap(df1, save_fn=join(out_dir,'heatmap_replicate_1.png'),
                        wildtype_seq=wt_seq, positions_per_row=100)

    out_pdb=join(out_dir, "replicate_1_bsplines_docked_on_struct.pdb")
    docked_on_struct_heatmap(pdb_file, heatmap_df1, output_pdb=out_pdb)
    save_new_pdb_image(join(out_dir, "replicate_1_bsplines_docked_on_struct.png"), out_pdb)

    df1.to_csv(join(out_dir,'replicate_1.csv'),index=True,index_label='variant')

    heatmap_df2=plot_weights_heatmap(df2, save_fn=join(out_dir, 'heatmap_replicate_2.png'),
                         wildtype_seq=wt_seq, positions_per_row=100)

    out_pdb=join(out_dir, "replicate_2_bspline_docked_on_struct.pdb")
    docked_on_struct_heatmap(pdb_file, heatmap_df2, output_pdb=out_pdb)
    save_new_pdb_image(join(out_dir,"replicate_2_bsplines_docked_on_struct.png"), out_pdb)

    df2.to_csv(join(out_dir, 'replicate_2.csv'), index=True,index_label='variant')


