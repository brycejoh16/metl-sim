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


    df1.to_csv(join(out_dir,'replicate_1.csv'),index=True,index_label='variant')

    heatmap_df2=plot_weights_heatmap(df2, save_fn=join(out_dir, 'heatmap_replicate_2.png'),
                         wildtype_seq=wt_seq, positions_per_row=100)

    df2.to_csv(join(out_dir, 'replicate_2.csv'), index=True,index_label='variant')


