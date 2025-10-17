import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

def read_csv_file_and_group(file_path, trial_filter=None):
    """Read CSV file, group by mutation info, and get interface_delta_X from row with minimum total_score for each group

    Args:
        file_path: Path to CSV file
        trial_filter: Dict with 'min_trial' and 'max_trial' to filter by trial index, or None for all trials
    """
    try:
        # Read the CSV file
        df = pd.read_csv(file_path)

        # Filter by trial index if specified
        if trial_filter is not None:
            min_trial = trial_filter.get('min_trial', 1)
            max_trial = trial_filter.get('max_trial', float('inf'))
            if 'trial' in df.columns:
                df = df[(df['trial'] >= min_trial) & (df['trial'] <= max_trial)]
                print(f"  - Filtered to trials {min_trial}-{max_trial}: {len(df)} rows")
            else:
                print("  - Warning: 'trial_idx' column not found, using all data")

        # Group by Mutation_point and Mutate_into, then get the row with minimum total_score for each group
        grouped = df.groupby(['Mutation_point', 'Mutate_into']).apply(
            lambda x: x.loc[x['interface_delta_X'].idxmin()]
        ).reset_index(drop=True)

        # Store both interface_delta_X and total_score for each mutation group (no filtering)
        data = {}
        total_count = len(grouped)

        for _, row in grouped.iterrows():
            description = f"{row['Mutation_point']}_{row['Mutate_into']}"
            interface_delta_x = row['interface_delta_X']
            total_score = row['total_score']

            data[description] = {
                'interface_delta_X': interface_delta_x,
                'total_score': total_score
            }

        print(f"  - Total mutation groups: {total_count}")
        print(f"  - All groups included (no filtering): {len(data)}")

        return data

    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return {}

def match_data(data1, data2):
    """Match data between two dictionaries based on description keys"""
    matched_data = {}
    for description in data1:
        if description in data2:
            matched_data[description] = {
                'file1': data1[description],
                'file2': data2[description]
            }
    return matched_data

def calculate_correlation_stats(x, y):
    """Calculate correlation statistics"""
    # Pearson correlation coefficient
    pearson_r, pearson_p = stats.pearsonr(x, y)
    
    # Spearman correlation coefficient
    spearman_r, spearman_p = stats.spearmanr(x, y)
    
    # R-squared value
    r_squared = pearson_r**2
    
    # Mean squared error from y=x line
    mse = np.mean((x - y)**2)
    rmse = np.sqrt(mse)
    
    # Mean absolute error from y=x line
    mae = np.mean(np.abs(x - y))
    
    return {
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'r_squared': r_squared,
        'rmse': rmse,
        'mae': mae
    }

def create_scatter_plot_and_stats(matched_data, output_prefix, data_type='interface_delta_X'):
    """Create scatter plot and calculate statistics for specified data type"""
    # Extract the matched values for the specified data type
    x_values = [matched_data[desc]['file1'][data_type] for desc in matched_data]
    y_values = [matched_data[desc]['file2'][data_type] for desc in matched_data]
    
    # Convert to numpy arrays for statistical calculations
    x_values = np.array(x_values)
    y_values = np.array(y_values)
    
    # Calculate correlation statistics
    stats_dict = calculate_correlation_stats(x_values, y_values)
    
    # Create the scatter plot
    plt.figure(figsize=(10, 8))
    plt.scatter(x_values, y_values, alpha=0.6, s=50)
    
    # Set labels based on data type
    if data_type == 'interface_delta_X':
        plt.xlabel('interface_delta_X (first trial)', fontsize=12)
        plt.ylabel('interface_delta_X (second trial)', fontsize=12)
        # Set fixed axis limits from -20 to 0 for interface_delta_X
        plt.xlim(-20, 0)
        plt.ylim(-20, 0)
        # Add y=x line in red from -20 to 0
        plt.plot([-20, 0], [-20, 0], 'r-', linewidth=2)
        # Position text at top left of plot area (fixed coordinates for -20 to 0 range)
        x_pos = -19.6
        y_pos = -0.4
    else:
        plt.xlabel('total_score (first trial)', fontsize=12)
        plt.ylabel('total_score (second trial)', fontsize=12)
        # Auto-scale for total_score
        min_val = min(np.min(x_values), np.min(y_values))
        max_val = max(np.max(x_values), np.max(y_values))
        margin = (max_val - min_val) * 0.05
        plt.xlim(min_val - margin, max_val + margin)
        plt.ylim(min_val - margin, max_val + margin)
        # Add y=x line for the auto-scaled range
        plt.plot([min_val, max_val], [min_val, max_val], 'r-', linewidth=2)
        # Position text at top left of auto-scaled plot area
        x_pos = min_val + (max_val - min_val) * 0.02
        y_pos = max_val - (max_val - min_val) * 0.02

    # Add statistics as text
    stats_text = f"""RMSE = {stats_dict['rmse']:.4f}
N = {len(matched_data)}"""

    plt.text(x_pos, y_pos, stats_text, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.grid(True, alpha=0.3)
    
    # Save the plot with high DPI
    plt.tight_layout()
    plot_filename = f'{output_prefix}_{data_type}_correlation_plot.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.close()  # Close to free memory
    
    return stats_dict, plot_filename

def save_detailed_stats(stats_dict, matched_data, output_prefix, data_type='interface_delta_X'):
    """Save detailed statistics to a text file"""
    stats_filename = f'{output_prefix}_{data_type}_detailed_stats.txt'
    
    with open(stats_filename, 'w') as f:
        f.write("Detailed Correlation Statistics\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Pearson correlation coefficient: {stats_dict['pearson_r']:.6f}\n")
        f.write(f"Pearson p-value: {stats_dict['pearson_p']:.2e}\n")
        f.write(f"Spearman correlation coefficient: {stats_dict['spearman_r']:.6f}\n")
        f.write(f"Spearman p-value: {stats_dict['spearman_p']:.2e}\n")
        f.write(f"R-squared: {stats_dict['r_squared']:.6f}\n")
        f.write(f"Root Mean Square Error (RMSE): {stats_dict['rmse']:.6f}\n")
        f.write(f"Mean Absolute Error (MAE): {stats_dict['mae']:.6f}\n")
        f.write(f"Number of matched data points: {len(matched_data)}\n")
        
        # Add some interpretation
        f.write("\nInterpretation:\n")
        f.write("-" * 20 + "\n")
        if abs(stats_dict['pearson_r']) > 0.8:
            f.write("Strong correlation detected.\n")
        elif abs(stats_dict['pearson_r']) > 0.5:
            f.write("Moderate correlation detected.\n")
        else:
            f.write("Weak correlation detected.\n")
            
        if stats_dict['pearson_p'] < 0.001:
            f.write("Correlation is highly significant (p < 0.001).\n")
        elif stats_dict['pearson_p'] < 0.05:
            f.write("Correlation is significant (p < 0.05).\n")
        else:
            f.write("Correlation is not statistically significant (p >= 0.05).\n")
    
    return stats_filename

def calculate_ddg_and_create_bar_graph(data, output_prefix, data_type='total_score'):
    """Calculate ddG (each position's value - lowest value) and create bar graph"""

    # Parse mutation data to extract position and mutation info
    position_data = {}
    for description, values in data.items():
        # Extract position and amino acid from description like "74_ALA"
        parts = description.split('_')
        position = int(parts[0])
        amino_acid = parts[1]

        if position not in position_data:
            position_data[position] = {}

        position_data[position][amino_acid] = values[data_type]

    # Calculate ddG for each position (value - minimum value across all mutations for that position)
    ddg_data = {}
    position_summary = {}

    for position in position_data:
        values = list(position_data[position].values())
        min_value = min(values)

        position_summary[position] = {
            'min_value': min_value,
            'ddg_values': []
        }

        for amino_acid, value in position_data[position].items():
            ddg = value - min_value
            ddg_key = f"{position}_{amino_acid}"
            ddg_data[ddg_key] = ddg
            position_summary[position]['ddg_values'].append(ddg)

    # Prepare data for bar graph - group by position
    positions = sorted(position_summary.keys())
    x_labels = []
    ddg_values = []
    colors = []

    # Color scheme: different colors for different amino acids
    aa_colors = {
        'ALA': 'red', 'ARG': 'blue', 'ASN': 'green', 'ASP': 'orange', 'CYS': 'purple',
        'GLN': 'brown', 'GLU': 'pink', 'GLY': 'gray', 'HIS': 'olive', 'ILE': 'cyan',
        'LEU': 'magenta', 'LYS': 'yellow', 'MET': 'lime', 'PHE': 'navy', 'PRO': 'teal',
        'SER': 'maroon', 'THR': 'silver', 'TRP': 'gold', 'TYR': 'indigo', 'VAL': 'coral'
    }

    for position in positions:
        for amino_acid in sorted(position_data[position].keys()):
            x_labels.append(f"{position}\n{amino_acid}")
            ddg_key = f"{position}_{amino_acid}"
            ddg_values.append(ddg_data[ddg_key])
            colors.append(aa_colors.get(amino_acid, 'black'))

    # Create bar graph
    plt.figure(figsize=(40, 10))
    bars = plt.bar(range(len(x_labels)), ddg_values, color=colors, alpha=0.7)

    # Customize the plot
    plt.xlabel('Mutation Position and Amino Acid', fontsize=14)
    if data_type == 'total_score':
        plt.ylabel('ddG (Total Score - Min Total Score)', fontsize=14)
        plot_title = 'ddG Values by Mutation Position (Total Score)'
    else:
        plt.ylabel('ddG (Interface Delta Score - Min Interface Delta Score)', fontsize=14)
        plot_title = 'ddG Values by Mutation Position (Interface Delta Score)'

    plt.title(plot_title, fontsize=16)
    plt.xticks(range(len(x_labels)), x_labels, rotation=45, ha='right', fontsize=8)
    plt.grid(True, alpha=0.3, axis='y')

    # Add horizontal line at y=0
    plt.axhline(y=0, color='black', linestyle='-', linewidth=1)

    # Save the plot
    plt.tight_layout()
    plot_filename = f'{output_prefix}_{data_type}_ddg_bar_graph.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.close()

    # Save ddG data to CSV
    ddg_df = pd.DataFrame([
        {
            'Position': int(key.split('_')[0]),
            'Amino_Acid': key.split('_')[1],
            'ddG': value,
            'Original_Value': data[key][data_type]
        }
        for key, value in ddg_data.items()
    ])

    csv_filename = f'{output_prefix}_{data_type}_ddg_data.csv'
    ddg_df.to_csv(csv_filename, index=False)

    print(f"ddG bar graph saved as: {plot_filename}")
    print(f"ddG data saved as: {csv_filename}")

    return plot_filename, csv_filename

def create_per_position_bar_graph(data, output_prefix, data_type='total_score', stat_type='mean'):
    """Create per-position bar graph using mean or median values across amino acids"""

    # Parse mutation data to extract position and mutation info
    position_data = {}
    for description, values in data.items():
        # Extract position and amino acid from description like "74_ALA"
        parts = description.split('_')
        position = int(parts[0])
        amino_acid = parts[1]

        if position not in position_data:
            position_data[position] = {}

        position_data[position][amino_acid] = values[data_type]

    # Calculate ddG for each position first (value - minimum value across all mutations for that position)
    position_ddg_data = {}
    position_stats = {}

    for position in position_data:
        values = list(position_data[position].values())
        min_value = min(values)

        # Calculate ddG values for this position
        ddg_values = [value - min_value for value in values]

        if stat_type == 'mean':
            stat_value = np.mean(ddg_values)
        else:  # median
            stat_value = np.median(ddg_values)

        position_stats[position] = {
            'stat_value': stat_value,
            'ddg_values': ddg_values,
            'amino_acids': list(position_data[position].keys())
        }

    # Prepare data for bar graph
    positions = sorted(position_stats.keys())
    stat_values = [position_stats[pos]['stat_value'] for pos in positions]

    # Create bar graph
    plt.figure(figsize=(20, 8))
    bars = plt.bar(positions, stat_values, alpha=0.7, color='steelblue')

    # Customize the plot
    plt.xlabel('Position', fontsize=14)
    if data_type == 'total_score':
        plt.ylabel(f'{stat_type.title()} ddG (Total Score - Min Total Score)', fontsize=14)
        plot_title = f'Per-Position {stat_type.title()} ddG Values (Total Score)'
    else:
        plt.ylabel(f'{stat_type.title()} ddG (Interface Delta Score - Min Interface Delta Score)', fontsize=14)
        plot_title = f'Per-Position {stat_type.title()} ddG Values (Interface Delta Score)'

    plt.title(plot_title, fontsize=16)
    plt.grid(True, alpha=0.3, axis='y')

    # Add horizontal line at y=0
    plt.axhline(y=0, color='black', linestyle='-', linewidth=1)

    # Save the plot
    plt.tight_layout()
    plot_filename = f'{output_prefix}_{data_type}_per_position_{stat_type}_bar_graph.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.close()

    # Save position statistics to CSV
    position_df = pd.DataFrame([
        {
            'Position': position,
            f'{stat_type.title()}_ddG': position_stats[position]['stat_value'],
            'Num_Mutations': len(position_stats[position]['ddg_values'])
        }
        for position in positions
    ])

    csv_filename = f'{output_prefix}_{data_type}_per_position_{stat_type}_data.csv'
    position_df.to_csv(csv_filename, index=False)

    print(f"Per-position {stat_type} ddG bar graph saved as: {plot_filename}")
    print(f"Per-position {stat_type} ddG data saved as: {csv_filename}")

    return plot_filename, csv_filename, position_stats

def create_boxplots_top_positions(data, output_prefix, data_type='total_score', n_top=20):
    """Create box plots for top N positions with highest mean ddG energies"""

    # Parse mutation data to extract position and mutation info
    position_data = {}
    for description, values in data.items():
        # Extract position and amino acid from description like "74_ALA"
        parts = description.split('_')
        position = int(parts[0])
        amino_acid = parts[1]

        if position not in position_data:
            position_data[position] = {}

        position_data[position][amino_acid] = values[data_type]

    # Calculate ddG for each position and get mean values
    position_stats = {}
    for position in position_data:
        values = list(position_data[position].values())
        min_value = min(values)

        # Calculate ddG values for this position
        ddg_values = [value - min_value for value in values]
        mean_ddg = np.mean(ddg_values)

        position_stats[position] = {
            'mean_ddg': mean_ddg,
            'ddg_values': ddg_values,
            'amino_acids': list(position_data[position].keys())
        }

    # Sort positions by mean ddG and get top N
    sorted_positions = sorted(position_stats.items(), key=lambda x: x[1]['mean_ddg'], reverse=True)
    top_positions = sorted_positions[:n_top]

    # Prepare data for box plots
    boxplot_data = []
    position_labels = []

    for position, stats in top_positions:
        boxplot_data.append(stats['ddg_values'])
        position_labels.append(str(position))

    # Create box plot
    plt.figure(figsize=(16, 10))
    box_plot = plt.boxplot(boxplot_data, labels=position_labels, patch_artist=True)

    # Color the boxes
    colors = plt.cm.viridis(np.linspace(0, 1, len(boxplot_data)))
    for patch, color in zip(box_plot['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Customize the plot
    plt.xlabel('Position', fontsize=14)
    if data_type == 'total_score':
        plt.ylabel('ddG (Total Score - Min Total Score)', fontsize=14)
        plot_title = f'Box Plots for Top {n_top} Positions with Highest Mean ddG (Total Score)'
    else:
        plt.ylabel('ddG (Interface Delta Score - Min Interface Delta Score)', fontsize=14)
        plot_title = f'Box Plots for Top {n_top} Positions with Highest Mean ddG (Interface Delta Score)'

    plt.title(plot_title, fontsize=16)
    plt.xticks(rotation=45)
    plt.grid(True, alpha=0.3, axis='y')

    # Add horizontal line at y=0
    plt.axhline(y=0, color='black', linestyle='-', linewidth=1)

    # Save the plot
    plt.tight_layout()
    plot_filename = f'{output_prefix}_{data_type}_top_{n_top}_positions_boxplot.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.close()

    # Save top positions data to CSV
    top_positions_df = pd.DataFrame([
        {
            'Rank': i+1,
            'Position': position,
            'Mean_ddG': stats['mean_ddg'],
            'Std_ddG': np.std(stats['ddg_values']),
            'Min_ddG': np.min(stats['ddg_values']),
            'Max_ddG': np.max(stats['ddg_values']),
            'Num_Mutations': len(stats['ddg_values'])
        }
        for i, (position, stats) in enumerate(top_positions)
    ])

    csv_filename = f'{output_prefix}_{data_type}_top_{n_top}_positions_data.csv'
    top_positions_df.to_csv(csv_filename, index=False)

    print(f"Top {n_top} positions box plot saved as: {plot_filename}")
    print(f"Top {n_top} positions data saved as: {csv_filename}")

    return plot_filename, csv_filename

def get_wild_type_mapping():
    """Get wild type amino acid mapping from PDB file"""
    wild_type_map = {}

    # Wild type sequence from the PDB file - position: amino acid
    wild_type_sequence = {
        1: 'MET', 2: 'GLN', 3: 'HIS', 4: 'THR', 5: 'TYR', 6: 'PRO', 7: 'ALA', 8: 'GLN', 9: 'LEU', 10: 'MET',
        11: 'ARG', 12: 'PHE', 13: 'GLY', 14: 'THR', 15: 'ALA', 16: 'ALA', 17: 'ARG', 18: 'ALA', 19: 'GLU', 20: 'HIS',
        21: 'MET', 22: 'THR', 23: 'ILE', 24: 'ALA', 25: 'ALA', 26: 'ALA', 27: 'ILE', 28: 'GLY', 29: 'ALA', 30: 'LEU',
        31: 'GLY', 32: 'ALA', 33: 'ASP', 34: 'GLU', 35: 'ALA', 36: 'ASP', 37: 'ALA', 38: 'ILE', 39: 'VAL', 40: 'MET',
        41: 'ASP', 42: 'ILE', 43: 'VAL', 44: 'PRO', 45: 'ASP', 46: 'GLY', 47: 'GLU', 48: 'ARG', 49: 'ASP', 50: 'ALA',
        51: 'TRP', 52: 'TRP', 53: 'ASP', 54: 'ASP', 55: 'GLU', 56: 'GLY', 57: 'PHE', 58: 'SER', 59: 'SER', 60: 'SER',
        61: 'PRO', 62: 'PHE', 63: 'THR', 64: 'LYS', 65: 'ASN', 66: 'ALA', 67: 'HIS', 68: 'HIS', 69: 'ALA', 70: 'GLY',
        71: 'ILE', 72: 'VAL', 73: 'ALA', 74: 'THR', 75: 'SER', 76: 'VAL', 77: 'THR', 78: 'LEU', 79: 'GLY', 80: 'GLN',
        81: 'LEU', 82: 'GLN', 83: 'ARG', 84: 'GLU', 85: 'GLN', 86: 'GLY', 87: 'ASP', 88: 'LYS', 89: 'LEU', 90: 'VAL',
        91: 'SER', 92: 'LYS', 93: 'ALA', 94: 'ALA', 95: 'GLU', 96: 'TYR', 97: 'PHE', 98: 'GLY', 99: 'ILE', 100: 'ALA',
        101: 'CYS', 102: 'ARG', 103: 'VAL', 104: 'ASN', 105: 'ASP', 106: 'GLY', 107: 'LEU', 108: 'ARG', 109: 'THR', 110: 'THR',
        111: 'ARG', 112: 'PHE', 113: 'VAL', 114: 'ARG', 115: 'LEU', 116: 'PHE', 117: 'SER', 118: 'ASP', 119: 'ALA', 120: 'LEU',
        121: 'ASP', 122: 'ALA', 123: 'LYS', 124: 'PRO', 125: 'LEU', 126: 'THR', 127: 'ILE', 128: 'GLY', 129: 'HIS', 130: 'ASP',
        131: 'TYR', 132: 'GLU', 133: 'VAL', 134: 'GLU', 135: 'PHE', 136: 'LEU', 137: 'LEU', 138: 'ALA', 139: 'THR', 140: 'ARG',
        141: 'ARG', 142: 'VAL', 143: 'TYR', 144: 'GLU', 145: 'PRO', 146: 'PHE', 147: 'GLU', 148: 'ALA', 149: 'PRO', 150: 'PHE',
        151: 'ASN', 152: 'PHE', 153: 'ALA', 154: 'PRO', 155: 'HIS', 156: 'CYS', 157: 'GLY', 158: 'ASP', 159: 'VAL', 160: 'SER',
        161: 'TYR', 162: 'GLY', 163: 'ARG', 164: 'ASP', 165: 'THR', 166: 'VAL', 167: 'ASN', 168: 'TRP', 169: 'PRO', 170: 'LEU',
        171: 'LYS', 172: 'ARG', 173: 'SER', 174: 'PHE', 175: 'PRO', 176: 'ARG', 177: 'GLN', 178: 'LEU', 179: 'GLY', 180: 'GLY',
        181: 'PHE', 182: 'LEU', 183: 'THR', 184: 'ILE', 185: 'GLN', 186: 'GLY', 187: 'ALA', 188: 'ASP', 189: 'ASN', 190: 'ASP',
        191: 'ALA', 192: 'GLY', 193: 'MET', 194: 'VAL', 195: 'MET', 196: 'TRP', 197: 'ASP', 198: 'ASN', 199: 'ARG', 200: 'PRO',
        201: 'GLU', 202: 'SER', 203: 'ARG', 204: 'ALA', 205: 'ALA', 206: 'LEU', 207: 'ASP', 208: 'GLU', 209: 'MET', 210: 'HIS',
        211: 'ALA', 212: 'GLU', 213: 'TYR', 214: 'ARG', 215: 'GLU', 216: 'THR', 217: 'GLY', 218: 'ALA', 219: 'ILE', 220: 'ALA',
        221: 'ALA', 222: 'LEU', 223: 'GLU', 224: 'ARG', 225: 'ALA', 226: 'ALA', 227: 'LYS', 228: 'ILE', 229: 'MET', 230: 'LEU',
        231: 'LYS', 232: 'PRO', 233: 'GLN', 234: 'PRO', 235: 'GLY', 236: 'GLN', 237: 'LEU', 238: 'THR', 239: 'LEU', 240: 'PHE',
        241: 'GLN', 242: 'SER', 243: 'LYS', 244: 'ASN', 245: 'LEU', 246: 'HIS', 247: 'ALA', 248: 'ILE', 249: 'GLU', 250: 'ARG',
        251: 'CYS', 252: 'THR', 253: 'SER', 254: 'THR', 255: 'ARG', 256: 'ARG', 257: 'THR', 258: 'MET', 259: 'GLY', 260: 'LEU',
        261: 'PHE', 262: 'LEU', 263: 'ILE', 264: 'HIS', 265: 'THR', 266: 'GLU', 267: 'ASP', 268: 'GLY', 269: 'TRP', 270: 'ARG',
        271: 'MET', 272: 'PHE', 273: 'ASP'
    }

    return wild_type_sequence

def find_wild_type_energies(data, data_type='total_score'):
    """Find wild type energies for each position"""
    wild_type_map = get_wild_type_mapping()
    wild_type_energies = {}

    for description, values in data.items():
        # Extract position and amino acid from description like "74_ALA"
        parts = description.split('_')
        position = int(parts[0])
        amino_acid = parts[1]

        # Check if this is wild type for this position
        if position in wild_type_map and wild_type_map[position] == amino_acid:
            wild_type_energies[position] = values[data_type]

    return wild_type_energies

def calculate_ddg_wt_and_create_bar_graph(data, output_prefix, data_type='total_score'):
    """Calculate ddG using wild type as reference and create bar graph"""

    # Get wild type energies
    wild_type_energies = find_wild_type_energies(data, data_type)

    # Parse mutation data to extract position and mutation info
    position_data = {}
    for description, values in data.items():
        # Extract position and amino acid from description like "74_ALA"
        parts = description.split('_')
        position = int(parts[0])
        amino_acid = parts[1]

        if position not in position_data:
            position_data[position] = {}

        position_data[position][amino_acid] = values[data_type]

    # Calculate ddG for each position (value - wild type value for that position)
    ddg_data = {}
    position_summary = {}

    for position in position_data:
        if position not in wild_type_energies:
            print(f"Warning: No wild type energy found for position {position}")
            continue

        wt_value = wild_type_energies[position]
        position_summary[position] = {
            'wt_value': wt_value,
            'ddg_values': []
        }

        for amino_acid, value in position_data[position].items():
            ddg = value - wt_value
            ddg_key = f"{position}_{amino_acid}"
            ddg_data[ddg_key] = ddg
            position_summary[position]['ddg_values'].append(ddg)

    # Prepare data for bar graph - group by position
    positions = sorted(position_summary.keys())
    x_labels = []
    ddg_values = []
    colors = []

    # Color scheme: different colors for different amino acids
    aa_colors = {
        'ALA': 'red', 'ARG': 'blue', 'ASN': 'green', 'ASP': 'orange', 'CYS': 'purple',
        'GLN': 'brown', 'GLU': 'pink', 'GLY': 'gray', 'HIS': 'olive', 'ILE': 'cyan',
        'LEU': 'magenta', 'LYS': 'yellow', 'MET': 'lime', 'PHE': 'navy', 'PRO': 'teal',
        'SER': 'maroon', 'THR': 'silver', 'TRP': 'gold', 'TYR': 'indigo', 'VAL': 'coral'
    }

    for position in positions:
        for amino_acid in sorted(position_data[position].keys()):
            x_labels.append(f"{position}\n{amino_acid}")
            ddg_key = f"{position}_{amino_acid}"
            ddg_values.append(ddg_data[ddg_key])
            colors.append(aa_colors.get(amino_acid, 'black'))

    # Create bar graph
    plt.figure(figsize=(40, 10))
    bars = plt.bar(range(len(x_labels)), ddg_values, color=colors, alpha=0.7)

    # Customize the plot
    plt.xlabel('Mutation Position and Amino Acid', fontsize=14)
    if data_type == 'total_score':
        plt.ylabel('ddG (Total Score - Wild Type Total Score)', fontsize=14)
        plot_title = 'ddG Values by Mutation Position vs Wild Type (Total Score)'
    else:
        plt.ylabel('ddG (Interface Delta Score - Wild Type Interface Delta Score)', fontsize=14)
        plot_title = 'ddG Values by Mutation Position vs Wild Type (Interface Delta Score)'

    plt.title(plot_title, fontsize=16)
    plt.xticks(range(len(x_labels)), x_labels, rotation=45, ha='right', fontsize=8)
    plt.grid(True, alpha=0.3, axis='y')

    # Add horizontal line at y=0
    plt.axhline(y=0, color='black', linestyle='-', linewidth=1)

    # Save the plot
    plt.tight_layout()
    plot_filename = f'{output_prefix}_{data_type}_ddg_wt_bar_graph.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.close()

    # Save ddG data to CSV
    ddg_df = pd.DataFrame([
        {
            'Position': int(key.split('_')[0]),
            'Amino_Acid': key.split('_')[1],
            'ddG_vs_WT': value,
            'Original_Value': data[key][data_type],
            'WT_Value': wild_type_energies.get(int(key.split('_')[0]), 'N/A')
        }
        for key, value in ddg_data.items()
    ])

    csv_filename = f'{output_prefix}_{data_type}_ddg_wt_data.csv'
    ddg_df.to_csv(csv_filename, index=False)

    print(f"Wild type ddG bar graph saved as: {plot_filename}")
    print(f"Wild type ddG data saved as: {csv_filename}")

    return plot_filename, csv_filename, ddg_data

def create_per_position_bar_graph_wt(data, output_prefix, data_type='total_score', stat_type='mean'):
    """Create per-position bar graph using mean or median values across amino acids (wild type reference)"""

    # Get wild type energies
    wild_type_energies = find_wild_type_energies(data, data_type)

    # Parse mutation data to extract position and mutation info
    position_data = {}
    for description, values in data.items():
        # Extract position and amino acid from description like "74_ALA"
        parts = description.split('_')
        position = int(parts[0])
        amino_acid = parts[1]

        if position not in position_data:
            position_data[position] = {}

        position_data[position][amino_acid] = values[data_type]

    # Calculate ddG for each position first (value - wild type value for that position)
    position_stats = {}

    for position in position_data:
        if position not in wild_type_energies:
            continue

        wt_value = wild_type_energies[position]

        # Calculate ddG values for this position
        ddg_values = [value - wt_value for value in position_data[position].values()]

        if stat_type == 'mean':
            stat_value = np.mean(ddg_values)
        else:  # median
            stat_value = np.median(ddg_values)

        position_stats[position] = {
            'stat_value': stat_value,
            'ddg_values': ddg_values,
            'amino_acids': list(position_data[position].keys())
        }

    # Prepare data for bar graph
    positions = sorted(position_stats.keys())
    stat_values = [position_stats[pos]['stat_value'] for pos in positions]

    # Create bar graph
    plt.figure(figsize=(20, 8))
    bars = plt.bar(positions, stat_values, alpha=0.7, color='steelblue')

    # Customize the plot
    plt.xlabel('Position', fontsize=14)
    if data_type == 'total_score':
        plt.ylabel(f'{stat_type.title()} ddG (Total Score - Wild Type Total Score)', fontsize=14)
        plot_title = f'Per-Position {stat_type.title()} ddG Values vs Wild Type (Total Score)'
    else:
        plt.ylabel(f'{stat_type.title()} ddG (Interface Delta Score - Wild Type Interface Delta Score)', fontsize=14)
        plot_title = f'Per-Position {stat_type.title()} ddG Values vs Wild Type (Interface Delta Score)'

    plt.title(plot_title, fontsize=16)
    plt.grid(True, alpha=0.3, axis='y')

    # Add horizontal line at y=0
    plt.axhline(y=0, color='black', linestyle='-', linewidth=1)

    # Save the plot
    plt.tight_layout()
    plot_filename = f'{output_prefix}_{data_type}_per_position_{stat_type}_wt_bar_graph.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.close()

    # Save position statistics to CSV
    position_df = pd.DataFrame([
        {
            'Position': position,
            f'{stat_type.title()}_ddG_vs_WT': position_stats[position]['stat_value'],
            'Num_Mutations': len(position_stats[position]['ddg_values'])
        }
        for position in positions
    ])

    csv_filename = f'{output_prefix}_{data_type}_per_position_{stat_type}_wt_data.csv'
    position_df.to_csv(csv_filename, index=False)

    print(f"Per-position {stat_type} ddG vs WT bar graph saved as: {plot_filename}")
    print(f"Per-position {stat_type} ddG vs WT data saved as: {csv_filename}")

    return plot_filename, csv_filename, position_stats

def create_boxplots_top_positions_wt(data, output_prefix, data_type='total_score', n_top=20, lowest=True):
    """Create box plots for top N positions with highest/lowest mean ddG energies vs wild type"""

    # Get wild type energies
    wild_type_energies = find_wild_type_energies(data, data_type)

    # Parse mutation data to extract position and mutation info
    position_data = {}
    for description, values in data.items():
        # Extract position and amino acid from description like "74_ALA"
        parts = description.split('_')
        position = int(parts[0])
        amino_acid = parts[1]

        if position not in position_data:
            position_data[position] = {}

        position_data[position][amino_acid] = values[data_type]

    # Calculate ddG for each position and get mean values
    position_stats = {}
    for position in position_data:
        if position not in wild_type_energies:
            continue

        wt_value = wild_type_energies[position]

        # Calculate ddG values for this position
        ddg_values = [value - wt_value for value in position_data[position].values()]
        mean_ddg = np.mean(ddg_values)

        position_stats[position] = {
            'mean_ddg': mean_ddg,
            'ddg_values': ddg_values,
            'amino_acids': list(position_data[position].keys())
        }

    # Sort positions by mean ddG and get top N
    sorted_positions = sorted(position_stats.items(), key=lambda x: x[1]['mean_ddg'], reverse=not lowest)
    top_positions = sorted_positions[:n_top]

    # Prepare data for box plots
    boxplot_data = []
    position_labels = []

    for position, stats in top_positions:
        boxplot_data.append(stats['ddg_values'])
        position_labels.append(str(position))

    # Create box plot
    plt.figure(figsize=(16, 10))
    box_plot = plt.boxplot(boxplot_data, labels=position_labels, patch_artist=True)

    # Color the boxes
    colors = plt.cm.viridis(np.linspace(0, 1, len(boxplot_data)))
    for patch, color in zip(box_plot['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Customize the plot
    plt.xlabel('Position', fontsize=14)
    if data_type == 'total_score':
        plt.ylabel('ddG (Total Score - Wild Type Total Score)', fontsize=14)
        title_type = "Lowest" if lowest else "Highest"
        plot_title = f'Box Plots for Top {n_top} Positions with {title_type} Mean ddG vs Wild Type (Total Score)'
    else:
        plt.ylabel('ddG (Interface Delta Score - Wild Type Interface Delta Score)', fontsize=14)
        title_type = "Lowest" if lowest else "Highest"
        plot_title = f'Box Plots for Top {n_top} Positions with {title_type} Mean ddG vs Wild Type (Interface Delta Score)'

    plt.title(plot_title, fontsize=16)
    plt.xticks(rotation=45)
    plt.grid(True, alpha=0.3, axis='y')

    # Add horizontal line at y=0
    plt.axhline(y=0, color='black', linestyle='-', linewidth=1)

    # Save the plot
    plt.tight_layout()
    suffix = "lowest" if lowest else "highest"
    plot_filename = f'{output_prefix}_{data_type}_{suffix}_{n_top}_positions_wt_boxplot.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    plt.close()

    # Save top positions data to CSV
    top_positions_df = pd.DataFrame([
        {
            'Rank': i+1,
            'Position': position,
            'Mean_ddG_vs_WT': stats['mean_ddg'],
            'Std_ddG': np.std(stats['ddg_values']),
            'Min_ddG': np.min(stats['ddg_values']),
            'Max_ddG': np.max(stats['ddg_values']),
            'Num_Mutations': len(stats['ddg_values'])
        }
        for i, (position, stats) in enumerate(top_positions)
    ])

    csv_filename = f'{output_prefix}_{data_type}_{suffix}_{n_top}_positions_wt_data.csv'
    top_positions_df.to_csv(csv_filename, index=False)

    print(f"Top {n_top} {suffix} positions vs WT box plot saved as: {plot_filename}")
    print(f"Top {n_top} {suffix} positions vs WT data saved as: {csv_filename}")

    return plot_filename, csv_filename

def print_mutations_better_than_wt(data, output_prefix):
    """Print all mutations that are better than wild type for both total_score and interface_delta_X"""

    # Get wild type energies
    wt_total_scores = find_wild_type_energies(data, 'total_score')
    wt_interface_deltas = find_wild_type_energies(data, 'interface_delta_X')

    better_total_score = []
    better_interface_delta = []

    for description, values in data.items():
        # Extract position and amino acid from description like "74_ALA"
        parts = description.split('_')
        position = int(parts[0])
        amino_acid = parts[1]

        # Skip if this is wild type
        wild_type_map = get_wild_type_mapping()
        if position in wild_type_map and wild_type_map[position] == amino_acid:
            continue

        # Check if better than wild type for total_score (lower is better)
        if position in wt_total_scores:
            wt_score = wt_total_scores[position]
            if values['total_score'] < wt_score:
                ddg = values['total_score'] - wt_score
                better_total_score.append({
                    'mutation': description,
                    'total_score': values['total_score'],
                    'wt_total_score': wt_score,
                    'ddG': ddg
                })

        # Check if better than wild type for interface_delta_X (lower is better)
        if position in wt_interface_deltas:
            wt_interface = wt_interface_deltas[position]
            if values['interface_delta_X'] < wt_interface:
                ddg = values['interface_delta_X'] - wt_interface
                better_interface_delta.append({
                    'mutation': description,
                    'interface_delta_X': values['interface_delta_X'],
                    'wt_interface_delta_X': wt_interface,
                    'ddG': ddg
                })

    # Sort by ddG (most negative first - most improvement)
    better_total_score = sorted(better_total_score, key=lambda x: x['ddG'])
    better_interface_delta = sorted(better_interface_delta, key=lambda x: x['ddG'])

    print(f"\n=== MUTATIONS BETTER THAN WILD TYPE ===")

    print(f"\nMutations with better (lower) total_score than wild type ({len(better_total_score)} found):")
    print("-" * 80)
    print(f"{'Mutation':<15} {'Total Score':<12} {'WT Score':<12} {'ddG':<12} {'Improvement':<12}")
    print("-" * 80)
    for mutation in better_total_score:
        improvement = abs(mutation['ddG'])
        print(f"{mutation['mutation']:<15} {mutation['total_score']:<12.4f} {mutation['wt_total_score']:<12.4f} {mutation['ddG']:<12.4f} {improvement:<12.4f}")

    print(f"\nMutations with better (lower) interface_delta_X than wild type ({len(better_interface_delta)} found):")
    print("-" * 80)
    print(f"{'Mutation':<15} {'Interface ΔX':<12} {'WT Interface':<12} {'ddG':<12} {'Improvement':<12}")
    print("-" * 80)
    for mutation in better_interface_delta:
        improvement = abs(mutation['ddG'])
        print(f"{mutation['mutation']:<15} {mutation['interface_delta_X']:<12.4f} {mutation['wt_interface_delta_X']:<12.4f} {mutation['ddG']:<12.4f} {improvement:<12.4f}")

    # Save to files
    if better_total_score:
        total_df = pd.DataFrame(better_total_score)
        total_csv = f'{output_prefix}_mutations_better_total_score.csv'
        total_df.to_csv(total_csv, index=False)
        print(f"\nBetter total_score mutations saved to: {total_csv}")

    if better_interface_delta:
        interface_df = pd.DataFrame(better_interface_delta)
        interface_csv = f'{output_prefix}_mutations_better_interface_delta.csv'
        interface_df.to_csv(interface_csv, index=False)
        print(f"Better interface_delta_X mutations saved to: {interface_csv}")

    return better_total_score, better_interface_delta

def print_top_worst_mutations(data, n=10):
    """Print top N mutations with worst (largest) total_score and interface_delta_X values"""

    # Convert data to list of tuples for sorting
    mutations_list = []
    for description, values in data.items():
        mutations_list.append({
            'mutation': description,
            'total_score': values['total_score'],
            'interface_delta_X': values['interface_delta_X']
        })

    # Sort by total_score (descending - worst first)
    worst_total_score = sorted(mutations_list, key=lambda x: x['total_score'], reverse=True)[:n]

    # Sort by interface_delta_X (descending - worst first)
    worst_interface_delta = sorted(mutations_list, key=lambda x: x['interface_delta_X'], reverse=True)[:n]

    print(f"\n=== TOP {n} WORST MUTATIONS ===")

    print(f"\nTop {n} mutations with worst (largest) total_score:")
    print("-" * 60)
    print(f"{'Rank':<4} {'Mutation':<15} {'Total Score':<15}")
    print("-" * 60)
    for i, mutation in enumerate(worst_total_score, 1):
        print(f"{i:<4} {mutation['mutation']:<15} {mutation['total_score']:<15.4f}")

    print(f"\nTop {n} mutations with worst (largest) interface_delta_X:")
    print("-" * 60)
    print(f"{'Rank':<4} {'Mutation':<15} {'Interface Delta X':<18}")
    print("-" * 60)
    for i, mutation in enumerate(worst_interface_delta, 1):
        print(f"{i:<4} {mutation['mutation']:<15} {mutation['interface_delta_X']:<18.4f}")

def main():
    # Read and process the Output_100st_deep_div_lig folder for bar graphs
    '''
    print("Reading and processing Output_500st_deep/docking_results_parallel.csv for bar graphs...")
    data_deep = read_csv_file_and_group('Output_500st_deep_2_best_structures/docking_results_parallel.csv')

    if len(data_deep) == 0:
        print("No data found in Output_100st_deep!")
        return
    
    # Create bar graphs for ddG analysis (original per-mutation graphs)
    print("Creating ddG bar graph for total_score...")
    total_plot_file, total_csv_file = calculate_ddg_and_create_bar_graph(data_deep, 'deep_ddg_analysis', 'total_score')

    print("Creating ddG bar graph for interface_delta_X...")
    interface_plot_file, interface_csv_file = calculate_ddg_and_create_bar_graph(data_deep, 'deep_ddg_analysis', 'interface_delta_X')

    # Create per-position bar graphs using mean values
    print("Creating per-position mean ddG bar graph for total_score...")
    total_mean_plot, total_mean_csv, _ = create_per_position_bar_graph(data_deep, 'deep_ddg_analysis', 'total_score', 'mean')

    print("Creating per-position median ddG bar graph for total_score...")
    total_median_plot, total_median_csv, _ = create_per_position_bar_graph(data_deep, 'deep_ddg_analysis', 'total_score', 'median')

    print("Creating per-position mean ddG bar graph for interface_delta_X...")
    interface_mean_plot, interface_mean_csv, _ = create_per_position_bar_graph(data_deep, 'deep_ddg_analysis', 'interface_delta_X', 'mean')

    print("Creating per-position median ddG bar graph for interface_delta_X...")
    interface_median_plot, interface_median_csv, _ = create_per_position_bar_graph(data_deep, 'deep_ddg_analysis', 'interface_delta_X', 'median')

    # Create box plots for top 20 positions with highest mean energies
    print("Creating box plot for top 20 positions (total_score)...")
    total_boxplot, total_boxplot_csv = create_boxplots_top_positions(data_deep, 'deep_ddg_analysis', 'total_score', 20)

    print("Creating box plot for top 20 positions (interface_delta_X)...")
    interface_boxplot, interface_boxplot_csv = create_boxplots_top_positions(data_deep, 'deep_ddg_analysis', 'interface_delta_X', 20)

    # Create wild type-based ddG analysis
    print("\n=== WILD TYPE-BASED ddG ANALYSIS ===")

    # Create bar graphs for ddG analysis vs wild type
    print("Creating wild type ddG bar graph for total_score...")
    wt_total_plot_file, wt_total_csv_file, _ = calculate_ddg_wt_and_create_bar_graph(data_deep, 'deep_wt_ddg_analysis', 'total_score')

    print("Creating wild type ddG bar graph for interface_delta_X...")
    wt_interface_plot_file, wt_interface_csv_file, _ = calculate_ddg_wt_and_create_bar_graph(data_deep, 'deep_wt_ddg_analysis', 'interface_delta_X')

    # Create per-position bar graphs using mean values vs wild type
    print("Creating per-position mean ddG bar graph vs wild type for total_score...")
    wt_total_mean_plot, wt_total_mean_csv, _ = create_per_position_bar_graph_wt(data_deep, 'deep_wt_ddg_analysis', 'total_score', 'mean')

    print("Creating per-position median ddG bar graph vs wild type for total_score...")
    wt_total_median_plot, wt_total_median_csv, _ = create_per_position_bar_graph_wt(data_deep, 'deep_wt_ddg_analysis', 'total_score', 'median')

    print("Creating per-position mean ddG bar graph vs wild type for interface_delta_X...")
    wt_interface_mean_plot, wt_interface_mean_csv, _ = create_per_position_bar_graph_wt(data_deep, 'deep_wt_ddg_analysis', 'interface_delta_X', 'mean')

    print("Creating per-position median ddG bar graph vs wild type for interface_delta_X...")
    wt_interface_median_plot, wt_interface_median_csv, _ = create_per_position_bar_graph_wt(data_deep, 'deep_wt_ddg_analysis', 'interface_delta_X', 'median')

    # Create box plots for top 20 positions with lowest mean energies vs wild type
    print("Creating box plot for top 20 lowest positions vs wild type (total_score)...")
    wt_total_boxplot, wt_total_boxplot_csv = create_boxplots_top_positions_wt(data_deep, 'deep_wt_ddg_analysis', 'total_score', 20, lowest=True)

    print("Creating box plot for top 20 lowest positions vs wild type (interface_delta_X)...")
    wt_interface_boxplot, wt_interface_boxplot_csv = create_boxplots_top_positions_wt(data_deep, 'deep_wt_ddg_analysis', 'interface_delta_X', 20, lowest=True)
    '''
    # Read and process the Output_500st_deep data split by trial index for correlation analysis
    print("\nReading and processing Output_500st_deep/docking_results_parallel.csv for correlation (trial 1-250)...")
    data1 = read_csv_file_and_group('docking_results_90st_temp0.csv', {'min_trial': 1, 'max_trial': 45})

    print("Reading and processing Output_500st_deep/docking_results_parallel.csv for correlation (trial 251-500)...")
    data2 = read_csv_file_and_group('docking_results_90st_temp0.csv', {'min_trial': 46, 'max_trial': 90})

    print(f"\nFiltering Summary:")
    print(f"Dataset 1 final count: {len(data1)} mutation groups")
    print(f"Dataset 2 final count: {len(data2)} mutation groups")

    if len(data1) == 0 or len(data2) == 0:
        print("One or both datasets have no data!")
        return

    # Match the data based on mutation description
    print("Matching mutation groups between datasets...")
    matched_data = match_data(data1, data2)

    print(f"Found {len(matched_data)} matching mutation groups")

    if len(matched_data) == 0:
        print("No matching mutation groups found!")
        print("Sample keys from dataset 1:", list(data1.keys())[:5])
        print("Sample keys from dataset 2:", list(data2.keys())[:5])
        return

    print("Creating interface_delta_X scatter plot and calculating statistics...")
    interface_stats_dict, interface_plot_filename = create_scatter_plot_and_stats(matched_data, 'deep_500st_split_analysis', 'interface_delta_X')

    print("Saving interface_delta_X detailed statistics...")
    interface_stats_filename = save_detailed_stats(interface_stats_dict, matched_data, 'deep_500st_split_analysis', 'interface_delta_X')

    print("Creating total_score scatter plot and calculating statistics...")
    total_stats_dict, total_plot_filename = create_scatter_plot_and_stats(matched_data, 'deep_500st_split_analysis', 'total_score')

    print("Saving total_score detailed statistics...")
    total_stats_filename = save_detailed_stats(total_stats_dict, matched_data, 'deep_500st_split_analysis', 'total_score')
    '''
    # Print summary to console
    print(f"\n=== ddG BAR GRAPH ANALYSIS COMPLETE ===")
    print(f"Original per-mutation bar graphs:")
    print(f"  Total score ddG bar graph saved as: {total_plot_file}")
    print(f"  Total score ddG data saved as: {total_csv_file}")
    print(f"  Interface delta ddG bar graph saved as: {interface_plot_file}")
    print(f"  Interface delta ddG data saved as: {interface_csv_file}")

    print(f"\nPer-position summary bar graphs:")
    print(f"  Total score mean ddG bar graph saved as: {total_mean_plot}")
    print(f"  Total score mean ddG data saved as: {total_mean_csv}")
    print(f"  Total score median ddG bar graph saved as: {total_median_plot}")
    print(f"  Total score median ddG data saved as: {total_median_csv}")
    print(f"  Interface delta mean ddG bar graph saved as: {interface_mean_plot}")
    print(f"  Interface delta mean ddG data saved as: {interface_mean_csv}")
    print(f"  Interface delta median ddG bar graph saved as: {interface_median_plot}")
    print(f"  Interface delta median ddG data saved as: {interface_median_csv}")

    print(f"\nBox plots for top 20 positions:")
    print(f"  Total score box plot saved as: {total_boxplot}")
    print(f"  Total score box plot data saved as: {total_boxplot_csv}")
    print(f"  Interface delta box plot saved as: {interface_boxplot}")
    print(f"  Interface delta box plot data saved as: {interface_boxplot_csv}")

    print(f"\n=== WILD TYPE ddG ANALYSIS COMPLETE ===")
    print(f"Wild type-based per-mutation bar graphs:")
    print(f"  Total score ddG vs WT bar graph saved as: {wt_total_plot_file}")
    print(f"  Total score ddG vs WT data saved as: {wt_total_csv_file}")
    print(f"  Interface delta ddG vs WT bar graph saved as: {wt_interface_plot_file}")
    print(f"  Interface delta ddG vs WT data saved as: {wt_interface_csv_file}")

    print(f"\nWild type-based per-position summary bar graphs:")
    print(f"  Total score mean ddG vs WT bar graph saved as: {wt_total_mean_plot}")
    print(f"  Total score mean ddG vs WT data saved as: {wt_total_mean_csv}")
    print(f"  Total score median ddG vs WT bar graph saved as: {wt_total_median_plot}")
    print(f"  Total score median ddG vs WT data saved as: {wt_total_median_csv}")
    print(f"  Interface delta mean ddG vs WT bar graph saved as: {wt_interface_mean_plot}")
    print(f"  Interface delta mean ddG vs WT data saved as: {wt_interface_mean_csv}")
    print(f"  Interface delta median ddG vs WT bar graph saved as: {wt_interface_median_plot}")
    print(f"  Interface delta median ddG vs WT data saved as: {wt_interface_median_csv}")

    print(f"\nBox plots for top 20 lowest positions vs wild type:")
    print(f"  Total score box plot saved as: {wt_total_boxplot}")
    print(f"  Total score box plot data saved as: {wt_total_boxplot_csv}")
    print(f"  Interface delta box plot saved as: {wt_interface_boxplot}")
    print(f"  Interface delta box plot data saved as: {wt_interface_boxplot_csv}")
    '''
    print(f"\n=== CORRELATION ANALYSIS COMPLETE ===")
    print(f"Interface_delta_X plot saved as: {interface_plot_filename}")
    print(f"Interface_delta_X statistics saved as: {interface_stats_filename}")
    print(f"Total_score plot saved as: {total_plot_filename}")
    print(f"Total_score statistics saved as: {total_stats_filename}")

    print(f"\nQuick Summary - Interface_delta_X:")
    print(f"Pearson correlation: {interface_stats_dict['pearson_r']:.4f} (p = {interface_stats_dict['pearson_p']:.2e})")
    print(f"Spearman correlation: {interface_stats_dict['spearman_r']:.4f} (p = {interface_stats_dict['spearman_p']:.2e})")
    print(f"RMSE: {interface_stats_dict['rmse']:.4f}")

    print(f"\nQuick Summary - Total_score:")
    print(f"Pearson correlation: {total_stats_dict['pearson_r']:.4f} (p = {total_stats_dict['pearson_p']:.2e})")
    print(f"Spearman correlation: {total_stats_dict['spearman_r']:.4f} (p = {total_stats_dict['spearman_p']:.2e})")
    print(f"RMSE: {total_stats_dict['rmse']:.4f}")

    print(f"\nNumber of data points: {len(matched_data)}")
    '''
    # Print top 10 mutations with largest deviations between trials
    print("\n=== TOP 10 MUTATIONS WITH LARGEST DEVIATION BETWEEN TRIALS ===")

    # Calculate deviations for all matched mutations
    deviations = []
    for description in matched_data:
        total_score_dev = abs(matched_data[description]['file1']['total_score'] -
                             matched_data[description]['file2']['total_score'])
        interface_dev = abs(matched_data[description]['file1']['interface_delta_X'] -
                           matched_data[description]['file2']['interface_delta_X'])

        deviations.append({
            'mutation': description,
            'total_score_trial1': matched_data[description]['file1']['total_score'],
            'total_score_trial2': matched_data[description]['file2']['total_score'],
            'total_score_deviation': total_score_dev,
            'interface_delta_X_trial1': matched_data[description]['file1']['interface_delta_X'],
            'interface_delta_X_trial2': matched_data[description]['file2']['interface_delta_X'],
            'interface_delta_X_deviation': interface_dev
        })

    # Sort by total_score deviation and get top 10
    top_total_dev = sorted(deviations, key=lambda x: x['total_score_deviation'], reverse=True)[:80]

    print("\nTop 10 mutations with largest total_score deviation:")
    print("-" * 100)
    print(f"{'Rank':<4} {'Mutation':<15} {'Trial 1':<15} {'Trial 2':<15} {'Deviation':<15}")
    print("-" * 100)
    for i, mut in enumerate(top_total_dev, 1):
        print(f"{i:<4} {mut['mutation']:<15} {mut['total_score_trial1']:<15.4f} {mut['total_score_trial2']:<15.4f} {mut['total_score_deviation']:<15.4f}")

    # Sort by interface_delta_X deviation and get top 10
    top_interface_dev = sorted(deviations, key=lambda x: x['interface_delta_X_deviation'], reverse=True)[:80]

    print("\nTop 10 mutations with largest interface_delta_X deviation:")
    print("-" * 100)
    print(f"{'Rank':<4} {'Mutation':<15} {'Trial 1':<15} {'Trial 2':<15} {'Deviation':<15}")
    print("-" * 100)
    for i, mut in enumerate(top_interface_dev, 1):
        print(f"{i:<4} {mut['mutation']:<15} {mut['interface_delta_X_trial1']:<15.4f} {mut['interface_delta_X_trial2']:<15.4f} {mut['interface_delta_X_deviation']:<15.4f}")

    # Save to CSV files
    total_dev_df = pd.DataFrame(top_total_dev)
    total_dev_csv = 'deep_500st_split_analysis_top_total_score_deviations.csv'
    total_dev_df.to_csv(total_dev_csv, index=False)
    print(f"\nTop total_score deviations saved to: {total_dev_csv}")

    interface_dev_df = pd.DataFrame(top_interface_dev)
    interface_dev_csv = 'deep_500st_split_analysis_top_interface_delta_X_deviations.csv'
    interface_dev_df.to_csv(interface_dev_csv, index=False)
    print(f"Top interface_delta_X deviations saved to: {interface_dev_csv}")

    # Create scatter plots for top deviation mutations
    print("\n=== CREATING SCATTER PLOTS FOR TOP DEVIATION MUTATIONS ===")

    # Filter matched_data to only include top deviation mutations
    top_total_dev_keys = set([mut['mutation'] for mut in top_total_dev])
    top_interface_dev_keys = set([mut['mutation'] for mut in top_interface_dev])

    matched_data_top_total = {k: v for k, v in matched_data.items() if k in top_total_dev_keys}
    matched_data_top_interface = {k: v for k, v in matched_data.items() if k in top_interface_dev_keys}

    print(f"Creating scatter plot for top {len(matched_data_top_total)} total_score deviation mutations...")
    total_dev_stats_dict, total_dev_plot_filename = create_scatter_plot_and_stats(
        matched_data_top_total, 'deep_500st_split_analysis_top_total_dev', 'total_score')
    total_dev_stats_filename = save_detailed_stats(
        total_dev_stats_dict, matched_data_top_total, 'deep_500st_split_analysis_top_total_dev', 'total_score')

    print(f"Creating scatter plot for top {len(matched_data_top_interface)} interface_delta_X deviation mutations...")
    interface_dev_stats_dict, interface_dev_plot_filename = create_scatter_plot_and_stats(
        matched_data_top_interface, 'deep_500st_split_analysis_top_interface_dev', 'interface_delta_X')
    interface_dev_stats_filename = save_detailed_stats(
        interface_dev_stats_dict, matched_data_top_interface, 'deep_500st_split_analysis_top_interface_dev', 'interface_delta_X')

    print(f"\nTop deviation scatter plots saved:")
    print(f"  Total score: {total_dev_plot_filename}")
    print(f"  Total score stats: {total_dev_stats_filename}")
    print(f"  Interface delta X: {interface_dev_plot_filename}")
    print(f"  Interface delta X stats: {interface_dev_stats_filename}")

    # Print mutations better than wild type
    print_mutations_better_than_wt(data_deep, 'deep_wt_analysis')

    # Print top 10 worst mutations from the deep dataset
    print_top_worst_mutations(data_deep, 10)
    '''
if __name__ == "__main__":
    main()
