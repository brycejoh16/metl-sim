import sys
sys.path.insert(0, "/Users/brycejohnson/Desktop/proteins/metl-sim/code")
from protocol_analysis import parity_plot
import pandas as pd

if __name__ == '__main__':
    df_500k = pd.read_csv('/Users/brycejohnson/Desktop/proteins/metl-sim/SadX/large_dock_mutants_v3_ligand+mutate_area_no_low_res/processed_run/energies_df.csv')
    df_500k_singles = df_500k.loc[df_500k.groupby(df_500k['variant'])["interface_delta_X"].idxmin()]
    df_low_res = pd.read_csv('/Users/brycejohnson/Desktop/proteins/metl-sim/SadX/v1_ligand+mutate_area_no_low_res/processed_run/energies_df.csv')
    df_low_res_singles = df_low_res.loc[df_low_res.groupby(df_low_res['variant'])["interface_delta_X"].idxmin()]

    df_low_res_singles.set_index('variant',inplace=True)
    df_500k_singles.set_index('variant', inplace=True)

    df_low_res_singles['interface_delta_X_500k']=  df_500k_singles['interface_delta_X']
    df_low_res_singles['run_time_500k'] = df_500k_singles['run_time']




    parity_plot(x=df_low_res_singles['interface_delta_X'],
                y=df_low_res_singles['interface_delta_X_500k'],
                x_label='WT (Normal Prepare)',
                y_label='WT (500,000 nstruct Prepare)',
               save_path='low_res_singles_vs_500k_singles.png')

    parity_plot(x=df_low_res_singles['run_time'],
                y=df_low_res_singles['run_time_500k'],
                x_label='Runtime (seconds) for WT (Normal Prepare)',
                y_label='Runtime (seconds) for WT (500,000 nstruct Prepare)',
                save_path='low_res_singles_vs_500k_singles_runtime.png')
