import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from utils import get_seq_from_pdb


def heatmap(input_pdb,input_fn_df,output_fn,field='interface_delta_X'):


    seq = get_seq_from_pdb(input_pdb)
    seq = list(seq)


    df = pd.read_csv(input_fn_df)  # ADD DATAFRAME HERE###
    # df= df[~df['interface_delta_X'].between(-15.5, -14.3, inclusive='both')]


    if field=='interface_delta_X':
        WT_energy = -15
        params= {'vmin':-5, 'vmax':5,'center':0}
    else:
        WT_energy=0
        params={'center':0}



    correct_aas=list('QENHDRKTSAGMLVIWYFPC')
    out = np.zeros((20, len(seq) + 1))
    for row in df.iterrows():
        variant =row[1]['variant']
        wt_aa, pos, mut_aa  =variant[0] , int(variant[1:-1]),variant[-1]

        assert wt_aa==SadX[pos-1]
        out[correct_aas.index(mut_aa)][pos-1] = row[1][field]-WT_energy
    plt.figure(figsize=(50, 5), facecolor="white")

    sns.set(rc={'figure.figsize': (20, 15)})

    sns.heatmap(out,cmap='bwr_r', yticklabels=[x for x in correct_aas], xticklabels=[
        (str(x) + seq[x - 1]) for x in range(len(seq) + 1)[1:]],**params)


    for i in range(0, len(seq)):
        # print(i)
        plt.scatter(i+0.5, correct_aas.index(seq[i])+0.5, s=30, color='black')


    plt.savefig(output_fn)



    # display(df)


    #df=df[~df['interface_delta_X'].between(-15.5, -14.3, inclusive='both')]

aa_map = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
    "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
    "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
    "W": "TRP", "Y": "TYR"
}

reversed_aa_map = {value: key for key, value in aa_map.items()}

SadX='MQHTYPAQLMRFGTAARAEHMTIAAAIHALDADEADAIVMDIVPDGERDAWWDDEGFSSSPFTKNAHHAGIVATSVTLGQLQREQGDKLVSKAAEYFGIACRVNDGLRTTRFVRLFSDALDAKPLTIGHDYEVEFLLATRRVYEPFEAPFNFAPHCGDVSYGRDTVNWPLKRSFPRQLGGFLTIQGADNDAGMVMWDNRPESRAALDEMHAEYRETGAIAALERAAKIMLKPQPGQLTLFQSKNLHAIERCTSTRRTMGLFLIHTEDGWRMFD'
global_list =[]
def return_Variant(x):
    x= x.split('_')[1].split('.')[0]

    mut_aa=reversed_aa_map[x[-3:]]

    idx=int(x[:-3])-1
    wt_aa= SadX[idx]
    if mut_aa==wt_aa:
        print(f'{wt_aa}{idx+1}{mut_aa}')
        # global_list.append(f'{wt_aa}{idx}{mut_aa}')
        # return '_wt'
    # assert mut_aa!=wt_aa
    return f'{wt_aa}{idx+1}{mut_aa}'

def process_score_file_andres():
    df=pd.read_csv('templates/sadA_docking_andres_template/score_2024_8_14_ver_2021.36+release.57ac713 (3).sc',sep=r'\s+',skiprows=1)

    df['variant'] = df['description'].apply(lambda x: return_Variant(x))
    df=df.drop(columns=['SCORE:','description'])
    df.set_index('variant',inplace=True)
    df.to_csv('output/energize_outputs/energize_andres_2024-09-24/energize.csv')


def diff_scoring_in_andres_docking_protocol():
    df1=pd.read_csv('output/energize_outputs/energize_andres_2024-09-24/energize.csv').set_index('variant')
    df2=pd.read_csv('output/htcondor_runs/condor_energize_2024-09-23_00-16-21_sadA_singles_andres_docking/processed_run/energies_df.csv').set_index('variant')

    # drop duplicates
    df1=df1.loc[~df1.index.duplicated(keep='first')]
    df2 = df2.loc[~df2.index.duplicated(keep='first')]


    # remove the wildtype variants that andres generated and i did not


    df1['interface_delta_X_metl_sim'] = df2['interface_delta_X']

    df1 = df1[~np.isnan(df1['interface_delta_X_metl_sim'])]


    df1['interface_delta_X-X_metl_sim']  = df1['interface_delta_X'] - df1['interface_delta_X_metl_sim']


    df1.reset_index(inplace=True)
    df = df1[['variant','interface_delta_X-X_metl_sim']]

    df.to_csv('output/energize_outputs/energize_andres_2024-09-25/energize_diff.csv')

    heatmap('pdb_files/prepared_pdb_files/SadA_NSLeu_Corrected_3701_best_structure_0044_correct_seq_2024_3_6_unrelaxed_ver_2021.36+release.57ac713_p.pdb',
        'output/energize_outputs/energize_andres_2024-09-25/energize_diff.csv',
        'output/energize_outputs/energize_andres_2024-09-25/heatmap_interface_delta_X-X_metl_sim.png',
        'interface_delta_X-X_metl_sim')



if __name__ == '__main__':
    # diff_scoring_in_andres_docking_protocol()



    # process_score_file_andres()

    # heatmap('pdb_files/prepared_pdb_files/SadX_NSLeu_Corrected_3701_best_structure_0044_correct_seq_100.pdb',
    #         'output/htcondor_runs/condor_energize_2024-05-29_16-48-56_sadA_singles_docking/processed_run/energies_df.csv',
    #         'output/htcondor_runs/condor_energize_2024-05-29_16-48-56_sadA_singles_docking/processed_run/heatmap_filtered.png')

    heatmap('pdb_files/prepared_pdb_files/SadA_NSLeu_Corrected_3701_best_structure_0044_correct_seq_2024_3_6_unrelaxed_ver_2021.36+release.57ac713_p.pdb',
            'output/htcondor_runs/condor_energize_2024-09-23_00-16-21_sadA_singles_andres_docking/processed_run/energies_df.csv',
            'output/htcondor_runs/condor_energize_2024-09-23_00-16-21_sadA_singles_andres_docking/heatmap_interface_delta_X.png',
            'interface_delta_X')

    # heatmap(
    #     'pdb_files/prepared_pdb_files/SadA_NSLeu_Corrected_3701_best_structure_0044_correct_seq_2024_3_6_unrelaxed_ver_2021.36+release.57ac713_p.pdb',
    #     'output/energize_outputs/energize_andres_2024-09-24/energize.csv',
    #     'output/energize_outputs/energize_andres_2024-09-24/heatmap_interface_delta_X.png',
    #     'interface_delta_X')

    # df=pd.read_csv('output/energize_outputs/energize_local_local_2024-09-23_04-49-41_CczCue7hcLPM/energies.csv')
    # print(df[['variant','interface_delta_X']])