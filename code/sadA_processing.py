


def preprocess_single_variants():
    with open('/Users/brycejohnson/Desktop/proteins/metl-sim/variant_lists/SadA_AlphaFold_p_all_NS-1.txt', 'r') as f:
        out= f.read()

    out= out.replace('SadA_AlphaFold_p.pdb','SadX_NSLeu_Corrected_3701_best_structure_0044_correct_seq_100.pdb')

    with open('/Users/brycejohnson/Desktop/proteins/metl-sim/variant_lists/SadX_NSLeu_Corrected_3701_best_structure_0044_correct_seq_100_NS-1.txt','w') as f:
        f.write(out)

if __name__ == '__main__':
    preprocess_single_variants()