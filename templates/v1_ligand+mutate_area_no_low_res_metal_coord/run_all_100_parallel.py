from subprocess import Popen, STDOUT, PIPE
import os
import sys
import pandas as pd
import glob
from multiprocessing import Pool
import time
import math

threetoone = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
            'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
onetothree = dict([(threetoone[x],x) for x in threetoone.keys()])

# Read the XML templates at the start
mutate_template = open('SadA_mutate.xml').read().strip()
dock_template = open('dock_template.xml').read().strip()

# Define custom folder names
base_mutants_folder = "Mutants_2st_temp1"
base_output_folder = "Output_2st_temp1"

def process_mutation_group(args):
    """Process a group of mutations in parallel with isolated files and folders"""
    group_id, mutations = args
    
    # Create isolated folders for this process group
    mutants_folder = f"{base_mutants_folder}_group_{group_id}"
    output_folder = f"{base_output_folder}_group_{group_id}"
    
    # Ensure directories exist
    os.makedirs(mutants_folder, exist_ok=True)
    os.makedirs(output_folder, exist_ok=True)
    
    group_results = []
    
    print(f"Group {group_id}: Starting {len(mutations)} mutations")
    
    for pos, aa in mutations:
        mutation_id = f"{pos}_{aa}"
        
        try:
            # --- 1. Mutation Step ---
            mut_file = mutate_template.replace(str("POSITION"), str(pos)).replace(str("AMINOACID"), aa)
            xml_filename = f'tmp_mut_file_group_{group_id}_{mutation_id}.xml'

            with open(xml_filename, 'w') as f:
                f.write(mut_file)

            print(f"Group {group_id}: Processing mutation {pos}{aa}")

            options_mutate_template = open('options_mutate_1.txt').read()
            temp_mutate_options_content = options_mutate_template.replace('tmp_mut_file_1.xml', xml_filename)
            temp_mutate_options_file = f'options_mutate_group_{group_id}_{mutation_id}.txt'

            with open(temp_mutate_options_file, 'w') as f:
                f.write(temp_mutate_options_content)

            cmd = f'/home/romeroroot/rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease @{temp_mutate_options_file} -out:path:pdb {mutants_folder}/'
            process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                print(f"Group {group_id}: Mutation failed for {pos}{aa}: {stderr.decode()}")
                continue
            
            mutated_pdb = f'{mutants_folder}/SadA_NSLeu_Corrected_3701_best_structure_0044_0001.pdb'
            temp_pdb = f'{mutants_folder}/temp.pdb'
            
            if os.path.isfile(mutated_pdb):
                os.rename(mutated_pdb, temp_pdb)
            else:
                print(f"Group {group_id}: Mutated structure not found for {pos}{aa}")
                continue
            
            # --- 2. Docking Step (with dynamic XML) ---

            # Create an isolated dock XML for this mutation with the correct residue number
            dock_xml_content = dock_template.replace('MUTATION_POSITION', str(pos))
            dock_xml_filename = f'tmp_dock_file_group_{group_id}_{mutation_id}.xml'
            with open(dock_xml_filename, 'w') as f:
                f.write(dock_xml_content)

            # Create a temporary options file for docking
            options_template = open('options_dock_temp_1.txt').read()
            temp_options_content = options_template.replace('Mutants/temp.pdb', f'{mutants_folder}/temp.pdb')
            
            # MODIFIED: Target the dynamically created dock XML
            temp_options_content = temp_options_content.replace('dock.xml', dock_xml_filename)
            
            temp_options_file = f'options_dock_temp_group_{group_id}_{mutation_id}.txt'
            
            with open(temp_options_file, 'w') as f:
                f.write(temp_options_content)
            
            score_file_path = f'{output_folder}/score.sc'
            if os.path.exists(score_file_path):
                os.remove(score_file_path)
            
            cmd = f'/home/romeroroot/rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease @{temp_options_file} -packing:ex1 -packing:ex2 -out:suffix -nstruct 2 -out:path:pdb {output_folder}/ -out:path:score {output_folder}/'
            process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            stdout, stderr = process.communicate()
            
            if process.returncode != 0:
                print(f"Group {group_id}: Docking failed for {pos}{aa}: {stderr.decode()}")
                continue
            
            # --- 3. Process Results ---
            if os.path.isfile(score_file_path):
                with open(score_file_path, 'r') as f:
                    score_content = f.read().strip()
                    score_lines = score_content.split('\n')
                
                header_line_idx = -1
                header = []
                for idx, line in enumerate(score_lines):
                    if line.startswith('SCORE:') and 'total_score' in line:
                        header = line.split()[1:]
                        header_line_idx = idx
                        break
                
                if header_line_idx != -1:
                    for i, line in enumerate(score_lines[header_line_idx + 1:], 1):
                        if line.strip() and line.startswith('SCORE:'):
                            parts = line.split()[1:]
                            if len(parts) >= len(header):
                                result = {'Mutation_point': pos, 'Mutate_into': aa, 'trial': i}
                                for j, col_name in enumerate(header):
                                    if j < len(parts):
                                        result[col_name] = parts[j]
                                group_results.append(result)
            
            # Clean up temporary files
            if os.path.exists(xml_filename):
                os.remove(xml_filename)
            if os.path.exists(temp_mutate_options_file):
                os.remove(temp_mutate_options_file)
            if os.path.exists(dock_xml_filename):
                os.remove(dock_xml_filename)
            if os.path.exists(temp_options_file):
                os.remove(temp_options_file)
                
            # Rename PDB files
            for trial_num in range(1, 3):
                old_pdb_name = f'{output_folder}/temp_{trial_num:04d}.pdb'
                new_pdb_name = f'{output_folder}/SadA_NSLeu_{pos}{aa}_trial_{trial_num:03d}.pdb'
                if os.path.isfile(old_pdb_name):
                    os.rename(old_pdb_name, new_pdb_name)

            mutation_structures = [r for r in group_results if r['Mutation_point']==pos and r['Mutate_into']==aa]
            print(f"Group {group_id}: Completed {pos}{aa} with {len(mutation_structures)} structures.")
            
        except Exception as e:
            print(f"Group {group_id}: Error processing mutation {pos}{aa}: {str(e)}")
            continue
    
    print(f"Group {group_id}: Completed {len(mutations)} mutations with {len(group_results)} total structures")
    return group_id, group_results, mutants_folder, output_folder

def main():
    # Example: Process mutations at positions 1 through 273
    specific_positions = range(1, 274)
    #specific_positions = range(1,2) # For testing
    
    
    metal_binding_sites = [155, 157, 246]


    mutations = []
    for pos in specific_positions:
        if pos in metal_binding_sites:
            print(f"Skipping metal binding site: {pos}")
            continue
        
        for aa in threetoone.keys():
        #for aa in ['ALA']:  # For testing
            mutations.append((pos,aa))
    
    total_mutations = len(mutations)
    num_cores = 1 # Adjust based on your machine
    mutations_per_core = math.ceil(total_mutations / num_cores)
    
    print(f"Total mutations to process: {total_mutations}")
    print(f"Using {num_cores} cores")
    print(f"~{mutations_per_core} mutations per core")
    
    # Split mutations into groups
    mutation_groups = [(i, mutations[i*mutations_per_core:min((i+1)*mutations_per_core, total_mutations)]) for i in range(num_cores)]
    
    os.makedirs(base_mutants_folder, exist_ok=True)
    os.makedirs(base_output_folder, exist_ok=True)
    
    start_time = time.time()
    
    with Pool(processes=num_cores) as pool:
        results = pool.map(process_mutation_group, mutation_groups)
    
    all_results = []
    for group_id, group_results, mutants_folder, output_folder in results:
        all_results.extend(group_results)
        pdb_files = glob.glob(f'{output_folder}/SadA_NSLeu_*.pdb')
        for pdb_file in pdb_files:
            filename = os.path.basename(pdb_file)
            dest_path = f'{base_output_folder}/{filename}'
            if not os.path.exists(dest_path):
                os.rename(pdb_file, dest_path)
        
        try:
            if os.path.exists(output_folder):
                os.system(f'rm -rf {output_folder}')
            if os.path.exists(mutants_folder):
                os.system(f'rm -rf {mutants_folder}')
        except Exception as e:
            print(f"Could not clean up folders for group {group_id}: {e}")
    
    if all_results:
        df = pd.DataFrame(all_results)
        df.to_csv(f'{base_output_folder}/docking_results_parallel.csv', index=False)
        print(f"Saved {len(all_results)} results to {base_output_folder}/docking_results_parallel.csv")
    
    end_time = time.time()
    total_time = end_time - start_time
    
    successful_mutations = len(set((r['Mutation_point'], r['Mutate_into']) for r in all_results))
    
    print(f"\n=== PARALLEL PROCESSING COMPLETE ===")
    print(f"Total processing time: {total_time:.2f} seconds ({total_time/60:.1f} minutes)")
    if successful_mutations > 0:
        print(f"Successful mutations: {successful_mutations}/{total_mutations}")
        print(f"Total structures generated: {len(all_results)}")
        print(f"Average time per mutation: {total_time/successful_mutations:.2f} seconds")

if __name__ == "__main__":
    main()