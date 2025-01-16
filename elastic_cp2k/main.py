import numpy as np
import shutil
from ase import io
import subprocess
import os
from pathlib import Path
import pandas as pd
import yaml
import argparse

from .utils import (
    check_existing_results,
    create_folders,
    copy_files_to_folders,
    run_cp2k_optimization,
    extract_last_frame,
    read_last_frame_cell_vectors,
    save_xyz_file,
    extract_last_stress_tensor,
    modify_cp2k_input,
    remap_xyz,
    symmetrize_and_save,
    process_epsilon
)

def load_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config

def main():
    parser = argparse.ArgumentParser(description="elastic_cp2k")
    parser.add_argument(
        '-c', '--config',
        type=str,
        default='config.yaml',
        help='Path to the YAML configuration file.'
    )
    args = parser.parse_args()
    
    if not os.path.exists(args.config):
        print(f"Error: Configuration file {args.config} not found.")
        return
    
    # Load configuration
    config = load_config(args.config)
    
    # Assign configuration variables
    num_cores = config['num_cores']
    cp2k_exec_path = config['cp2k_exec_path']
    cellopt_inp = config['cellopt_inp']
    cellopt_out = config['cellopt_out']
    cellopt_jobname = config['cellopt_jobname']
    geoopt_inp = config['geoopt_inp']
    geoopt_out = config['geoopt_out']
    geoopt_jobname = config['geoopt_jobname']
    potential_file = config['potential_file']
    basis_set_file = config['basis_set_file']
    dftd3_dat = config['dftd3_dat']
    up = config['up']
    
    # Define base_folders within the code
    base_folders = [
        "./epsilonxx", 
        "./epsilonyy", 
        "./epsilonzz", 
        "./epsilonyz", 
        "./epsilonxz", 
        "./epsilonxy"
    ]
    
    # Define epsilon values
    epsilon_values = {
        "./epsilonxx": np.array([[up, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
        "./epsilonyy": np.array([[0.0, 0.0, 0.0], [0.0, up, 0.0], [0.0, 0.0, 0.0]]),
        "./epsilonzz": np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, up]]),
        "./epsilonyz": np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, up, 0.0]]),
        "./epsilonxz": np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [up, 0.0, 0.0]]),
        "./epsilonxy": np.array([[0.0, 0.0, 0.0], [up, 0.0, 0.0], [0.0, 0.0, 0.0]]),
    }
    
    base_path = Path.cwd()
    
    print(f"{'*' * 10} CELL OPTIMIZATION STARTED {'*' * 10}")
    if os.path.exists(cellopt_inp):
        run_cp2k_optimization(cellopt_inp, cellopt_out, cellopt_jobname, num_cores, cp2k_exec_path)
    else:
        print(f"Error: The input file {cellopt_inp} does not exist.")
    print(f"{'*' * 10} CELL OPTIMIZATION ENDED {'*' * 10}\n")
    
    print("Generating the optimized structure...")
    cellopt_file_path = f'./{cellopt_jobname}-pos-1.xyz'
    last_frame_atoms = extract_last_frame(cellopt_file_path)
    
    file_path = f'./{cellopt_jobname}-1.cell'
    last_frame_cell_vectors = read_last_frame_cell_vectors(file_path)
    
    cellopt_xyz_output = './coord_cellopt.xyz'
    save_xyz_file(last_frame_atoms, last_frame_cell_vectors, cellopt_xyz_output)
    print("Optimized structure saved\n")
    
    print("Creating folders for cell stretching...")
    create_folders(base_folders)
    files_to_copy = [potential_file, basis_set_file, dftd3_dat, cellopt_xyz_output]
    copy_files_to_folders(files_to_copy, base_folders)
    print("Done!\n")
    
    stress_tensor_cellopt = extract_last_stress_tensor(cellopt_out)
    cell_vectors_cellopt = last_frame_cell_vectors
    
    print("Applying strain to the cell....")
    input_file_path = f"./{geoopt_inp}"
    for folder, epsilon in epsilon_values.items():
        process_epsilon(folder, epsilon, cell_vectors_cellopt, input_file_path, num_cores, cp2k_exec_path)
    print("Done!\n")
    
    print(f"{'*' * 10} GEO OPTIMIZATION STARTED {'*' * 10}")
    for folder in base_folders:
        pos_path = os.path.join(base_path, folder, "pos")
        neg_path = os.path.join(base_path, folder, "neg")
        
        os.chdir(pos_path)
        print(f"Switch to the {Path.cwd()} directory for geometry optimization.")
        run_cp2k_optimization(geoopt_inp, geoopt_out, geoopt_jobname, num_cores, cp2k_exec_path)
        
        os.chdir(neg_path)
        print(f"Switch to the {Path.cwd()} directory for geometry optimization.")
        run_cp2k_optimization(geoopt_inp, geoopt_out, geoopt_jobname, num_cores, cp2k_exec_path)
    os.chdir(base_path)
    print(f"{'*' * 10} GEO OPTIMIZATION ENDED {'*' * 10}\n")
    
    print(f"{'*' * 10} CALCULATING ELASTIC CONSTANTS {'*' * 10}")
    Cij = np.zeros((6, 6))
    for index, folder in enumerate(base_folders):
        os.chdir(base_path)
        stress_tensor_cellopt = extract_last_stress_tensor(cellopt_out)
        
        os.chdir(os.path.join(base_path, folder, "pos"))
        stress_tensor = extract_last_stress_tensor(geoopt_out)
        Cn1_pos = - (stress_tensor[0,0] - stress_tensor_cellopt[0,0]) / (up)
        Cn2_pos = - (stress_tensor[1,1] - stress_tensor_cellopt[1,1]) / (up)
        Cn3_pos = - (stress_tensor[2,2] - stress_tensor_cellopt[2,2]) / (up)
        Cn4_pos = - (stress_tensor[1,2] - stress_tensor_cellopt[1,2]) / (up) # yz
        Cn5_pos = - (stress_tensor[0,2] - stress_tensor_cellopt[0,2]) / (up) # xz
        Cn6_pos = - (stress_tensor[0,1] - stress_tensor_cellopt[0,1]) / (up) # xy
        
        os.chdir(os.path.join(base_path, folder, "neg"))
        stress_tensor = extract_last_stress_tensor(geoopt_out)
        Cn1_neg = (stress_tensor[0,0] - stress_tensor_cellopt[0,0]) / (up)
        Cn2_neg = (stress_tensor[1,1] - stress_tensor_cellopt[1,1]) / (up)
        Cn3_neg = (stress_tensor[2,2] - stress_tensor_cellopt[2,2]) / (up)
        Cn4_neg = (stress_tensor[1,2] - stress_tensor_cellopt[1,2]) / (up) # yz
        Cn5_neg = (stress_tensor[0,2] - stress_tensor_cellopt[0,2]) / (up) # xz
        Cn6_neg = (stress_tensor[0,1] - stress_tensor_cellopt[0,1]) / (up) # xy
    
        Cn1 =  np.round((Cn1_pos + Cn1_neg) / 2, 3)
        Cn2 =  np.round((Cn2_pos + Cn2_neg) / 2, 3)
        Cn3 =  np.round((Cn3_pos + Cn3_neg) / 2, 3)
        Cn4 =  np.round((Cn4_pos + Cn4_neg) / 2, 3)
        Cn5 =  np.round((Cn5_pos + Cn5_neg) / 2, 3)
        Cn6 =  np.round((Cn6_pos + Cn6_neg) / 2, 3)
        Cij[:, index] = [Cn1, Cn2, Cn3, Cn4, Cn5, Cn6]

        os.chdir(base_path)
    
    Cij = symmetrize_and_save(Cij)
    df = pd.DataFrame(Cij)
    print("\nCij:")
    print(df.to_string(index=False, header=False))
    print(f"{'*' * 10} DONE! {'*' * 10}")

if __name__ == "__main__":
    main()