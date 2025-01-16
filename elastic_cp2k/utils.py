import os
import shutil
import subprocess
import numpy as np
from ase import io

def check_existing_results(cp2k_out, jobname):
    if os.path.exists(cp2k_out) and os.path.exists(f"{jobname}-1.cell") and os.path.exists(f"{jobname}-pos-1.xyz"):
        return True
    else:
        try:
            if os.path.exists(cp2k_out):
                os.remove(cp2k_out)
            if os.path.exists(f"{jobname}-1.cell"):
                os.remove(f"{jobname}-1.cell")
            if os.path.exists(f"{jobname}-pos-1.xyz"):
                os.remove(f"{jobname}-pos-1.xyz")
        except FileNotFoundError:
            pass
        return False

def create_folders(base_folders):
    for base_folder in base_folders:
        pos_folder = os.path.join(base_folder, "pos")
        neg_folder = os.path.join(base_folder, "neg")
        
        if os.path.exists(base_folder):
            shutil.rmtree(base_folder)
        
        os.makedirs(pos_folder, exist_ok=True)
        os.makedirs(neg_folder, exist_ok=True)

def copy_files_to_folders(files, base_folders):
    for folder in base_folders:
        for subfolder in ["pos", "neg"]:
            target_folder = os.path.join(folder, subfolder)
            
            if os.path.exists(target_folder):
                for file in files:
                    if os.path.exists(file):
                        shutil.copy(file, target_folder)
                    else:
                        print(f"Warning: {file} not found. Skipping file.")
            else:
                print(f"Warning: Target folder {target_folder} not found. Skipping.")

def run_cp2k_optimization(cp2k_inp, cp2k_out, jobname, num_cores, cp2k_exec_path):
    if check_existing_results(cp2k_out, jobname):
        print("This folder already contains the optimization results.")
        return
        
    command = [
        "mpirun", 
        "-n", str(num_cores), 
        cp2k_exec_path, 
        "-o", cp2k_out, 
        cp2k_inp
    ]
    
    try:
        print(f"Running CP2K with {num_cores} cores...")
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        print(f"CP2K run completed successfully. Output is saved to {cp2k_out}.")
        print("Standard Output:", result.stdout.decode())
        print("Standard Error:", result.stderr.decode())
    
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during CP2K execution: {e}")
        if e.stdout:
            print("Standard Output:", e.stdout.decode())
        if e.stderr:
            print("Standard Error:", e.stderr.decode())

def extract_last_frame(file_path):
    with open(file_path, 'r') as file:
        file_content = file.readlines()

    frame_start_index = None
    for i in range(len(file_content) - 1, -1, -1):
        if file_content[i].startswith(" i ="):
            frame_start_index = i
            break

    if frame_start_index is None:
        print(f"Warning: No frame found in {file_path}.")
        return []

    atoms_in_last_frame = file_content[frame_start_index + 1:]

    frame_end_index = None
    for i in range(len(atoms_in_last_frame)):
        if atoms_in_last_frame[i].startswith(" i ="):
            frame_end_index = i
            break
            
    if frame_end_index is not None:
        atoms_in_last_frame = atoms_in_last_frame[:frame_end_index]

    atom_data = []
    for line in atoms_in_last_frame:
        atom_info = line.split()
        if len(atom_info) < 4:
            continue  # Skip incomplete lines
        atom_name = atom_info[0]
        coordinates = [float(coord) for coord in atom_info[1:4]]
        atom_data.append((atom_name, coordinates))
    return atom_data

def read_last_frame_cell_vectors(file_path):
    with open(file_path, 'r') as file:
        content = [line for line in file.readlines() if line.strip()]

    if not content:
        print(f"Warning: {file_path} is empty.")
        return np.identity(3)

    last_frame_line = content[-1]
    parts = last_frame_line.split()

    if len(parts) < 11:
        print(f"Warning: Unexpected format in {file_path}.")
        return np.identity(3)

    vector_a = [np.round(float(parts[2]), 10), np.round(float(parts[3]), 10), np.round(float(parts[4]), 10)]
    vector_b = [np.round(float(parts[5]), 10), np.round(float(parts[6]), 10), np.round(float(parts[7]), 10)]
    vector_c = [np.round(float(parts[8]), 10), np.round(float(parts[9]), 10), np.round(float(parts[10]), 10)]

    cell_vectors = np.array([vector_a, vector_b, vector_c])
    return cell_vectors

def save_xyz_file(atom_data, cell_vectors, output_path):
    num_atoms = len(atom_data)
    with open(output_path, 'w') as file:
        file.write(f"{num_atoms}\n")
        file.write(f'Lattice="{" ".join(map(str, cell_vectors.flatten()))}" Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        
        for atom in atom_data:
            atom_name = atom[0]
            coordinates = atom[1]
            formatted_coords = " ".join([f"{coord:.10f}" for coord in coordinates])
            file.write(f"{atom_name} {formatted_coords}\n")

def extract_last_stress_tensor(file_path):
    stress_tensors = []
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        # Traverse each line in the file
        for i, line in enumerate(lines):
            if 'STRESS| Analytical stress tensor [GPa]' in line:
                # Found the stress tensor section
                stress_matrix = []
                for j in range(i + 2, i + 5):  # Next 3 lines contain the tensor
                    if j >= len(lines):
                        break
                    row_data = lines[j].split()[2:]  # Skip first two columns
                    if len(row_data) < 3:
                        continue  # Skip incomplete lines
                    stress_matrix.append([np.round(float(val), 11) for val in row_data[:3]])
                if len(stress_matrix) == 3:
                    stress_tensors.append(np.array(stress_matrix))  # Add to list
    if stress_tensors:
        return stress_tensors[-1]
    else:
        print(f"Warning: No stress tensors found in {file_path}.")
        return np.zeros((3, 3))

def modify_cp2k_input(input_file_path, output_file_path, cell_vectors):
    with open(input_file_path, 'r') as file:
        input_lines = file.readlines()

    # Flag to check if inside &CELL section
    in_cell_section = False
    new_lines = []

    # Record indentation for &CELL
    cell_indent = None

    for line in input_lines:
        if "&CELL" in line:
            in_cell_section = True
            cell_indent = len(line) - len(line.lstrip())
            new_lines.append(line)
            continue

        if in_cell_section:
            if "&END CELL" in line:
                in_cell_section = False
                new_lines.append(line)
            else:
                stripped_line = line.strip()
                if stripped_line.startswith("A"):
                    new_lines.append(f"{' ' * (cell_indent + 2)}A  {cell_vectors[0,0]:.10f} {cell_vectors[0,1]:.10f} {cell_vectors[0,2]:.10f}\n")
                elif stripped_line.startswith("B"):
                    new_lines.append(f"{' ' * (cell_indent + 2)}B  {cell_vectors[1,0]:.10f} {cell_vectors[1,1]:.10f} {cell_vectors[1,2]:.10f}\n")
                elif stripped_line.startswith("C"):
                    new_lines.append(f"{' ' * (cell_indent + 2)}C  {cell_vectors[2,0]:.10f} {cell_vectors[2,1]:.10f} {cell_vectors[2,2]:.10f}\n")
                else:
                    new_lines.append(line)
        else:
            new_lines.append(line)

    # Write modified input to new file
    with open(output_file_path, 'w') as output_file:
        output_file.writelines(new_lines)

def remap_xyz(cell_vectors, coord_file, output_file):
    try:
        atoms = io.read(coord_file)
        atoms.set_cell(cell_vectors, scale_atoms=False)
        atoms.wrap()
        atoms.positions = np.round(atoms.positions, 10)
        atoms.write(output_file)
    except Exception as e:
        print(f"Error in remapping XYZ: {e}")

def symmetrize_and_save(array, filename="Cij.dat"):
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError("Cij must be a 2D square array!")
    symm_array = (array + array.T) / 2
    np.fill_diagonal(symm_array, np.diag(array))
    np.savetxt(filename, symm_array, fmt='%.3f')
    return symm_array

def process_epsilon(folder, epsilon, cell_vectors_cellopt, input_file_path, num_cores, cp2k_exec_path):
    cell_vectors_disp_pos = np.round(np.dot(cell_vectors_cellopt, (np.eye(3) + epsilon)), 10)
    cell_vectors_disp_neg = np.round(np.dot(cell_vectors_cellopt, (np.eye(3) - epsilon)), 10)
    
    pos_input_path = os.path.join(folder, "pos", input_file_path)
    neg_input_path = os.path.join(folder, "neg", input_file_path)
    
    # Modify input files with displaced cell vectors
    modify_cp2k_input(input_file_path, pos_input_path, cell_vectors_disp_pos)
    modify_cp2k_input(input_file_path, neg_input_path, cell_vectors_disp_neg)

    # Remap coordinates
    output_xyz_pos = os.path.join(folder, "pos", "coord.xyz")
    output_xyz_neg = os.path.join(folder, "neg", "coord.xyz")
    remap_xyz(cell_vectors_disp_pos, "./coord_cellopt.xyz", output_xyz_pos)
    remap_xyz(cell_vectors_disp_neg, "./coord_cellopt.xyz", output_xyz_neg)
