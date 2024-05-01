#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 14:46:58 2023

@author: akshat
"""

import os 
import uuid
import time
import subprocess
import itertools
import argparse
import multiprocessing    

def read_config_file(filename):
    params = {}
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line and not line.startswith("#"):  # skip comments and empty lines
                key, value = line.split("=", 1)
                if value.isdigit():
                    params[key] = int(value)
                else:
                    params[key] = value
    return params

# Read the all.ctrl file to get parameters
config_params = read_config_file("all.ctrl")

RECEPTOR_LOCATION  = config_params.get("RECEPTOR_LOCATION", "./DATA/docking_receptor.pdbqt")
EXHAUSTIVENESS     = str( config_params.get("EXHAUSTIVENESS", "1"))
CENTER_X           = str( config_params.get("CENTER_X"))
CENTER_Y           = str( config_params.get("CENTER_Y"))
CENTER_Z           = str( config_params.get("CENTER_Z"))
SIZE_X             = str( config_params.get("SIZE_X"))
SIZE_Y             = str( config_params.get("SIZE_Y"))
SIZE_Z             = str( config_params.get("SIZE_Z"))
SMILES_FILE        = str( config_params.get("SMILES_FILES"))
MAX_NUM_JOBS       = config_params.get("MAX_NUM_JOBS")
MAX_NUM_MOLS       = config_params.get("NUM_MOLS")
DOCKING_SCORE_THRS = float(config_params.get("DOCKING_SCORE_THRESHOLD"))

def generate_unique_file_name(base_name, extension):
    timestamp = int(time.time() * 1000)
    unique_id = uuid.uuid4().hex
    file_name = f"{base_name}_{timestamp}_{unique_id}.{extension}"
    return file_name


def check_energy(lig_): 
    """
    Check the quality of a generated structure by computing its total energy using the Open Babel obenergy tool.
    Parameters:
        lig_ (str): the name of the ligand file in PDBQT format.
    Returns:
        total_energy (float): the computed total energy of the ligand in Kcal/mol.
    """
    # Check the quality of generated structure (some post-processing quality control):
    try: 
        ob_cmd = ['obenergy', lig_]
        command_obabel_check = subprocess.run(ob_cmd, capture_output=True)
        command_obabel_check = command_obabel_check.stdout.decode("utf-8").split('\n')[-2]
        total_energy         = float(command_obabel_check.split(' ')[-2])
    except: 
        
        total_energy = 10000 # Calculation has failed. 
        
    return total_energy


def run_docking(lig_location, out_location, method='qvina'): 
    """
    Perform molecular docking with a specific methods (QuickVina/Smina) on the 6Y2F protein. 
    An exhaustiveness of 10 is used for the QuickVina calculations, while an 
    exhaustivesness of 100 is used for a smina calculation 

    Parameters
    ----------
    method : str, The calculation type to be run qvina/smina 

    Returns
    -------
    (float) Docking score.
    """
    if method == 'qvina': 
        command_run = subprocess.run(["./DATA/qvina", "--receptor", RECEPTOR_LOCATION, "--ligand", lig_location, "--center_x", CENTER_X, "--center_y", CENTER_Y, "--center_z", CENTER_Z, "--size_x", SIZE_X, "--size_y", SIZE_Y, "--size_z", SIZE_Z, "--exhaustiveness", EXHAUSTIVENESS, "--out", out_location], capture_output=True)
    elif method == 'smina': 
        command_run = subprocess.run(["./DATA/smina", "--receptor", RECEPTOR_LOCATION, "--ligand", lig_location, "--center_x", CENTER_X, "--center_y", CENTER_Y, "--center_z", CENTER_Z, "--size_x", SIZE_X, "--size_y", SIZE_Y, "--size_z", SIZE_Z, "--exhaustiveness", EXHAUSTIVENESS, "--out", out_location], capture_output=True)
    
    # Note: TO BE SPECIFIED BY USER: Any new docking method can be specified by the user
    # Ensure to update the conditional block accordingly.
    # If a new method is added, it should be handled explicitly, otherwise, an exception will be raised, indicating the available options.
    
    else: 
        raise Exception('Possible docking softwares: qvina/smina')

    # Ensure the pose of the output molecule is not broken: 
    pose_energy = check_energy(out_location)
    if pose_energy == 10000: # broken molecule generated (docking failed)
        return 10000
        
    # Obtain the docking score: 
    command_run = command_run.stdout.decode("utf-8").split('\n')

    docking_score = []
    for item in command_run: 
        line_split = item.split(' ')
        line_split = [x for x in line_split if x != '']
        if len(line_split) == 4: 
            try: 
                _ = float(line_split[0])
                vr_2 = float(line_split[1])
                _ = float(line_split[2])
                _ = float(line_split[3])
                docking_score.append(vr_2)
            except: 
                continue
    docking_score = min(docking_score)

    return docking_score




def perform_calc_single(args): 
    
    out_location    = generate_unique_file_name('pose', 'pdbqt') # For the docking pose
    output_filename = generate_unique_file_name('lig', 'pdbqt')  # For the 3D ligand (obabel converted smi)
    
    try: 
        smi, chunk_1, chunk_2, enamine_id = args
        
        # print('smi: {} fname: {}'.format(smi, output_filename))
        cmd = ["obabel", "-ismi","-:" + smi,"-O", output_filename, "--gen3d"]
        # subprocess.run(cmd, timeout=120)
        with open(os.devnull, 'w') as devnull:
            subprocess.run(cmd, stdout=devnull, stderr=devnull, timeout=120)

        # Ensure a stable molecule: 
        lig_energy = check_energy(output_filename)
    
        # Specifying docking input file & output file: 
        lig_location = output_filename
        
        # Perform docking: 
        if lig_energy < 10000: 
            score_3 = run_docking(lig_location, out_location, method='qvina')
                
        if score_3 > DOCKING_SCORE_THRS: 
            if os.path.exists(out_location):
                os.system('rm {}'.format(out_location))
            if os.path.exists(output_filename):
                os.system('rm {}'.format(output_filename))
            with open('./OUTPUT_{}_{}.txt'.format(chunk_1, chunk_2), 'a+') as f: 
                f.writelines(['{}, {}, {}\n'.format(smi, score_3, enamine_id)]) 
        else: 
            with open('./OUTPUT_{}_{}.txt'.format(chunk_1, chunk_2), 'a+') as f: 
                f.writelines(['{}, {}, {}, {}, {}\n'.format(smi, score_3, enamine_id, lig_location, out_location)]) 
            
            # Move the files lig_location, out_location to the directory OUTPUTS/ there is no need to keep the file in the original location anymore
            os.system('mv {} OUTPUTS/'.format(lig_location))
            os.system('mv {} OUTPUTS/'.format(out_location))

    except: 

        if os.path.exists(out_location):
            os.system('rm {}'.format(out_location))
        if os.path.exists(output_filename):
            os.system('rm {}'.format(output_filename))
        
        with open('./OUTPUT_{}_{}.txt'.format(chunk_1, chunk_2), 'a+') as f: 
            f.writelines(['{}, {}, {}\n'.format(smi, 10000, 10000)]) 


def main(filename, chunk_1, chunk_2):
    
    smiles_all     = []
    enamine_id_all = []
    with open(filename, 'r') as f:
        for line_number, line in enumerate(f):
            # Skip the header row
            if line_number == 0:
                continue
    
            # Read only the lines within the desired range
            if chunk_1 <= line_number < chunk_2:
                smiles = line.split('\t')[0]
                smiles_all.append(smiles)
                enamine_id_all.append( line.split('\t')[1] )
            elif line_number >= chunk_2:
                # Break the loop once you've read up to chunk_2
                break
    
    print('Num smiles:', len(smiles_all))
    data = [(smiles, chunk_1, chunk_2, enamine_id_all[i]) for i,smiles in enumerate(smiles_all)]
    
    # pool object with number of element
    pool = multiprocessing.Pool()
    
    # Parallel time: 
    start_time = time.time()
    pool.map(perform_calc_single, data)
    total_time = time.time() - start_time
    print('Total Time: ', total_time)
    
        
    
    
parser = argparse.ArgumentParser()
parser.add_argument("job_id", help="Array ID: SLURM_ARRAY_TASK_ID")
args = parser.parse_args()
job_idx = int(args.job_id)    


ratio = (MAX_NUM_MOLS+100000) // MAX_NUM_JOBS # Number of molecules that will be processed by each subjob

start_idx = (job_idx-1) * ratio
end_idx   = start_idx + ratio

print('Start idx: {} End idx: {}'.format(start_idx, end_idx))

main(SMILES_FILE, start_idx, end_idx)

    
    
    
