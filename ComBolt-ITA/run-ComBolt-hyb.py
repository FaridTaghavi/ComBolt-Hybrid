# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
#
# This file is part of a DFG-funded research project (grant no. 517518417).
# See the top-level LICENSE file for license details.

#!/usr/bin/env python3

import numpy as np
import frzout
import ROOT
ROOT.gSystem.Load("build/src/libsrc.so")  
from ROOT import particle_frzout, particle_urqmd, std
import array
import os
import subprocess
import glob
import shutil
from input_parameters import *


if INIT_MODE == 1:
    path_of_IS = "/home/farid/MyRepositories/KineticTheory/trento-master/events/event_set9/"
elif INIT_MODE == 2:
    path_of_IS = "/home/farid/MyRepositories/KineticTheory/kineticTheory/initial_state_AMPT/IS/"



#  ---------------------- Functions ----------------
def build_TrENTo_cmd(folder_path):

    # folder_path = temp_path+"trento_IS" 
    # Check if the directory exists
    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        shutil.rmtree(folder_path)
        print(f"Deleted folder: {folder_path}")
    # else:
    #     print(f"Folder does not exist: {folder_path}")
    
    cmd = [
        "../trento-master/build/src/trento",
        target, projectile,
        "--normalization", str(1),      # The normalization is set in the ComBolt initiale_state header
        "--grid-max", str(grid_max),
        "--grid-step", str(grid_step),
        "--output", folder_path,
        "--reduced-thickness", str(reduced_thickness_function),
        "--fluctuation", str(fluctuation),
        "--nucleon-width", str(nucleon_width),
        "--constit-width", str(constituent_width),
        "--constit-number", str(constituent_number),
        "--nucleon-min-dist", str(nucleon_minimum_distance),
        "--cross-section", str(cross_section),
        "--no-header",
        "--quiet"
    ]
    if not minimum_bias:
        cmd.extend([
                    "--b-min", str(minimum_impact_parameter),
                    "--b-max", str(maximum_impact_parameter)
                ])
    return cmd
def build_AMPT_cmd(folder_path):
    src = "../Ampt-v1.26t9b-v2.26t9b"
    dst = os.path.join(folder_path, "Ampt-v1.26t9b-v2.26t9b")
    shutil.copytree(src, dst, dirs_exist_ok=True)

    exec_path = os.path.join(dst, "exec")
    cmd = [exec_path]
    return cmd, dst  
    # Check if the directory exists
    # if os.path.exists(folder_path) and os.path.isdir(folder_path):
    #     shutil.rmtree(folder_path)
    #     print(f"Deleted folder: {folder_path}")
    
    
    
    # src = "../Ampt-v1.26t9b-v2.26t9b"
    # dst = os.path.join(folder_path, "Ampt-v1.26t9b-v2.26t9b")
    # shutil.copytree(src, dst, dirs_exist_ok=True)
    # # os.chdir( os.path.join(folder_path, "Ampt-v1.26t9b-v2.26t9b")  )
    # # exec_path = os.path.join(dst, "exec")
    # exec_path = os.path.join(folder_path, "Ampt-v1.26t9b-v2.26t9b", "exec")
    # cmd = [exec_path]
    # return cmd

def load_surface(filename):
    data = np.loadtxt(filename)

    if data.shape[1] != 12:
        raise ValueError("Expected 12 columns per line")

    x = data[:, 0:3]
    sigma = data[:, 3:6]
    v = data[:, 6:8]
    pi = {
        'xx': data[:, 8],
        'xy': data[:, 9],
        'yy': data[:, 10]
    }
    Pi = data[:, 11]

    return x, sigma, v, pi, Pi

def build_ComBolt_cmd(IS_FILE, init_mode):
    
    trento_file = ""
    ampt_file = ""
 
    if init_mode == 0:
        print("Non-fluctuating Gaussian mode!")
    elif init_mode == 1:
        trento_file = IS_FILE
        # print("TrENTo mode!")
    elif init_mode == 2:
        ampt_file = IS_FILE
        # print("AMPT mode!")
    
    cmd = [
        "./build/ComBolt-ITA",
        "--output", OUTPUT_DIR,
        "--quiet", str(QUITE_MODE),
        "--surface_output", SURFACE_DIR,
        "--trento_file", trento_file,
        "--save_BD", str(SAVE_BD),
        "--save_hydro", str(SAVE_HYDRO),
        "--save_Tmunu", str(SAVE_TMUNU),
        "--save_observables", str(SAVE_OBSERVABLES),
        "--nx", str(NX),
        "--xmax", str(XMAX),
        "--nvp", str(NVP),
        "--nvphi", str(NVPHI),
        "--Z2_symmetry", str(Z2_SYMMETRY),
        "--h", str(H),
        "--tau_max", str(TAU_MAX),
        "--tau0", str(TAU0),
        "--initialization_mode", str(INIT_MODE),
        "--momentum_isotropicity", str(MOM_ISO),
        "--energy_tot_tau0", str(ENERGY_TAU0),
        "--R0", str(R0),
        "--ellipticity", str(ELLIPTICITY),
        "--C0", str(C0),
        "--epsFO", str(EPSFO),
        "--dissipation_mode", str(DISS_MODE),
        "--gamma_hat", str(GAMMA_HAT),
        "--eta_over_s", str(ETA_S),
        "--grid_max", str(grid_max),
        "--grid_step", str(grid_step),
        "--norm", str(NORM),
        "--smear_sigmaV", str(SMEAR_SIGMA_V),
        "--smear_lambda", str(SMEAR_LAMBDA),
        "--smear_sigmaPhi", str(SMEAR_SIGMA_PHI),
        "--AMPT_norm", str(AMPT_NORM),
        "--AMPT_etas_cut", str(ETAS_CUT),
        "--AMPT_file", ampt_file,
        "--insert_threshold", str(INSERT_THRESHOLD),
        "--number_of_divistion", str(NUM_DIVISION),
        "--remove_threshold", str(REMOVE_THRESHOLD),
        "--number_of_cores_BoltzEq_solver", str(CORES_BOLTZ),
        "--number_of_cores_integration", str(CORES_INT),
        "--number_of_cores_interpolation", str(CORES_INTERP),
        "--number_of_cores_calc_obs", str(CORES_OBS),
        "--end_evol_last_frozen_cell", str(END_EVOL_LAST_FROZEN_CELL),
        "--calculate_freeze_out", str(CALCULATE_FREEZE_OUT)
    ]
    return cmd

# Function to find all files in a directory matching a specific pattern
def find_IS_files(directory):
    # Use glob to find files with a pattern like 000...00.dat, 000...01.dat, etc.
    pattern = os.path.join(directory, "*.dat")
    files = glob.glob(pattern)

    # Sort files numerically by their names (assuming the numbers are at the end of the filenames)
    files.sort(key=lambda x: int(os.path.basename(x).split('.')[0]))  # Extract the number part for sorting

    return files

def main():

    # Create a root file and a tree
    if particlization or afterburner_stage:
        if save_frzout_particles_ROOT_format or save_urqmd_particles_ROOT_format:
            file = ROOT.TFile(save_particles_path + "ComBolt-events.root", "RECREATE")
            file.SetCompressionLevel(1)  # faster I/O (less CPU usage, larger file)
            tree = ROOT.TTree("Events", "ComBolt_events")
            
            # Keep trace of the Combolt event ID and the sample ID
            eventID = np.zeros(1, dtype=int)
            tree.Branch("eventID", eventID, "eventID/I")
            
            # Marks events with same ComBolt evolution origin
            ComBolt_eventID = np.zeros(1, dtype=int)
            tree.Branch("ComBolt_eventID", ComBolt_eventID, "ComBolt_eventID/I")
            
            # Marks different samples of the same ComBolt event
            sampleID = np.zeros(1, dtype=int)
            tree.Branch("sampleID", sampleID, "sampleID/I")
    
        if save_frzout_particles_ROOT_format:
            particles_frzout = std.vector(ROOT.particle_frzout)()    
            tree.Branch("particles_frzout", particles_frzout)

        if save_urqmd_particles_ROOT_format:
            particles_urqmd = std.vector(ROOT.particle_urqmd)()    
            tree.Branch("particles_urqmd", particles_urqmd)

    # We generate many events
    event_id = 0
    for nevent in range(0, number_of_events):

        # => Generate initial state

        IS_FILE = ""
        # ---> TrENTo initil state
        if INIT_MODE == 1:
            print(f"\nGenerating TrENTo event {nevent} ...")
            try:
                trento_cmd = build_TrENTo_cmd(TEMP_PATH + "IS_data")
                subprocess.run(trento_cmd, check=True)
                IS_FILE = TEMP_PATH + "IS_data/0.dat"
                print("TrENTo event is generated.")
            except subprocess.CalledProcessError as e:
                print("TrENTo execution failed:", e)
        # ---> AMPT initial state
        elif INIT_MODE == 2:
            print(f"\nGenerating AMPT event {nevent} ...")
            AMPT_cmd, AMPT_folder = build_AMPT_cmd(TEMP_PATH + "IS_data")
            IS_FILE = TEMP_PATH + "IS_data/Ampt-v1.26t9b-v2.26t9b/ana/parton-initial-afterPropagation.dat"
            AMPT_FILE = IS_FILE
            subprocess.run(AMPT_cmd, cwd=AMPT_folder, check=True)
        else:
            print("Initial state mode is invalid.")

        # => Running ComBolt-ITA
        try:
            print("ComBolt-ITA is running ... (it takes a while)")
            combolt_cmd = build_ComBolt_cmd(IS_FILE, INIT_MODE)
            subprocess.run(combolt_cmd, check=True)
            print("ComBolt execution completed.")
        except subprocess.CalledProcessError as e:
            print("ComBolt execution failed:", e)


        # => Particlization and afterburner
        if particlization or afterburner_stage:
            if afterburner_stage:
                print("Duke particlization and UrQMD afterburner is running ...")
            else:    
                print("Duke particlization is running ...")
            
            SURFACE_FILE = f"{SURFACE_DIR}surface.dat"
            if not os.path.exists(SURFACE_FILE):
                print(f"Surface file not found at: {SURFACE_FILE}")
                print("Aborting event due to missing freeze-out surface file.")
                break 
    
            x, sigma, v, pi, Pi = load_surface(SURFACE_FILE)

            surface = frzout.Surface(x, sigma, v, pi=pi, Pi=Pi, ymax=ymax)
            hrg = frzout.HRG(swiching_temp, species='urqmd', res_width=True)
        
            urqmd_input      = os.path.join(TEMP_PATH, "urqmd_input.dat")
            osc2u_executable = os.path.join(afterburner_path, "build/osc2u/osc2u")
            urqmd_executable = os.path.join(afterburner_path, "build/urqmd/urqmd")
            urqmd_conf       = os.path.join(TEMP_PATH, "urqmd.conf")

            ComBolt_eventID[0] = nevent
            # Sampling particles
            for sample_itr in range(0, overSampling):
                # CLEAR VECTORS BEFORE ADDING NEW PARTICLES
                if save_frzout_particles_ROOT_format:
                    particles_frzout.clear()
    
                if save_urqmd_particles_ROOT_format:
                    particles_urqmd.clear()
                
                eventID[0] = event_id
                sampleID[0] = sample_itr


                if save_frzout_particles_plain_format:
                    input_path_filename = os.path.join(save_particles_path, f"p_in_{nevent}_{sample_itr}.dat")
                else:
                    input_path_filename = os.path.join(TEMP_PATH, f"p_in_{nevent}_{sample_itr}.dat")

                if save_urqmd_particles_plain_format:
                    output_path_filename = os.path.join(save_particles_path, f"p_out_{nevent}_{sample_itr}.dat")
                else:
                    output_path_filename = os.path.join(TEMP_PATH, f"p_out_{nevent}_{sample_itr}.dat")
        
                particles = frzout.sample(surface, hrg)
                if particles.size == 0:
                    print("No particles sampled, skipping this iteration.")
                    continue    

                with open(input_path_filename, 'w') as f:
                    print('#', particles.size, file=f)
                    for p in particles:
                        print(p['ID'], *p['x'], *p['p'], file=f)

                        if save_frzout_particles_ROOT_format:
                            particle_frzout = ROOT.particle_frzout(
                                                            int(p['ID']),
                                                            float(p['x'][0]),
                                                            float(p['x'][1]),
                                                            float(p['x'][2]),
                                                            float(p['x'][3]),
                                                            float(p['p'][0]),
                                                            float(p['p'][1]),
                                                            float(p['p'][2]),
                                                            float(p['p'][3])
                                                        )
                            particles_frzout.push_back(particle_frzout)
                
                if afterburner_stage:
                    with open(urqmd_input, 'w') as out_file, open(input_path_filename, 'r') as in_file:
                        subprocess.run([osc2u_executable], 
                                       stdin=in_file, 
                                       stdout=out_file, 
                                       stderr=subprocess.DEVNULL, 
                                       check=True,
                                       cwd=TEMP_PATH  # Run the command *inside* temp_dir
                                       )

                    # Create urqmd.conf if not exists
                    if not os.path.exists(urqmd_conf):
                        with open(urqmd_conf, 'w') as f:
                            f.write("tim 8000 8000\ncto 40 2\nf13\nf14\nf15\nf16\nf18\nf19\nf20\nxxx")
                    # Set environment and run URQMD silently
                    env = os.environ.copy()
                    env['ftn09'] = urqmd_conf
                    env['ftn10'] = urqmd_input
                    env['ftn30'] = output_path_filename

                    # Run UeQMD
                    result = subprocess.run([urqmd_executable], 
                                   env=env, 
                                   stdout=subprocess.DEVNULL, 
                                   stderr=subprocess.DEVNULL,
                                   cwd=TEMP_PATH  # Run the command *inside* temp_dir
                                   ) # subprocess.DEVNULL: makes the UrQMD output silent           

                    if result.returncode != 0:
                        print(" UrQMD failed with return code:", result.returncode)
                    # Save UrQMD output in ROOT format
                    if save_urqmd_particles_ROOT_format:
                
                        with open(output_path_filename, 'r') as f:
                            for line in f:
                                if line.startswith("#"):
                                    continue
                                parts = line.split()
                                if len(parts) < 8:
                                    continue
                                # The output format is: 0. pid 1. charge 2. pT 3. ET 4. mT 5. phi 6. y 7. eta
                                particle_urqmd = ROOT.particle_urqmd(
                                                                    int(parts[0]),   
                                                                    int(parts[1]),  
                                                                    float(parts[2]), 
                                                                    float(parts[3]),
                                                                    float(parts[4]),
                                                                    float(parts[5]),
                                                                    float(parts[6]),
                                                                    float(parts[7])
                                                                )
                                particles_urqmd.push_back(particle_urqmd)
                
                # Fill the TTree
                event_id += 1
                if particlization or afterburner_stage:
                    if save_frzout_particles_ROOT_format or save_urqmd_particles_ROOT_format:
                        tree.Fill()
                        # Save the result after all over samplings to avoid loosing data due to intruption.
                        if event_id % overSampling == 0:
                            tree.AutoSave("SaveSelf")
            
            if afterburner_stage:
                print("Duke particlization and UrQMD afterburner execution completed.")
            else:    
                print("Duke particlization is done.")
        

        # Clean up the temporary folder
        if os.path.exists(TEMP_PATH) and os.path.isdir(TEMP_PATH):
            shutil.rmtree(TEMP_PATH)

    
    if particlization or afterburner_stage:
        if save_frzout_particles_ROOT_format or save_urqmd_particles_ROOT_format:
            file.Write()
            file.Close()


# Run the script
if __name__ == "__main__":
    main()






