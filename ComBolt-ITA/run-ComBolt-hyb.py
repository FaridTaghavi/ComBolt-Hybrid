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
from ROOT import particle_frzout, particle_urqmd, initial_info, std
import array
import os
import subprocess
import glob
import shutil
import matplotlib.pyplot as plt

import importlib.util
import sys






if len(sys.argv) < 2:
    print("Usage: python3 script.py path/to/input_parameters.py [TEMP_PATH] [OUTPUT_DIR] [OUTPUT_FILE]")
    sys.exit(1)

param_file_path = sys.argv[1]
temp_path = sys.argv[2] if len(sys.argv) > 2 else None
out_path = sys.argv[3] if len(sys.argv) > 3 else None
out_file = sys.argv[4] if len(sys.argv) > 4 else None

# Dynamically import input_parameters.py
spec = importlib.util.spec_from_file_location("input_parameters", param_file_path)
input_parameters = importlib.util.module_from_spec(spec)
spec.loader.exec_module(input_parameters)

# Use TEMP_PATH and OUT_PATH from command line if provided, otherwise default
if temp_path is None:
    print("The TEMP_PATH is set from input parameter file.")
    temp_path = input_parameters.TEMP_PATH
if out_path is None:
    print("The OUTPUT_PATH is set from input parameter file.")
    out_path = input_parameters.OUTPUT_DIR
if out_file is None:
    print("The OUTPUT_FILE is ComBolt-events.root.")
    out_file = input_parameters.ROOTFILE_NAME


# # Get path to input_parameters.py from command-line
# if len(sys.argv) < 2:
#     print("Usage: python3 script.py path/to/input_parameters.py")
#     sys.exit(1)
# 
# param_file_path = sys.argv[1]
# 
# # Dynamically import the module
# spec = importlib.util.spec_from_file_location("input_parameters", param_file_path)
# input_parameters = importlib.util.module_from_spec(spec)
# spec.loader.exec_module(input_parameters)



if input_parameters.INIT_MODE == 1:
    path_of_IS = "/home/farid/MyRepositories/KineticTheory/trento-master/events/event_set9/"
elif input_parameters.INIT_MODE == 2:
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
   
    # cmd = [
    #     "../trento-master/build/src/trento",
    # ]
    # 
    # # Decide whether to use target/projectile names or input files
    # if input_parameters.use_files:  # you can add this flag in your class
    #     cmd.extend([input_parameters.file1, input_parameters.file2])
    # else:
    #     cmd.extend([input_parameters.target, input_parameters.projectile])
    # 
    # # Now add the rest of the options
    # cmd.extend([
    #     "--normalization", str(1),
    #     "--grid-max", str(input_parameters.grid_max),
    #     "--grid-step", str(input_parameters.grid_step),
    #     "--output", folder_path,
    #     "--reduced-thickness", str(input_parameters.reduced_thickness_function),
    #     "--fluctuation", str(input_parameters.fluctuation),
    #     "--nucleon-width", str(input_parameters.nucleon_width),
    #     "--constit-width", str(input_parameters.constituent_width),
    #     "--constit-number", str(input_parameters.constituent_number),
    #     "--nucleon-min-dist", str(input_parameters.nucleon_minimum_distance),
    #     "--cross-section", str(input_parameters.cross_section),
    #     "--no-header",
    #     "--quiet"
    # ])

 
    cmd = [
        "../trento-master/build/src/trento",
        input_parameters.target, input_parameters.projectile,
        # "--normalization", str(1),      # The normalization is set in the ComBolt initiale_state header
        "--normalization", str(input_parameters.trento_norm),  
        "--grid-max", str(input_parameters.grid_max),
        "--grid-step", str(input_parameters.grid_step),
        "--output", folder_path,
        "--reduced-thickness", str(input_parameters.reduced_thickness_function),
        "--fluctuation", str(input_parameters.fluctuation),
        "--nucleon-width", str(input_parameters.nucleon_width),
        "--constit-width", str(input_parameters.constituent_width),
        "--constit-number", str(input_parameters.constituent_number),
        "--nucleon-min-dist", str(input_parameters.nucleon_minimum_distance),
        "--cross-section", str(input_parameters.cross_section),
        "--no-header",
        "--quiet"
    ]
    if not input_parameters.minimum_bias:
        cmd.extend([
                    "--b-min", str(input_parameters.minimum_impact_parameter),
                    "--b-max", str(input_parameters.maximum_impact_parameter)
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
        "--output", out_path,
        "--quiet", str(input_parameters.QUITE_MODE),
        "--surface_output", temp_path,
        "--trento_file", trento_file,
        "--save_BD", str(input_parameters.SAVE_BD),
        "--save_hydro", str(input_parameters.SAVE_HYDRO),
        "--save_Tmunu", str(input_parameters.SAVE_TMUNU),
        "--save_observables", str(input_parameters.SAVE_OBSERVABLES),
        "--save_init_info", str(input_parameters.SAVE_INIT_INFO),
        "--nx", str(input_parameters.NX),
        "--xmax", str(input_parameters.XMAX),
        "--nvp", str(input_parameters.NVP),
        "--nvphi", str(input_parameters.NVPHI),
        "--Z2_symmetry", str(input_parameters.Z2_SYMMETRY),
        "--h", str(input_parameters.H),
        "--tau_max", str(input_parameters.TAU_MAX),
        "--tau0", str(input_parameters.TAU0),
        "--initialization_mode", str(input_parameters.INIT_MODE),
        "--momentum_isotropicity", str(input_parameters.MOM_ISO),
        "--energy_tot_tau0", str(input_parameters.ENERGY_TAU0),
        "--R0", str(input_parameters.R0),
        "--ellipticity", str(input_parameters.ELLIPTICITY),
        "--C0", str(input_parameters.C0),
        "--epsFO", str(input_parameters.EPSFO),
        "--TFO", str(input_parameters.TEMP_FO),
        "--EOS_mode", str(input_parameters.EOS_MODE),
        "--dissipation_mode", str(input_parameters.DISS_MODE),
        "--eta_over_s_min", str(input_parameters.ETA_OVER_S_MIN),
        "--eta_over_s_slope", str(input_parameters.ETA_OVER_S_SLOPE),
        "--eta_over_s_pow", str(input_parameters.ETA_OVER_S_POW),
        "--eta_over_s_Tc", str(input_parameters.ETA_OVER_S_TC),
        "--grid_max", str(input_parameters.grid_max),
        "--grid_step", str(input_parameters.grid_step),
        "--norm", str(input_parameters.NORM),
        "--smear_sigmaV", str(input_parameters.SMEAR_SIGMA_V),
        "--smear_lambda", str(input_parameters.SMEAR_LAMBDA),
        "--smear_sigmaPhi", str(input_parameters.SMEAR_SIGMA_PHI),
        "--AMPT_norm", str(input_parameters.AMPT_NORM),
        "--AMPT_etas_cut", str(input_parameters.ETAS_CUT),
        "--AMPT_file", ampt_file,
        "--insert_threshold", str(input_parameters.INSERT_THRESHOLD),
        "--number_of_divistion", str(input_parameters.NUM_DIVISION),
        "--remove_threshold", str(input_parameters.REMOVE_THRESHOLD),
        "--number_of_cores_BoltzEq_solver", str(input_parameters.CORES_BOLTZ),
        "--number_of_cores_integration", str(input_parameters.CORES_INT),
        "--number_of_cores_interpolation", str(input_parameters.CORES_INTERP),
        "--number_of_cores_calc_obs", str(input_parameters.CORES_OBS),
        "--end_evol_last_frozen_cell", str(input_parameters.END_EVOL_LAST_FROZEN_CELL),
        "--calculate_freeze_out", str(input_parameters.CALCULATE_FREEZE_OUT)
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

def main(temp_path, out_path, out_file):
    print(temp_path, out_path, out_file)
    # Create a root file and a tree
    if input_parameters.particlization or input_parameters.afterburner_stage:
        if input_parameters.save_frzout_particles_ROOT_format or input_parameters.save_urqmd_particles_ROOT_format:



            os.makedirs(out_path, exist_ok=True)  # Create dir if it doesn't exist

            filename = os.path.join(out_path, out_file + ".root")

            print(f"Opening ROOT file: '{filename}'")
            file = ROOT.TFile(filename, "RECREATE")


            # file = ROOT.TFile(input_parameters.save_particles_path + "ComBolt-events.root", "RECREATE")
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
    
        if input_parameters.save_init_info_ROOT_format:
            initial_information = std.vector(ROOT.initial_info)()    
            tree.Branch("initial_information", initial_information)
        
        if input_parameters.save_frzout_particles_ROOT_format:
            particles_frzout = std.vector(ROOT.particle_frzout)()    
            tree.Branch("particles_frzout", particles_frzout)

        if input_parameters.save_urqmd_particles_ROOT_format:
            particles_urqmd = std.vector(ROOT.particle_urqmd)()    
            tree.Branch("particles_urqmd", particles_urqmd)

    # We generate many events
    event_id = 0
    for nevent in range(0, input_parameters.number_of_events):

        # => Generate initial state

        IS_FILE = ""
        # ---> TrENTo initil state
       

        if input_parameters.INIT_MODE == 1:
             print(f"\nGenerating TrENTo event {nevent} ...")
             try:
                 trento_cmd = build_TrENTo_cmd(temp_path + "IS_data")
                 print(trento_cmd)
                 subprocess.run(trento_cmd, check=True)
                 IS_FILE = temp_path + "IS_data/0.dat"
                 print("TrENTo event is generated.")
             except subprocess.CalledProcessError as e:
                 print("TrENTo execution failed:", e)
        
        # if input_parameters.INIT_MODE == 1:
        #     print(f"\nGenerating TrENTo event {nevent} ...")
        #     try:
        #        
        #         is_data_dir = os.path.join(temp_path, "IS_data")
        #         #  print(f"Creating directory: {is_data_dir}")
        #         os.makedirs(is_data_dir, exist_ok=True)
        #         # print(f"Directory exists now? {os.path.isdir(is_data_dir)}")
 
        #         trento_cmd = build_TrENTo_cmd(temp_path + "IS_data")
        #         print("Running command:", ' '.join(trento_cmd))
        # 
        #         # Optionally: open log file
        #         # with open("trento_output.log", "w") as logfile:
        #         proc = subprocess.Popen(
        #             trento_cmd,
        #             stdout=subprocess.PIPE,
        #             stderr=subprocess.STDOUT,
        #             text=True
        #         )


        #         
        #         for line in proc.stdout:
        #             print(line, end='')  # Print TrENTo output live
        #             # logfile.write(line)  # Optionally save to file
        # 
        #         proc.wait()
        #         if proc.returncode != 0:
        #             raise subprocess.CalledProcessError(proc.returncode, trento_cmd)
        # 
        #         IS_FILE = temp_path + "IS_data/0.dat"
        #         print("TrENTo event is generated.")
        # 
        #     except subprocess.CalledProcessError as e:
        #         print("\n[ERROR] TrENTo execution failed with exit code:", e.returncode)
        #     except Exception as e:
        #         print("\n[ERROR] Unexpected error running TrENTo:", e)
        


        #---> AMPT initial state
        elif input_parameters.INIT_MODE == 2:
            print(f"\nGenerating AMPT event {nevent} ...")
            AMPT_cmd, AMPT_folder = build_AMPT_cmd(temp_path + "IS_data")
            IS_FILE = temp_path + "IS_data/Ampt-v1.26t9b-v2.26t9b/ana/parton-initial-afterPropagation.dat"
            AMPT_FILE = IS_FILE
            subprocess.run(AMPT_cmd, cwd=AMPT_folder, check=True)
        else:
            print("Initial state mode is invalid.")

        # => Running ComBolt-ITA
        try:
            print("ComBolt-ITA is running ... (it takes a while)")
            combolt_cmd = build_ComBolt_cmd(IS_FILE, input_parameters.INIT_MODE)
            subprocess.run(combolt_cmd, check=True)
            print("ComBolt execution completed.")
        except subprocess.CalledProcessError as e:
            print("ComBolt execution failed:", e)
        
        # => Save initial information in root

        if input_parameters.save_init_info_ROOT_format:
            initial_info_path_filename = os.path.join(temp_path, f"init_inf.dat")
            with open(initial_info_path_filename, 'r') as f:
                for line in f:
                    parts = line.split()
                    init_inf = ROOT.initial_info(
                                                     float(parts[0]),
                                                     float(parts[1]),
                                                     float(parts[2]),
                                                     float(parts[3]),
                                                     float(parts[4]),
                                                     float(parts[5]),
                                                     float(parts[6]),
                                                     float(parts[7]),
                                                     float(parts[8]),
                                                     float(parts[9]),
                                                     float(parts[10]),
                                                     float(parts[11]),
                                                     float(parts[12]),
                                                     float(parts[13])
                                                 )
                    initial_information.push_back(init_inf)

        # => Particlization and afterburner
        if input_parameters.particlization or input_parameters.afterburner_stage:
            if input_parameters.afterburner_stage:
                print("Duke particlization and UrQMD afterburner is running ...")
            else:    
                print("Duke particlization is running ...")
            
            SURFACE_FILE = f"{temp_path}surface.dat"
            if not os.path.exists(SURFACE_FILE):
                print(f"Surface file not found at: {SURFACE_FILE}")
                print("Aborting event due to missing freeze-out surface file.")
                break 
    
            x, sigma, v, pi, Pi = load_surface(SURFACE_FILE)

            surface = frzout.Surface(x, sigma, v, pi=pi, Pi=Pi, ymax=input_parameters.ymax)
            hrg = frzout.HRG(input_parameters.swiching_temp, species='urqmd', res_width=True, decay_f500=True)
            # hrg = frzout.HRG(input_parameters.swiching_temp, species='urqmd', res_width=True)
        
            urqmd_input      = os.path.abspath( os.path.join(temp_path, "urqmd_input.dat")  )
            osc2u_executable = os.path.abspath( os.path.join(input_parameters.afterburner_path, "build/osc2u/osc2u")  )
            urqmd_executable = os.path.abspath( os.path.join(input_parameters.afterburner_path, "build/urqmd/urqmd")  )
            urqmd_conf       = os.path.abspath( os.path.join(temp_path, "urqmd.conf")  )

            ComBolt_eventID[0] = nevent
            # Sampling particles
            for sample_itr in range(0, input_parameters.overSampling):
                # CLEAR VECTORS BEFORE ADDING NEW PARTICLES
                if input_parameters.save_frzout_particles_ROOT_format:
                    particles_frzout.clear()
    
                if input_parameters.save_urqmd_particles_ROOT_format:
                    particles_urqmd.clear()
                
                eventID[0] = event_id
                sampleID[0] = sample_itr


                if input_parameters.save_frzout_particles_plain_format:
                    input_path_filename = os.path.join(save_particles_path, f"p_in_{nevent}_{sample_itr}.dat")
                else:
                    input_path_filename = os.path.join(temp_path, f"p_in_{nevent}_{sample_itr}.dat")

                if input_parameters.save_urqmd_particles_plain_format:
                    output_path_filename = os.path.join(save_particles_path, f"p_out_{nevent}_{sample_itr}.dat")
                else:
                    output_path_filename = os.path.join(temp_path, f"p_out_{nevent}_{sample_itr}.dat")
        
                particles = frzout.sample(surface, hrg)
                if particles.size == 0:
                    print("No particles sampled, skipping this iteration.")
                    continue    

                with open(input_path_filename, 'w') as f:
                    print('#', particles.size, file=f)
                    for p in particles:
                        print(p['ID'], *p['x'], *p['p'], file=f)

                        if  input_parameters.save_frzout_particles_ROOT_format:
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


                            # phiphi = np.arctan2( p['p'][2], p['p'][1] ) + np.pi
                            # # make histogram
                            # plt.hist(phiphi, bins=36, range=(0, 2*np.pi), density=True, histtype='step', linewidth=2)
                            # 
                            # # nice labels
                            # plt.xlabel(r'$\phi$ [rad]')
                            # plt.ylabel('Normalized counts')
                            # plt.title('Azimuthal distribution of particles')
                            # 
                            # plt.xlim(0, 2*np.pi)
                            # print(phiphi)
                            # # save figure
                            # plt.savefig("/home/ktas/ge57vag/02-04-Aug-2025-ComBolt-ITA/ROOT-macro/01-pTSpectrum/phi_distribution_fout.png", dpi=300, bbox_inches="tight")
                            # plt.close() 
                
                if input_parameters.afterburner_stage:
                    with open(urqmd_input, 'w') as out_file, open(input_path_filename, 'r') as in_file:
                        subprocess.run([osc2u_executable], 
                                       stdin=in_file, 
                                       stdout=out_file, 
                                       stderr=subprocess.DEVNULL, 
                                       check=True,
                                       cwd=temp_path  # Run the command *inside* temp_dir
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
                                   cwd=temp_path  # Run the command *inside* temp_dir
                                   ) # subprocess.DEVNULL: makes the UrQMD output silent           

                    if result.returncode != 0:
                        print(" UrQMD failed with return code:", result.returncode)
                    # Save UrQMD output in ROOT format
                    if input_parameters.save_urqmd_particles_ROOT_format:
                
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
                                # # make histogram
                                # plt.hist(parts[5], bins=36, range=(0, 2*np.pi), density=True, histtype='step', linewidth=2)
                                # 
                                # # nice labels
                                # plt.xlabel(r'$\phi$ [rad]')
                                # plt.ylabel('Normalized counts')
                                # plt.title('Azimuthal distribution of particles')
                                # 
                                # plt.xlim(0, 2*np.pi)
                                # print(parts[5])
                                # # save figure
                                # plt.savefig("/home/ktas/ge57vag/02-04-Aug-2025-ComBolt-ITA/ROOT-macro/01-pTSpectrum/phi_distribution_pyth.png", dpi=300, bbox_inches="tight")
                                # plt.close() 
                # Fill the TTree
                event_id += 1
                if input_parameters.particlization or input_parameters.afterburner_stage:
                    if input_parameters.save_frzout_particles_ROOT_format or input_parameters.save_urqmd_particles_ROOT_format:
                        tree.Fill()
                        # Save the result after all over samplings to avoid loosing data due to intruption.
                        if event_id % input_parameters.overSampling == 0:
                            tree.AutoSave("SaveSelf")
            
            if input_parameters.afterburner_stage:
                print("Duke particlization and UrQMD afterburner execution completed.")
            else:    
                print("Duke particlization is done.")
        

        # Clean up the temporary folder
        if os.path.exists(temp_path) and os.path.isdir(temp_path):
           shutil.rmtree(temp_path)

    
    if input_parameters.particlization or input_parameters.afterburner_stage:
        if input_parameters.save_frzout_particles_ROOT_format or input_parameters.save_urqmd_particles_ROOT_format:
            # file.Write()    #AutoSave should handel this
            file.Close()


# Run the script
if __name__ == "__main__":
    # Run main with out_file
    main(temp_path,out_path, out_file)






