# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
#
# This file is part of a DFG-funded research project (grant no. 517518417).
# See the top-level LICENSE file for license details.

# input_parameters.py

# =============  ComBolt-ITA parameters =================

# Spatial lattice setup (Gaussian mode and AMPT mode)
NX = 150
XMAX = 15.0  # [fm]

# Momentum lattice setup
NVP = 15
NVPHI = 20
Z2_SYMMETRY = True # The Boltzmann distribution is symmetric for pozitive and negative vz.

# Numerical method parameters
H = 0.002       # time step [fm]
ALPHA_H = 0.003 # time steps stiffness
TAU_MAX = 5.0  # [fm] maximum time

# Initial conditions
TAU0 = 1.0e-1        # [fm] initial time 
INIT_MODE = 1          # Initial condition mode: 0. Gaussian 1. External energy density (from TrENTo), 2. external Boltz Dist From AMPT
MOM_ISO = 0.1         # Initial momentum anisotropy: MOM_ISO -> 0 highly anisotropic, MOM_ISO -> infinity close to isotropy


# Gaussian spatial distribution
ENERGY_TAU0 = 1280     # [GeV / fm] initial total transverse energy * initial time (Gaussian mode)
R0 = 1.97              # [fm] initial system size (Gaussian mode)
ELLIPTICITY = 0.0; # 0.416    # initial ellipticity (Gaussian mode)

# Trento settings
TRENTO_FILE = "../trento-master/events/event_set7/0.dat"
# GRID_MAX = 15.0     # TrENTo
# GRID_STEP = 0.2
# NORM = 50.0         # Scaling the initial TrENTo energy density (The nomalization in TrENTois set  to 1)
NORM = 1.0         # Scaling the initial TrENTo energy density (The nomalization in TrENTois set  to 1)

# AMPT settings
SMEAR_SIGMA_V = 0.5
SMEAR_LAMBDA = 10.5
SMEAR_SIGMA_PHI = 10.5
AMPT_NORM = 0.5
ETAS_CUT = 2.0
AMPT_FILE = "initial_state_AMPT/0002.dat"



# Medium properties  
EOS_MODE = 1    # 0: conformal EOS, 1: lattice QCD EOS
C0 = 4.27                                    # Conformal equation of state: energy = C0 T^4
TEMP_FO = 0.155                              # [GeV] Freeze-out temperature
EPSFO = C0 * TEMP_FO**4 / 0.197**3           # [GeV/fm^3] Freeze-out energy density 
# EPSFO = 0.322213                           # [GeV/fm^3] Freeze-out energy density 
# TEMP_FO = (EPSFO * 0.197**3 / C0)**(1/4)   # [GeV] Freeze-out temperature


# Dissipation
DISS_MODE = 0    # 0: eta/s constant, 1: opacity constant, 2:  constant eta / s with tau_relax(ed) from lattice.
GAMMA_HAT = 0.08
ETA_OVER_S_MIN = 0.18
ETA_OVER_S_SLOPE = 0.0    # in GeV^-1
ETA_OVER_S_POW = 0.0
ETA_OVER_S_TC = TEMP_FO     # in GeV


# Interpolation
INSERT_THRESHOLD = 0.15    # Dividie segments in the vz directon if it is greater thatn INSERT_THRESHOLD
NUM_DIVISION = 2           # Number of divisions
REMOVE_THRESHOLD = 5e-7    # Remove the vz node if its distane to vz = 0 is REMOVE_THRESHOLD

# Saving options
# OUTPUT_DIR = "/scratch8/ge57vag/26-10Sep2025-ComBolt-OO/evolution/"
OUTPUT_DIR = "/scratch8/ge57vag/26-10Sep2025-ComBolt-OO/"
TEMP_PATH =  "/scratch8/ge57vag/26-10Sep2025-ComBolt-OO/temp/"
ROOTFILE_NAME = "ComBolt-events"
SURFACE_DIR = TEMP_PATH
SAVE_BD          = False
SAVE_HYDRO       = False
SAVE_TMUNU       = False
SAVE_OBSERVABLES = False  
SAVE_INIT_INFO   = True   # if you want to save the initial info in ROOT, keep it true.
 

# Multithreading
CORES_BOLTZ = 40
CORES_INT = 40
CORES_INTERP = 40
CORES_OBS = 40

# Freeze-out parameters
END_EVOL_LAST_FROZEN_CELL = True # True: end evolution of the last frozen cell, False: continue evolution until tau_max
CALCULATE_FREEZE_OUT = True 

# Shell message
QUITE_MODE = False # True: no evolution progress message

# ---------------  TrENTo parameters  ----------------

""" The normalization is initiated in the ComBolt-ITA section."""


use_target_projectile_from_folder = True



Oxygen_1 = "/home/ktas/ge57vag/02-04-Aug-2025-ComBolt-ITA/O_and_Ne_samples/NLEFT_dmin_0.5fm_negativeweights_O.h5";
Oxygen_2 = "/home/ktas/ge57vag/02-04-Aug-2025-ComBolt-ITA/O_and_Ne_samples/NLEFT_dmin_0.5fm_positiveweights_O.h5";
Oxygen_3 = "/home/ktas/ge57vag/02-04-Aug-2025-ComBolt-ITA/O_and_Ne_samples/PGCM_clustered_dmin0_O.h5";
Oxygen_4 = "/home/ktas/ge57vag/02-04-Aug-2025-ComBolt-ITA/O_and_Ne_samples/PGCM_uniform_dmin0_O.h5";

Neon_1   = "/home/ktas/ge57vag/02-04-Aug-2025-ComBolt-ITA/O_and_Ne_samples/NLEFT_dmin_0.5fm_negativeweights_Ne.h5";

target     = Oxygen_2
projectile = Oxygen_2
trento_norm = 70
reduced_thickness_function = 0.0
fluctuation = 1.92
nucleon_minimum_distance = 0.0   # [fm]
nucleon_width = 0.6      # [fm]
constituent_width = 0.46   # [fm]
constituent_number = 4    # For constituent number = 1 the constituent width should be set equal to the nucleon width
cross_section = 6.8      # [fm^2]
grid_max = 10.0    # [fm]
grid_step = 0.2    # [fm]

minimum_bias = True
minimum_impact_parameter = 0.0    # [fm]
maximum_impact_parameter = 0.0  # [fm]

# -------------  Save the initial state info in ROOT -------

save_init_info_ROOT_format = True

# -------------  Particlization parameters  ----------------

particlization = True
ymax = 2.5
swiching_temp = TEMP_FO
overSampling = 1000 
save_particles_path = "/scratch8/ge57vag/15-14Aug2025-ComBolt-pp/evolution/ "
save_frzout_particles_plain_format = False
save_frzout_particles_ROOT_format = True

# ---------------  UrQMD parameters  ----------------
afterburner_stage = True
afterburner_path = "/home/ktas/ge57vag/02-04-Aug-2025-ComBolt-ITA/urqmd-afterburner/"
save_urqmd_particles_plain_format = False
save_urqmd_particles_ROOT_format = True

# --------------- generating many events ------------------ 

number_of_events = 1
run_IS_automatically = True   # True: run IS automatically, False: use pre-generated IS
