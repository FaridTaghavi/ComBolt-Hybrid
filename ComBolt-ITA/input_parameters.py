# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
#
# This file is part of a DFG-funded research project (grant no. 517518417).
# See the top-level LICENSE file for license details.

# input_parameters.py

# =============  ComBolt-ITA parameters =================

# Spatial lattice setup (Gaussian mode and AMPT mode)
NX = 150
XMAX = 2.0  # [fm]

# Momentum lattice setup
NVP = 15
NVPHI = 15
Z2_SYMMETRY = True # The Boltzmann distribution is symmetric for pozitive and negative vz.

# Numerical method parameters
H = 0.002       # time step [fm]
ALPHA_H = 0.003 # time steps stiffness
TAU_MAX = .12  # [fm] maximum time

# Initial conditions
TAU0 = 5.0e-4        # [fm] initial time 
INIT_MODE = 1          # Initial condition mode: 0. Gaussian 1. External energy density (from TrENTo), 2. external Boltz Dist From AMPT
MOM_ISO = 1000.0         # Initial momentum anisotropy: MOM_ISO -> 0 highly anisotropic, MOM_ISO -> infinity close to isotropy
# Gaussian spatial distribution
ENERGY_TAU0 = 1280     # [GeV / fm] initial total transverse energy * initial time (Gaussian mode)
R0 = 1.97              # [fm] initial system size (Gaussian mode)
ELLIPTICITY = 0.416    # initial ellipticity (Gaussian mode)

# Trento settings
TRENTO_FILE = "../trento-master/events/event_set7/0.dat"
# GRID_MAX = 15.0     # TrENTo
# GRID_STEP = 0.2
NORM = 25.0         # Scaling the initial TrENTo energy density (we set the nomaliatin in TrENTo equal to 1)

# AMPT settings
SMEAR_SIGMA_V = 0.5
SMEAR_LAMBDA = 10.5
SMEAR_SIGMA_PHI = 10.5
AMPT_NORM = 0.5
ETAS_CUT = 2.0
AMPT_FILE = "initial_state_AMPT/0002.dat"

# Medium properties
# C0 = 13.8997     # Conformal equation of state: energy = C0 T^4
C0 = 4.6

TEMP_FO = 0.155                        # [GeV] Freeze-out temperature
EPSFO = C0 * TEMP_FO**4 / 0.197**3      # [GeV/fm^3] Freeze-out energy density 

# EPSFO = 0.3                           # [GeV/fm^3] Freeze-out energy density 
# TEMP_FO = (EPSFO * 0.197**3 / C0)**(1/4) # [GeV] Freeze-out temperature


# Dissipation
DISS_MODE = 0
GAMMA_HAT = 20.0
ETA_S = 0.010

# Interpolation
INSERT_THRESHOLD = 0.15    # Dividie segments in the vz directon if it is greater thatn INSERT_THRESHOLD
NUM_DIVISION = 2           # Number of divisions
REMOVE_THRESHOLD = 5e-7    # Remove the vz node if its distane to vz = 0 is REMOVE_THRESHOLD

# Saving options
OUTPUT_DIR = "../evolution/"
TEMP_PATH = "/home/farid/MyRepositories/KineticTheory/evolution/temp/"
SURFACE_DIR = TEMP_PATH
SAVE_BD          = False
SAVE_HYDRO       = False
SAVE_TMUNU       = False
SAVE_OBSERVABLES = False

# Multithreading
CORES_BOLTZ = 13
CORES_INT = 13
CORES_INTERP = 13
CORES_OBS = 13

# Freeze-out parameters

END_EVOL_LAST_FROZEN_CELL = True # True: end evolution of the last frozen cell, False: continue evolution until tau_max
CALCULATE_FREEZE_OUT = True 

# Shell message
QUITE_MODE = False # True: no evolution progress message

# ---------------  TrENTo parameters  ----------------

""" The normalization is initiated in the ComBolt-ITA section."""

target = "p"
projectile = "p"
reduced_thickness_function = 0.0
fluctuation = 10.
nucleon_minimum_distance = 0.0   # [fm]
nucleon_width = 1.0      # [fm]
constituent_width = 0.8   # [fm]
constituent_number = 5     # For constituent number = 1 the constituent width shoul be set equal to the nucleon width
cross_section = 7.32      # [fm^2]
grid_max = 4.0    # [fm]
grid_step = 0.1    # [fm]

minimum_bias = False
minimum_impact_parameter = 0.0    # [fm]
maximum_impact_parameter = 0.0   # [fm]

# -------------  Particlization parameters  ----------------

particlization = True
ymax = 2.5
swiching_temp = TEMP_FO
overSampling = 3000
save_particles_path = "/home/farid/MyRepositories/KineticTheory/evolution/"
save_frzout_particles_plain_format = False
save_frzout_particles_ROOT_format = True

# ---------------  UrQMD parameters  ----------------
afterburner_stage = True 
afterburner_path = "/home/farid/MyRepositories/KineticTheory/urqmd-afterburner/"
save_urqmd_particles_plain_format = False
save_urqmd_particles_ROOT_format = True

# --------------- generating many events ------------------ 

number_of_events = 100
run_IS_automatically = True   # True: run IS automatically, False: use pre-generated IS
