# ComBolt-Hybrid

A modular hybrid heavy-ion collision event generator combining initial-state modeling, kinetic pre-equilibrium evolution, and hadronic afterburner stages within a unified framework.

This project implements the **CoMBolt-ITA hybrid approach**, based on the relativistic Boltzmann equation in the Isotropization Time Approximation (ITA), and connects it to realistic initial condition generators and hadronic cascade simulations. The framework is designed to study collective dynamics in small and large collision systems directly from kinetic theory.

---

## 🚀 Overview

ComBolt-Hybrid provides an end-to-end simulation pipeline for relativistic nuclear collisions:

1. Initial condition generation  
2. Pre-equilibrium kinetic evolution (CoMBolt-ITA)  
3. Particlization / freeze-out  
4. Hadronic afterburner evolution  

The framework bridges microscopic transport theory and macroscopic collective flow observables in a controlled and modular way.

---

## 🧠 Physics Background

### CoMBolt-ITA

The CoMBolt-ITA model evolves the phase-space distribution function using the relativistic Boltzmann equation in the Isotropization Time Approximation (ITA). This approach:

- Interpolates naturally between free streaming and hydrodynamic behavior  
- Dynamically generates collective flow  
- Allows controlled variation of effective opacity  
- Connects kinetic theory to viscous hydrodynamics in the appropriate limit  

The model has been developed and applied in:

1. S. F. Taghavi and S. M. A. Tabatabaee Mehr,  
   *CoMBolt-ITA, a Collective Model via relativistic Boltzmann equation in Isotropization Time Approximation*,  
   arXiv:2504.17707  

2. S. F. Taghavi and S. M. A. Tabatabaee Mehr,  
   *Opacity estimation of OO collision from CoMBolt-ITA hybrid*,  
   arXiv:2512.05009  

These works demonstrate that CoMBolt-ITA reproduces viscous hydrodynamic behavior in the fluid regime while providing controlled deviations at larger viscosities or anisotropies. It has been applied to small systems and OO collisions at the LHC to estimate effective medium opacity.

---

## 📦 Repository Structure

ComBolt-Hybrid/
├── ComBolt-ITA/              # Core kinetic evolution module  
├── trento-master/            # Initial condition generator (TrENTo)  
├── urqmd-afterburner/        # Hadronic cascade (UrQMD)  
├── evolution/                # Evolution drivers and configuration  
├── frzout/                   # Freeze-out and particlization tools  
├── Ampt-v1.26t9b-v2.26t9b/    # AMPT snapshots (reference / comparison)  
└── .vscode/                  # Editor configuration files  

Main components:

- TrENTo – Parametric initial condition generator  
- CoMBolt-ITA – Kinetic theory evolution  
- Freeze-out tools (frzout) – Conversion to particle distributions  
- UrQMD afterburner – Hadronic rescattering stage  

---

## 🛠️ Requirements

- C++ compiler with C++17 support  
- CMake (≥ 3.10 recommended)  
- Fortran compiler (for UrQMD)  
- Python (optional, for analysis scripts)  

---

## 🔧 Build Instructions

git clone https://github.com/FaridTaghavi/ComBolt-Hybrid.git  
cd ComBolt-Hybrid  

mkdir build  
cd build  

cmake ..  
make -j$(nproc)  

Adjust compiler flags or toolchain settings if needed for your system.

---

## ▶️ Typical Workflow

1. Generate initial conditions

   ./trento [options]

2. Run CoMBolt-ITA evolution

   ./ComBoltITA --config config.yaml

3. Particlization / freeze-out

   Convert kinetic output into particle distributions using the provided freeze-out utilities.

4. Hadronic afterburner

   afterburner input.dat output.dat

5. Analysis

   Post-process final particle data to compute flow harmonics (v_n), multiplicity distributions, transverse momentum spectra, and correlation observables.

---

## 🔬 Scientific Scope

This framework is suited for:

- Small-system collectivity studies (pp, pPb, OO)  
- Opacity and transport coefficient estimation  
- Pre-equilibrium dynamics  
- Comparison between kinetic theory and hydrodynamics  
- Hybrid modeling across different collision systems  

---

## 📚 Citation

If you use this framework, please cite:

@article{Taghavi:2025ComBoltITA,  
  author  = {Taghavi, S. F. and Tabatabaee Mehr, S. M. A.},  
  title   = {CoMBolt-ITA, a Collective Model via relativistic Boltzmann equation in Isotropization Time Approximation},  
  eprint  = {2504.17707},  
  archivePrefix = {arXiv},  
  primaryClass  = {nucl-th}  
}

@article{Taghavi:2025OpacityOO,  
  author  = {Taghavi, S. F. and Tabatabaee Mehr, S. M. A.},  
  title   = {Opacity estimation of OO collision from CoMBolt-ITA hybrid},  
  eprint  = {2512.05009},  
  archivePrefix = {arXiv},  
  primaryClass  = {nucl-th}  
}

---

## 🤝 Contributing

Contributions are welcome. Possible improvements include:

- Extended documentation and tutorials  
- Validation benchmarks  
- Additional collision systems  
- Performance optimization  
- Extended transport kernels or equation-of-state options  

Please open an issue or submit a pull request.

---

## 📄 License

Specify the project license here (e.g., MIT, GPL, etc.).

---

## Acknowledgments

This project builds upon widely used community tools such as TrENTo and UrQMD and contributes to ongoing efforts in kinetic theory modeling of relativistic nuclear collisions.
