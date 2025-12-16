// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.


// collision_kernel.h
#ifndef collision_kernel_H
#define collision_kernel_H

// #include "interpolation.h"
#include "lattice_EOS.h"

inline void collision_kernel(const field &BD, const hydrodynamic_field &hyd, const double &gamma, double &kernel)
{
    double vp {BD.coordinate.vp};
    double vphi {BD.coordinate.vphi};
    
    double v_perp {(std::abs(vp) > 1 - 1.e-16) ? 0.0 : std::sqrt(1 - vp * vp)};
    
    double vx {v_perp * std::cos(vphi )};
    double vy {v_perp * std::sin(vphi )};

    double utau {hyd.utau   }; 
    double ed {hyd.ed};
    double ux {hyd.ux}; 
    double uy {hyd.uy}; 
    double ueta {hyd.ueta   }; 
    
    double VmuUmu {1. * utau - vx * ux - vy * uy - BD.tau  * vp * ueta};  
    
    kernel =   -4 * vp * vp / BD.tau * BD.value  +   gamma *  ( - VmuUmu * std::pow(ed, 0.25) * BD.value +  std::pow(ed, 1.25) / std::pow(VmuUmu, 3.)  );

};


inline void eta_over_s_temperature_dependent(
    const double &T_eff,
    double &eta_over_s,
    const double min,
    const double Slope,
    const double Pow,
    const double Tc
){
    // eta/s = min + Slope * (T - Tc) * (T/Tc)^Pow
    eta_over_s = min + Slope * (T_eff - Tc) * std::pow(T_eff / Tc, Pow);
    // std::cout << "eta / s    " << eta_over_s << "    " << min  << "    " << Slope  << "    " << Pow  << "    " << Tc  << "\n";
}

inline void pseudo_lattice_gamma(const double &EtaOverS, const double &ed, const LatticeEOSInverse &eos_inv, double &gamma)
{
    double T_eff; 
    // double ed_min = 1e-3 / 0.197; // 1/fm^4. The freezeout energy density should be in the order of 0.3 GeV/fm^3.
    
    eos_inv.get_T_from_ed(ed, T_eff);
    double T_min = 54.43 / 197.3; // T = 54.43 MeV  is the minimum temperature in the lattice interpolation function that gives positive energy density, which is well smaller that 
                                    // switching temperature 155 MeV. 
    if (T_eff <= T_min) {
        gamma = 0.0;
    } else {
        double tau_relax = 5 * EtaOverS / T_eff; 
        gamma = 1.0 / ( tau_relax * std::pow(ed, 0.25) );
    }
    // std::cout << "T_eff = " << T_eff << "\n"; 
};
#endif // collision_kernel
