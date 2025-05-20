// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.


// collision_kernel.h
#ifndef collision_kernel_H
#define collision_kernel_H

// #include "interpolation.h"

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
#endif // collision_kernel