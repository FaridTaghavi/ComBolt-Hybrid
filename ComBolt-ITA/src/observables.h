// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// observables.h
#ifndef observables_H
#define observablea_H

#include "initial_state.h"


// One should be able to combine this loop into other loops
// For the moment, we keep it separate for readability.
struct calculate_observables
{
    // COMMENT! the distribution is not centered!!!!! move it to the center first!
    const lattice &lattice_;
    calculate_observables(const lattice &latt ) : lattice_(latt) {}

    double Delta_x {lattice_.axis.x[1] - lattice_.axis.x[0]};
    double Delta_y {lattice_.axis.y[1] - lattice_.axis.y[0]};
    
    void operator()( const parameters_class::multicore_parameters &ncore,
                    const Tmunu_grid &discret_Tmunu, const Hydro_grid &discret_hydro, observables &obs )
    {

        obs.tau = discret_Tmunu[0].tau;
        // double eps_total {0.};
        double av_x2 {0.};
        double av_y2 {0.};
        double av_xy {0.};
        
        int num_threads = ncore.number_of_cores_calc_obs;
        std::vector<double> Etr_local(num_threads, 0.0);
        std::vector<double> eps_p_x_local(num_threads, 0.0);
        std::vector<double> eps_p_y_local(num_threads, 0.0);
        std::vector<double> u_perp_local(num_threads, 0.0);
        std::vector<double> av_invReyn_local(num_threads, 0.0);
        std::vector<double> eps_tot_local(num_threads, 0.0);

        #pragma omp parallel num_threads(num_threads)
        {
            int tid = omp_get_thread_num();
            #pragma omp for 
            for (size_t ind = 0; ind < discret_Tmunu.size(); ++ind) 
            {
                Etr_local[tid]      += (discret_Tmunu[ind].Txx + discret_Tmunu[ind].Tyy);
                eps_p_x_local[tid]  += (discret_Tmunu[ind].Txx - discret_Tmunu[ind].Tyy);
                eps_p_y_local[tid]  += (2 * discret_Tmunu[ind].Txy);
                u_perp_local[tid]   += (std::sqrt(discret_hydro[ind].ux * discret_hydro[ind].ux + discret_hydro[ind].uy * discret_hydro[ind].uy)) * discret_hydro[ind].ed;
                av_invReyn_local[tid] += (discret_hydro[ind].inverseReyn) * discret_hydro[ind].ed;
                eps_tot_local[tid]  += discret_hydro[ind].ed;
            }
        }

        // Merge thread-local results outside the loop
        for (int i = 0; i < num_threads; i++) {
            obs.Etr        += Etr_local[i];
            obs.eps_p_x    += eps_p_x_local[i];
            obs.eps_p_y    += eps_p_y_local[i];
            obs.u_perp     += u_perp_local[i];
            obs.av_invReyn += av_invReyn_local[i];
            obs.eps_tot    += eps_tot_local[i];
        }


        av_x2 *= Delta_x * Delta_y;
        av_y2 *= Delta_x * Delta_y;
        av_xy *= Delta_x * Delta_y;

        obs.eps_2_x = (av_y2 - av_x2 ) / (av_y2 + av_x2);
        obs.eps_2_y = 2 * av_xy / (av_y2 + av_x2);
        obs.R = std::sqrt( av_x2 + av_y2 );
        
        obs.eps_p_x     /= obs.Etr;
        obs.eps_p_y     /= obs.Etr; 
        obs.u_perp      /= obs.eps_tot;
        obs.av_invReyn  /= obs.eps_tot;
        
        obs.eps_tot *= Delta_x * Delta_y;
        obs.Etr     *= Delta_x * Delta_y;
    }

 
};

#endif // observables.h