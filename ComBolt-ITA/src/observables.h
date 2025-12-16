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
                     const parameters_class::dissipation_info &disp_inf,
                     const Tmunu_grid &discret_Tmunu, const Hydro_grid &discret_hydro, observables &obs )
    {
        int num_threads = ncore.number_of_cores_calc_obs;

        // First find the center of the initial state
        double x_bar {0.0};
        double y_bar {0.0};
        std::vector<double> eps_tot_local(num_threads, 0.0);
        std::vector<double> x_bar_local(num_threads, 0.0);
        std::vector<double> y_bar_local(num_threads, 0.0);
        

        #pragma omp parallel num_threads(num_threads)
        {
            int tid = omp_get_thread_num();
            #pragma omp for 
            for (size_t ind = 0; ind < discret_Tmunu.size(); ++ind) 
            {
                eps_tot_local[tid]  += discret_hydro[ind].ed;
                x_bar_local[tid] += discret_hydro[ind].x * discret_hydro[ind].ed;
                y_bar_local[tid] += discret_hydro[ind].y * discret_hydro[ind].ed;
            }
        }

        
        // Merge thread-local results outside the loop
        for (int i = 0; i < num_threads; i++) {
            obs.eps_tot    += eps_tot_local[i];
            x_bar      += x_bar_local[i];
            y_bar      += y_bar_local[i];
        }

        x_bar /= obs.eps_tot;
        y_bar /= obs.eps_tot;

        // Calculate the initial ecc and momentim, ...
        
        obs.tau = discret_Tmunu[0].tau;
        // double eps_total {0.};
        double av_x2 {0.};
        double av_y2 {0.};
        double av_xy {0.};
        double av_r2 {0.0};

        double av_x3 {0.};
        double av_y3 {0.};
        double av_x2y {0.};
        double av_xy2 {0.};
        double av_r3 {0.0};


        double av_gamma {0.};
        size_t n_total_nonzero_Teff {0};
        
        std::vector<double> Etr_local(num_threads, 0.0);
        std::vector<double> eps_p_x_local(num_threads, 0.0);
        std::vector<double> eps_p_y_local(num_threads, 0.0);
        std::vector<double> u_perp_local(num_threads, 0.0);
        std::vector<double> av_invReyn_local(num_threads, 0.0);

        std::vector<double> av_x2_local(num_threads, 0.0);
        std::vector<double> av_y2_local(num_threads, 0.0);
        std::vector<double> av_xy_local(num_threads, 0.0);
        std::vector<double> av_r2_local(num_threads, 0.0);
        
        std::vector<double> av_x3_local(num_threads, 0.0);
        std::vector<double> av_y3_local(num_threads, 0.0);
        std::vector<double> av_x2y_local(num_threads, 0.0);
        std::vector<double> av_xy2_local(num_threads, 0.0);
        std::vector<double> av_r3_local(num_threads, 0.0);
        
        std::vector<double> av_gamma_local(num_threads, 0.0);
        std::vector<size_t> n_nonzero_Teff(num_threads, 0);
        
        #pragma omp parallel num_threads(num_threads)
        {
            int tid = omp_get_thread_num();
            #pragma omp for 
            for (size_t ind = 0; ind < discret_Tmunu.size(); ++ind) 
            {

                double x = discret_hydro[ind].x - x_bar;
                double y = discret_hydro[ind].y - y_bar;
                

                Etr_local[tid]      += (discret_Tmunu[ind].Txx + discret_Tmunu[ind].Tyy);
                eps_p_x_local[tid]  += (discret_Tmunu[ind].Txx - discret_Tmunu[ind].Tyy);
                eps_p_y_local[tid]  += (2 * discret_Tmunu[ind].Txy);
                u_perp_local[tid]   += (std::sqrt(discret_hydro[ind].ux * discret_hydro[ind].ux + discret_hydro[ind].uy * discret_hydro[ind].uy)) * discret_hydro[ind].ed;
                av_invReyn_local[tid] += (discret_hydro[ind].inverseReyn) * discret_hydro[ind].ed;

                
                av_x2_local[tid] += x * x * discret_hydro[ind].ed;
                av_y2_local[tid] += y * y * discret_hydro[ind].ed;
                av_xy_local[tid] += x * y * discret_hydro[ind].ed;
                av_r2_local[tid] += (x * x + y * y) * discret_hydro[ind].ed;

                av_x3_local[tid]  += x * x * x * discret_hydro[ind].ed;
                av_y3_local[tid]  += y * y * y * discret_hydro[ind].ed;
                av_x2y_local[tid] += x * x * y * discret_hydro[ind].ed;
                av_xy2_local[tid] += x * y * y * discret_hydro[ind].ed;
                av_r3_local[tid] += std::pow( x * x + y * y, 1.5 ) * discret_hydro[ind].ed;

                // if (discret_hydro[ind].T_eff != 0)
                //     av_gamma_local[tid] +=  ( 1.0 / (5.0 * disp_inf.eta_over_s)  ) *  discret_hydro[ind].T_eff / std::pow( discret_hydro[ind].ed, 0.25);
                // if (discret_hydro[ind].T_eff != 0)
                // {
                //     av_gamma_local[tid] +=  ( 1.0 / (5.0 * disp_inf.eta_over_s)  ) *  discret_hydro[ind].T_eff / std::pow( discret_hydro[ind].ed, 0.25);
                //     n_nonzero_Teff[tid] ++;
                // }  

                if (discret_hydro[ind].T_eff != 0)
                {
                    double local_eta_over_s;
                    eta_over_s_temperature_dependent(discret_hydro[ind].T_eff, local_eta_over_s, disp_inf.eta_over_s_min, disp_inf.eta_over_s_slope, disp_inf.eta_over_s_pow, disp_inf.eta_over_s_Tc);
                    av_gamma_local[tid] +=  ( 1.0 / (5.0 * local_eta_over_s)  ) *  discret_hydro[ind].T_eff / std::pow( discret_hydro[ind].ed, 0.25);
                    n_nonzero_Teff[tid] ++;
                } 

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
            av_x2      += av_x2_local[i];
            av_y2      += av_y2_local[i];
            av_xy      += av_xy_local[i];
            av_r2      += av_r2_local[i];

            av_x3      += av_x3_local[i];
            av_y3      += av_y3_local[i];
            av_x2y     += av_x2y_local[i];
            av_xy2     += av_xy2_local[i];
            av_r3      += av_r3_local[i];

            av_gamma   += av_gamma_local[i];
            n_total_nonzero_Teff += n_nonzero_Teff[i];
        }

        av_gamma /= n_total_nonzero_Teff;


        obs.eps_2_x = -(av_x2 - av_y2 ) / av_r2;
        obs.eps_2_y = -2 * av_xy / av_r2;

        obs.eps_3_x = -( av_x3 - 3 * av_xy2 ) / av_r3;
        obs.eps_3_y = -( 3 * av_x2y - av_y3 ) / av_r3;

        obs.Rsq = av_r2 / obs.eps_tot;
        
        // av_x2 *= Delta_x * Delta_y;
        // av_y2 *= Delta_x * Delta_y;
        // av_xy *= Delta_x * Delta_y;

        // av_x3  *= Delta_x * Delta_y;
        // av_y3  *= Delta_x * Delta_y;
        // av_x2y *= Delta_x * Delta_y;
        // av_xy2 *= Delta_x * Delta_y;


        // obs.Rsq = av_x2 + av_y2;

        // obs.eps_2_x = -(av_x2 - av_y2 ) / obs.Rsq;
        // obs.eps_2_y = -2 * av_xy / obs.Rsq;

        // obs.eps_3_x = -( av_x3 - 3 * av_xy2 ) / std::pow(obs.Rsq, 1.5);
        // obs.eps_3_y = -( 3 * av_x2y - av_y3 ) / std::pow(obs.Rsq, 1.5);
        
        obs.eps_p_x     /= obs.Etr;
        obs.eps_p_y     /= obs.Etr; 
        obs.u_perp      /= obs.eps_tot;
        obs.av_invReyn  /= obs.eps_tot;
        
        obs.eps_tot *= Delta_x * Delta_y;
        obs.Etr     *= Delta_x * Delta_y;

        // obs.Rsq /= obs.eps_tot;

        obs.gamma_hat = av_gamma * std::pow(obs.Rsq, 0.125 ) *  std::pow(obs.eps_tot * obs.tau / M_PI, 0.25);
    }

 
};

#endif // observables.h
