// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// BoltzEQ_solver.h

#ifndef BoltzEQ_solver_H
#define BoltzEQ_solver_H

#include <omp.h>
#include <filesystem>
#include "Tmunu_and_hyd.h"
#include "observables.h"
#include "interpolation.h"
#include "freeze-out.h"
#include "cornelius.h"

inline void coordinates_step_forward(const double &tau0, const double &tau, const coord &coord0, coord &coord)
{
    double norm_shift;
    if ( std::abs(coord0.vp) < 1 - 1e-16)
    {
        norm_shift = ((-tau0 + tau * std::sqrt(1 - (1 - (tau0 * tau0) / (tau * tau)) * (coord0.vp * coord0.vp)))    ) / std::sqrt(1 - (coord0.vp * coord0.vp));
    } else
    {
        norm_shift = 0;
    }

    coord.x = coord0.x + norm_shift * std::cos(coord0.vphi); 
    coord.y = coord0.y + norm_shift * std::sin(coord0.vphi);
    coord.vp = (tau0 * coord0.vp) / (tau * std::sqrt(1 - (1 - (tau0 * tau0) / (tau * tau)) * (coord0.vp * coord0.vp)));
    coord.vphi = coord0.vphi;

};

struct euler_step_forward
{
    void operator()(const parameters_class::multicore_parameters &ncore,
                    const parameters_class::interpolation_parameters &interp_param, 
                    const double h, const double &gamma, lattice &lattice_I, const Field_grid &BD_I, const Hydro_grid &hyd_grd, lattice &lattice_O, Field_grid &BD_O)
    {

        const size_t totalDim { lattice_I.axis.x.size() * lattice_I.axis.y.size() * lattice_I.axis.vp.size() * lattice_I.axis.vphi.size() };
        Field_grid BD_interm(totalDim);

        #pragma omp parallel for num_threads(ncore.number_of_cores_BoltzEq_solver)
        for (size_t point_ind = 0; point_ind < totalDim; ++point_ind)
        {
            size_t hyd_ind;
            lattice_I.project_to_poistion_grid_index(point_ind, hyd_ind);
            
            hydrodynamic_field hyd;
            hyd = hyd_grd[hyd_ind];
            
            double CF;
            collision_kernel(BD_I[point_ind], hyd, gamma, CF);
        
            BD_interm[point_ind].tau   = BD_I[point_ind].tau   + h;
            BD_interm[point_ind].value = BD_I[point_ind].value + h * CF;
            coordinates_step_forward(BD_I[point_ind].tau, BD_interm[point_ind].tau, BD_I[point_ind].coordinate, BD_interm[point_ind].coordinate);

        }
        BoltzDist_interpolation boltzInterp;
        boltzInterp(ncore, interp_param, lattice_I, BD_interm, lattice_O, BD_O);
        
    }
};

struct solve_BoltzmannEQ
{

    axes_grid axis;
    Field_grid discret_BD_O;
    Tmunu_grid discrete_Tmunu;
    Hydro_grid discrete_hyd;
    Administration admin;
    euler_step_forward euler_step_frwrd;
    dissipation dissipation_;
    parameters_class params;
    size_t i_tau {0};
    double h; 
    double gamma; 

    solve_BoltzmannEQ(const parameters_class &params_) :
        params(params_)
    {}
    
    void operator()()
    {
        // Initiate the lattice
        lattice lattice_I(axis), lattice_O(axis);
        // lattice_info l_inf;
        lattice_I.initiate_lattice(params.lattice_info_, params.initial_values_, params.trento_info_, axis);
        
        Field_grid discret_BD_I;
        switch (params.initial_values_.initialization) {
            case 0: {
                // Setup the discretized initial state
                discretize_F discF;
                discF(params.initial_values_, lattice_I, discret_BD_I);
                break;
            }
            case 1 :{
                construct_BD_from_Trento trento_IS;
                trento_IS(params.initial_values_, params.trento_info_, lattice_I, discret_BD_I);
                break;
            }
            case 2 :{
                construct_BD_from_AMPT discF_AMPT;
                discF_AMPT(params.ampt_info_, lattice_I, discret_BD_I);
                break;
            }
            default:
                std::cerr << "--> The initial value is not defined." << std::endl;
}
        std::string folder_path = params.saving_info_.save_folder;
        std::string folder_surface = params.saving_info_.surface_folder;
        
        std::string file_OB = "observables.dat";
        if (std::filesystem::exists(folder_path + file_OB)) {
            std::filesystem::remove(folder_path + file_OB);
        }
        
        if (!std::filesystem::exists(folder_surface)) {
            std::filesystem::create_directory(folder_surface); // or fs::create_directories(folder) for nested paths
        std::cout << "--> Created directory: " << folder_surface << std::endl;
        }
        // else {
        //     std::cout << "--> Directory already exists: " << folder_surface << std::endl;
        // }

        std::string file_surface = "surface.dat";
        if (std::filesystem::exists(folder_surface + file_surface)) {
            std::filesystem::remove(folder_surface + file_surface);
        }

        
        dissipation_(params.dissipation_info_,
                     params.trento_info_,
                     params.medium_info_,
                     params.initial_values_,
                     gamma);

        
        //=========================
        // double epsilon_FO = 0.3 / 0.197; // fm^-4
        Hydro_grid discret_hyd_prev_step; // Keep the previous step of hydrodynamic fields for the cornelius code
        bool first_step = true;         
        bool found_freeze_out = true;    
        bool continue_loop = true;
        double last_tau = 0.0;
        //=========================
        
        std::cout << "--> Starting the Boltzmann equation solver ...\n";
        // while (discret_BD_I[0].tau < params.numerical_calc_param_.tau_max)
        while (continue_loop)
        {  
            // Check if the evolution stops after all parts are frozen out or continue until tau_max
            last_tau = discret_BD_I[0].tau;
            if (params.freeze_out_info_.end_evol_last_frozen_cell)
            {
                if (params.freeze_out_info_.calculate_freeze_out)
                {
                    continue_loop = found_freeze_out;
                } else {
                    std::cout << "   * No freeze out calculation is demanded, continue until tau_max.\n";
                    continue_loop = (discret_BD_I[0].tau < params.numerical_calc_param_.tau_max);
                }
            
            } else {

                continue_loop = (discret_BD_I[0].tau < params.numerical_calc_param_.tau_max);
            }
            
            // double tau_0 = params.initial_values_.tau0;
            double h_inf = params.numerical_calc_param_.h;
            double alpha = params.numerical_calc_param_.alpha;
			h = h_inf * ( ( discret_BD_I[0].tau )/ (  discret_BD_I[0].tau  +  h_inf / alpha) );
            
            // if (!params.shell_message_.quiet)
            //     std::cout << "\n ======== step number " << i_tau << " at tau = " << discret_BD_I[0].tau << " ========\n\n";
            if (!params.shell_message_.quiet)
                std::cout << "   * Step number " << i_tau << " at tau = " << discret_BD_I[0].tau << "\n";

            Tmunu_and_hyd::Tmunu_and_hydro_from_BoltzmannDist tmunu_hyd(lattice_I);
            tmunu_hyd(  params.lattice_info_, 
                        params.multicore_parameters_,
                        discret_BD_I, discrete_Tmunu, discrete_hyd);

            //=================== Freeze out: cornelius on action ===================
            if (params.freeze_out_info_.calculate_freeze_out)
            {
                if (first_step)
                {
                    discret_hyd_prev_step = discrete_hyd;
                    first_step = false;

                } else 
                {
                    cornelius_on_action cornelius;
                    cornelius(params, h, lattice_I, discrete_hyd, discret_hyd_prev_step, folder_surface + file_surface, found_freeze_out );
                    discret_hyd_prev_step = discrete_hyd;
                }
            }
            //=======================================================================
            
            save_fields_new(params.saving_info_,
                            discret_BD_I, 
                            discrete_hyd, 
                            discrete_Tmunu, 
                            folder_path,  
                            i_tau);

            if (params.saving_info_.save_observables)
            {
                observables obs;
                calculate_observables calc_obs(lattice_I);
                calc_obs(params.multicore_parameters_, discrete_Tmunu, discrete_hyd, obs);
                admin.save_obs(obs, folder_path + file_OB);
            }
   

            // Timer timer;
            // timer.start();
            euler_step_frwrd(params.multicore_parameters_, 
                             params.interpolation_parameters_,  
                             h, 
                             gamma, 
                             lattice_I, 
                             discret_BD_I, 
                             discrete_hyd, 
                             lattice_O, 
                             discret_BD_O);
            // timer.stop();
            
            discret_BD_I = discret_BD_O;
            lattice_I.axis = lattice_O.axis;

            i_tau++;
        }
        std::cout << "--> End of the Boltzmann equation solver at tau = " << last_tau << " [fm/c].\n";
    }
};

#endif // BoltzEQ_solver_H
