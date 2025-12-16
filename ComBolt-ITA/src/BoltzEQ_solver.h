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

#include "lattice_EOS.h"

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
                    const double h, const parameters_class::dissipation_info &dissip_inf, 
                    // const latticeEOS_Table &EOS_table,
                    // const double &gamma,
                    lattice &lattice_I, const Field_grid &BD_I, const Hydro_grid &hyd_grd, lattice &lattice_O, Field_grid &BD_O)
    {

        const size_t totalDim { lattice_I.axis.x.size() * lattice_I.axis.y.size() * lattice_I.axis.vp.size() * lattice_I.axis.vphi.size() };
        Field_grid BD_interm(totalDim);
        double gamma;
        double local_eta_over_s;

        // LatticeEOSInverse eos_inv(EOS_table);
        #pragma omp parallel for num_threads(ncore.number_of_cores_BoltzEq_solver)
        for (size_t point_ind = 0; point_ind < totalDim; ++point_ind)
        {
            size_t hyd_ind;
            lattice_I.project_to_poistion_grid_index(point_ind, hyd_ind);
            
            hydrodynamic_field hyd;
            hyd = hyd_grd[hyd_ind];


            if (hyd.T_eff != 0.0){
                eta_over_s_temperature_dependent(hyd.T_eff, local_eta_over_s, dissip_inf.eta_over_s_min, dissip_inf.eta_over_s_slope, dissip_inf.eta_over_s_pow, dissip_inf.eta_over_s_Tc);
                gamma = ( 1.0 / (5.0 * local_eta_over_s)  )  *  ( hyd.T_eff / std::pow(hyd.ed, 0.25) );
            } else {
                gamma = 0.0;
            } 

            // if (hyd.T_eff != 0.0)
            //     gamma = ( 1.0 / (5.0 * dissip_inf.eta_over_s)  )  *  ( hyd.T_eff / std::pow(hyd.ed, 0.25) );
            // else
            //     gamma = 0.0;  
            // // Determine gamma_internal based on disspation_mode
            // // 0: eta/s constant, 1: opacity constant, 2: eta/s constant with tau_relax(ed) from lattice.
            // if (  dissip_inf.dissipation_mode == 0 || dissip_inf.dissipation_mode == 1) {
            //     gamma_internal = gamma;
            // } else if (dissip_inf.dissipation_mode == 2) {
            //     // std::cout << "Before pseudo_lattice_gamma: " << gamma_internal << "\n";
            //     pseudo_lattice_gamma(dissip_inf.eta_over_s, hyd.ed, eos_inv, gamma_internal);
            //     // std::cout << "After pseudo_lattice_gamma: " << gamma_internal << "\n";
            // } else {
            //     std::cerr << "The dissipation_mode is not defined.\n";
            //     throw std::runtime_error("The dissipation_mode is not defined.");
            // }

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
    // dissipation dissipation_;
    parameters_class params;
    size_t i_tau {0};
    double h; 
    double gamma; 


    // Massive
    axes_grid_massive axis_massive;

    solve_BoltzmannEQ(const parameters_class &params_) :
        params(params_)
    {}
    
    void operator()()
    {
        // initialize the axes grid massive
        lattice_massive lattice_massive_I;
        lattice_massive_I.initiate_lattice(params.lattice_info_massive_, params.initial_values_, params.trento_info_);
        size_t ind;
        lattice_massive_I.grid_index(15, 10, 3, 8, 12, ind);
        std::cout << "Check grid function: " << ind << "\n";
        std::array<size_t, 5> ijklm ;
        lattice_massive_I.inverse_grid_index(ind, ijklm);
        std::cout << "Check inverse grid function: " << ijklm[0] << " " << ijklm[1] << " " << ijklm[2] << " " << ijklm[3] << " " << ijklm[4] << "\n";


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
        
        std::string file_init_inf = "init_inf.dat";
        if (std::filesystem::exists(folder_surface + file_init_inf)) {
            std::filesystem::remove(folder_surface + file_init_inf);
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

        
        // dissipation_(params.dissipation_info_,
        //              params.trento_info_,
        //              params.medium_info_,
        //              params.initial_values_,
        //              gamma);

        
        //=========================
        // double epsilon_FO = 0.3 / 0.197; // fm^-4
        Hydro_grid discret_hyd_prev_step; // Keep the previous step of hydrodynamic fields for the cornelius code
        bool first_step = true;         
        bool first_step_calculate_initial = true;         
        bool found_freeze_out = true;    
        bool continue_loop = true;
        double last_tau = 0.0;
        double save_interval = 0.2; // Every 0.5 fm/c calculate and save the observables
        double next_save_time = 0.0;
        double min_hadrton_form_time = 0.0; // fm
        //=========================

           
        // Tablulate the EOS for later use

        // Find the maximum value of energy density and corresponding temperature
        max_value_ed max_ed;
        double ed_max;
        max_ed(discret_BD_I, ed_max);

        double C0_comformal = 1.0; // For High temp QCD C0_comformal is around 13.9 but we take small value get the high temperature limit.
        double T_max = std::pow(ed_max / C0_comformal, 0.25); // Initial guess for T_max
        std::cout << "--> Maximum energy density at initial time is ed_max = " << ed_max << " [GeV/fm^3].\n";
        std::cout << "--> Corresponding maximum temperature is T_max = " << T_max << " [GeV].\n";
        double T_step = 0.1 / 197.3; // 0.1 MeV in GeV

        latticeEOS_Table EOS_table;
        make_lattice_table make_table;
        make_table(T_max, T_step, EOS_table);

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
            // h = 0.025;      
            // if (!params.shell_message_.quiet)
            //     std::cout << "\n ======== step number " << i_tau << " at tau = " << discret_BD_I[0].tau << " ========\n\n";
            if (!params.shell_message_.quiet)
                std::cout << "   * Step number " << i_tau << " at tau = " << discret_BD_I[0].tau << "\n";

            Tmunu_and_hyd::Tmunu_and_hydro_from_BoltzmannDist tmunu_hyd(lattice_I);
            // tmunu_hyd(  params.lattice_info_, 
            //             params.multicore_parameters_,
            //             discret_BD_I, discrete_Tmunu, discrete_hyd);
            tmunu_hyd(  params,
                        EOS_table,
                        discret_BD_I, discrete_Tmunu, discrete_hyd);
            //=================== Freeze out: cornelius on action ===================
                if (params.freeze_out_info_.calculate_freeze_out)
                {
                    if (discret_BD_I[0].tau > min_hadrton_form_time)
                    {
                        if (first_step)
                        {
                            discret_hyd_prev_step = discrete_hyd;
                            first_step = false;

                        } else 
                        {
                            {
                                cornelius_on_action cornelius;
                                cornelius(params, h, lattice_I, discrete_hyd, discret_hyd_prev_step, folder_surface + file_surface, found_freeze_out );
                                discret_hyd_prev_step = discrete_hyd;
                            }
                        }
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
                calc_obs(params.multicore_parameters_, 
                         params.dissipation_info_,
                         discrete_Tmunu, discrete_hyd, obs);
                admin.save_obs(obs, folder_path + file_OB);
            }
   
            // // I just commented out the previous lines saving only the initial state later properly
            // if (params.saving_info_.save_init_inf && first_step_calculate_initial )
            // // if (params.saving_info_.save_observables )
            // {
            //     first_step_calculate_initial = false;
            //     observables obs;
            //     calculate_observables calc_obs(lattice_I);
            //     calc_obs(params.multicore_parameters_, 
            //              params.dissipation_info_,
            //              discrete_Tmunu, discrete_hyd, obs);
            //     // admin.save_obs(obs, folder_path + file_OB);
            //     admin.save_obs(obs, folder_surface + file_init_inf);
            // }
   
            if (params.saving_info_.save_init_inf && discret_BD_I[0].tau >= next_save_time - 0.001 )
            // if (params.saving_info_.save_observables )
            {
                first_step_calculate_initial = false;
                observables obs;
                calculate_observables calc_obs(lattice_I);
                calc_obs(params.multicore_parameters_, 
                         params.dissipation_info_,
                         discrete_Tmunu, discrete_hyd, obs);
                // admin.save_obs(obs, folder_path + file_OB);
                obs.life_time = last_tau;
                admin.save_obs(obs, folder_surface + file_init_inf);
                next_save_time += save_interval;
                // std::cout << "  save interval ------------------------->>>>  " << save_interval << "\n";
            }
            // Timer timer;
            // timer.start();
            euler_step_frwrd(params.multicore_parameters_, 
                             params.interpolation_parameters_,  
                             h, 
                             params.dissipation_info_,
                             // EOS_table,
                             // gamma, 
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
