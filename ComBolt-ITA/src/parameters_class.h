// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// parameters.h
#ifndef parameters_class_H
#define parameters_class_H

#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <string>
#include <vector>

class parameters_class
{
private:
    
        
public:

    static constexpr double hbarC = 0.19733; // GeV fm
    struct saving_info
    {
        std::string save_folder, surface_folder;
        bool save_BD;
        bool save_hydro;
        bool save_Tmunu;
        bool save_observables; 

        saving_info( 
            const std::string save_folder_ = "/home/farid/MyRepositories/KineticTheory/evolution/",
            const std::string surface_folder_ = "/home/farid/MyRepositories/KineticTheory/surface/",
            const bool save_BD_ = false,
            bool save_hydro_ = true,
            bool save_Tmunu_ = false,
            bool save_observables_ = true
          ) :
          save_folder(save_folder_),
          surface_folder(surface_folder_),
          save_BD(save_BD_),
          save_hydro(save_hydro_),
          save_Tmunu(save_Tmunu_),
          save_observables(save_observables_)
          {};
    };
    struct lattice_info
    {
        const size_t  n_x;
        const double  x_max;
        const size_t  n_v_p;
        const size_t  n_v_phi;
        const bool    Z2_symmetry;
    
        const size_t  n_y;
        const double  x_min;
        const double y_min, y_max;
        const double v_phi_min, v_phi_max;

        lattice_info( 
            const size_t n_x_ = 101,
            const double x_max_ = 19.5,
            const size_t n_v_p_ = 15,
            const size_t n_v_phi_ = 40,
            const bool Z2_symmetry_ = true
          ) :
          n_x(n_x_),
          x_max(x_max_),
          n_v_p(n_v_p_),
          n_v_phi(n_v_phi_),
          Z2_symmetry(Z2_symmetry_),
          n_y(n_x_),
          x_min(-x_max_),
          y_min(-x_max_),
          y_max(x_max_),
          v_phi_min(0.),
          v_phi_max(2. * M_PI - 2 * M_PI / n_v_phi_) 
          {};

    };
    struct numerical_calc_param
    {
        const double h;
        const double tau_max;
        const double alpha;

        numerical_calc_param(
            const double h_ = 0.0025,
            const double tau_max_ = 15.5,
            const double alpha_ = 0.02
          ) :
          h(h_),
          tau_max(tau_max_),
          alpha(alpha_)
          {};
    };
    struct initial_values 
    {
        const double tau0;
        const int initialization; 
        const double momentum_isotropicity;

        // For Gaussian input
        const double energy_tot_tau0;
        const double R0;
        const double ellipticity;

        initial_values(
            const double tau0_ = 5e-5,
            const int initialization_ = 2, // 0. Gaussian, 1. From Trento, 2. From AMPT
            const double momentum_isotropicity_ = 0.05,
            const double energy_tot_tau0_ = 1280 / hbarC,
            const double R0_ = 1.97,
            const double ellipticity_ = 0.416
          ) :
          tau0(tau0_),
          initialization(initialization_),
          momentum_isotropicity(momentum_isotropicity_),
          energy_tot_tau0(energy_tot_tau0_),
          R0(R0_),
          ellipticity(ellipticity_)
          {};

    };
    struct medium_info
    {
        const double C0; //Confomtal equation of state: energy = C0 * T^4
        const double epsFO; 

        medium_info(
            const double C0_ = 13.8997,
            const double epsFO_ = 0.30 / hbarC // 1/fm^4
          ) :
          C0(C0_),
          epsFO(epsFO_)
          {};

    };
    struct dissipation_info
    {
        const unsigned int dissipation_mode; 
        const double eta_over_s; 
        const double gamma_hat;

        dissipation_info(
            const unsigned int dissipation_mode_ = 0,// 0. eta / s, 1. opacity   
            const double eta_over_s_ = 1.0 / (4. * M_PI),
            const double gamma_hat_ = 200.0
          ) :
          dissipation_mode(dissipation_mode_),
          eta_over_s(eta_over_s_),
          gamma_hat(gamma_hat_)
          {};

    };
    struct trento_info
    {
        const double grid_max;
        const double grid_step;
        const double norm;
        std::string save_folder;

        trento_info(
            const double grid_max_ = 15.0,
            const double grid_step_ = 0.2,
            const double norm_ = 15.0 / hbarC, // 1/fm   
            // NOTE: We assumed Trento normalization has dimension [fm].
            // We assumed Reduced Thickness Function T_R = tau_0 * energy_density, whihc has dimension GeV/fm^2.
            // Therefore, we consider the GeV to fm conversion in the above normalization to find energy_density = T_R / tau_0 / hbarC, which has dimensionl 1/fm^4.
            const std::string save_folder_ = "/home/farid/MyRepositories/KineticTheory/trento-master/events/event_set6/8.dat"
          ) :
          grid_max(grid_max_),
          grid_step(grid_step_),
          norm(norm_),
          save_folder(save_folder_)
          {};
    };
    struct AMPT_info
    {
        const double sigma_v, lambda, sigma_phi;
        const double ampt_norm;
        const double etas_cut;
        std::string save_folder;
        AMPT_info(
            const double sigma_v_ = 0.05,
            const double lambda_ = 0.5,
            const double sigma_phi_ = 0.5,
            const double ampt_norm_ = 1.0,
            const double etas_cut_ = 0.5,
            const std::string save_folder_ = "/home/farid/MyRepositories/KineticTheory/kineticTheory/initial_state_AMPT/parton-initial-afterPropagation.dat"
          ) :
          sigma_v(sigma_v_),
          lambda(lambda_),
          sigma_phi(sigma_phi_),
          ampt_norm(ampt_norm_),
          etas_cut(etas_cut_),
          save_folder(save_folder_)
          {};
    };
    struct interpolation_parameters
    {
        const double insert_threshold; 
        const unsigned int number_of_divistion;
        const double remove_threshold; 

        interpolation_parameters(
            const double insert_threshold_ = 0.2, // For eta/s = 1.0/4pi , 2.0/4pi  0.2 is good
            const unsigned int number_of_divistion_ = 2,
            const double remove_threshold_ = 5e-7 // For eta/s = 1.0/4pi , 2.0/4pi  5e-7 is good
          ) :
          insert_threshold(insert_threshold_),
          number_of_divistion(number_of_divistion_),
          remove_threshold(remove_threshold_)
          {};
    };
    struct multicore_parameters
    {
        const unsigned int number_of_cores_BoltzEq_solver;
        const unsigned int number_of_cores_integration;
        const unsigned int number_of_cores_interpolation;
        const unsigned int number_of_cores_calc_obs;

        multicore_parameters(
            const unsigned int number_of_cores_BoltzEq_solver_ = 13,
            const unsigned int number_of_cores_integration_ = 13,
            const unsigned int number_of_cores_interpolation_ = 13,
            const unsigned int number_of_cores_calc_obs_ = 13
          ) :
          number_of_cores_BoltzEq_solver(number_of_cores_BoltzEq_solver_),
          number_of_cores_integration(number_of_cores_integration_),
          number_of_cores_interpolation(number_of_cores_interpolation_),
          number_of_cores_calc_obs(number_of_cores_calc_obs_)
          {};
    };

    struct freeze_out_info
    {
        const bool calculate_freeze_out = true;
        const bool end_evol_last_frozen_cell = true;

        freeze_out_info(
            const bool calculate_freeze_out_ = true,
            const bool end_evol_last_frozen_cell_ = true
          ) :
          
            calculate_freeze_out(calculate_freeze_out_),
            end_evol_last_frozen_cell(end_evol_last_frozen_cell_)
          {};
    };
    
    struct shell_message
    {
        const bool quiet = false;
        shell_message(
            const bool quiet_ = false
          ) :
          quiet(quiet_)
          {};
    };

    saving_info saving_info_;
    lattice_info lattice_info_;
    numerical_calc_param numerical_calc_param_;
    initial_values initial_values_;
    medium_info medium_info_;
    dissipation_info dissipation_info_;
    trento_info trento_info_;
    AMPT_info ampt_info_;
    interpolation_parameters interpolation_parameters_;
    multicore_parameters multicore_parameters_;
    freeze_out_info freeze_out_info_;
    shell_message shell_message_;

    parameters_class(
        const saving_info& saving_info_,
        const lattice_info& lattice_info_,
        const numerical_calc_param& numerical_calc_param_,
        const initial_values& initial_values_,
        const medium_info& medium_info_,
        const dissipation_info& dissipation_info_,
        const trento_info& trento_info_,
        const AMPT_info& ampt_info_,
        const interpolation_parameters& interpolation_parameters_,
        const multicore_parameters& multicore_parameters_,
        const freeze_out_info& freeze_out_info_,
        const shell_message& shell_message_ 
    ) 
        : saving_info_(saving_info_), 
          lattice_info_(lattice_info_), 
          numerical_calc_param_(numerical_calc_param_), 
          initial_values_(initial_values_), 
          medium_info_(medium_info_), 
          dissipation_info_(dissipation_info_), 
          trento_info_(trento_info_), 
          ampt_info_(ampt_info_), 
          interpolation_parameters_(interpolation_parameters_), 
          multicore_parameters_(multicore_parameters_),
          freeze_out_info_(freeze_out_info_),
          shell_message_(shell_message_)
    {}

};


#endif //parameters
