// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// parameters.h
#ifndef parameters_H
#define parameters_H

#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <string>
#include <vector>

const double hbarC = 0.19733; // GeV fm

struct final_parameters
{
    std::string save_folder = "/home/farid/MyRepositories/KineticTheory/evolution/";
    bool save_BD          = false;
    bool save_hydro       = true;
    bool save_Tmunu       = false;
    bool save_observables = true; 
};
struct lattice_info
{
    const size_t  n_x     {101 };
    const double  x_max   {19.5};
    const size_t  n_v_p   {15};
    const size_t  n_v_phi {15};
    const bool    Z2_symmetry = true;
   
    const double  x_min { -x_max };
    const size_t  n_y {n_x };
    const double y_min  {-x_max}, y_max  {x_max};
    const double v_phi_min {0.} , v_phi_max {2. * M_PI - 2 * M_PI / n_v_phi};

};
struct Numerical_calc_param
{
    const double h {0.005};
    const double tau_max {15.5};
};
struct initial_values 
{
    const double tau0 {5e-5};
    const int initialization {2}; // 0. Gaussian, 1. From folder initial state, 2. From Trento

    const double momentum_isotropicity {0.05};
    
    // For Gaussian input
    const double energy_tot_tau0 {1280 / hbarC};
    const double R0 {1.97};
    const double ellipticity {0.416};

};
struct medium_info
{
    const double C0 {13.8997}; //Confomtal equation of state: energy = C0 * T^4
    const double epsFO {0.30 / hbarC}; // 1/fm^4

};
struct dissipation_info
{
    const unsigned int dissipation_mode {0}; // 0. eta / s, 1. opacity   
    
    const double eta_over_s {1.0 / (4. * M_PI)}; 
    const double gamma_hat {200.0};

};
struct trento_info
{
    const double grid_max  {15.0};
    const double grid_step {0.2};
    const double norm      {15.0 / hbarC}; // 1/fm
    // const double energy_tot_tau0 {1 / hbarC}; // 1/fm
    std::string save_folder = "/home/farid/MyRepositories/KineticTheory/trento-master/events/event_set6/8.dat";
    // std::string save_folder = "/home/farid/MyRepositories/KineticTheory/trento-master/events/Gaussian/0.dat";
};
struct interpolation_parameters
{
    const double insert_threshold = 0.15; // For eta/s = 1.0/4pi , 2.0/4pi  0.2 is good
    const unsigned int number_of_divistion = 2;
    const double remove_threshold = 5e-7; // For eta/s = 1.0/4pi , 2.0/4pi  5e-7 is good
};
struct multicore_parameters
{
    const unsigned int number_of_cores_BoltzEq_solver = 13;
    const unsigned int number_of_cores_integration = 13;
    const unsigned int number_of_cores_interpolation = 13;
    const unsigned int number_of_cores_calc_obs = 13;
};

#endif //parameters
