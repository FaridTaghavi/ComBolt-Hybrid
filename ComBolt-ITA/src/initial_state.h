// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.


// initial_state.h
#ifndef initial_state_H
#define initial_state_H

#include <omp.h>

#include "collision_kernel.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

#include "Administration.h"
#include "grid.h"
#include "parameters_class.h"
#include "Tmunu_and_hyd.h"

#include <boost/math/special_functions/jacobi_theta.hpp>
// #include "interpolation.h"

inline double sech(double x) {
    return 1.0 / std::cosh(x);
}

void F0_smooth( const parameters_class::initial_values &init,
                const double &x, const double &y, const double &vp, const double &vphi, double &result)
{
    double tyn{init.momentum_isotropicity};
    double R0 = init.R0;
    double eps2{init.ellipticity};
    double epsilon0 = init.energy_tot_tau0 / init.tau0;

    // result = (std::exp((std::pow(x, 2.) / (-1 + eps2) - std::pow(y, 2.) / (1 + eps2)) / (2. * std::pow(R0, 2.))) * epsilon0 * (1 / std::cosh(vp / tyn))) / (4. * std::sqrt(1 - std::pow(eps2, 2.)) * M_PI * std::pow(R0, 2.) * tyn * std::atan(std::sinh(1 / tyn)));
    result = 2.0 * (std::exp((std::pow(x, 2) / (-1.0 + eps2) - std::pow(y, 2) / (1.0 + eps2)) / 
                (2.0 * std::pow(R0, 2))) * epsilon0 * sech(vp / tyn)) / 
               (4.0 * std::sqrt(1.0 - std::pow(eps2, 2)) * M_PI * std::pow(R0, 2) * 
                tyn * std::atan(std::sinh(1.0 / tyn)));
};

void parton_smearing_function(const double &x,  const double &y,  const double &vz,  const double &vphi, 
                              const double &xi, const double &yi, const double &vzi, const double &vphi_i,
                              const double &sigma_v, const double &lambda, const double &sigma_phi, 
                              double &result)
{
    double smearX   {std::exp((-std::pow(x - xi,2) - std::pow(y - yi,2))/(2.*std::pow(sigma_v,2)))}; 
    double smearVz  { 1  /  (  (lambda*(std::atan(std::sinh((1 - vzi)/lambda)) + std::atan(std::sinh((1 + vzi)/lambda)))) * std::cosh((vz - vzi)/lambda))};
    double smearPhi {(1 / (2 * M_PI) ) * boost::math::jacobi_theta3(  (vphi-vphi_i) / 2.  , std::exp( - sigma_phi * sigma_phi / 2 )  )};
    
    result = smearX * smearVz * smearPhi; 
};

struct discretize_F
{
    void operator()(const parameters_class::initial_values &init, 
                    const lattice &lattice_, 
                    Field_grid &fg)
    {
        std::cout << "--> Discretizing the initial state\n";
        // const double dim { lattice_.axis.x.size()  };
        for (double x_point : lattice_.axis.x)
            for (double y_point : lattice_.axis.y)
                for (double vp_point : lattice_.axis.vp)
                    for (double vphi_point : lattice_.axis.vphi)
                    {
                        field fl;
                        fl.tau = init.tau0;
                        fl.coordinate.x = x_point;
                        fl.coordinate.y = y_point;
                        fl.coordinate.vp = vp_point;
                        fl.coordinate.vphi = vphi_point;
                        F0_smooth(init, x_point, y_point, vp_point, vphi_point, fl.value);
                        fg.push_back(fl);
                    };
        std::cout << "--> End of discretization.\n";
    }
};

struct construct_BD_from_Trento
{
    double tyn;
    double N; 

    void operator()(const parameters_class::initial_values &init, 
                    const parameters_class::trento_info &trento_inf, 
                    const lattice &lattice_, 
                    Field_grid &initial_BD)
    {
        tyn = init.momentum_isotropicity;
        N = trento_inf.norm / init.tau0;
        std::cout << "--> Reading initial state from TrENTo output ...\n";

        std::ifstream file(trento_inf.save_folder);         

        if (!file.is_open())
        {
            std::cerr << "--> Error opening file: " << trento_inf.save_folder << std::endl;
        }

        std::string line;
        double F0_xy;
        size_t ix{0}, iy{0};
        const size_t dim = lattice_.axis.x.size() * lattice_.axis.y.size() * lattice_.axis.vp.size() * lattice_.axis.vphi.size();
        initial_BD.resize(dim);
        while (std::getline(file, line))
        {
            ix = 0;
            std::istringstream iss(line);
            while (iss >> F0_xy) 
            { 
                for (size_t ivp = 0; ivp < lattice_.axis.vp.size(); ivp++)
                {
                    double vp_point = lattice_.axis.vp[ivp];
                    for (size_t ivphi = 0; ivphi < lattice_.axis.vphi.size(); ivphi++)
                    {
                        double vphi_point = lattice_.axis.vphi[ivphi];
                        double F0_smooth_momentum = (1 / std::cosh(vp_point / tyn)) / (2. * tyn * std::atan(std::sinh(1 / tyn)));
                        size_t ind;
                        lattice_.grid_index(ix, iy, ivp, ivphi, ind);
                        field fl;
                        fl.tau = init.tau0;
                        fl.coordinate.x = lattice_.axis.x[ix];
                        fl.coordinate.y = lattice_.axis.y[iy];
                        fl.coordinate.vp = vp_point;
                        fl.coordinate.vphi = vphi_point;
                        fl.value = 2.0 * N * F0_xy * F0_smooth_momentum;

                        initial_BD[ind] = fl;
                    };
                };
                ix++;
            };
            iy++;
        };
        std::cout << "--> End of reading initial state.\n";
    }
};

struct construct_BD_from_AMPT
{

    // We construct the initial state from AMPT. 
    // Since the AMPT initial state contains parton information in 3D, we need to project the 
    // longitudinal space into the transverse space by demanding partons are inside |delta_eta| < eta_cut. 
    // We also assume the initial time is the average time of initiated partons.
    // We also assume an overall normalization to fix the probable deviation in the final multiplicity.
    double hbarC = 0.19733; // GeV fm
    struct PartonData {
        double px, py, pz;
        double mass;
        double x, y, z, t;
    };
    std::vector<PartonData> partons;

    void operator()(const parameters_class::AMPT_info &AMPT_inf,
                    const lattice &lattice_, 
                    Field_grid &fg)
    {
        std::cout << "--> Reading AMPT file ...\n";

        std::ifstream file(AMPT_inf.save_folder);         

        if (!file.is_open())
        {
            std::cerr << "--> Error opening file: " << AMPT_inf.save_folder << std::endl;
        }
        
        std::string line; 
        // Skip the first line
        std::getline(file, line);

        double average_t = 0;
        double acceped_partons = 0;
        
        // Here, we read the AMPT file and store the parton information
        // We also calculate the average time of partons
        // We also filter the partons based on the eta cut and z > t condition
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::vector<double> entries;
            double value;

            // Read all 12 entries
            while (iss >> value) {
                entries.push_back(value);
            }

            if (entries.size() >= 9) {

                double etas = std::atanh(entries[7] / entries[8]);

                if (std::abs(etas) <= AMPT_inf.etas_cut)  // We keep partons in the midrapidity range
                {
                    if (entries[8] >= entries[7]) // Sometimes, AMPT makes partons at z > t, we ignore them
                    {
                        acceped_partons++;
                        average_t += entries[8];

                        PartonData pdata = {
                            entries[1] / hbarC, // px
                            entries[2] / hbarC, // py
                            entries[3] / hbarC, // pz
                            entries[4] / hbarC, // mass
                            entries[5], // x
                            entries[6], // y
                            entries[7], // z
                            entries[8]  // t
                        };
                        
                        partons.push_back(pdata);

                    }

                }
               
            } else {
                std::cerr << "--> Warning: Skipping malformed line.\n";
            }
        }
        std::cout << "--> Accepted partons: " << acceped_partons << "\n";
        if (acceped_partons > 0)
        {
            average_t /= acceped_partons;
            std::cout << "--> Average of partons' creation time: " << average_t << " [fm/c] \n";
        }
        else
        {
            std::cerr << "--> No partons accepted.\n";
            return;
        }
        std::cout << "--> End of reading AMPT file.\n";

        std::cout << "--> Constructing discretized initial state from AMPT ...\n";
        const size_t totalDim { lattice_.axis.x.size() * lattice_.axis.y.size() * lattice_.axis.vp.size() * lattice_.axis.vphi.size() };
        fg.resize(totalDim);


        #pragma omp parallel for collapse(2) num_threads(13)
        for (std::size_t i = 0; i < lattice_.axis.x.size(); ++i)
            for (std::size_t j = 0; j < lattice_.axis.y.size(); ++j)
                for (std::size_t k = 0; k < lattice_.axis.vp.size(); ++k)
                    for (std::size_t l = 0; l < lattice_.axis.vphi.size(); ++l) 
        // for (double x_point : lattice_.axis.x)
        //     for (double y_point : lattice_.axis.y)
        //         for (double vp_point : lattice_.axis.vp)
        //             for (double vphi_point : lattice_.axis.vphi)
                    {
                        size_t index;
                        lattice_.grid_index(i, j, k, l, index);
                        // std::cout << "i,j,k,l: " << i << " " << j << " " << k << " " << l << "    ->  "  << index <<  "\n";
                        field fl;
                        fl.tau = average_t;
                        double x_point = lattice_.axis.x[i];
                        double y_point = lattice_.axis.y[j];
                        double vp_point = lattice_.axis.vp[k];
                        double vphi_point = lattice_.axis.vphi[l];

                        fl.coordinate.x = x_point;
                        fl.coordinate.y = y_point;
                        fl.coordinate.vp = vp_point;
                        fl.coordinate.vphi = vphi_point;

                        double Fi;
                        fl.value = 0;
                        for (auto part : partons)
                        {
                            double tau_i = part.t;
                            double x_i = part.x;
                            double y_i = part.y;
                            double E_i = std::sqrt(part.px * part.px + part.py * part.py + part.pz * part.pz + part.mass * part.mass);
                            double vp_i = part.pz / E_i;
                            double vphi_i = std::atan2(part.py, part.px);

                            parton_smearing_function(x_point, y_point, vp_point,        vphi_point,
                                                     x_i,     y_i,     vp_i,            vphi_i,
                                                     AMPT_inf.sigma_v, AMPT_inf.lambda, AMPT_inf.sigma_phi, 
                                                    Fi);

                            fl.value += AMPT_inf.ampt_norm * E_i * Fi / tau_i;
                        }
                        // fg.push_back(  fl);
                        fg[index] = fl;
                        // std::cout << "fg " << index << " : " << fg[index].value << "\n";
                        
                    };
        std::cout << "--> End of construction.\n";

    }
};
struct dissipation
{
    void operator()(const parameters_class::dissipation_info &dissipation_inf,
                    const parameters_class::trento_info &trento_inf,
                    const parameters_class::medium_info &medium_inf,
                    const parameters_class::initial_values &init,
                    double &gamma)
    {
        if (dissipation_inf.dissipation_mode == 0)
        {
            gamma = 1 / (std::pow(medium_inf.C0, 0.25) * 5 * dissipation_inf.eta_over_s) ;
        }
        else if (dissipation_inf.dissipation_mode == 1 && init.initialization == 0)
        {
            gamma = dissipation_inf.gamma_hat * std::pow( init.R0 * init.energy_tot_tau0, -0.25 ) / M_PI;
        }
        else if (dissipation_inf.dissipation_mode == 1 && init.initialization == 1)
        {
            double R0;
            std::ifstream file(trento_inf.save_folder);          // Open the file
            if (!file.is_open())
            {
                std::cerr << "--> Error opening file: " << trento_inf.save_folder << std::endl;
            }

            std::string line;
            double F0_xy, x, y, e_tot, R02;
            double xmax = trento_inf.grid_max;
            double dx = trento_inf.grid_step;
            size_t ix{0}, iy{0};
            while (std::getline(file, line))
            {
                ix = 0;
                std::istringstream iss(line);
                while (iss >> F0_xy) 
                { 
                    x = -xmax + dx / 2 + ix * dx;
                    y = -xmax + dx / 2 + iy * dx;
                    e_tot += F0_xy;
                    R02 += F0_xy * (x * x + y * y);
                    ix++;
                }
            iy++;
            }
            R0 = std::sqrt(R02 / e_tot);
            gamma = dissipation_inf.gamma_hat * ( R0 * init.energy_tot_tau0 ) / M_PI;
        }
        else
        {
            std::cerr << "--> Dissipation mode is not defined.\n";
        }
    }
};

#endif // initial_state