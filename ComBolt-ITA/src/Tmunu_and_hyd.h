// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// Tmunu_and_hyd.h
#ifndef Tmunu_and_hyd_H
#define Tmunu_and_hyd_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include "Administration.h"
#include "grid.h"
#include "initial_state.h"
#include "interpolation.h"
#include "error_handler.h"
#include <omp.h>

class Tmunu_and_hyd     
{
private:

public:
    
    Tmunu_and_hyd() {}

    struct integral_abscissas_weights
    {
        integral_abscissas_weights() {};

        std::vector<double> trapezodial_weights(const std::vector<double> &axis)
        {
            std::vector<double> weight;
            size_t n = axis.size();
            for (size_t ind = 0; ind < n; ++ind){
                
                if(ind == 0) {
                    weight.push_back( 0.5 * (axis[ind + 1] - axis[ind] )  );
                } else if(ind == n - 1) {
                    weight.push_back( 0.5 * (axis[ind]- axis[ind - 1] )  );
                } else {
                    weight.push_back( 0.5 * (axis[ind + 1] - axis[ind - 1] ) );
                };
                
            };
            return weight;

        };
            
    };
    
    struct Tmunu_integrand
    {
        void operator()(const parameters_class::lattice_info &lat_inf,
            const double &F, const double &vp, const double &vphi, const double &tau, energy_momentum_field &Tmunu_integrnd_O)
        {
            double Z2;
            lat_inf.Z2_symmetry ? Z2 = 2. : Z2 = 1.;
            
            Tmunu_integrnd_O.Ttautau  = 1 / (4 * M_PI) * F * Z2;
            Tmunu_integrnd_O.Ttaueta  = 1 / (4 * M_PI) *  vp  / tau * F  * (2 - Z2); //Zero because of  Z2 symm and odd terrm with respect to vp
            Tmunu_integrnd_O.Txx      = 1 / (4 * M_PI) * (1 - vp * vp) * std::cos( vphi) * std::cos( vphi) * F * Z2;
            Tmunu_integrnd_O.Txy      = 1 / (4 * M_PI) * (1 - vp * vp) * std::cos( vphi) * std::sin( vphi) * F * Z2;
            Tmunu_integrnd_O.Tetaeta  = 1 / (4 * M_PI) * (vp * vp) / (tau * tau) * F * Z2;
            Tmunu_integrnd_O.Tyy      = 1 / (4 * M_PI) * (1 - vp * vp ) * std::sin( vphi) * std::sin( vphi) * F * Z2;
            if (std::abs(vp) >= 1 - 1e-16) 
            {
                Tmunu_integrnd_O.Ttaux = 0;
                Tmunu_integrnd_O.Ttauy = 0;
                Tmunu_integrnd_O.Txeta = 0;
                Tmunu_integrnd_O.Tyeta = 0;
            } else {
                Tmunu_integrnd_O.Ttaux    = 1 / (4 * M_PI) * std::sqrt( 1 - vp * vp ) * std::cos( vphi) * F * Z2;
                Tmunu_integrnd_O.Ttauy    = 1 / (4 * M_PI) * std::sqrt( 1 - vp * vp ) * std::sin( vphi) * F * Z2;
                Tmunu_integrnd_O.Txeta    = 1 / (4 * M_PI) * (vp * std::sqrt( 1 - vp * vp) * std::cos( vphi)) / tau *  F * (2 - Z2); 
                Tmunu_integrnd_O.Tyeta    = 1 / (4 * M_PI) * (vp * std::sqrt(1 - vp * vp) * std::sin( vphi)) / tau * F * (2 - Z2);
            }

        };
    };

    template<typename T>
    struct Landau_matching
    {

        // First index up second index down.
        void operator()(const T &em, hydrodynamic_field &hyd)
        { 
            hyd.tau = em.tau;
            hyd.x = em.x;
            hyd.y = em.y;
        
            // First index up second index down.
            // In this case the Laudau matching condition  
            // T^{munu} u_mu =  lambda g^{munu}u_mu  leads to
            // |T^mu_nu - lambda delta^mu_nu| = 0, meaning finding eigen system of
            // T^mu_nu
            gsl_matrix* Tmu_nu = gsl_matrix_alloc(4, 4);
            gsl_matrix_set(Tmu_nu, 0, 0,  em.Ttautau);
            gsl_matrix_set(Tmu_nu, 0, 1, -em.Ttaux);
            gsl_matrix_set(Tmu_nu, 0, 2, -em.Ttauy);
            gsl_matrix_set(Tmu_nu, 0, 3, -em.Ttaueta * em.tau * em.tau);
            gsl_matrix_set(Tmu_nu, 1, 0,  em.Ttaux);
            gsl_matrix_set(Tmu_nu, 1, 1, -em.Txx);
            gsl_matrix_set(Tmu_nu, 1, 2, -em.Txy);
            gsl_matrix_set(Tmu_nu, 1, 3, -em.Txeta * em.tau * em.tau);
            gsl_matrix_set(Tmu_nu, 2, 0,  em.Ttauy);
            gsl_matrix_set(Tmu_nu, 2, 1, -em.Txy);
            gsl_matrix_set(Tmu_nu, 2, 2, -em.Tyy);
            gsl_matrix_set(Tmu_nu, 2, 3, -em.Tyeta * em.tau * em.tau);
            gsl_matrix_set(Tmu_nu, 3, 0,  em.Ttaueta);
            gsl_matrix_set(Tmu_nu, 3, 1, -em.Txeta);
            gsl_matrix_set(Tmu_nu, 3, 2, -em.Tyeta);
            gsl_matrix_set(Tmu_nu, 3, 3, -em.Tetaeta * em.tau * em.tau);
            
            // Allocate space for eigenvalues and eigenvectors
            gsl_vector_complex* eval = gsl_vector_complex_alloc(4);
            gsl_matrix_complex* evec = gsl_matrix_complex_alloc(4, 4);
            gsl_eigen_nonsymmv_workspace* workspace = gsl_eigen_nonsymmv_alloc(4);
                 // Perform the eigenvalue decomposition
            gsl_set_error_handler_off();
            // std::cout << "Tmunu matrix " << gsl_matrix_get( Tmu_nu,0 ,0) << "\n";
            int status = gsl_eigen_nonsymmv(Tmu_nu, eval, evec, workspace);
            // if (status != GSL_SUCCESS) {
                // std::cout << "The landau matching has been failed at tau = " << em.tau << " and x, y " << em.x << " , " << em.y << "\n"; 
            try{
                
                landau_matching_eigenvalues(status);
            } catch(const std::exception& e) {
                std::cerr << e.what() << '\n';
                // hyd.ed   = 0.;
                // hyd.utau = 1.;
                // hyd.ux   = 0.;
                // hyd.uy   = 0.;
                // hyd.ueta = 0.;
                exit(EXIT_FAILURE);
            }
                 // u is an upper index vector
            double utau {1.}, ux {0.}, uy {0.}, ueta {0.};
            double energy {-1.};
            for (size_t i = 0; i < 4; i++) 
            {
                gsl_complex eval_i = gsl_vector_complex_get(eval, i);
                // std::cout << " inside Landau mathcing  " << GSL_REAL(eval_i) <<  "   " << GSL_IMAG(eval_i) << "\n";
                if (std::abs(GSL_IMAG(eval_i))  < 1.e-10  && GSL_REAL(eval_i) >= 0 )
                {
                    // a is an upper index vector
                    double a[4]; 
                    for (size_t j = 0; j < 4; j++) {
                        gsl_complex evec_ij = gsl_matrix_complex_get(evec, j, i);
                        a[j]   = GSL_REAL(evec_ij);
                    }
                    double normSq = a[0] * a[0] - a[1] * a[1] - a[2] * a[2] - a[3] * a[3] * (em.tau * em.tau);
                    // std::cout << "normSq  " << normSq << "\n";
                    
                    // Check if the vector is time-like
                    if (normSq > 0 ) 
                    {
                        energy = GSL_REAL(eval_i);
                        
                        if (a[0] < 0) for (int q = 0; q < 4; q++) a[q] = -a[q]; // The first component of four velocity should be always positive. 
                        
                        utau = a[0] / std::sqrt(normSq);
                        ux   = a[1] / std::sqrt(normSq);
                        uy   = a[2] / std::sqrt(normSq);
                        ueta = a[3] / std::sqrt(normSq);
                        break; // If the energy is found, stop checking other eigenvalues
                    } 
                }
            }
            // try {
            //     landau_matching_energy(energy);
            // } catch(const std::exception& e) {
            //     std::cerr << e.what() << '\n';
            //     exit(EXIT_FAILURE);
            // }
            
            if (energy == -1. ) {
                energy = 0;
                utau = 1;
                ux = 0;
                uy = 0;
                ueta = 0;
                // throw std::runtime_error("Landau matching: reals and positive eignevalue could not be found.");
            }

            // Free the allocated memory
            gsl_vector_complex_free(eval);
            gsl_matrix_complex_free(evec);
            gsl_matrix_free(Tmu_nu);
            gsl_eigen_nonsymmv_free(workspace);
        
            hyd.ed   = energy;
            hyd.utau = utau;
            hyd.ux   = ux;
            hyd.uy   = uy;
            hyd.ueta = ueta;
            
            double Px, Py, Pz; //, eps_l;

            // eps_l = std::pow(hyd.tau,4)*em.Tetaeta*std::pow(ueta,2) + em.Ttautau*std::pow(utau,2) - 2*em.Ttaux*utau*ux + em.Txx*std::pow(ux,2) - 2*em.Ttauy*utau*uy + 2*em.Txy*ux*uy + em.Tyy*std::pow(uy,2) + 
            //     2*std::pow(hyd.tau,2)*ueta*(-(em.Ttaueta*utau) + em.Txeta*ux + em.Tyeta*uy);
            // eps_l = em.Ttautau*std::pow(utau,2) - 2*em.Ttaux*utau*ux + em.Txx*std::pow(ux,2) - 2*em.Ttauy*utau*uy + 2*em.Txy*ux*uy + em.Tyy*std::pow(uy,2);
            Px = (em.Txx*std::pow(1 + utau + std::pow(ux,2),2) + ux*(std::pow(hyd.tau,4)*em.Tetaeta*std::pow(ueta,2)*ux + em.Ttautau*std::pow(1 + utau,2)*ux - 2*em.Ttaux*(1 + utau)*(1 + utau + std::pow(ux,2)) + 
                2*(-(em.Ttauy*(1 + utau)*ux) + em.Txy*(1 + utau + std::pow(ux,2)))*uy + em.Tyy*ux*std::pow(uy,2) + 
                2*std::pow(hyd.tau,2)*ueta*(-(em.Ttaueta*(1 + utau)*ux) + em.Txeta*(1 + utau + std::pow(ux,2)) + em.Tyeta*ux*uy)))/std::pow(1 + utau,2);

            Py = (em.Tyy*std::pow(1 + utau + std::pow(uy,2),2) + uy*(2*em.Txy*(1 + utau)*ux + std::pow(hyd.tau,4)*em.Tetaeta*std::pow(ueta,2)*uy + 
                (em.Ttautau*std::pow(1 + utau,2) + ux*(-2*em.Ttaux*(1 + utau) + em.Txx*ux))*uy + 2*em.Txy*ux*std::pow(uy,2) - 2*em.Ttauy*(1 + utau)*(1 + utau + std::pow(uy,2)) + 
                2*std::pow(hyd.tau,2)*ueta*(-((em.Ttaueta + em.Ttaueta*utau - em.Txeta*ux)*uy) + em.Tyeta*(1 + utau + std::pow(uy,2)))))/std::pow(1 + utau,2); 
            
            Pz = (std::pow(hyd.tau,2)*(em.Tetaeta*std::pow(1 + std::pow(hyd.tau,2)*std::pow(ueta,2) + utau,2) + 
               ueta*(em.Ttautau*ueta*std::pow(1 + utau,2) - 2*em.Ttaueta*(1 + utau)*(1 + std::pow(hyd.tau,2)*std::pow(ueta,2) + utau) + 
                  ux*(-2*em.Ttaux*ueta*(1 + utau) + 2*em.Txeta*(1 + std::pow(hyd.tau,2)*std::pow(ueta,2) + utau) + em.Txx*ueta*ux) + 
                  2*(-(em.Ttauy*ueta*(1 + utau)) + em.Tyeta*(1 + std::pow(hyd.tau,2)*std::pow(ueta,2) + utau) + em.Txy*ueta*ux)*uy + em.Tyy*ueta*std::pow(uy,2))))/std::pow(1 + utau,2);
            
            // ====================
            // hyd.BulkTensorXX = eps_l;
            hyd.BulkTensorXX = 0.0;

            hyd.Px = Px;
            hyd.Py = Py;
            hyd.Pz = Pz;

            hyd.pi_tautau = em.Ttautau  - energy / 3 * ( 4 * utau * utau - 1);
            hyd.pi_taux   = em.Ttaux    - energy / 3 * ( 4 * utau * ux      );
            hyd.pi_tauy   = em.Ttauy    - energy / 3 * ( 4 * utau * uy      );
            hyd.pi_taueta = em.Ttaueta  - energy / 3 * ( 4 * utau * ueta    );
            hyd.pi_xx     = em.Txx      - energy / 3 * ( 4 * ux   * ux   + 1);
            hyd.pi_xy     = em.Txy      - energy / 3 * ( 4 * ux   * uy      );
            hyd.pi_xeta   = em.Txeta    - energy / 3 * ( 4 * ux   * ueta    );
            hyd.pi_yy     = em.Tyy      - energy / 3 * ( 4 * uy   * uy   + 1);
            hyd.pi_yeta   = em.Tyeta    - energy / 3 * ( 4 * uy   * ueta    );
            hyd.pi_etaeta = em.Tetaeta  - energy / 3 * ( 4 * ueta * ueta + 1 / ( hyd.tau * hyd.tau) );

            double pimunu_pimunu {   std::pow(hyd.pi_tautau,2) - 2*std::pow(hyd.pi_taux,2) - 2*std::pow(hyd.pi_tauy,2) + std::pow(hyd.pi_xx,2) + 2*std::pow(hyd.pi_xy,2) + std::pow(hyd.pi_yy,2) 
                              - 2*std::pow(hyd.pi_taueta,2)*std::pow(hyd.tau,2) +  2*std::pow(hyd.pi_xeta,2)*std::pow(hyd.tau,2) + 2*std::pow(hyd.pi_yeta,2)*std::pow(hyd.tau,2)
                              + std::pow(hyd.pi_etaeta,2)*std::pow(hyd.tau,4)};
            if ( pimunu_pimunu < 0)
            {
                pimunu_pimunu = -pimunu_pimunu;
            }               
            if (energy != 0)
            {
                hyd.inverseReyn = std::sqrt( 6 * pimunu_pimunu   /   (  energy * energy  )  )  ;
            } else {

                hyd.inverseReyn = 0;
            }
            if (std::isnan(hyd.inverseReyn)) {
                std::cout << "The value of inverseReyn is nan at point: tau " << hyd.tau << " x " << hyd.x << " y " << hyd.y << " energy " << energy << " pimunupimunu " << pimunu_pimunu << "\n";
            }

            hyd.isotropicity = ( Px + Py - 2 * Pz ) / energy;
        
        };
    };

    struct Tmunu_from_BoltzmannDist
    {
        const lattice &lattice_;
        Tmunu_integrand Tmunu_integ;
        energy_momentum_field Tmunu_result;
        integral_abscissas_weights integ_abs_weight;
        energy_momentum_field Tmunu_integrnd;
        
        Tmunu_from_BoltzmannDist(const lattice &latt ) : lattice_(latt) {}
        
        void operator()(const parameters_class::lattice_info &l_inf,
                        const parameters_class::multicore_parameters &ncore,
                        const Field_grid &BD,  const size_t &ix, const size_t &iy, energy_momentum_field &Tmunu_O)
        {
            // lattice_info l_inf;
            // initial_axis_grid iniAxGrid;
            // initial_state init_state;
            // iniAxGrid = init_state.initiate_axis(l_inf);


            Tmunu_result.tau = BD[0].tau;
            Tmunu_result.x   = lattice_.axis.x[ix];
            Tmunu_result.y   = lattice_.axis.y[iy]; 

            size_t index;
            std::vector<double> vp_grid_at_this_point;
            vp_grid_at_this_point.reserve(lattice_.axis.vp.size());
            for (size_t i_vp = 0; i_vp < lattice_.axis.vp.size(); ++i_vp)
            {
                lattice_.grid_index(ix, iy, i_vp, 0, index);
                vp_grid_at_this_point.push_back (BD[index].coordinate.vp);
            }  
            
            std::vector<double> vp_weight = integ_abs_weight.trapezodial_weights(   vp_grid_at_this_point );
            
            // Important! The grid in the vphi direction does not change! Unlike vp 
            std::vector<double> vphi_axis_2pi =  lattice_.axis.vphi;
            vphi_axis_2pi.push_back(2 * M_PI); //Add 2*pi to the end of vphi_grid
            std::vector<double> vphi_weight = integ_abs_weight.trapezodial_weights(  vphi_axis_2pi );
            
            /* Loop over vp */
            size_t vp_abs_size { lattice_.axis.vp.size()};
            
            // #pragma omp parallel for collapse(2) num_threads(ncore.number_of_cores_integration)
            for (size_t i_vp = 0.; i_vp < vp_abs_size; ++i_vp)
            {
                energy_momentum_field Tmunu_over_vphi;
                double vp   { vp_grid_at_this_point[i_vp] };
                double w_vp { vp_weight[i_vp] };
                double Fi_times_wi;
                
                /* Loop over vphi ring */
                size_t ind;
                size_t vphi_abs_size { vphi_axis_2pi.size() };
                for (size_t i_vphi = 0. ; i_vphi  < vphi_abs_size; ++i_vphi)
                {
                    size_t ii_vphi = i_vphi;
                    if (i_vphi == vphi_abs_size -1) ii_vphi = 0; // For vphi = 2*pi we set F at vphi = 0
                    
                    double vphi   { vphi_axis_2pi[i_vphi]};
                    double w_vphi { vphi_weight[i_vphi] };

                    double F;
                    lattice_.grid_index(ix, iy, i_vp, ii_vphi, ind);
                    F = BD[ind].value;
                    Fi_times_wi = F * w_vphi;
                    // energy_momentum_field Tmunu_integrnd;
                    Tmunu_integ(l_inf, Fi_times_wi, vp, vphi, BD[0].tau, Tmunu_integrnd);
                    
                    Tmunu_over_vphi += Tmunu_integrnd;

                }

                Tmunu_result += Tmunu_over_vphi * w_vp; //The scalar should be on the right-hand side
            }
            Tmunu_O = Tmunu_result;
        };
    };

    struct Tmunu_and_hydro_from_BoltzmannDist
    {
        const lattice &lattice_;
        Tmunu_and_hydro_from_BoltzmannDist(const lattice &latt ) : lattice_(latt) {}
        
        void operator()(const parameters_class::lattice_info &l_inf,
                        const parameters_class::multicore_parameters &ncore,
                        const Field_grid &BD, Tmunu_grid &discrete_Tmunu, Hydro_grid &discrete_hyd)
        {
            // iniAxGrid = init_state.initiate_axis(l_inf);
            size_t position_dim = lattice_.axis.x.size() * lattice_.axis.y.size();
            
            Tmunu_grid Tmunu_Grd(position_dim);
            Hydro_grid Hydro_Grd(position_dim);

            #pragma omp parallel for num_threads(ncore.number_of_cores_integration)
            for(size_t ind = 0; ind < position_dim; ++ind)
            {   
                std::array<size_t, 2> ij;
                lattice_.inverse_grid_index_position(ind, ij);
                

                Tmunu_from_BoltzmannDist Tmunu(lattice_);
                Landau_matching<energy_momentum_field > Landau_match;
                energy_momentum_field Tmunu_result;
                Tmunu(l_inf, ncore, BD, ij[0], ij[1], Tmunu_result);
                hydrodynamic_field hydField;

                Landau_match( Tmunu_result, hydField );

                Tmunu_Grd[ind] = Tmunu_result;
                Hydro_Grd[ind] = hydField;

            }
            discrete_Tmunu = Tmunu_Grd;
            discrete_hyd = Hydro_Grd;
            // discrete_Tmunu = std::move(Tmunu_Grd);
            // discrete_hyd = std::move(Hydro_Grd);
        }
    };

};


#endif //Tmunu_and_hyd_H