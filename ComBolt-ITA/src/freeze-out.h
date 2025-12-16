// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// M.A. Tabatabee, S.F. Taghavi
#ifndef freeze_out_H
#define freeze_out_H

#include "parameters_class.h"
#include "Administration.h"
#include <iomanip>
#include "grid.h"
#include "cornelius.h"


struct surface_element
{
    double c00, c01, c10, c11;
};

// void trilinear_interpolation( const double &tau, const double &x, const double &y,
//                               const double &tau0, const double &tau1, 
//                               const double &x0, const double &x1, 
//                               const double &y0, const double &y1,
//                               const double &f000, const double &f001, const double &f010, const double &f011,
//                               const double &f100, const double &f101, const double &f110, const double &f111,
//                               double &f_interped   
//                             )
// {
//     // Wikipedia: https://en.wikipedia.org/wiki/Trilinear_interpolation
//     double dt = tau1 - tau0;
//     double dx = x1 - x0;
//     double dy = y1 - y0;
// 
//     double xd = (x - x0) / dx;
//     double yd = (y - y0) / dy;
//     double td = (tau - tau0) / dt;
// 
//     double c00 = f000 * (1 - td) + f100 * td;
//     double c01 = f001 * (1 - td) + f101 * td;
//     double c10 = f010 * (1 - td) + f110 * td;
//     double c11 = f011 * (1 - td) + f111 * td;
// 
//     double c0 = c00 * (1 - xd) + c10 * xd;
//     double c1 = c01 * (1 - xd) + c11 * xd;
// 
//     f_interped = c0 * (1 - yd) + c1 * yd;
// };

void trilinear_interpolation_new( const double &tau, const double &x, const double &y,
                              const double &tau0, const double &tau1, 
                              const double &x0, const double &x1, 
                              const double &y0, const double &y1,
                              const std::array<surface_element, 2> & vi_tau0, const std::array<surface_element, 2> & vi_tau1,
                              const std::array<surface_element, 3> & pi_ij_tau0, const std::array<surface_element, 3> & pi_ij_tau1,
                              const surface_element & Pi_bulk_tau0, const surface_element & Pi_bulk_tau1,
                              std::array<double, 2> & vi_interped,
                              std::array<double, 3> & pi_ij_interped,
                              double &Pi_bulk_interped
                            )
{
    //The function performances need to be validated.
    // Wikipedia: https://en.wikipedia.org/wiki/Trilinear_interpolation
    double dt = tau1 - tau0;
    double dx = x1 - x0;
    double dy = y1 - y0;

    double xd = (x - x0) / dx;
    double yd = (y - y0) / dy;
    double td = (tau - tau0) / dt;

    std::array<surface_element, 2> vi_;
    std::array<double, 2> vi_c0, vi_c1;

    std::array<surface_element, 3> pi_ij_;
    std::array<double, 3> pi_ij_c0, pi_ij_c1;
    
    surface_element Pi_bulk_;
    double Pi_bulk_c0, Pi_bulk_c1;
    // std::cout << "---------------------------------------------\n";
    // std::cout << "tau0: " << tau0 << " tau1: " << tau1 << " x0 " << x0 << " x1: " << x1 << " y0: " << y0 << " y1: " << y1 << "\n";
    // std::cout << "xd: " << xd << " yd: " << yd << " td: " << td << "\n"; 
    for (int i = 0; i < 2; ++i)
    {
        vi_[i].c00 = vi_tau0[i].c00 * (1 - td) + vi_tau1[i].c00 * td;
        vi_[i].c01 = vi_tau0[i].c01 * (1 - td) + vi_tau1[i].c01 * td;
        vi_[i].c10 = vi_tau0[i].c10 * (1 - td) + vi_tau1[i].c10 * td;
        vi_[i].c11 = vi_tau0[i].c11 * (1 - td) + vi_tau1[i].c11 * td;
        
        vi_c0[i] = vi_[i].c00 * (1 - xd) + vi_[i].c10 * xd;
        vi_c1[i] = vi_[i].c01 * (1 - xd) + vi_[i].c11 * xd;
        
        vi_interped[i] = vi_c0[i] * (1 - yd) + vi_c1[i] * yd;
    }


    for (int i = 0; i < 3; ++i)
    {
        pi_ij_[i].c00 = pi_ij_tau0[i].c00 * (1 - td) + pi_ij_tau1[i].c00 * td;
        pi_ij_[i].c01 = pi_ij_tau0[i].c01 * (1 - td) + pi_ij_tau1[i].c01 * td;
        pi_ij_[i].c10 = pi_ij_tau0[i].c10 * (1 - td) + pi_ij_tau1[i].c10 * td;
        pi_ij_[i].c11 = pi_ij_tau0[i].c11 * (1 - td) + pi_ij_tau1[i].c11 * td;
        
        pi_ij_c0[i] = pi_ij_[i].c00 * (1 - xd) + pi_ij_[i].c10 * xd;
        pi_ij_c1[i] = pi_ij_[i].c01 * (1 - xd) + pi_ij_[i].c11 * xd;
        
        pi_ij_interped[i] = pi_ij_c0[i] * (1 - yd) + pi_ij_c1[i] * yd;
    }
    
    Pi_bulk_.c00 = Pi_bulk_tau0.c00 * (1 - td) + Pi_bulk_tau1.c00 * td;
    Pi_bulk_.c01 = Pi_bulk_tau0.c01 * (1 - td) + Pi_bulk_tau1.c01 * td;
    Pi_bulk_.c10 = Pi_bulk_tau0.c10 * (1 - td) + Pi_bulk_tau1.c10 * td;
    Pi_bulk_.c11 = Pi_bulk_tau0.c11 * (1 - td) + Pi_bulk_tau1.c11 * td;

    Pi_bulk_c0 = Pi_bulk_.c00 * (1 - xd) + Pi_bulk_.c10 * xd;
    Pi_bulk_c1 = Pi_bulk_.c01 * (1 - xd) + Pi_bulk_.c11 * xd;

    Pi_bulk_interped = Pi_bulk_c0 * (1 - yd) + Pi_bulk_c1 * yd;
};

struct cornelius_on_action
{
    const int dimension {3};


    void operator()(const parameters_class  &param,
                    const double &h,
                    const lattice &lattice_I, 
                    const Hydro_grid &discrete_hyd, 
                    const Hydro_grid &discret_hyd_prev_step,
                    const std::string path_surface,
                    bool &found_freeze_out
                )    
    {
                bool intersect = true;
                double epsilon_FO, dummuy_P;

                if (param.dissipation_info_.dissipation_mode == 0 || param.dissipation_info_.dissipation_mode == 1)
                {
                    epsilon_FO = param.medium_info_.epsFO; // 1/fm^4
                } else if (param.dissipation_info_.dissipation_mode == 2) {
                    lattice_ed_pressure( param.medium_info_.TFO, 
                                        epsilon_FO, 
                                        dummuy_P);
                    // std::cout << "Freeze-out energy " << epsilon_FO << "\n";
                } else {
                    std::cerr << "The dissipation mode is not defined.\n";
                    return;
                }
                
                const int step = 1;
                const double DTau = 1.0 * h; //Time spacing, we keep as the evolution time step, since it is dynamically changing
                const double DX   = 1.0 * (lattice_I.axis.x[step] - lattice_I.axis.x[0]);
                const double DY   = 1.0 * (lattice_I.axis.y[step] - lattice_I.axis.y[0]);
                
                double dx[3] = {DTau, DX, DY};
                Cornelius* cornelius_ptr = new Cornelius();
                cornelius_ptr->init(dimension, epsilon_FO, dx);
            
                double ***cube = new double ** [2];
                for (int i = 0; i < 2; i++){
                    cube[i] = new double * [2];
                    for (int j = 0; j < 2; j++){
                        cube[i][j] = new double[2];
                        for (int k = 0; k < 2; k++){
                            cube[i][j][k] = 0.0;
                        }
                    }
                }
                found_freeze_out = false;

                for (size_t ix=0; ix < lattice_I.axis.x.size() - 1; ix += step){
                    for (size_t iy=0; iy < lattice_I.axis.y.size() - 1; iy += step){
                        // ------------------- AGAIN in hyd the ij -> ji problem!! -------------------
                        size_t hyd_ind_00, hyd_ind_01, hyd_ind_10, hyd_ind_11;
                        lattice_I.grid_index_position(ix, iy, hyd_ind_00);
                        lattice_I.grid_index_position(ix, iy+1, hyd_ind_10); // Fix the bug here
                        lattice_I.grid_index_position(ix+1, iy, hyd_ind_01);  // Fix the bug here
                        lattice_I.grid_index_position(ix+1, iy+1, hyd_ind_11); 
                        intersect=true;
                        if ((discrete_hyd[hyd_ind_11].ed - epsilon_FO)*(discret_hyd_prev_step[hyd_ind_00].ed - epsilon_FO) > 0.)
                            if ((discrete_hyd[hyd_ind_10].ed - epsilon_FO)*(discret_hyd_prev_step[hyd_ind_01].ed - epsilon_FO) > 0.)
                                if ((discrete_hyd[hyd_ind_01].ed - epsilon_FO)*(discret_hyd_prev_step[hyd_ind_10].ed - epsilon_FO) > 0.)
                                    if ((discrete_hyd[hyd_ind_00].ed - epsilon_FO)*(discret_hyd_prev_step[hyd_ind_11].ed - epsilon_FO) > 0.)
                                        intersect = 0;
                        if (intersect == false) continue;
                        found_freeze_out = true;

                        cube[0][0][0] = discret_hyd_prev_step[hyd_ind_00].ed;
                        cube[0][0][1] = discret_hyd_prev_step[hyd_ind_01].ed;
                        cube[0][1][0] = discret_hyd_prev_step[hyd_ind_10].ed;
                        cube[0][1][1] = discret_hyd_prev_step[hyd_ind_11].ed;
                        cube[1][0][0] = discrete_hyd[hyd_ind_00].ed;
                        cube[1][0][1] = discrete_hyd[hyd_ind_01].ed;
                        cube[1][1][0] = discrete_hyd[hyd_ind_10].ed;
                        cube[1][1][1] = discrete_hyd[hyd_ind_11].ed;

                        cornelius_ptr->find_surface_3d(cube);
                        
                        // Sometimes there are more that one surface element exists for a given cell. kk is the index of the surface element
                        for (int kk = 0; kk < cornelius_ptr->get_Nelements(); kk++){
                            freeze_out_vars surface;
                            // std::cout << "Number of elements: " << cornelius_ptr->get_Nelements() << std::endl;
                            for (int ii = 0; ii < dimension; ii++){
                                // vol_element[ii] = cornelius_ptr->get_normal_elem(kk, ii);
                                // int sign = (ii == 0) ? 1 : -1;
                                int sign = 1;
                                surface.sigma_mu[ii] = sign * cornelius_ptr->get_normal_elem(kk, ii);
                            }
                            
                            // The time of the previous step shifted to centrois
                            surface.xmu[0] = discret_hyd_prev_step[hyd_ind_00].tau + cornelius_ptr->get_centroid_elem(kk, 0);  // Tau
                            // The location of lower left corner is shifted to centroid
                            surface.xmu[1] = discrete_hyd[hyd_ind_00].x + cornelius_ptr->get_centroid_elem(kk, 1);          // x
                            surface.xmu[2] = discrete_hyd[hyd_ind_00].y + cornelius_ptr->get_centroid_elem(kk, 2);          // y 

                            // ============ Finding the interpolated vatiables over the surface ============
                            // ------------------- AGAIN in hyd the ij -> ji problem!! -------------------
                            const double tau0 = discret_hyd_prev_step[hyd_ind_11].tau;
                            const double tau1 = discrete_hyd[hyd_ind_11].tau;
                            const double x0 = discrete_hyd[hyd_ind_00].x;  
                            // const double x1 = discrete_hyd[hyd_ind_01].x; // Fix the bug here
                            const double x1 = discrete_hyd[hyd_ind_10].x; 
                            const double y0 = discrete_hyd[hyd_ind_00].y;
                            // const double y1 = discrete_hyd[hyd_ind_10].y; // Fix the bug here
                            const double y1 = discrete_hyd[hyd_ind_01].y; 
                            
                            const surface_element vx_tau0 = { 
                                discret_hyd_prev_step[hyd_ind_00].ux / discret_hyd_prev_step[hyd_ind_00].utau,
                                discret_hyd_prev_step[hyd_ind_01].ux / discret_hyd_prev_step[hyd_ind_01].utau,
                                discret_hyd_prev_step[hyd_ind_10].ux / discret_hyd_prev_step[hyd_ind_10].utau,
                                discret_hyd_prev_step[hyd_ind_11].ux / discret_hyd_prev_step[hyd_ind_11].utau
                            };
                            const surface_element vy_tau0 = { 
                                discret_hyd_prev_step[hyd_ind_00].uy / discret_hyd_prev_step[hyd_ind_00].utau,
                                discret_hyd_prev_step[hyd_ind_01].uy / discret_hyd_prev_step[hyd_ind_01].utau,
                                discret_hyd_prev_step[hyd_ind_10].uy / discret_hyd_prev_step[hyd_ind_10].utau,
                                discret_hyd_prev_step[hyd_ind_11].uy / discret_hyd_prev_step[hyd_ind_11].utau
                            };
                            const std::array<surface_element, 2> vi_tau0 = {vx_tau0, vy_tau0};
                            
                            const surface_element vx_tau1 = { 
                                discrete_hyd[hyd_ind_00].ux / discrete_hyd[hyd_ind_00].utau,
                                discrete_hyd[hyd_ind_01].ux / discrete_hyd[hyd_ind_01].utau,
                                discrete_hyd[hyd_ind_10].ux / discrete_hyd[hyd_ind_10].utau,
                                discrete_hyd[hyd_ind_11].ux / discrete_hyd[hyd_ind_11].utau
                            };
                            const surface_element vy_tau1 = { 
                                discrete_hyd[hyd_ind_00].uy / discrete_hyd[hyd_ind_00].utau,
                                discrete_hyd[hyd_ind_01].uy / discrete_hyd[hyd_ind_01].utau,
                                discrete_hyd[hyd_ind_10].uy / discrete_hyd[hyd_ind_10].utau,
                                discrete_hyd[hyd_ind_11].uy / discrete_hyd[hyd_ind_11].utau
                            };
                            const std::array<surface_element, 2> vi_tau1 = {vx_tau1, vy_tau1};

                            const surface_element pi_xx_tau0 = { discret_hyd_prev_step[hyd_ind_00].pi_xx, discret_hyd_prev_step[hyd_ind_01].pi_xx,
                                                                 discret_hyd_prev_step[hyd_ind_10].pi_xx, discret_hyd_prev_step[hyd_ind_11].pi_xx };
                            
                            const surface_element pi_xy_tau0 = { discret_hyd_prev_step[hyd_ind_00].pi_xy, discret_hyd_prev_step[hyd_ind_01].pi_xy,
                                                                 discret_hyd_prev_step[hyd_ind_10].pi_xy, discret_hyd_prev_step[hyd_ind_11].pi_xy };                      

                            const surface_element pi_yy_tau0 = { discret_hyd_prev_step[hyd_ind_00].pi_yy, discret_hyd_prev_step[hyd_ind_01].pi_yy,
                                                                 discret_hyd_prev_step[hyd_ind_10].pi_yy, discret_hyd_prev_step[hyd_ind_11].pi_yy };
                            const std::array<surface_element, 3> pi_ij_tau0 = {pi_xx_tau0, pi_xy_tau0, pi_yy_tau0};
                            
                            const surface_element pi_xx_tau1 = { discrete_hyd[hyd_ind_00].pi_xx, discrete_hyd[hyd_ind_01].pi_xx ,
                                                                 discrete_hyd[hyd_ind_10].pi_xx, discrete_hyd[hyd_ind_11].pi_xx };

                            const surface_element pi_xy_tau1 = { discrete_hyd[hyd_ind_00].pi_xy, discrete_hyd[hyd_ind_01].pi_xy,
                                                                 discrete_hyd[hyd_ind_10].pi_xy, discrete_hyd[hyd_ind_11].pi_xy };

                            const surface_element pi_yy_tau1 = { discrete_hyd[hyd_ind_00].pi_yy, discrete_hyd[hyd_ind_01].pi_yy,
                                                                 discrete_hyd[hyd_ind_10].pi_yy, discrete_hyd[hyd_ind_11].pi_yy };

                            const std::array<surface_element, 3> pi_ij_tau1 = {pi_xx_tau1, pi_xy_tau1, pi_yy_tau1};
                            
                            const surface_element Pi_bulk_tau0 = { discret_hyd_prev_step[hyd_ind_00].BulkTensorXX, discret_hyd_prev_step[hyd_ind_01].BulkTensorXX,
                                                                   discret_hyd_prev_step[hyd_ind_10].BulkTensorXX, discret_hyd_prev_step[hyd_ind_11].BulkTensorXX };
                            
                            const surface_element Pi_bulk_tau1 = { discrete_hyd[hyd_ind_00].BulkTensorXX, discrete_hyd[hyd_ind_01].BulkTensorXX,
                                                                   discrete_hyd[hyd_ind_10].BulkTensorXX, discrete_hyd[hyd_ind_11].BulkTensorXX };

                            // std::cout << " BulkTensorXX: " << discret_hyd_prev_step[hyd_ind_11].BulkTensorXX << "  "
                            //           << discrete_hyd[hyd_ind_11].BulkTensorXX << std::endl;

                            trilinear_interpolation_new( surface.xmu[0], surface.xmu[1], surface.xmu[2],
                                                            tau0, tau1, 
                                                            x0, x1, 
                                                            y0, y1,
                                                            vi_tau0, vi_tau1,
                                                            pi_ij_tau0, pi_ij_tau1,
                                                            Pi_bulk_tau0, Pi_bulk_tau1,
                                                            surface.vi,
                                                            surface.pi_ij,
                                                            surface.Pi_bulk
                                                            );  

                            // std::cout << "vx_tau0: " << vx_tau0.c00 << "  " << vx_tau0.c01 << "  " << vx_tau0.c10 << "  " << vx_tau0.c11 << std::endl;
                            // std::cout << "vx_tau1: " << vx_tau1.c00 << "  " << vx_tau1.c01 << "  " << vx_tau1.c10 << "  " << vx_tau1.c11 << std::endl;
                            // std::cout << "vx_interped: " << surface.vi[0] << "   " << surface.vi[1] << std::endl;
                            //double interpedF; 
                            //trilinear_interpolation(0.3, 0.42, 0.7,
                            //                        0.11, 0.36,
                            //                        0.3 , 0.43,
                            //                        0.01, 7.4,
                            //                        0.09801,4.86873e-6, 0.413675, 0.0000205493,
                            //                        1.04976, 0.0000520514, 4.43074, 0.000219264,
                            //                        interpedF);

                            //std::cout << " triginterp:  "<< std::setprecision(10)  << interpedF << std::endl;
                            
                            save_surface(surface, path_surface);

                        }
                    }
    
                }               
                for (int i = 0; i < 2; i++){
                    for (int j = 0; j < 2; j++){
                        delete [] cube[i][j];
                    }
                    delete [] cube[i];
                }
                delete [] cube;
                delete cornelius_ptr; 
                // std::cout << "Number of elements: " << cornelius_ptr->get_Nelements() << " intersect  " << found_freeze_out << std::endl;
             }
    };

#endif // freeze_out_H 
