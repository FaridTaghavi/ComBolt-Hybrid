// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.


// grid.h
#ifndef grid_H
#define grid_H


#include "parameters_class.h"
#include <vector>
#include <array>
#include <cmath>
#include <iostream>

struct Boltzmann_dist_field 
{
    double tau {0.}, x {0.}, y {0.}, vp {0.}, vphi {0.}, F {0.};
};
struct coord
{
    double x {0.}, y{0.}, vp {0.}, vphi {0.};
};


/* We can extend the following structure  */
struct hydrodynamic_field
{
    double tau, x, y;
    double ed {0.};
    double Px {0.};
    double Py {0.};
    double Pz {0.};
    //u must consider with upper index
    double utau {0.};
    double ux   {0.};
    double uy   {0.};
    double ueta {0.};
    //Extend if needed
    // double SheatTensorXX {0.};
    double BulkTensorXX  {0.};
    // double pi_tx {0.};
    // double pi_ty {0.};
    // double pi_tz {0.};
    // double pi_xy {0.};
    // double pi_xz {0.};
    // double pi_yz {0.};
    double pi_tautau {0.};
    double pi_taux   {0.};
    double pi_tauy   {0.};
    double pi_taueta {0.};
    double pi_xx     {0.};
    double pi_xy     {0.};
    double pi_xeta   {0.};
    double pi_yy     {0.};
    double pi_yeta   {0.};
    double pi_etaeta {0.};

    double inverseReyn  {0.}; 
    double isotropicity {0.};
    // bool   freeze_hyper {false};
    // bool   isotro_hyper {false};
};

struct energy_momentum_field
{
    double tau, x, y;
    double Ttautau {0.};
    double Ttaux   {0.};
    double Ttauy   {0.};
    double Txx {0.};
    double Txy {0.};
    double Tyy {0.};
    double Ttaueta {0.};
    double Txeta   {0.};
    double Tyeta   {0.};
    double Tetaeta {0.};

    energy_momentum_field operator+(const energy_momentum_field& other) const {
        return {
            tau, x, y, // Assume tau, x, y are unchanged or handled separately
            Ttautau + other.Ttautau,
            Ttaux + other.Ttaux,
            Ttauy + other.Ttauy,
            Txx + other.Txx,
            Txy + other.Txy,
            Tyy + other.Tyy,
            Ttaueta + other.Ttaueta,
            Txeta + other.Txeta,
            Tyeta + other.Tyeta,
            Tetaeta + other.Tetaeta
        };
    }
    energy_momentum_field& operator+=(const energy_momentum_field& other) {
        Ttautau += other.Ttautau;
        Ttaux += other.Ttaux;
        Ttauy += other.Ttauy;
        Txx += other.Txx;
        Txy += other.Txy;
        Tyy += other.Tyy;
        Ttaueta += other.Ttaueta;
        Txeta += other.Txeta;
        Tyeta += other.Tyeta;
        Tetaeta += other.Tetaeta;
        return *this;
    }
    energy_momentum_field operator*(const double &scalar) const {
        return {
            tau, x, y, // tau, x, y are not scaled
            Ttautau * scalar,
            Ttaux * scalar,
            Ttauy * scalar,
            Txx * scalar,
            Txy * scalar,
            Tyy * scalar,
            Ttaueta * scalar,
            Txeta * scalar,
            Tyeta * scalar,
            Tetaeta * scalar
        };

    }

};   

struct field
{
    double tau {0.};//, x {0.}, y {0.}, vp {0.}, vphi {0.};
    coord coordinate = {0., 0., 0., 0.};
    double value {0.};

};


struct observables
{
    double tau;
    double eps_tot    {0};
    double R          {0};
    double eps_2_x    {0};
    double eps_2_y    {0};
    double eps_p_x    {0};
    double eps_p_y    {0};
    double u_perp     {0};
    double av_invReyn {0};
    double Etr {0};
};

struct freeze_out_vars
{
    std::array<double, 3> xmu      {{0., 0., 0.}}; // (t, x, y) in fm
    std::array<double, 3> sigma_mu {{0., 0., 0.}}; // covariant normal vector fm^2
    std::array<double, 2> vi       {{0., 0.}};     // (vx,vy)
    std::array<double, 3> pi_ij    {{0., 0., 0.}}; // shear tensor (pi_xx, pi_xy, pi_yy)
    double             Pi_bulk     {0.};     // bulk tensor 
};

/* The grid is a 1D array each element contains the information of the field.   */
using Field_grid = std::vector< field >;
using Hydro_grid = std::vector< hydrodynamic_field >;
using Tmunu_grid = std::vector< energy_momentum_field >;
using obs_evol = std::vector< observables >;
using freeze_out_surface = std::vector< freeze_out_vars >;

struct axes_grid
{
    std::vector<double> x, y, vp, vphi;
};

struct lattice
{
    axes_grid axis;
    
    lattice(axes_grid axis_ ) : axis(axis_) {}
        
    void initiate_lattice(  const parameters_class::lattice_info &l_inf,
                            const parameters_class::initial_values &init_val,
                            const parameters_class::trento_info &trento_inf, 
                            axes_grid &ax_grid)
    {
        
        double dx, dy, nx, ny, xmin, ymin; 
        
        if (init_val.initialization == 1){
            
            double grid_max = trento_inf.grid_max;
            double grid_step = trento_inf.grid_step;
            
            dx = grid_step;
            dy = grid_step;
            
			nx = static_cast<size_t>(2*grid_max / dx) ;
            ny = static_cast<size_t>(2*grid_max / dy) ;
            
			xmin = -grid_max + dx / 2;
            ymin = -grid_max + dy / 2;
        
        } else if (init_val.initialization == 0 || init_val.initialization == 2) {
            dx = (l_inf.x_max - l_inf.x_min) / (l_inf.n_x - 1);
            dy = (l_inf.y_max - l_inf.y_min) / (l_inf.n_y - 1);
            nx = l_inf.n_x;
            ny = l_inf.n_y;
            xmin = l_inf.x_min;
            ymin = l_inf.y_min;
        } else {
            std::cerr << "The lattice is not defined.\n";
            return;
        }
        
        double dv_phi = (l_inf.v_phi_max + 2 * M_PI / l_inf.n_v_phi - l_inf.v_phi_min ) / (l_inf.n_v_phi);
        for (size_t ix = 0; ix < nx; ++ix)
        {
            axis.x.push_back(xmin + ix * dx); 
        }
        for (size_t iy = 0; iy <ny; ++iy)
        {
            axis.y.push_back(ymin + iy * dy); 
        }
        for (size_t ivphi = 0; ivphi < l_inf.n_v_phi; ++ivphi)
        {
            axis.vphi.push_back(l_inf.v_phi_min + ivphi * dv_phi); 
        }
        for (size_t ivp = 0; ivp < l_inf.n_v_p; ++ivp)
        {
            size_t n_of_intervals {l_inf.n_v_p - 1};
            double q =  static_cast<double>(ivp) / n_of_intervals;
            axis.vp.push_back( q ); 
        }
        if (!l_inf.Z2_symmetry)
        {
            std::vector <double> vp_copy = axis.vp;
            for (auto vp : axis.vp) {
                if (vp != 0.0) {  // Avoid duplicating zero
                    vp_copy.insert(vp_copy.begin(),  -vp);
                }
            }
            axis.vp = vp_copy;
        }
   
    };
    
    void grid_index(const size_t &i_x, const size_t &i_y, const size_t &i_v_p, const size_t &i_v_phi, size_t &index) const 
    {
        size_t Nx = axis.x.size();
        size_t Ny = axis.y.size();
        size_t Nvp = axis.vp.size();
        size_t Nvphi = axis.vphi.size();

        // std::cout << "---------------------------size of lattice: " << Nx << "    " << Ny << "    " << Nvp << "    " << Nvphi << "\n";
        
        if (i_x >= Nx || i_y >= Ny || i_v_p >= Nvp || i_v_phi >= Nvphi){
            std::cerr << "grid_index: Index overflow\n";
            throw std::runtime_error("grid_index: Index overflow");
            // std::cout << "Index overflow\n";
            // return;
        };

        index = i_x * Ny * Nvp * Nvphi + i_y * Nvp * Nvphi + i_v_p * Nvphi + i_v_phi;
    };

    void grid_index_position(const size_t &i_x, const size_t &i_y, size_t &index) const 
    {
        size_t Nx = axis.x.size();
        size_t Ny = axis.y.size();
        // std::cout << "---------------------------size of lattice: " << Nx << "    " << Ny << "\n";
        
        if (i_x >= Nx || i_y >= Ny){
            std::cerr << "grid_index: Index overflow\n";
            throw std::runtime_error("grid_index: Index overflow");
            // std::cout << "Index overflow\n";
            // return;
        };

        index = i_x * Ny + i_y;
    };

    void inverse_grid_index_position (const size_t &ind, std::array<size_t, 2> &output) const
    {
        size_t i,j;
        size_t dim2;

        dim2 = axis.y.size();
        j = ind % dim2;
        i = (ind - j) / dim2;

        // output = {i,j}; 
        output = {j,i}; // IMPORTANT!!!! Why TRANSPOSE is required? DEBUG IT LATER!!
    };
    
    void inverse_grid_index(const size_t &ind_, std::array<size_t, 4> &output) const
    {
        size_t i,j,k,l;
        size_t dim2, dim3, dim4;
        size_t ind = ind_;
        dim2 = axis.y.size();
        dim3 = axis.vp.size();
        dim4 = axis.vphi.size();

        l = ind % dim4;
        ind = (ind - l) / dim4;
        k = ind % dim3;
        ind = (ind - k) / dim3;
        j = ind % dim2;
        ind = (ind - j) / dim2;
        i = ind;

        output = {i,j,k,l};
    };

    void project_to_poistion_grid_index(const size_t &ind, size_t &output) const
    {
        size_t dim2 = axis.y.size();
        std::array<size_t, 4> ijkl;
        inverse_grid_index(ind, ijkl);
        output = ijkl[0] + dim2 * ijkl[1];
    }

};



#endif //grid
