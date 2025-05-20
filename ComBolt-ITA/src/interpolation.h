// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// interpolation.h
#ifndef interpolation_H
#define interpolation_H

#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline.h>

#include "initial_state.h"

struct interpolation1D
{
    void operator()(const std::vector<double> &x_base, const std::vector<double> &x_input, const std::vector<double> &f_input, const gsl_interp_type *gsl_interp_method, std::vector<double> &f_interped)
    {
        size_t n = x_input.size();
        double xmin = x_input.front();
        double xmax = x_input.back();

        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        // gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, n);
        gsl_spline *spline = gsl_spline_alloc(gsl_interp_method, n);

        gsl_spline_init(spline, x_input.data(), f_input.data(), n);

        for (double xi : x_base)
        {
            double yi;
            if (xi < xmin || xi > xmax)
            {
                yi = 0;
            }
            else
            {
                yi = gsl_spline_eval(spline, xi, acc);
            }
            f_interped.push_back(yi);
        }

        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
        
        // if (  isnan(  f_interped[0] ) ) std::cout << " min  " << xmin << " max  " << xmax  << " n "  << n << "  " << f_interped[0] << "  x input size  " << x_input[3]  << "\n"; 
    }
};

struct interpolationBiX
{
    void operator()(const std::vector<double> &x_base, const std::vector<double> &y_base, const std::vector<double> &x_input, const std::vector<double> &y_input, const std::vector<double> &f_input, const gsl_interp2d_type *gsl_interp_method, std::vector<double> &f_interped)
    {
        double xmin = x_input.front();
        double xmax = x_input.back();
        double ymin = y_input.front();
        double ymax = y_input.back();
        
        size_t nx = x_input.size();
        size_t ny = y_input.size();
        // gsl_spline2d *spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, nx, ny);
        gsl_spline2d *spline = gsl_spline2d_alloc(gsl_interp_method, nx, ny);

        gsl_interp_accel *x_acc = gsl_interp_accel_alloc();
        gsl_interp_accel *y_acc = gsl_interp_accel_alloc();

        gsl_spline2d_init(spline, x_input.data(), y_input.data(), f_input.data(), nx, ny);

        for (auto xi : x_base)
        {
            for (auto yi : y_base)
            {
                double zi;
                if (xi < xmin || xi > xmax || yi < ymin || yi > ymax)
                {
                    zi = 0;
                }
                else
                {
                    zi = gsl_spline2d_eval(spline, xi, yi, x_acc, y_acc);
                }
                f_interped.push_back(zi);
            }
        }

        gsl_spline2d_free(spline);
        gsl_interp_accel_free(x_acc);
        gsl_interp_accel_free(y_acc);
    }
};

struct interpolation2D
{
    void operator()(const std::vector<double> &x_base, const std::vector<double> &y_base, const std::vector<double> &x_input, const std::vector<double> &y_input, const std::vector<double> &f_input, const gsl_interp_type *gsl_interp_method, std::vector<double> &f_interped)
    {
        std::vector<std::vector<double>> list_of_f_along_x;
        for (size_t i_y = 0; i_y < y_input.size(); ++i_y)
        {
            std::vector<double> f_input_along_x;
            std::vector<double> f_interped_along_x;
            for (size_t i_x = 0; i_x < x_input.size(); ++i_x)
            {
                f_input_along_x.push_back(f_input[y_input.size() * i_x + i_y]);
            }
            interpolation1D interp1D;
            interp1D(x_base, x_input, f_input_along_x, gsl_interp_method, f_interped_along_x);
            list_of_f_along_x.push_back(f_interped_along_x);
        }
        for (size_t i_x = 0; i_x < x_base.size(); ++i_x)
        {
            std::vector<double> f_input_along_y;
            for (size_t i_y = 0; i_y < y_input.size(); ++i_y)
            {
                f_input_along_y.push_back(list_of_f_along_x[i_y][i_x]);
            }
            interpolation1D interp1D;
            std::vector<double> f_interped_along_y;
            interp1D(y_base, y_input, f_input_along_y, gsl_interp_method, f_interped_along_y);
            f_interped.insert(f_interped.end(), f_interped_along_y.begin(), f_interped_along_y.end());
        }
    }
};




// Function to divide the distance and insert intermediate points
std::vector<double> refineVector(const std::vector<double>& input, const double h, const int number_of_divistion) {
    std::vector<double> refined;

    for (size_t i = 0; i < input.size() - 1; ++i) {
        refined.push_back(input[i]);

        double distance = std::abs(input[i + 1] - input[i]);

        if (distance > h) {
            // Calculate the number of intermediate points
            // int n = std::ceil(distance / h);

            // Insert intermediate points
            double step = distance / number_of_divistion;
            
            for (int j = 1; j < number_of_divistion; ++j) {
                refined.push_back(input[i] + j * step  );
            }
        }
    }

    // Add the last element
    refined.push_back(input.back());
    return refined;
}

void removeCloseIntervals(std::vector<double>& vec, double threshold) {
    if (std::abs(vec[1] - vec[0]) < threshold && vec.size() >= 5) 
        vec.erase(vec.begin() + 1); 
   
}
struct BoltzDist_interpolation
{
    const gsl_interp2d_type *interp2d_method_along_position = gsl_interp2d_bicubic;
    const gsl_interp_type *interp_method_along_momentum     = gsl_interp_cspline;

    void operator()(const parameters_class::multicore_parameters &ncore, 
                    const parameters_class::interpolation_parameters &interp_param,
                    const lattice &lattice_I, const Field_grid &discrBD, lattice &lattice_O, Field_grid &interpedBD_out)
    {

        std::vector<double> x_base = lattice_I.axis.x;
        std::vector<double> y_base = lattice_I.axis.y;
        size_t Nx = lattice_I.axis.x.size();
        size_t Ny = lattice_I.axis.y.size();
        size_t Nvp = lattice_I.axis.vp.size();
        size_t Nvphi = lattice_I.axis.vphi.size();
        // std::vector<double>  vp_base  = lattice_I.axis.vp;

        Field_grid interpedBD(discrBD.size());
        #pragma omp parallel for collapse(2) num_threads(ncore.number_of_cores_interpolation)
        for (size_t ivp = 0; ivp < Nvp; ++ivp)
        {
            for (size_t ivphi = 0; ivphi < Nvphi; ++ivphi)
            {

                std::vector<double> x_input;
                std::vector<double> y_input;

                size_t x_index, y_index;

                for (size_t ix = 0; ix < Nx; ++ix)
                {
                    lattice_I.grid_index(ix, 0, ivp, ivphi, x_index);
                    x_input.push_back(discrBD[x_index].coordinate.x);

                }
                for (size_t iy = 0; iy < Ny; ++iy)
                {
                    lattice_I.grid_index(0, iy, ivp, ivphi, y_index);
                    y_input.push_back(discrBD[y_index].coordinate.y);
                }

                std::vector<double> f_input;

                // BUG IN THE (x,y) -> (y,x) in the following line
                // I used an IF condition to fix the problem manyually.
                for (size_t iy = 0; iy < Ny; ++iy)
                {
                    for (size_t ix = 0; ix < Nx; ++ix)
                    {
                        size_t ind;
                        lattice_I.grid_index(ix, iy, ivp, ivphi, ind);
                        interpedBD[ind].tau = discrBD[ind].tau;
                        interpedBD[ind].coordinate.x = x_base[ix];
                        interpedBD[ind].coordinate.y = y_base[iy];
                        interpedBD[ind].coordinate.vp = discrBD[ind].coordinate.vp;
                        interpedBD[ind].coordinate.vphi = discrBD[ind].coordinate.vphi;

                        f_input.push_back(discrBD[ind].value);
                    }
                }

                std::vector<double> f_interped;

                interpolationBiX interp2D;
                interp2D(x_base, y_base, x_input, y_input, f_input, interp2d_method_along_position, f_interped);

                for (size_t ix = 0; ix < Nx; ++ix)
                {
                    for (size_t iy = 0; iy < Ny; ++iy)
                    {
                        size_t ind;
                        lattice_I.grid_index(ix, iy, ivp, ivphi, ind);
                        interpedBD[ind].value = f_interped[Ny * ix + iy];
                    }
                }
            }
        }
        lattice_O.axis    = lattice_I.axis;
        // interpedBD_out = interpedBD;
		
		// We match the new lattice with the shifted vp grid (We should be able to that in the previous loop)
        lattice_O.axis = lattice_I.axis;
        lattice_O.axis.vp.clear();
        for (size_t ivp = 0; ivp < lattice_I.axis.vp.size(); ++ivp)
        {
            size_t ind;
            lattice_I.grid_index(0, 0, ivp, 0, ind);
            lattice_O.axis.vp.push_back(interpedBD[ind].coordinate.vp);
        }


        lattice_O.axis.vp = refineVector(lattice_O.axis.vp, interp_param.insert_threshold, interp_param.number_of_divistion);
        // removeCloseIntervals(lattice_O.axis.vp, interp_param.remove_threshold);

        // parameters_class::shell_message shell_msg;
        // if (!shell_msg.quiet)
        //     std::cout << "Size of vp grid:  " << lattice_O.axis.vp.size() << " \n";
        
        interpedBD_out.resize(Nx * Ny * lattice_O.axis.vp.size() * Nvphi);
        std::vector<double> vp_base = lattice_O.axis.vp;
        
        #pragma omp parallel for collapse(3) num_threads(ncore.number_of_cores_interpolation)
        for (size_t ix = 0; ix < Nx; ++ix)
        {
            for (size_t iy = 0; iy < Ny; ++iy)
            {
                for (size_t ivphi = 0; ivphi < Nvphi; ++ivphi)
                {
                    std::vector<double> vp_input;
                    std::vector<double> f_input;

                    for (size_t ivp = 0; ivp < Nvp; ++ivp)
                    {
                        size_t ind;
                        lattice_I.grid_index(ix, iy, ivp, ivphi, ind);
                        vp_input.push_back(discrBD[ind].coordinate.vp);

                        interpedBD[ind].coordinate.vp = vp_base[ivp];
                        f_input.push_back(interpedBD[ind].value);
                    }
                    std::vector<double> f_interped;
                    interpolation1D interp1D;
                    interp1D(vp_base, vp_input, f_input, interp_method_along_momentum, f_interped);
                    for (size_t ivp = 0; ivp < lattice_O.axis.vp.size(); ++ivp)
                    {
                        size_t ind;
                        lattice_O.grid_index(ix, iy, ivp, ivphi, ind);
                        interpedBD_out[ind].tau = interpedBD[0].tau;
                        interpedBD_out[ind].coordinate.x = lattice_O.axis.x[ix];
                        interpedBD_out[ind].coordinate.y = lattice_O.axis.y[iy];
                        interpedBD_out[ind].coordinate.vp = lattice_O.axis.vp[ivp];
                        interpedBD_out[ind].coordinate.vphi = lattice_O.axis.vphi[ivphi];
                        interpedBD_out[ind].value = f_interped[ivp];
                    }
                }
            }
        }

    }
};

#endif // interpolation_H
