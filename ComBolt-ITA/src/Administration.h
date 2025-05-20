// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// administration.h
#ifndef ADMINISTRATION_H
#define ADMINISTRATION_H


#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <functional>
#include <filesystem>
#include <chrono>
#include "grid.h"

 //include <H5Cpp.h> 

namespace fs = std::filesystem;

// Enum to specify operation type
enum Operation { ADD, SUBTRACT, MULTIPLE, DIVIDE };

struct path_struct {
    fs::path save_path = fs::current_path();
    fs::path load_path = fs::current_path();
};

void save_BD(const Field_grid&  grid, const std::string& filename)
{

    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    const int fieldWidth = 25;

    for (const auto& f : grid) {
            
        outfile  
        << std::left << std::setw(fieldWidth) << f.tau 
        << std::left << std::setw(fieldWidth) << f.coordinate.x 
        << std::left << std::setw(fieldWidth) << f.coordinate.y 
        << std::left << std::setw(fieldWidth) << std::setprecision(15) << f.coordinate.vp 
        << std::left << std::setw(fieldWidth) << f.coordinate.vphi 
        << std::left << std::setw(fieldWidth) << f.value 
        << std::endl;
    }

    outfile.close();
    if (!outfile) {
        std::cerr << "Error writing to file: " << filename << std::endl;
    }
};

void save_hydro(const Hydro_grid& grid, const std::string& filename)
{
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
        
    const int fieldWidth = 17;

    for (const auto& f : grid) {
            
        outfile  
        << std::left << std::setw(fieldWidth) << f.tau 
        << std::left << std::setw(fieldWidth) << f.x 
        << std::left << std::setw(fieldWidth) << f.y 
        << std::left << std::setw(fieldWidth) << f.ed 
        << std::left << std::setw(fieldWidth) << f.Px 
        << std::left << std::setw(fieldWidth) << f.Py 
        << std::left << std::setw(fieldWidth) << f.Pz 
        << std::left << std::setw(fieldWidth) << f.utau 
        << std::left << std::setw(fieldWidth) << f.ux 
        << std::left << std::setw(fieldWidth) << f.uy 
        << std::left << std::setw(fieldWidth) << f.ueta 
        << std::left << std::setw(fieldWidth) << f.pi_tautau 
        << std::left << std::setw(fieldWidth) << f.pi_taux 
        << std::left << std::setw(fieldWidth) << f.pi_tauy 
        << std::left << std::setw(fieldWidth) << f.pi_taueta 
        << std::left << std::setw(fieldWidth) << f.pi_xx 
        << std::left << std::setw(fieldWidth) << f.pi_xy 
        << std::left << std::setw(fieldWidth) << f.pi_xeta 
        << std::left << std::setw(fieldWidth) << f.pi_yy 
        << std::left << std::setw(fieldWidth) << f.pi_yeta 
        << std::left << std::setw(fieldWidth) << f.pi_etaeta 

        << std::left << std::setw(fieldWidth) << f.inverseReyn 
        << std::left << std::setw(fieldWidth) << f.isotropicity 
        // << std::left << std::setw(fieldWidth) << f.pi_tx 
        // << std::left << std::setw(fieldWidth) << f.pi_ty 
        // << std::left << std::setw(fieldWidth) << f.pi_tz 
        // << std::left << std::setw(fieldWidth) << f.pi_xy 
        // << std::left << std::setw(fieldWidth) << f.pi_xz 
        // << std::left << std::setw(fieldWidth) << f.pi_yz 
        // << std::left << std::setw(fieldWidth) << f.BulkTensorXX 
        << std::endl;
    }


    outfile.close();
    if (!outfile) {
        std::cerr << "Error writing to file: " << filename << std::endl;
    }
};
void save_Tmunu(const Tmunu_grid& grid, const std::string& filename)
{
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
        
    const int fieldWidth = 17;

    for (const auto& f : grid) {
            
        outfile  
        << std::left << std::setw(fieldWidth) << f.tau 
        << std::left << std::setw(fieldWidth) << f.x 
        << std::left << std::setw(fieldWidth) << f.y 
        << std::left << std::setw(fieldWidth) << f.Ttautau 
        << std::left << std::setw(fieldWidth) << f.Ttaux 
        << std::left << std::setw(fieldWidth) << f.Ttauy 
        << std::left << std::setw(fieldWidth) << f.Ttaueta 
        << std::left << std::setw(fieldWidth) << f.Txx 
        << std::left << std::setw(fieldWidth) << f.Txy 
        << std::left << std::setw(fieldWidth) << f.Txeta 
        << std::left << std::setw(fieldWidth) << f.Tyy 
        << std::left << std::setw(fieldWidth) << f.Tyeta 
        << std::left << std::setw(fieldWidth) << f.Tetaeta 
        << std::endl;
    }


    outfile.close();
    if (!outfile) {
        std::cerr << "Error writing to file: " << filename << std::endl;
    }
};

void save_fields_new(const parameters_class::saving_info &save_info,
                     const Field_grid &BD, 
                     const Hydro_grid Hydro, 
                     const Tmunu_grid Tmunu, 
                     const std::string &folder_path, 
                     const size_t index)
{

    if (save_info.save_BD)
    {
        std::string file = "New_Boltzmann_Evol_" + std::to_string(  index ) + ".dat";
        save_BD(BD, folder_path + file);
    }
            

    if (save_info.save_hydro)
    {
        std::string file_hyd = "New_Hydro_Evol_" + std::to_string( index ) + ".dat";
        save_hydro(Hydro, folder_path + file_hyd);
    }

    if (save_info.save_Tmunu)
    {
        std::string file_tmunu = "New_Tmunu_Evol_" + std::to_string( index ) + ".dat";
        save_Tmunu(Tmunu, folder_path + file_tmunu);
    }

};
void save_observables(const observables& obs, const std::string& filename)
{
    std::ofstream outfile(filename, std::ios::app); 

    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    const int fieldWidth = 17;

            
    outfile      
    << std::left << std::setw(fieldWidth) << obs.tau 
    << std::left << std::setw(fieldWidth) << obs.eps_tot 
    << std::left << std::setw(fieldWidth) << obs.eps_2_x
    << std::left << std::setw(fieldWidth) << obs.eps_2_y 
    << std::left << std::setw(fieldWidth) << obs.eps_p_x 
    << std::left << std::setw(fieldWidth) << obs.eps_p_y
    << std::left << std::setw(fieldWidth) << obs.u_perp
    << std::left << std::setw(fieldWidth) << obs.av_invReyn
    << std::left << std::setw(fieldWidth) << obs.Etr 
    << std::endl;

    outfile.close();
    if (!outfile) {
        std::cerr << "Error writing to file: " << filename << std::endl;
    }
};

void save_surface(const freeze_out_vars& surface, const std::string& filename)
{
    std::ofstream outfile(filename, std::ios::app); 

    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    const int fieldWidth = 17;

            
    outfile      
    << std::left << std::setw(fieldWidth) << surface.xmu[0]
    << std::left << std::setw(fieldWidth) << surface.xmu[1]
    << std::left << std::setw(fieldWidth) << surface.xmu[2]
    << std::left << std::setw(fieldWidth) << surface.sigma_mu[0]
    << std::left << std::setw(fieldWidth) << surface.sigma_mu[1]
    << std::left << std::setw(fieldWidth) << surface.sigma_mu[2]
    << std::left << std::setw(fieldWidth) << surface.vi[0]
    << std::left << std::setw(fieldWidth) << surface.vi[1]
    << std::left << std::setw(fieldWidth) << surface.pi_ij[0]
    << std::left << std::setw(fieldWidth) << surface.pi_ij[1]
    << std::left << std::setw(fieldWidth) << surface.pi_ij[2]
    << std::left << std::setw(fieldWidth) << surface.Pi_bulk 
    << std::endl;

    outfile.close();
    if (!outfile) {
        std::cerr << "Error writing to file: " << filename << std::endl;
    }
};

// using lattice = std::vector<std::vector<std::vector<double>>>;

class Timer {
    private:
        std::chrono::high_resolution_clock::time_point start_time;
        std::chrono::high_resolution_clock::time_point end_time;
        bool running = false;
    
    public:
        // Start the timer
        void start() {
            start_time = std::chrono::high_resolution_clock::now();
            running = true;
        }
    
        // Stop the timer and print the duration
        void stop() {
            if (!running) {
                std::cerr << "Timer was not started!" << std::endl;
                return;
            }
    
            end_time = std::chrono::high_resolution_clock::now();
            running = false;
    
            // double duration_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
            double duration_seconds = std::chrono::duration<double>(end_time - start_time).count();


            // parameters_class::shell_message shell_msg;
            // if (!shell_msg.quiet)
            //     std::cout << "Execution time: " << duration_seconds << " seconds" << std::endl;
        }
    };
class Administration// : public TObject
{
public:

    void save_hydro(const Hydro_grid& grid, const std::string& filename);
    void save_Tmunu(const Tmunu_grid& grid, const std::string& filename);
    void save_obs(const observables& obs, const std::string& filename);
    // void save_surface(const freeze_out_vars& surface, const std::string& filename);
    void delete_file_if_exists(const std::string& filename);
    void clear_file_contents(const std::string& filename);
    // void save_hydro_root(TFile);

    // ClassDef(Administration, 1);
};


#endif