// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

  
//#include <H5Cpp.h>  
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>
//#include <H5Cpp.h> 

#include "Administration.h"

void Administration::save_hydro(const Hydro_grid& grid, const std::string& filename)
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
void Administration::save_Tmunu(const Tmunu_grid& grid, const std::string& filename)
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

void Administration::save_obs(const observables& obs, const std::string& filename)
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
    << std::left << std::setw(fieldWidth) << obs.R 
    << std::endl;

    outfile.close();
    if (!outfile) {
        std::cerr << "Error writing to file: " << filename << std::endl;
    }
};
void Administration::delete_file_if_exists(const std::string& filename) {
    namespace fs = std::filesystem;

    if (fs::exists(filename)) {
        if (fs::remove(filename)) {
            std::cout << "File " << filename << " deleted successfully." << std::endl;
        } else {
            std::cerr << "Failed to delete file: " << filename << std::endl;
        }
    } else {
        std::cout << "File " << filename << " does not exist." << std::endl;
    }
};
void Administration::clear_file_contents(const std::string& filename) {
    // Open the file in truncation mode to clear its contents
    std::ofstream outfile(filename, std::ios::trunc);

    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
    } else {
        std::cout << "File contents cleared: " << filename << std::endl;
    }

    // Closing the file automatically happens when `outfile` goes out of scope
}