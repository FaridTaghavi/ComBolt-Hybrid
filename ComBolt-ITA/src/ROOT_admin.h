// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

#ifndef ROOT_ADMIN_H
#define ROOT_ADMIN_H


#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <functional>
#include <filesystem>

#include "grid.h"
#include "TObject.h"

namespace fs = std::filesystem;

// Enum to specify operation type
// enum Operation { ADD, SUBTRACT, MULTIPLE, DIVIDE };
// 
// struct path_struct {
//     fs::path save_path = fs::current_path();
//     fs::path load_path = fs::current_path();
// };
// 
// using lattice = std::vector<std::vector<std::vector<double>>>;


class ROOT_admin : public TObject
{
public:
    ROOT_admin(){
        
    };
    // Destructor
    ~ROOT_admin() {

    };

    void save_hydro_in_ROOT(const Hydro_grid& grid, const std::string& filename, const int & time_step_index);

  
    ClassDef(ROOT_admin, 1);
};





#endif