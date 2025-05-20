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
//#include <H5Cpp.h> 
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"

#include "ROOT_admin.h"


void ROOT_admin::save_hydro_in_ROOT(const Hydro_grid& grid, const std::string& filename, const int & time_step_index){

    TFile *outputFile = new TFile(filename.c_str(), "UPDATE");

    std::string tree_name {"At_time_step_" + std::to_string(time_step_index)};
    TTree *tree = new TTree(tree_name.c_str(), tree_name.c_str());

    Hydro_grid hyd_grd {grid};

    for (size_t i = 0; i < grid.size(); ++i) {
        std::string branch_name {"At_point_" + std::to_string(i)};
        tree->Branch(branch_name.c_str(), &hyd_grd[i], "tau/D:x/D:y/D:ed/D:Px/D:Py/D:Pz/D:utau/D:ux/D:uy/D:ueta/D:ShearTensorXX/D:BulkTensorXX/D");
    }
    
    tree->Fill();
    tree->Write("", TObject::kOverwrite); //The overwrite helps if I do not remove the old file, it overwrite the time step tree and keep the rest intact.
    delete tree;

    outputFile->Close();
    delete outputFile;

};