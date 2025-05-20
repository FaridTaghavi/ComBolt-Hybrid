// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// error_handler.h
#ifndef error_handler_H
#define error_handler_H

#include <iostream>
#include <stdexcept>

inline void landau_matching_energy(const double& energy) {
    if (energy < 0) {
        throw std::runtime_error("Landau matching failed: No real and positive eigenvalue found.");
    }
}
inline void landau_matching_eigenvalues(int STATUS) {
    if (STATUS) {
        throw std::runtime_error("In GSL: Eigenvalue computation failed. ");
    }

}
#endif // error_handler_H
