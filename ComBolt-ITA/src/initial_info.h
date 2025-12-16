// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// This class is used to define a particle in a PyROOT

#ifndef initial_info_H
#define initial_info_H

#include "Rtypes.h"
#include "TObject.h"

// class initial_info : public TObject {
// public:
//     Float_t tau0, Etr, eps_tot;
//     Float_t eps_2_x, eps_2_y, Rsq;
// 
//     initial_info() {};
//     
//     initial_info(Float_t tau0_, Float_t Etr_, Float_t eps_tot_, Float_t eps_2_x_, Float_t eps_2_y_,
//         Float_t Rsq_)
//    : tau0(tau0_), Etr(Etr_), eps_tot(eps_tot_), eps_2_x(eps_2_x_), eps_2_y(eps_2_y_), Rsq(Rsq_) {}
// 
//     virtual ~initial_info() {};
// 
//     ClassDef(initial_info, 1);
// };


class initial_info : public TObject {
public:
    // Matches the text output column-by-column
    Float_t tau;         // obs.tau
    Float_t eps_tot;     // obs.eps_tot
    Float_t eps_2_x;     // obs.eps_2_x
    Float_t eps_2_y;     // obs.eps_2_y
    Float_t eps_3_x;     // obs.eps_3_x
    Float_t eps_3_y;     // obs.eps_3_y
    Float_t eps_p_x;     // obs.eps_p_x
    Float_t eps_p_y;     // obs.eps_p_y
    Float_t u_perp;      // obs.u_perp
    Float_t av_invReyn;  // obs.av_invReyn
    Float_t Etr;         // obs.Etr
    Float_t Rsq;         // obs.Rsq
    Float_t gamma_hat;   // obs.gamma_hat
    Float_t life_time;   // obs.life_time

    initial_info() {};

    initial_info(Float_t tau_, Float_t eps_tot_, Float_t eps_2_x_, Float_t eps_2_y_,
                 Float_t eps_3_x_, Float_t eps_3_y_, Float_t eps_p_x_, Float_t eps_p_y_,
                 Float_t u_perp_, Float_t av_invReyn_, Float_t Etr_, Float_t Rsq_,
                 Float_t gamma_hat_, Float_t life_time_)
        : tau(tau_), eps_tot(eps_tot_), eps_2_x(eps_2_x_), eps_2_y(eps_2_y_),
          eps_3_x(eps_3_x_), eps_3_y(eps_3_y_), eps_p_x(eps_p_x_), eps_p_y(eps_p_y_),
          u_perp(u_perp_), av_invReyn(av_invReyn_), Etr(Etr_), Rsq(Rsq_),
          gamma_hat(gamma_hat_), life_time(life_time_) {}

    virtual ~initial_info() {};

    ClassDef(initial_info, 1);
};


#endif // initial_info_H
