// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// This class is used to define a particle in a PyROOT

#ifndef particle_urqmd_H
#define particle_urqmd_H

#include "Rtypes.h"
#include "TObject.h"

class particle_urqmd  : public TObject {
public:
    Int_t pid;        // Particle ID
    Int_t charge;     // Charge of the particle

    // Kinematic quantities
    Float_t pT;       // Transverse momentum
    Float_t ET;       // Transverse energy
    Float_t mT;       // Transverse mass
    Float_t phi;      // Azimuthal angle
    Float_t y_rap;    // Rapidity
    Float_t eta;      // Pseudorapidity

    particle_urqmd() {}; // Default constructor
    // Constructor
    particle_urqmd(Int_t pid_, Int_t charge_, Float_t pT_, Float_t ET_, Float_t mT_,
                   Float_t phi_, Float_t y_rap_, Float_t eta_)
        : pid(pid_), charge(charge_), pT(pT_), ET(ET_), mT(mT_), phi(phi_), y_rap(y_rap_), eta(eta_) {}

    virtual ~particle_urqmd() {}

    ClassDef(particle_urqmd, 1);
};

#endif // particle_frzout_H
