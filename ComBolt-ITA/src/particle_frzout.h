// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// This class is used to define a particle in a PyROOT

#ifndef particle_frzout_H
#define particle_frzout_H

#include "Rtypes.h"
#include "TObject.h"

class particle_frzout : public TObject {
public:
    Int_t pid;
    Float_t t, x, y, z;
    Float_t E, px, py, pz;

    particle_frzout() {};
    
    particle_frzout(Int_t pid_, Float_t t_, Float_t x_, Float_t y_, Float_t z_,
        Float_t E_, Float_t px_, Float_t py_, Float_t pz_)
   : pid(pid_), t(t_), x(x_), y(y_), z(z_),
     E(E_), px(px_), py(py_), pz(pz_) {}

    virtual ~particle_frzout() {};

    ClassDef(particle_frzout, 1);
};

#endif // particle_frzout_H
