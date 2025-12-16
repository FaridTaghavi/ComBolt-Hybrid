// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.

// The parametrization of lattice QCD EOS from Jonah E. Bernhard arXiv:1804.06469 [nucl-th].

// lattice_EOS.h
#ifndef lattice_EOS_H
#define lattice_EOS_H

#include <cmath>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

double Sech(double x) {
    return 1.0 / std::cosh(x);
}

// The input temperature is in 1/fm and the output energy density and pressure are in 1/fm^4.
inline void lattice_ed_pressure(const double &T, double &ed, double &P)
{
    double hbarC = 197.3; // MeV fm
        // Constants
    const double an = -8.7704;
    const double bn = 3.9200;
    const double cn = 0.0;
    const double dn = 0.3419;
    const double ad = -1.2600;
    const double bd = 0.8425;
    const double cd = 0.0;
    const double dd = -0.0475;
    const double ct = 3.8706;
    const double t0 = 0.9761;
    const double Tc = 154.0 / hbarC; // 154 MeV in 1/fm;
    const double pid = (95.0 * M_PI * M_PI) / 180.0;

    double t = T / Tc;
    P = ((pid + dn/std::pow(t,4) + cn/std::pow(t,3)
         + bn/std::pow(t,2) + an/t)*std::pow(T,4)*(1 + std::tanh(ct*(t - t0))))/(2.*(1 + dd/std::pow(t,4) + cd/std::pow(t,3) + bd/std::pow(t,2) + ad/t));

    double edMinus3P = (std::pow(T,4)*(ct*t*Tc*(dd*std::pow(Tc,4) + cd*t*std::pow(Tc,4) + bd*std::pow(t,2)*std::pow(Tc,4) + ad*std::pow(t,3)*std::pow(Tc,4) + std::pow(t,4)*std::pow(Tc,4))*
        (pid*std::pow(t,4)*std::pow(Tc,4) + Tc*(an*std::pow(t,3)*std::pow(Tc,3) + Tc*(dn*std::pow(Tc,2) + cn*t*std::pow(Tc,2) + bn*std::pow(t,2)*std::pow(Tc,2))))*std::pow(Sech(ct*(t - t0)),2) - 
       std::pow(Tc,2)*(dd*std::pow(Tc,4) + cd*t*std::pow(Tc,4) + bd*std::pow(t,2)*std::pow(Tc,4) + ad*std::pow(t,3)*std::pow(Tc,4) + std::pow(t,4)*std::pow(Tc,4))*
        (an*std::pow(t,3)*std::pow(Tc,3) + Tc*(2*bn*std::pow(t,2)*std::pow(Tc,2) + Tc*(4*dn*Tc + 3*cn*t*Tc)))*(1 + std::tanh(ct*(t - t0))) + 
       std::pow(Tc,2)*(ad*std::pow(t,3)*std::pow(Tc,3) + Tc*(2*bd*std::pow(t,2)*std::pow(Tc,2) + Tc*(4*dd*Tc + 3*cd*t*Tc)))*
        (pid*std::pow(t,4)*std::pow(Tc,4) + Tc*(an*std::pow(t,3)*std::pow(Tc,3) + Tc*(dn*std::pow(Tc,2) + cn*t*std::pow(Tc,2) + bn*std::pow(t,2)*std::pow(Tc,2))))*(1 + std::tanh(ct*(t - t0)))))/
   (2.*Tc*std::pow(dd*std::pow(Tc,4) + cd*t*std::pow(Tc,4) + bd*std::pow(t,2)*std::pow(Tc,4) + ad*std::pow(t,3)*std::pow(Tc,4) + std::pow(t,4)*std::pow(Tc,4),2)); 

    ed = 3 * P + edMinus3P;

};


// =============================
// Inverse mapping: ed → T (via spline)
// =============================
class LatticeEOSInverse {
private:
    gsl_spline* spline;
    gsl_interp_accel* acc;
    double log_ed_min;
    double log_ed_max;
public:
    LatticeEOSInverse(const latticeEOS_Table &table) {
        
        acc = gsl_interp_accel_alloc();
        spline = gsl_spline_alloc(gsl_interp_cspline, table.Ts.size());
        gsl_spline_init(spline, table.log_eds.data(), table.Ts.data(), table.Ts.size());

        log_ed_min = table.log_eds.front();
        log_ed_max = table.log_eds.back();
    }

    ~LatticeEOSInverse() {
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }

    void get_T_from_ed(const double &ed, double &T) const {
        double log_ed = std::log(ed);
        if (log_ed < log_ed_min) {
            T = 0.0;
            return;
        } else if (log_ed > log_ed_max)
        {
            std::cerr << "--> Warning: Energy density " << ed << " is out of interpolation range, " << log_ed_min << " < " << " log(ed) " << " < " << log_ed_max << ".\n";
            return;
        }
        
        T = gsl_spline_eval(spline, std::log(ed), acc);
    //  return std::pow( ed / 13.9997, 0.25); // for conformal EOS);
    }

};

#endif // lattice_EOS_H