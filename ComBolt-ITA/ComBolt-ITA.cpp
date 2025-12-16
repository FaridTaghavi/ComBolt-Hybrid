// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.


#include <boost/program_options.hpp>
#include "Administration.h"
#include "BoltzEQ_solver.h"
#include "parameters_class.h"

using namespace boost::program_options;
int main(int argc, const char* argv[]) {

    try {
        namespace po = boost::program_options;
        
        po::options_description desc{"combolt-ITA"};
        desc.add_options()
            ("output,o", po::value<std::string>(), "Folder path to save the output files")
            ("surface_output", po::value<std::string>(), "Folder path to save freeze-out surface files")
            ("quiet", po::value<bool>(), "No evolution progress")
            ("trento_file", po::value<std::string>(), "Trento file path")
            ("AMPT_file", po::value<std::string>(), "AMPT file path")
            ("smear_sigmaV", po::value<double>(), "AMPT smearing partons in spatail dimensions")
            ("smear_lambda", po::value<double>(), "AMPT smearing vz")
            ("smear_sigmaPhi", po::value<double>(), "AMPT smearing phi_p")
            ("AMPT_norm", po::value<double>(), "AMPT normalization factor")
            ("AMPT_etas_cut", po::value<double>(), "AMPT pseudorapidity cut")
            ("save_BD", po::value<bool>(), "Save the Boltzmann distribution")
            ("save_hydro", po::value<bool>(), "Save the hydrodynamic fields")
            ("save_Tmunu", po::value<bool>(), "Save the energy-momentum tensor")
            ("save_observables", po::value<bool>(), "Save the observables")
            ("save_init_info", po::value<bool>(), "Save the initial conditions information")
            ("nx", po::value<size_t>(), "Number of grid points in x direction")
            ("xmax", po::value<double>(), "Maximum value of x")
            ("nvp", po::value<size_t>(), "Number of grid points in v_p direction")
            ("nvphi", po::value<size_t>(), "Number of grid points in v_phi direction")
            ("Z2_symmetry", po::value<bool>(), "Z2 symmetry")
            ("h", po::value<double>(), "Time step")
            ("alpha_h", po::value<double>(), "Time step stiffness")
            ("tau_max", po::value<double>(), "Maximum value of tau")
            ("tau0", po::value<double>(), "Initial time")
            ("momentum_isotropicity", po::value<double>(), "Momentum isotropicity")
            ("initialization_mode", po::value<int>(), "initialization_mode")
            ("energy_tot_tau0", po::value<double>(), "Total energy at tau0")
            ("R0", po::value<double>(), "R0")
            ("ellipticity", po::value<double>(), "Ellipticity")
            ("grid_max", po::value<double>(), "Grid max")
            ("grid_step", po::value<double>(), "Grid step")
            ("norm", po::value<double>(), "Norm")
            ("insert_threshold", po::value<double>(), "Insert threshold")
            ("number_of_divistion", po::value<unsigned int>(), "Number of division")
            ("remove_threshold", po::value<double>(), "Remove threshold")
            ("number_of_cores_BoltzEq_solver", po::value<unsigned int>(), "Number of cores for Boltzmann equation solver")
            ("number_of_cores_integration", po::value<unsigned int>(), "Number of cores for integration")
            ("number_of_cores_interpolation", po::value<unsigned int>(), "Number of cores for interpolation")
            ("number_of_cores_calc_obs", po::value<unsigned int>(), "Number of cores for calculating observables")
            ("dissipation_mode", po::value<unsigned int>(), "Dissipation mode")
            ("gamma_hat", po::value<double>(), "Gamma hat")
            ("eta_over_s_min", po::value<double>(), "Eta over s minimum")
            ("eta_over_s_slope", po::value<double>(), "Eta over s slope")
            ("eta_over_s_pow", po::value<double>(), "Eta over s power")
            ("eta_over_s_Tc", po::value<double>(), "Eta over s Tc")
            ("C0", po::value<double>(), "conformal equation of state: energy = C0 * T^4")
            ("epsFO", po::value<double>(), "freeze-out energy density")
            ("TFO", po::value<double>(), "freeze-out temperature")
            ("EOS_mode", po::value<unsigned int>(), "Equation of state mode")
            ("calculate_freeze_out", po::value<bool>(), "Calculate freeze-out")
            ("end_evol_last_frozen_cell", po::value<bool>(), "End evolution at last frozen cell")
            ("help,h", "Help screen");

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        // Extract values from `vm` with default fallbacks
        parameters_class::saving_info save_info{
            vm.count("output") ? vm["output"].as<std::string>() : "../../evolution/",
            vm.count("surface_output") ? vm["surface_output"].as<std::string>() : "../../evolution/",
            vm.count("save_BD") ? vm["save_BD"].as<bool>() : false,
            vm.count("save_hydro") ? vm["save_hydro"].as<bool>() : true,
            vm.count("save_Tmunu") ? vm["save_Tmunu"].as<bool>() : false,
            vm.count("save_observables") ? vm["save_observables"].as<bool>() : true,
            vm.count("save_init_info") ? vm["save_init_info"].as<bool>() : true
        };

        parameters_class::lattice_info lattice_info{
            vm.count("nx") ? vm["nx"].as<size_t>() : 101,
            vm.count("xmax") ? vm["xmax"].as<double>() : 19.5,
            vm.count("nvp") ? vm["nvp"].as<size_t>() : 15,
            vm.count("nvphi") ? vm["nvphi"].as<size_t>() : 40,
            vm.count("Z2_symmetry") ? vm["Z2_symmetry"].as<bool>() : true
        };

        parameters_class::numerical_calc_param numerical_param{
            vm.count("h") ? vm["h"].as<double>() : 0.0025,
            vm.count("tau_max") ? vm["tau_max"].as<double>() : 15.5,
            vm.count("alpha_h") ? vm["alpha_h"].as<double>() : 0.02
        };

        parameters_class::initial_values initial_vals{
            vm.count("tau0") ? vm["tau0"].as<double>() : 5e-5,
            vm.count("initialization_mode") ?  vm["initialization_mode"].as<int>() : 0 ,
            vm.count("momentum_isotropicity") ? vm["momentum_isotropicity"].as<double>() : 0.05,
            (vm.count("energy_tot_tau0") ? vm["energy_tot_tau0"].as<double>() : 1280 ) / parameters_class::hbarC,
            vm.count("R0") ? vm["R0"].as<double>() : 1.97,
            vm.count("ellipticity") ? vm["ellipticity"].as<double>() : 0.416
        };
        
        parameters_class::medium_info medium{
            vm.count("C0") ? vm["C0"].as<double>() : 13.8997, 
           ( vm.count("epsFO") ? vm["epsFO"].as<double>() : 0.30 ) / parameters_class::hbarC,
           ( vm.count("TFO") ? vm["TFO"].as<double>() : 0.155 ) / parameters_class::hbarC,
            vm.count("EOS_mode") ? vm["EOS_mode"].as<unsigned int>() : 0
        };
        
        // parameters_class::dissipation_info dissipation{
        //     vm.count("dissipation_mode") ? vm["dissipation_mode"].as<unsigned int>() : 0, 
        //     vm.count("eta_over_s") ? vm["eta_over_s"].as<double>() : 0.1 / (4. * M_PI),
        //     vm.count("gamma_hat") ? vm["gamma_hat"].as<double>() : 20.0
        // };

        parameters_class::dissipation_info dissipation{
            vm.count("dissipation_mode") ? vm["dissipation_mode"].as<unsigned int>() : 0, 
            vm.count("eta_over_s_min") ? vm["eta_over_s_min"].as<double>() : 0.1 / (4. * M_PI),
            (vm.count("eta_over_s_slope") ? vm["eta_over_s_slope"].as<double>() : 0.0) * parameters_class::hbarC,
            vm.count("eta_over_s_pow") ? vm["eta_over_s_pow"].as<double>() : 0.0,
            (vm.count("eta_over_s_Tc") ? vm["eta_over_s_Tc"].as<double>() : 0.155) / parameters_class::hbarC,
            vm.count("gamma_hat") ? vm["gamma_hat"].as<double>() : 20.0
        };

        parameters_class::trento_info trento{
            vm.count("grid_max") ? vm["grid_max"].as<double>() : 15.0,
            vm.count("grid_step") ? vm["grid_step"].as<double>() : 0.2,
            (vm.count("norm") ? vm["norm"].as<double>() : 15.0 )/ parameters_class::hbarC,
            vm.count("trento_file") ? vm["trento_file"].as<std::string>() : "../../trento-master/events/event_set6/8.dat"
        };

        parameters_class::AMPT_info ampt{
            vm.count("smear_sigmaV") ? vm["smear_sigmaV"].as<double>() : 0.05,
            vm.count("smear_lambda") ? vm["smear_lambda"].as<double>() : 0.5,
            vm.count("smear_sigmaPhi") ? vm["smear_sigmaPhi"].as<double>() : 0.5,
            vm.count("AMPT_norm") ? vm["AMPT_norm"].as<double>() : 1.0,
            vm.count("AMPT_etas_cut") ? vm["AMPT_etas_cut"].as<double>() : 0.5,
            vm.count("AMPT_file") ? vm["AMPT_file"].as<std::string>() : "../../kineticTheory/initial_state_AMPT/parton-initial-afterPropagation8.dat"
        };

        parameters_class::interpolation_parameters interpolation{
            vm.count("insert_threshold") ? vm["insert_threshold"].as<double>() : 0.2,
            vm.count("number_of_divistion") ? vm["number_of_divistion"].as<unsigned int>() : 2,
            vm.count("remove_threshold") ? vm["remove_threshold"].as<double>() : 5e-7
        };

        parameters_class::multicore_parameters multicore{
            vm.count("number_of_cores_BoltzEq_solver") ? vm["number_of_cores_BoltzEq_solver"].as<unsigned int>() : 13,
            vm.count("number_of_cores_integration") ? vm["number_of_cores_integration"].as<unsigned int>() : 13,
            vm.count("number_of_cores_interpolation") ? vm["number_of_cores_interpolation"].as<unsigned int>() : 13,
            vm.count("number_of_cores_calc_obs") ? vm["number_of_cores_calc_obs"].as<unsigned int>() : 13
        };
        parameters_class::freeze_out_info freeze_out{
            vm.count("calculate_freeze_out") ? vm["calculate_freeze_out"].as<bool>() : true,
            vm.count("end_evol_last_frozen_cell") ? vm["end_evol_last_frozen_cell"].as<bool>() : true
        };
        parameters_class::shell_message shell_message{
            vm.count("quiet") ? vm["quiet"].as<bool>() : false
        };
        parameters_class::lattice_info_massive lattice_info_massive{
            vm.count("nx") ? vm["nx"].as<size_t>() : 101,
            vm.count("xmax") ? vm["xmax"].as<double>() : 19.5,
            vm.count("npT") ? vm["npT"].as<size_t>() : 15,
            vm.count("nphi_p") ? vm["nphi_p"].as<size_t>() : 15,
            vm.count("npz") ? vm["npz"].as<size_t>() : 15,
            vm.count("pT_max") ? vm["pT_max"].as<double>() : 10.0,
            vm.count("pz_max") ? vm["pz_max"].as<double>() : 10.0,
            vm.count("Z2_symmetry") ? vm["Z2_symmetry"].as<bool>() : true
        };

        parameters_class params(
            save_info, lattice_info, lattice_info_massive, numerical_param, initial_vals,
            medium, dissipation, trento, ampt, interpolation, multicore,
            freeze_out, shell_message
        );
        
        std::cout << "--> Parameters initialized successfully!\n";
    
        solve_BoltzmannEQ solve(params);
        solve();

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }    

       
    return 0;
};
