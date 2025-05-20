
struct hydrodynamic_field
{
    double tau, x, y;
    double ed {0.};
    double Px {0.};
    double Py {0.};
    double Pz {0.};
    //u must consider with upper index
    double utau {0.};
    double ux {0.};
    double uy {0.};
    double ueta {0.};
    //Extend if needed
    // double SheatTensorXX  {0.};
    // double BulkTensoXX  {0.};
};

void plotDensity() {



    TFile *outputFile = new TFile("Hydrodynamic_Evolution.root", "RECREATE");

	size_t number_of_timeSteps = 2000;
	for (size_t time_step = 0; time_step < number_of_timeSteps; ++time_step){
		std::cout << "Time " << time_step << " is saved.\n";
    	std::ifstream inputFile("../evolution_gamma0.0/HydroEvol_" + std::to_string(time_step) + ".dat");
    	
    	if (!inputFile.is_open()) {
    	    std::cerr << "Error: Could not open the file!" << std::endl;
    	    return;
    	}
		std::string time_step_name {"time_step" + std::to_string(time_step) };
    	TTree *tree = new TTree(time_step_name.c_str(), time_step_name.c_str());

		std::vector<hydrodynamic_field> hyd_grid;
		std::string line;
		while (std::getline(inputFile, line)){
			std::istringstream iss(line);
			hydrodynamic_field hyd; 
    		if (iss >> hyd.tau >> hyd.x >> hyd.y >> hyd.ed >> hyd.Px >> hyd.Py >> hyd.Pz >> hyd.utau >> hyd.ux >> hyd.uy >> hyd.ueta ) {
				hyd_grid.push_back(hyd);
    		}
		}

    	for (size_t i = 0; i < hyd_grid.size(); ++i) {
		    std::string branchnumber {"At_point_" + std::to_string(i)};
			tree->Branch(branchnumber.c_str(), &hyd_grid[i]);
		 }

   		tree->Fill();
    	tree->Write();
    	delete tree;
		inputFile.close();
	}

    outputFile->Close();
	delete outputFile;

    std::cout << "Data successfully read from file and stored in TTree!" << std::endl;



    // Create a histogram with 100 bins, ranging from -5 to 5
    
	// TH1F *h = new TH1F("h", "Example Histogram;X axis;Y axis", 100, -5, 5);

    // // Fill the histogram with 10000 random numbers from a Gaussian distribution
    // for (int i = 0; i < 10000; ++i) {
    //     double x = gRandom->Gaus(0, 1); // Mean 0, standard deviation 1
    //     h->Fill(x);
    // }

    // // Draw the histogram
    // h->Draw();
}
