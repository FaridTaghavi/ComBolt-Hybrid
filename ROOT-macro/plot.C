
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TList.h>
#include <iostream>
#include "Math/IFunction.h"



//============================= Functions and structures ====================================


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
    double SheatTensorXX  {0.};
    double BulkTensorXX  {0.};
};

double BesselI(unsigned int n, double x ){
	return ROOT::Math::cyl_bessel_i(n, x);
};

double free_streaming_energy_density(double tau, double r,  double theta, double tau0, double R0, double eps0){
	
	return 

		(eps0*tau0*(std::pow(R0,2)*BesselI(1,(r*(tau - tau0))/std::pow(R0,2)) + std::sqrt(4*std::pow(r,2)*std::pow(tau - tau0,2)*std::pow(BesselI(0,(r*(tau - tau0))/std::pow(R0,2)),2) + 
	    4*r*std::pow(R0,2)*(-tau + tau0)*BesselI(0,(r*(tau - tau0))/std::pow(R0,2))*BesselI(1,(r*(tau - tau0))/std::pow(R0,2)) + 
	    (std::pow(R0,4) - 4*std::pow(r,2)*std::pow(tau - tau0,2))*std::pow(BesselI(1,(r*(tau - tau0))/std::pow(R0,2)),2))))/
	    (2.*std::exp((std::pow(r,2) + std::pow(tau - tau0,2))/(2.*std::pow(R0,2)))*r*tau*(tau - tau0))
	 ;

};



size_t nx {101}, ny {101};
size_t midnx {50}, midny {50};
size_t time_step_1 {15}, time_step_2 {500}, time_step_3 {1000};


size_t grid_index_poistion(size_t i_x, size_t i_y)
{
	size_t n_x {nx}, n_y {ny};
	/*  A condition to avoid overflow index */
	if (i_x >= n_x || i_y >= n_y ){
		std::cout << "Index overflow\n";
		return 0;
	};

	return i_x * n_y + i_y;
}
// ==================================================================================================================================

void readAllTreesFromRootFile(const char* fileName) {

    // Open the ROOT file in READ mode
    TFile* file = TFile::Open(fileName, "READ");
    
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return;
    }
	
	//== Defining the variables =================
	//
	//-------------------------------------------
	//
	std::vector< double > time, analy_free_stream_time,  tau4Over3Times_energy_dens, tau4Over3Times_energy_dens_anal_freeStream;
	std::array < std::vector< double >, 2 > energy_vs_r_at_t1, energy_vs_r_at_t2, energy_vs_r_at_t3, energy_vs_r_at_y1, energy_vs_r_at_y2;
	
	//==========================================
	//
	//------------------------------------------
	
    TList* keys = file->GetListOfKeys();
	size_t num_of_time_step = keys->GetSize();
	// num_of_time_step = 1020;
	// Loop over TTrees
	for (size_t time_step = 0; time_step < num_of_time_step;  time_step += 30	){
		std::string tree_name {"At_time_step_" + std::to_string(time_step)};
		TTree *tree = (TTree*)file->Get(tree_name.c_str());
		 if (!tree) {
		 	    std::cerr << "Error: TTree not found!" << std::endl;
		 		    return;
	    }

		hydrodynamic_field hyd;
		
		
		size_t branch_number { grid_index_poistion(midnx, midny)};
		std::string branch_name {"At_point_" + std::to_string(branch_number) };
    	tree->SetBranchAddress(branch_name.c_str(), &hyd);
    	
		
		// Loop over the entries in the TTree
		//
		//
		// 
		//
		// ---------------------------------
    	// Long64_t nEntries = tree->GetEntries();
    	// for (Long64_t entry = 0; entry < nEntries; ++entry) {
    	// 	tree->GetEntry(entry);
		// 	time.push_back(hyd.tau);
		// 	tau4Over3Times_energy_dens.push_back(10 * std::pow( hyd.tau ,4/3.) *  hyd.ed);


    	// }
		// size_t entry {grid_index_poistion(15, 15)   }; 
    	tree->GetEntry(0);
	 	time.push_back(hyd.tau);
	 	tau4Over3Times_energy_dens.push_back(10 * std::pow( hyd.tau ,4/3.) *  hyd.ed);
		//=================================
		//
		//
		//
		//
		//--------------------------------
		delete tree;
	}
	// Find cross-sectional plot at constant time -----
	//
	// ================================================
	std::string tree_name {"At_time_step_" + std::to_string(time_step_1)};
	TTree *tree = (TTree*)file->Get(tree_name.c_str());
    if (!tree) {
	    std::cerr << "Error: TTree not found!" << std::endl;
	    return;
    }
	
	hydrodynamic_field hyd;
	double tau_at_time_step;
	for (int ix = 0; ix < nx; ++ix){
		std::string branch_name {"At_point_" + std::to_string(grid_index_poistion(ix, midnx + 0 )  ) };
    	tree->SetBranchAddress(branch_name.c_str(), &hyd);
    	tree->GetEntry(0);
		tau_at_time_step = hyd.tau;
		energy_vs_r_at_t1[0].push_back( hyd.x );
		energy_vs_r_at_t1[1].push_back( hyd.ed);

		std::string branch_name2 {"At_point_" + std::to_string(grid_index_poistion(ix, midnx + 3  )  ) };
    	tree->SetBranchAddress(branch_name2.c_str(), &hyd);
    	tree->GetEntry(0);
		tau_at_time_step = hyd.tau;
		energy_vs_r_at_y1[0].push_back( hyd.x );
		energy_vs_r_at_y1[1].push_back( hyd.ed);

		std::string branch_name3 {"At_point_" + std::to_string(grid_index_poistion(ix, midnx + 7 )  ) };
    	tree->SetBranchAddress(branch_name3.c_str(), &hyd);
    	tree->GetEntry(0);
		tau_at_time_step = hyd.tau;
		energy_vs_r_at_y2[0].push_back( hyd.x );
		energy_vs_r_at_y2[1].push_back( hyd.ed);
	}

	delete tree;

	file->Close();
	delete file;




	// Analytical free-streaming ================
	for (int itime = 0; itime < 1000; ++itime){
		double tau0_fs = 0.1;
		double tau_fs = itime * 0.01 + tau0_fs + 1e-12; 

		analy_free_stream_time.push_back(tau_fs);
		double freeStreamED { free_streaming_energy_density(tau_fs, 0.0001, 0., tau0_fs, 1., 1.  )   };
		tau4Over3Times_energy_dens_anal_freeStream.push_back( 0.5 *  std::pow( tau_fs ,4/3.) *  freeStreamED  ); // Check the 0.5 prefactor! 
	
	}
	
	std::array < std::vector <double>, 2 > fs_analy_energy_vs_r_at_t1;
	for (int ix = 0; ix < 10000; ++ix){
		double tau0_fs = 0.1;
		double tau_fs = tau_at_time_step; 
		double x = ix * 0.001 + 1e-12; 

		fs_analy_energy_vs_r_at_t1[0].push_back( x );
		fs_analy_energy_vs_r_at_t1[1].push_back (0.1 * 0.5 *  free_streaming_energy_density(tau_fs, x, 0., tau0_fs, 1., 1.  )   );
	}
	
	// ==========================================

    std::cout << "Finished reading from " << fileName << std::endl;
	// Create a TGraph
    TGraph *graph = new TGraph(time.size(), time.data(), tau4Over3Times_energy_dens.data());
    graph->SetTitle(";#tau [fm]; #tau^{3/4} #epsilon(0,0)"); // Set the title and axis labels
	   
	// Set marker style and line style
    graph->SetMarkerStyle(21);
    graph->SetLineColor(kBlue);
    graph->SetLineWidth(4);

	TGraph *graph2 = new TGraph(analy_free_stream_time.size(), analy_free_stream_time.data(), tau4Over3Times_energy_dens_anal_freeStream.data());
    graph2->SetLineColor(kRed);
    graph2->SetMarkerStyle(22);
    graph2->SetLineWidth(4);
	graph2->SetLineStyle(2);								  



	graph->GetXaxis()->SetLimits(0, 4.);	
	graph->GetYaxis()->SetRangeUser(0, 0.051);	
    // Create a canvas to draw the graph
    TCanvas *c1 = new TCanvas("c1", "Energy density from Kinetic Theory", 800, 600);
    // c1->SetGrid(); // Optional: draw a grid on the canvas
	   
	c1->SetMargin(0.12, 0.05, 0.12, 0.05);
    graph->Draw("AL"); // "A" for axis, "P" for points, "L" for line
    graph2->Draw("L");	

    TLegend *legend = new TLegend(0.65, 0.7, 0.9, 0.9);
    // legend->SetHeader("Legend", "C"); // Set header
    legend->AddEntry(graph, "Simulation", "l"); // "l" for line, "p" for points
    legend->AddEntry(graph2, "Analytic", "l"); // "l" for line, "p" for points
	legend->SetBorderSize(0);
	legend->SetTextSize(0.04);
    legend->Draw(); // Draw the legend on the canvas


    c1->Modified();
    c1->Update();

	c1->SaveAs("plots/p1.pdf");
	// c1->SaveAs("plots/p1.eps");
	// c1->SaveAs("plots/p1.svg");

	// delete graph, graph2, c1;	
	
	// --------------------------------------------------------------------
	
	// Create a TGraph
    TGraph *graph3 = new TGraph(energy_vs_r_at_t1[0].size(), energy_vs_r_at_t1[0].data(), energy_vs_r_at_t1[1].data());
    graph3->SetTitle(";x [fm]; #epsilon(x,0)"); // Set the title and axis labels
	   
	// Set marker style and line style
    graph3->SetMarkerStyle(21);
    graph3->SetLineColor(kBlue);
    graph3->SetLineWidth(4);

	graph3->GetXaxis()->SetLimits(-5.3, 5.3);	
	// graph3->GetYaxis()->SetRangeUser(0, 0.00051);	



	TGraph *graph4 = new TGraph(fs_analy_energy_vs_r_at_t1[0].size(), fs_analy_energy_vs_r_at_t1[0].data(), fs_analy_energy_vs_r_at_t1[1].data());
    graph4->SetLineColor(kRed);
    graph4->SetMarkerStyle(22);
    graph4->SetLineWidth(4);
	graph4->SetLineStyle(2);								  
   
   
	TGraph *graph5 = new TGraph(energy_vs_r_at_y1[0].size(), energy_vs_r_at_y1[0].data(), energy_vs_r_at_y1[1].data());
    graph5->SetMarkerStyle(21);
    graph5->SetLineColor(kCyan);
    graph5->SetLineWidth(4);

	graph5->GetXaxis()->SetLimits(0, 5.3);	
	// graph5->GetYaxis()->SetRangeUser(0, 0.00051);	

	TGraph *graph7 = new TGraph(energy_vs_r_at_y2[0].size(), energy_vs_r_at_y2[0].data(), energy_vs_r_at_y2[1].data());
    graph7->SetMarkerStyle(21);
    graph7->SetLineColor(kGreen);
    graph7->SetLineWidth(4);

	graph5->GetXaxis()->SetLimits(0, 5.3);	
	// graph5->GetYaxis()->SetRangeUser(0, 0.00051);	
	
	
	// Create a canvas to draw the graph
    TCanvas *c2 = new TCanvas("c2", "Energy density from Kinetic Theory", 800, 600);
    // c2->SetGrid(); // Optional: draw a grid on the canvas
	   
	c2->SetMargin(0.12, 0.05, 0.12, 0.05);
    graph3->Draw("AL"); // "A" for axis, "P" for points, "L" for line
    // graph4->Draw("L");	
    // graph5->Draw("L");	
    // graph7->Draw("L");	

    TLegend *legend2 = new TLegend(0.65, 0.7, 0.9, 0.9);
    // legend->SetHeader("Legend", "C"); // Set header
    legend2->AddEntry(graph, "Simulation", "l"); // "l" for line, "p" for points
    legend2->AddEntry(graph2, "Analytic", "l"); // "l" for line, "p" for points
	legend2->SetBorderSize(0);
	legend2->SetTextSize(0.04);
    legend2->Draw(); // Draw the legend on the canvas
   
    std::stringstream stream;
	stream << std::fixed << std::setprecision(2) << tau_at_time_step;
	std::string textInsidePlot = "#tau =  " + stream.str();

   	TLatex latex;
    latex.SetTextSize(0.04);  // Set the size of the text
    latex.SetTextColor(kBlack); // Set the color of the text
	// std::string textInsidePlot = "#tau = " + std::to_string(tau_at_time_step);
    latex.DrawLatexNDC(0.71, 0.65, textInsidePlot.c_str()	); // (x, y, text)

    c2->Modified();
    c2->Update();

	c2->SaveAs("plots/p2.pdf");
	// c1->SaveAs("plots/p1.eps");
	// c1->SaveAs("plots/p1.svg");

	// delete graph, graph2, c1;	
//====================================================================================	
	
	// file->ls();


    // Get the list of keys (objects) in the file
    // TList* keys = file->GetListOfKeys();
	// size_t num_of_time_step = keys->GetSize();

	// std::cout << "How many trees? " << num_of_time_step << "\n";
/* 
    // for (int i = 0; i <  num_of_time_step; ++i) {
    for (int i = 0; i <  10; ++i) {
        TKey* key = static_cast<TKey*>(keys->At(i));
    	if (key->GetClassName() == TString("TTree")) {
    	    
			TTree* tree = static_cast<TTree*>(key->ReadObj());
			size_t nEntries = tree->GetEntries();
			
			std::cout << "how many points? " << nEntries << "\n";

			     hydrodynamic_field hyd;
   			// Loop over the entries
	        for (Long64_t j = 0; j < 1; ++j) {
				
				 std::string branch_name {"At_point_" + std::to_string(j)};
		   	     tree->SetBranchAddress(branch_name.c_str(), &hyd ); 
			// size_t pos_ind { grid_index_poistion(16, 16) };
			// for (size_t j = pos_ind; j < pos_ind + 1; ++j  ){
			}		
	   		// Load the entry
	        for (Long64_t j = 0; j < nEntries; ++j) {
	  	  	    tree->GetEntry(j);
	
	       		 //  Now 'event' contains the data for the current entry
	       		 std::cout << "Entry " << j << ": "
	       		   << "tau = " << hyd.tau << ", "
	       		   << "x = " << hyd.x << ", "
	       		   << "ed = " << hyd.ed << std::endl;
	        }
	        // Close the file

		delete tree;
		}
	}
	file->Close();
	delete file;
	*/
    // Loop through all keys
    // for (int i = 0; i < keys->GetSize(); ++i) {
    //     // Get each key
    //     TKey* key = static_cast<TKey*>(keys->At(i));
    //     // Check if the key corresponds to a TTree
    //     if (key->GetClassName() == TString("TTree")) {
    //         // Get the TTree
    //         TTree* tree = static_cast<TTree*>(key->ReadObj());
    //         std::cout << "Found TTree: " << tree->GetName() << std::endl;

    //         // Loop through branches
    //         TObjArray* branches = tree->GetListOfBranches();
    //         for (int j = 0; j < branches->GetSize(); ++j) {
    //             TBranch* branch = static_cast<TBranch*>(branches->At(j));
    //             std::cout << "  Branch: " << branch->GetName() << std::endl;
    //         }

    //         // Print the number of entries in the TTree
    //         std::cout << "  Number of entries: " << tree->GetEntries() << std::endl;
    //         
    //         // Optional: Reading entries
    //         // For demonstration, you can uncomment the following block to read and print entries
    //         /*
    //         for (int entry = 0; entry < tree->GetEntries(); ++entry) {
    //             tree->GetEntry(entry); // Read the entry
    //             // Print or process the entry as needed
    //         }
    //         */
    //     }
    // }

    // Close the file
    // file->Close();
    // delete file;

}
void plot() {
    readAllTreesFromRootFile("/home/farid/MyRepositories/KineticTheory/kineticTheory/evolution/HydroEvol_ROOT.root ");
}









