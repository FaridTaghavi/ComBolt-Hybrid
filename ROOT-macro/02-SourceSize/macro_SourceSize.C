#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <TDirectory.h>
#include <TKey.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TList.h>
#include <TROOT.h>
#include <algorithm>  // for std::upper_bound

#include "/home/farid/MyRepositories/KineticTheory/kineticTheory/src/particle_frzout.h"
#include "/home/farid/MyRepositories/KineticTheory/kineticTheory/src/particle_urqmd.h"


// std::string ROOT_file = "/home/farid/MyRepositories/KineticTheory/evolution_root/ComBolt-events_etaOS=0.08.root";
std::string ROOT_file = "/home/farid/MyRepositories/KineticTheory/evolution_root/ComBolt-events.root";
std::string library_src = "/home/farid/MyRepositories/KineticTheory/kineticTheory/build/src/libsrc.so";


int findBin(double B, const std::vector<double>& bin_edges) {
    auto it = std::upper_bound(bin_edges.begin(), bin_edges.end(), B);
    int idx = static_cast<int>(it - bin_edges.begin()) - 1;
    if (idx >= 0 && idx < static_cast<int>(bin_edges.size() - 1)) {
        return idx;
    }
    return -1;  // out of range
}
void read_particles() {
    // Load your dictionary library
    gSystem->Load(library_src.c_str() );

    // Open ROOT file
    TFile *f = TFile::Open(ROOT_file.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "Failed to open ROOT file.\n";
        return;
    }

    // Retrieve tree
    TTree *tree = (TTree*)f->Get("Events");
    if (!tree) {
        std::cerr << "Failed to get TTree 'Events'.\n";
        f->Close();
        return;
    }

    // Declare variables to read branches
    int eventID = -1, sampleID = -1, ComBolt_eventID = -1;
    std::vector<particle_frzout> *particles_frzout = nullptr;
    std::vector<particle_urqmd>  *particles_urqmd  = nullptr;

    // Set branch addresses
    tree->SetBranchAddress("eventID", &eventID);
    tree->SetBranchAddress("sampleID", &sampleID);
    tree->SetBranchAddress("ComBolt_eventID", &ComBolt_eventID);
    tree->SetBranchAddress("particles_frzout", &particles_frzout);
    tree->SetBranchAddress("particles_urqmd",  &particles_urqmd);

    size_t nevents = tree->GetEntries();
    std::cout << "Number of samples (tree entries): " << nevents << "\n";

	// Set cuts:
	double eta_cut = 0.8;
	double y_cut = 1.0;
	double pT_min_pi = 0.14, pT_max_pi = 4.0;

	// Average momentum and pair center frame
	struct TwoBodySource{
		double r2_star, M_T;
	};
	std::vector<TwoBodySource> S2_pi_pi, S2_pr_pr; 
	
	// pion mass
	const double pion_mass = 0.139;
	// proton mass
	const double proton_mass = 0.938;
	
	// >>>>>>>>>>>>>>>>>>>>>>>  Handel the bins for pions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	// std::vector<double> bin_edges_pi_pT = {0.15, 0.30, 0.50, 0.70, 0.90, 1.50};
	std::vector<double> bin_edges_pi_pT = {0.15, 0.30, 0.48, 0.68, 0.84, 1.30}; // ALICE paper
	// std::vector<double> bin_edges_pi_pT = {0.15,0.175, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0};
	// std::vector<double> bin_edges_pi_pT = {0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05};
	std::vector<double> bin_edges_pi_MT;
	for (double b : bin_edges_pi_pT){
		double mt =  std::sqrt(  pion_mass * pion_mass +  b * b ) ; 
		bin_edges_pi_MT.push_back( mt    );
		// std::cout << "mT  " << mt << "\n";
	}
	

	std::vector<double> bin_cent_pi;
	std::vector<double> bin_width_pi;
	for (int i = 0; i < bin_edges_pi_MT.size() - 1; ++i){
		double center =   ( bin_edges_pi_MT[i+1] + bin_edges_pi_MT[i] ) / 2;   
		double width =   ( bin_edges_pi_MT[i+1] - bin_edges_pi_MT[i] ) / 2;   
		bin_cent_pi.push_back( center   ); 
		bin_width_pi.push_back( width   ); 
		// std::cout << "mT  center " << center << "    " << width  << "\n";
	}

	// std::vector<double> bin_edges = {0.204, 0.33, 0.52, 0.71, 0.91, 1.51};
	// std::vector<double> bin_cent = {0.267, 0.425, 0.615, 0.81, 1.21};
	// std::vector<double> bin_width = {0.063 / 1, 0.095 / 1, 0.095 / 1, 0.1 / 1, 0.3 / 1};
	std::vector<double> rs_pipi_error;//  = {0.02, 0.02, 0.02, 0.02, 0.02};
	// std::vector<double> bin_edges = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};

	const size_t n_bins_pi = bin_edges_pi_MT.size() - 1;
    std::vector<double> sum_r2_pipi_star(n_bins_pi, 0.0);
    std::vector<int> count_pi(n_bins_pi, 0);
	std::vector<double> avg_r_pipi_source(n_bins_pi, 0.0);
	std::vector<double> var_r_pipi_source(n_bins_pi, 0.0);
	std::vector<std::vector< double>>  rs_pipi; 	
	// ======================================================================================
	

	// >>>>>>>>>>>>>>>>>>>>>>>  Handel the bins for protons <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	// std::vector<double> bin_edges_pro_pT = {0.52, .67, .70}; // ALICE paper
	std::vector<double> bin_edges_pro_MT = { 1.0774, 1.138, 1.1986, 1.2592, 1.37, 1.555, 1.82, 2.60    }  ; // ALICE paper
	// for (double b : bin_edges_pro_pT){
	// 	double mt =  std::sqrt(  proton_mass * proton_mass +  b * b ) ; 
	// 	bin_edges_pro_MT.push_back( mt    );
	// 	std::cout << "mT  " << mt << "\n";
	// }
	

	std::vector<double> bin_cent_pro;
	std::vector<double> bin_width_pro;
	for (int i = 0; i < bin_edges_pro_MT.size() - 1; ++i){
		double center =   ( bin_edges_pro_MT[i+1] + bin_edges_pro_MT[i] ) / 2;   
		double width =   ( bin_edges_pro_MT[i+1] - bin_edges_pro_MT[i] ) / 2;   
		bin_cent_pro.push_back( center   ); 
		bin_width_pro.push_back( width   ); 
		// std::cout << "mT  center " << center << "    " << width  << "\n";
	}

	// std::vector<double> bin_edges = {0.204, 0.33, 0.52, 0.71, 0.91, 1.51};
	// std::vector<double> bin_cent = {0.267, 0.425, 0.615, 0.81, 1.21};
	// std::vector<double> bin_width = {0.063 / 1, 0.095 / 1, 0.095 / 1, 0.1 / 1, 0.3 / 1};
	std::vector<double> rs_prpr_error;//  = {0.02, 0.02, 0.02, 0.02, 0.02};
	// std::vector<double> bin_edges = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};

	const size_t n_bins_pro = bin_edges_pro_MT.size() - 1;
    std::vector<double> sum_r2_prpr_star(n_bins_pro, 0.0);
    std::vector<int> count_pro(n_bins_pro, 0);
	std::vector<double> avg_r_prpr_source(n_bins_pro, 0.0);
	std::vector<double> var_r_prpr_source(n_bins_pro, 0.0);
	std::vector<std::vector< double>>  rs_prpr; 	
	// ======================================================================================



    // Loop over entries, each corresponds to one oversample
	for (size_t i = 0; i < nevents; ++i) {
        tree->GetEntry(i);

        if (particles_frzout && !particles_frzout->empty()) {
			
			
			// We check all pair of particles
            for (size_t j = 0; j < particles_frzout->size(); ++j) {
            	
				const auto& particle_1 = particles_frzout->at(j);
				double y1 = 0.5 * std::log(( particle_1.E + particle_1.pz) / (particle_1.E - particle_1.pz));
				
				double p1 = std::sqrt(particle_1.px * particle_1.px +
                     particle_1.py * particle_1.py +
                     particle_1.pz * particle_1.pz);

				double eta1 = 0.5 * std::log((p1 + particle_1.pz) / (p1 - particle_1.pz));

				double pT1 = std::sqrt( particle_1.px * particle_1.px +  particle_1.py * particle_1.py );

				// if (std::abs(y1)> y_cut || pT1 < pT_min_pi || pT1 > pT_max_pi    ) continue; // rapidity cut
				if (std::abs(eta1)> eta_cut || pT1 < pT_min_pi || pT1 > pT_max_pi    ) continue; // pseudorapidity cut
				
				for (size_t k = j + 1; k < particles_frzout->size(); ++k) {
            	    
					const auto& particle_2 = particles_frzout->at(k);
					double y2 = 0.5 * std::log(( particle_2.E + particle_2.pz) / (particle_2.E - particle_2.pz));
					double pT2 = std::sqrt( particle_2.px * particle_2.px +  particle_2.py * particle_2.py );

					double p2 = std::sqrt(particle_2.px * particle_2.px +
                	     particle_2.py * particle_2.py +
                	     particle_2.pz * particle_2.pz);

					double eta2 = 0.5 * std::log((p2 + particle_2.pz) / (p2 - particle_2.pz));
				
					// if (std::abs(y2) > y_cut || pT2 < pT_min_pi || pT2 > pT_max_pi   ) continue; // rapidity cut
					if (std::abs(eta2) > eta_cut || pT2 < pT_min_pi || pT2 > pT_max_pi   ) continue;// pseudorapidity cut
					
					TwoBodySource  TwoS_pipi, TwoS_prpr;
					
					// ============================= Pair of pions ================================
					if ( std::abs(  particle_1.pid) == 221 && std::abs( particle_2.pid  ) == 221  ){
						double xs =  (particle_2.x - particle_1.x);	
						double ys =  (particle_2.y - particle_1.y);	
						double zs =  (particle_2.z - particle_1.z);

				
						double k_T = 0.5 * std::sqrt( ( particle_2.px + particle_1.px  ) * (  particle_2.px + particle_1.px  ) +   ( particle_2.py + particle_1.py  ) * (  particle_2.py + particle_1.py )  ); 
						TwoS_pipi.M_T =  std::sqrt(pion_mass * pion_mass + k_T * k_T  );
						// double measure = 1.0 / TwoS.M_T ; 	// / TwoS.M_T is the trem comes from the measure; 
						
						
						TwoS_pipi.r2_star = ( xs * xs + ys * ys + zs * zs);
					    // measure /= TwoS.r2_star;

						// TwoS.r2_star = xs * xs + ys * ys; 	

						S2_pi_pi.push_back( TwoS_pipi);
					}
					
					// ============================= Pair of protons ================================
					if ( std::abs(  particle_1.pid) == 2212 && std::abs( particle_2.pid  ) == 2212  ){
						double xs =  (particle_2.x - particle_1.x);	
						double ys =  (particle_2.y - particle_1.y);	
						double zs =  (particle_2.z - particle_1.z);

				
						double k_T = 0.5 * std::sqrt( ( particle_2.px + particle_1.px  ) * (  particle_2.px + particle_1.px  ) +   ( particle_2.py + particle_1.py  ) * (  particle_2.py + particle_1.py )  ); 
						TwoS_prpr.M_T =  std::sqrt(proton_mass * proton_mass + k_T * k_T  );
						// double measure = 1.0 / TwoS.M_T ; 	// / TwoS.M_T is the trem comes from the measure; 
						
						
						TwoS_prpr.r2_star = ( xs * xs + ys * ys + zs * zs);
					    // measure /= TwoS.r2_star;

						// TwoS.r2_star = xs * xs + ys * ys; 	

						S2_pr_pr.push_back( TwoS_prpr);
					}
					// Apply cuts
					// if ( particle.pid != 221)
					// 	continue;
            	    
            	    // std::cout << "    pid: " << particle_1.pid
            	    //           << " | x = (" << particle_1.x << ", " << particle_1.y << ", " << particle_1.z << ", " << particle_1.t
            	    //           << ") | p = (" << particle_1.px << ", " << particle_1.py << ", " << particle_1.pz << ", " << particle_1.E << ")\n";
           			
					// pT_list_ch.push_back(particle.pT);
           			// if ( std::abs(  particle.pid ) == 211 )   pT_list_pi_pm.push_back(particle.pT);  // pion+-
           			// if ( std::abs(  particle.pid ) == 321 )   pT_list_K_pm.push_back(particle.pT);   // K+-
           			// if ( std::abs(  particle.pid ) == 2212 )   pT_list_p_pBar.push_back(particle.pT);// proton antiproton 
		   		}
			}
        }
		// std::cout << "numner of pion pairs:  " << S2_pi_pi.size() << "\n";

	// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Pion Pion souece size <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<	
		// Loop over two particle source and bin
    	for (const auto& dp : S2_pi_pi) {
    	    int bin = findBin(dp.M_T, bin_edges_pi_MT);
    	    if (bin != -1) {
    	        sum_r2_pipi_star[bin] +=  dp.r2_star / dp.M_T; // 1 / dp.M_T is the trem comes from the measure; 
    	        count_pi[bin] += 1;
    	    }
    	}
		
		std::vector<double> rs_pipi_event;
		for (size_t i = 0; i < n_bins_pi; ++i) {
    	    if (count_pi[i] > 0) {
				rs_pipi_event.push_back( std::sqrt(  sum_r2_pipi_star[i] / count_pi[i] / 6.0  ) ); // 1/6 comes from the standard definition of r0.
    	    } else {
		    	rs_pipi_event.push_back( 0.0);
    	    }
    	}
		rs_pipi.push_back(rs_pipi_event);
		S2_pi_pi.clear();	
		
	// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Proton Proton souece size <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<	
		// Loop over two particle source and bin
    	for (const auto& dp : S2_pr_pr) {
    	    int bin = findBin(dp.M_T, bin_edges_pro_MT);
    	    if (bin != -1) {
    	        sum_r2_prpr_star[bin] +=  dp.r2_star / dp.M_T; // 1 / dp.M_T is the trem comes from the measure; 
    	        count_pro[bin] += 1;
    	    }
    	}
		
		std::vector<double> rs_prpr_event;
		for (size_t i = 0; i < n_bins_pro; ++i) {
    	    if (count_pro[i] > 0) {
				rs_prpr_event.push_back( std::sqrt(  sum_r2_prpr_star[i] / count_pro[i] / 6.0  ) ); // 1/6 comes from the standard definition of r0.
    	    } else {
		    	rs_prpr_event.push_back( 0.0);
    	    }
    	}
		// std::cout << "proton?   " << rs_prpr_event[0] << "\n";
		rs_prpr.push_back(rs_prpr_event);
		S2_pr_pr.clear();	
    }

    // Find average of r pion pion  source at each mT bin

	for (size_t i = 0; i < n_bins_pi; ++i) {
		size_t ent = 0;
		for (auto r  : rs_pipi){
			avg_r_pipi_source[i] += r[i];
			ent++;
		}
		avg_r_pipi_source[i] /= ent;
	}
	
	std::cout << "Pion Pion \n"; 	
	for (size_t i = 0; i < n_bins_pi; ++i) {
		size_t ent = 0;
		for (auto r  : rs_pipi){
			double diff = r[i] - avg_r_pipi_source[i];
        	var_r_pipi_source[i] += diff * diff;
			ent++;
		}
		var_r_pipi_source[i] /= ent;
		
		double stddev = std::sqrt(var_r_pipi_source[i]);      // standard deviation
    	double stderr = stddev / std::sqrt(ent);               // standard error of the mean
    	
		rs_pipi_error.push_back( stderr);
		std::cout << "Bin " << i << ": m_T = " << bin_cent_pi[i]
        	      << ", mean = " << avg_r_pipi_source[i]
            	  << ", stderr = " << stderr << '\n';
	}



		

    // Find average of r proton proton  source at each mT bin

	for (size_t i = 0; i < n_bins_pro; ++i) {
		size_t ent = 0;
		for (auto r  : rs_prpr){
			avg_r_prpr_source[i] += r[i];
			ent++;
		}
		avg_r_prpr_source[i] /= ent;
	}

	std::cout << "Proton Proton \n"; 	
	for (size_t i = 0; i < n_bins_pro; ++i) {
		size_t ent = 0;
		for (auto r  : rs_prpr){
			double diff = r[i] - avg_r_prpr_source[i];
        	var_r_prpr_source[i] += diff * diff;
			ent++;
		}
		var_r_prpr_source[i] /= ent;
		
		double stddev = std::sqrt(var_r_prpr_source[i]);      // standard deviation
    	double stderr = stddev / std::sqrt(ent);               // standard error of the mean
    	
		rs_prpr_error.push_back( stderr);
		std::cout << "Bin " << i << ": m_T = " << bin_cent_pro[i]
        	      << ", mean = " << avg_r_prpr_source[i]
            	  << ", stderr = " << stderr << '\n';
	}
	// std::cout << "pt list size: " << pT_list.size() << "\n";
	// for (auto p : pT_list){
	// 	std::cout << "elements "  <<  p << "\n";
	// }
	// delete particles_frzout;
	
//-------------------------  Analysis ------------------------------------------

	
	
    
	
	TCanvas* c1 = new TCanvas("c1", "Source size", 600, 600);
    // TGraph* graph = new TGraph(n_bins, bin_centers.data(), avg_rs.data());
	TGraphErrors* graph = new TGraphErrors(n_bins_pi, bin_cent_pi.data(), avg_r_pipi_source.data(), bin_width_pi.data(), rs_pipi_error.data());
    graph->SetTitle(";#LT m_{T} #GT [GeV / #it{c}^{2}] ; #it{r}_{core} [fm]");
    
	   
	// graph->SetFillColor(2);
    // graph->SetFillStyle(3001);
    // graph->Draw("2A");
    // graph->Draw("P same");
	
	graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.2);
	graph->SetFillColorAlpha(kRed, 0.3);  // semi-transparent error boxes
	graph->SetMarkerStyle(20);
	graph->SetMarkerSize(1);
	graph->SetLineColor(kBlue);
    graph->Draw("2A");
    graph->Draw("P same");
	graph->GetXaxis()->SetLimits(0.2, 2.5);
	graph->GetYaxis()->SetRangeUser(0.5, 3.0);  // Set Y range
	
	
	
	// --- Add second graph ---
	TGraphErrors* graph2 = new TGraphErrors(
	    n_bins_pro,
	    bin_cent_pro.data(),
	    avg_r_prpr_source.data(),
	    bin_width_pro.data(),
	    rs_prpr_error.data()
	);
	graph2->SetMarkerStyle(21);
	graph2->SetMarkerSize(1.2);
	graph2->SetFillColorAlpha(kGreen + 2, 0.3);
	graph2->SetLineColor(kGreen + 2);
	graph2->Draw("2 same");
	graph2->Draw("P same");
	
	
	
	
	gPad->Update();
    c1->Update();
	c1->SaveAs("source_size.png");

/*
//>>>>>>>>>>>>>>>>   Read data from HEPDATA <<<<<<<<<<<<<<<<

 
    // 1) Open your HEPData ROOT file
    TFile *f_ALICE = TFile::Open("/home/farid/MyRepositories/KineticTheory/evolution_root/HEPData-ins1357424-v1-Table_1.root"); // replace with your filename
    if (!f_ALICE || f_ALICE->IsZombie()) {
        std::cerr << "Error: Cannot open HEPData ROOT file\n";
        return;
    }

    // 2) Navigate to the “Table 1;1” directory
    //    HEPData often encodes semicolons in the name—TFile::GetDirectory handles it.
    TDirectory *dir = f_ALICE->GetDirectory("Table 1");
    if (!dir) {
        std::cerr << "Error: Directory 'Table 1' not found\n";
        f_ALICE->Close();
        return;
    }

    // 3) Retrieve the graphs by their base names (omit the ";1" cycle)
    TGraph *g1 = (TGraph*)dir->Get("Graph1D_y1");
    TGraph *g2 = (TGraph*)dir->Get("Graph1D_y2");
    TGraph *g3 = (TGraph*)dir->Get("Graph1D_y3");

    if (!g1 || !g2 || !g3) {
      std::cerr << "One of the graphs was not found\n";
      f_ALICE->Close();
      return;
    }

    g1->SetLineColor(kBlue);    g1->SetMarkerColor(kBlue);    g1->SetMarkerStyle(24);
    g2->SetLineColor(kRed);     g2->SetMarkerColor(kRed);     g2->SetMarkerStyle(25);
    g3->SetLineColor(kGreen+2); g3->SetMarkerColor(kGreen+2); g3->SetMarkerStyle(26);

	// >>>>>>>>>>>>>>>>>>>>  Make histogram from ComBolt-ITA <<<<<<<<<<<<<<<<<<<<<<<<<< 


	double pT_bins[] = {0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, \
									0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.1, 1.2, 1.3, \
									1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, \
									2.8, 2.9, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.2, 4.4, 4.6, 4.8, 5.}; 
	int npTbins = sizeof(pT_bins) / sizeof(double) - 1;
	
	TH1F *ptHist_ch = new TH1F("ptHist_ch", "pT spectrum of charged particles", npTbins, pT_bins);
	TH1F *ptHist_pion = new TH1F("ptHist_pion", "pT spectrum of pions", npTbins, pT_bins);
	TH1F *ptHist_K = new TH1F("ptHist_K", "pT spectrum of K", npTbins, pT_bins);
	TH1F *ptHist_proton = new TH1F("ptHist_proton", "pT spectrum of proton", npTbins, pT_bins);

	ptHist_ch   ->SetStats(0);
	ptHist_pion ->SetStats(0);
    ptHist_K    ->SetStats(0);
    ptHist_proton->SetStats(0);
	
	ptHist_ch   ->SetTitle(0);
	ptHist_pion ->SetTitle(0);
    ptHist_K    ->SetTitle(0);
    ptHist_proton->SetTitle(0);

	ptHist_ch->Sumw2();                  
	ptHist_pion->Sumw2();                 
	ptHist_K->Sumw2();                  
	ptHist_proton->Sumw2();           
	
	// pT_list = {0.15, 0.22, 0.35, 0.45, 0.8, 1.2, 2.5, 3.3};  // Example
	for (auto p : pT_list_ch)     ptHist_ch->Fill(p);
	for (auto p : pT_list_pi_pm)  ptHist_pion->Fill(p);
	for (auto p : pT_list_K_pm)   ptHist_K->Fill(p);
	for (auto p : pT_list_p_pBar) ptHist_proton->Fill(p);
	
	
	// for (auto p : pT_list_K_pm){
	// 	// std::cout << "-->  " << p << "\n";
	// 	ptHist->Fill(p);
	// }
	ptHist_ch->Scale(1.0, "width");      // divide each bin by its width
	ptHist_pion->Scale(1.0, "width");    
	ptHist_K->Scale(1.0, "width");     
	ptHist_proton->Scale(1.0, "width");

	double Nevn = tree->GetEntries();
	Nevn *= 1.0; // A manual normalization  correction by hand
				 
	ptHist_ch->Scale(1.0 / Nevn);
	ptHist_pion->Scale(1.0 / Nevn);
	ptHist_K->Scale(1.0 / Nevn);
	ptHist_proton->Scale(1.0 / Nevn);
	
	// >>>>>>>>>>>>>>>>>>>>  Calculate the integrated dN/dy >>>>>>>>>>>>>>>>>>>> 
	double dNdy_pion = 0., dNdy_K = 0., dNdy_proton = 0.;
	double dNdy_pion_err2 = 0., dNdy_K_err2 = 0., dNdy_proton_err2 = 0.; 
	
	for (int i = 1; i <= ptHist_pion->GetNbinsX(); ++i) {
	    double binContent = ptHist_pion->GetBinContent(i);
	    double binError   = ptHist_pion->GetBinError(i);
	    double binWidth   = ptHist_pion->GetBinWidth(i);
	
	    dNdy_pion += binContent * binWidth;
	    dNdy_pion_err2 += std::pow(binError * binWidth, 2);
	}

	double dNdy_pion_err = std::sqrt(dNdy_pion_err2);
	
	for (int i = 1; i <= ptHist_K->GetNbinsX(); ++i) {
	    double binContent = ptHist_K->GetBinContent(i);
	    double binError   = ptHist_K->GetBinError(i);
	    double binWidth   = ptHist_K->GetBinWidth(i);
	
	    dNdy_K += binContent * binWidth;
	    dNdy_K_err2 += std::pow(binError * binWidth, 2);
	}

	double dNdy_K_err = std::sqrt(dNdy_K_err2);
	
	for (int i = 1; i <= ptHist_proton->GetNbinsX(); ++i) {
	    double binContent = ptHist_proton->GetBinContent(i);
	    double binError   = ptHist_proton->GetBinError(i);
	    double binWidth   = ptHist_proton->GetBinWidth(i);
	
	    dNdy_proton += binContent * binWidth;
	    dNdy_proton_err2 += std::pow(binError * binWidth, 2);
	}

	double dNdy_proton_err = std::sqrt(dNdy_proton_err2);
	
	std::cout << "The dN/dy of pion:   " << dNdy_pion << " +/- " << dNdy_pion_err  << "\n";	
	std::cout << "The dN/dy of K:      " << dNdy_K << " +/- " << dNdy_K_err  << "\n";	
	std::cout << "The dN/dy of proton: " << dNdy_proton << " +/- " << dNdy_proton_err  << "\n";	
	
	// =========================================================================

	TCanvas *c = new TCanvas("c", "Custom Histogram", 800, 600);
    c->SetLogy(1);
   
    ptHist_pion->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    ptHist_pion->GetYaxis()->SetTitle("(1/N)  d^{2}N/(dp_{T}dy) [(GeV/c)^{-1}]");

	ptHist_pion->SetLineColor(kBlue + 2);
	ptHist_pion->SetMarkerColor(kBlue + 2);
    ptHist_pion->SetMarkerStyle(20);
    ptHist_pion->SetMarkerSize(1.3);

	ptHist_K->SetLineColor(kRed + 2);
	ptHist_K->SetMarkerColor(kRed +2 );
    ptHist_K->SetMarkerStyle(21);
	ptHist_K->SetMarkerSize(1.3);
	
	ptHist_proton->SetLineColor(kGreen + 4);
	ptHist_proton->SetMarkerColor(kGreen + 4);
    ptHist_proton->SetMarkerStyle(22);
	ptHist_proton->SetMarkerSize(1.3);
    
	
	
	// Draw and update canvas
    // ptHist->Draw("E1");  // "E1" shows error bars with markers
    // gPad->Update();      // Ensure the canvas refreshes
    // c->Modified();
    // c->Update();         // Make sure canvas is updated

    // Optional: Set the axis range to ensure the plot is visible
	ptHist_pion->SetMinimum(0.5e-2);
	ptHist_pion->GetXaxis()->SetRangeUser(0, 3);
	
	// ptHist_pion->SetFillColorAlpha(kBlue-2, 0.2); 
	// ptHist_pion->SetFillStyle(1001); 
    ptHist_pion->Draw("E");
	
	// ptHist_K->SetFillColorAlpha(kRed-2, 0.2); 
	// ptHist_K->SetFillStyle(1001); 
    ptHist_K->Draw("E same");
	
	// ptHist_proton->SetFillColorAlpha(kGreen, 0.2); 
	// ptHist_proton->SetFillStyle(1001); 
    ptHist_proton->Draw("E same");
	
	g1->Draw("P same");
    g2->Draw("P same");
    g3->Draw("P same");	
	// gPad->Update();  // Update the canvas again

    // Run ROOT event loop
    // gApplication->Run(kTRUE);

	// 1) Create your legend (x1,y1,x2,y2 in NDC [0–1] coordinates)
	auto leg = new TLegend(0.60, 0.65, 1.5, 0.85);
	
	// 2) Layout
	leg->SetNColumns(2);          // e.g. 3 columns (one per species pair)
	leg->SetTextFont(42);         // Helvetica
	leg->SetTextSize(0.03);       // relative to pad height
	leg->SetColumnSeparation(0.02);
	leg->SetEntrySeparation(0.01);
	
	// 3) Style
	leg->SetBorderSize(0);        // no border
	leg->SetFillStyle(0);         // transparent
	leg->SetMargin(0.2);          // padding inside the box
	
	// 4) Add entries in the desired order
	// leg->AddEntry(ptHist_pion,    "ComBolt-ITA #pi^{#pm}",   "lep");
	// leg->AddEntry(g1,      "ALICE #pi^{#pm}",   "lep");
	// leg->AddEntry(ptHist_K,  "ComBolt-ITA K^{#pm}",     "lep");
	// leg->AddEntry(g2,      "ALICE K^{#pm}",     "lep");
	// leg->AddEntry(ptHist_proton,     "ComBolt-ITA p-#bar{p}",     "lep");
	// leg->AddEntry(g3,      "ALICE p--#bar{p}",     "lep");
	leg->AddEntry(ptHist_pion,    "",   "lep");
	leg->AddEntry(g1,      "#pi^{+} + #pi^{-}",   "lep");
	leg->AddEntry(ptHist_K,  "",     "lep");
	leg->AddEntry(g2,      "K^{+} + K^{-}",     "lep");
	leg->AddEntry(ptHist_proton,     "",     "lep");
	leg->AddEntry(g3,      "p + #bar{p}",     "lep");
	
	// 5) Draw it *after* all histos/graphs
	leg->Draw();

	TLatex *latex = new TLatex();
	latex->SetNDC();               // use normalized coordinates (0–1)
	latex->SetTextSize(0.04);      // adjust text size
								   //
	// latex->SetTextColor(kRed+2);
	latex->SetTextFont(42);         
	latex->SetTextAngle(0);         // 0 = horizontal

	latex->SetTextSize(0.03);
	latex->DrawLatex(0.725, 0.85, "ALICE");
	latex->SetTextSize(0.03);
	latex->DrawLatex(0.58, 0.85, " CoMBolt-ITA");
	latex->SetTextSize(0.04);
	latex->DrawLatex(0.20, 0.23, "pp #sqrt{s} = 7 TeV");
	latex->SetTextSize(0.04);
	latex->DrawLatex(0.20, 0.16, " |y| < 0.5");

    c->SaveAs("pt_spectrum.png");
	// delete particles_urqmd;
	delete g1, g2, g3;
    f_ALICE->Close();
  */
  	f->Close();
}

