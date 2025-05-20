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

#include "/home/farid/MyRepositories/KineticTheory/kineticTheory/src/particle_frzout.h"
#include "/home/farid/MyRepositories/KineticTheory/kineticTheory/src/particle_urqmd.h"

void read_particles() {
    // Load your dictionary library
    gSystem->Load("/home/farid/MyRepositories/KineticTheory/kineticTheory/build/src/libsrc.so");

    // Open ROOT file
    TFile *f = TFile::Open("/home/farid/MyRepositories/KineticTheory/evolution_root/ComBolt-events.root");
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
	double eta_cut = 0.5;
	double y_cut = 0.5;

	// Vector of partilce
	std::vector<double> pT_list_ch;
	std::vector<double> pT_list_pi_pm;
	std::vector<double> pT_list_K_pm;
	std::vector<double> pT_list_p_pBar;

    // Loop over entries, each corresponds to one oversample
    for (size_t i = 0; i < nevents; ++i) {
        tree->GetEntry(i);

        // std::cout << "\n At event " << eventID 
        //           << ", ComBolt event  " << ComBolt_eventID
        //           << ", particlization sample " << sampleID << "\n";

        // Show freeze-out particles
        // if (particles_frzout && !particles_frzout->empty()) {
        //     std::cout << "  [particles_frzout] count = " << particles_frzout->size() << "\n";
        //     for (size_t j = 0; j < particles_frzout->size(); ++j) {
        //         const auto& p = particles_frzout->at(j);
        //         std::cout << "    pid: " << p.pid
        //                   << " | x = (" << p.x << ", " << p.y << ", " << p.z << ", " << p.t
        //                   << ") | p = (" << p.px << ", " << p.py << ", " << p.pz << ", " << p.E << ")\n";
        //     }
        // }

        // Show URQMD output particles
        if (particles_urqmd && !particles_urqmd->empty()) {
            // std::cout << "  [particles_urqmd] count = " << particles_urqmd->size() << "\n";
            for (size_t j = 0; j < particles_urqmd->size(); ++j) {
                const auto& particle = particles_urqmd->at(j);
				// Apply cuts
				// if ( particle.charge == 0 || std::abs(particle.eta) > eta_cut  )
				if ( particle.charge == 0 || std::abs(particle.y_rap) > y_cut  )
					continue;
                // std::cout << "    pid: " << particle.pid
                //           << " | charge: " << particle.charge
                //           << " | pT: " << particle.pT << " | ET: " << particle.ET
                //           << " | mT: " << particle.mT
                //           << " | phi: " << particle.phi
                //           << " | y rapidity: " << particle.y_rap
                //           << " | eta: " << particle.eta << "\n";
           		pT_list_ch.push_back(particle.pT);
           		if ( std::abs(  particle.pid ) == 211 )   pT_list_pi_pm.push_back(particle.pT);  // pion+-
           		if ( std::abs(  particle.pid ) == 321 )   pT_list_K_pm.push_back(particle.pT);   // K+-
           		if ( std::abs(  particle.pid ) == 2212 )   pT_list_p_pBar.push_back(particle.pT);// proton antiproton 
		   	}
        }
    }
	// std::cout << "pt list size: " << pT_list.size() << "\n";
	// for (auto p : pT_list){
	// 	std::cout << "elements "  <<  p << "\n";
	// }
	// delete particles_frzout;
	
//-------------------------  Analysis ------------------------------------------

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
	Nevn *= 0.9; // A manual normalization  correction by hand
	ptHist_ch->Scale(1.0 / Nevn);
	ptHist_pion->Scale(1.0 / Nevn);
	ptHist_K->Scale(1.0 / Nevn);
	ptHist_proton->Scale(1.0 / Nevn);
	
	TCanvas *c = new TCanvas("c", "Custom Histogram", 800, 600);
    c->SetLogy(1);
   
    ptHist_pion->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    ptHist_pion->GetYaxis()->SetTitle("(1/N)  d^2N/(dp_{T}dy) [(GeV/c)^{-1}]");

	ptHist_pion->SetLineColor(kBlue);
	ptHist_pion->SetMarkerColor(kBlue);
    ptHist_pion->SetMarkerStyle(20);
    // ptHist_pion->SetMarkerSize(1.0);

	ptHist_K->SetLineColor(kRed);
	ptHist_K->SetMarkerColor(kRed);
    ptHist_K->SetMarkerStyle(21);
	// ptHist_K->SetMarkerSize(1.2);
	
	ptHist_proton->SetLineColor(kGreen+2);
	ptHist_proton->SetMarkerColor(kGreen+2);
    ptHist_proton->SetMarkerStyle(22);
	// ptHist_proton->SetMarkerSize(1.3);
    
	
	
	// Draw and update canvas
    // ptHist->Draw("E1");  // "E1" shows error bars with markers
    // gPad->Update();      // Ensure the canvas refreshes
    // c->Modified();
    // c->Update();         // Make sure canvas is updated

    // Optional: Set the axis range to ensure the plot is visible
	ptHist_pion->SetMinimum(0.5e-2);
	ptHist_pion->GetXaxis()->SetRangeUser(0, 3);
    ptHist_pion->Draw("E1");
    ptHist_K->Draw("E1 same");
    ptHist_proton->Draw("E1 same");
	
	g1->Draw("P same");
    g2->Draw("P same");
    g3->Draw("P same");	
	// gPad->Update();  // Update the canvas again

    // Run ROOT event loop
    // gApplication->Run(kTRUE);

	// 1) Create your legend (x1,y1,x2,y2 in NDC [0–1] coordinates)
	auto leg = new TLegend(0.60, 0.45, 1.5, 0.85);
	
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
	latex->DrawLatex(0.35, 0.78, "pp #sqrt{s} = 7 TeV");
	latex->SetTextSize(0.04);
	latex->DrawLatex(0.35, 0.71, " |y| < 0.5");

    c->SaveAs("pt_spectrum.png");
	// delete particles_urqmd;
	delete g1, g2, g3;
    f_ALICE->Close();
    f->Close();
}

