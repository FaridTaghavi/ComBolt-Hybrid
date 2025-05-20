

void read_particles() {
    // Load your shared library containing the dictionary
    // gSystem->Load("build/src/libsrc.so");
    // gSystem->AddIncludePath("-I/home/farid/MyRepositories/KineticTheory/kineticTheory/src");
    // #include "particle_frzout.h"
    #include "/home/farid/MyRepositories/KineticTheory/kineticTheory/src/particle_frzout.h"    
    #include "/home/farid/MyRepositories/KineticTheory/kineticTheory/src/particle_urqmd.h"    

    // Open the ROOT file
    TFile *file = TFile::Open("ComBolt-events.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file." << std::endl;
        return;
    }

    // Access the tree
    TTree *tree = (TTree*)file->Get("Events");
    if (!tree) {
        std::cerr << "Cannot find TTree 'Events'." << std::endl;
        return;
    }

    // Create vectors and set branch addresses
    std::vector<particle_frzout> *particles_frzout = nullptr;
    std::vector<particle_urqmd> *particles_urqmd = nullptr;

    tree->SetBranchAddress("particles_frzout", &particles_frzout);
    tree->SetBranchAddress("particles_urqmd", &particles_urqmd);

    Long64_t nentries = tree->GetEntries();
    std::cout << "Number of events: " << nentries << std::endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        std::cout << "Event " << i << std::endl;

        std::cout << "  [Frzout] Number of particles: " << particles_frzout->size() << std::endl;
        for (const auto& p : *particles_frzout) {
            std::cout << "    pid: " << p.pid << ", px: " << p.px << ", py: " << p.py << std::endl;
        }

        std::cout << "  [UrQMD] Number of particles: " << particles_urqmd->size() << std::endl;
        for (const auto& p : *particles_urqmd) {
            std::cout << "    pid: " << p.pid << ", px: " << p.px << ", py: " << p.py << std::endl;
        }

        std::cout << std::endl;
    }

    file->Close();
}

