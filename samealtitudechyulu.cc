#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLegend.h"

void samealtitudechyulu() {
    // Open file
    TFile *f = new TFile("DAT000006samealtchyulu.root");

    // Get the event-level tree
    TTree *t1 = (TTree*)f->Get("Event");
    
    // Connect the branches with their member variables.
    Double_t Phi, Theta, Ntot, Nmu, Eprim;
    t1->SetBranchAddress("Theta", &Theta);
    t1->SetBranchAddress("Phi", &Phi);
    t1->SetBranchAddress("Ntot", &Ntot);
    t1->SetBranchAddress("Nmu", &Nmu);
    t1->SetBranchAddress("Eprim", &Eprim);
    
    // Define event-level histograms (cos(Theta) instead of Theta)
    TH1F* hcosthe = new TH1F("hcosthe", "cos(Theta)", 100, 0, 1);
    TH1F* hntot = new TH1F("hntot", "Ntot", 100, 0, 30);
    TH1F* hnmu = new TH1F("hnmu", "Nmu", 100, 0, 30);
    TH1F* heprim = new TH1F("heprim", "Eprim", 100, 0, 30);

    // Read all entries 
    Int_t nentries = t1->GetEntries();
    for (Int_t i = 0; i < nentries; i++) {
        t1->GetEntry(i);
        hcosthe->Fill(std::cos(Theta));
        hntot->Fill(Ntot);
        hnmu->Fill(Nmu);
        heprim->Fill(Eprim);
    }

    // Create a canvas for event-level plots
    TCanvas *cevent = new TCanvas("cevent", "Event level distributions samealtitudechyulu", 800, 600);
    cevent->Divide(2,2);
    cevent->cd(1); 
    hcosthe->GetXaxis()->SetTitle("cos(Theta)");
    hcosthe->GetXaxis()->CenterTitle();
    hcosthe->GetYaxis()->SetTitle("Counts");
    hcosthe->GetYaxis()->CenterTitle();
    hcosthe->SetLineColor(kGreen+2);
    hcosthe->Draw();
    
    cevent->cd(2); 
    hntot->GetXaxis()->SetTitle("Total Number of Particles");
    hntot->GetXaxis()->CenterTitle();
    hntot->GetYaxis()->SetTitle("Counts");
    hntot->GetYaxis()->CenterTitle();
    hntot->SetLineColor(kRed+2);
    hntot->Draw();
    
    cevent->cd(3); 
    hnmu->GetXaxis()->SetTitle("Number of Muons");
    hnmu->GetXaxis()->CenterTitle();
    hnmu->GetYaxis()->SetTitle("Counts");
    hnmu->GetYaxis()->CenterTitle();
    hnmu->SetLineColor(kBlue+2);
    hnmu->Draw();

    cevent->cd(4);  
    heprim->GetXaxis()->SetTitle("Primary Energy (GeV)");
    heprim->GetXaxis()->CenterTitle();
    heprim->GetYaxis()->SetTitle("Counts");
    heprim->GetYaxis()->CenterTitle();
    heprim->SetLineColor(kMagenta+2);
    heprim->Draw();
    cevent->SaveAs("event_distributions_samealtitudechyulu.jpeg");

    // Particle-level variables
    TTree *d;
    Double_t id, theta, ek;

    // Histogram for muon cos(theta) and energy
    TH1F* hmu_costheta = new TH1F("hmu_costheta", "Muon cos(Theta) Distribution samealtitudechyulu", 100, 0, 1);
    TH1F* hmu_energy = new TH1F("hmu_energy", "Muon Energy Distribution samealtitudechyulu", 100, 0, 30); // energy [GeV]

    // Loop over the tuples (one per shower with particles)
    Int_t nshowers = t1->GetEntries("Ntot>0");
    for (Int_t j = 1; j <= nshowers; j++) {
        std::ostringstream tname;
        tname << "data_1;" << j;
        d = (TTree*)f->Get(tname.str().c_str());

        if (!d) continue;

        d->SetBranchAddress("id", &id);
        d->SetBranchAddress("theta", &theta); 
        d->SetBranchAddress("ek", &ek); 
        
        Int_t npart = d->GetEntries();
        for (Int_t i = 0; i < npart; i++) {
            d->GetEntry(i);
            if (id == 5 || id == 6) {
                hmu_costheta->Fill(std::cos(theta));
                hmu_energy->Fill(ek);
            }
        }
    }

    std::cout << "Number of muons: " << hmu_costheta->GetEntries() << std::endl;

    // Draw muon particle distributions
    TCanvas *cparticle = new TCanvas("cparticle", "Particle level distributions samealtitudechyulu", 800, 600);
    cparticle->Divide(1,2);

    cparticle->cd(1);
    hmu_costheta->GetXaxis()->SetTitle("cos(Theta)");
    hmu_costheta->GetXaxis()->CenterTitle();
    hmu_costheta->GetYaxis()->SetTitle("Counts");
    hmu_costheta->GetYaxis()->CenterTitle();
    hmu_costheta->SetLineColor(kGreen+2);
    hmu_costheta->Draw();

    cparticle->cd(2);
    hmu_energy->GetXaxis()->SetTitle("Energy (GeV)");
    hmu_energy->GetXaxis()->CenterTitle();
    hmu_energy->GetYaxis()->SetTitle("Counts");
    hmu_energy->GetYaxis()->CenterTitle();
    hmu_energy->SetLineColor(kBlue+1);
    hmu_energy->Draw();

    cparticle->SaveAs("particle_distributions_samealtitudechyulu.jpeg");
}

