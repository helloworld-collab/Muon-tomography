#include <iostream>
#include <cmath>
#include <fstream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"

void basic_config() {
    // Open the ROOT file and retrieve the TTree named "Event"
    TFile *f = new TFile("/home/samantha/Documents/corsika-78000/src/utils/coast/CorsikaRead/DAT000002.root");
    TTree *t1 = (TTree*)f->Get("Event");

    // Define variables for the branches.
    double Nmu;
    t1->SetBranchAddress("Nmu", &Nmu);

    // Create a histogram for the number of muons.
    // Adjust the number of bins and range to handle the larger values of muons in the events.
    TH1F *hMuon = new TH1F("hMuon", "Number of Muons;N_{muons};Counts", 50, 0, 2500);  // Increased the bin range

    Long64_t nentries = t1->GetEntries();

    for (Long64_t i = 0; i < nentries; i++) {
        t1->GetEntry(i);

        // Only fill the histogram if the value is non-negative and valid
        if (Nmu >= 0) {
            hMuon->Fill(Nmu);
        } else {
            std::cerr << "Warning: Invalid value for Nmu in event " << i + 1 << " (" << Nmu << ")\n";
        }
    }

    // Create a canvas to display the histogram.
    TCanvas *canvas = new TCanvas("canvas", "Muon Histogram", 800, 600);
    //canvas->SetLogy();  // Use logarithmic scale on the y-axis

    // Draw the histogram.
    hMuon->Draw();

    // Save the canvas to a file.
    canvas->SaveAs("MuonHisto20events.png");

    // Print the statistics (mean, std dev, etc.)
    std::cout << "Mean number of muons: " << hMuon->GetMean() << std::endl;
    std::cout << "Standard deviation: " << hMuon->GetStdDev() << std::endl;
    std::cout << "Number of entries: " << hMuon->GetEntries() << std::endl;
}

