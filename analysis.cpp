#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include <iomanip>
#include <vector>
#include <iostream>

// Cosmetics
void set_style()
{
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(57);
    gStyle->SetOptStat(11220);
    gStyle->SetOptFit(111);
}

// Histogram cosmetic
void histo_cosmetic(TH1F *h)
{
    int fillColor = 426;
    h->SetFillColor(fillColor);
    h->SetLineWidth(1);
    h->SetLineColor(kBlack);
    h->GetYaxis()->SetTitle("Entries");
    h->GetYaxis()->SetTitleOffset(1.5);
    h->DrawCopy("H");
    h->DrawCopy("E,SAME");
}

// Function cosmetic
void function_cosmetic(TF1 *f)
{
    f->SetLineStyle(2);
    f->SetLineWidth(3);
    f->SetLineColor(kRed);
}

void analysis()
{
    // Cosmetics
    set_style();

    // Reading from simulation.root and creating analysis.root
    TFile *file1 = new TFile("simulation.root", "READ");
    TFile *file2 = new TFile("analysis.root", "RECREATE");
    file1->ls();

    // Copying histograms
    TH1F *h0 = (TH1F *)file1->Get("h0");
    TH1F *h1 = (TH1F *)file1->Get("h1");
    TH1F *h2 = (TH1F *)file1->Get("h2");
    TH1F *h3 = (TH1F *)file1->Get("h3");
    TH1F *h4 = (TH1F *)file1->Get("h4");
    TH1F *h5 = (TH1F *)file1->Get("h5");
    TH1F *h6 = (TH1F *)file1->Get("h6");
    TH1F *h7 = (TH1F *)file1->Get("h7");
    TH1F *h8 = (TH1F *)file1->Get("h8");
    TH1F *h9 = (TH1F *)file1->Get("h9");
    TH1F *h10 = (TH1F *)file1->Get("h10");
    TH1F *h11 = (TH1F *)file1->Get("h11");
    TH3F *h12 = (TH3F *)file1->Get("h12");

    // Creating canvas
    TCanvas *c = new TCanvas("c", "Various histograms", 200, 10, 950, 700);
    c->Divide(2, 2);
    TCanvas *cEnergy = new TCanvas("cEnergy", "Particle energy", 200, 10, 950, 700);
    TCanvas *cImpulse = new TCanvas("cImpulse", "Impulse and transverse impulse", 200, 10, 950, 700);
    cImpulse->Divide(2, 1);
    TCanvas *c3DAngles = new TCanvas("c3DAngles", "3D angles distribution", 200, 10, 950, 700);
    TCanvas *cInvMass = new TCanvas("cInvMass", "Invariant mass", 200, 10, 950, 700);
    cInvMass->Divide(3, 2);
    TCanvas *cResonance = new TCanvas("cResonance", "Resonance signal", 200, 10, 950, 700);
    cResonance->Divide(2, 1);
    TCanvas *cvuota = new TCanvas("cvuota", "cvuota signal", 20, 10, 90, 50);

    // Creating functions for fitting
    TF1 *f1 = new TF1("f1", "expo", 0, 10);
    TF1 *f4 = new TF1("f4", "pol0", 0, 2 * TMath::Pi());
    TF1 *f5 = new TF1("f5", "pol0", 0, TMath::Pi());
    TF1 *f11 = new TF1("f11", "gaus", 0.5, 1.3);

    // Creating histograms for resonance signal
    TH1F *hDiff1 = new TH1F(*h6);
    TH1F *hDiff2 = new TH1F(*h8);

    hDiff1->Add(h6, h7, 1, -1);
    hDiff2->Add(h8, h9, 1, -1);

    // Creating functions for resonance signal histograms fitting
    TF1 *fDiff1 = new TF1("fDiff1", "gaus", 0, 3.5);
    TF1 *fDiff2 = new TF1("fDiff2", "gaus", 0, 3.5);

    // Setting functions parameters
    double k_mass = 0.89166;
    double k_width = 0.050;
    f1->SetParameter(0, -1);
    f4->SetParameter(0, 1E5);
    f5->SetParameter(0, 1E5);
    f11->SetParameter(k_mass, k_width);
    fDiff1->SetParameters(k_mass, k_width);
    fDiff2->SetParameters(k_mass, k_width);

    // Setting functions cosmetic
    function_cosmetic(f1);
    function_cosmetic(f4);
    function_cosmetic(f5);
    function_cosmetic(f11);
    function_cosmetic(fDiff1);
    function_cosmetic(fDiff2);

    // Fitting histograms
    h1->Fit("f1");
    h4->Fit("f4");
    h5->Fit("f5");
    h11->Fit("f11");
    hDiff1->Fit("fDiff1");
    hDiff2->Fit("fDiff2");

    vector<TH1F *> histo{h0, h1, h4, h5, h3, h2, h6, h7, h8, h9, h10, h11, hDiff1, hDiff2};
    vector<const char *> bin_names{"#pi+", "#pi-", "K+", "K-", "p+", "p-", "K*"};

    //  Drawing on c
    for (int i = 0; i < 4; ++i)
    {
        c->cd(i + 1);
        if (i == 0)
        {
            histo[i]->GetXaxis()->SetTitle("Particle");
            for (int j = 1; j < 8; ++j)
            {
                histo[i]->GetXaxis()->SetBinLabel(j, bin_names[j - 1]);
            }
            histo[i]->SetMarkerStyle(20);
            histo[i]->SetMarkerSize(0.3);
        }
        if (i == 1)
        {
            histo[i]->GetXaxis()->SetTitle("Impulse modulus (GeV/c)");
        }
        if (i == 2)
        {
            histo[i]->GetXaxis()->SetTitle("#theta (rad)");
            histo[i]->SetMinimum(8000);
            histo[i]->SetMaximum(12000);
        }
        if (i == 3)
        {
            histo[i]->GetXaxis()->SetTitle("#phi (rad)");
            histo[i]->SetMinimum(8000);
            histo[i]->SetMaximum(12000);
        }
        histo_cosmetic(histo[i]);
    }

    // Drawing on cEnergy
    cEnergy->cd();
    histo[4]->GetXaxis()->SetTitle("Energy (GeV)");
    histo_cosmetic(histo[4]);

    // Drawing on cImpulse
    cImpulse->cd(1);
    histo[1]->GetXaxis()->SetTitle("Impulse modulus (GeV/c)");
    histo_cosmetic(histo[1]);

    cImpulse->cd(2);
    histo[5]->GetXaxis()->SetTitle("Transverse impulse modulus (GeV/c)");
    histo_cosmetic(histo[5]);

    // Drawing on c3DAngles
    c3DAngles->cd();
    h12->GetXaxis()->SetTitle("sin(#theta)*cos(#phi)");
    h12->GetYaxis()->SetTitle("sin(#theta)*sin(#phi)");
    h12->GetZaxis()->SetTitle("cos(#theta)");
    h12->GetXaxis()->SetTitleOffset(2);
    h12->GetYaxis()->SetTitleOffset(2);
    h12->GetZaxis()->SetTitleOffset(1);
    int fillColor = 426;
    h12->SetMarkerColor(fillColor);
    h12->SetLineWidth(2);
    h12->DrawCopy("H");
    h12->DrawCopy("E, SAME");

    // Drawing on cInvMass
    for (int i = 6; i < 12; ++i)
    {
        cInvMass->cd(i - 5);
        histo[i]->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
        histo_cosmetic(histo[i]);
    }

    // Drawing on cResonance
    for (int i = 12; i < 14; ++i)
    {
        cResonance->cd(i - 11);
        histo[i]->SetMinimum(-6000);
        histo[i]->SetMaximum(16000);
        histo[i]->SetAxisRange(0., 2.5, "X");
        histo[i]->GetYaxis()->SetTitleOffset(2.5);
        if (i == 12)
        {
            histo[i]->SetTitle("Resonant signal (all particles)");
        }
        if (i == 13)
        {
            histo[i]->SetTitle("Resonant signal (#pi/K)");
        }
        histo_cosmetic(histo[i]);
    }

    // Creating legends
    TLegend *leg1 = new TLegend(.2, .7, .6, .9);
    TLegend *leg4 = new TLegend(.1, .7, .5, .9);
    TLegend *leg5 = new TLegend(.1, .7, .5, .9);
    TLegend *leg11 = new TLegend(.1, .7, .4, .9);
    TLegend *legDiff1 = new TLegend(.1, .7, .3, .9);
    TLegend *legDiff2 = new TLegend(.1, .7, .3, .9);

    // Drawing legends
    leg1->SetFillColor(0);
    leg1->AddEntry(h1, "Impulse distribution");
    leg1->AddEntry(f1, "Exponential distribution");
    c->cd(2);
    leg1->Draw("SAME");

    leg4->SetFillColor(0);
    leg4->AddEntry(h4, "#theta distribution");
    leg4->AddEntry(f4, "Uniform distribution");
    c->cd(3);
    leg4->Draw("SAME");

    leg5->SetFillColor(0);
    leg5->AddEntry(h5, "#phi distribution");
    leg5->AddEntry(f5, "Uniform distribution");
    c->cd(4);
    leg5->Draw("SAME");

    leg11->SetFillColor(0);
    leg11->AddEntry(h11, "Invariant mass distribution");
    leg11->AddEntry(f11, "Gaussian distribution");
    cInvMass->cd(6);
    leg11->Draw("SAME");

    legDiff1->SetFillColor(0);
    legDiff1->AddEntry(h1, "Resonance signal distribution");
    legDiff1->AddEntry(f1, "Gaussian distribution");
    cResonance->cd(1);
    legDiff1->Draw("SAME");

    legDiff2->SetFillColor(0);
    legDiff2->AddEntry(h1, "Resonance signal distribution");
    legDiff2->AddEntry(f1, "Gaussian distribution");
    cResonance->cd(2);
    legDiff2->Draw("SAME");

    // Writing statistics on shell
    std::cout << "------------------------------------\n";
    std::cout << "Number of \u03C0+ generated: " << left << setw(10) << h0->GetBinContent(1) << " ± " << h0->GetBinError(1)
              << " percentage: " << (h0->GetBinContent(1) / h0->GetEntries()) * 100 << "%\n";
    std::cout << "Number of \u03C0- generated: " << left << setw(10) << h0->GetBinContent(2) << " ± " << h0->GetBinError(2)
              << " percentage: " << (h0->GetBinContent(2) / h0->GetEntries()) * 100 << "%\n";
    for (int i{3}; i < 8; ++i)
    {
        std::cout << "Number of " << bin_names[i - 1] << left << setw(10)
                  << " generated: " << h0->GetBinContent(i) << " ± " << h0->GetBinError(i)
                  << " percentage: " << (h0->GetBinContent(i) / h0->GetEntries()) * 100 << "%\n";
    }
    std::cout << "------------------------------------\n";
    std::cout << "Azimuthal angle distribution statistics:\n"
              << "Parameter = " << f4->GetParameter(0) << " ± " << f4->GetParError(0) << '\n'
              << "Chisquare = " << f4->GetChisquare() << '\n'
              << "DOF = " << f4->GetNDF() << '\n'
              << "Chisquare/DOF = " << f4->GetChisquare() / f4->GetNDF() << '\n'
              << "Probability = " << f4->GetProb() << '\n';
    std::cout << "------------------------------------\n";
    std::cout << "Polar angle distribution statistics:\n"
              << "Parameter = " << f5->GetParameter(0) << " ± " << f5->GetParError(0) << '\n'
              << "Chisquare = " << f5->GetChisquare() << '\n'
              << "DOF = " << f5->GetNDF() << '\n'
              << "Chisquare/DOF = " << f5->GetChisquare() / f5->GetNDF() << '\n'
              << "Probability = " << f5->GetProb() << '\n';
    std::cout << "------------------------------------\n";
    std::cout << "Impulse modulus distribution statistics:\n"
              << "First parameter = " << f1->GetParameter(0) << " ± " << f1->GetParError(0) << '\n'
              << "Second parameter = " << f1->GetParameter(1) << " ± " << f1->GetParError(1) << '\n'
              << "Chisquare = " << f1->GetChisquare() << '\n'
              << "DOF = " << f1->GetNDF() << '\n'
              << "Chisquare/DOF = " << f1->GetChisquare() / f1->GetNDF() << '\n'
              << "Probability = " << f1->GetProb() << '\n';
    std::cout << "------------------------------------\n";
    std::cout << "Invariant mass from daugther particles statistics:\n"
              << "Mean = " << f11->GetParameter(1) << " ± " << f11->GetParError(1) << '\n'
              << "Sigma = " << f11->GetParameter(2) << " ± " << f11->GetParError(2) << '\n'
              << "Amplitude = " << f11->GetParameter(0) << " ± " << f11->GetParError(0) << '\n'
              << "Chisquare/DOF = " << f11->GetChisquare() / f11->GetNDF() << '\n'
              << "Probability = " << f11->GetProb() << '\n';
    std::cout << "------------------------------------\n";
    std::cout << "Resonance signal (all particles) statistics:\n"
              << "Mean = " << fDiff1->GetParameter(1) << " ± " << fDiff1->GetParError(1) << '\n'
              << "Sigma = " << fDiff1->GetParameter(2) << " ± " << fDiff1->GetParError(2) << '\n'
              << "Amplitude = " << fDiff1->GetParameter(0) << " ± " << fDiff1->GetParError(0) << '\n' // Non so cosa sia!!!!!
              << "Chisquare/DOF = " << fDiff1->GetChisquare() / fDiff1->GetNDF() << '\n'
              << "Probability = " << fDiff1->GetProb() << '\n';
    std::cout << "------------------------------------\n";
    std::cout << "Resonance signal (\u03C0/K) statistics:\n"
              << "Mean = " << fDiff2->GetParameter(1) << " ± " << fDiff2->GetParError(1) << '\n'
              << "Sigma = " << fDiff2->GetParameter(2) << " ± " << fDiff2->GetParError(2) << '\n'
              << "Amplitude = " << fDiff2->GetParameter(0) << " ± " << fDiff2->GetParError(0) << '\n'
              << "Chisquare/DOF = " << fDiff2->GetChisquare() / fDiff2->GetNDF() << '\n'
              << "Probability = " << fDiff2->GetProb() << '\n';
    std::cout << "------------------------------------\n";

    // Drawing histograms on file2
    file2->cd();
    c->Write();
    cEnergy->Write();
    cImpulse->Write();
    c3DAngles->Write();
    cInvMass->Write();
    cResonance->Write();

    // Saving Canvas as .pdf, .C, .root and .jpeg
    c->Print("VariousDistributions.pdf");
    c->Print("VariousDistributions.C");
    c->Print("VariousDistributions.root");
    c->Print("VariousDistributions.jpeg");
    cEnergy->Print("EnergyDistribution.pdf");
    cEnergy->Print("EnergyDistribution.C");
    cEnergy->Print("EnergyDistribution.root");
    cEnergy->Print("EnergyDistribution.jpeg");
    cImpulse->Print("Impulse.pdf");
    cImpulse->Print("Impulse.C");
    cImpulse->Print("Impulse.root");
    cImpulse->Print("Impulse.jpeg");
    c3DAngles->Print("3DAngles.pdf");
    c3DAngles->Print("3DAngles.C");
    c3DAngles->Print("3DAngles.root");
    c3DAngles->Print("3DAngles.jpeg");
    cInvMass->Print("InvMass.pdf");
    cInvMass->Print("InvMass.C");
    cInvMass->Print("InvMass.root");
    cInvMass->Print("InvMass.jpeg");
    cResonance->Print("ResonanceSignal.pdf");
    cResonance->Print("ResonanceSignal.C");
    cResonance->Print("ResonanceSignal.root");
    cResonance->Print("ResonanceSignal.jpeg");

    file2->Close();
    file1->Close();
}
