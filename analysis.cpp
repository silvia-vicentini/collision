#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include <iomanip>
#include <vector>
#include <iostream>

void setStyle()
{
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(57);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(1);
}

void histo_cosmetic(TH1F *h)
{
    int fillColor = 426;
    h->SetFillColor(fillColor);
    h->SetLineWidth(1);
    h->SetLineColor(kBlack);
    h->GetYaxis()->SetTitle("Entries");
    h->DrawCopy("H");
    h->DrawCopy("E,SAME");
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.5);
}

void function_cosmetic(TF1 *f)
{
    f->SetLineStyle(2);
    f->SetLineWidth(3);
    f->SetLineColor(kRed);
}

void analysis()
{
    using namespace std;

    TFile *file1 = new TFile("simulation.root", "READ");
    TFile *file2 = new TFile("analysis.root", "RECREATE");
    file1->ls();

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

    TCanvas *c = new TCanvas("c", "c", 200, 10, 600, 400);
    c->Divide(2, 2);
    TCanvas *cAngles = new TCanvas("cAngles", "cAngles", 200, 10, 600, 400);
    cAngles->Divide(2, 2);
    TCanvas *cInvMass = new TCanvas("cInvMass", "cInvMass", 200, 10, 600, 400);
    cInvMass->Divide(3, 2);
    TCanvas *cDiff = new TCanvas("cDiff", "cDiff", 200, 10, 600, 400);
    cDiff->Divide(2, 1);
    TCanvas *cvuota = new TCanvas("cvuota", "cvuota", 200, 10, 600, 400);

    TF1 *f1 = new TF1("f1", "expo", 0, 10);
    TF1 *f4 = new TF1("f4", "pol0", 0, 2 * TMath::Pi());
    TF1 *f5 = new TF1("f5", "pol0", 0, TMath::Pi());

    TH1F *hDiff1 = new TH1F("hDiff1", "Diff1", 200, 0, 2);
    TH1F *hDiff2 = new TH1F("hDiff2", "Diff2", 200, 0, 2);

    hDiff1->Add(h6, h7, 1, -1);
    hDiff2->Add(h8, h9, 1, -1);

    TF1 *fDiff1 = new TF1("fDiff1", "gaus", 0, 2);
    TF1 *fDiff2 = new TF1("fDiff2", "gaus", 0, 2);

    double k_mass = 0.89166;
    double k_width = 0.050;
    f1->SetParameter(0, -1);
    f4->SetParameter(0, 1E5);
    f5->SetParameter(0, 1E5);
    fDiff1->SetParameters(k_mass, k_width);
    fDiff2->SetParameters(k_mass, k_width);

    function_cosmetic(f1);
    function_cosmetic(f4);
    function_cosmetic(f5);
    function_cosmetic(fDiff1);
    function_cosmetic(fDiff2);

    h1->Fit("f1");
    h4->Fit("f4");
    h5->Fit("f5");
    hDiff1->Fit("fDiff1");
    hDiff2->Fit("fDiff2");

    vector<TH1F *> histo{h0, h1, h3, h2, h4, h5, h6, h7, h8, h9, h10, h11, hDiff1, hDiff2};
    vector<const char *> bin_names{"#pi+", "#pi-", "K+", "K-", "p+", "p-", "K*"};
    // bisogna cambiare la legenda
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
        }
        if (i == 1)
        {
            histo[i]->GetXaxis()->SetTitle("Impulse modulus (GeV/c)");
        }
        if (i == 2)
        {
            histo[i]->GetXaxis()->SetTitle("Energy (GeV)");
        }
        if (i == 3)
        {
            histo[i]->GetXaxis()->SetTitle("Trasvers impulse modulus (GeV/c)");
        }
        histo_cosmetic(histo[i]);
    }

    // Drawing on cAngles
    for (int i = 4; i < 6; ++i)
    {
        cAngles->cd(i - 3);
        if (i == 4)
        {
            histo[i]->GetXaxis()->SetTitle("#theta (rad)");
        }
        if (i == 5)
        {
            histo[i]->GetXaxis()->SetTitle("#phi (rad)");
        }
        histo[i]->SetMinimum(85 * 1E3);
        histo[i]->SetMaximum(115 * 1E3);
        histo_cosmetic(histo[i]);
    }

    cAngles->cd(3);
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
    // non disegna l'h11 perché?
    for (int i = 6; i < 12; ++i)
    {
        cInvMass->cd(i - 5);
        histo[i]->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
        histo_cosmetic(histo[i]);
    }

    // non li crea o non li disegna perché?
    for (int i = 12; i < 14; ++i)
    {
        cDiff->cd(i - 11);
        histo[i]->SetMinimum(-4000);
        histo[i]->SetMaximum(5000);
        if (i == 12)
        {
            histo[i]->SetTitle("Difference of invariant mass between discordant/concordant charges");
        }
        if (i == 13)
        {
            histo[i]->SetTitle("Difference of invariant mass between discordant/concordant types");
        }
        histo_cosmetic(histo[i]);
    }

    /*
        cDiff->cd(1);
        hDiff1->SetTitle("Difference of invariant mass between discordant/concordant charges");
        hDiff1->SetMinimum(-4000);
        hDiff1->SetMaximum(5000);
        hDiff1->SetAxisRange(0., 2.5, "X");
        histo_cosmetic(hDiff1);*/
    /*
        TLegend *leg1 = new TLegend(0.2, 0.2, .8, .8);
        TLegend *leg2 = new TLegend(0.2, 0.2, .8, .8);
        TLegend *leg3 = new TLegend(0.2, 0.2, .8, .8);
        TLegend *leg4 = new TLegend(0.2, 0.2, .8, .8);
        TLegend *leg5 = new TLegend(0.2, 0.2, .8, .8);
        TLegend *leg6 = new TLegend(0.2, 0.2, .8, .8);*/

    c->Print("VariousDistributions.pdf");
    cAngles->Print("Angles.pdf");
    cInvMass->Print("InvMass.pdf");
    cDiff->Print("InvMassDiff.pdf");
    c->Print("VariousDistributions.C");
    cAngles->Print("Angles.C");
    cInvMass->Print("InvMass.C");
    cDiff->Print("InvMassDiff.C");
    c->Print("VariousDistributions.root");
    cAngles->Print("Angles.root");
    cInvMass->Print("InvMass.root");
    cDiff->Print("InvMassDiff.root");

    file2->Close();
    file1->Close();

    /*


        cDiff->cd(2);
        hDiff2->SetTitle("Difference of invariant mass between discordant/concordant types");
        histo_cosmetic(hDiff2);*/

    // scrivi meglio

    cout << "------------------------------------------------------------------------\n";
    cout << "Number of \u03C0+ generated: " << left << setw(10) << h0->GetBinContent(1) << " ± " << h0->GetBinError(1)
         << " percentage: " << (h0->GetBinContent(1) / h0->GetBinError(1)) * 100 << "%\n";
    cout << "Number of \u03C0- generated: " << left << setw(10) << h0->GetBinContent(2) << " ± " << h0->GetBinError(2)
         << " percentage: " << (h0->GetBinContent(2) / h0->GetBinError(2)) * 100 << "%\n";
    for (int i{3}; i < 8; ++i)
    {
        cout << "Number of " << bin_names[i - 1] << left << setw(10)
             << " generated: " << h0->GetBinContent(i) << " ± " << h0->GetBinError(i)
             << " percentage: " << (h0->GetBinContent(i) / h0->GetBinError(i)) * 100 << "%\n";
    }
    cout << "------------------------------------------------------------------------\n";
    cout << "Azimuthal angle distribution parameter = " << f4->GetParameter(0) << " +/- " << f4->GetParError(0) << '\n';
    cout << "Azimutal angle distribution chisquare = " << f4->GetChisquare() / f4->GetNDF() << '\n';
    cout << "Azimutal angle distribution probability = " << f4->GetProb() << '\n';
    cout << "------------------------------------------------------------------------\n";
    cout
        << "Polar angle distribution parameter = " << f5->GetParameter(0) << " +/- " << f5->GetParError(0) << '\n';
    cout << "Polar angle distribution chisquare = " << f5->GetChisquare() / f5->GetNDF() << '\n';
    cout << "Polar angle distribution probability = " << f5->GetProb() << '\n';
    cout << "------------------------------------------------------------------------\n";
    cout << "Impulse modulus distribution parameter = " << f1->GetParameter(0) << " +/- " << f1->GetParError(0) << '\n';
    cout << "Impulse modulus distribution mean = " << h1->GetMean() << " +/- " << h1->GetMeanError() << '\n';
    cout << "Impulse modulus distribution mean set in phase of generation = 1\n";
    cout << "Impulse modulus distribution chisquare = " << f1->GetChisquare() / f1->GetNDF() << '\n';
    cout << "Impulse modulus distribution probability = " << f1->GetProb() << '\n';
    cout << "------------------------------------------------------------------------\n";
    /* cout << "Difference of invariant mass between discordant/concordant charges chisquare = " << fDiff1->GetChisquare() / fDiff1->GetNDF() << '\n';
     cout << "Difference of invariant mass between discordant/concordant charges probability = " << fDiff1->GetProb() << '\n';
     cout << "Difference of invariant mass between discordant/concordant types chisquare = " << fDiff2->GetChisquare() / fDiff2->GetNDF() << '\n';
     cout << "Difference of invariant mass between discordant/concordant types probability = " << fDiff2->GetProb() << '\n';
     cout << "K* mass from fDiff1 = " << fDiff1->GetParameter(0) << '\n';
     cout << "K* mass from fDiff2 = " << fDiff2->GetParameter(0) << '\n';
     cout << "K* width from fDiff1 = " << fDiff1->GetParameter(1) << '\n';
     cout << "K* width from fDiff2 = " << fDiff2->GetParameter(1) << '\n';*/
    // non mi crea il file 2

    // scrivi alla fine void main per poter compilare da terminale
}
// Add main in order to compile from SHELL
int main()
{
    analysis();

    return EXIT_SUCCESS;
}