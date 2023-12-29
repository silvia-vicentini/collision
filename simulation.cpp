#include "particle.hpp"

#include "TFile.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TStyle.h"

/*void setStyle()
{
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(57);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(1111);
}*/

void simulation()
{
    Particle::AddParticleType("pi+", 0.13957, +1);
    Particle::AddParticleType("pi-", 0.13957, -1);
    Particle::AddParticleType("K+", 0.49367, +1);
    Particle::AddParticleType("K-", 0.49367, -1);
    Particle::AddParticleType("p+", 0.93827, +1);
    Particle::AddParticleType("p-", 0.93827, -1);
    Particle::AddParticleType("K*", 0.89166, 0, 0.050);

    gRandom->GetSeed();

    TFile *file = new TFile("simulation.root", "RECREATE");

    TH1F *h0 = new TH1F("h0", "Paricle types", 7, 0, 7);
    TH1F *h1 = new TH1F("h1", "Impulse modulus", 500, 0, 10);
    TH1F *h2 = new TH1F("h2", "Trasvers impulse modulus", 500, 0, 10);
    TH1F *h3 = new TH1F("h3", "Particel energy", 500, 0, 10);
    TH1F *h4 = new TH1F("h4", "Azimuthal angle", 100, 0, 2 * TMath::Pi());
    TH1F *h5 = new TH1F("h5", "Polar angle", 100, 0, TMath::Pi());
    TH3F *h12 = new TH3F("h12", "3D azimuthal and polar angles", 100, -1, 1, 100, -1, 1, 100, -1, 1);
    TH1F *h6 = new TH1F("h6", "Invariant mass between discordant charges", 200, 0, 2);
    TH1F *h7 = new TH1F("h7", "Invariant mass between concordant charges", 200, 0, 2);
    TH1F *h8 = new TH1F("h8", "Invariant mass between #pi+/K- and #pi-/K+", 200, 0, 2);
    TH1F *h9 = new TH1F("h9", "Invariant mass between #pi+/K+ and #pi-/K-", 200, 0, 2);
    TH1F *h10 = new TH1F("h10", "Invariant mass between all particels", 200, 0, 2);
    TH1F *h11 = new TH1F("h11", "Invariant mass between particels from decayment", 200, 0, 2);

    h6->Sumw2();
    h7->Sumw2();
    h8->Sumw2();
    h9->Sumw2();
    h10->Sumw2();
    h11->Sumw2();

    int nGen = 1E5;
    for (int j = 0; j < nGen; ++j)
    {
        int nParticels = 100;
        std::vector<Particle> EventParticle;
        double phi, theta, p;
        double x = 0;
        for (int i = 0; i < nParticels; ++i)
        {
            Particle particle;
            x = gRandom->Rndm();
            if (x < 0.4)
            {
                particle.SetIndex(0);
            }
            else if (x < 0.8)
            {
                particle.SetIndex(1);
            }
            else if (x < 0.85)
            {
                particle.SetIndex(2);
            }
            else if (x < 0.9)
            {
                particle.SetIndex(3);
            }
            else if (x < 0.945)
            {
                particle.SetIndex(4);
            }
            else if (x < 0.99)
            {
                particle.SetIndex(5);
            }
            else
            {
                particle.SetIndex(6);
            }
            EventParticle.push_back(particle);

            h0->Fill(EventParticle[i].GetIndex());
            phi = gRandom->Rndm() * TMath::Pi();
            theta = gRandom->Rndm() * 2 * TMath::Pi();
            p = gRandom->Exp(1);
            EventParticle[i].SetP(p * sin(theta) * cos(phi), p * sin(theta) * sin(phi), p * cos(theta));
            h12->Fill(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
            h1->Fill(p);
            h2->Fill(p * sin(theta));
            h3->Fill(EventParticle[i].Energy());
            h4->Fill(theta);
            h5->Fill(phi);

            if (EventParticle[i].GetIndex() == 6)
            {
                double y = 0;
                y = gRandom->Rndm();
                if (y < 0.5)
                {
                    Particle dau1("pi-");
                    Particle dau2("K+");
                    EventParticle[i].Decay2Body(dau1, dau2);
                    // EventParticle.pop_back();
                    EventParticle.push_back(dau1);
                    EventParticle.push_back(dau2);
                    h11->Fill(dau1.InvMass(dau2));
                }
                else
                {
                    Particle dau1("pi+");
                    Particle dau2("K-");
                    EventParticle[i].Decay2Body(dau1, dau2);
                    // EventParticle.pop_back();
                    EventParticle.push_back(dau1);
                    EventParticle.push_back(dau2);
                    h11->Fill(dau1.InvMass(dau2));
                }
            }
        }

        int size = EventParticle.size();
        /* if (size > 100)
         {
             std::cout << size << '\n';
         }*/

        for (int a = 0; a < size; ++a)
        {
            for (int b = a + 1; b < size; ++b)
            {
                h10->Fill(EventParticle[a].InvMass(EventParticle[b]));
                if (EventParticle[a].GetCharge() * EventParticle[b].GetCharge() < 0)
                {
                    h6->Fill(EventParticle[a].InvMass(EventParticle[b]));
                }
                if (EventParticle[a].GetCharge() * EventParticle[b].GetCharge() > 0)
                {
                    h7->Fill(EventParticle[a].InvMass(EventParticle[b]));
                }
                if ((EventParticle[a].GetIndex() == 0 && EventParticle[b].GetIndex() == 3) || (EventParticle[a].GetIndex() == 3 && EventParticle[b].GetIndex() == 0) || (EventParticle[a].GetIndex() == 1 && EventParticle[b].GetIndex() == 2) || (EventParticle[a].GetIndex() == 2 && EventParticle[b].GetIndex() == 1))
                {
                    h8->Fill(EventParticle[a].InvMass(EventParticle[b]));
                }
                if ((EventParticle[a].GetIndex() == 0 && EventParticle[b].GetIndex() == 2) || (EventParticle[a].GetIndex() == 2 && EventParticle[b].GetIndex() == 0) || (EventParticle[a].GetIndex() == 1 && EventParticle[b].GetIndex() == 3) || (EventParticle[a].GetIndex() == 3 && EventParticle[b].GetIndex() == 1))
                {
                    h9->Fill(EventParticle[a].InvMass(EventParticle[b]));
                }
            }
        }
    }

    file->Write();
    file->Close();
}

// Add main in order to compile from SHELL
int main()
{
    setStyle();

    simulation();

    return EXIT_SUCCESS;
}