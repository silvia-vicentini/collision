#include "particle.hpp"

#include <iomanip>
#include <cmath>
#include <cassert>

Particle::Particle()
{
    fIndex = -1;
    fPx = 0;
    fPy = 0;
    fPz = 0;
}

Particle::Particle(const char *name, double Px, double Py, double Pz) : fName{name}, fPx{Px}, fPy{Py}, fPz{Pz}
{
    fIndex = FindParticle(name);
    if (fIndex == -1)
    {
        std::cerr << '\n'
                  << "Error : " << name << " is not one of the types considered.\n";
    }
}

std::vector<ParticleType *> Particle::fParticleType{};

int Particle::FindParticle(const char *name)
{
    for (int i = 0; i < Particle::GetSize(); ++i)
    {
        if (name == fParticleType[i]->GetName())
        {
            return i;
        }
    }
    return -1;
};

double Particle::GetPx() const { return fPx; };
double Particle::GetPy() const { return fPy; };
double Particle::GetPz() const { return fPz; };
int Particle::GetIndex() const { return fIndex; };
double Particle::GetMass() const { return fParticleType[fIndex]->GetMass(); };
int Particle::GetCharge() const { return fParticleType[fIndex]->GetCharge(); };
int Particle::GetSize() { return fParticleType.size(); }

double Particle::InvMass(Particle &p)
{
    double p_sum_norm2 = pow(p.GetPx() + GetPx(), 2) + pow(p.GetPy() + GetPy(), 2) + pow(p.GetPz() + GetPz(), 2);
    return sqrt(pow((Energy() + p.Energy()), 2) - p_sum_norm2);
}

double Particle::Energy() const
{
    double norm2 = fPx * fPx + fPy * fPy + fPz * fPz;
    return sqrt(GetMass() * GetMass() + norm2);
}

void Particle::AddParticleType(const char *name, const double mass, const int charge, double width)
{
    if (GetSize() < fMaxNumParticleType)
    {
        if (FindParticle(name) != -1)
        {
            std::cerr << '\n'
                      << name << " is alseady included.\n";
        }
        else
        {
            if (width == 0)
            {
                ParticleType *particle = new ParticleType(name, mass, charge);
                fParticleType.push_back(particle);
            }
            else
            {
                ResonanceType *particle = new ResonanceType(name, mass, charge, width);
                fParticleType.push_back(particle);
            }
        }
    }
    else
    {
        std::cerr << "Can't include more than " << fMaxNumParticleType << " types.\n";
    }
}

void Particle::SetIndex(const int index)
{
    if (index >= 0 && index <= GetSize())
    {
        fIndex = index;
        fName = fParticleType[fIndex]->GetName();
    }
    else
    {
        std::cerr << "Particle can't be one of the possible types.\n";
    }
};

void Particle::SetIndex(const char *name)
{
    SetIndex(FindParticle(name));
};

void Particle::PrintfParticleType()
{
    std::cout << "Particle Types\n";
    std::cout << std::setw(12) << std::left << "Name" << std::setw(12)
              << std::left << "Mass" << std::setw(12) << std::left
              << "Charge" << std::setw(6) << std::left << "Width" << '\n';
    for (ParticleType *particle : fParticleType)
    {
        std::cout << std::setw(12) << std::left << particle->GetName() << std::setw(12)
                  << std::left << particle->GetMass() << std::setw(12) << std::left
                  << particle->GetCharge() << std::setw(6) << std::left << particle->GetWidth() << '\n';
    }
};

void Particle::PrintParticle() const
{
    if (fIndex != -1)
    {
        std::cout << "-----------------------------------\n";
        std::cout << "Particle " << fIndex << " = " << fName << "\n";
        std::cout << "Px = " << GetPx() << "\n";
        std::cout << "Py = " << GetPy() << "\n";
        std::cout << "Pz = " << GetPz() << "\n";
    }
    else
    {
        std::cout << "-----------------------------------\n";
        std::cout << "Can't return particle index because it's not one of the types considered.\n";
        std::cout << "Particle name = " << fName << "\n";
        std::cout << "Px = " << GetPx() << "\n";
        std::cout << "Py = " << GetPy() << "\n";
        std::cout << "Pz = " << GetPz() << "\n";
    }
}

void Particle::SetP(double Px, double Py, double Pz)
{
    fPx = Px;
    fPy = Py;
    fPz = Pz;
}

void Particle::Boost(double bx, double by, double bz)
{
    double energy = Energy();

    // Boost this Lorentz vector
    double b2 = bx * bx + by * by + bz * bz;
    double gamma = 1.0 / sqrt(1.0 - b2);
    double bp = bx * fPx + by * fPy + bz * fPz;
    double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

    fPx += gamma2 * bp * bx + gamma * bx * energy;
    fPy += gamma2 * bp * by + gamma * by * energy;
    fPz += gamma2 * bp * bz + gamma * bz * energy;
}

int Particle::Decay2Body(Particle &dau1, Particle &dau2) const
{
    if (GetMass() == 0.0)
    {
        printf("Decayment cannot be preformed if mass is zero\n");
        return 1;
    }

    double massMot = GetMass();
    double massDau1 = dau1.GetMass();
    double massDau2 = dau2.GetMass();

    if (fIndex > -1)
    { // add width effect

        // gaussian random numbers

        float x1, x2, w, y1;

        double invnum = 1. / RAND_MAX;
        do
        {
            x1 = 2.0 * rand() * invnum - 1.0;
            x2 = 2.0 * rand() * invnum - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);

        w = sqrt((-2.0 * log(w)) / w);
        y1 = x1 * w;

        massMot += fParticleType[fIndex]->GetWidth() * y1;
    }

    if (massMot < massDau1 + massDau2)
    {
        printf("Decayment cannot be preformed because mass is too low in this channel\n");
        return 2;
    }

    double pout = sqrt((massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) * (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) / massMot * 0.5;

    double norm = 2 * M_PI / RAND_MAX;

    double phi = rand() * norm;
    double theta = rand() * norm * 0.5 - M_PI / 2.;
    dau1.SetP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi), pout * cos(theta));
    dau2.SetP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi), -pout * cos(theta));

    double energy = sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

    double bx = fPx / energy;
    double by = fPy / energy;
    double bz = fPz / energy;

    dau1.Boost(bx, by, bz);
    dau2.Boost(bx, by, bz);

    return 0;
}
