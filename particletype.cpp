#include "particletype.hpp"

// constructor
ParticleType::ParticleType() : fMass(0), fCharge(0) {}

ParticleType::ParticleType(const char *name,
                           const double mass,
                           const int charge) : fName(name), fMass(mass), fCharge(charge)
{
}
const char *ParticleType::GetName() const
{
    return fName;
}
double ParticleType::GetMass() const
{
    return fMass;
}
int ParticleType::GetCharge() const
{
    return fCharge;
}
double ParticleType::GetWidth() const { return 0; }

void ParticleType::Print() const
{
    std::cout << "-----------------------------------\n";
    std::cout << "Particle name = " << fName << '\n';
    std::cout << "Particle mass = " << fMass << '\n';
    if (fCharge > 0)
    {
        std::cout << "Particle charge = +" << fCharge << '\n';
    }
    else
    {
        std::cout << "Particle charge = " << fCharge << '\n';
    }
}
