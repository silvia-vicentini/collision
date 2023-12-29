#include "resonancetype.hpp"

// constructor
ResonanceType::ResonanceType(const char *name, const double mass, const int charge, double const width) : ParticleType(name, mass, charge), fWidth(width)
{
}

double ResonanceType::GetWidth() const
{
    return fWidth;
}

void ResonanceType::Print() const
{
    ParticleType::Print();
    std::cout << "Particle width = " << fWidth << '\n';
}
