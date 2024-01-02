#include "resonancetype.hpp"

ResonanceType::ResonanceType(const char *name, const double mass, const int charge, const double width) : ParticleType(name, mass, charge), fWidth(width)
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
