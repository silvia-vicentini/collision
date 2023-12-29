#ifndef RESONANCETYPE_HPP
#define RESONANCETYPE_HPP
#include "particletype.hpp"

class ResonanceType : public ParticleType
{
private:
    const double fWidth;

public:
    // costruttore
    ResonanceType(const char *, const double, const int, const double);
    double GetWidth() const;
    void Print() const;
};

#endif