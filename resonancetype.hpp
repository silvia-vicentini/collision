#ifndef RESONANCETYPE_HPP
#define RESONANCETYPE_HPP
#include "particletype.hpp"

class ResonanceType : public ParticleType
{
private:
    const double fWidth;

public:
    // Constructor
    ResonanceType(const char *, const double, const int, const double);

    // Getter method
    double GetWidth() const override;

    // Printer method
    void Print() const override;
};

#endif
