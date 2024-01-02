#ifndef PARTICLETYPE_HPP
#define PARTICLETYPE_HPP
#include <iostream>
#include <string>

class ParticleType
{
private:
    const char *fName;
    const double fMass;
    const int fCharge;

public:
    // Constructors
    ParticleType(const char *, const double, const int);
    ParticleType();

    // Getter methods
    const char *GetName() const;
    double GetMass() const;
    int GetCharge() const;
    virtual double GetWidth() const;

    // Printer method
    virtual void Print() const;
};
#endif
