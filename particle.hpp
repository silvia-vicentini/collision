#ifndef PARTICLE_HPP
#define PARTICLE_HPP
#include "resonancetype.hpp"
#include <vector>

class Particle
{
private:
    // Container of the differet types of particles generated
    static std::vector<ParticleType *> fParticleType;

    static const int fMaxNumParticleType = 10;
    int fIndex;
    const char *fName;
    double fPx;
    double fPy;
    double fPz;

    // Funtion to define the corresponding type of the particle above those incuded in fParticleType
    static int FindParticle(const char *);

    // Function used in Decay2Body
    void Boost(double, double, double);

public:
    // Constructor
    Particle(const char *, double Px = 0, double Py = 0, double Pz = 0);

    // Default constructor
    Particle();

    // Getter methods
    double GetPx() const;
    double GetPy() const;
    double GetPz() const;
    int GetIndex() const;
    double GetMass() const;
    int GetCharge() const;
    static int GetSize();

    // Function to calculate the energy of a particle
    double Energy() const;

    // Function to calculate the invariant mass of a particle fron the products of decayment
    double InvMass(Particle &) const;

    // Function to add a type of particle to fParticleType
    static void AddParticleType(const char *, const double, const int, double width = 0);

    // Setter methods
    void SetIndex(const int);
    void SetIndex(const char *);
    void SetP(double, double, double);

    // Printer methods
    static void PrintfParticleType();
    void PrintParticle() const;

    // Function to set the products of decayment
    int Decay2Body(Particle &, Particle &) const;
};

#endif
