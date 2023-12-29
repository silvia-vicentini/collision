#ifndef PARTICLE_HPP
#define PARTICLE_HPP
#include "resonancetype.hpp"
#include <vector>

// scarica root su ubuntu!

class Particle
{
private:
    static std::vector<ParticleType *> fParticleType;
    static const int fMaxNumParticleType = 10;
    int fIndex;
    const char *fName;
    double fPx;
    double fPy;
    double fPz;
    static int FindParticle(const char *);
    void Boost(double, double, double);

public:
    // costruttore
    Particle(const char *, double Px = 0, double Py = 0, double Pz = 0);

    // Default constructor
    Particle();

    double GetPx() const;
    double GetPy() const;
    double GetPz() const;
    int GetIndex() const;
    double GetMass() const;
    int GetCharge() const;
    static int GetSize();

    double Energy() const;
    double InvMass(Particle &);

    static void AddParticleType(const char *, const double, const int, double width = 0);

    void SetIndex(const int);
    void SetIndex(const char *);
    void SetP(double, double, double);

    static void PrintfParticleType();
    void PrintParticle() const;

    int Decay2Body(Particle &, Particle &) const;
};

#endif