#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "particle.hpp"

#include "doctest.h"

#include <vector>

// va bene la width usata?

TEST_CASE("Testing ParticleType and ResonanceType's methods")
{
    SUBCASE("ParticleType")
    {
        ParticleType particle1("e", 0.511, -1);
        CHECK(particle1.GetName() == "e");
        CHECK(particle1.GetMass() == doctest::Approx(0.511));
        CHECK(particle1.GetCharge() == -1);

        particle1.Print();
    }

    SUBCASE("ResonanceType")
    {
        ResonanceType particle2("K*", 497.648, 0, 0.6);
        CHECK(particle2.GetName() == "K*");
        CHECK(particle2.GetMass() == doctest::Approx(497.648));
        CHECK(particle2.GetCharge() == 0);
        CHECK(particle2.GetWidth() == 0.6);

        particle2.Print();
    }

    SUBCASE("Testing the Print method")
    {
        ParticleType particle3("p", 938.272, +1);
        ResonanceType particle4("K+", 493.677, +1, 0.9);
        std::vector<ParticleType *> test;
        test.push_back(&particle3);
        test.push_back(&particle4);
        for (ParticleType *particle : test)
        {
            particle->Print();
        }
    }

    SUBCASE("Testing the const methods for ParticleType")
    {
        const ParticleType particle5("π+", 139.6, +1);
        CHECK(particle5.GetName() == "π+");
        CHECK(particle5.GetMass() == doctest::Approx(139.6));
        CHECK(particle5.GetCharge() == 1);
        particle5.Print();
    }

    SUBCASE("Testing the const methods for ResonanceType")
    {
        const ResonanceType particle6("π-", 139.6, -1, 0.3);
        CHECK(particle6.GetName() == "π-");
        CHECK(particle6.GetMass() == doctest::Approx(139.6));
        CHECK(particle6.GetCharge() == -1);
        CHECK(particle6.GetWidth() == 0.3);
        particle6.Print();
    }
}

TEST_CASE("Testing Particle's methods")
{
    // m is expressed in GeV/c
    Particle::AddParticleType("e", 0.511, -1);
    Particle::AddParticleType("K*", 497.648, 0, 0.6);
    Particle::AddParticleType("p", 938.272, +1);
    Particle::AddParticleType("π-", 139.6, -1, 0.3);

    // p is expressed in GeV
    Particle particle1("e", 0.4599, 0.44457, 0.36792);
    Particle particle2("K*");
    Particle particle3("p", -938.272);
    const Particle particle4("π-", 110.284, -115.868, 114.472);
    Particle particle5("p-");
    Particle particle6("p", 798.272);

    SUBCASE("Testing FindParticle function")
    {
        CHECK(particle1.GetIndex() == 0);
        CHECK(particle2.GetIndex() == 1);
        CHECK(particle3.GetIndex() == 2);
        CHECK(particle4.GetIndex() == 3);
        CHECK(particle5.GetIndex() == -1);
        CHECK(particle6.GetIndex() == 2);
    }

    SUBCASE("Testing getter methods")
    {
        CHECK(Particle::GetSize() == 4);

        CHECK(particle1.GetPx() == doctest::Approx(0.4599));
        CHECK(particle1.GetPy() == doctest::Approx(0.44457));
        CHECK(particle1.GetPz() == doctest::Approx(0.36792));
        CHECK(particle2.GetPx() == 0);
        CHECK(particle2.GetPy() == 0);
        CHECK(particle2.GetPz() == 0);
        CHECK(particle3.GetPx() == doctest::Approx(-938.272));
        CHECK(particle3.GetPy() == 0);

        CHECK(particle1.GetMass() == 0.511);
        CHECK(particle1.GetCharge() == -1);
    }

    SUBCASE("Testing Energy and InvMass methods")
    {

        CHECK(particle1.Energy() == doctest::Approx(0.897572627));
        CHECK(particle2.Energy() == doctest::Approx(497.648));

        CHECK(particle1.InvMass(particle2) == doctest::Approx(498.5450265));
    }

    Particle::PrintfParticleType();

    particle1.PrintParticle();
    particle2.PrintParticle();
    particle3.PrintParticle();
    particle4.PrintParticle();
    particle5.PrintParticle();

    SUBCASE("Testing setter methods")
    {
        Particle particle7;
        particle7.SetIndex(3);
        CHECK(particle7.GetIndex() == 3);
        CHECK(particle7.GetPx() == 0);
        CHECK(particle7.GetPy() == 0);
        CHECK(particle7.GetPz() == 0);

        particle7.SetIndex("e");
        CHECK(particle7.GetIndex() == 0);

        particle7.SetP(323.47, 353.33, 467.80);
        CHECK(particle7.GetPx() == 323.47);
        CHECK(particle7.GetPy() == 353.33);
        CHECK(particle7.GetPz() == 467.80);
    }
}
