#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillatorNum : public Hamiltonian{
public:
    HarmonicOscillatorNum(System* system, double omega_ho, double omega_r);
    double computePotentialEnergy(std::vector<Particle*> particles);

protected:
    double m_omega_ho = 0;
    double m_omega_r = 0;
};

class HarmonicOscillator : public HarmonicOscillatorNum{
public:
    HarmonicOscillator(System* system, double omega_ho, double omega_r);
    double computeKineticEnergy(std::vector<Particle*> particles);
};
