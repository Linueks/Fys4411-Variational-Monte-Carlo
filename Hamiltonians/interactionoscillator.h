#pragma once
#include "hamiltonian.h"
//#include "harmonicoscillator.h"
#include <vector>

class InteractionOscillatorNum : public Hamiltonian{
public:
    InteractionOscillatorNum(System* system, double omega_ho, double omega_r);
    double computePotentialEnergy(std::vector<Particle*> particles);

protected:
    double m_omega_ho = 0;
    double m_omega_r = 0;
};

class InteractionOscillator : public InteractionOscillatorNum{
public:
    InteractionOscillator(System* system, double omega_ho, double omega_r);
    double computeKineticEnergy(std::vector<Particle*> particles);
};
