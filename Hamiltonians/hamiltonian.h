#pragma once
#include <vector>
#include "../WaveFunctions/wavefunction.h"
#include "../particle.h"



class Hamiltonian {
public:
    Hamiltonian(class System* system);
    double computeLocalEnergy(std::vector<class Particle*> particles);
    virtual double computeKineticEnergy(std::vector<Particle *> particles);
    virtual double computePotentialEnergy(std::vector<class Particle*> particles) = 0;

protected:
    class System* m_system = nullptr;

private:
    double derivativeStepSize = 0;
    double derivativeSquaredInverse = 0;
    double kineticEnergy = 0;
};

