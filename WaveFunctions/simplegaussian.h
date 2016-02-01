#pragma once
#include "wavefunction.h"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha);
    double evaluate(class Particle* particle);
    double computeKineticEnergy(class Particle* particles);
};
