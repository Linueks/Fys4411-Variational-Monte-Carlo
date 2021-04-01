#pragma once
#include "wavefunction.h"

class Interacting : public WaveFunction {
public:
    Interacting(System* system, double alpha, double beta, double a);
    double evaluate(std::vector<Particle*> particles);
    double alphaDerivative(std::vector<Particle*> particles);
    double computeDoubleDerivative(std::vector<Particle*> particles);
    double uPrimeR(std::vector<double> ithParticle, std::vector<double> jthParticle);
    double uDoublePrimeR(std::vector<double> ithParticle, std::vector<double> kthParticle);
};
