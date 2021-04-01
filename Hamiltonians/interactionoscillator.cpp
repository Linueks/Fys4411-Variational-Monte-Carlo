#include "interactionoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;


InteractionOscillatorNum::InteractionOscillatorNum(System* system, double omega_ho, double omega_r) :
    Hamiltonian(system){
        assert(omega_ho > 0);
        assert(omega_r > 0);
        m_omega_ho = omega_ho;
        m_omega_r = omega_r;
}

InteractionOscillator::InteractionOscillator(System* system, double omega_ho, double omega_r) :
    InteractionOscillatorNum(system, omega_ho, omega_r){
        assert(omega_ho > 0);
        assert(omega_r > 0);
        m_omega_ho = omega_ho;
        m_omega_r = omega_r;     
}

double InteractionOscillator::computeKineticEnergy(std::vector<Particle*> particles){
    WaveFunction* wavefunction = m_system->getWaveFunction();
    double kineticEnergy = 0;
    kineticEnergy = -0.5 * wavefunction->computeDoubleDerivative(particles);

    return kineticEnergy / wavefunction->evaluate(particles);
}

double InteractionOscillatorNum::computePotentialEnergy(std::vector<Particle*> particles){
    double potentialEnergy = 0;
    double internalEnergy = 0;
    double rxySquared = 0;
    double rzSquared = 0;
    double ijNorm = 0;
    const int numberOfParticles = m_system->getNumberOfParticles();
    const int numberOfDimensions = m_system->getNumberOfDimensions();
    double m_a = m_system->getCutoffRadius();

    for (int i=0; i < numberOfParticles; i++){
        std::vector<double> ithParticle = particles[i]->getPosition();
        for (int d=0; d < numberOfDimensions; d++){
            if (d==2){
                rzSquared += ithParticle[d]*ithParticle[d];
            }
            else{
                rxySquared += ithParticle[d]*ithParticle[d];
            }
            //cout << "Hamiltonian interaction, " <<"i:"<< i <<", d:"<< d << ", term:" << rxySquared << endl;
  
        }

        for (int j=i+1; j < numberOfParticles; j++){
            std::vector<double> jthParticle = particles[i]->getPosition();
            for (int d=0; d < numberOfDimensions; d++){
                ijNorm += (jthParticle[d] - ithParticle[d]) * (jthParticle[d] - ithParticle[d]);
            }
            ijNorm = sqrt(ijNorm);
            if (ijNorm < m_a){
                // some big number, technically is infinite
                potentialEnergy += 1e6;
            }
        }
    }
    return potentialEnergy += 0.5 * (m_omega_ho*m_omega_ho*rxySquared + m_omega_r*m_omega_r*rzSquared);
}