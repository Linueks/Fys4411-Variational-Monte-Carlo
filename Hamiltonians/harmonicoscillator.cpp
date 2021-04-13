#pragma once
#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillatorNum::HarmonicOscillatorNum(System* system, double omega_ho, double omega_r) :
    Hamiltonian(system){
        assert(omega_ho > 0);
        assert(omega_r > 0);
        m_omega_ho = omega_ho;
        m_omega_r = omega_r;
}

HarmonicOscillator::HarmonicOscillator(System* system, double omega_ho, double omega_r) :
    HarmonicOscillatorNum(system, omega_ho, omega_r){
        assert(omega_ho > 0);
        assert(omega_r > 0);
        m_omega_ho = omega_ho;
        m_omega_r = omega_r;
}

double HarmonicOscillator::computeKineticEnergy(std::vector<Particle*> particles){
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */
    WaveFunction* waveFunction = m_system->getWaveFunction();

    double kineticEnergy = 0;

    kineticEnergy = -0.5 * waveFunction->computeDoubleDerivative(particles);

    return kineticEnergy;
}

double HarmonicOscillatorNum::computePotentialEnergy(std::vector<Particle*> particles){
    double potentialEnergy = 0;
    const int numberOfParticles = m_system->getNumberOfParticles();
    const int numberOfDimensions = m_system->getNumberOfDimensions();
    bool interactionOrNot = m_system->getInteractionOrNot();


    // last change here, moved the true code below from interactionsoscillator.cpp to try to fix it'
    // ended up with two seemingly identical codes below that produce different values for 1 particle 1 dimension. 
    // they should be equal unless im completely wrong. 
/*     if (interactionOrNot){
        double internalEnergy = 0;
        double rxySquared = 0;
        double rzSquared = 0;
        double ijNorm = 0;
        double m_a = m_system->getCutoffRadius();

        cout << "got in here" << endl;

        for (int i=0; i < numberOfParticles; i++){
            std::vector<double> ithParticle = particles[i]->getPosition();
            for (int d=0; d < numberOfDimensions; d++){
                if (d==2){
                    rzSquared += ithParticle[d]*ithParticle[d];
                }
                else{
                    rxySquared += ithParticle[d]*ithParticle[d];
                }
                cout << "Hamiltonian interaction, " <<"i:"<< i <<", d:"<< d << ", term:" << rxySquared << endl;
    
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
    } */

    //else{
        double rSquared = 0;
        for (int i=0; i < numberOfParticles; i++){
            std::vector<double> particlePosition = particles[i]->getPosition();
            for (int d=0; d < numberOfDimensions; d++){
                    rSquared += particlePosition[d] * particlePosition[d];
                //cout << "Hamiltonian no interaction, " <<"i:"<< i <<", d:"<< d << ", term:" << rSquared << endl;
                }       
            }
        return potentialEnergy = 0.5 * m_omega_ho*m_omega_ho * rSquared;
    //}
}