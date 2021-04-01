#include "hamiltonian.h"
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include "../particle.h"
#include <iostream>

using namespace std;



Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeKineticEnergy(std::vector<Particle *> particles) {
    //std::cout << "hei" << std::endl;
    double derivativeStepSize = m_system->getDerivativeStepSize();
    double derivativeSquaredInverse = 1 / (derivativeStepSize*derivativeStepSize);
    const int numberOfParticles = m_system->getNumberOfParticles();
    const int numberOfDimensions = m_system->getNumberOfDimensions();

    WaveFunction* waveFunction = m_system->getWaveFunction();

    double rSquared = 0;
    double kineticEnergy = 0;

    for (int i=0; i<numberOfParticles; i++){
        std::vector<double> particlePosition = particles[i]->getPosition();
        for (int j=0; j<numberOfDimensions; j++){
                rSquared += particlePosition[j] * particlePosition[j];
                double waveFunctionOld = waveFunction->evaluate(particles);
    
                particles[i]->adjustPosition(-derivativeStepSize, j);
                double waveFunctionMinus = waveFunction->evaluate(particles);
    
                particles[i]->adjustPosition(2*derivativeStepSize, j);
                double waveFunctionPlus = waveFunction->evaluate(particles);
                particles[i]->adjustPosition(-derivativeStepSize, j);

                kineticEnergy -= (waveFunctionMinus+waveFunctionPlus - 2*waveFunctionOld);
        }
    }
    return 0.5 * derivativeSquaredInverse * kineticEnergy / waveFunction->evaluate(particles);
}

double Hamiltonian::computeLocalEnergy(std::vector<Particle*> particles) {
    double kineticEnergy = computeKineticEnergy(particles);
    double potentialEnergy = computePotentialEnergy(particles);
    //cout <<"kinetic: "<< kineticEnergy<<", potential: "<<potentialEnergy<<endl;
    return kineticEnergy + potentialEnergy;
    //return computePotentialEnergy(particles);
}

