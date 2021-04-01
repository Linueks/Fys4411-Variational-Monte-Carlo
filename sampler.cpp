#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;
using std::fstream;
using std::ios;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
    m_acceptedSteps = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep) {
    double localEnergy = 0;
    double alphaDerivative = 0;
    bool writeToFile = m_system->getWriteToFile();

    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
        m_deltaPsi = 0;
        m_psiEnergyDerivative = 0;
        m_derivative = 0;
    }

    if (writeToFile && m_stepNumber == 0) {
        m_energies = std::vector<double>();

    }

    m_acceptedSteps += acceptedStep;
    

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */

    localEnergy = m_system->getHamiltonian()->computeLocalEnergy(m_system->getParticles());
    alphaDerivative = m_system->getWaveFunction()->alphaDerivative(m_system->getParticles());

    m_deltaPsi += alphaDerivative;
    m_cumulativeEnergy += localEnergy;
    m_psiEnergyDerivative += alphaDerivative*localEnergy;
    m_stepNumber++;

    if (m_stepNumber % 1000 == 0) {
        cout << "i: " << m_stepNumber << endl;
    }


    if (writeToFile) {
        double energy = m_cumulativeEnergy / (m_stepNumber+1);
        m_energies.push_back(energy);
    }
}

void Sampler::printEnergiesToFile() {
        fstream outfile;
        //cout << m_system->getOutfileName() << endl;
        outfile.open(m_system->getOutfileName(), ios::app);
        if (!outfile) {
            cout << "Error in opening file..";
        }
        for (auto i : m_energies) {
            outfile << i << endl;
        }

        outfile.close();
}


void Sampler::printOutputToTerminal() {
    int np = m_system->getNumberOfParticles();
    int nd = m_system->getNumberOfDimensions();
    int ms = m_system->getNumberOfMetropolisSteps();
    int p = m_system->getWaveFunction()->getNumberOfParameters();
    double ef = m_system->getEquilibrationFraction();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();

    cout << endl;
    cout << "N particles: " << np;
    cout << ", N dimensons: " << nd;
    cout << ", Metrosteps: " << ms;
    cout << ", Equilsteps: " << ms*ef;
    cout << ", Acceptratio: " << m_acceptedRatio;
    cout << ", N parameters: " << p;
    for (int i=0; i < p; i++) {
        cout << ", Parameter " << i+1 << ": " << pa.at(i);
    }
    cout << ", Energy : " << m_energy;
}

void Sampler::computeAverages() {
    // Compute the averages of the sampled quantities.

    m_energy = m_cumulativeEnergy / m_stepNumber;
    m_psiAlphaDerivative = m_deltaPsi / m_stepNumber;
    m_derivative = m_psiEnergyDerivative / m_stepNumber;
    m_acceptedRatio = m_acceptedSteps / (double) m_numberOfMetropolisSteps;
}

double Sampler::getEnergyGradient() {
    return 2 * (m_derivative - m_psiAlphaDerivative*m_energy);
}

void Sampler::resetValues()
{
    m_stepNumber = 0;
    m_acceptedRatio = 0;
    m_cumulativeEnergy = 0;
    m_acceptedRatio = 0;
    m_acceptedSteps = 0;
    m_energy = 0;
    m_psiAlphaDerivative = 0;
    m_psiEnergyDerivative = 0;
}