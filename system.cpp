#include "system.h"
#include <cassert>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include "iostream"
#include "fstream"
#include "chrono"

using namespace std;


System::System() {
    m_random = new Random();
}

System::System(int seed) {
    m_random = new Random(seed);
}

void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
    //cout << "hamiltonian sat: " << m_hamiltonian << endl;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
    //cout << "wavefunction sat: " << m_waveFunction <<endl;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
    //cout << "initial state sat: " << m_initialState << endl;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
    //cout << "Equilibration Fraction sat: " << m_equilibrationFraction << endl;
}

void System::setInteraction(bool interactionOrNot) {
    m_interactionOrNot= interactionOrNot;
}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
    //cout << "metropolis step length sat: " << m_stepLength << endl;
}

void System::setCutoffRadius(double a){
    assert(a > 0);
    double m_a = a;
}

/*
void System::setVariationalStepLength(double variationStep) {
    assert(variationStep >= 0);
    m_variationStep = variationStep;
    cout << "variation step length sat: " << m_variationStep << endl;
}

void System::setAlphaStart(double alphaStart) {
    assert(alphaStart >= 0);
    m_alphaStart = alphaStart;
    cout << "alphaStart sat: " << m_alphaStart << endl;
}
*/

void System::setNumberOfMetropolisSteps(int numberOfSteps) {
    assert(numberOfSteps >= 0);
    m_numberOfMetropolisSteps = numberOfSteps;
    //cout << "Number of metropolis steps sat: " << m_numberOfMetropolisSteps << endl;
}

void System::setDerivativeStepSize(double derivativeStep) {
    assert(derivativeStep >= 0);
    m_derivativeStepSize = derivativeStep;
    //cout << "Numerical derivative step size sat: " << m_derivativeStepSize << endl;
}

void System::activateImportanceSampling(bool activate, double dt) {
    m_dt = dt;
    m_importanceSampling = activate;
    //cout << "Importance Sampling flag: " << m_importanceSampling << endl;
}

void System::writeToFile(bool writeToFile, string filename) {
    m_writeToFile = writeToFile;
    m_filename = filename;

    if (m_writeToFile) {
        fstream file;
        file.open(m_filename, ios::out | ios::trunc);
        if(!file) {
            cout << "Error creating file.." << endl;
        }
        file.close();
    }
}

void System::computeAnalytical(bool computeAnalytical) {
    m_computeAnalytical = computeAnalytical;
}

double System::quantumForce(int particleIndex, int dimensionIndex) {
    double quantumForce = -4 * m_waveFunction->getParameters()[0] * m_particles[particleIndex]->getPosition()[dimensionIndex];
    return quantumForce;
}

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */


    double stepProbability = 0;
    //double dx = 0;

    int randomParticle = m_random->nextInt(m_numberOfParticles-1); 
    int randomDirection = m_random->nextInt(m_numberOfDimensions-1);

    if (m_importanceSampling) {
        double oldWaveFunction = m_waveFunction->evaluate(m_particles);
        double qForceOld = quantumForce(randomParticle, randomDirection);
        double dx = m_random->nextGaussian(0, 1)*sqrt(m_dt) + 0.5*qForceOld*m_dt;
        m_particles[randomParticle]->adjustPosition(dx, randomDirection);

        double newWaveFunction = m_waveFunction->evaluate(m_particles);

        double qForceNew = quantumForce(randomParticle, randomDirection);

        double greensFunctionInner = 0.5*(qForceOld+qForceNew) * ((0.5*m_dt*0.5)*(qForceNew-qForceOld) - dx);                       //the stupid term here 0.5*m_dt*0.5 is because i want to make one of the 0.5's into m_D later i think
        stepProbability = exp(greensFunctionInner) * newWaveFunction*newWaveFunction / (oldWaveFunction*oldWaveFunction);
        double acceptStep = m_random->nextDouble();
        if (acceptStep < stepProbability){
            //cout << "hei" << endl; 
            return true;
        }
        else {
            m_particles[randomParticle]->adjustPosition(-dx, randomDirection);
            return false;
        }
    }
    else {
        double oldWaveFunction = m_waveFunction->evaluate(m_particles);
        //cout << "oldwave: " << oldWaveFunction << endl;
        double dx = m_stepLength * (2 * m_random->nextDouble() - 1);    //changed from nextdouble - 0.5
        m_particles[randomParticle]->adjustPosition(dx, randomDirection);
        double newWaveFunction = m_waveFunction->evaluate(m_particles);
        //cout << "newwave: " << newWaveFunction << endl;
        stepProbability = newWaveFunction*newWaveFunction / (oldWaveFunction*oldWaveFunction);
        double acceptStep = m_random->nextDouble();

        if (acceptStep < stepProbability){
            return true;
        }
        else {
            m_particles[randomParticle]->adjustPosition(-dx, randomDirection);
            return false;
        }
    }

    //cout << "got to the bottom of metro step" << endl;
    return false;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    //cout << "--------starting runMetropolisSteps--------" << endl;
    auto start = chrono::high_resolution_clock::now();
    m_particles = m_initialState->getParticles();

    //cout << "Size of m_particles: " <<static_cast<int>(m_particles.size()) << endl;

    m_sampler = new Sampler(this);
    m_numberOfMetropolisSteps = numberOfMetropolisSteps;
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);
    //cout << "Number of metropolis steps: " << m_numberOfMetropolisSteps << endl;

    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();
        //cout << acceptedStep << endl;

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */


        if (i > m_equilibrationFraction * m_numberOfMetropolisSteps) {
            //cout << "metrostep: "<<i<<",--------------"<< endl;
            m_sampler->sample(acceptedStep);
        }
        
    }

    auto stop = chrono::high_resolution_clock::now();
    auto simulationDuration = chrono::duration_cast<chrono::milliseconds>(stop - start);

    cout << "Simulation time: " << simulationDuration.count() << "ms, ";
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
    if (m_writeToFile){
        m_sampler->printEnergiesToFile();
    }
    //double energyGradient = m_sampler->getEnergyGradient();
    //cout << ", Energy Gradient: " <<energyGradient << endl;
}

