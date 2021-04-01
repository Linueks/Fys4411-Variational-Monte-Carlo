#include <iostream>
#include "system.h"
#include "particle.h"
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/interacting.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/interactionoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"

//#include <iostream>
#include <fstream>
#include <string>

using namespace std;


int main() {

    ifstream myFile_Handler;
    string myLine;

    // File Open in the Read Mode
/*     myFile_Handler.open("..\\input_variables.txt");
    if(myFile_Handler.is_open())
    {
        // Keep reading the file
        while(getline(myFile_Handler, myLine))
        {
            // print the line on the standard output
            cout << myLine[0] << endl;
        }
    // File Close
    myFile_Handler.close();
    }
    else
    {
        cout << "Unable to open the file!";
    }
    return 0; */



    // Seed for the random number generator
    int seed = 2021;

    int numberOfDimensions = 1;
    int numberOfParticles = 1;
    int numberOfSteps = 1 << 6;
    double omega_ho = 1.0;                          // Oscillator xy frequency.
    double omega_r = 1.0;                           // Oscillator z frequency  
    double stepLength = 0.1;                        // Metropolis step length.
    double equilibration = 0.1;                     // Amount of the total steps used
    // for equilibration.

    double alpha = 0.5;                             // Variational parameter.
    double beta = 2.82843;                                // Trap ellipticity
    double a = 0.0043;                                 // cutoff parameter for interaction

    // for when using HarmonicOscillatorNum class
    double derivativeStepSize = 1e-3;

    //bool computeAnalytical = false;
    bool activateImportanceSampling = true;
    double importanceDt = 0.05;

    bool writeToFile = true;
    string filename = "..\\..\\data\\energy.txt";


    // different method flags, bad solution but better than comment / uncommenting
    bool calculateNumericDerivative = false;
    bool simpleHarmonicOscillator = false;
    bool interactingHarmonicOscillator = true;



    System* system = new System(seed);

    if (simpleHarmonicOscillator){
        if (calculateNumericDerivative){
            system->setHamiltonian(new HarmonicOscillatorNum(system, omega_ho, omega_r)); //uncomment to run with numeric differentiation
        }
        else{
            system->setHamiltonian(new HarmonicOscillator(system, omega_ho, omega_r));
        }
        system->setWaveFunction(new SimpleGaussian(system, alpha));
    }


    if (interactingHarmonicOscillator){
        if (calculateNumericDerivative){
            system->setHamiltonian(new InteractionOscillatorNum(system, omega_ho, omega_r));
        }
        else{
            system->setHamiltonian(new InteractionOscillator(system, omega_ho, omega_r));
        }
        system->setWaveFunction(new Interacting(system, alpha, beta, a));
    }


    system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setCutoffRadius(a);
    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setNumberOfMetropolisSteps(numberOfSteps);
    system->activateImportanceSampling(activateImportanceSampling, importanceDt);
    system->setDerivativeStepSize(derivativeStepSize);
    system->writeToFile(writeToFile, filename);
    system->runMetropolisSteps(numberOfSteps);



/*     double variationalStepLength = 0.05;            // step length for alpha and beta values
    double variationalIterations = 5;              // alpha + variationalStepLength*variationalIterations will be the end alpha
    for (int i=0; i<variationalIterations; i++) {   
        System* system = new System(seed);
        //system->setHamiltonian(new HarmonicOscillatorNum(system, omega_ho, omega_r)); //uncomment to run with numeric differentiation
        system->setHamiltonian(new HarmonicOscillator(system, omega_ho, omega_r));
        system->setWaveFunction(new SimpleGaussian(system, alpha));
        system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction(equilibration);
        system->setStepLength(stepLength);
        system->setNumberOfMetropolisSteps(numberOfSteps);
        system->activateImportanceSampling(activateImportanceSampling, importanceDt);        
        system->setDerivativeStepSize(derivativeStepSize);
        system->runMetropolisSteps(numberOfSteps);
        alpha += variationalStepLength;
    } */



/*     double nGradientDescentRuns = 50;
    double learningRate = 0.35;

    for (int i=0; i<nGradientDescentRuns; i++) {
        System* system = new System(seed);
        //system->setHamiltonian(new HarmonicOscillatorNum(system, omega_ho, omega_r)); //uncomment to run with numeric differentiation
        system->setHamiltonian(new HarmonicOscillator(system, omega_ho, omega_r));
        system->setWaveFunction(new SimpleGaussian(system, alpha));
        system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
        system->setEquilibrationFraction(equilibration);
        system->setStepLength(stepLength);
        system->setNumberOfMetropolisSteps(numberOfSteps);
        system->activateImportanceSampling(activateImportanceSampling, importanceDt);        
        system->setDerivativeStepSize(derivativeStepSize);
        system->writeToFile(writeToFile, filename);
        system->runMetropolisSteps(numberOfSteps);
        
        Sampler* sampler = system->getSampler();
        double energyGradient = sampler->getEnergyGradient();
        alpha -= learningRate*energyGradient;
    }
 */

    return 0;
}
