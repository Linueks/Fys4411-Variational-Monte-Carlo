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
    int numberOfParticles = 100;
    int numberOfSteps = 1 << 23;
    double omega_ho = 1.0;                          // Oscillator xy frequency.
    double omega_r = 1.0;                           // Oscillator z frequency  
    double stepLength = 1.0;                        // Metropolis step length.
    double equilibration = 0.0;                     // Amount of the total steps used
    double alpha = 0.5;                             // Variational parameter.
    double beta = 1;//2.82843;                      // Trap ellipticity
    double a = 0.0043;                              // Cutoff parameter for particle interaction
    double importanceDt = 0.001;                     // Timestep in importance sampling
    double derivativeStepSize = 1e-3;               // For when using numeric differentiation


    // different method flags, bad solution but better than comment / uncommenting
    bool calculateNumericDerivative = false;
    bool interactionOrNot = false;
    
    
    bool activateImportanceSampling = true;
    string filename = "..\\..\\data\\energy.txt";  
    bool writeToFile = true;

    // set up system with everything defined above
    
    System* system = new System();
    if (interactionOrNot){
        if (calculateNumericDerivative){
            system->setHamiltonian(new InteractionOscillatorNum(system, omega_ho, omega_r));
        }
        else{
            system->setHamiltonian(new InteractionOscillator(system, omega_ho, omega_r));
        }
        system->setWaveFunction(new Interacting(system, alpha, beta, a));
    }
    else{
        if (calculateNumericDerivative){
            system->setHamiltonian(new HarmonicOscillatorNum(system, omega_ho, omega_r)); //uncomment to run with numeric differentiation
        }
        else{
            system->setHamiltonian(new HarmonicOscillator(system, omega_ho, omega_r));
        }
        system->setWaveFunction(new SimpleGaussian(system, alpha));
    }

    system->setInteraction(interactionOrNot);
    system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles));
    system->setCutoffRadius(a);
    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setNumberOfMetropolisSteps(numberOfSteps);
    system->activateImportanceSampling(activateImportanceSampling, importanceDt);
    system->setDerivativeStepSize(derivativeStepSize);
    system->writeToFile(writeToFile, filename);
    system->runMetropolisSteps(numberOfSteps);
    
    
    std::vector<int> nParticles = {1, 10, 50};    
/*     for (int n : nParticles){
        for (int r=0; r<5; r++){
            System* system = new System();
            if (interactionOrNot){
                if (calculateNumericDerivative){
                    system->setHamiltonian(new InteractionOscillatorNum(system, omega_ho, omega_r));
                }
                else{
                    system->setHamiltonian(new InteractionOscillator(system, omega_ho, omega_r));
                }
                system->setWaveFunction(new Interacting(system, alpha, beta, a));
            }
            else{
                if (calculateNumericDerivative){
                    system->setHamiltonian(new HarmonicOscillatorNum(system, omega_ho, omega_r)); //uncomment to run with numeric differentiation
                }
                else{
                    system->setHamiltonian(new HarmonicOscillator(system, omega_ho, omega_r));
                }
                system->setWaveFunction(new SimpleGaussian(system, alpha));
            }

            system->setInteraction(interactionOrNot);
            system->setInitialState(new RandomUniform(system, numberOfDimensions, n));
            system->setCutoffRadius(a);
            system->setEquilibrationFraction(equilibration);
            system->setStepLength(stepLength);
            system->setNumberOfMetropolisSteps(numberOfSteps);
            system->activateImportanceSampling(activateImportanceSampling, importanceDt);
            system->setDerivativeStepSize(derivativeStepSize);
            system->writeToFile(writeToFile, filename);
            system->runMetropolisSteps(numberOfSteps);
        }
    } */





    // brute force search of the alpha parameter space
    double variationalStepLength = 0.05;            // step length for alpha and beta values
    double variationalIterations = 30;              // alpha + variationalStepLength*variationalIterations will be the end alpha
    

/*     for (int d=1; d<4; d++){
        for (int n : nParticles){
            alpha = 0.35;
            cout << "---------- d: " << d << "---------- n: " << n << endl;
            for (int i=0; i<variationalIterations; i++) {   
                System* system = new System(seed);
                //system->setHamiltonian(new HarmonicOscillatorNum(system, omega_ho, omega_r)); //uncomment to run with numeric differentiation
                system->setHamiltonian(new HarmonicOscillator(system, omega_ho, omega_r));
                system->setWaveFunction(new SimpleGaussian(system, alpha));
                system->setInitialState(new RandomUniform(system, d, n));
                system->setEquilibrationFraction(equilibration);
                system->setStepLength(stepLength);
                system->setNumberOfMetropolisSteps(numberOfSteps);
                system->activateImportanceSampling(activateImportanceSampling, importanceDt);        
                system->setDerivativeStepSize(derivativeStepSize);
                system->runMetropolisSteps(numberOfSteps);
                alpha += variationalStepLength;
            }
        }
    } */
    




    // simplest gradient descent search of the alpha parameter space

    // for 1D 1P learning rate=0.35 and 10000 metrosteps converges nicely
    // 1D 10P 0.05, 1D 100P 0.3, 1D 500 0.3

    // for 2D 1P learning rate=0.35 and 10000 metrosteps converges nicely
    // 2D 10P 0.3, 2D 100P 0.05 2D 500P 0.1 

    // for 3D 1P learning rate=0.15 and 10000 metrosteps converges nicely
    // 3D 10P 0.05, 3D 100P 0.05, 3D 500P 0.15


    double nGradientDescentRuns = 50;
    double learningRate = 0.15;

/*     for (int i=0; i<nGradientDescentRuns; i++) {
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
    } */


    return 0;
}
