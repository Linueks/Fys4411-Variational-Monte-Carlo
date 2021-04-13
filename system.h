#pragma once
#include <vector>
#include <Math/random.h>

class System {
public:
    System();
    System(int seed);
    bool metropolisStep();
    void runMetropolisSteps(int numberOfMetropolisSteps);
    double quantumForce(int particleIndex, int dimensionIndex);
    void iterateVariational(int maxVariations);
    void setNumberOfParticles(int numberOfParticles);
    void setNumberOfDimensions(int numberOfDimensions);
    void setStepLength(double stepLength);
    void setCutoffRadius(double a);
    void setInteraction(bool interactionOrNoInteraction);
    void setNumberOfMetropolisSteps(int numberOfSteps);
    void setEquilibrationFraction(double equilibrationFraction);
    
    
    //void setVariationalStepLength(double variationStep);
    //void setAlphaStart(double alphaStart);
    
    
    void setDerivativeStepSize(double derivativeStep);
    void activateImportanceSampling(bool activate, double dt);
    void writeToFile(bool writeToFile, std::string filename);
    void setHamiltonian(class Hamiltonian* hamiltonian);
    void setWaveFunction(class WaveFunction* waveFunction);
    void setInitialState(class InitialState* initialState);
    
    void computeAnalytical(bool computeAnalytical);
    
    class WaveFunction* getWaveFunction() {return m_waveFunction;}
    class Hamiltonian* getHamiltonian() {return m_hamiltonian;}
    class Sampler* getSampler() {return m_sampler;}
    std::vector<class Particle*> getParticles() {return m_particles;}
    class Random* getRandomEngine() {return m_random;}
    int getNumberOfParticles() {return m_numberOfParticles;}
    int getNumberOfDimensions() {return m_numberOfDimensions;}
    int getNumberOfMetropolisSteps() {return m_numberOfMetropolisSteps;}
    double getDerivativeStepSize() {return m_derivativeStepSize;}
    double getEquilibrationFraction() {return m_equilibrationFraction;}
    double getCutoffRadius() {return m_a;}
    bool getWriteToFile() {return m_writeToFile;}
    bool getInteractionOrNot() {return m_interactionOrNot;}
    std::string getOutfileName() {return m_filename;}

    bool m_computeAnalytical = false;

private:
    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;
    int m_numberOfMetropolisSteps = 0;
    double m_equilibrationFraction = 0.0;
    double m_stepLength = 0.1;
    double m_variationStep = 0.0;
    double m_alphaStart = 0.4;
    double m_derivativeStepSize = 0;
    double m_dt = 0;
    double m_a = 0;
    bool m_importanceSampling = false;
    bool m_writeToFile = false;
    bool m_interactionOrNot = false;
    std::string m_filename = "";

    
    
    std::vector<double> m_alphaValues = std::vector<double>();
    class WaveFunction* m_waveFunction = nullptr;
    class Hamiltonian* m_hamiltonian = nullptr;
    class InitialState* m_initialState = nullptr;
    class Sampler* m_sampler = nullptr;
    std::vector<class Particle*> m_particles = std::vector<class Particle*>();
    class Random* m_random = nullptr;
};

