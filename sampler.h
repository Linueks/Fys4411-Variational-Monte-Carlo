#pragma once
#include <vector>

class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();
    void printEnergiesToFile();
    void computeAverages();
    void resetValues();
    double getEnergy() { return m_energy; }
    double getEnergyGradient();

private:
    int m_numberOfMetropolisSteps = 0;
    int m_stepNumber = 0;
    int m_acceptedSteps = 0;
    double m_energy = 0;
    double m_cumulativeEnergy = 0;
    std::vector<double> m_energies = std::vector<double>();
    double m_deltaPsi = 0;
    double m_psiAlphaDerivative = 0;
    double m_psiEnergyDerivative = 0;
    double m_derivative = 0;
    double m_acceptedRatio = 0;
    class System* m_system = nullptr;
};
