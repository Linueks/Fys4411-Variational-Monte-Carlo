#include "interacting.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>
//#include <math.h>
using namespace std;


Interacting::Interacting(System* system, double alpha, double beta, double a) :
    WaveFunction(system){
    assert(alpha >= 0);
        m_numberOfParameters = 3;
        m_parameters.reserve(m_numberOfParameters);
        m_parameters.push_back(alpha);
        m_parameters.push_back(beta);
        m_parameters.push_back(a);
}

double Interacting::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
    double m_alpha = m_parameters[0];
    double m_beta = m_parameters[1];
    double m_a = m_parameters[2];
    double rSquared = 0;
    const int numberOfParticles = m_system->getNumberOfParticles();
    const int numberOfDimensions = m_system->getNumberOfDimensions();
    double absoluteDiff = 0;
    double correlationFactor = 1;

    for (int i=0; i < numberOfParticles; i++){
        std::vector<double> iParticle = particles[i]->getPosition();
        for (int j=0; j < numberOfDimensions; j++){
            if (j==2){
                rSquared += iParticle[j]*iParticle[j]*m_beta;
            }
            else{
                rSquared += iParticle[j]*iParticle[j];
            }
        }
        for (int j=i+1; j < numberOfParticles; j++){
            std::vector<double> jParticle = particles[j]->getPosition();
            double diff = 0;
            for (int d=0; d < numberOfDimensions; d++){
                diff += iParticle[d]*iParticle[d] - jParticle[d]*jParticle[d];
            }
            absoluteDiff = sqrt(diff);
            if (absoluteDiff <= m_a){
                correlationFactor *= 1;
            }
            else {
                correlationFactor *= 1 - m_a / absoluteDiff;
            }
        }
    }
    return exp(m_alpha * rSquared) * correlationFactor;
}

double Interacting::computeDoubleDerivative(std::vector<class Particle*> particles) {
    double m_alpha = m_parameters[0];
    double m_beta = m_parameters[1];
    double m_a = m_parameters[2];    
    const int numberOfParticles = m_system->getNumberOfParticles();
    const int numberOfDimensions = m_system->getNumberOfDimensions();
    double mainTerm1 = 0;
    double mainTerm2 = 0;
    double mainTerm3 = 0;
    double mainTerm4 = 0;
    double mainTerm5 = 0;
      
    
    for (int i=0; i < numberOfParticles; i++){
        double temp = 0;
        std::vector<double> iParticle = particles[i]->getPosition();
        for (int d=0; d < numberOfDimensions; d++){
            if (d==2){
                mainTerm1 += iParticle[d]*iParticle[d]*m_beta;
                cout << "interaction, " <<"i:"<< i <<", d:"<< d << ", term:" << mainTerm1 << endl;
            }
            else{
                mainTerm1 += iParticle[d]*iParticle[d];
                cout << "interaction, " <<"i:"<< i <<", d:"<< d << ", term:" << mainTerm1 << endl;
            }
/*             mainTerm1 += iParticle[d]*iParticle[d];
            cout << "interaction, " <<"i:"<< i <<", d:"<< d << ", term:" << mainTerm1 << endl; */
        }
        for (int j=i+1; j < numberOfParticles; j++){
            double temp = 0;
            std::vector<double> jParticle = particles[j]->getPosition();
            for (int h=0; h < numberOfDimensions; h++){
                cout << "kom inn i 2" << endl;
                if (h==2){
                    mainTerm2 += iParticle[h] * (iParticle[h] - jParticle[h]) *m_beta;
                }
                else{
                    mainTerm2 += iParticle[h] * (iParticle[h] - jParticle[h]);
                }
                
            }
            temp *= uPrimeR(iParticle, jParticle);
            mainTerm2 += temp;

            for (int k=j+1; k < numberOfParticles; k++){
                cout << "kom inn i 3" << endl;
                double temp = 0;
                std::vector<double> kParticle = particles[k]->getPosition();
                for (int u=0; u < numberOfDimensions; u++){
                    temp += (iParticle[u] - jParticle[u]) * (iParticle[u] - kParticle[u]);
                }
                temp *= uPrimeR(iParticle, kParticle) * uPrimeR(jParticle, iParticle);
                mainTerm3 += temp;
            }
        mainTerm4 += 2*uPrimeR(jParticle, iParticle);
        mainTerm5 += uDoublePrimeR(jParticle, iParticle);
        }   
    }
    //std::cout << "1: "<<mainTerm1<<", 2: "<<mainTerm2<<", 3: "<<mainTerm3<<", 4: "<<mainTerm4<<", 5: "<<mainTerm5 << endl;
    return 4*m_alpha*mainTerm1 - 2*m_alpha*(2 + m_beta) - 4*m_alpha*mainTerm2 + mainTerm3 + mainTerm4 + mainTerm5;    
}

double Interacting::uPrimeR(std::vector<double> ithParticle, std::vector<double> jthParticle) {
    double positionDifference = 0;
    double norm = 0;
    double uPrime = 0;
    double m_a = m_parameters[2];    
    const int numberOfDimensions = m_system->getNumberOfDimensions();

    for (int d=0; d < numberOfDimensions; d++){
        positionDifference += (ithParticle[d] - jthParticle[d]) * (ithParticle[d] - jthParticle[d]);
    }
    // taking the norm of the difference
    norm = sqrt(positionDifference);
    uPrime = m_a / (norm*norm - m_a*norm);

    return uPrime / norm;
}

double Interacting::uDoublePrimeR(std::vector<double> ithParticle, std::vector<double> kthParticle) {
    double positionDifference = 0;
    double norm = 0;
    double uDoublePrime = 0;
    double m_a = m_parameters[2];    
    const int numberOfDimensions = m_system->getNumberOfDimensions();

    for (int d=0; d < numberOfDimensions; d++){
        positionDifference += (ithParticle[d] - kthParticle[d]) * (ithParticle[d] - kthParticle[d]);
    }
    norm = sqrt(positionDifference);
    uDoublePrime = -m_a*(2*norm - m_a) / ((norm*norm - m_a*norm)*(norm*norm - m_a*norm));
    return uDoublePrime;
}


double Interacting::alphaDerivative(std::vector<class Particle*> particles) {
    return 0;
}