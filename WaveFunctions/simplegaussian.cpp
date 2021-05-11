#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>

using namespace std;

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
    WaveFunction(system){
        assert(alpha >= 0);
        m_numberOfParameters = 1;
        m_parameters.reserve(m_numberOfParameters);
        m_parameters.push_back(alpha);
        //m_parameters.push_back(beta);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
    double m_alpha = m_parameters[0];
    double rSquared = 0;
    const int numberOfParticles = m_system->getNumberOfParticles();
    const int numberOfDimensions = m_system->getNumberOfDimensions();

    for (int i=0; i < numberOfParticles; i++){
        std::vector<double> particlePosition = particles[i]->getPosition();
        for (int j=0; j < numberOfDimensions; j++){
            // check if we have a z-dimension as this is different than x and y with parameter beta varying the trap shape
            rSquared += particlePosition[j] * particlePosition[j];
        }
    }

    //cout << exp(-m_alpha * r_squared) << endl;
    //cout << "r squared: " << r_squared << endl;
    return exp(-m_alpha * rSquared);
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    // using analytical derivative for n particles with option to have it elliptical if beta is different from one
    // m_omega_ho is the frequency in xy direction, m_omega_r is the frequency in z direction, beta is the strechedness of the potential in the z direction, 
    double m_alpha = m_parameters[0];
    const int numberOfParticles = m_system->getNumberOfParticles();
    const int numberOfDimensions = m_system->getNumberOfDimensions();

    // Below is the calculation of the kinetic portion of the local energy
    double rSquared = 0;

    for (int i=0; i < numberOfParticles; i++){
        std::vector<double> particlePosition = particles[i]->getPosition();
        for (int d=0; d < numberOfDimensions; d++){
            //cout << d << endl;
            rSquared += particlePosition[d] * particlePosition[d];
            //cout << "no interaction, " <<"i:"<< i <<", d:"<< d << ", term:" << rSquared << endl;

        }
    }
    //cout << rSquared << endl;
    return (4 * m_alpha*m_alpha * rSquared - 2 * numberOfDimensions * numberOfParticles * m_alpha);
}

double SimpleGaussian::alphaDerivative(std::vector<class Particle*> particles) {
    const int numberOfParticles = m_system->getNumberOfParticles();
    const int numberOfDimensions = m_system->getNumberOfDimensions();

    double rSquared = 0;

    //cout << "hei" << endl;
    for (int i=0; i<numberOfParticles; i++){
        std::vector<double> particlePosition = particles[i]->getPosition();
        for (int d=0; d<numberOfDimensions; d++){
            rSquared += particlePosition[d] * particlePosition[d];
        }
    }
    return -rSquared;
}
