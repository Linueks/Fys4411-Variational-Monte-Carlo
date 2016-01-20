#pragma once

class Particle {
public:
    Particle();
    void setPosition(double* position);
    void adjustPosition(double change, int dimension);
    void setNumberOfDimensions(int numberOfDimensions);
    double* getPosition() { return m_position; }

private:
    int     m_numberOfDimensions = 0;
    double* m_position = nullptr;
};

