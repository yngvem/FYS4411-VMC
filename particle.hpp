#ifndef PARTICLES
#define PARTICLES

#include <vector>

using namespace std;

class Particle {
    int num_dimensions;
    vector<double> r;
    double m;

    void initiate_r();

public:
    Particle(int num_dimensions, double mass);
    Particle(int num_dimensions);

    // Obtaining information about the position.
    const vector<double>& position() const;
    const double norm() const;
    const double& operator[](int i); // TODO: Should this be const and reference?

    // Arithmetic operations between multiple particles.
    Particle operator+(const Particle& particle);
    Particle operator-(const Particle& particle);
    Particle operator*(const Particle& particle);
    
    // Arithmetic operations between a particle and a number.
    Particle operator+(double a);
    Particle operator*(double a);
    Particle operator-(double a);
    Particle operator/(double a);
};

#endif