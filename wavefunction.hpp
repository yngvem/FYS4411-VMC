#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include "particles.hpp"

struct WaveFunctionParameters {
    double a;
    double alpha;
    double beta;
    double omega_ho;
    double omega_z;
    double hbar;
}

class SingleParticleFunction{
    WaveFunctionParameters* parameters;

public:
    double evaluate(Particles particles, int particle_idx);
    double evaluate_gradient(Particles particles, int particle_idx);
    double evaluate_laplace(Particles particles, int particle_idx);

    SingleParticleFunction (& parameters_)
}

class WaveFunction {
    WaveFunctionParameters* parameters;

public:
    double local_energy(Particles particles);
    double evaluate(Particles particles);
    double evaluate_PDF(Particles particle);

};


#endif // WAVEFUNCTION_HPP
