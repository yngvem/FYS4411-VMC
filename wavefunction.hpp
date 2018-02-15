#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include <cmath>
#include "particles.hpp"

struct WaveFunctionParameters {
    double a;
    double alpha;
    double beta;
    double omega_ho;
    double omega_z;
    double hbar; 
};

class SingleParticleFunction {
    WaveFunctionParameters* parameters;

public:
    double evaluate(Particles particles, int particle_idx);
    vector<double> evaluate_gradient(Particles particles, int particle_idx);
    double evaluate_laplace(Particles particles, int particle_idx);

    SingleParticleFunction ();
};

class WaveFunction {
    WaveFunctionParameters* parameters;

    double diff_function_u(Particles particles, int particle_idx_1,
                           int particle_idx_2);
    double first_deriv_u(Particles particles, int particle_idx_1,
                           int particle_idx_2);
    double second_deriv_u(Particles particles, int particle_idx_1,
                           int particle_idx_2);

public:

    double local_energy(Particles particles);
    double evaluate(Particles particles);
    double evaluate_PDF(Particles particle);


    WaveFunction ();
};


#endif // WAVEFUNCTION_HPP
