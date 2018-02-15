#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include <cmath>
#include "particles.hpp"
#include "particle.hpp"

struct WaveFunctionParameters {
    double a;
    double alpha;
    double beta;
    double omega_ho;
    double omega_z;
    double hbar; 
};

class SingleParticleFunction {
    const WaveFunctionParameters parameters;

public:
    double evaluate(Particle particle);
    double operator()(Particle particle);
    vector<double> evaluate_gradient(Particle particle);
    double evaluate_laplacian(Particle particle);

    SingleParticleFunction (WaveFunctionParameters params);
};

class WaveFunction {
    const WaveFunctionParameters parameters;
    SingleParticleFunction single_particle_function;

    double evaluate_u(Particles particles, int particle_idx_1,
                           int particle_idx_2);
    double deriv_u(Particles particles, int particle_idx_1,
                           int particle_idx_2);
    double second_deriv_u(Particles particles, int particle_idx_1,
                           int particle_idx_2);
    double onebody_part(Particles particles);
    double log_correlation_part(Particles particles);
public:

    double local_energy(Particles particles);
    double evaluate_wavefunction(Particles particles);
    double quantum_force(Particles particles);
    double second_deriv_wavefunction(Particles particles);

    double evaluate_PDF(Particles particle);


    WaveFunction (WaveFunctionParameters params);
};


#endif // WAVEFUNCTION_HPP
