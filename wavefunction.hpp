#ifndef WAVEFUNCTION_HPP
#define WAVEFUNCTION_HPP

#include <cmath>
#include <limits>
#include "particles.hpp"
#include "particle.hpp"

struct WaveFunctionParameters {
    double a;
    double alpha;
    double beta;
    double omega_ho;
    double omega_z;
    double hbar;

    WaveFunctionParameters();
    void set_default_parameters();
};

class SingleParticleFunction {
    WaveFunctionParameters parameters;

public:
    double evaluate(const Particle& particle);
    double operator()(const Particle& particle);
    vector<double> evaluate_gradient(const Particle& particle);
    double evaluate_laplacian(const Particle& particle);

    SingleParticleFunction (WaveFunctionParameters params);
    SingleParticleFunction ();
    void set_parameters (WaveFunctionParameters params);
};

class WaveFunction {
    WaveFunctionParameters parameters;
    SingleParticleFunction single_particle_function;

    double evaluate_u (const Particles& particles, int particle_idx_1,
                           int particle_idx_2);
    double deriv_u (const Particles& particles, int particle_idx_1,
                           int particle_idx_2);
    double second_deriv_u (const Particles& particles, int particle_idx_1,
                           int particle_idx_2);
    double onebody_part (const Particles& particles);
    double log_correlation_part (const Particles& particles);
public:

    double local_energy (const Particles& particles);
    double evaluate_wavefunction (const Particles& particles);
    vector<vector <double> > quantum_force (const Particles& particles);
    double kinetic_energy (const Particles& particles);
    double ext_potential (const Particles& particles);
    double int_potential (const Particles& particles);

    double evaluate_PDF (const Particles& particle);


    WaveFunction (WaveFunctionParameters params);
    WaveFunction ();

    void set_parameters (WaveFunctionParameters params);
    void set_default_parameters ();
};


#endif // WAVEFUNCTION_HPP
