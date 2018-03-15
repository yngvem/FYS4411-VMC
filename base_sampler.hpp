#ifndef SAMPLER_BASE
#define SAMPLER_BASE

#include <iostream>
#include <cmath>
#include <random>
#include "wavefunction.hpp"
#include "particles.hpp"
#include "particle.hpp"

using namespace std;

class BaseSampler {

protected:
    Particles particles;
    WaveFunction wave_function;
    double probability_value;
    mt19937 generator;
    uniform_real_distribution<double> uniform;

    virtual void propose_pertubation() = 0;
    virtual void reject_perturbation() = 0;
    virtual double rejction_criteria() = 0;

    double compute_probability();
    double local_energy();

    void MC_step();

public:
    void warm_up(int num_steps);
    vector<double> perform_iterations(int num_steps, int memory_frequency);
    vector<double> compute_local_energy(int num_steps);

    BaseSampler (ParticlesParams particles_params,
                 WaveFunctionParameters wavefunction_params,
                 int seed);
};

#endif
