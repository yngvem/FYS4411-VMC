#ifndef SAMPLER_BASE
#define SAMPLER_BASE

#include <iostream>
#include <cmath>
#include "wavefunction.hpp"
#include "particles.hpp"
#include "particle.hpp"

using namespace std;

class SamplerBase {

protected:
    Particles particles;
    WaveFunction wave_function;
    double probability_value;

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

    SamplerBase (ParticlesParams particles_params,
                 WaveFunctionParameters wavefunction_params);
};

#endif
