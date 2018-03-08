#include <iostream>
#include <cmath>
#include "wavefunction.hpp"
#include "particles.hpp"
#include "particle.hpp"
#include "sampler_base.hpp"

using namespace std;

WaveFunctionParameters params;

SamplerBase::SamplerBase (ParticlesParams particles_params,
                          WaveFunctionParameters wavefunction_params) :
    particles(Particles(particles_params)),
    wave_function(WaveFunction(wavefunction_params)),
    probability_value(wave_function.evaluate_PDF(particles))
{}

double SamplerBase::compute_probability(){
    return wave_function.evaluate_PDF(particles);
}

double SamplerBase::local_energy(){
    return wave_function.local_energy(particles);
}

/*
void SamplerBase::MC_step () {
    propose_pertubation();
    double potential_probability = compute_probability();

    if (r > rejction_criteria())
        probability_value = potential_probability;
    else
        reject_perturbation();
}
*/

void SamplerBase::warm_up (int num_steps) {
    for (int i = 0; i < num_steps; ++i)
        MC_step();
}


vector<double> SamplerBase::perform_iterations(int num_steps, int memory_frequency) {
    int num_save_steps = num_steps/memory_frequency;
    vector<double> energies(num_save_steps);
    for (int i = 0; i < num_steps; ++i) {
        MC_step();

        if (i % memory_frequency == 0)
            energies[i/memory_frequency] = local_energy();
    }

    return energies;
}

vector<double> SamplerBase::compute_local_energy(int num_steps) {
    double mean_energy;
    double mean_squared_energy;

    for (int i = 0; i < num_steps; ++i) {
        MC_step();

        mean_energy += local_energy();
        mean_squared_energy += pow(mean_energy, 2);
    }

    mean_energy /= num_steps;
    mean_squared_energy /= num_steps;

    vector<double> energy_vector(2);
    energy_vector[0] = mean_energy;
    energy_vector[1] = mean_squared_energy - pow(mean_energy, 2);

    return energy_vector;
}

