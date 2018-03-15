#include <iostream>
#include <cmath>
#include <random>
#include "wavefunction.hpp"
#include "particles.hpp"
#include "particle.hpp"
#include "metropolis_sampler.hpp"

using namespace std;

MetropolisSampler::MetropolisSampler (double step_size_,
                                        ParticlesParams particles_params,
                                        WaveFunctionParameters wavefunction_params) :
    BaseSampler(particles_params, wavefunction_params),
    step_size(step_size_)
{}

void MetropolisSampler::propose_pertubation(){
    current_particle_idx = (current_particle_idx + 1) % particles.num_particles;
    //Propose perturbation
    vector<double> perturbation(particles.num_dimensions);
    for (int i = 0; i< particles.num_dimensions; ++i)
        perturbation[i] = uniform(generator);

    particles.perturb_particle(current_particle_idx,perturbation);
}

void MetropolisSampler::reject_perturbation(){
    particles.reject_perturbation(current_particle_idx);
}

double MetropolisSampler::rejction_criteria(){
    double ratio;
    ration = probability_value/wave_function.evaluate_PDF(particles);
}

