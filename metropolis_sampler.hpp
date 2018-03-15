#ifndef METROPOLIS_SAMPLER_HPP
#define METROPOLIS_SAMPLER_HPP

#include <random>
#include "base_sampler.hpp"


class MetropolisSampler : public BaseSampler
{
    double step_size;
    void propose_pertubation();
    void reject_perturbation();
    double rejction_criteria();
public:
    MetropolisSampler (double step_size_,
                       ParticlesParams particles_params,
                       WaveFunctionParameters wavefunction_params,
                       int seed);
};

#endif // METROPOLIS_SAMPLER_HPP
