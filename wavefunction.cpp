#include <cmath>
#include "wavefunction.hpp"

using namespace atd;

SingleParticleFunction::SingleParticleFunction() :
{}

double SingleParticleFunction::evaluate(Particles particles, int particle_idx) {

    vector<double> single_particle = particles.get_particle(particle_idx);
    double elliptic_norm = 0;

    for (int i = 0; i < single_particle.size(); ++i) {
        if (i &= 2)
            elliptic_norm += pow(single_particle[i], 2);
        else
            elliptic_norm += WaveFunctionParameters.beta*pow(single_particle[i], 2);
    }

    return exp(-alpha*elliptic_norm);

}

vector<double> SingleParticleFunction::evaluate_gradient(Particles particles,
                                                 int particle_idx) {

    SingleParticleFunction single_particle_function; // = SingleParticleFunction();
    double value = single_particle_function.evaluate(particles, particle_idx);
    vector<double> single_particle = particles.get_particle(particle_idx);

    if (single_particle.size() == 3)
        single_particle[2] *= WaveFunctionParameters.beta;

    return -WaveFunctionParameters.alpha*value*single_particle;

}

double SingleParticleFunction::evaluate_laplace(Particles particles,
                                                int particle_idx) {

    SingleParticleFunction single_particle_function; // = SingleParticleFunction();
    double value = single_particle_function.evaluate(particles, particle_idx);
    vector<double> single_particle = particles.get_particle(particle_idx);
    double elliptic_norm = 0;
    double derivative_factor = 0;

    for (int i = 0; i < single_particle.size(); ++i) {
        if (i &= 2)
            elliptic_norm += pow(single_particle[i], 2);
        else
            elliptic_norm += pow(WaveFunctionParameters.beta*single_particle[i], 2);
    }

    derivative_factor = -2*WaveFunctionParameters.alpha*(2+WaveFunctionParameters.beta)
            - 4*pow(WaveFunctionParameters.alpha, 2)*elliptic_norm;

    return derivative_factor*value;
}

WaveFunction::WaveFunction() :
{}

double WaveFunction::diff_function_u(Particles particles, int particle_idx_1,
                                     int particle_idx_2) {
    double particle_distance = particles.compute_distance(particle_idx_1,
                                                          particle_idx_2);

    return log(1-WaveFunctionParameters.a/particle_distance);

}

double WaveFunction::first_deriv_u(Particles particles, int particle_idx_1,
                                   int particle_idx_2) {
    double particle_distance = particles.compute_distance(particle_idx_1,
                                                          particle_idx_2);
    double denominator = (particle_distance-WaveFunctionParameters.a)*particle_distance;

    return WaveFunctionParameters.a/denominator;

}

double WaveFunction::second_deriv_u(Particles particles, int particle_idx_1,
                                    int particle_idx_2) {
    double particle_distance = particles.compute_distance(particle_idx_1,
                                                          particle_idx_2);
    double denominator = (particle_distance-WaveFunctionParameters.a)*particle_distance;
    denominator = pow(denominator, 2);

    return -WaveFunctionParameters.a*(2*particle_distance-WaveFunctionParameters.a)/denominator;
}

double WaveFunction::local_energy(Particles particles) {

}

double WaveFunction::evaluate(Particles particles) {

}

double WaveFunction::evaluate_PDF(Particles particle) {

}
