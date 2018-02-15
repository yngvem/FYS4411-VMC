#include <cmath>
#include "wavefunction.hpp"
#include "particle.hpp"

using namespace std;

SingleParticleFunction::SingleParticleFunction (WaveFunctionParameters params) :
    parameters(params)
{}

double SingleParticleFunction::operator() (Particle particle) {
    return evaluate(particle);
}

double SingleParticleFunction::evaluate (Particle particle) {

    double elliptic_norm = 0;
    for (int i = 0; i< particle.num_dimensions; ++i) {
        if (i == 2)
            elliptic_norm += parameters.beta*pow(particle[i], 2);
        else
            elliptic_norm += pow(particle[i], 2);
    }

    return exp(-parameters.alpha*elliptic_norm);

}

vector<double> SingleParticleFunction::evaluate_gradient (Particle particle) {
    double value = evaluate(particle);
    vector<double> gradient(particle.num_dimensions, -2*parameters.alpha*value);

    for (int i = 0; i < particle.num_dimensions; ++i) {
        if (i == 2)
            gradient[i] *= parameters.beta * particle[i];
        else
            gradient[i] *= particle[i];
    }

    return gradient;
}

double SingleParticleFunction::evaluate_laplacian (Particle particle) {
    double value = evaluate(particle);
    double elliptic_norm = 0;
    double derivative_factor = 0;

    for (int i = 0; i < particle.num_dimensions; ++i) {
        if (i != 2)
            elliptic_norm += pow(particle[i], 2);
        else
            elliptic_norm += pow(parameters.beta*particle[i], 2);
    }

    derivative_factor = -2*parameters.alpha*(2+parameters.beta)
                        -4*pow(parameters.alpha, 2)*elliptic_norm;

    return derivative_factor*value;
}

WaveFunction::WaveFunction (WaveFunctionParameters params) :
    parameters(params)
{
    SingleParticleFunction single_particle_function(params);
}

double WaveFunction::evaluate_u (Particles particles, int particle_idx_1,
                                     int particle_idx_2) {
    double particle_distance = particles.compute_distance(particle_idx_1,
                                                          particle_idx_2);

    return log(1-parameters.a/particle_distance);

}

double WaveFunction::deriv_u (Particles particles, int particle_idx_1,
                                   int particle_idx_2) {
    double particle_distance = particles.compute_distance(particle_idx_1,
                                                          particle_idx_2);
    double denominator = (particle_distance-parameters.a)*particle_distance;

    return parameters.a/denominator;

}

double WaveFunction::second_deriv_u (Particles particles, int particle_idx_1,
                                    int particle_idx_2) {
    double particle_distance = particles.compute_distance(particle_idx_1,
                                                          particle_idx_2);
    double denominator = (particle_distance-parameters.a)*particle_distance;
    denominator = pow(denominator, 2);

    return -parameters.a*(2*particle_distance-parameters.a)/denominator;
}

double WaveFunction::onebody_part (Particles particles) {
    double g = 1;
    for (int i = 0; i < particles.num_particles; ++i) 
        g *= single_particle_function(particles.get_particle(i));
    
    return g;
}

double WaveFunction::local_energy (Particles particles) {

}

double WaveFunction::log_correlation_part (Particles particles) {
    double u = 0;
    for (int i = 0; i < particles.num_particles-1; ++i)
        for (int j = i+1; j < particles.num_particles; ++j)
            u += evaluate_u(particles, i, j);
}

double WaveFunction::evaluate_wavefunction (Particles particles) {
    return onebody_part(particles) * exp(log_correlation_part(particles));
}

double WaveFunction::quantum_force (Particles particles) {
    
}

double WaveFunction::second_deriv_wavefunction (Particles particles) {

}

double WaveFunction::evaluate_PDF (Particles particle) {

}
