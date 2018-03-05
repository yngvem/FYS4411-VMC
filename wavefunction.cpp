#include <cmath>
#include <limits>
#include "wavefunction.hpp"
#include "particle.hpp"
#include "particles.hpp"

using namespace std;

//----------------Parameters------------------//
WaveFunctionParameters::WaveFunctionParameters () :
    a(0),
    alpha(1),
    beta(1),
    omega_ho(1),
    omega_z(1),
    hbar(1)
{}
void WaveFunctionParameters::set_default_parameters () {
    a = 0;
    alpha = 1;
    beta = 1;
    omega_ho = 1;
    omega_z = 1;
    hbar = 1;
}

//----------------Single particle function-------------------//
SingleParticleFunction::SingleParticleFunction (WaveFunctionParameters params) :
    parameters(params)
{}


SingleParticleFunction::SingleParticleFunction () {};


void SingleParticleFunction::set_parameters (WaveFunctionParameters params) {
    parameters = params;
}


double SingleParticleFunction::operator() (const Particle& particle) {
    return evaluate(particle);
}


double SingleParticleFunction::evaluate (const Particle& particle) {

    double elliptic_norm = 0;
    for (int i = 0; i< particle.num_dimensions; ++i) {
        if (i == 2)
            elliptic_norm += parameters.beta*pow(particle[i], 2);
        else
            elliptic_norm += pow(particle[i], 2);
    }

    return exp(-parameters.alpha*elliptic_norm);

}


vector<double> SingleParticleFunction::evaluate_gradient (const Particle& particle) {
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


double SingleParticleFunction::evaluate_laplacian (const Particle& particle) {
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


//----------------Wavefunction-------------------//
WaveFunction::WaveFunction (WaveFunctionParameters params) :
    parameters(params)
{
    single_particle_function.set_parameters(params);
}

WaveFunction::WaveFunction () {}


double WaveFunction::evaluate_u (const Particles& particles, int particle_idx_1,
                                     int particle_idx_2) {
    double particle_distance = particles.compute_distance(particle_idx_1,
                                                          particle_idx_2);

    return log(1-parameters.a/particle_distance);

}


double WaveFunction::deriv_u (const Particles& particles, int particle_idx_1,
                                   int particle_idx_2) {
    double particle_distance = particles.compute_distance(particle_idx_1,
                                                          particle_idx_2);
    double denominator = (particle_distance-parameters.a)*particle_distance;

    return parameters.a/denominator;

}


double WaveFunction::second_deriv_u (const Particles& particles, int particle_idx_1,
                                    int particle_idx_2) {
    double particle_distance = particles.compute_distance(particle_idx_1,
                                                          particle_idx_2);
    double denominator = (particle_distance-parameters.a)*particle_distance;
    denominator = pow(denominator, 2);

    return -parameters.a*(2*particle_distance-parameters.a)/denominator;
}


double WaveFunction::onebody_part (const Particles& particles) {
    double g = 1;
    for (int i = 0; i < particles.num_particles; ++i) 
        g *= single_particle_function(particles.get_particle(i));  // We get error because we have not fixed get_particle yet
    
    return g;
}


double WaveFunction::local_energy (const Particles& particles) {

    return second_deriv_wavefunction_quotient(particles) + ext_potential(particles)
          +int_potential(particles);

}


double WaveFunction::ext_potential(const Particles& particles){

    vector<double> weights = {parameters.omega_ho, parameters.omega_ho, parameters.omega_z};
    double external_potential = 0;

    for (int k = 0; k < particles.num_particles; ++k){
        Particle particle = particles.get_particle(k);

        external_potential += 0.5*particle.m*pow(particle.weighted_distance_from_origin(weights), 2);
    }

    return external_potential;
}

double WaveFunction::int_potential (const Particles& particles) {
    double a = parameters.a;

    for(int k = 0; k < particles.num_particles; ++k){
        for (int i = 0; i < k; ++i){
         if (particles.compute_distance(k, i) <= a)
            return numeric_limits<double>::infinity();
        }
    }

    return 0;
}


double WaveFunction::log_correlation_part (const Particles& particles) {
    double u = 0;
    for (int i = 0; i < particles.num_particles-1; ++i)
        for (int j = i+1; j < particles.num_particles; ++j)
            u += evaluate_u(particles, i, j);
}


double WaveFunction::evaluate_wavefunction (const Particles& particles) {
    return onebody_part(particles) * exp(log_correlation_part(particles));
}


// Quantum force = drift force
// TODO: Test quantum force
vector<double> WaveFunction::quantum_force (const Particles& particles) {
    vector<double> force(particles.num_dimensions*particles.num_particles);

    // Compute quantum force per particle
    for (int i = 0; i < particles.num_particles; ++i) {
        // Unary terms
        vector<double> single_particle_grad = single_particle_function.evaluate_gradient(
            particles.get_particle(i)
        );
        vector<double> u_gradient(particles.num_dimensions, 0);
        // Pairwise interactions
        for (int j = 0; j < particles.num_particles; ++j) {
            if (j == i) continue;

            // TODO: Wrap this stuff in a function
            double u_diff = deriv_u(particles, i, j);
            vector<double> distance = particles.compute_normalised_distance_vector(i, j);

            for (int k = 0; k < particles.num_dimensions; ++k)
                u_gradient[k] += u_diff*distance[k];
            // End wrap
        }

        // Update force
        for (int j = 0; j < particles.num_dimensions; ++j)
            force[i*particles.num_dimensions+j] = single_particle_grad[j] + u_gradient[j];
    }

    return force;
}


double WaveFunction::second_deriv_wavefunction_quotient (const Particles& particles) {
    double wave_function_quotient;

    // Differentiate wrt all particles
    for (int k = 0; k < particles.num_particles; ++k){
        const Particle current_particle = particles.get_particle(k);
        double single_particle_value = single_particle_function.evaluate(current_particle);

        // Single particle relative laplacian
        double relative_laplacian = single_particle_function.evaluate_laplacian(current_particle)
                                   /single_particle_value;
        wave_function_quotient += relative_laplacian;

        // Single particle relative gradient
        vector<double> rel_gradient = single_particle_function.evaluate_gradient(current_particle);
        for (int l = 0; l < particles.num_dimensions; ++l)
            rel_gradient[l] /= single_particle_value;

        // First summation
        for (int j = 0; j < particles.num_particles; ++j){
            if (j == k) continue;

            vector<double> normalised_distance = particles.compute_normalised_distance_vector(j,k);
            double u_diff = deriv_u(particles, j, k);

            for (int l= 0; l < particles.num_dimensions; ++l)
                wave_function_quotient += rel_gradient[k]*normalised_distance[l]*u_diff;
        }

        // Double summation
        for (int i = 0; i < particles.num_particles; ++i) {
            if (i == k) continue;
            vector<double> rel_distance_i = particles.compute_normalised_distance_vector(k, i);
            double u_diff_i = deriv_u(particles, k, i);

            for (int j = 0; j < particles.num_particles; ++j) {
                if (j == k) continue;
                vector<double> rel_distance_j = particles.compute_normalised_distance_vector(k, j);
                double u_diff_j = deriv_u(particles, k, j);

                for (int l = 0; l < particles.num_dimensions; ++l){
                    double i_term = u_diff_i*rel_distance_i[l];
                    double j_term = u_diff_j*rel_distance_j[l];

                    wave_function_quotient += i_term*j_term;
                }
            }
        }

        // Final Summation
        for (int j = 0; j < particles.num_particles; ++j){
            if (j == k) continue;
            wave_function_quotient += second_deriv_u(particles, j, k);
            wave_function_quotient += 2*deriv_u(particles, j, k)/particles.compute_distance(j, k);
        }

    }

    return wave_function_quotient;
}


double WaveFunction::evaluate_PDF (const Particles& particles) {
    return pow(evaluate_wavefunction(particles), 2);
}

