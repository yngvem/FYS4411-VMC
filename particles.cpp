#include <vector>
#include <cmath>
#include <stdexcept>
#include "particles.hpp"

using namespace std;

Particles::Particles (ParticlesParams params):
    num_particles(params.num_particles),
    num_dimensions(params.num_dimensions)
{
for (int i = 0; i < num_particles; ++i)
    particles.push_back(Particle(num_dimensions, params.mass));
}

Particles::Particles (int num_particles_, int num_dimensions_, double mass) :
    num_particles(num_particles_),
    num_dimensions(num_dimensions_)
{
    for (int i = 0; i < num_particles; ++i) 
        particles.push_back(Particle(num_dimensions, mass));
}


Particles::Particles (int num_particles_, int num_dimensions_, 
                      vector<double> mass) :
    num_particles(num_particles_),
    num_dimensions(num_dimensions_)
{   
    for (int i = 0; i < num_particles; ++i)
        particles.push_back(Particle(num_dimensions, mass[i]));
}

Particles::Particles (vector<vector <double>> positions, double mass) :
    num_particles(positions.size()),
    num_dimensions(positions[0].size())
{
    for (int i = 0; i < num_particles; ++i)
        particles.push_back(Particle(positions[i], mass));
}

Particles::Particles (vector<vector <double>> positions, vector<double> mass):
    num_particles(positions.size()),
    num_dimensions(positions[0].size())
{
    for (int i = 0; i < num_particles; ++i)
        particles.push_back(Particle(positions[i], mass[i]));
}

const Particle& Particles::get_particle(int particle_idx) const {
    return particles[particle_idx];
}


vector<double> Particles::compute_normalised_distance_vector (int first_particle_idx,
                                                              int second_particle_idx) const {
    vector<double> distance_vector = compute_distance_vector(first_particle_idx,
                                                             second_particle_idx);
    double distance = compute_distance(first_particle_idx, second_particle_idx);
    for (int i = 0; i < num_dimensions; ++i) 
        distance_vector[i] /= distance;
    
    return distance_vector;
}


vector<double> Particles::compute_distance_vector (int first_particle_idx,
                                                   int second_particle_idx) const {
    return particles[first_particle_idx].distance_vector(
        particles[second_particle_idx]
    );
}


double Particles::compute_distance(int first_particle_idx , int second_particle_idx) const {
    return particles[first_particle_idx].distance_between(
        particles[second_particle_idx]
    );
}

void Particles::perturb_particle(int particle_idx, vector<double> perturbation) {
    particles[particle_idx].perturb(perturbation);
}

void Particles::reject_perturbation(int particle_idx) {
    particles[particle_idx].reject_perturbation();
}
