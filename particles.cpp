#include <vector>
#include <cmath>
#include <stdexcept>
#include "particles.hpp"

using namespace std;

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
