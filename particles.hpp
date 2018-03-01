#ifndef PARTICLES
#define PARTICLES

#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>
#include "particle.hpp"

using namespace std;

// TODO: Create single particle class

class Particles {
    vector<double> R;
    
    vector<Particle> particles;

    // Initialisation funcitons
    vector<Particle> init_particles(double mass);
    vector<Particle> init_particles(vector<double> mass);
    // vector<double> initial_R();
    // vector<double> initial_r();
    // vector<double> initial_m(double mass_);

    // Arithmetic operations
    double norm(vector<double>& v) const;
    vector<double> difference(vector<double>& v_1, vector<double>& v_2) const;
    vector<double> multiply_by_scalar(vector<double>& v, double a) const;
    double inner_product(vector<double>& v_1, vector<double>& v_2) const;

public:
    const int num_particles;
    const int num_dimensions;
    // const vector<double> mass;

    Particles (int num_particles_, int num_dimensions_, vector<double> mass);
    Particles (int num_particles_, int num_dimensions_, double mass);
    
    const Particle& get_particle (int particle_number) const;
    vector<double> compute_distance_vector (int first_particle_idx,
                                            int second_particle_idx) const;
    vector<double> compute_normalised_distance_vector (int first_particle_idx,
                                                       int second_particle_idx) const;
    double compute_distance  (int first_particle_idx , int second_particle_idx) const;
    double compute_R_squared () const;
};


#endif
