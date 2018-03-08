#ifndef PARTICLES
#define PARTICLES

#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>
#include "particle.hpp"

using namespace std;

struct ParticlesParams {
    int num_particles = 1;
    int num_dimensions = 3;
    double mass = 1;
};

class Particles {
    vector<Particle> particles;

    // Initialisation funcitons
    vector<Particle> init_particles(double mass);
    vector<Particle> init_particles(vector<double> mass);

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
    Particles (vector<vector <double>> positions, double mass);
    Particles (vector<vector <double>> positions, vector<double> mass);
    Particles (ParticlesParams params);
    
    const Particle& get_particle (int particle_number) const;
    vector<double> compute_distance_vector (int first_particle_idx,
                                            int second_particle_idx) const;
    vector<double> compute_normalised_distance_vector (int first_particle_idx,
                                                       int second_particle_idx) const;
    double compute_distance  (int first_particle_idx , int second_particle_idx) const;
    double compute_R_squared () const;

    void perturb_particle(int particle_idx, vector<double> perturbation);
    void reject_perturbation(int particle_idx);
};


#endif
