#ifndef PARTICLES
#define PARTICLES

#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>

using namespace std;

// TODO: Create single particle class

class Particles {
    int num_particles;
    int num_dimensions;
    vector<double> mass;
    vector<double> R;

    // Initialisation funcitons
    vector<double> initial_R();
    vector<double> initial_r();
    vector<double> initial_m(double mass_);

    // Arithmetic operations
    double norm(vector<double>& v);
    vector<double> difference(vector<double>& v_1, vector<double>& v_2);
    vector<double> multiply_by_scalar(vector<double>& v, double a);
    double inner_product(vector<double>& v_1, vector<double>& v_2);

public:
    Particles (int num_particles_, int num_dimensions_, vector<double> mass_);
    Particles (int num_particles_, int num_dimensions_, vector<double> mass_,
               vector<double>& R_);
    Particles (int num_particles_, int num_dimensions_, double mass_);
    Particles (int num_particles_, int num_dimensions_, double mass_,
               vector<double>& R_);
    

    vector<double> get_particle (int particle_number);
    vector<double> compute_distance_vector (int first_particle_idx,
                                            int second_particle_idx);
    double compute_distance  (int first_particle_idx , int second_particle_idx);
    double compute_R_squared ();
};


#endif
