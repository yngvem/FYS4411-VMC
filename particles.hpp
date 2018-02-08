#ifndef PARTICLES
#define PARTICLES

#include <vector>
#include <cmath>
#include <random>
#include <stdexcept>
#include "vec3/vec3.h"

using namespace std;

class Particles {
    int num_particles;
    int num_dimensions;
    vector<double> R;

    // Initialisation funcitons
    vector<double> initial_R();
    vector<double> initial_r();

    // Arithmetic operations
    double norm(vector<double>& v);
    vector<double>& difference(vector<double>& v_1, vector<double>& v_2);
    vector<double>& multiply_by_scalar(vector<double>& v, double a);
    double inner_product(vector<double>& v_1, vector<double>& v_2);

public:
    Particles (int num_particles_, int num_dimensions_);
    Particles (int num_particles_, int num_dimensions_, vector<double>& R_);

    const vector<double>& get_particle (int particle_number) const;
    vector<double> compute_distance_vector (int first_particle_idx,
                                            int second_particle_idx);
    double compute_distance  (int first_particle_idx , int second_particle_idx);
};


#endif
