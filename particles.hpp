#ifndef PARTICLES
#define PARTICLES

#include <vector>
#include <cmath>
#include <random>

using namespace std;

class Particles {
    int num_particles;
    int num_dimensions;
    vector<double> R;

    vector<double> initial_R();
    vector<double> initial_r();

    double compute_norm(vector<double>& v);

public:
    Particles (int num_particles_, int num_dimensions_);
    Particles (int num_particles_, int num_dimensions_, vector<double>& R_);

    vector<double> get_particle (int particle_number);
    vector<double> compute_distance_vector (int first_particle , int second_particle);
    double compute_distance  (int first_particle , int second_particle);
};


#endif
