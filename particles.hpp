#ifndef PARTICLES
#define PARTICLES

#include <vector>
#include <random>

using namespace std;

class Particles {
    int n_particles;
    int n_dimensions;
    vector<double> R;

    vector<double> initial_R(int num_dimensions, int num_particles);
    vector<double> initial_r(int num_dimensions);

public:
    vector<double> get_particle (int particle_number);
    Particles (int num_particles, int num_dimensions);
    Particles (int num_particles, int num_dimensions, int R);
};


#endif