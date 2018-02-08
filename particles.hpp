#ifndef PARTICLES
#define PARTICLES

#include <vector>
#include <random>

using namespace std;

class Particles {
    int num_particles;
    int num_dimensions;
    vector<double> R;

    vector<double> initial_R();
    vector<double> initial_r();

public:
    vector<double> get_particle (int particle_number);
    Particles (int num_particles_, int num_dimensions_);
    Particles (int num_particles_, int num_dimensions_, vector<double>& R_);
};


#endif
