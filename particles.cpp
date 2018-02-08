#include <vector>
#include <cmath>
#include "particles.hpp"

using namespace std;

Particles::Particles (int num_particles_, int num_dimensions_) :
    num_particles(num_particles_),
    num_dimensions(num_dimensions_)
{
    R = initial_R();
}

Particles::Particles (int num_particles_, int num_dimensions_, 
                      vector<double>& R_) :
    num_particles(num_particles_),
    num_dimensions(num_dimensions_),
    R(R_)
{}

vector<double> Particles::initial_R() {
    int position_length = num_particles*num_dimensions;
    vector<double> R0(position_length);

    for (int i = 0; i < position_length; ++i){
        vector<double> single_particle_r0 = initial_r();

        for (int j = 0; j < num_dimensions; ++j)
            R0[i*num_dimensions + j] = single_particle_r0[j];
    }
}

vector<double> Particles::initial_r() {
    vector<double> r(num_dimensions, 0);
    return r;
}

double Particles::compute_norm(vector<double> &v) {
    double v_sq = 0;
    for(auto* r : v)
        v_sq += v*v;

    return sqrt(v_sq);
}

vector<double> Particles::compute_distance_vector(int first_particle, int second_particle) {

}
