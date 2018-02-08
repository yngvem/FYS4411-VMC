#include <vector>
#include <cmath>
#include <stdexcept>
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

void assert_same_size (vector<double>& v_1, vector<double>& v_2) {
    if (v_1.size() != v_2.size())
        throw invalid_argument(
            "Can't take inner product of vectors of different size"
        );
}

double Particles::norm(vector<double> &v) {
    return sqrt(inner_product(v, v));
}

vector<double> Particles::difference(vector<double> &v_1, vector<double> &v_2) {
    assert_same_size(v_1, v_2);
    vector<double> difference_vector(v_1.size());
    for (int i = 0; i < v_1.size(); ++i)
        difference_vector[i] = v_1-v_2;

    return difference_vector;
}

vector<double> Particles::multiply_by_scalar(vector<double>& v, double a) {
    vector<double> result(v.size());
    for(int i = 0; i < v.size() ; ++i)
        result[i] = a*v[i];
    return result;
}

double Particles::inner_product (vector<double>& v_1, vector<double>& v_2) {
    assert_same_size(v_1, v_2);
    double product = 0;
    for (int i = 0; i < v_1.size(); ++i)
        product += v_1[i]*v_2[i];

    return product;
}

vector<double> Particles::get_particle(int particle_number) {
    int start_idx = particle_number*num_dimensions;
    vector<double> single_particle(num_dimensions);
    for(int i = 0; i < num_dimensions; ++i)
        single_particle[i] = R[start_idx+i];

    return single_particle;
}

vector<double> Particles::compute_distance_vector (int first_particle_idx,
                                                   int second_particle_idx) {
    vector<double> first_particle = get_particle(first_particle_idx);
    vector<double> second_particle = get_particle(second_particle_idx);
    return difference(first_particle, second_particle);
}

double compute_distance(int first_particle_idx , int second_particle_idx) {
    return norm(compute_distance_vector(first_particle_idx, second_particle_idx));
}

