#include <vector>
#include <cmath>
#include "particle.hpp"

using namespace std;


Particle::Particle (int num_dimensions_, double mass) :
    num_dimensions(num_dimensions_),
    m(mass)
{
    r = initiate_r();
}

Particle::Particle (double mass, const vector<double> r_) :
    num_dimensions(r_.size()),
    m(mass),
    r(r_)
{}

vector<double> Particle::initiate_r () {
    vector <double> r_(num_dimensions, 0);
    return r_;
}

double Particle::weighted_norm (const vector<double>& vec,
                                const vector<double>& weights) {
    // TODO: Assert equal size
    double sum = 0;
    for (int i = 0; i < vec.size(); ++i) {
        sum += vec[i]*vec[i]*weights[i];
    }
    return sqrt(sum);
}

double Particle::norm (const vector<double>& vec) {
    vector<double> weights(vec.size(), 1);  // TODO: Precompute and store as private.
    return weighted_norm(vec, weights);
}

const vector<double>& Particle::get_position () const {
    return r;
}

const double Particle::distance_from_origin () {
    return norm(r);
}

const vector<double> Particle::distance_vector (const Particle& particle) {
    // TODO: Assert equal size
    vector<double> distance_vector_(0, num_dimensions);
    for (int i = 0; i < num_dimensions; ++i) {
        distance_vector_[i] = r[i] - particle.get_position()[i];
    }

    return distance_vector_;
}

const double Particle::distance_between (const Particle& particle) {
    return norm(distance_vector(particle));
}

double Particle::operator[] (int i) {
    return r[i];
}