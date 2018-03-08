#include <vector>
#include <cmath>
#include <iostream>
#include "particle.hpp"

using namespace std;


Particle::Particle (int num_dimensions_, double mass) :
    num_dimensions(num_dimensions_),
    m(mass)
{
    initiate_r();
}

Particle::Particle (const vector<double> r_, double mass) :
    num_dimensions(r_.size()),
    m(mass),
    r(r_),
    r_old(r_)
{
}

void Particle::initiate_r () {
    for (int i = 0; i < num_dimensions; ++i) {
        r.push_back(0);
        r_old.push_back(r[i]);
    }
}

double Particle::weighted_norm (const vector<double>& vec,
                                const vector<double>& weights) const {
    // TODO: Assert equal size
    double sum = 0;
    for (int i = 0; i < vec.size(); ++i) {
        sum += vec[i]*vec[i]*weights[i];
    }
    return sqrt(sum);
}

double Particle::norm (const vector<double>& vec) const {
    vector<double> weights(vec.size(), 1);  // TODO: Precompute and store as private.
    return weighted_norm(vec, weights);
}

const vector<double>& Particle::get_position () const {
    return r;
}

double Particle::distance_from_origin () const {
    return norm(r);
}

double Particle::weighted_distance_from_origin (vector<double> weights) const {
    return weighted_norm(r, weights);
}

vector<double> Particle::distance_vector (const Particle& particle) const {
    // TODO: Assert equal size
    vector<double> distance_vector_(num_dimensions, 0);
    for (int i = 0; i < num_dimensions; ++i) {
        distance_vector_[i] = r[i] - particle.get_position()[i];
    }

    return distance_vector_;
}

double Particle::distance_between (const Particle& particle) const {
    return norm(distance_vector(particle));
}

double Particle::operator[] (int i) const {
    return r[i];
}


void Particle::perturb(vector<double> pertubation){
    for (int i = 0; i < num_dimensions; ++i){
        r_old[i] = r[i];
        r[i] += pertubation[i];
    }
}
void Particle::reject_perturbation(){
    for (int i = 0; i < num_dimensions; ++i)
        r[i] = r_old[i];
}
