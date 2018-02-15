#ifndef PARTICLE
#define PARTICLE

#include <vector>

using namespace std;

class Particle {
    vector<double> r;
    

    vector<double> initiate_r();

    // Linear algebra functions.
    double norm(const vector<double>& vec);
    double weighted_norm(const vector<double>& vec,
                         const vector<double>& weights);

public:
    const int num_dimensions;
    const double m;

    Particle(int num_dimensions_, double mass);
    Particle(double mass, vector<double> r);

    // Obtaining information about the position.
    const vector<double>& get_position() const;
    double weighted_distance_from_origin(vector<double> weights);
    double distance_from_origin();
    const vector<double> distance_vector(const Particle& particle);
    const double distance_between(const Particle& particle);
    double operator[](int i); // TODO: Should this be const and reference?
};

#endif