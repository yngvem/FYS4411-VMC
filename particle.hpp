#ifndef PARTICLE
#define PARTICLE

#include <vector>

using namespace std;

class Particle {
    vector<double> r;
    vector<double> r_old;

    void initiate_r();

    // Linear algebra functions.
    double norm(const vector<double>& vec) const;
    double weighted_norm(const vector<double>& vec,
                         const vector<double>& weights) const;

public:
    const int num_dimensions;
    const double m;

    Particle(int num_dimensions_, double mass);
    Particle(vector<double> r, double mass);

    // Obtaining information about the position.
    const vector<double>& get_position() const;
    double weighted_distance_from_origin(vector<double> weights) const;
    double distance_from_origin() const;
    vector<double> distance_vector(const Particle& particle) const;
    double distance_between(const Particle& particle) const;
    double operator[](int i) const; // TODO: Should this be const and reference?
    void perturb(vector<double> pertubation);
    void reject_perturbation();
};

#endif
