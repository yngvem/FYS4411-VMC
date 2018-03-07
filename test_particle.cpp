#include "particle.hpp"

using namespace std;


void test_distance_from_origin() {
    Particle particle(3, 1);

    double x = particle.distance_from_origin();
}


int main() {
    test_distance_from_origin();
}