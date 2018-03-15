#include <iostream>
#include "particles.hpp"
#include "test_wavefunction.hpp"

using namespace std;

vector <vector <double>> create_positions_vector(
        int num_particles,
        int num_dimensions,
        double positions_)
{
    vector <vector <double>> positions(num_particles, vector<double>(num_dimensions, positions_));
    return positions;
}

int main()
{
    WaveFunctionParameters params;


    cout << "Testing laplacian of single particle wf evaluation:" << endl;
    cout << "Alpha equal to 1/2 and beta equal to 1" << endl;
    params.alpha = 0.5;
    params.beta = 1;

    std::cout << "Alpha: " << params.alpha << endl
              << "Beta: " << params.beta << endl;

    cout << "Particle placed at the origin" << endl;
    vector <vector <double>> position = create_positions_vector(1, 1, 0);

    cout << "Mass equal to 1" << endl;
    double mass = 1;

    // Start test
    Particles particles(position, mass);
    SingleParticleFunction single_particle_function(params);

    const Particle current_particle = particles.get_particle(0);
    double single_particle_laplacian = single_particle_function.evaluate_laplacian(current_particle);

    cout << "laplacian of wf is equal to: " <<single_particle_laplacian << endl;
    assert(abs(single_particle_laplacian+ 1)< 0.00001 );
    cout << "Test passed." << endl;

}
