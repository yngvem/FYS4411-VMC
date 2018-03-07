#define CATCH_CONFIG_MAIN
#include <iostream>
#include <assert.h>
#include "test_wavefunction.hpp"
#include "wavefunction.hpp"
#include "catch.hpp"

using namespace std;

vector <vector <double>> get_pos_vectors(int num_particles, int num_dimensions, double position_) {
    return vector<vector <double>>(num_particles, vector<double>(num_dimensions, position_));
}

TEST_CASE("SINGLE PARTICLE FUNCTION") {
    // Set wave function parameters
    INFO("alpha = 0.5");
    INFO("beta = 1");
    WaveFunctionParameters params;
    params.alpha = 0.5;
    params.beta = 1;

    // Set global particle parameters
    INFO("Particles have mass = 1");
    double mass = 1;

    // Create single particle function
    SingleParticleFunction single_particle_function(params);

    SECTION("SINGLE PARTICLE") {
        int num_particles = 1;
        SECTION("ONE DIMENSION") {
            int num_dimensions = 1;
            SECTION("PLACED AT THE ORIGIN") {
                vector<double> position(1, 0);
                Particle particle(position, mass);

                SECTION("TESTING EVALUATION") {
                    REQUIRE(single_particle_function.evaluate(particle) == Approx(1));
                }

                SECTION("TESTING GRADIENT") {
                    vector<double> grad = single_particle_function.evaluate_gradient(particle);
                    for (auto x : grad)
                       REQUIRE(x == Approx(0));
                }
/*
                SECTION("TESTING LAPLACIAN") {
                    REQUIRE(single_particle_function.evaluate_laplacian(particle) == Approx(-1));
                }*/
            }

            SECTION("PLACED AT x=1") {
                vector<double> position(1, 1);
                Particle particle(position, mass);

                SECTION("TESTING EVALUATION") {
                    REQUIRE(single_particle_function.evaluate(particle) == Approx(exp(-0.5)));
                }

                SECTION("TESTING GRADIENT") {
                    vector<double> grad = single_particle_function.evaluate_gradient(particle);
                    for (auto x : grad)
                        REQUIRE(x == Approx(-exp(-0.5)));
                }

                SECTION("TESTING LAPLACIAN") {
                    REQUIRE(single_particle_function.evaluate_laplacian(particle) == Approx(0));
                }
            }
        }
    }
}

































TestSingleParticleFunction::TestSingleParticleFunction(){
    test_evaluate();
    test_laplace();
    test_gradient();
}

vector <vector <double>> TestSingleParticleFunction::create_positions_vector(
        int num_particles,
        int num_dimensions,
        double positions_)
{
    vector <vector <double>> positions(num_particles, vector<double>(num_dimensions, positions_));
    return positions;
}

void TestSingleParticleFunction::test_evaluate () {
     WaveFunctionParameters params;

     cout << "Testing single particle wf evaluation:" << endl;
     cout << "Alpha equal to 1/2 and beta equal to 1" << endl;
     params.alpha = 1/2;
     params.beta = 1;


     cout << "Particle placed at the origin" << endl;
     vector <vector <double>> position = create_positions_vector(1, 1, 0);

     cout << "Mass equal to 1" << endl;
     double mass = 1;

     // Start test
     Particles particles(position, mass);
     SingleParticleFunction single_particle_function(params);

     const Particle current_particle = particles.get_particle(0);
     double single_particle_value = single_particle_function.evaluate(current_particle);

     cout << "evaluation of wf is equal to: " <<single_particle_value << endl;
     assert(abs(single_particle_value- 1)< 0.00001 );
     cout << "Test passed." << endl;
}

void TestSingleParticleFunction::test_laplace() {
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

void TestSingleParticleFunction::test_gradient() {

}


// Test the wave function.
TestWavefunction::TestWavefunction ()
{
    test_second_deriv_quotient();
    test_local_energy();
}

// Creates a vector of vectors with all values equal to `positions`.
vector <vector <double>> TestWavefunction::create_positions_vector(
        int num_particles,
        int num_dimensions,
        double positions_)
{
    vector <vector <double>> positions(num_particles, vector<double>(num_dimensions, positions_));
    return positions;
}

void TestWavefunction::test_second_deriv_quotient () {
    WaveFunctionParameters params;

    cout << "Testing second derivative quotient:" << endl;
    cout << "One particle in one dimension:" << endl;
    int num_dimensions = 1;
    int num_particles = 1;

    cout << "hbar equal to 1" << endl;
    double hbar = 1;
    params.hbar = hbar;

    cout << "Alpha equal to 1/2 and beta equal to 1" << endl;
    params.alpha = 1/2;
    params.beta = 1;


    cout << "Particle placed at the origin" << endl;
    vector <vector <double>> position = create_positions_vector(num_particles, num_dimensions, 10);

    cout << "Mass equal to 1" << endl;
    double mass = 1;

    // Start test
    WaveFunction wave_function(params);
    Particles particles(position, mass);

    cout << "Wave function quotient: " <<wave_function.local_energy(particles) << endl;
    assert(abs(wave_function.second_deriv_wavefunction_quotient(particles) - 0.5) < 0.00001);
    cout << "Test passed." << endl;
}

void TestWavefunction::test_local_energy () {
    WaveFunctionParameters params;

    cout << "Testing local energy:" << endl;
    cout << "One particle in one dimension:" << endl;
    int num_dimensions = 1;
    int num_particles = 1;

    cout << "Omega equal to 1" << endl;
    double omega_ho = 1;
    params.omega_ho = omega_ho;

    cout << "hbar equal to 1" << endl;
    double hbar = 1;
    params.hbar = hbar;

    cout << "Alpha equal to 1/2 and beta equal to 1" << endl;
    params.alpha = 1/2;
    params.beta = 1;


    cout << "Particle placed at the origin" << endl;
    vector <vector <double>> position = create_positions_vector(num_particles, num_dimensions, 1);

    cout << "Mass equal to 1" << endl;
    double mass = 1;

    // Start test
    WaveFunction wave_function(params);
    Particles particles(position, mass);

    cout << "Local energy is equal to: " <<wave_function.local_energy(particles) << endl;
    assert(wave_function.local_energy(particles) == num_dimensions*hbar*omega_ho/2);
    cout << "Test passed." << endl;
/*
    // Second test:
    cout << "One particle in two dimensions:" << endl;
    num_dimensions = 2;

    cout << "Particle placed at the origin" << endl;
    position = create_positions_vector(num_particles, num_dimensions, 0);

    // Start test
    wave_function(params);
    particles(position, mass);

    assert(wave_function.local_energy(particles) == num_dimensions*hbar*omega_ho/2);
    cout << "Test passed." << endl;
    */
}

