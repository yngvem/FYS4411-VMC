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
    INFO("beta = 0.5");
    WaveFunctionParameters params;
    params.alpha = 0.5;
    params.beta = 0.5;

    // Set global particle parameters
    INFO("Particles have mass = 1");
    double mass = 1;

    // Create single particle function
    SingleParticleFunction single_particle_function(params);

    SECTION("ONE DIMENSION") {
        int num_dimensions = 1;
        SECTION("PLACED AT THE ORIGIN") {
            vector<double> position(num_dimensions, 0);
            Particle particle(position, mass);

            SECTION("TESTING EVALUATION") {
                REQUIRE(single_particle_function.evaluate(particle) == Approx(1));
            }

            SECTION("TESTING GRADIENT") {
                vector<double> grad = single_particle_function.evaluate_gradient(particle);
                for (auto x : grad)
                    REQUIRE(x == Approx(0));
            }

            SECTION("TESTING LAPLACIAN") {
                REQUIRE(single_particle_function.evaluate_laplacian(particle) == Approx(-1));
            }
        }

        SECTION("PLACED AT x=1") {
            vector<double> position(num_dimensions, 1);
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

    SECTION("TWO DIMENSIONS") {
        int num_dimensions = 2;
        SECTION("PLACED AT THE ORIGIN") {
            vector<double> position (num_dimensions, 0);
            Particle particle(position, mass);

            SECTION("TESTING EVALUATION") {
                REQUIRE(single_particle_function.evaluate(particle) == Approx(1));
            }

            SECTION("TESTING GRADIENT") {
                vector<double> grad = single_particle_function.evaluate_gradient(particle);
                for (auto x : grad)
                    REQUIRE(x == Approx(0));
            }

            SECTION("TESTING LAPLACIAN") {
                REQUIRE(single_particle_function.evaluate_laplacian(particle) == Approx(-2));
            }
        }

        SECTION("PLACED AT x=1") {
            vector<double> position(num_dimensions, 1);
            Particle particle(position, mass);

            SECTION("TESTING EVALUATION") {
                REQUIRE(single_particle_function.evaluate(particle) == Approx(exp(-1)));
            }

            SECTION("TESTING GRADIENT") {
                vector<double> grad = single_particle_function.evaluate_gradient(particle);
                for (auto x : grad)
                   REQUIRE(x == Approx(-exp(-1)));
            }

            SECTION("TESTING LAPLACIAN") {
                REQUIRE(single_particle_function.evaluate_laplacian(particle) == Approx(0));
            }
        }
    }

    SECTION("Three DIMENSIONS") {
        int num_dimensions = 3;
        SECTION("PLACED AT THE ORIGIN") {
            vector<double> position (num_dimensions, 0);
            Particle particle(position, mass);

            SECTION("TESTING EVALUATION") {
                REQUIRE(single_particle_function.evaluate(particle) == Approx(1));
            }

            SECTION("TESTING GRADIENT") {
                vector<double> grad = single_particle_function.evaluate_gradient(particle);
                for (auto x : grad)
                    REQUIRE(x == Approx(0));
            }

            SECTION("TESTING LAPLACIAN") {
                REQUIRE(single_particle_function.evaluate_laplacian(particle) == Approx(-2.5));
            }
        }

        SECTION("PLACED AT x=1") {
            vector<double> position(num_dimensions, 1);
            Particle particle(position, mass);

            SECTION("TESTING EVALUATION") {
                REQUIRE(single_particle_function.evaluate(particle) == Approx(exp(-1.25)));
            }

            SECTION("TESTING GRADIENT") {
                vector<double> grad = single_particle_function.evaluate_gradient(particle);
                   REQUIRE(grad[0] == Approx(-exp(-1.25)));
                   REQUIRE(grad[1] == Approx(-exp(-1.25)));
                   REQUIRE(grad[2] == Approx(-0.5*exp(-1.25)));
            }

            SECTION("TESTING LAPLACIAN") {
                REQUIRE(single_particle_function.evaluate_laplacian(particle) == Approx(-0.25*exp(-1.25)));
            }
        }
    }

}



TEST_CASE("WAVE FUNCTION") {
    // Set wave function parameters
    INFO("alpha = 0.5");
    INFO("beta = 1");
    INFO("hbar = 1");
    INFO("omega_ho = 1");
    INFO("omega_z = 0.5");
    INFO("a = 1");
    WaveFunctionParameters params;
    params.alpha = 0.5;
    params.beta = 1;
    params.hbar = 1;
    params.omega_ho = 1;
    params.omega_z = 0.5;
    params.a = 1;

    // Set global particle parameters
    INFO("Particles have mass = 1");
    double mass = 1;

    // Create single particle function
    WaveFunction wave_function(params);

    SECTION("SINGLE PARTICLE") {
        int num_particles = 1;
        SECTION("ONE DIMENSION") {
            int num_dimensions = 1;
            SECTION("PLACED AT THE ORIGIN") {
                vector<vector<double>> position(num_particles, vector<double>(num_dimensions, 0));
                Particles particles(position, mass);

                SECTION("TESTING PDF") {
                    REQUIRE(wave_function.evaluate_PDF(particles) == Approx(1));
                }

                SECTION("TESTING LOCAL ENERGY") {
                    REQUIRE(wave_function.local_energy(particles) == Approx(0.5));                    
                }

                SECTION("TESTING QUANTUM FORCE") {
                    vector<double> force = wave_function.quantum_force(particles);
                    for (auto x : force) 
                        REQUIRE(x == Approx(0));
                }
            }

            SECTION("PLACED AT x=1") {
                vector<vector<double>> position = get_pos_vectors(num_particles, num_dimensions, 1);
                Particles particles(position, mass);

                SECTION("TESTING PDF") {
                    REQUIRE(wave_function.evaluate_PDF(particles) == Approx(pow(exp(-0.5), 2)));
                }

                SECTION("TESTING LOCAL ENERGY") {
                    REQUIRE(wave_function.local_energy(particles) == Approx(0.5));
                }

                SECTION("TESTING QUANTUM FORCE") {
                    vector<double> force = wave_function.quantum_force(particles);
                    for (auto x: force)
                        REQUIRE(x == Approx(-2));
                    INFO("Dpsi = -psi");
                    INFO("2Dpsi/psi = -2");
                }
            }
        }
    }
}


// void TestWavefunction::test_second_deriv_quotient () {
//     cout << "Particle placed at the origin" << endl;
//     vector <vector <double>> position = create_positions_vector(num_particles, num_dimensions, 10);
// 
//     cout << "Mass equal to 1" << endl;
//     double mass = 1;
// 
//     // Start test
//     WaveFunction wave_function(params);
//     Particles particles(position, mass);
// 
//     cout << "Wave function quotient: " <<wave_function.local_energy(particles) << endl;
//     assert(abs(wave_function.second_deriv_wavefunction_quotient(particles) - 0.5) < 0.00001);
//     cout << "Test passed." << endl;
// }
// 
// void TestWavefunction::test_local_energy () {
//     WaveFunctionParameters params;
// 
//     cout << "Testing local energy:" << endl;
//     cout << "One particle in one dimension:" << endl;
//     int num_dimensions = 1;
//     int num_particles = 1;
// 
//     cout << "Omega equal to 1" << endl;
//     double omega_ho = 1;
//     params.omega_ho = omega_ho;
// 
//     cout << "hbar equal to 1" << endl;
//     double hbar = 1;
//     params.hbar = hbar;
// 
//     cout << "Alpha equal to 1/2 and beta equal to 1" << endl;
//     params.alpha = 1/2;
//     params.beta = 1;
// 
// 
//     cout << "Particle placed at the origin" << endl;
//     vector <vector <double>> position = create_positions_vector(num_particles, num_dimensions, 1);
// 
//     cout << "Mass equal to 1" << endl;
//     double mass = 1;
// 
//     // Start test
//     WaveFunction wave_function(params);
//     Particles particles(position, mass);
// 
//     cout << "Local energy is equal to: " <<wave_function.local_energy(particles) << endl;
//     assert(wave_function.local_energy(particles) == num_dimensions*hbar*omega_ho/2);
//     cout << "Test passed." << endl;
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
//}

// Creates a vector of vectors with all values equal to `positions`.
vector <vector <double>> TestWavefunction::create_positions_vector(
        int num_particles,
        int num_dimensions,
        double positions)
{
    return vector<vector <double>> (num_particles, vector<double>(num_dimensions, positions));
}
