#ifndef TEST_WAVEFUNCTION_HPP
#define TEST_WAVEFUNCTION_HPP

#include <iostream>
#include <assert.h>
#include "wavefunction.hpp"

using namespace std;

class TestWavefunction
{
    vector <vector <double>> create_positions_vector(int num_particles, int num_dimensions,
                                                     double positions_);
public:
    TestWavefunction();
    void test_local_energy();
    void test_second_deriv_quotient();
    void test_pdf();
    void test_quantum_force();
    void test_wavefunction();
};

class TestSingleParticleFunction{

    vector <vector <double>> create_positions_vector(int num_particles, int num_dimensions,
                                                     double positions_);
public:
    TestSingleParticleFunction();
    void test_evaluate();
    void test_gradient();
    void test_laplace();
};

#endif // TEST_WAVEFUNCTION_HPP
