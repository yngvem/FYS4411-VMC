#ifndef TEST_WAVEFUNCTION_HPP
#define TEST_WAVEFUNCTION_HPP

#include <iostream>
#include <assert.h>
#include "wavefunction.hpp"

using namespace std;

class TestWavefunction
{
public:
    vector <vector <double>> position;
    Particles particles;
    WaveFunction wave_function;

    TestWavefunction();
    void test_local_energy();
    void test_pdf();
    void test_quantum_force();
    void test_wavefunction();
};

#endif // TEST_WAVEFUNCTION_HPP
