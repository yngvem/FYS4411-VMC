#include <iostream>
#include <assert.h>
#include "test_wavefunction.hpp"

using namespace std;

TestWavefunction::TestWavefunction () :
    position(vector<vector <double>>(1, vector<double>(1, 1))),
    particles(Particles(position, 1))
{
    test_local_energy();
}

void TestWavefunction::test_local_energy () {
    cout << particles.get_particle(0).get_position()[0] << endl;
    cout << wave_function.ext_potential(particles) << endl;
    cout << wave_function.int_potential(particles) << endl;
    cout << wave_function.second_deriv_wavefunction_quotient(particles) << endl;
    cout << wave_function.local_energy(particles) << endl;
}
