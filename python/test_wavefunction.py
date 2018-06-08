"""

"""


__author__ = 'Yngve Mardal Moe'
__email__ = 'yngve.m.moe@gmail.com'


import numpy as np
from pytest import fixture
from particles import Particles
from wave_function import SingleParticleFunction, WaveFunction


EPS = 1e-5
np.random.seed(1)

def move_to_origin(particles):
    particles._positions *= 0

def move_to_1(particles):
    move_to_origin(particles)
    particles._positions += 1


@fixture
def single_point_particle_1d():
    return Particles(num_particles=1, num_dimensions=1, mass=1, diameter=0)


@fixture
def single_point_particle_1d_origin(single_point_particle_1d):
    move_to_origin(single_point_particle_1d)
    return single_point_particle_1d


@fixture
def single_point_particle_1d_at_1(single_point_particle_1d):
    move_to_1(single_point_particle_1d)
    return single_point_particle_1d


class TestSingleParticleFunction:
    @fixture
    def single_particle_function_1d_at_0(self, single_point_particle_1d_origin):
        return SingleParticleFunction(particles=single_point_particle_1d_origin,
                                      alpha=1)

    @fixture
    def single_particle_function_1d_at_1(self, single_point_particle_1d_at_1):
        return SingleParticleFunction(particles=single_point_particle_1d_at_1,
                                      alpha=1)

    def test_evaluate_1d_origin(self, single_particle_function_1d_at_0):
        assert np.abs(single_particle_function_1d_at_0()-1).sum() < EPS

    def test_gradient_1d_origin(self, single_particle_function_1d_at_0):
        gradient = single_particle_function_1d_at_0.gradient
        assert np.abs(gradient).sum() < EPS

    def test_laplacian_1d_origin(self, single_particle_function_1d_at_0):
        laplacian = single_particle_function_1d_at_0.laplacian
        assert np.abs(laplacian + 2).sum() < EPS

    def test_evaluate_1d_at_1(self, single_particle_function_1d_at_1):
        assert np.abs(single_particle_function_1d_at_1() - 1/np.e).sum() < EPS

    def test_gradient_1d_at_1(self, single_particle_function_1d_at_1):
        gradient = single_particle_function_1d_at_1.gradient
        assert np.abs(gradient + 2/np.e).sum() < EPS

    def test_laplacian_1d_at_1(self, single_particle_function_1d_at_1):
        laplacian = single_particle_function_1d_at_1.laplacian
        assert np.abs(laplacian - 2/np.e).sum() < EPS
        

class TestWaveFunction:
    @fixture
    def two_particles_1d_random(self):
        return Particles(num_particles=2, num_dimensions=1, mass=1, diameter=0.1)

    @fixture
    def wf_two_particles_1d_random(self, two_particles_1d_random):
        return WaveFunction(two_particles_1d_random, alpha=1)

    @fixture
    def five_particles_3d_random(self):
        return Particles(num_particles=5, num_dimensions=3, mass=1, diameter=0.1)

    @fixture
    def wf_5_particles_3d_random(self, five_particles_3d_random):
        return WaveFunction(five_particles_3d_random)

    @fixture
    def wf_one_particle_1d_origin(self, single_point_particle_1d_origin):
        return WaveFunction(single_point_particle_1d_origin, alpha=1)

    def test_local_energy_one_particle_1d(self, wf_one_particle_1d_origin):
        assert np.abs(wf_one_particle_1d_origin.local_energy - 1) < EPS

    def test_local_energy_two_particles_1d(self, wf_two_particles_1d_random):
        spf = wf_two_particles_1d_random.single_particle_function
        parts = wf_two_particles_1d_random.particles
        difference = parts.positions[0] - parts.positions[1]
        u_deriv = wf_two_particles_1d_random.u_derivatives(np.abs(parts.differences))[0, 1]
        u_second_deriv = wf_two_particles_1d_random.u_second_deriv(np.abs(parts.differences))[0, 1]

        single_particle_part = spf.relative_laplacian.sum()
        gradient_part = np.sign(difference)*u_deriv
        gradient_part *= spf.relative_gradient[0] - spf.relative_gradient[1]
        
        pairwise_part1 = 2*u_deriv**2
        pairwise_part = u_deriv * (2/np.abs(difference))
        pairwise_part += u_second_deriv
        pairwise_part *= 2

        print('True values:')
        print(f'Single particle parts {single_particle_part}')
        print(f'Gradient part {single_particle_part + gradient_part}')
        print(f'First pairwise part {single_particle_part + gradient_part + pairwise_part1}')

        ke = single_particle_part + gradient_part + pairwise_part + pairwise_part1
        print(f'All parts {ke}')
        ke *= -wf_two_particles_1d_random.hbar/(2*parts.mass)
        print(f'All parts after constant factor {ke}')

        print('\nComputed values:')
        assert wf_two_particles_1d_random.kinetic_energy - ke < 1e-5


if __name__ == '__main__':
    pass

