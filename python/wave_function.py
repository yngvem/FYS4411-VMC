"""

"""

__author__ = 'Yngve Mardal Moe'
__email__ = 'yngve.m.moe@gmail.com'

import numpy as np
from numba import jit, float64, void, int64
from particles import Particles
import pickle


class SingleParticleFunction:
    def __init__(self, particles, alpha=0.5, beta=1):
        self.particles = particles

        self.alpha = alpha
        self.beta = beta
        self._weights = np.ones([1, self.particles.num_dimensions])
        if particles.num_dimensions == 3:
            self._weights[0, 2] = beta
        self._weights_sq = self._weights**2

    def __call__(self):
        return np.exp(-self.alpha*self.elliptic_distances**2)

    @property
    def gradient(self):
        return self.relative_gradient*self()[..., np.newaxis]

    @property
    def laplacian(self):
        return self.relative_laplacian*self()

    @property
    def relative_gradient(self):
        positions = self.particles.positions
        return -2*self.alpha*self._weights*positions
        
    @property
    def relative_laplacian(self):
        positions_sq = self.particles.positions**2
        weighted_dot = 4*(self.alpha**2)*positions_sq@self._weights_sq.T
        constant_term = -2*self.alpha*self._weights.sum()
        return weighted_dot+constant_term
        
    @property
    def elliptic_distances(self):
        return np.linalg.norm(self._weights*self.particles.positions, axis=-1)
    
    def set_particles(self, particles):
        self.__init__(particles, alpha=self.alpha, beta=self.beta)
        

class WaveFunction:
    def __init__(self, particles, alpha=0.5, beta=1, hbar=1, omega_ho=1,
                 omega_z=1):
        self.single_particle_function = SingleParticleFunction(
            particles=particles,
            alpha=alpha,
            beta=beta
        )

        self.alpha = alpha
        self.beta = beta
        self.hbar = hbar
        self.a = particles.diameter
        self.particles = particles
        self.num_particles = particles.num_particles

        self.omegas = omega_ho + np.zeros(particles.num_dimensions)
        if particles.num_dimensions == 3:
            self.omegas[2] = omega_z

    @property
    def f(self):
        f = np.zeros((self.num_particles, self.num_particles), dtype='float128')
        nonzero = np.logical_not(self.particles.overlapping)
        f[nonzero] = 1 - self.a/np.abs(self.particles.distances[nonzero])
        return f

    @property
    def u(self):
        # TODO: Think about how to deal with rij<a
        u = np.log(self.f)
        return u

    @property
    def u_derivatives(self):
        return self.a/self._u_deriv_denominator
    
    @property
    def u_second_deriv(self):
        distances = self.particles.distances
        return -self.a*(2*distances - self.a)/(self._u_deriv_denominator**2)

    @property
    def _u_deriv_denominator(self):
        distances = self.particles.distances
        overlapping = self.particles.overlapping
        denom = 1/(self.a*distances - distances**2)
        denom[overlapping] = np.inf

        return denom

    @property
    def onebody_part(self):
        return np.prod(single_particle_function())

    @property
    def ext_potential(self):
        elliptic_distances = np.linalg.norm(self.omegas*self.particles.positions,
                                            axis=-1)
        return 0.5*self.particles.mass*np.sum(elliptic_distances**2)

    @property
    def int_potential(self):
        overlapping = self.particles.overlapping.copy()
        overlapping[range(self.num_particles), range(self.num_particles)] = False
        return np.inf if np.any(overlapping) else 0

    @property
    def kinetic_energy(self):
        particles = self.particles
        n = self.particles.num_particles
        d = self.particles.num_dimensions

        # Unary terms
        kinetic = self.single_particle_function.relative_laplacian.sum(1)

        if self.a == 0:
            return (-0.5*self.hbar/self.particles.mass)*kinetic.sum()
        
        # Intermediate terms used for the pairwise terms
        relative_gradient = self.single_particle_function.relative_gradient
        rec_distances = 1/particles.distances
        rec_distances[range(n), range(n)] = 0

        normalised_distance = particles.differences*rec_distances[..., np.newaxis]
        weighted_distance = self.u_derivatives[..., np.newaxis]*normalised_distance

        # Pairwise terms
        kinetic += np.einsum('ikj,ij->i', weighted_distance, relative_gradient)
        
        difference_sums = np.zeros(n)
        self._difference_sums(weighted_distance, n, difference_sums)
        kinetic += difference_sums
        # This could be 
        # kinetic += np.einsum('ikl,kjl->k', weighted_distance, weighted_distance)
        # but that i slower.

        u_primes = self.u_derivatives
        double_prime_term = self.u_second_deriv + 2*u_primes*rec_distances
        kinetic += double_prime_term.sum(axis=1)
        
        kinetic *= -0.5*self.hbar/self.particles.mass

        return kinetic.sum()

    @staticmethod
    @jit
    def _difference_sums(differences, n, output):
        for i in range(n):
            output[i] = np.sum(differences[i]@differences[i].T)

    @property
    def local_energy(self):
        return self.kinetic_energy + self.int_potential + self.ext_potential

    @property
    def quantum_force(self):
        u_derivs = self.u_derivatives[..., np.newaxis]
        rel_diffs = self.particles.relative_differences

        return 2*self.single_particle_function.relative_gradient + \
                np.sum((u_derivs*rel_diffs), axis=1)

    @property
    def log_correlation_part(self):
        return 0.5*self.u_matrix.sum()

    def __call__(self):
        single_particle_part = self.single_particle_function().astype('float128').prod()
        if self.a == 0:
            return single_particle_part
        f = self.f
        f[range(self.num_particles), range(self.num_particles)] = 1
        interactions = f.prod()
        return single_particle_part*interactions
    
    def __getstate__(self):
        particles_pickle = pickle.dumps(self.particles)
        single_particle_dict = {
            k: v for k, v in self.single_particle_function.__dict__.items()
                     if k != 'particles'
        }
        self_dict = {
            k: v for k, v in self.__dict__.items()
                if k not in ('single_particle_function', 'particles')
        }
        return particles_pickle, single_particle_dict, self_dict
    
    def __setstate__(self, state):
        self.__dict__ = state[2]
        self.particles = pickle.loads(state[0])
        self.single_particle_function = SingleParticleFunction(self.particles)
        self.single_particle_function.__dict__ = state[1]
        self.single_particle_function = SingleParticleFunction(self.particles)


if __name__ == '__main__':
    pass


