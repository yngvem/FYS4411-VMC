"""

"""

__author__ = 'Yngve Mardal Moe'
__email__ = 'yngve.m.moe@gmail.com'

import numpy as np
from numba import jit, float64, void, int64
from particles import Particles


def compute_wf_differences(wf, h=1e-5):
    n = wf.particles.num_particles
    d = wf.particles.num_dimensions
    forward_differences = np.zeros((n, d))
    backward_differences = np.zeros((n, d))

    for i in range(n):
        for j in range(d):
            perturbation = np.zeros(d)
            perturbation[j] = h

            wf.particles.perturb_particle(perturbation, i)
            forward_differences = wf()
            wf.particles.reset_positions()

            wf.particles.perturb_particle(-perturbation, i)
            backward_differences = wf()
            wf.particles.reset_positions

    return forward_differences, backward_differences


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
        return self.relative_laplacian*self()[..., np.newaxis]

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

    def numerical_gradient(self, h=1e-8):
        fd, bd = compute_wf_differences(self, h)
        return (fd - bd)/(2*h)

    def numerical_laplacian(self, h=1e-5):
        fd, bd = compute_wf_differences(self, h)
        return (fd - 2*self() + bd)/(2*h)
        

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

    def f(self, r):
        return 1-(self.a/r)

    def u(self, r):
        return np.log(self.f(r))

    def u_derivatives(self, r):
        return self.a/(r*(r-self.a))
    
    def u_second_deriv(self, r):
        return -self.a*(2*r - self.a)/((r*(r - self.a))**2)

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

        # Unary terms of Laplacian
        kinetic = self.single_particle_function.relative_laplacian.sum()

        # print(f'Single particle part: {kinetic}')
        if self.a == 0:
            return (-0.5*self.hbar/self.particles.mass)*kinetic
        
        # Intermediate terms regarding single particle function 
        relative_spf_gradient = self.single_particle_function.relative_gradient

        # Intermediate terms regarding the u-functions
        u_derivs = self.u_derivatives(self.particles.distances)
        u_derivs[range(self.num_particles), range(self.num_particles)] = 0
        u_second_deriv = self.u_second_deriv(self.particles.distances)
        u_second_deriv[range(self.num_particles), range(self.num_particles)] = 0

        # Intermediate terms regarding distances
        rec_distances = 1/particles.distances
        rec_distances[range(n), range(n)] = 0
        normalised_distance = particles.differences*rec_distances[..., np.newaxis]
        weighted_distance = u_derivs[..., np.newaxis]*normalised_distance

        # Pairwise terms of Laplacian
        # from pdb import set_trace; set_trace()
        kinetic -= np.einsum('ikj,ij->i', weighted_distance, relative_spf_gradient).sum()
        # print(f'Gradient pairwise part: {kinetic}')
        
        kinetic -= self._triple_sums()
        # print(f'First pairwise part: {kinetic}')

        double_prime_term = u_second_deriv + 2*u_derivs*rec_distances
        kinetic += double_prime_term.sum()
        # print(f'Second pairwise part: {kinetic}')
        
        kinetic *= -0.5*self.hbar/self.particles.mass

        return kinetic

    def _triple_sums(self):
        particles = self.particles
        n = self.particles.num_particles

        # Intermediate terms regarding the u-functions
        u_derivs = self.u_derivatives(self.particles.distances)
        u_derivs[range(self.num_particles), range(self.num_particles)] = 0
        
        # Intermediate terms regarding positions
        rec_distances = 1/particles.distances
        rec_distances[range(n), range(n)] = 0
        normalised_distance = particles.differences*rec_distances[..., np.newaxis]
        weighted_distance = u_derivs[..., np.newaxis]*normalised_distance

        # Compute the difference_sums
        difference_sums = np.zeros(n)
        self._compute_triple_sums(weighted_distance, n, difference_sums)
        # This could be 
        # kinetic += np.einsum('ikl,kjl->k', weighted_distance, weighted_distance)
        # but that i slower.
        return difference_sums.sum()

    @staticmethod
    @jit
    def _compute_triple_sums(differences, n, output):
        for i in range(n):
            output[i] = np.sum(differences[i]@differences[i].T)

    @property
    def local_energy(self):
        return self.kinetic_energy + self.int_potential + self.ext_potential

    @property
    def quantum_force(self):
        u_derivs = self.u_derivatives(self.particles.distances)[..., np.newaxis]
        u_derivs[range(self.num_particles), range(self.num_particles)] = 0
        rel_diffs = self.particles.relative_differences
        spf_grad = self.single_particle_function.relative_gradient

        return 2*spf_grad + 4*np.sum((u_derivs*rel_diffs), axis=1)

    @property
    def log_correlation_part(self):
        return 0.5*self.u_matrix.sum()

    def __call__(self):
        single_particle_part = self.single_particle_function().astype('float128').prod()
        
        if self.a == 0:
            return single_particle_part
        
        f = self.f(self.particles.distances)
        f[self.particles.overlapping] = 0
        f[range(self.num_particles), range(self.num_particles)] = 1
        interactions = f.prod()
        return single_particle_part*interactions
    

    def numerical_quantum_force(self, h=1e-8):
        fd, bd = compute_wf_differences(self, h)
        gradient = (fd - bd)/(2*h)
        return 2*gradient/self()
    
    def numerical_kinetic_energy(self, h=1e-5):
        fd, bd = compute_wf_differences(self, h)
        laplacian = np.sum((fd - 2*self() + bd)/(h**2))
        relative_laplacian = laplacian/self()
        energy_factor = -0.5*self.hbar/self.particles.mass

        return relative_laplacian*energy_factor

    def numerical_u_derivatives(self, r, h=1e-8):
        fd = self.u(r + h)
        bd = self.u(r - h)
        return (fd - bd)/(2*h)

    def numerical_u_second_derivs(self, r, h=1e-5):
        u = self.u(r)
        fd = self.u(r + h)
        bd = self.u(r - h)
        return (fd - 2*u + bd)/(h**2)


if __name__ == '__main__':
    pass


