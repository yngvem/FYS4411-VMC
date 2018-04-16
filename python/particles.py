"""
"""


__author__ = 'Yngve Mardal Moe'
__email__ = 'yngve.m.moe@gmail.com'


import numpy as np
from scipy.spatial.distance import pdist, squareform
from numba import jit


class Particles:
    """Wrapper for the particle positions.
    """
    def __init__(self, num_particles, num_dimensions, mass, diameter):
        self._num_particles = num_particles
        self._num_dimensions = num_dimensions
        self.mass = mass
        self.diameter = diameter

        self._positions = np.random.randn(num_particles, num_dimensions)
        self._old_positions = self._positions.copy()

        # Difference matrices
        self._distances = squareform(pdist(self._positions))
        self._rec_distances = 1/self._distances
        self._rec_distances[range(num_particles), range(num_particles)] = 0
        self._differences = np.zeros((num_particles, num_particles, num_dimensions))
        self._rel_differences = self._differences*self._rec_distances[..., np.newaxis]
        self._overlapping = self._distances < diameter

        self._resolve_overlapping_particles()
        self._positions_changed = True

    def _resolve_overlapping_particles(self):
        overlapping = self._overlapping.copy()
        overlapping[range(self.num_particles), range(self.num_particles)] = False
        if not np.any(overlapping):
            return
        
        self._positions *= 0
        var = 1
        for i in range(self.num_particles):
            self._positions[i] = var*np.random.randn(self.num_dimensions)
            while np.any(pdist(self._positions[:i])):
                var += 1
                self._positions[i] = var*np.random.randn(self.num_dimensions)

    def set_position(self, new_position, particle_num=None):
        """Change the position of one or all of the particles.
        Parameters
        ----------
        new_position : Numpy array
            The new position vector (or matrix)
        particle_num : int or None
            The index of the particle to move. If None, all the particle
            positions are changed.
        """
        self._positions_changed = True
        self._num_pos_changes += 1
        self.old_positions[:, :] = self._positions

        if particle_num is None:
            self._positions[:, :] = new_position[:, :]
        else:
            self._positions[particle_num] = new_position

    def perturb_particle(self, perturbation, particle_num):
        """Perturb the position of one or all of the particles.
        Parameters
        ----------
        perturbation : Numpy array
            The perturbation to add to the position vector (or matrix).
        particle_num : int or None
            The index of the particle to perturb. If None, all the particles
            are perturbed.
        """
        self._positions_changed = True
        self._old_positions[:, :] = self._positions
        
        if particle_num is None:
            self._positions += perturbation
        else:
            self._positions[particle_num] += perturbation

    def reset_positions(self):
        """Reset positions to previous timestep.
        """
        self._positions_changed = True 
        self._positions[:, :] = self._old_positions

    def _update_distances(self):
        """Update all differences and distances.
        """
        self._positions_changed = False

        self._distances = squareform(pdist(self._positions))
        self._rec_distances = 1/self._distances
        self._rec_distances[range(self.num_particles), range(self.num_particles)] = 0
        self._overlapping = self._distances <= self.diameter
        self._differences *= 0
        self._compute_differences(self._positions, self._differences)
        self._rel_differences = self._differences*self._rec_distances[..., np.newaxis]

    @staticmethod
    @jit
    def _compute_differences(x, diff):
        """Compute all pairwise position differences.
        """
        for i in range(x.shape[0]):
            diff[i, i+1:] = x[i] - x[i+1:]
        diff -= np.rollaxis(diff, 1)

    @property
    def positions(self):
        """The position matrix of the particles.
        """
        return self._positions

    @property
    def distances(self):
        """The pairwise distances of all particles.
        """
        if self._positions_changed:
            self._update_distances()
        return self._distances
    
    @property
    def differences(self):
        """The pairwise position differences (as vectors) of all particles.
        """
        if self._positions_changed:
            self._update_distances()
        return self._differences

    @property
    def relative_differences(self):
        return self._rel_differences

    @property
    def reciprocal_distances(self):
        return self._rec_distances

    @property
    def overlapping(self):
        """Boolean matrix indicating whether any particles overlap.
        """
        if self._positions_changed:
            self._update_distances()
        return self._overlapping

    @property
    def num_dimensions(self):
        return self._num_dimensions

    @property
    def num_particles(self):
        return self._num_particles

if __name__ == '__main__':
    pass
