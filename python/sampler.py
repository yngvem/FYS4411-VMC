"""

"""


__author__ = 'Yngve Mardal Moe'
__email__ = 'yngve.m.moe@gmail.com'


import numpy as np


class BaseSampler:
    def __init__(self, wave_function, debug=False):
        self.wave_function = wave_function
        self._P = wave_function()**2
        self._old_P = self._P
        self.energies = np.array([])
        self.current_energy = wave_function.local_energy 
        self.debug = debug 
        self.debug_list = []

    def _rejection_criteria(self): 
        """Defined in subclasses. Returns probability of keeping current position 
        """
        raise NotImplementedError('The rejection criteria should be implemented '
                                  'before using the sampler.')

    def _propose_perturbation(self):
        """Defined in subclasses. Proposes a new position of the particles.
        """ 
        raise NotImplementedError('The perturbations should be implemented '
                                  'before using the sampler.')

    def single_step(self):
        """Perform a single MCMC iteration. 
        """
        self._propose_perturbation() 
        self.P = self.wave_function()**2 
        if self._keep_position: 
            self.current_energy = self.wave_function.local_energy 
        else: 
            self._reject_perturbation() 
        
        return self.current_energy

    @property
    def _keep_position(self):
        """Whether the current position should be kept or not.
        """
        if self.debug:
            r = np.random.random() < self._rejection_criteria()
            self.debug_list.append((self._rejection_criteria(), r))
            return r

        return np.random.random() < self._rejection_criteria()

    def _reject_perturbation(self):
        """Reject proposed perturbation.
        """
        self.wave_function.particles.reset_positions()
        self._P = self._old_P

    def warm_up(self, num_steps=1000):
        """Perform a given set of iterations without storing anything.
        """
        for _ in num_steps:
            self.single_step()

    def compute_local_energy(self, num_steps):
        """Perform MCMC iterations and return the estimated energy and variance.
        """
        energies = self.perform_mc_iterations(num_steps)
        return energies.mean(), energies.std()**2

    def perform_mc_iterations(self, num_steps, save_frequency=1):
        """Perform the given amount of MC cycles and store the energies.
        """
        energies = np.zeros(num_steps//save_frequency)
        for i in range(num_steps):
            self.single_step()
            if i % save_frequency == 0:
                energies[i//save_frequency] = self.current_energy
        self.energies = np.concatenate((self.energies, energies))
        return energies

    @property
    def P(self):
        return self._P

    @P.setter
    def P(self, value):
        self._old_P = self._P
        self._P = value


class MetropolisSampler(BaseSampler):
    def __init__(self, wave_function, step_size, debug=False):
        self.step_size = step_size
        self._current_particle = 0
        self._num_particles = wave_function.particles.num_particles
        self._num_dimensions = wave_function.particles.num_dimensions
        
        super().__init__(wave_function, debug=debug)

    def _rejection_criteria(self):
        return self.P/self._old_P

    def _propose_perturbation(self):
        perturbation = self.step_size*np.random.randn(self._num_dimensions)
        self.wave_function.particles.perturb_particle(perturbation,
                                                      self._current_particle)
        self._current_particle = (self._current_particle+1)%self._num_particles


class ImportanceSampler(BaseSampler):
    def __init__(self, wave_function, timestep, debug=False):
        self.D = wave_function.hbar/(2*wave_function.particles.mass)
        self.timestep = timestep
        self._F = wave_function.quantum_force
        self._old_F = self._F.copy()

        self._greens = 0

        self._current_particle = 0
        self._num_particles = wave_function.particles.num_particles
        self._num_dimensions = wave_function.particles.num_dimensions
        super().__init__(wave_function, debug=debug)

    def _rejection_criteria(self):
        return self.greens*self.P/(self._old_P)

    def _propose_perturbation(self):
        random_part  = np.random.randn(self._num_dimensions)*np.sqrt(self.timestep)
        

        quantum_force = self.wave_function.quantum_force[self._current_particle]
        quantum_force *= self.D*self.timestep
        
        perturbation = random_part + quantum_force
        self.wave_function.particles.perturb_particle(perturbation,
                                                      self._current_particle)
        self._update_greens(perturbation)
        self._current_particle = (self._current_particle+1)%self._num_particles

    def _reject_perturbation(self):
        self._greens = self._old_greens
        self._F = self._old_F
        super()._reject_perturbation()

    @property
    def greens(self):
        return self._greens

    @greens.setter
    def greens(self, value):
        self._old_greens = self._greens
        self._greens = value
    
    @property
    def F(self):
        return self._F

    @F.setter
    def F(self, value):
        self._old_F = self._F
        self._F = value

    def _update_greens(self, perturbation):
        self.F = self.wave_function.quantum_force
        delta_F = self.F[self._current_particle] - self._old_F[self._current_particle]
        F_sum = self.F[self._current_particle] + self._old_F[self._current_particle]

        greens = 0.5*F_sum * (0.5*self.D*self.timestep*delta_F + perturbation)
        self.greens = np.exp(greens.sum())

if __name__ == '__main__':
    pass
