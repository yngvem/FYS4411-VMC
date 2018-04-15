"""

"""


__author__ = 'Yngve Mardal Moe'
__email__ = 'yngve.m.moe@gmail.com'


import numpy as np


class BaseSampler:
    def __init__(self, wave_function, step_size, seed=None):
        np.random.seed(seed)
        self.step_size = step_size
        self.wave_function = wave_function
        self._P = wave_function()**2
        self._old_P = self._P
        self.energies = np.array([])
        self.current_energy = wave_function.local_energy 
        self.acceptance_list = []

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
        r = np.random.random() < self._rejection_criteria()
        self.acceptance_list.append(r)
        return r

    def _reject_perturbation(self):
        """Reject proposed perturbation.
        """
        self.wave_function.particles.reset_positions()
        self._P = self._old_P

    def warm_up(self, num_steps=1000):
        """Perform a given set of iterations without storing anything.
        """
        for _ in range(num_steps):
            self.single_step()

    def compute_local_energy(self, num_steps):
        """Perform MCMC iterations and return the estimated energy and variance.
        """
        energies = self.perform_mc_iterations(num_steps)
        return energies.mean(), energies.std()**2

    def perform_mc_iterations(self, num_steps, save_frequency=1,
                              return_positions=False, return_wf=False):
        """Perform the given amount of MC cycles and store the energies.
        """
        energies = np.zeros(num_steps//save_frequency)
        positions = []
        wf_evals = []
        for i in range(num_steps):
            self.single_step()
            if i % save_frequency == 0:
                energies[i//save_frequency] = self.current_energy

                if return_positions:
                    positions.append(self.wave_function.particles.positions.copy())
                if return_wf:
                    wf_evals.append(self.wave_function())

        positions = np.array(positions)
        wf_evals = np.array(wf_evals)
        self.energies = np.concatenate((self.energies, energies))

        if return_positions and return_wf:
            return energies, positions, wf_evals
        elif return_positions:
            return energies, positions
        elif return_wf:
            return energies, wf_evals
        else:
            return energies

    @property
    def P(self):
        return self._P

    @P.setter
    def P(self, value):
        self._old_P = self._P
        self._P = value

    @property
    def acceptance_ratio(self):
        return np.mean(self.acceptance_list)
    
    def compute_acceptance_ratio(self, num_iterations=100):
        self.acceptance_list = []
        self.warm_up(num_iterations)
        return self.acceptance_ratio
    
    
    def _ideal_step_size(self, delta, ideal):
        return np.abs(self.acceptance_ratio-ideal) < delta
    
    def find_ideal_step_size(self, delta=0.1, ideal=0.5, num_iterations=100, 
                             max_bisection_steps=20, verbose=False):
        """Find step size which has the specified acceptance ratio.
        
        Assumes that acceptance ratio decreases monotonically with step size
        and use the bisection algorithm to find a step size with an acceptable
        acceptance ratio.
        
        A maximum of 20 splits is performed and if no acceptable ratios are 
        computed a runtime error is raised.
        """
        step_size1 = self.step_size
        step_size2 = step_size1
        
        if self._ideal_step_size(delta, ideal):
            return self.step_size, self.acceptance_ratio
        
        while self.compute_acceptance_ratio(num_iterations) < ideal:
            self.step_size = step_size1
            step_size1 *= 0.1
            
        ar1 = self.compute_acceptance_ratio(num_iterations)
        
        if self._ideal_step_size(delta, ideal):
            return self.step_size, self.acceptance_ratio
        
        if verbose:
            print(f'Found left bound, {step_size1}')
            print(f'  With ar {self.acceptance_ratio}')
        
        self.step_size = step_size2
        while self.compute_acceptance_ratio(num_iterations) > ideal:
            self.step_size = step_size2
            step_size2 *= 10
            ar2 = self.compute_acceptance_ratio(num_iterations)
        
        if verbose:
            print(f'Found right bound, {step_size2}')
            print(f'  With ar {self.acceptance_ratio}')
            print('\nStarting iterations now')
       
        n = 0
        while not self._ideal_step_size(delta, ideal):
            self.step_size = (step_size2 + step_size1)/2
            ar = self.compute_acceptance_ratio(num_iterations)
                
            if np.sign((ar-ideal)*(ar1-ideal)) < 0:
                step_size2 = self.step_size
            else:
                ar1 = ar
                step_size1 = self.step_size
            
            n += 1
            if n > max_bisection_steps:
                raise RuntimeError('No acceptable step size could be found')
            if verbose: print(f'    Iteration number {n} completed')

        if verbose:
            print(f'Found acceptable step size after {n} iterations')

        return self.step_size, self.acceptance_ratio

    
class MetropolisSampler(BaseSampler):
    def __init__(self, wave_function, step_size, seed=None):
        self._current_particle = 0
        self._num_particles = wave_function.particles.num_particles
        self._num_dimensions = wave_function.particles.num_dimensions
        
        super().__init__(wave_function, step_size=step_size, seed=seed)

    def _rejection_criteria(self):
        return self.P/self._old_P

    def _propose_perturbation(self):
        perturbation = self.step_size*np.random.randn(self._num_dimensions)
        self.wave_function.particles.perturb_particle(perturbation,
                                                      self._current_particle)
        self._current_particle = (self._current_particle+1)%self._num_particles
        
            


class ImportanceSampler(BaseSampler):
    def __init__(self, wave_function, step_size, seed=None):
        self.D = wave_function.hbar/(2*wave_function.particles.mass)
        self._F = wave_function.quantum_force
        self._old_F = self._F.copy()

        self._greens = 0

        self._current_particle = 0
        self._num_particles = wave_function.particles.num_particles
        self._num_dimensions = wave_function.particles.num_dimensions
        super().__init__(wave_function, step_size=step_size, seed=seed)

    def _rejection_criteria(self):
        return self.greens_ratio*self.P/(self._old_P)

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
        self._F = self._old_F
        super()._reject_perturbation()

    @property
    def greens_ratio(self):
        return self._greens

    @greens_ratio.setter
    def greens_ratio(self, value):
        self._greens = value
    
    @property
    def F(self):
        return self._F

    @F.setter
    def F(self, value):
        self._old_F = self._F
        self._F = value
    
    @property
    def timestep(self):
        return self.step_size
    
    @timestep.setter
    def timestep(self, value):
        self.step_size = value

    def _update_greens(self, perturbation):
        self.F = self.wave_function.quantum_force
        delta_F = self.F[self._current_particle] - self._old_F[self._current_particle]
        F_sum = self.F[self._current_particle] + self._old_F[self._current_particle]

        greens = -0.5*F_sum * (0.5*self.D*self.timestep*delta_F + perturbation)
        self.greens_ratio = np.exp(greens.sum())

if __name__ == '__main__':
    pass
