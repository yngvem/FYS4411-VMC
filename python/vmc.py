"""

"""

__author__ = 'Yngve Mardal Moe'
__email__ = 'yngve.m.moe@gmail.com'

import numpy as np
from pickle import dumps, loads
from copy import deepcopy
from particles import Particles
from wave_function import WaveFunction
import sampler


class VMC:
    def __init__(self, sampler_type, wf_params, particles_params):
        if type(sampler_type) is str:
            sampler_type = getattr(sampler, sampler_type)
        self.sampler_type = sampler_type
        self.particles_params = particles_params
        self.wave_function_params = wf_params
        self.ideal_step_sizes = {}
    
    def _create_sampler(self, alpha):
        parts = Particles(**self.particles_params)
        wave_function = WaveFunction(
            parts,
            **self.wave_function_params,
            alpha=alpha
        )
        return self.sampler_type(step_size=0.1, wave_function=wave_function)
    
    def find_ideal_step_size(self, alpha, ideal_ar, ar_error=0.05,
                             max_ar_iterations=20, ar_int_iters=500):
        """Performs MCMC iterations with the specified alpha.

        Parameters:
        -----------
        alpha : float
            Variational parameter.
        ideal_ar : float
            Ideal acceptance ratio for MCMC computations
        ar_error : float (=0.05)
            Acceptable deviation between the actual acceptance ratio and the
            ideal acceptance ratio.
        max_ar_iterations : int (=20)
            Maximum number of bisection iterations to perform when computing
            the ideal step size.
        ar_int_iters : int (=500)
            Number of iterations used to compute the acceptance ratio
        seed : int (optional)
            Random seed.
        """
        if (alpha, ideal_ar) in self.ideal_step_sizes:
            return self.ideal_step_sizes[alpha, ideal_ar]
        sampler = self._create_sampler(alpha)
        sampler.find_ideal_step_size(
            ideal=ideal_ar,
            delta=ar_error,
            num_iterations=ar_int_iters,
            max_bisection_steps=max_ar_iterations
        )
        self.ideal_step_sizes[alpha, ideal_ar] = sampler.step_size
        return sampler.step_size
    
    def perform_mc_iterations(self, alpha, num_iterations, ideal_ar,
                              ar_error=0.05, max_ar_iterations=20,
                              ar_int_iters=500, warm_up=1000,
                              return_positions=False, return_wf=False,
                              seed=None):
        """Performs MCMC iterations with the specified alpha.

        Parameters:
        -----------
        alpha : float
            Variational parameter.
        num_iterations : int
            Number of MCMC iterations to perform.
        ideal_ar : float
            Ideal acceptance ratio for MCMC computations
        ar_error : float (=0.05)
            Acceptable deviation between the actual acceptance ratio and the
            ideal acceptance ratio.
        max_ar_iterations : int (=20)
            Maximum number of bisection iterations to perform when computing
            the ideal step size.
        ar_int_iters : int (=500)
            Number of iterations used to compute the acceptance ratio
        seed : int (optional)
            Random seed.
        """
        np.random.seed(seed)
        sampler = self._create_sampler(alpha)
        sampler.step_size = self.find_ideal_step_size(
            alpha=alpha,
            ideal_ar=ideal_ar,
            ar_error=ar_error,
            ar_int_iters=ar_int_iters,
            max_ar_iterations=max_ar_iterations
            )
        sampler.warm_up(warm_up)
        return sampler.perform_mc_iterations(
            num_iterations,
            return_positions=return_positions,
            return_wf=return_wf
        )
    
    def compute_gradient_and_energy(self, alpha, num_iterations, ideal_ar,
                                    ar_error=0.5, max_ar_iterations=20,
                                    ar_int_iters=500, warm_up=1000,
                                    seed=None):
        """Compute the derivative of the energy with respect to alpha.

        Parameters:
        -----------
        alpha : float
            Variational parameter.
        num_iterations : int
            Number of MCMC iterations to perform.
        ideal_ar : float
            Ideal acceptance ratio for MCMC computations
        ar_error : float (=0.05)
            Acceptable deviation between the actual acceptance ratio and the
            ideal acceptance ratio.
        max_ar_iterations : int (=20)
            Maximum number of bisection iterations to perform when computing
            the ideal step size.
        ar_int_iters : int (=500)
            Number of iterations used to compute the acceptance ratio
        seed : int (optional)
            Random seed.
        """
        energies, positions, wfs = self.perform_mc_iterations(
            alpha=alpha,
            num_iterations=num_iterations,
            return_positions=True,
            return_wf=True,
            ideal_ar=ideal_ar,
            ar_error=ar_error,
            max_ar_iterations=max_ar_iterations,
            ar_int_iters=ar_int_iters
        )

        if positions.shape[1] == 3:
            beta_vector = np.array([[1, 1, self.wave_function_params['beta']]])
            positions_sum = np.array([(p*beta_vector).sum() for p in positions])
        else:
            positions_sum = np.array([p.sum() for p in positions])

        energy = energies.mean()
        psi_deriv = -2*alpha*positions_sum*wfs
        

        first_part = np.mean(psi_deriv*energies)
        second_part = np.mean(psi_deriv)*energy
        derivative = 2*(first_part-second_part)

        return derivative, energy
        

    def __getstate__(self):
        self_dict = deepcopy(self.__dict__)
        self_dict['sampler_type'] = dumps(self.sampler_type)
        return self_dict
        

    def __setstate__(self, state):
        self.__dict__ = state
        self.sampler_type = loads(self.sampler_type)




if __name__ == '__main__':
    pass


