import numpy as np
import pymc as pm
import pytensor.tensor as pt
from pymc.ode import DifferentialEquation
from scipy.integrate import solve_ivp

class SEIRModelEngine:
    def __init__(self, N=1_000_000.0):
        self.N = N
        self.gamma_base = 1/7.0  # Recovery rate 
        self.sigma = 1/12.0      # Incubation rate 

    def seir_ode_pymc(self, y, t, p):
        """ODE system for PyMC HMC samplers."""
        S, E, I, R = y
        beta_0, eta = p[0], p[1]
        # Seasonal forcing: beta fluctuates to model mosquito density peaks 
        beta_t = beta_0 * (1 + eta * pt.cos(2 * np.pi * (t - 105) / 365))
        return [-beta_t*S*I/self.N, beta_t*S*I/self.N - self.sigma*E, 
                self.sigma*E - self.gamma_base*I, self.gamma_base*I]

    def seir_scipy(self, t, y, b0, e, g):
        """Standard ODE system for deterministic validation."""
        bt = b0 * (1 + e * np.cos(2 * np.pi * (t - 105) / 365))
        return [-bt*y[0]*y[2]/self.N, bt*y[0]*y[2]/self.N - self.sigma*y[1], 
                self.sigma*y[1] - g*y[2], g*y[2]]

    def run_inference(self, observed_times, observed_cases, y0):
        """Performs Bayesian inference using NUTS."""
        seir_model_wrapper = DifferentialEquation(
            func=self.seir_ode_pymc, times=observed_times, n_states=4, n_theta=2, t0=0
        )
        
        with pm.Model() as malaria_model:
            beta_0 = pm.Lognormal("beta_0", mu=np.log(0.4), sigma=0.3)
            eta = pm.Beta("eta", alpha=2, beta=5)
            phi = pm.HalfNormal("phi", sigma=10)

            sol = seir_model_wrapper(y0=y0, theta=[beta_0, eta])
            pm.NegativeBinomial("Y_obs", mu=sol[:, 2], alpha=phi, observed=observed_cases)
            
            trace = pm.sample(1000, tune=1000, target_accept=0.9, chains=2)
            pm.sample_posterior_predictive(trace, extend_inferencedata=True)
            
        return trace