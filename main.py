import numpy as np
from seir_engine import SEIRModelEngine
from Visualizer import SEIRVisualizer
from scipy.integrate import solve_ivp

# 1. Initialize the Engine and Visualizer
engine = SEIRModelEngine(N=1_000_000.0)
viz = SEIRVisualizer(output_dir='figures')

# 2. Set Initial Conditions (S, E, I, R) 
I0 = 8760.0
E0 = 2.5 * I0  # Modeling the latent parasite reservoir 
y0 = [1_000_000.0 - E0 - I0, E0, I0, 0.0]
observed_times = np.array([7, 14, 21, 28])
observed_cases = np.array([1800, 2200, 2450, 2580])

# 3. Run Bayesian Inference 
print("Starting NUTS Sampling...")
trace = engine.run_inference(observed_times, observed_cases, y0)

# 4. Extract Results and Run Deterministic Simulations 
post_beta = trace.posterior["beta_0"].mean().item()
post_eta = trace.posterior["eta"].mean().item()
t_eval = np.linspace(0, 210, 210)

res = solve_ivp(engine.seir_scipy, [0, 210], y0, 
                args=(post_beta, post_eta, engine.gamma_base), t_eval=t_eval)

# Simulate 30% reduction via intervention (e.g., bed-nets) 
res_nets = solve_ivp(engine.seir_scipy, [0, 210], y0, 
                     args=(0.7 * post_beta, post_eta, engine.gamma_base), t_eval=t_eval)

# 5. Generate and Save Figures 
print("Generating figures...")
viz.plot_all(t_eval, res, res_nets, observed_times, observed_cases, 
             trace, post_beta, post_eta, engine.gamma_base, engine, y0)
print("Done! Check the 'figures/' folder.")