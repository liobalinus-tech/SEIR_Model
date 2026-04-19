import matplotlib.pyplot as plt
import seaborn as sns
import arviz as az
import numpy as np
import os
from scipy.integrate import solve_ivp  # Essential import for sensitivity analysis 

class SEIRVisualizer:
    def __init__(self, output_dir='figures'):
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        plt.style.use('seaborn-v0_8-whitegrid')

    def plot_all(self, t_eval, res, res_nets, observed_times, observed_cases, trace, post_beta, post_eta, gamma_base, engine, y0, output_dir=None):
        """Generates and saves all clinical visuals for the SEIR model."""
        
        # FIG 1: Mass Action Balance 
        plt.figure(figsize=(8, 5))
        plt.stackplot(t_eval, res.y[2], res.y[1], res.y[0], 
                      labels=['I', 'E', 'S'], 
                      colors=['#e74c3c', '#f39c12', '#3498db'], alpha=0.7)
        plt.title("FIG 1: Full SEIR State Transitions", fontweight='bold')
        plt.xlabel("Days")
        plt.ylabel("Individuals")
        plt.legend()
        plt.savefig(f'{self.output_dir}/fig1_seir_states.png', dpi=300)
        plt.close()

        # FIG 2: Uncertainty Quantification 
        plt.figure(figsize=(8, 5))
        az.plot_hdi(observed_times, trace.posterior_predictive["Y_obs"], 
                    color="gray", fill_kwargs={"alpha": 0.3, "label": "94% HDI"})
        plt.plot(t_eval, res.y[2], color='black', label="Posterior Mean")
        plt.scatter(observed_times, observed_cases, color='red', label="Observed Data", zorder=5)
        plt.title("FIG 2: Model Calibration & Uncertainty", fontweight='bold')
        plt.legend()
        plt.savefig(f'{self.output_dir}/fig2_calibration.png', dpi=300)
        plt.close()

        # FIG 3: Reproduction Number Dynamics 
        plt.figure(figsize=(8, 5))
        Rt = (post_beta * (1 + post_eta * np.cos(2*np.pi*(t_eval-105)/365))) / gamma_base
        plt.plot(t_eval, Rt, color='#8e44ad', linewidth=2)
        plt.axhline(1, color='black', linestyle='--')
        plt.title("FIG 3: Seasonal Effective Reproduction Number R(t)", fontweight='bold')
        plt.ylabel("R(t)")
        plt.xlabel("Days")
        plt.savefig(f'{self.output_dir}/fig3_rt_dynamics.png', dpi=300)
        plt.close()

        # FIG 4: Phase Portrait 
        plt.figure(figsize=(8, 5))
        plt.plot(res.y[0], res.y[2], color='#2c3e50', lw=2)
        plt.title("FIG 4: Epidemiological Phase Portrait", fontweight='bold')
        plt.xlabel("Susceptible (S)")
        plt.ylabel("Infectious (I)")
        plt.savefig(f'{self.output_dir}/fig4_phase_portrait.png', dpi=300)
        plt.close()

        # FIG 5: Strategic Sensitivity Analysis 
        print("Calculating sensitivity matrix (this may take a moment)...")
        plt.figure(figsize=(8, 5))
        b_range = np.linspace(0.2, 0.6, 15)
        g_range = np.linspace(0.1, 0.2, 15)
        peaks = np.zeros((len(b_range), len(g_range)))
        
        for i, b in enumerate(b_range):
            for j, g in enumerate(g_range):
                # Fixed: ensure solve_ivp has access to engine and y0 passed in plot_all
                sol = solve_ivp(engine.seir_scipy, [0, 210], y0, 
                                args=(b, post_eta, g), t_eval=t_eval)
                peaks[i, j] = np.max(sol.y[2])

        plt.imshow(peaks, origin='lower', extent=[0.1, 0.2, 0.2, 0.6], 
                   aspect='auto', cmap='YlOrRd')
        plt.colorbar(label="Peak Infection Count")
        plt.title("FIG 5: Strategic Sensitivity Analysis", fontweight='bold')
        plt.xlabel("Recovery Rate (Gamma)")
        plt.ylabel("Transmission Rate (Beta)")
        plt.savefig(f'{self.output_dir}/fig5_sensitivity.png', dpi=300, bbox_inches='tight')
        plt.close()

        # FIG 6: Policy Simulation 
        plt.figure(figsize=(8, 5))
        plt.fill_between(t_eval, res.y[2], res_nets.y[2], color='green', alpha=0.2, label="Averted Cases")
        plt.plot(t_eval, res.y[2], 'k--', alpha=0.4, label="Baseline")
        plt.plot(t_eval, res_nets.y[2], color='green', label="Bed-Net Strategy (30% reduction)")
        plt.title("FIG 6: Mitigation Strategy Impact", fontweight='bold')
        plt.legend()
        plt.savefig(f'{self.output_dir}/fig6_intervention.png', dpi=300)
        plt.close()