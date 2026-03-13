# SEIR_Model
# Bayesian Probabilistic Modeling of Seasonal Malaria Dynamics: An SEIR-Inference Framework for Public Health Policy in Northern Namibia

This repository presents a rigorous computational study of malaria epidemiology. The project utilizes a deterministic **SEIR (Susceptible-Exposed-Infectious-Recovered)** framework, augmented by **Bayesian Markov Chain Monte Carlo (MCMC)** methods to calibrate parameters against stochastic real-world data.

## 🔬 Research Overview
Conventional epidemiological models often rely on fixed parameters that fail to account for the stochastic nature of disease transmission. This project addresses this by implementing a Bayesian hierarchical model to estimate the posterior distributions of the transmission rate ($\beta$) and the incubation rate ($\eta$). The result is a probabilistic forecast that quantifies uncertainty in public health outcomes.

### Core Objectives:
1.  **Deterministic Modeling:** Solving the non-linear system of Ordinary Differential Equations (ODEs) defining the SEIR transitions.
2.  **Stochastic Calibration:** Implementing the **No-U-Turn Sampler (NUTS)** to derive latent parameters from observed case counts.
3.  **Stability Analysis:** Utilizing the **Jacobian Matrix** and eigenvalue analysis to evaluate the local asymptotic stability of the disease-free equilibrium.
4.  **Counterfactual Policy Simulation:** Quantifying the impact of a 30% reduction in $\beta$ (simulating Insecticide-Treated Net interventions) on the effective reproduction number.

## 🧮 Mathematical Formulation

### 1. The SEIR System Dynamics
The population flow is governed by the following system of non-linear ODEs:

$$
\begin{aligned}
\frac{dS}{dt} &= -\frac{\beta SI}{N} \\
\frac{dE}{dt} &= \frac{\beta SI}{N} - \eta E \\
\frac{dI}{dt} &= \eta E - \gamma I \\
\frac{dR}{dt} &= \gamma I
\end{aligned}
$$

**Where:**
* $\beta$ (**Transmission Rate**): The product of the contact rate and the probability of transmission per contact.
* $\eta$ (**Incubation Rate**): The inverse of the average latent period.
* $\gamma$ (**Recovery Rate**): The inverse of the average infectious period.

### 2. Bayesian Hierarchical Likelihood
To account for overdispersion in epidemiological reporting, the observed daily incidence $Y_t$ is modeled via a **Negative Binomial distribution**:

$$Y_t \sim \text{NB}(\mu_t, \phi)$$

Where $\mu_t$ is the model-predicted incidence and $\phi$ is the concentration parameter. We assign weakly informative **Half-Normal priors** to the rates to ensure positivity and model convergence.

## 📊 Computational Implementation

### Data Pipeline
* **Solver:** `SciPy.integrate.solve_ivp` using the Runge-Kutta 45 method.
* **Inference Engine:** `PyMC` for high-dimensional MCMC sampling.
* **Diagnostics:** `ArviZ` for evaluating R-hat convergence and Effective Sample Size (ESS).

### Results and Figures
* **Trace Plots (FIG 1-4):** Confirmation of chain stationarity and posterior density for $\beta$ and $\eta$.
* **Sensitivity Analysis (FIG 5):** A bivariate analysis illustrating the relationship between recovery kinetics and transmission intensity.
* **Intervention Efficacy (FIG 6):** Comparative visualization of the "Averted Cases" metric, demonstrating the non-linear relationship between transmission reduction and total disease burden.

## 🚀 Execution
```bash
# Clone the repository
git clone [https://github.com/liobalinus-tech/SEIR_Modelling.git](https://github.com/liobalinus-tech/SEIR_Modelling.git)

# Install dependencies
pip install numpy scipy matplotlib pymc arviz

# Execute analysis
jupyter notebook SEIR_Modelling.ipynball numpy scipy matplotlib pymc arvizy:**
   ```bash
   git clone [https://github.com/](https://github.com/)[liobalinus-tech]/[SEIR_Modelling].git
