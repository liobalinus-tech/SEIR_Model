# Bayesian SEIR Modeling of Malaria Outbreaks in Northern Namibia
This repository contains a modular Python framework for modeling malaria transmission dynamics using Bayesian Inference. The project captures seasonal fluctuations in infection rates and evaluates the statistical impact of health interventions.

### Project Overview
The model utilizes a compartmental SEIR (Susceptible-Exposed-Infectious-Recovered) structure to simulate disease spread. By leveraging No-U-Turn Sampling (NUTS), the system calibrates epidemiological parameters against observed clinical data, providing high-fidelity uncertainty quantification.

### Technical Architecture
The codebase is decoupled into three primary modules for scalability and testing:
* seir_engine.py: The mathematical core, containing the ODE systems and the PyMC Bayesian model configuration.
* Visualizer.py: A comprehensive visualization suite that generates clinical figures, including uncertainty intervals and phase portraits.
* main.py: The execution controller that orchestrates the inference and generates results.
### Key Results
The model generates six critical clinical visuals (located in the /figures directory):
* Mass Action Balance: State transitions through the SEIR compartments.
* Uncertainty Quantification: Model calibration with 94% Highest Density Intervals (HDI).
* Reproduction Number Dynamics: Seasonal oscillations of $R(t)$.
* Phase Portrait: Geometric stability of the $S-I$ relationship.
* Strategic Sensitivity Analysis: Heatmap identifying peak clinical loads under varying transmission/recovery rates.
* Policy Simulation: Quantitative impact of intervention strategies, such as bed-net distribution.
### Getting Started
**Prerequisites**
Ensure you have the following libraries installed:\
*Bash*\
pip install numpy pymc arviz matplotlib seaborn scipy pytensor\
**Execution**\
To run the full simulation and generate all figures, execute:\
*Bash*\
python main.py\
**Research Documentation**\
For a deep dive into the mathematical proofs, clinical parameter selection, and socio-economic analysis, please refer to the Mini-Thesis PDF included in this repository.
