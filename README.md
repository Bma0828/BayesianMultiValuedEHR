# BayesianMultiValuedEHR: Bayesian Feature Selection for Multi-valued Treatment Comparisons

This repository contains the implementation of a Bayesian feature selection framework for analyzing multi-valued treatment comparisons in electronic health records, focusing on vasopressor effectiveness.

## Overview

This code implements a novel framework combining instrumental variable (IV) analysis with Bayesian feature selection methods and neural networks to estimate causal effects in multi-valued treatment settings. The methodology addresses three key challenges:
- Handling multiple treatment comparisons simultaneously
- Comparing Bayesian feature selection methods
- Selecting relevant features while capturing complex nonlinear relationships

## Features

- Instrumental Variable (IV) analysis using physician prescribing preferences
- Implementation of multiple feature selection methods:
  - Spike-and-Slab priors
  - Bayesian LASSO
  - Standard LASSO
- Neural network models for outcome prediction
- Comprehensive causal effect estimation

## Application

The framework was applied to compare three commonly used vasopressors (norepinephrine, vasopressin, and phenylephrine) using the MIMIC-IV database, revealing a clear hierarchical pattern in treatment effectiveness.

## Citation

If you use this code in your research, please cite our paper:

```
Qian, Yunzhe, and Bowen Ma. "Bayesian Feature Selection for Multi-valued Treatment Comparisons: 
An Electronic Health Records Study of Vasopressor Effectiveness." medRxiv (2024): 2024-12.
```

## Contributors

- Yunzhe Qian - Harvard T.H. Chan School of Public Health
- Bowen Ma - Harvard T.H. Chan School of Public Health
