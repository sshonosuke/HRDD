# Hierarchical regression disconinuity design

This repository provides R code implementing Hierarchical regression disconinuity design (HRDD), as proposed by the following paper.

Sugasawa, S., Ishihara, T. and Kurisu, D. (2023). Hierarchical regression disconinuity design: pursuing subgroup treatment effects. *arXiv*.

The repository includes the following files.

- `HRDD-Continuous.R` : Implementation of HRDD under continuous response (quadratic loss)  
- `HRDD-Binary.R` : Implementation of HRDD under binary response (negative logistic log-likelihood loss)  
- `Sim-oneshot-continuous.R`: One-shot simulation study with continuous response
- `Sim-oneshot-binary.R`: One-shot simulation study with binary response
- `Dataset.RData`: Columbia scholarship data 
- `Application.R`: Example of applying HRDD to Columbia scholarship data 
