# melas-heteroplasmy

This repository contains code supporting the paper ["Mitochondrial DNA density homeostasis accounts for a threshold effect in a cybrid model of a human mitochondrial disease."](https://doi.org/10.1042/BCJ20170651)

This repo is a work in progress as we curate the old code for public release. We are still verifying the full, paper-scale simulation run. If you get stuck using this project, please [raise a git issue](https://github.com/StochasticBiology/melas-heteroplasmy/issues) and we will do our best to help.

The code in `src` runs parameter inference using MCMC to estimate values and uncertainties for the parameters of a simple model relating several biological quantities related to a disease-causing mtDNA mutation. The MCMC inference is performed using custom C code, and the results are plotted using Python. The required data is included in `data`. The model itself is contained in the function `Eval_model` in `src/simulate/met_functions.h`.

Requirements: C, Python with `numpy`, `pandas`, `matplotlib`.

## Wrapper script
The script `run.sh --demo` executes a fast demo run of the pipeline. `run.sh --full` performs the full analysis for the paper, using more computer time. Raw output goes to `src/simulate`; plots in PDF and PNG format go to `src/plot`.

## Run MCMC
MCMC and model simulation is done with `met_hast.c` in `src/simulate`. This takes command-line arguments: `met_hast.ce [demo] [regulariser scale] [random seed]`. `[demo]` is 1 for the fast demo version and 0 for the full simulation; `[regulariser scale]` is 0.1 by default; `[random seed]` is 42 by default.

The `[demo]` argument controls internal parameters setting the size and reporting of the MCMC simulation.

## Plotting
Plotting is done with the `.py` scripts in `src/plot`. 

## Limitations
Please check `TODO`s wherever you see them in the code for major limitations.
