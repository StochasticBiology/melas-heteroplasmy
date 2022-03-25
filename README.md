# melas-heteroplasmy

This repository contains code supporting the paper ["Mitochondrial DNA density homeostasis accounts for a threshold effect in a cybrid model of a human mitochondrial disease."](https://doi.org/10.1042/BCJ20170651)

This project was carried out roughly between 2014-16 (back before the first author knew what version control was :) ), and its first upload to GitHub was made in 2022. Whilst we have made a modest effort to curate this project for ease of use, the amount of time that has passed and the complexity of the underlying project makes it difficult for us to guarantee that all of the code here is in working order.

If you get stuck using this project, please [raise a git issue](https://github.com/StochasticBiology/melas-heteroplasmy/issues) and we will do our best to help.

# Run MCMC

To run the MCMC:
1. `cd src/simulate`
2. Set the parameters `NUM_STEPS`, `REPORT` and `REPORT_BUFFER_LEN` appropriately. The code is currently set up to complete quickly in a demo mode, but the comments give values used in the paper.
3. Run `sh run.sh` in any linux environment with gcc set up

# Plotting

Install the conda environment
```bash
conda env create -f environment.yml
```
then run the code under `src/plot`. You may need to have LaTeX installed to make the plots, but if you can live without then just remove the lines
```python 
plt.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
```
from plotting scripts wherever you see them.

# Limitations
Please check `TODO`s wherever you see them in the code for major limitations.
