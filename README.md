# GasSteadySim.jl 

This package implements a nonlinear solver to determine the steady state gas flow for a pipeline network with ideal and non-ideal equations of state.
``GasSteadySim.jl`` is not a registered Julia package. Hence installation of the package should be done as follows:

```julia 
using Pkg
Pkg.add("https://github.com/kaarthiksundar/GasSteadySim.jl.git")
```

For the API usage, users are referred to the ``test/`` and the ``examples/`` directories.

## Citation
The details of the algorithm are described in the following paper:
([doi link](https://dx.doi.org/10.1109/TCNS.2022.3232524)): 

```bibtex
@article{nrsolver,
 author = {Shriram Srinivasan and Kaarthik Sundar and Vitaliy Gyrya and Anatoly Zlotnik},
 date = {2023},
 doi = {10.1109/TCNS.2022.3232524},
 journal = {IEEE Transactions on Control of Network Systems},
 number = {3},
 pages = {1449--1461},
 title = {Numerical Solution of the Steady-State Network Flow Equations for a Non-Ideal Gas},
 volume = {10}
}
```