# LES of Channnel flow with OpenFOAM

## Introduction

This repository contains OpenFOAM v2306 case setup files and results for channel
flow at $Re_\tau = 1000$. 

In the data M1, M2, and M3 correspond to 3 different grids, increasingly fine.
The M2_rk and M2_sigma cases correspond to simulations on the M2 grid using the
[RKSymFOAM](https://github.com/janneshopman/RKSymFoam) solver, and the Sigma LES
model.

## Software

All the results were generated using OpenFOAM v2306. Additional software
includes:

- [swak4Foam](https://sourceforge.net/p/openfoam-extend/swak4Foam) A large
toolkit for OpenFOAM. In this project only `funkySetFields` is used.
- [RKSymFOAM](https://github.com/janneshopman/RKSymFoam) Low-dissipation solver
  with Runge-Kutta time steppers.
- [runTimeChannelBudgets](https://github.com/janneshopman/runTimeChannelBudgets)
  Tool for computing TKE budget terms.
- [ofreaders](https://gitlab.com/chalmers-marine-technology/ofreaders/) Tool for
handling various OpenFOAM-generated files in Python. Here used for converting
probe data to HDF5.
- Standard Python packages for scientific computing: `numpy`, `matplotlib`,
`h5py`.

## Directory structure

- The folder `scripts` contains Python scripts that can be used to generate
figures and compute quantities like spectra.
- The folder `dns` contains a bash script for downloading the DNS data by Lee
and Moser (2015).


## Data for each case

For each of the cases the following is provided:

- A dummy `0.org` directory with constant initial conditions.
- A `constant` directory, with no mesh.
- A `system` directory with all the necessary configuration files, including a
  `blockMeshDict` to generate the grid.
- A `postProcessing` directory with subdirectory `collapsedFields` containing
  the statistics at time 750 and (for selected cases) subdirectory
  `averageYPlus` containing the development of spatially-averaged $y^+$ in time.
- A `graphs/750` directory containing the profiles of the TKE budget terms.  
- A `plane_stats.txt` file containing single-point statistics computed using
  probe data. These can be compared to those in the `postProcessing` directory.
- `z_spectrum_*.txt` files containing the spanwise spectra for each wall-normal
  value.
- `coherence_*.txt` files containing the linear coherence spectra with respect
  to $y^+ = 2$.


## Probe data

Due to the size of the produced probe data files (> 100 Gb for the M3 case) we
are unable to easily share them at scale.
The data can be considered available upon request.
However, we provide a portion of the data for the M1 grid via the following
links.
In these files only every 4th time-step is retained.

[streamwise component](https://u.pcloud.link/publink/show?code=XZg7AF5ZVzj67WwapIu4kM4GRPwY859RflfX)

[wall-normal component](https://u.pcloud.link/publink/show?code=XZqJAF5ZmctdvgfsfdX5ny99q4RFyFnesxTk)

[spanwise component](https://u.pcloud.link/publink/show?code=XZAJAF5Zrkd0e60GXdY2DOpMsrppn7UfnLh7)



