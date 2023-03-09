EOS
===

Equation of State from Du, Steiner, and Holt, Phys. Rev. C (2019).

Full documentation at https://neutronstars.utk.edu/code/eos .

This C++ code constructs the equation of state of homogeneous nucleonic matter for use in simulations of core-collapse supernovae and neutron star mergers.

This code was originally described in a paper on arxiv.org[https://arxiv.org/abs/1802.09710].

The source code is available on github[https://github.com/awsteiner/eos].

The installation of Boost, GSL (versions 1.16 and later), HDF5 (versions 1.8.14 and later), and the most current version of O2scl is required in order to compile the code to generate and analyze EOS tables. You do not need to compile the code to use the EOS tables - they can be read by any application which reads HDF5 files.

To run the code in a docker container you need to run
```
docker build - < docker
```
currently the code reads the EOS table and depending on the baryon density $n_B$, electron fraction $Y_e$ and temperature $T$ computes various nuclear properties.

The table is an HDF5 file and the datasets are given below. Several quantities are stored in tensor_grid objects, which are stored as HDF5 groups with the name of the group given below. The contents of the rank 3 tensor are stored in the dataset named data in each group.

Grid
$n_{n_B}$: number of points in baryon density grid

$n_{Y_e}$: number of points in electron fraction grid

$n_T$: number of points in temperature grid

$n_{B_{grid}}$: array containing baryon density grid (in fm $^{-3}$)

$Y_{e_{grid}}$: array containing electron fraction grid

$T_{grid}$: array containing temperature grid (in MeV )