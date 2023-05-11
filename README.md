EOS
===

Equation of State from Du, Steiner, and Holt, Phys. Rev. C (2019).

Full documentation at https://neutronstars.utk.edu/code/eos .

This C++ code constructs the equation of state of homogeneous nucleonic matter for use in simulations of core-collapse supernovae and neutron star mergers.

This code was originally described in a paper on arxiv.org[https://arxiv.org/abs/1802.09710].

The source code is available on github[https://github.com/awsteiner/eos].

The installation of Boost, GSL (versions 1.16 and later), HDF5 (versions 1.8.14 and later), and the most current version of O2scl is required in order to compile the code to generate and analyze EOS tables. You do not need to compile the code to use the EOS tables - they can be read by any application which reads HDF5 files.

To simply download the source code and run it, use

```
git clone https://github.com/awsteiner/eos
```
The newest changes with compatibility of Docker is avaiable
in the v2 branch and this version has the most recent 
updates with Pions.
The dockerfile `docker` contains all the commands necessary to install all the dependencies and compile the executable `eos_nuclei`.
To do this in a docker container you need to run
```
docker build - < docker
```
This process is lengthy and will take a while to finish. To run the code, read an EOS table and see all the different outputs, one needs to run
```
docker run docker-compose.yml
```
This downloads a fiducial EOS table from https://neutronstars.utk.edu/code/eos/download.html and depending on the baryon density $n_B$, electron fraction $Y_e$ and temperature $T$ computes various nuclear properties and outputs to terminal.

One can edit the `docker-compose.yml` file to change the desired $n_B$, $Y_e$ or $T$ or change the EOS model if they so wish.

The table is an HDF5 file and the datasets are given below. Several quantities are stored in tensor_grid objects, which are stored as HDF5 groups with the name of the group given below. The contents of the rank 3 tensor are stored in the dataset named data in each group.

Grid

$n_{n_B}$: number of points in baryon density grid

$n_{Y_e}$: number of points in electron fraction grid

$n_T$: number of points in temperature grid

$n_{B_{grid}}$: array containing baryon density grid (in fm $^{-3}$)

$Y_{e_{grid}}$: array containing electron fraction grid

$T_{grid}$: array containing temperature grid (in MeV )