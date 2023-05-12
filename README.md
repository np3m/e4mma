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
The dockerfile appropriately named `Dockerfile` contains all the commands necessary to create an ubuntu image, install all the dependencies and compile the executable `eos_nuclei`. To do this you need to run
```
docker build -t eosv2 - < Dockerfile
```
This process is lengthy and will take a while to finish. Afterwards, the user can create a docker container from the`eosv2` image and get into its terminal using

```
 docker run -it --name=utk_eos eosv2 bash
```
To run the code, read an EOS table and see all the different outputs, one needs to first download an eos table from https://neutronstars.utk.edu/code/eos/download.html
```
curl https://isospin.roam.utk.edu/public_data/eos_tables/du21/fid_3_5_22.o2 --output data/fid_3_5_22.o2; 
```
and use somethng like
```
./eos_nuclei 
        -set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load data/fid_3_5_22.o2 \
		-set recompute 1 \
        -point-nuclei 0.1 0.4 30
```
Using this fiducial EOS table and depending on the baryon density $n_B$, electron fraction $Y_e$ and temperature $T$, the code computes various nuclear properties and outputs to terminal.

One can edit the arguments of the `point-nuclei` function for now to change the $n_B$, $Y_e$ or $T$ respectively. If they wish to use a different EOS model, they need to download a different table from https://neutronstars.utk.edu/code/eos/download.html.

The table is an HDF5 file and the datasets are given below. Several quantities are stored in tensor_grid objects, which are stored as HDF5 groups with the name of the group given below. The contents of the rank 3 tensor are stored in the dataset named data in each group.

Grid

$n_{n_B}$: number of points in baryon density grid

$n_{Y_e}$: number of points in electron fraction grid

$n_T$: number of points in temperature grid

$n_{B_{grid}}$: array containing baryon density grid (in fm $^{-3}$)

$Y_{e_{grid}}$: array containing electron fraction grid

$T_{grid}$: array containing temperature grid (in MeV )