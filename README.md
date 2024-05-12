UTK-EoS
===

Equation of State from Du, Steiner, and Holt, Phys. Rev. C (2019).

Full documentation at https://neutronstars.utk.edu/code/eos .

## Introduction

Our EoS describes the relationship between energy density and pressure in a large three-dimensional ($\mathrm{nB}$, $\mathrm{Ye}$, $\mathrm{T}$) space,

There are several different physical regimes each constrained by different observables and theoretical approaches. 

First the zero temperature nuclear
matter at nuclear saturation density is closely connected
to nuclear masses, charge radii, giant resonances, and
other laboratory observables. Global fits to experimental data have been performed with Skyrme and covariant mean-field models. 

Secondly, cold neutron matter below nuclear saturation density, is difficult to probe experimentally but is well-constrained by theoretical calculations based on
semi-phenomenological nuclear forces or microscopic chiral effective field theory-based interactions. 

The third regime, strongly-interacting high-temperature matter, is best described by interactions and many-body approaches similar to those applied to cold neutron matter near saturation density. 

The fourth regime, low-density and high-temperature matter that is nearly
non-degenerate, is best described by the virial expansion.
The equation of state in this regime is determined from
nucleon-nucleon scattering phase shifts. 

Finally, neutron-rich matter at densities above twice saturation
density is most strongly constrained by observations of
neutron star masses and radii, particularly the observation of neutron stars with $M \simeq 2M_\odot$.

In this work, we construct a phenomenological free energy density that is consistent with observational and
theoretical constraints in the five aforementioned physical regimes. This is in contrast to works which attempt to
describe matter over the entire density and temperature
range with a single detailed model of the nucleon-nucleon
interaction.

Our second advance is in the treatment of uncertainties. The most relevant parameters which describe
the uncertainties in different density and temperature
regimes are not clearly related. In this work, through the
construction of a phenomenological model one can vary
uncertainties in different regimes independently, without
spoiling agreement elsewhere.

The module produces an EoS table in a 3 dimensional grid of $\mathrm{nB}$, $\mathrm{Ye}$ and $\mathrm{T}$ and stores various physical and thermodynamics properties in the data files (detailed documentation [here](https://neutronstars.utk.edu/code/eos/table_format.html#thermodynamic-quantities)).

## Physics
First, we define the symmetry energy to include a zero
temperature contribution which combines the QMC EOS
 near saturation density, the neutron star fit at higher den
sities, and the Skyrme interaction for isospin-symmetric
 matter
$$\begin{equation}
\epsilon_{sym}(n_B) = h(n_B)\epsilon_{QMC}(nB) + [1-h(n_B)]\epsilon_{NS}(n_B) - f_{Skyrme}(nB,x_p = 1/2, T=0)
\end{equation}$$
Defining the isospin asymmetry $ \delta = 1-2x_p$, we can
 combine this with the model described in [Du et al](https://arxiv.org/pdf/1802.09710) to obtain
 the free energy density of degenerate matter

$$\begin{equation}
f_{deg}(n_B,x_p,T) = f_{Skyrme}(nB,x_p = 1/2, T=0) + \delta^2\epsilon_{sym}(n_B) + \delta^2\Delta f_{hot}(nB,x_p = 0, T) + (1-\delta^2)\Delta f_{hot}(nB,x_p = 1/2, T)
\end{equation}$$

Finally,we ensure that the total nucleonic free energy
 gives the result from the virial expansion at high tem
peratures using 

$$\begin{equation}
f_{np}(n_B,x_p,T) = f_{virial}(n_n,x_p,T)g+f_{deg}(n_B,x_p,T)(1-g)
\end{equation}$$

When we need to include the
 electrons, positrons, and photons, we define the free en
ergy density

$$\begin{equation}
f_{npe\gamma} \equiv f_{np} + f_{e^-}+f_{e^+}+f_\gamma
\end{equation}$$

Using this formalism, the chemical potentials and entropy can be computed directly (eq. 28-32 in [Du et al](https://arxiv.org/pdf/1802.09710)).

We enforce causality at high densities.

## Docker
This code was originally described in [Du et al](https://arxiv.org/pdf/1802.09710). The source code is available on github. The installation of Boost, GSL (versions 1.16 and later), HDF5 (versions 1.8.14 and later), and a more current version of O2scl (version 0.928 or later) is required in order to compile the code to generate and analyze EOS tables. You do not need to compile the code to use the EOS tables - they can be read by any application which reads HDF5 files. The EoS tables are availlable to download at our [website](https://neutronstars.utk.edu/code/eos/download.html)

### Build the docker image
To build the EOS code inside a docker container:

clone the Github repository. Particularly the `V2` branch.

```
git clone https://github.com/awsteiner/eos && \
    cd eos && \
    git checkout v2 && git checkout 93ca543
``` 

Download the EOS table and copy it to that `eos/data/` folder. This is done so the calculations are much faster. Since the code reads the table and creates an output with the MUSES standard.
```
curl https://isospin.roam.utk.edu/public_data/eos_tables/du21/fid_3_5_22.o2 --output data/fid_3_5_22.o2
```
The dockerfile `Dockerfile` uses a multistage docker build to build a minimal image with the executable `eos_nuclei` in it.
The user can build the docker image inside the eos folder themselves using
`
docker build . -t utk
`
However the process is lengthy and it is advised to just download the already built image from dockerhub.

### Download the docker image
To download the built image from dockerhub the user can use
`
docker pull nostrad1/utk-eos:v2
`
### Running the module
After either building or downloding the image the user can just run `docker_run_mount.sh` script locally inside the `test` folder to mount the local input, output and data
folders inside the eos folder to the container and execute the function `utk_for_lepton` inside the container with default configuration that creates the eos output 
for lepton module in the `output` folder in `csv` format.

```
bash docker_run_mount.sh
```
This grabs the default `config.yaml` file, validates it, runs the eos code with the validated configuration and afterwards postprocesses the output using `muses-porter`.

Now if they want to use another configuration, they need to run the `yaml_generator.py` in the `src` folder to create a user specific `config.yaml` like:
```
cd ../src
python3 yaml_generator.py --load data/fid_3_5_22.o2 \
	--output_format HDF5 \
    --nB_grid_spec '150,10^(i*0.04-12)*2.0' \
	--Ye_grid_spec '30,0.01*(i+1)' \
```
before running the previous command.

### Possible inputs for the module:
1. `load` : Loads the EoS table. Currently only `"../data/fid_3_5_22.o2"` is supported.
2. `output format`: Format of the output files for Lepton module (either `csv` or `hdf5`).
3. `verbose`: verbosity parameter for the code.(either 0,1,2)
4. `nB_grid_spec`: The function for default baryon density grid. `'N,func(i)'`, 
                    i takes values from 0-N 
                    and func(i) fills up the grid . The user can change the grid length N and the 
                    desired function (default: `'301,10^(i*0.04-12)*2.0'`  ).
                    `nB_grid` ranges from $2.0\times10^{-12}-2~\mathrm{fm^{-3}}$. Values outside this range will be ignored for now.
5. `Ye_grid_spec`: The function for default electron fraction grid. `'N,func(i)'`, 
                    i takes values from 0-N 
                    and func(i) fills up the grid. The user can change the grid length N and the 
                    desired function (default: `'70,0.01*(i+1)'`).
                    `Ye_grid` ranges from $1.0\times10^{-2}-0.7$. Values outside this range will be ignored for now.

More functions will be added later.

### Use EoS inside docker
If the user wants to get into the container and run the code from inside, use
```
docker run -it --rm --name utk -u 0:0 \
  -v "${PWD}/input:/opt/eos/input" \
  -v "${PWD}/output:/opt/eos/output" \
  -v "${PWD}/data:/opt/eos/data" \
  nostrad1/utk-eos:v2 /bin/bash
```
 in the `eos` folder to get into the container. 
Creating a user specific config.yaml is similar inside the container as well. Finally run `run_utk_for_lepton.sh` script inside the `test` folder using
```
bash run_utk_for_lepton.sh
```
to validate the `config.yaml` generate the eos output file from the user-specified configuration and post-process the file in the specified format in the output directory.

### Experimental
To compute EoS at a certain $\mathrm{nB~(fm^{-3})}$ , $\mathrm{Ye}$ and $\mathrm{T~(MeV)}$, Generate a configuration, using `point_generator.py` in the `src` folder to create a user specific `point.yaml` like:
```
python3 point_generator.py \
	--select_model "470 738 0.5 13.0 62.4 32.8 0.9" \
	--a_virial 10 --b_virial 10 \
	--load ../data/fid_3_5_22.o2 \
	--point_nuclei "0.16 0.4 30" 
```

#### Possible inputs
##### Options:
Select as many as you like.

`select_model`: Select an EOS model. The possible inputs are 
                P_FIDUCIAL=`"470 738 0.5 13.0 62.4 32.8 0.9"`
                P_LARGE_MMAX=`"783 738 0.5 13.0 62.4 32.8 0.9"`
                P_SMALL_R=`"214 738 0.5 13.0 62.4 32.8 0.9"`
                P_SMALLER_R=`"256 738 0.5 13.0 62.4 32.8 0.9"`
                P_LARGE_R=`"0 738 0.5 13.0 62.4 32.8 0.9"`
                P_SMALL_SL=`"470 738 0.5 13.0 23.7 29.5 0.9"`
                P_LARGE_SL=`"470 738 0.5 13.0 100.0 36.0 0.9"`
                default is `"470 738 0.5 13.0 62.4 32.8 0.9"`
                Select only one

`load`: Loads an EOS table in to memory. default: `"../data/fid_3_5_22.o2"`

`select_high_T`: Select 0 for the original DSH fit, 1 for NRAPR, 2 for Sk chi 414, 3 for
                Skchi450, 4 for Skchi500, 5 for ?, "+ and 6 for Sk chi m* (the default).
              default: 6

`eos_deriv`: Compute derivatives numerically.
              default: `'0'`


`a_virial`: Coefficient for modulation of virial EOS. default: `10.0`

`b_virial`: Coefficient for modulation of virial EOS. default: `10.0`

`include_muons`: If true, include muons. default: `false`

`max_ratio`: The maximum value of N/Z or Z/N. default: `7.0`

`mh_tol_rel`: Relative tolerance for the solver in the `eos_fixed_dist()` function. default: `1.0e-06`

`recompute`: If true, recompute all points, irrespective of the value of the convergence flag. default: `false`

#### commands:
Select only one

`get`: This command gets the value of a parameter. default: `'a_virial' `
        
`point`: Evaluate the EOS at one (nB,Ye,T) point.
              default: `"0.16 0.4 0.1"`

`point_nuclei`: Compute and/or show EOS results at one `(n_B,Y_e,T)` point. default: `"0.16 0.4 0.1"`
              
`random`: Select a random EOS, checking several physical constraints and re-selecting a
                new random EOS until all the constraints are met.
              default: `'0'`

Run `test_conf.sh` script inside the `test` folder using
```
bash test_conf.sh
```
To see a terminal output of the computed quantities. However functionality of this script is currently limited and more options will be added later. To use the eos fully, use makefile and CLI which are avilable.
