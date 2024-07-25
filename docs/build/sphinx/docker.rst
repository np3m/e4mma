EOS with Docker
====================

This code was originally described in [Du et al.
2019](https://arxiv.org/pdf/1802.09710) and improved upon [Du et al.
2022](https://arxiv.org/pdf/2107.06697) with nuclei. The source code
is available on github. The installation of Boost, GSL (versions 1.16
and later), HDF5 (versions 1.8.14 and later), and a more current
version of O2scl (version 0.928 or later) is required in order to
compile the code to generate and analyze EOS tables. You do not need
to compile the code to use the EOS tables - they can be read by any
application which reads HDF5 files. The EoS tables are availlable to
download at our
[website](https://neutronstars.utk.edu/code/eos/download.html)

### Build the docker image
To build the EOS code inside a docker container:

clone the Github repository. Particularly the `V2` branch.

```
git clone https://github.com/awsteiner/eos && \
    cd eos && \
    git checkout v2 && git checkout 5955e74
``` 

Download the EOS table and copy it to that `eos/data/` folder. This is
done so the calculations are much faster. Since the code reads the
table and creates an output with the MUSES standard.
```
curl https://isospin.roam.utk.edu/public_data/eos_tables/du21/fid_3_5_22.o2 --output data/fid_3_5_22.o2
```
The dockerfile `Dockerfile` uses a multistage docker build to build a
minimal image with the executable `eos_nuclei` in it. The user can
build the docker image inside the eos folder themselves using
`
docker build . -t utk
`
However the process is lengthy and it is advised to just download the
already built image from dockerhub.

### Download the docker image
To download the built image from dockerhub the user can use
`
docker pull nostrad1/utk-eos:v2
`
### Running the module
After either building or downloding the image the user can just run
`docker_run_mount.sh` script locally inside the `test` folder to mount
the local input, output and data folders inside the eos folder to the
container and execute the function `utk_for_lepton` inside the
container with default configuration that creates the eos output for
lepton module in the `output` folder in `csv` format.

```
bash docker_run_mount.sh
```
This grabs the default `config.yaml` file, validates it, runs the eos
code with the validated configuration and afterwards postprocesses the
output using `muses-porter`.

Now if they want to use another configuration, they need to run the
`yaml_generator.py` in the `src` folder to create a user specific
`config.yaml` like:
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
