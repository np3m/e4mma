EOS
===

Equation of State from Du, Steiner, and Holt, Phys. Rev. C (2019).

Full documentation at https://neutronstars.utk.edu/code/eos .

To run the EOS code in a docker container:

1.  clone the Github repository. Particularly the ```V2``` branch.

```
git clone https://github.com/awsteiner/eos && \
    cd eos && \
    git switch v2
```  

2. Download the EOS table and copy it to that ```eos/data/``` folder.
```
curl https://isospin.roam.utk.edu/public_data/eos_tables/du21/fid_3_5_22.o2 --output data/fid_3_5_22.o2
```
3. The dockerfile ```dockerfilev3``` uses a multistage docker build to build a minimal image with the executable ```eos_nuclei``` in it.
The user can build the docker image inside the eos folder themselves using
```
docker build . -t utk -f dockerfilev3
```
However the process is lengthy and it is advised to just download the already built image from dockerhub.
```
docker pull nostrad1/utk-eos:v2
```
4. After step 3 the user can just run ```docker_run_mount.sh``` script locally to mount the local input, output and data
folders inside the eos folder to the container and execute the function ```utk_for_lepton``` inside the container with default configuration that creates the eos output 
for lepton module in the ```api/output``` folder in ```csv``` format.

```
bash docker_run_mount.sh
```
This grabs the default config.yaml file, validates it, runs the eos code with the validated configuration and afterwards postprocesses the output using muses-porter.

Now if they want to use another configuration, they need to run the ```yaml_generator.py``` to create a user specific ```config.yaml``` like:
```
python3 yaml_generator.py --select_model "470 738 0.5 13.0 62.4 32.8 0.9" \
	--a_virial 10.0 --b_virial 10.0 \
	--load data/fid_3_5_22.o2 \
	--output_format HDF5
```
before running the previous command.

5. If the user wants to get into the container and run the code from inside, use
```
docker run -it --rm --name utk -u 0:0 \
  -v "${PWD}/api/input:/opt/eos/api/input" \
  -v "${PWD}/api/output:/opt/eos/api/output" \
  -v "${PWD}/data:/opt/eos/data" \
  nostrad1/utk-eos:v2 /bin/bash
```
to get into the container. 
Creating a user specific config.yaml is similar inside the container as well. Finally run ```run_utk_for_lepton.sh``` script using
```
bash run_utk_for_lepton.sh
```
to validate the ```config.yaml``` generate the eos output file from the user-specified configuration and post-process the file in the specified format in the output directory.