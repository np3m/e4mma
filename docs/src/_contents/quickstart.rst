Quick Start Guide
====================
Docker
--------------------
Download the docker image
~~~~~~~~~~~~~~~~~~~~~~~~~~
To download the built image from dockerhub the user can use

.. code-block:: bash

    docker pull nostrad1/utk-eos:v1.9.3

Build the docker image
~~~~~~~~~~~~~~~~~~~~~~
To build the docker image themselves, the user can use:

Clone the GitHub repository. Particularly the ``muses`` branch.

.. code-block:: bash

    git clone https://github.com/np3m/e4mma && \
    cd e4mma && \
    git checkout muses 

The dockerfile ``Dockerfile`` uses a multistage docker build to build a
minimal image with the executable ``eos_nuclei`` in it. The user can
build the docker image inside the ``e4mma`` folder themselves using

.. code-block:: bash

    docker build . -t nostrad1/utk-eos:v1.9.3

Run the docker image
~~~~~~~~~~~~~~~~~~~~
To mount local directories and run the docker image as a container run

.. code-block:: bash

    docker run -it --rm --name crust-dft -u 0:0 \
    -v "${PWD}/input:/opt/eos/input" \
    -v "${PWD}/output:/opt/eos/output" \
    -v "${PWD}/data:/opt/eos/data" \
    nostrad1/utk-eos:v1.9.3 /bin/bash

This will put the user in a terminal inside the src folder where the 
main executable ``eos_nuclei`` is. To exit the container type ``exit`` 
and press enter.

Use Crust-DFT inside docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once the user is in the container, they can run the main executable 
``eos_nuclei`` followed by ``-help`` to get accustomed to the 
command-line and the available commands

.. code-block:: bash

    ./eos_nuclei -set data_dir "../data" -help

or

.. code-block:: bash

    ./eos_nuclei -set data_dir "../data" -help point-nuclei


Evaluate the EOS (without nuclei) at a particular point.
Note that the code needs to be read a few data files, so we begin
by setting the ``data_dir`` parameter.


This is the ``fiducial`` model from Du et al. (2022)

.. code-block:: bash

    ./eos -set data_dir "../data" \
    -select-model 470 738 0.5 13.0 62.4 32.8 0.9 -point 0.08 0.5 5.0

Pure Skyrme model

.. code-block:: bash

    ./eos -set data_dir "../data" -alt-model Skyrme NRAPR \
    -point 0.08 0.5 5.0

RMF support is still experimental

.. code-block:: bash

    ./eos -set data_dir "../data" -alt-model RMF SFHo \
    -point 0.08 0.5 5.0

For non-RMF models, the code without nuclei also works at T=0

.. code-block:: bash

    ./eos -set data_dir "../data" \
    -select-model 470 738 0.5 13.0 62.4 32.8 0.9 -point 0.08 0.5 5.0 -point 0.08 0.5 0.0
    ./eos -set data_dir "../data" -alt-model Skyrme NRAPR \
    -point 0.08 0.5 0.0


Create a small table with derivatives based on an initial guess


Download the initial guess. The file is compared with the SHA256
hash and only downloaded if the current file doesn't match the hash.
The `acol` command is part of O2scl (one of the e4mma dependencies).
Instead of acol, you can just use, e.g. 'curl' to download the file
and `openssl dgst -sha256` to obtain the hash.

.. code-block:: bash

    acol -download ../output/fid_3_5_22.o2 \
    https://isospin.roam.utk.edu/public_data/eos_tables/du21/fid_3_5_22.o2 \
    840f6f171f05081deed53fd8bf50bad1b16a865418c37b1b630817ae10ad6736

Select a random EOS parameterization, create the table, and then
compute derivatives and store it in E_table_deriv.o2. This table
does not include leptons.

.. code-block:: bash

    ./eos_nuclei -set data_dir "../data" -random \
    -set nB_grid_spec "5,0.01*(i+1)" -set Ye_grid_spec "3,0.4+0.01*i" \
    -set T_grid_spec "3,5+i" -generate-table \
    "ext_guess=../data/fid_3_5_22.o2" -eos-deriv \
    -output ../output/E_table_deriv.o2

For more examples see the script files in examples directory.

Calculation Engine
--------------------

Running the module
~~~~~~~~~~~~~~~~~~
After either building or downloading the image, the user needs a configuration file ``config.yaml`` 
in the ``input`` folder and an EOS table file in the ``data`` folder.

To generate ``config.yaml``, run the
``yaml_generator.py`` in the ``src`` folder like:

.. code-block:: bash

    cd ../src
    python3 yaml_generator.py \
	    --output_format HDF5 \
        --nB_grid_spec '150,10^(i*0.04-12)*2.0' \
	    --Ye_grid_spec '30,0.01*(i+1)' \
        --inc_lepton false

Download an EOS table and copy it to the ``data/`` folder as ``EOS_table.o2``. This is
done, so the calculations are much faster. Since the code reads the
table and creates an output with the MUSES standard. The tables and their contents 
are explained in developer guide.

.. code-block:: bash

    curl https://isospin.roam.utk.edu/public_data/eos_tables/du21/fid_3_5_22.o2 --output data/EOS_table.o2

The user can run ``docker_run_mount.sh`` script locally inside the ``test`` folder to mount
the local input, output and data folders inside the ``e4mma`` folder to the
container and execute the function ``utk_for_lepton`` inside the
container with default configuration that creates the crust-dft output for
lepton module in the ``output`` folder in ``csv`` format.

.. code-block:: bash

    bash docker_run_mount.sh

This grabs the ``config.yaml`` file, validates it, runs the crust-dft
code with the validated configuration and afterwards post-processes the
output using ``muses-porter``.

Possible inputs for the module:

- ``output format``: format of the output files for Lepton module (either ``csv`` or ``hdf5``)
- ``verbose``: verbosity parameter for the code.(either 0,1,2)
- ``nB_grid_spec``: the function for default baryon density grid. ``'N,func(i)'``, 
                    i takes values from 0-N-1 
                    and func(i) fills up the grid . The user can change the grid length N and the 
                    desired function (default: ``'301,10^(i*0.04-12)*2.0'``)
                    ``nB_grid`` ranges from in :math:`2.0\times10^{-12}-2~\mathrm{fm^{-3}}`. Values outside this range will be ignored for now
- ``Ye_grid_spec``: The function for default electron fraction grid. ``'N,func(i)'``, 
                    i takes values from 0-N-1 
                    and func(i) fills up the grid. The user can change the grid length N and the 
                    desired function (default: ``'70,0.01*(i+1)'``).
                    ``Ye_grid`` ranges from in :math:`1.0\times10^{-2}-0.7`. Values outside this range will be ignored for now

- ``inc_lepton``: whether to include leptons or not (boolean, default: ``False``)
More functions will be added later.



in the ``e4mma`` folder to get into the container. 
Creating a user specific ``config.yaml`` is similar inside the container as well. Finally, run ``run_utk_for_lepton.sh`` script inside the ``test`` folder using

.. code-block:: bash

    bash run_utk_for_lepton.sh

to validate the ``config.yaml`` generate the crust-dft output file from the user-specified configuration and post-process the file in the specified format in the output directory.
