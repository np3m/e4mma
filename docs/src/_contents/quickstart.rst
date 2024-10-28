Quick Start Guide
====================
Docker
--------------------
Download the docker image
~~~~~~~~~~~~~~~~~~~~~~~~~~
To download the built image from dockerhub the user can use

.. code-block:: bash

    docker pull nostrad1/utk-eos:v1.9.3

Run the docker image
~~~~~~~~~~~~~~~~~~~~
To mount local directories and run the docker image as a container run

.. code-block:: bash

    docker run -it --rm --name crust-dft -u 0:0 \
    -v "${PWD}/input:/opt/e4mma/input" \
    -v "${PWD}/output:/opt/e4mma/output" \
    -v "${PWD}/data:/opt/e4mma/data" \
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

Build the docker image
~~~~~~~~~~~~~~~~~~~~~~
To build the docker image themselves, the user can use:

Clone the GitHub repository. Particularly the ``muses`` branch.

.. code-block:: bash

    git clone https://github.com/np3m/e4mma && \
    cd e4mma && \
    git checkout muses 

The dockerfile ``Dockerfile`` uses a multistage docker build to build a
minimal image with the executables ``eos`` and ` ``eos_nuclei`` in it. The user can
build the docker image inside the ``e4mma`` folder themselves using

.. code-block:: bash

    docker build . -t nostrad1/utk-eos:v1.9.3

Calculation Engine
--------------------

Running the module
~~~~~~~~~~~~~~~~~~
To run Crust-DFT module and generate an EOS table without leptons, after either building or downloading the image, the user needs a configuration file ``config.yaml`` 
in the ``input`` folder and a precomputed EOS table file in the ``data`` folder depending on the configuration.

Possible inputs for the module:

- ``output format``: format of the output files for Lepton module (either ``csv`` or ``hdf5``)
- ``verbose``: verbosity parameter for the code.(either 0,1,2)
- ``nB_grid_spec``: the function for default baryon density grid. ``'N,func(i)'``, 
                    i takes values from 0-N-1 
                    and func(i) fills up the grid . The user can change the grid length N and the 
                    desired function (default: ``'301,10^(i*0.04-12)*2.0'``)
                    ``nB_grid`` ranges from in :math:`2.0\times10^{-12}-2~\mathrm{fm^{-3}}`. Values outside this range will be ignored for now. At least ``N=3`` is needed to compute derivatives.
- ``Ye_grid_spec``: The function for default electron fraction grid. ``'N,func(i)'``, 
                    i takes values from 0-N-1 
                    and func(i) fills up the grid. The user can change the grid length N and the 
                    desired function (default: ``'70,0.01*(i+1)'``).
                    ``Ye_grid`` ranges from in :math:`1.0\times10^{-2}-0.7`. Values outside this range will be ignored for now. At least ``N=3`` is needed to compute derivatives.
- ``select_model``: EOS parameter sets from Du et al. (2022)
                    ``P_FIDUCIAL="470 738 0.5 13.0 62.4 32.8 0.9"``
                    ``P_LARGE_MMAX="783 738 0.5 13.0 62.4 32.8 0.9"``
                    ``P_SMALL_R="214 738 0.5 13.0 62.4 32.8 0.9"``
                    ``P_SMALLER_R="256 738 0.5 13.0 62.4 32.8 0.9"``
                    ``P_LARGE_R="0 738 0.5 13.0 62.4 32.8 0.9"``
                    ``P_SMALL_SL="470 738 0.5 13.0 23.7 29.5 0.9"``
                    ``P_LARGE_SL="470 738 0.5 13.0 100.0 36.0 0.9"``
                    Use one of the EOS parameterization e.g. ``'470 738 0.5 13.0 62.4 32.8 0.9'``.
- ``generate_table``: Boolean, if ``true``, generate EOS tables instead of interpolating precalculated ones.
- ``ext_guess``: Boolean, only valid when ``generate_table`` is ``true``. Use a downloaded EOS table as 
                 initial guess when producing new EOS table. This makes the process of generating new tables faster. The user can choose to not use an external guess but the time to generate new EOS table with a dense ``nB`` and ``Ye`` grid can take up a lot of time.

To generate a ``config.yaml``, run the
``yaml_generator.py`` in the ``src`` folder like:

.. code-block:: bash

    python3 yaml_generator.py \
		--generate_table false \
		--ext_guess false \
		--output_format CSV \
		--nB_grid_spec "3,0.01*(i+1)" \
		--Ye_grid_spec "3,0.4+0.01*i"

Download an EOS table and copy it to the ``data/`` folder as ``EOS_table.o2``. This is
necessary if ``generate_table`` is false or ``ext_guess`` is true. Since the code reads the
table and creates an output with the MUSES standard. The tables and their contents 
are explained in developer guide. The user is encouraged to match the model and the EOS table.

.. code-block:: bash

    curl https://isospin.roam.utk.edu/public_data/eos_tables/du21/fid_3_5_22.o2 --output data/EOS_table.o2

The user can run ``docker_run_mount.sh`` script locally inside the ``test`` folder to mount
the local input, output and data folders inside the ``e4mma`` folder to the
container and execute the function ``utk_for_lepton`` inside the
container with default configuration that creates the crust-dft output for
lepton module in the ``output`` folder in ``csv`` format.

.. code-block:: bash

    bash ../test/docker_run_mount.sh

This grabs the ``config.yaml`` file, validates it, runs the crust-dft
code with the validated configuration and afterwards post-processes the
output using ``muses-porter``.

More functions will be added later.

in the ``e4mma`` folder to get into the container. 
Creating a user specific ``config.yaml`` is similar inside the container as well. Finally, run ``run_utk_for_lepton.sh`` script inside the ``test`` folder using

.. code-block:: bash

    cd ../test \
    bash run_utk_for_lepton.sh

to validate the ``config.yaml`` generate the crust-dft output file from the user-specified configuration and post-process the file in the specified format in the output directory.
