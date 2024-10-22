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

Use Crust-DFT inside docker
~~~~~~~~~~~~~~~~~~~~~
To mount local directories and run the docker image as a container run

.. code-block:: bash

    docker run -it --rm --name crust-dft -u 0:0 \
    -v "${PWD}/input:/opt/eos/input" \
    -v "${PWD}/output:/opt/eos/output" \
    -v "${PWD}/data:/opt/eos/data" \
    nostrad1/utk-eos:v1.9.3 /bin/bash

This will put the user in a terminal inside the src folder where the 
main executable ``eos_nuclei`` is.


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
