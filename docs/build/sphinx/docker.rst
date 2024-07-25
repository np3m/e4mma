EOS with Docker
====================

To run the code in a docker container use::

        docker build -t eosv2 - < dockerfilev2

This creates an ubuntu image, install all the dependencies and o2scl and compiles the 
executable ``eos_nuclei``. The process is lengthy and will take some time to finish. 
The user can change the number of processors to use during ``make`` command according 
to their machine. 
Afterwards, the user can create a docker container from the ``eosv2`` image and 
get into its terminal using::

        docker run -it --name=utk_eos eosv2 bash

To run the code, read an EOS table and see all the different outputs, one needs to 
first download an eos table from https://neutronstars.utk.edu/code/eos/download.html::

        curl https://isospin.roam.utk.edu/public_data/eos_tables/du21/fid_3_5_22.o2 
        --output data/fid_3_5_22.o2 

Currently the code reads the EOS table and depending on the baryon density :math:`n_B`, 
electron fraction :math:`Y_e` and temperature :math:`T` computes various nuclear properties.
To do this use::

        make mbnew