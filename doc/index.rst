Equations of state for Multi-Messenger Astronomy (E4MMA)
========================================================

This C++ code constructs the equation of state
for use in simulations of core-collapse supernovae
and neutron star mergers. This code is based on `Du et al. (2019)
<https://doi.org/10.1103/PhysRevC.99.025803>`_, and `Du et al. (2022)
<https://doi.org/10.1103/PhysRevC.105.035803>`_. The source code is
available on `github <https://github.com/np3m/e4mma>`_.

You may either use one of the Docker images (see :ref:`Using docker`)
or compile the code yourself. You do not need to compile the code to
use the equation of state (EOS) tables -- they can be read by any
application which can read HDF5 files. For more information on
compiling the code, see :ref:`Compiling E4MMA`.

More documentation will be added as time permits.

.. toctree:: 
   :maxdepth: 2

   download
   table_format
   docker
   python
   examples
   compile
   rest_mass
   chem_pot
   trans
   cs2
   eos_nuclei
   class_doc
   todos
   issues

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
