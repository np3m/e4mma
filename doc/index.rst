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
application which can read HDF5 files.

Compiling the code requires the installation of `Boost
<http://www.boost.org>`_, `GSL <http://www.gnu.org/software/gsl>`_
(versions 2.0 and later), `HDF5 <http://www.hdfgroup.org>`_ (versions
1.8.14 and later), and `O2scl
<https://awsteiner.org/code/o2scl/index.html>`_. You will need to
manually edit the makefile to work with your system and then compile
``eos_nuclei`` in order to generate an EOS.

More documentation will be added as time permits.

.. toctree:: 
   :maxdepth: 2

   download	      
   table_format
   docker
   python
   examples
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
