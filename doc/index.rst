Equations of state for Multi-Messenger Astronomy (E4MMA)
========================================================

This C++ code constructs the equation of state of homogeneous
nucleonic matter for use in simulations of core-collapse supernovae
and neutron star mergers. This code was originally described in a
paper on `arxiv.org <https://arxiv.org/abs/1802.09710>`_. The source
code is available on `github <https://github.com/awsteiner/eos>`_. The
installation of `Boost <http://www.boost.org>`_, `GSL
<http://www.gnu.org/software/gsl>`_ (versions 1.16 and later), `HDF5
<http://www.hdfgroup.org>`_ (versions 1.8.14 and later), and the most
current version of `O2scl
<https://awsteiner.org/code/o2scl/index.html>`_ is required in
order to compile the code to generate and analyze EOS tables. You do
not need to compile the code to use the EOS tables - they can be read
by any application which reads HDF5 files.

You will need to manually edit the makefile to work with your system
and then compile ``eos_nuclei`` in order to generate an EOS. The
homogeneous matter EOS from Du et al. (2019) has a separate executable
``eos``, which can also be compiled.

More documentation will be added as time permits.

.. toctree:: 
   :maxdepth: 2

   download	      
   table_format
   examples
   docker
   eos_nuclei
   class_doc
   todos
   rest_mass
   chem_pot
   trans
   cs2
   plotting

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
