Crust-DFT EOSs for Astrophysical Simulations
======================================

This c++ code was originally described in a
paper on `arxiv.org <https://arxiv.org/abs/1802.09710>`_. The source
code is available on `github <https://github.com/np3m/e4mma>`_. The
installation of `Boost <http://www.boost.org>`_, `GSL
<http://www.gnu.org/software/gsl>`_ (versions 1.16 and later), `HDF5
<http://www.hdfgroup.org>`_ (versions 1.8.14 and later), and at least 
version `eee9fd83` of `O2scl <https://neutronstars.utk.edu/code/o2scl/index.html>`_ 
is required in order to compile the code to generate and analyze EOS tables. 

You do not need to compile the code to use the EOS tables - they can be read
by any application which reads HDF5 files.

More documentation will be added as time permits.

.. toctree:: 
   :maxdepth: 3

   quickstart
   developer
   physics

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
