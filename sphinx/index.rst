Du, Steiner, and Holt (2018) Phenomenological EOS
=================================================

This code constructs the equation of state of homogeneous nucleonic
matter for use in simulations of core-collapse supernovae and neutron
star mergers. This code was originally described in a paper on
`arxiv.org <https://arxiv.org/abs/1802.09710>`_. The source code is
available on `github <https://github.com/padiac/eos>`_. The
installation of `Boost <http://www.boost.org>`_, `GSL
<http://www.gnu.org/software/gsl>`_ (versions 1.16 and later), `HDF5
<http://www.hdfgroup.org>`_ (versions 1.8.14 and later), the most
current version of `O2scl <http://web.utk.edu/~asteine1/o2scl>`_ is
required.

Class eos
---------
	     
.. doxygenclass:: eos
   :members:
   :protected-members:
   :undoc-members:

Class eos_crust_virial_v2
-------------------------
	     
.. doxygenclass:: eos_crust_virial_v2
   :members:
   :protected-members:
   :undoc-members:

Class virial_solver
-------------------
	     
.. doxygenclass:: virial_solver
   :members:
   :protected-members:
   :undoc-members:

.. toctree::
   :maxdepth: 2

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
