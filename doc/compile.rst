Compiling E4MMA
===============

Compiling the code requires the installation of `Boost
<http://www.boost.org>`_, `GSL <http://www.gnu.org/software/gsl>`_
(versions 2.0 and later), `HDF5 <http://www.hdfgroup.org>`_, and
`O2scl <https://awsteiner.org/code/o2scl/index.html>`_. You will need
to manually edit the makefile to work with your system and then
compile ``eos_nuclei`` in order to generate an EOS.

Alpha versions of E4MMA are currently built around specific versions
of O2scl/O2sclpy and commits from the ``andrew`` branch of E4MMA.

E4MMA Alpha 5
-------------

- O2scl/py v0.930a4 (based on docker file `v0.930a4_u24.04_py
  <https://github.com/awsteiner/o2scl/blob/dev/docker/v0.930a4_u24.04_py>`_)
  and E4MMA commit `b4ea4a7c
  <https://github.com/np3m/e4mma/commit/b4ea4a7c3d3abca43dee51d3763f14c4ff465796>`_.
  See E4MMA docker file `alpha5_ju_o0930a4_u24.04 <https://github.com/np3m/e4mma/blob/andrew/docker/alpha5_ju_o930a4_u24.04>`_.
