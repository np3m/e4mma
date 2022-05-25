Table Format
============

The table is an HDF5 file and the datasets are given below. Several
quantities are stored in :ref:`o2scl:tensor_grid` objects, which are
stored as HDF5 groups with the name of the group given below. The
contents of the rank 3 tensor are stored in the dataset named data in
each group.

Grid
----

- ``n_nB``: number of points in baryon density grid
- ``n_Ye``: number of points in electron fraction grid
- ``n_T``: number of points in temperature grid
- ``nB_grid``: array containing baryon density grid
  (in :math:`\mathrm{fm}^{-3}`)
- ``Ye_grid``: array containing electron fraction grid
- ``T_grid``: array containing temperature grid
  (in :math:`\mathrm{MeV}`)

Physical quantities
-------------------

These are stored as double-precision numbers.

- ``hc``: :math:`\hbar c` (in :math:`\mathrm{MeV~fm}`)
- ``alpha_em``: :math:`\alpha_{\mathrm{EM}}`, the fine structure constant
- ``m_neut``: the mass of the neutron
  (in :math:`\mathrm{MeV}`)
- ``m_prot``: the mass of the proton in MeV
  (in :math:`\mathrm{MeV}`)

Flags
-----

These are stored as integers.

- ``baryons_only``: 1 (true) if the thermodynamic quantities with
  baryons only are provided (always 1)
- ``with_leptons``: 1 (true) if the thermodynamic quantities with
  leptons are provided, and 0 (false) otherwise.
- ``derivs_computed``: 1 (true) if the baryon part of the pressure
  and the nucleon chemical potentials are included.
- ``alg_mode``: an integer representing the algorithm used to
  compute the table
- ``include_muons``: 0 (false) if muons are not included
  
Quantities for the solver
-------------------------

These are all :ref:`o2scl:tensor_grid` objects. The contents of the
tensor are stored in an HDF5 dataset named "data" in row-major
order (the first index is the baryon density, the second is the
electron fraction, and the third is the temperature).

- ``flag``: A flag indicating if each point is complete, empty, or
  if a guess has been stored (10 means the point is complete).
- ``log_xn``: :math:`\log_{10}(x_n)` where
  :math:`x_n\equiv n_n^{\prime}/n_B`.
- ``log_xp``: :math:`\log_{10}(x_p)` where
  :math:`x_p\equiv n_p^{\prime}/n_B`.
	
- ``A``: The average nuclear mass number	
- ``Z``: The average nuclear charge number

Also, if ``alg_mode`` is either 2, 3, or 4, then the following
tensors are included:
  
- ``A_min``: the smallest value of A in the distribution, not
  including n, p, d, t, :math:`^{3}\mathrm{He}`,
  :math:`^{4}\mathrm{Li}`, and :math:`\alpha`.
- ``A_max``: largest value of A in the distribution, not
  including n, p, d, t, :math:`^{3}\mathrm{He}`,
  :math:`^{4}\mathrm{Li}`, and :math:`\alpha`.
- ``NmZ_min``: the smallest value of :math:`N-Z` in the
  distribution, not including n, p, d, t, :math:`^{3}\mathrm{He}`,
  :math:`^{4}\mathrm{Li}`, and :math:`\alpha`.
- ``NmZ_max``: the largest value of :math:`N-Z` in the
  distribution, not including n, p, d, t, :math:`^{3}\mathrm{He}`,
  :math:`^{4}\mathrm{Li}`, and :math:`\alpha`.

Composition
-----------
	
These are all :ref:`o2scl:tensor_grid` objects and included for
all tables.

- ``Xn``: the baryon number fraction of neutrons (this includes only
  the neutrons in the low-density phase, i.e. the neutrons not inside
  nuclei)
- ``Xp``: the baryon number fraction of protons (this includes only
  the protons in the low-density phase, i.e. the protons not inside
  nuclei)
- ``Xalpha``: the baryon number fraction of alpha particles
- ``Xd``: the baryon number fraction of deuterons
- ``Xt``: the baryon number fraction of tritons
- ``XHe3``: the baryon number fraction of :math:`^{3}\mathrm{He}`,
- ``XLi4``: the baryon number fraction of :math:`^{4}\mathrm{Li}`
- ``Xnuclei``: the baryon number fraction of nuclei

The sum of these eight quantities should be close to 1 for all points
in the table. Note that ``Xn`` and ``Xp`` do not store the local
density of neutrons in the gas phase, those quantities are given in Du
et al. (2022): :math:`n_n^{\prime}\equiv n_n/\xi` and
:math:`n_n^{\prime}\equiv n_p/\xi`. These densities are not stored in
the tables, but are accessible, for example, in the ``point-nuclei``
command. At high densities or temperatures when nuclei are not
present, :math:`\xi=1`.

Thermodynamic quantities
------------------------

In this section, all quantities are stored as
:ref:`o2scl:tensor_grid` objects.

Three quantities are included for all tables:

- ``Fint``: the baryon part of the free energy per baryon
  (in :math:`\mathrm{MeV}`)
- ``Sint``: the baryon part of the entropy per baryon
- ``Eint``: the baryon part of the internal energy per baryon
  (in :math:`\mathrm{MeV}`)

The relation :math:`F_{\mathrm{int}}=E_{\mathrm{int}}-T
S_{\mathrm{int}}` should hold at all points in the table.
  
If ``include_muons`` is 1, then ``Ymu``, the muon fraction,
is also included. If either ``include_muons`` or ``with_leptons``
is 1, then ``mue``, the electron chemical potential is included.
The electron chemical potential includes the electron rest mass
and is in :math:`\mathrm{MeV}`. 

If ``derivs_computed`` is 1, then the following quantities are
also included:

- ``Pint``: the baryon part of the pressure 
  (in :math:`\mathrm{MeV}/\mathrm{fm}^3`)
- ``mun``: the neutron chemical potential
  (in :math:`\mathrm{MeV}`)
- ``mup``: the proton chemical potential
  (in :math:`\mathrm{MeV}`)

The thermodynamic identity,

.. math::

   E_{\mathrm{int}} n_B = - P_{\mathrm{int}} + S_{\mathrm{int}} n_B T
   + \bar{n}_n \mu_{n} + \bar{n}_p \mu_{p}

where :math:`\bar{n}_n\equiv n_B(1-Y_e)` and :math:`\bar{n}_p\equiv
n_B Y_e` should be satisfied for all of the points in the table.

The rest mass energy :math:`m_n (1-Y_e) + m_p Y_e` has
been subtracted out from ``Fint`` and ``Eint``. The neutron and proton
rest mass have been subtracted out from the neutron and proton
chemical potentials (indepedent of whether or not the model implies a
relativistic dispersion relation for the nucleons). See also
:ref:`Chemical Potentials` for some helpful information on chemical
potentials and the associated thermodynamic identity. If
``with_leptons`` is 1, then the electron chemical potential is
included (as described above) and the following four quantities are
also included:

- ``F``: the total free energy per baryon
  (in :math:`\mathrm{MeV}`)
- ``S``: the total entropy per baryon
- ``E``: the total internal energy per baryon
  (in :math:`\mathrm{MeV}`)
- ``P``: the total pressure 
  (in :math:`\mathrm{MeV}/\mathrm{fm}^3`)

The thermodynamic identity,

.. math::

   E n_B = - P + S n_B T
   + \bar{n}_n \mu_{n} + \bar{n}_p (\mu_{p}+\mu_{e})

should be satisfied for all of the points in the table. Note that
``F`` and ``E`` do not include the baryonic part of the rest mass
energy but they *do* include the electron rest mass energy, :math:`m_e
Y_e`.
  
String arrays
-------------

For compatibility with O\ :sub:`2`\ scl, a set of two string arrays is
also included. The first, ``oth_names`` contains the list: ``Xd, Xt,
XHe3, XLi4, flag, log_xn, and log_xp``. If ``alg_mode`` is 2 or
larger, ``oth_names`` also contains ``A_min, A_max, NmZ_min,
NmZ_max``. The second, ``oth_units``, contains a set of empty strings
because none of the tensors referred to in the ``oth_names`` list have
any units. The unsigned integer ``n_oth`` contains the size of the
``oth_names`` array.

Electron and photon table
-------------------------

The electron and photon table, contains five :ref:`o2scl:tensor_grid`
objects which includes electrons, positrons, and photons, 

- ``F``: the free energy per baryon
  (in :math:`\mathrm{MeV}`)
- ``S``: the entropy per baryon
- ``E``: the internal energy per baryon
  (in :math:`\mathrm{MeV}`)
- ``P``: the pressure 
  (in :math:`\mathrm{MeV}/\mathrm{fm}^3`)
- ``mue``: the electron chemical potential
  (in :math:`\mathrm{MeV}`)

Nuclear masses table
--------------------

The nuclear massses table, contains a :ref:`o2scl:table`
object with the columns

- ``Z``: the proton number,
- ``N``: the neutron number,
- ``g``: the spin degeneracy,
- ``m``: the total mass
  (in :math:`\mathrm{MeV}`),
- ``be``: the binding energy
  (in :math:`\mathrm{MeV}`),
- ``Sn``: the neutron separation energy
  (in :math:`\mathrm{MeV}`),
- ``Sp``: the proton separation energy
  (in :math:`\mathrm{MeV}`),
- ``mass_type``: 1 for light nucleus, 2 for AME, 3 for FRDM, and
  4 for extrapolated FRDM results, and
- ``spin_type``: 1 for light nucleus, 2 for Jexp from HFB fit, 3
  for Jth from HFB fit, 4 for simple ansatz.

