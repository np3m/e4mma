Rest mass contribution
======================

The Helmholtz free energy of homogeneous nucleonic matter is denoted
:math:`f_{\mathrm{Hom}}(n_n^{\prime},n_p^{\prime},T)` in Du et
al. (2022) where the primes refer to the local nucleonic densities in
the gaseous or low-density phase. The rest mass energy density
corresponding to this part of the full free energy density is then
:math:`m_n n_n^{\prime} + m_p n_p^{\prime}`. The rest mass energy
density is also omitted from the free energy density of nuclei,
referred to in Du et al. (2022) as :math:`\sum_i f_i`. The rest mass
energy density associated with the nuclear contribution is then
:math:`\sum_i N_i n_i m_n + \sum_i Z_i n_i m_p`. Multiplying
:math:`f_{\mathrm{Hom}}` by :math:`\xi` and then combining these two
contributions to the free energy, we find that the total rest mass
energy density (which is not included in the published tables) is:

.. math::

   f_{\mathrm{rest}} \equiv \xi n_n^{\prime} m_n + n_p^{\prime} \xi m_p + 
   \sum_i N_i n_i m_n + \sum_i Z_i n_i m_p

then by Eq. 2 in Du et al. (2022) this is equal to

.. math::

   f_{\mathrm{rest}} = n_B (1-Y_e) m_n + n_B Y_e m_p \, .

Dividing this by :math:`n_B` gives the contribution which has
been subtracted from ``Fint`` as described in :ref:`Table Format`.
