Chemical potentials
===================

Here, we explain the ambiguities associated with the definition of
the chemical potentials and how they relate to the code.

Denote the number density of neutrons in the low-density phase
:math:`n_n`, the number density of protons in the low-density phase
:math:`n_p`, the number density of electrons :math:`n_e`, and the
number density of nuclei as :math:`\{n_i\}`. Denote the total number
of neutrons and protons in both phases as :math:`\bar{n}_n` and
:math:`\bar{n}_p`. Using these definitions, one can write the free
energy for hot and dense matter in (at least) four different ways,
:math:`f_1(\bar{n}_n,\bar{n}_p,T)`,
:math:`f_2(\bar{n}_n,\bar{n}_p,n_e,T)`,
:math:`f_3(n_n,n_p,\{n_i\},T)`, :math:`f_4(n_n,n_p,\{n_i\},n_e,T)`. In
the first form, the Saha equations have been solved to determine
:math:`\{n_i\}` and charge neutrality has been used to determine
:math:`n_e`. In the second form, the Saha equations have been solved
but charge neutrality has not been used. In the third form, the Saha
equations have not been solved but charge neutrality has been used.
The electron contribution to the free energy is included in all four
free energies, but in the case of :math:`f_1` and :math:`f_3`, the
electron density is not independent of the other densities. For these
four free energies, there are four corresponding proton chemical
potentials, :math:`\partial f_1/\partial \bar{n}_p`, :math:`\partial
f_2/\partial \bar{n}_p`, :math:`\partial f_3/\partial n_p`, and
:math:`\partial f_4/\partial n_p`. *None of these four proton chemical
potentials are the same.* This documentation attempts to explain how
this complication relates to the code. In Du et al. (2022), we use a
confusing notation because we do not clearly distinguish
:math:`\bar{n}_n` and :math:`n_n`. The function :math:`f_1` is most
directly related to the tables which are generated and one can simply
identify :math:`\bar{n}_n=n_B(1-Y_e)` and :math:`\bar{n}_p=n_B Y_e`.
      
The comparison between :math:`f_1` and :math:`f_2` is the simplest
(now being a bit more careful about what is held constant)

.. math::

   \left(\frac{\partial f_1}{\partial \bar{n}_n}\right)_{\bar{n}_p,T} =
   \left(\frac{\partial f_2}{\partial \bar{n}_n}\right)_{\bar{n}_p,n_e,T}
   \quad \mathrm{and} \quad
   \left(\frac{\partial f_1}{\partial \bar{n}_p}\right)_{\bar{n}_n,T} =
   \left(\frac{\partial f_2}{\partial
   \bar{n}_p}\right)_{\bar{n}_n,n_e,T} +
   \left(\frac{\partial f_2}{\partial
   n_e}\right)_{\bar{n}_n,\bar{n}_p,T}

To simplify the discussion we use the following notation:

.. math::

   \mu_{p,i} \equiv \left( \frac{\partial f_i}{\partial \bar{n}_p}
   \right)

where all of the other densities are held constant, including either
:math:`n_n` or :math:`\bar{n}_n` as appropriate. Thus :math:`f_1` and
:math:`f_2` imply two thermodynamic identies

.. math::

   \varepsilon_1 &=& - P_1 + T s_1 + \bar{n}_n \mu_{n,1} +
   \bar{n}_p \mu_{p,1} \nonumber \\
   \varepsilon_2 &=& - P_2 + T s_2 + \bar{n}_n \mu_{n,2} +
   \bar{n}_p \mu_{p,2} + n_e \mu_e

When :math:`n_e=\bar{n}_p`, we have :math:`P_1=P_2`,
:math:`\varepsilon_1=\varepsilon_2`, and :math:`s_1=s_2`. In the EOS
literature, it has become standard to store :math:`\mu_{n,2}` and
refer to it as the "neutron chemical potential" and refer to
:math:`\mu_{p,2}` as the "proton chemical potential" even though
charge neutrality has been assumed so the electron density is not
independent. The tables generated at this website use the same
notation.

The distinction between :math:`\mu_{n,1}` and :math:`\mu_{n,3}` is
more complicated, see Eq. 36 of Du et al. (2022).

The neutron fraction ``Xn`` stored in the
table refers only to neutrons *outside* of nuclei, i.e. :math:`X_n
\equiv n_n/n_B \neq \bar{n}_n/n_B`.
