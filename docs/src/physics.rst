Physics Overview
=======================
Our EoS describes the relationship between energy density and pressure
in a large three-dimensional (:math:`\mathrm{nB}`, :math:`\mathrm{Ye}`,
:math:`\mathrm{T}`) space,

There are several physical regimes each constrained by
different observables and theoretical approaches.

First the zero temperature nuclear matter at nuclear saturation
density is closely connected to nuclear masses, charge radii, giant
resonances, and other laboratory observables. Global fits to
experimental data have been performed with Skyrme and covariant
mean-field models.

Secondly, cold neutron matter below nuclear saturation density, is
difficult to probe experimentally but is well-constrained by
theoretical calculations based on semi-phenomenological nuclear forces
or microscopic chiral effective field theory-based interactions.

The third regime, strongly-interacting high-temperature matter, is
best described by interactions and many-body approaches similar to
those applied to cold neutron matter near saturation density.

The fourth regime, low-density and high-temperature matter that is
nearly non-degenerate, is best described by the virial expansion. The
equation of state in this regime is determined from nucleon-nucleon
scattering phase shifts.

Finally, neutron-rich matter at densities above twice saturation
density is most strongly constrained by observations of neutron star
masses and radii, particularly the observation of neutron stars with
:math:`M \simeq 2M_\odot`.

In this work, we construct a phenomenological free energy density that
is consistent with observational and theoretical constraints in the
five aforementioned physical regimes. This is in contrast to works
which attempt to describe matter over the entire density and
temperature range with a single detailed model of the nucleon-nucleon
interaction.

Our second advance is in the treatment of uncertainties. The most
relevant parameters which describe the uncertainties in different
density and temperature regimes are not clearly related. In this work,
through the construction of a phenomenological model one can vary
uncertainties in different regimes independently, without spoiling
agreement elsewhere.

The module produces an EoS table in a 3 dimensional grid of
:math:`\mathrm{nB}`, :math:`\mathrm{Ye}` and
:math:`\mathrm{T}` and stores various
physical and thermodynamics properties in the data files (detailed
documentation
[here](https://neutronstars.utk.edu/code/eos/table_format.html#thermodynamic-quantities)).

Rest mass contribution
-----------------------

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

Variable Transformations
------------------------

It is useful to be able to convert derivative operators between the
various sets of composition variables. In the relations below, we omit
the "bars" and simply write :math:`n_n,n_p` for
:math:`\bar{n}_n,\bar{n}_p`. In other words, all of the nucleon
densities below are presumed to include nucleons both inside and
outside of nuclei.

Converting between :math:`(n_n,n_p)` and :math:`(n_B,n_e)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since :math:`n_p=n_e` and
:math:`n_n=n_B-n_e`,

.. math::
   
   \left(\frac{\partial }{\partial n_B}\right)_{n_e} &=& 
   \left(\frac{\partial n_n}{\partial n_B}\right)_{n_e}
   \left(\frac{\partial }{\partial n_n}\right)_{n_p} +
   \left(\frac{\partial n_p}{\partial n_B}\right)_{n_e}
   \left(\frac{\partial }{\partial n_p}\right)_{n_n} =
   \left(\frac{\partial }{\partial n_n}\right)_{n_p}
   \nonumber \\
   \left(\frac{\partial }{\partial n_e}\right)_{n_B} &=& 
   \left(\frac{\partial n_n}{\partial n_e}\right)_{n_B}
   \left(\frac{\partial }{\partial n_n}\right)_{n_p} +
   \left(\frac{\partial n_p}{\partial n_e}\right)_{n_B}
   \left(\frac{\partial }{\partial n_p}\right)_{n_n} =
   \left(\frac{\partial }{\partial n_p}\right)_{n_n} -
   \left(\frac{\partial }{\partial n_n}\right)_{n_p}

For second derivatives

.. math::
   
   \left(\frac{\partial^2 }{\partial n_B^2}\right)_{n_e} &=& 
   \left(\frac{\partial^2 }{\partial n_n^2}\right)_{n_p}
   \nonumber \\
   \left(\frac{\partial^2 }{\partial n_e\partial n_B}\right) &=& 
   \left(\frac{\partial^2 }{\partial n_p \partial n_n}\right) -
   \left(\frac{\partial^2 }{\partial n_n^2}\right)_{n_p}
   \nonumber \\
   \left(\frac{\partial^2 }{\partial n_e^2}\right)_{n_B} &=& 
   \left(\frac{\partial^2 }{\partial n_p^2}\right)_{n_n} -
   2\left(\frac{\partial^2 }{\partial n_p \partial n_n}\right) +
   \left(\frac{\partial^2 }{\partial n_n^2}\right)_{n_p}
   
Converting between :math:`(n_n,n_p)` and :math:`(n_B,Y_e)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since :math:`n_p=n_B Y_e` and :math:`n_n=n_B(1-Y_e)`,

.. math::
   
   \left(\frac{\partial }{\partial n_B}\right)_{Y_e} &=& 
   \left(\frac{\partial n_n}{\partial n_B}\right)_{Y_e}
   \left(\frac{\partial }{\partial n_n}\right)_{n_p} +
   \left(\frac{\partial n_p}{\partial n_B}\right)_{Y_e}
   \left(\frac{\partial }{\partial n_p}\right)_{n_n} =
   (1-Y_e) \left(\frac{\partial }{\partial n_n}\right)_{n_p} +
   Y_e \left(\frac{\partial }{\partial n_p}\right)_{n_n}
   \nonumber \\
   \left(\frac{\partial }{\partial Y_e}\right)_{n_B} &=& 
   \left(\frac{\partial n_n}{\partial Y_e}\right)_{n_B}
   \left(\frac{\partial }{\partial n_n}\right)_{n_p} +
   \left(\frac{\partial n_p}{\partial Y_e}\right)_{n_B}
   \left(\frac{\partial }{\partial n_p}\right)_{n_n} =
   n_B \left[\left(\frac{\partial }{\partial n_p}\right)_{n_n} -
   \left(\frac{\partial }{\partial n_n}\right)_{n_p} \right]

The inverse transformation is:

.. math::

   \left(\frac{\partial }{\partial n_n}\right)_{n_p} =
   \left(\frac{\partial }{\partial n_B}\right)_{Y_e}
   - \frac{Y_e}{n_B}
   \left(\frac{\partial }{\partial Y_e}\right)_{n_B}
   \nonumber \\
   \left(\frac{\partial }{\partial n_p}\right)_{n_n} =
   \left(\frac{\partial }{\partial n_B}\right)_{Y_e}
   + \frac{(1-Y_e)}{n_B}
   \left(\frac{\partial }{\partial Y_e}\right)_{n_B}

This transformation is used in ``stability()`` in ``eos_nuclei.cpp``.
There is a Maxwell relation:

.. math::

   \frac{\partial f^2}{\partial n_n \partial n_p} = 
   \frac{\partial f^2}{\partial n_p \partial n_n}

which implies    

.. math::

   \left(\frac{\partial \mu_n}{\partial n_p}\right) = 
   \left(\frac{\partial \mu_p}{\partial n_n}\right)
   \left(\frac{\partial \mu_e}{\partial n_n}\right)

or    

.. math::

   \left(\frac{\partial \mu_p}{\partial n_B}\right)_{Y_e}
   - \frac{Y_e}{n_B}
   \left(\frac{\partial \mu_p}{\partial Y_e}\right)_{n_B}
   =
   \left(\frac{\partial \mu_n}{\partial n_B}\right)_{Y_e}
   + \frac{(1-Y_e)}{n_B}
   \left(\frac{\partial \mu_n}{\partial Y_e}\right)_{n_B}

thus   

.. math::

   \left(\frac{\partial \mu_p}{\partial n_B}\right)_{Y_e}
   = 
   \frac{Y_e}{n_B}
   \left(\frac{\partial \mu_p}{\partial Y_e}\right)_{n_B}
   + \left(\frac{\partial \mu_n}{\partial n_B}\right)_{Y_e}
   + \frac{(1-Y_e)}{n_B}
   \left(\frac{\partial \mu_n}{\partial Y_e}\right)_{n_B}

This equality is also used in ``stability()`` in ``eos_nuclei.cpp``.

Converting between :math:`(n_n,n_p)` and :math:`(n_B,n_e)` with muons
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When muons are included, the expressions change, since :math:`n_p =
n_e + n_{\mu}(n_e)` and :math:`n_n = n_B - n_e - n_{\mu}(n_e)`,

.. math::

   \left(\frac{\partial }{\partial n_B}\right)_{n_e} &=& 
   \left(\frac{\partial n_n}{\partial n_B}\right)_{n_e}
   \left(\frac{\partial }{\partial n_n}\right)_{n_p} +
   \left(\frac{\partial n_p}{\partial n_B}\right)_{n_e}
   \left(\frac{\partial }{\partial n_p}\right)_{n_n} =
   \left(\frac{\partial }{\partial n_n}\right)_{n_p}
   \nonumber \\
   \left(\frac{\partial }{\partial n_e}\right)_{n_B} &=& 
   \left(\frac{\partial n_n}{\partial n_e}\right)_{n_B}
   \left(\frac{\partial }{\partial n_n}\right)_{n_p} +
   \left(\frac{\partial n_p}{\partial n_e}\right)_{n_B}
   \left(\frac{\partial }{\partial n_p}\right)_{n_n} =
   (1+\chi) \left[
   \left(\frac{\partial }{\partial n_p}\right)_{n_n} -
   \left(\frac{\partial }{\partial n_n}\right)_{n_p}\right]

where

.. math::
   
   \chi = \frac{\partial n_{\mu}}{\partial n_e} =
   \frac{\partial n_{\mu}}{\partial {\mu}_{\mu}}
   \frac{\partial {\mu}_{\mu}}{\partial {\mu}_e}
   \frac{\partial {\mu}_{e}}{\partial n_e} +
   \frac{\partial {\mu}_{e}}{\partial n_e} = 
   \frac{\partial n_{\mu}}{\partial {\mu}_{\mu}}
   \left(\frac{\partial n_e}{\partial {\mu}_{e}}\right)^{-1}

For second derivatives

.. math::
   
   \left(\frac{\partial^2 }{\partial n_B^2}\right)_{n_e} &=& 
   \left(\frac{\partial^2 }{\partial n_n^2}\right)_{n_p}
   \nonumber \\
   \left(\frac{\partial^2 }{\partial n_e\partial n_B}\right) &=& 
   (1+\chi)\left[\left(\frac{\partial^2 }{\partial n_p \partial n_n}\right) -
   \left(\frac{\partial^2 }{\partial n_n^2}\right)_{n_p}\right]
   \nonumber \\
   \left(\frac{\partial^2 }{\partial n_e^2}\right)_{n_B} &=&
   \left(1+\chi\right)^2 \left[
   \left(\frac{\partial^2 }{\partial n_p^2}\right)_{n_n} -
   2\left(\frac{\partial^2 }{\partial n_p \partial n_n}\right) +
   \left(\frac{\partial^2 }{\partial n_n^2}\right)_{n_p}\right]
   

Chemical Potentials
-------------------

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

Speed of sound in a multicomponent system
-----------------------------------------

Using :math:`\varepsilon` for energy density :math:`S` for entropy,
:math:`s` for entropy density, and :math:`\tilde{s}` for entropy per
baryon, the speed of sound is

.. math::
   
   c_s^2 = \left( \frac{\partial P}{\partial \varepsilon}
   \right)_{\tilde{s},\{ N_i \}}
   \, .

The energy density in the denominator must *include the rest mass
contribution to the energy density*. In infinite matter, it is useful
to rewrite this derivative in terms of fixed volume rather than fixed
number.

.. math::
   
   c_s^2 = \left( \frac{\partial P}{\partial \varepsilon}
   \right)_{S,\{ N_i \}} =
   \left( \frac{\partial P}{\partial V} \right)_{S,\{ N_i \}}
   \left( \frac{\partial \varepsilon}{\partial V} \right)_{S,\{ N_i \}}^{-1}
 
The second derivative is

.. math::
   
   \left( \frac{\partial \varepsilon}{\partial V} \right)_{S,\{ N_i \}} = 
   \left[ \frac{\partial  (E/V)}{\partial V} \right]_{S,\{ N_i \}} =
   -\frac{1}{V} P - \frac{E}{V^2} = - \frac{P+\varepsilon}{V}
   = - \frac{T s + \sum_i \mu_i n_i}{V}
 
and first derivative is

.. math::
   
   \left( \frac{\partial P}{\partial V} \right)_{S,\{ N_j \}} &=& -
   \left( \frac{\partial \varepsilon}{\partial V} \right)_{S,\{ N_j\}} +
   S \left[ \frac{\partial (T/V)}{\partial V} \right]_{S,\{ N_j \}} +
   \sum_i 
   N_i \left[ \frac{\partial  (\mu_i/V)}{\partial V} \right]_{S,\{ N_j \}}
   \nonumber \\ &=& -
   \left( \frac{\partial \varepsilon}{\partial V} \right)_{S,\{ N_j \}} +
   S \left[ -\frac{T}{V^2} + \left( \frac{\partial T}{\partial V}
   \right)_{S,\{ N_j \}}\right] +
   \sum_i 
   N_i \left[ -\frac{\mu_i}{V^2} +
   \left( \frac{\partial \mu_i}{\partial V} \right)_{S,\{ N_j \}}\right]
   \nonumber \\ &=& \frac{P + \varepsilon}{V} +
   S \left[ -\frac{T}{V^2} - \left( \frac{\partial P}{\partial S}
   \right)_{\{N_j\},V}\right] +
   \sum_i N_i \left[ -\frac{\mu_i}{V^2} -
   \left( \frac{\partial P}{\partial N_i}
   \right)_{S,\{N_{j\neq i}\},V}\right] \nonumber \\
   &=& - S \left( \frac{\partial P}{\partial S}\right)_{\{n_j\},V}
   - \sum_i N_i \left( \frac{\partial P}{\partial N_i}
   \right)_{S,\{n_{j\neq i}\},V}
     
Putting these two results together gives

.. math::
   
   c_s^2 = \left[s \left( \frac{\partial P}{\partial s}
   \right)_{\{n_j\},V} +
   \sum_i n_i \left( \frac{\partial P}
   {\partial n_i} \right)_{S,\{n_{j\neq i}\},V}\right] \left(
   T s + \sum_i \mu_i n_i \right)^{-1}
 
To re-express this in terms of derivatives of the free energy
(which again must include the rest mass contribution),

.. math::
   
   c_s^2 = \left\{s \left[ \frac{\partial (\sum_i \mu_i n_i - f)}
   {\partial s} \right]_{\{n_j\},V} +
   \sum_i 
   n_i\left[ \frac{\partial  ( \sum_k \mu_k n_k - f)}{\partial n_i}
   \right]_{s,\{n_{j\neq i}\},V}\right\} \left(
   T s + \sum_i \mu_i n_i \right)^{-1}
   
For the sum over :math:`k`,
all densities are constant except for :math:`n_i`, thus

.. math::
   
   \sum_i 
   n_i \frac{\partial}{\partial n_i}
   \left( \sum_k \mu_k n_k - f \right)_{s,\{n_{j\neq i}\},V}
   &=& \sum_i n_i \frac{\partial}{\partial n_i}
   \left( \sum_{k\neq i} \mu_k n_k + \mu_i n_i -f
   \right)_{s,\{n_{j\neq i}\},V} \nonumber \\
   &=& 
   \sum_i \left[ \sum_k n_k \left(\frac{\partial \mu_k }
   {\partial n_i}\right)_{s,\{n_{j\neq i}\},V} + \mu_i -
   \left(\frac{\partial f}{\partial n_i}\right)_{s,\{n_{j\neq i}\},V}
   \right]
 
To compute this we need

.. math::
   
   \left(\frac{\partial f}{\partial n_i}\right)_{s,\{n_{j\neq i}\},V} &=&
   \left(\frac{\partial f}{\partial n_i}\right)_{\{n_{j\neq i}\},T,V} +
   \left(\frac{\partial f}{\partial T}\right)_{n_B,\{n_{j\neq i}\},V}
   \left(\frac{\partial T}{\partial n_i}\right)_{\{n_{j\neq i}\},s,V}
   \nonumber \\
   &=& \mu_i - s \left(\frac{\partial T}{\partial n_i}
   \right)_{\{n_{j\neq i}\},s,V}
   \nonumber \\
   \left(\frac{\partial \mu_k}{\partial n_i}\right)_{s,\{n_{j\neq i}\},V} &=&
   \left(\frac{\partial \mu_k}{\partial n_i}\right)_{\{n_{j\neq i}\},T,V} +
   \left(\frac{\partial \mu_k}{\partial T}\right)_{n_i,\{n_{j\neq i}\},V}
   \left(\frac{\partial T}{\partial n_i}\right)_{\{n_{j\neq i}\},s,V} 
   \nonumber \\
   &=& f_{n_i n_k} + f_{n_k T}
   \left(\frac{\partial T}{\partial n_i}\right)_{\{n_{j\neq i}\},s,V}
 
which requires

.. math::
   
   \left(\frac{\partial T}{\partial n_i}\right)_{\{n_{j\neq i}\},s,V}
   = -\left(\frac{\partial s}{\partial n_i}\right)_{\{n_{j\neq i}\},T,V}
   \left(\frac{\partial s}{\partial T}\right)_{\{n\},V}^{-1}
   = -f_{n_i T}/f_{TT}
 
Finally, we get

.. math::
   
   c_s^2 = \left\{
   - \left(\frac{s}{f_{TT}}\right) \left( \sum_i n_i f_{n_i T}+s \right)
   + \sum_i n_i \left[ \sum_k n_k \left(f_{n_i n_k}- f_{n_k T}
   f_{n_i T} f_{TT}^{-1}\right) 
   - s f_{n_i T} f_{TT}^{-1}\right]
   \right\} \\
   \times \left(
   T s + \sum_i \mu_i n_i \right)^{-1}

and   
   
.. math::
   
   c_s^2 = \left[
   \sum_i \sum_k n_i n_k \left(f_{n_i n_k}- f_{n_k T}
   f_{n_i T} f_{TT}^{-1}\right)
   - 2\sum_i s n_i f_{n_i T} f_{TT}^{-1}
   - s^2 f_{TT}^{-1} \right] \left(
   T s + \sum_i \mu_i n_i \right)^{-1}
 
Note that, when applying this expression, one must be consistent about
the free energy which one differentiates and the densities and
chemical potentials which are used. See :ref:`Chemical Potentials` for
more information regarding this issue. 

WIP
----
First, we define the symmetry energy to include a zero
temperature contribution which combines the QMC EOS 
near saturation density, the neutron star fit at higher densities, 
and the Skyrme interaction for isospin-symmetric matter

.. math::

   \epsilon_{sym}(n_B) = h(n_B)\epsilon_{QMC}(nB) + [1-h(n_B)]\epsilon_{NS}(n_B) - f_{Skyrme}(nB,x_p = 1/2, T=0)


Defining the isospin asymmetry $ \delta = 1-2x_p$, we can combine this
 with the model described in `Du et al. 2019 <https://arxiv.org/pdf/1802.09710>`_ to obtain the free energy
 density of degenerate matter

.. math::

   f_{deg}(n_B,x_p,T) = f_{Skyrme}(nB,x_p = 1/2, T=0) + \delta^2\epsilon_{sym}(n_B) \\
   + \delta^2\Delta f_{hot}(nB,x_p = 0, T) + (1-\delta^2)\Delta f_{hot}(nB,x_p = 1/2, T)


Finally, we ensure that the total nucleonic free energy gives the
 result from the virial expansion at high tem peratures using

.. math::

   f_{np}(n_B,x_p,T) = f_{virial}(n_n,x_p,T)g+f_{deg}(n_B,x_p,T)(1-g)


When we need to include the
 electrons, positrons, and photons, we define the free en
ergy density

.. math::

   f_{npe\gamma} \equiv f_{np} + f_{e^-}+f_{e^+}+f_\gamma


Using this formalism, the chemical potentials and entropy can be computed directly (eq. 28-32 in `Du et al. 2019 <https://arxiv.org/pdf/1802.09710>`_).

We enforce causality at high densities.