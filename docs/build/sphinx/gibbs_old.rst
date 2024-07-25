Conserved charges and composition
=================================

In a core-collapse or merger, we can write the Gibbs energy
density as

.. math::
   
   g= \sum_i \mu_i n_i = \mu_n n_n + \mu_p n_p + \mu_e n_e
   + \mu_{\mu} n_{\mu} 
   + \mu_{\nu_e} (n_{\nu_e} - n_{\bar{\nu}_e}) 
   + \mu_{\nu_{\mu}} (n_{\nu_{\mu}} - n_{\bar{\nu}_{\mu}}) \, .

where :math:`n_e \equiv n_{e^{-}} - n_{e^{+}}` and :math:`n_{\mu}
\equiv n_{{\mu}^{-}} - n_{{\mu}^{+}}`. If we ignore neutrino mixing,
then electron lepton number and muon lepton number are conserved, so
there are associated chemical potential, :math:`\mu_{Le}` and
:math:`\mu_{L\mu}`. Thus :math:`\mu_n = \mu_B`, :math:`\mu_p = \mu_B +
\mu_Q`, :math:`\mu_e = - \mu_Q + \mu_{Le}`, :math:`\mu_{\nu_e} =
\mu_{Le}`, :math:`\mu_{\mu} = - \mu_Q + \mu_{L{\mu}}` and
:math:`\mu_{\nu_{\mu}} = \mu_{L{\mu}}`. This gives

.. math::
   
   g= \sum_i \mu_i n_i = \mu_B n_n + (\mu_B+\mu_Q) n_p +
   (-\mu_{Q}+\mu_{Le}) n_e
   + \mu_{Le} (n_{\nu_e} - n_{\bar{\nu}_e}) +
   (-\mu_{Q}+\mu_{L{\mu}}) n_{\mu}
   + \mu_{L{\mu}} (n_{\nu_{\mu}} - n_{\bar{\nu}_{\mu}}) \, .

Presuming that the system is locally charge neutral, :math:`n_p =
n_e + n_{\mu}`, then

.. math::
   
   g= \sum_i \mu_i n_i = \mu_B n_B + 
   \mu_{Le} (n_e + n_{\nu_e} - n_{\bar{\nu}_e}) +
   \mu_{L{\mu}} (n_{\mu} + n_{\nu_{\mu}} - n_{\bar{\nu}_{\mu}})
   = \mu_B n_B + \mu_{Le} n_{Le} + \mu_{L{\mu}} n_{L{\mu}}

where :math:`n_{Li}` is the number density of :math:`i`-type leptons.
Since there are two conserved charges, then we need two composition
variables plus the temperature to specify the EOS, and these are
typically chosen to be :math:`n_B`, :math:`Y_e\equiv n_e/n_B` and
:math:`T`. When neutrinos are trapped, the Gibbs energy density cannot
be simplified further. When neutrinos are not trapped, then
:math:`\mu_{Le} \neq \mu_{\nu_e}` and :math:`\mu_{L\mu} \neq
\mu_{\nu_\mu}` because electron lepton number is no longer conserved.
In this case :math:`\mu_{\nu_e} = \mu_{\nu_{\mu}} = 0` and
:math:`\mu_e = \mu_{\mu}`, thus :math:`g=n_B \mu_B + n_e \mu_{Le}+
n_{\mu} \mu_{L{\mu}} = n_B \mu_B + (n_e+n_{\mu}) ( \mu_p + \mu_e -
\mu_n)`. Alternatively, one can write :math:`g=n_n \mu_n + n_p
\hat{\mu}_p` where :math:`\hat{\mu}_p \equiv \mu_p + \mu_e`. Treating
the neutrinos as massless, then
:math:`k_{F,\nu_e}=k_{F,\bar{\nu}_e}=0` and thus
:math:`\mu_{\nu_e}=0`. When weak reactions are fast enough, beta
equilibrium holds, which gives
:math:`\mu_n + \mu_{\nu_e} = \mu_p + \mu_e` whether or not neutrinos
are trapped.
