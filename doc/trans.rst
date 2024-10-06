Variable transformations
========================

It is useful to be able to convert derivative operators between the
various sets of composition variables. In the relations below, we omit
the "bars" and simply write :math:`n_n,n_p` for
:math:`\bar{n}_n,\bar{n}_p`. In other words, all of the nucleon
densities below are presumed to include nucleons both inside and
outside of nuclei.

Converting between (nn,np) and (nB,ne)
--------------------------------------

Since :math:`n_p=n_e` and :math:`n_n=n_B-n_e`,

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
   
Converting between (nn,np) and (nB,Ye)
--------------------------------------

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

Converting between (nn,np) and (nB,Ye) with muons
-------------------------------------------------

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
   

