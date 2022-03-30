Speed of sound in a multicomponent system
=========================================

Using :math:`\varepsilon` for energy density (including rest mass
energy density), :math:`S` for entropy, :math:`s` for entropy density,
and :math:`\tilde{s}` for entropy per baryon, the speed of sound is

.. math::
   
   c_s^2 = \left( \frac{\partial P}{\partial \varepsilon}
   \right)_{\tilde{s},\{ N_i \}}
   \, .
 
In infinite matter, it is useful to rewrite this derivative in
terms of fixed volume rather than fixed number.

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
 
To re-express this in terms of derivatives of the free energy,

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
   = \mu_i - s \left(\frac{\partial T}{\partial n_i}
   \right)_{\{n_{j\neq i}\},s,V}
   \nonumber \\
   \left(\frac{\partial \mu_k}{\partial n_i}\right)_{s,\{n_{j\neq i}\},V} &=&
   \left(\frac{\partial \mu_k}{\partial n_i}\right)_{\{n_{j\neq i}\},T,V} +
   \left(\frac{\partial \mu_k}{\partial T}\right)_{n_i,\{n_{j\neq i}\},V}
   \left(\frac{\partial T}{\partial n_i}\right)_{\{n_{j\neq i}\},s,V} 
   = f_{n_i n_k} + f_{n_k T}
   \left(\frac{\partial T}{\partial n_i}\right)_{\{n_{j\neq i}\},s,V}
 
which requires

.. math::
   
   \left(\frac{\partial T}{\partial n_i}\right)_{\{n_{j\neq i}\},s,V}
   = -\left(\frac{\partial s}{\partial n_i}\right)_{\{n_{j\neq i}\},T,V}
   \left(\frac{\partial s}{\partial T}\right)_{\{n\},V}^{-1}
   = -f_{n_i T}/f_{TT}
 
Finally, we get

.. math::
   
   c_s^2 &=& \left\{
   - \left(\frac{s}{f_{TT}}\right) \left( \sum_i n_i f_{n_i T}+s \right)
   + \sum_i n_i \left[ \sum_k n_k \left(f_{n_i n_k}- f_{n_k T}
   f_{n_i T} f_{TT}^{-1}\right) 
   - s f_{n_i T} f_{TT}^{-1}\right]
   \right\} \left(
   T s + \sum_i \mu_i n_i \right)^{-1} \nonumber \\
   &=& \left[
   \sum_i \sum_k n_i n_k \left(f_{n_i n_k}- f_{n_k T}
   f_{n_i T} f_{TT}^{-1}\right)
   - 2\sum_i s n_i f_{n_i T} f_{TT}^{-1}
   - s^2 f_{TT}^{-1} \right] \left(
   T s + \sum_i \mu_i n_i \right)^{-1}
 
