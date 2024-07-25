# Documentation
First, we define the symmetry energy to include a zero
temperature contribution which combines the QMC EOS
 near saturation density, the neutron star fit at higher den
sities, and the Skyrme interaction for isospin-symmetric
 matter

$$\begin{equation}
\epsilon_{sym}(n_B) = h(n_B)\epsilon_{QMC}(nB) + [1-h(n_B)]\epsilon_{NS}(n_B) - f_{Skyrme}(nB,x_p = 1/2, T=0)
\end{equation}$$

Defining the isospin asymmetry $ \delta = 1-2x_p$, we can combine this
 with the model described in [Du et al.
 2019](https://arxiv.org/pdf/1802.09710) to obtain the free energy
 density of degenerate matter

$$\begin{equation}
f_{deg}(n_B,x_p,T) = f_{Skyrme}(nB,x_p = 1/2, T=0) + \delta^2\epsilon_{sym}(n_B) + \delta^2\Delta f_{hot}(nB,x_p = 0, T) + (1-\delta^2)\Delta f_{hot}(nB,x_p = 1/2, T)
\end{equation}$$

Finally, we ensure that the total nucleonic free energy gives the
 result from the virial expansion at high tem peratures using

$$\begin{equation}
f_{np}(n_B,x_p,T) = f_{virial}(n_n,x_p,T)g+f_{deg}(n_B,x_p,T)(1-g)
\end{equation}$$

When we need to include the
 electrons, positrons, and photons, we define the free en
ergy density

$$\begin{equation}
f_{npe\gamma} \equiv f_{np} + f_{e^-}+f_{e^+}+f_\gamma
\end{equation}$$

Using this formalism, the chemical potentials and entropy can be computed directly (eq. 28-32 in [Du et al. 2019](https://arxiv.org/pdf/1802.09710)).

We enforce causality at high densities.

## Local Sphinx build and website preview

To build and preview your documentation using Sphinx, run the build script to build and launch a Docker container like so:

```bash
bash docs/build/build.sh
```

Visit http://localhost:4000/docs/ in your browser to view the rendered website.
