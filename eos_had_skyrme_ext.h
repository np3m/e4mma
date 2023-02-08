/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2023, Xingfu Du, Zidu Lin, and Andrew W. Steiner
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
/** \file eos_had_skyrme.h
    \brief File defining \ref o2scl::eos_had_skyrme
*/
#ifndef O2SCL_EOS_HAD_SKYRME_H
#define O2SCL_EOS_HAD_SKYRME_H

#include <iostream>
#include <string>
#include <cmath>

#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/part.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/fermion_deriv_nr.h>

  /** \brief Extended Skyrme hadronic equation of state 

      This is the modified Skyrme model proposed by Holt and Lim.
   */
  class eos_had_skyrme_ext : public o2scl::eos_had_skyrme {
    
  protected:

    /** \brief Compute the base thermodynamic quantities

	This function computes the energy density, pressure,
	entropy, and chemical potentials.
     */
    template<class fermion_t>
      void base_thermo
      (fermion_t &ne, fermion_t &pr, double ltemper, o2scl::thermo &locth,
       double term, double term2, double ham1, double ham2,
       double ham3, double ham4, double ham5, double ham6) {

      double nb=ne.n+pr.n;
      double nba1=pow(nb,alpha);
      double nba2=pow(nb,alpha2);
      double nba3=pow(nb,alpha3);
      double nna1=pow(ne.n,alpha);
      double nna2=pow(ne.n,alpha2);
      double nna3=pow(ne.n,alpha3);
      double npa1=pow(pr.n,alpha);
      double npa2=pow(pr.n,alpha2);
      double npa3=pow(pr.n,alpha3);
      
      double ham32=a*t4/6.0*(1.0+0.5*x4);
      double ham42=a*t4*pow(2.0,alpha2-2.0)/6.0*(1.0-x4);
      double ham52=b*t4/12.0*(1.0+0.5*x4);
      double ham62=-b*t4/12.0*(0.5+x4);

      double ham33=a*t5/6.0*(1.0+0.5*x5);
      double ham43=a*t5*pow(2.0,alpha3-2.0)/6.0*(1.0-x5);
      double ham53=b*t5/12.0*(1.0+0.5*x5);
      double ham63=-b*t5/12.0*(0.5+x5);
      
      double ham=ne.ed+pr.ed+ham1*nb*nb+ham2*(ne.n*ne.n+pr.n*pr.n)+
	ham3*nba1*ne.n*pr.n+ham4*(nna1*ne.n*ne.n+npa1*pr.n*pr.n)+
	ham5*nb*nb*nba1+ham6*(ne.n*ne.n+pr.n*pr.n)*nba1+
	ham32*nba2*ne.n*pr.n+ham42*(nna2*ne.n*ne.n+npa2*pr.n*pr.n)+
	ham52*nb*nb*nba2+ham62*(ne.n*ne.n+pr.n*pr.n)*nba2+
	ham33*nba3*ne.n*pr.n+ham43*(ne.n*ne.n*ne.n+pr.n*pr.n*pr.n)+
	ham53*nb*nb*nba3+ham63*(ne.n*ne.n+pr.n*pr.n)*nba3;
      
      double gn, gp;
      if (ne.inc_rest_mass) {
	gn=2.0*ne.ms*(ne.ed-ne.n*ne.m);
      } else {
	gn=2.0*ne.ms*ne.ed;
      }
      if (pr.inc_rest_mass) {
	gp=2.0*pr.ms*(pr.ed-pr.n*pr.m);
      } else {
	gp=2.0*pr.ms*pr.ed;
      }
      
      // Variables dhdn{n,p} are the partial derivatives of the
      // Hamiltonian wrt the neutron and proton densities
      double common=2.0*ham1*nb+ham5*(2.0+alpha)*nb*nba1+
	ham52*(2.0+alpha2)*nb*nba2+ham53*(2.0+alpha3)*nb*nba3;
      
      double dhdnn=common+2.0*ham2*ne.n+
	ham3*nba1*pr.n*(alpha*ne.n/nb+1.0)+ham4*(nna1*ne.n*(2.0+alpha))+
	ham6*(2.0*ne.n*nba1+(ne.n*ne.n+pr.n*pr.n)*alpha*nba1/nb)+
	ham32*nba2*pr.n*(alpha2*ne.n/nb+1.0)+ham42*(nna2*ne.n*(2.0+alpha2))+
	ham62*(2.0*ne.n*nba2+(ne.n*ne.n+pr.n*pr.n)*alpha2*nba2/nb)+
	ham33*nba3*pr.n*(alpha3*ne.n/nb+1.0)+ham43*(nna3*ne.n*(2.0+alpha3))+
	ham63*(2.0*ne.n*nba3+(ne.n*ne.n+pr.n*pr.n)*alpha3*nba3/nb);
      
      double dhdnp=common+2.0*ham2*pr.n+
	ham3*nba1*ne.n*(alpha*pr.n/nb+1.0)+ham4*(npa1*pr.n*(2.0+alpha))+
	ham6*(2.0*pr.n*nba1+(ne.n*ne.n+pr.n*pr.n)*alpha*nba1/nb)+
	ham32*nba2*ne.n*(alpha2*pr.n/nb+1.0)+ham42*(npa2*pr.n*(2.0+alpha2))+
	ham62*(2.0*pr.n*nba2+(ne.n*ne.n+pr.n*pr.n)*alpha2*nba2/nb)+
	ham33*nba3*ne.n*(alpha3*pr.n/nb+1.0)+ham43*(npa3*pr.n*(2.0+alpha3))+
	ham63*(2.0*pr.n*nba3+(ne.n*ne.n+pr.n*pr.n)*alpha3*nba3/nb);

      // Compute the chemical potentials
      ne.mu=ne.nu+dhdnn+(gn+gp)*term+gn*term2;
      pr.mu=pr.nu+dhdnp+(gn+gp)*term+gp*term2;
      
      // Thermodynamics
      locth.ed=ham;
      locth.en=ne.en+pr.en;
      locth.pr=ltemper*locth.en+ne.mu*ne.n+pr.mu*pr.n-locth.ed;
      
      return;
    }
    
  public:

    /// \name Basic usage
    //@{
    /// Create a blank Skyrme EOS
    eos_had_skyrme_ext();

    /// Destructor
    virtual ~eos_had_skyrme_ext() {};

    /** \brief Equation of state as a function of densities
	at finite temperature
    */
    virtual int calc_temp_e(o2scl::fermion &ne, o2scl::fermion &pr,
			    double temper, o2scl::thermo &th);

    /** \brief Equation of state as a function of densities at 
	zero temperature
    */
    virtual int calc_e(o2scl::fermion &ne, o2scl::fermion &pr,
		       o2scl::thermo &lt);

    /// Return string denoting type ("eos_had_skyrme_ext")
    virtual const char *type() { return "eos_had_skyrme_ext"; }
    //@}

    /// \name The new parameters
    //@{
    double t4, x4, t5, x5;
    double alpha2, alpha3;
    //@}

  };

  /** \brief Extended Skyrme hadronic equation of state 

      This is the modified Skyrme model proposed by Holt and Lim.

      The Hamiltonian is 
      \f[
      {\cal H} = \frac{\hbar^2 \tau_n}{2 m_n^{*}}
      + \frac{\hbar^2 \tau_p}{2 m_p^{*}} + 
      {\cal H}_{\mathrm{pot}}
      \f]
      where \f$ \tau_i = k_{F,i}^5/(5 \pi^2) \f$,
      \f[
      \frac{1}{2 m_n^{*}} = \frac{1}{2 m_n} + f_n(n_n,n_p)
      \f]
      and 
      \f[
      \frac{1}{2 m_p^{*}} = \frac{1}{2 m_p} + f_p(n_n,n_p) \, .
      \f]
      Then, the chemical potentials are
      \f[
      \mu_n = \nu_n + \tau_n 
      \frac{\partial f_n}{\partial n_n} + 
      \tau_p 
      \frac{\partial f_p}{\partial n_n} 
      \f]
      and 
      \f[
      \mu_p = \nu_p + \tau_p 
      \frac{\partial f_p}{\partial n_p} + 
      \tau_n 
      \frac{\partial f_n}{\partial n_p} \, .
      \f]
   */
  class eos_had_lim_holt : public o2scl::eos_had_temp_eden_base {
    
  protected:

    /// Nonrelativistic fermion object
    o2scl::fermion_nonrel nrf;
    
  public:

    /// \name Parameters
    //@{
    double betaL, betaU, theta, thetaL, sigma;
    
    double alphaL, alphaU, etaL, etaU, zetaL, zetaU;

    double gamma, gamma2;
    //@}
    
  protected:

    /** \brief Compute the base thermodynamic quantities

	This function computes the energy density, pressure,
	entropy, and chemical potentials.
     */
    template<class fermion_t>
      void base_thermo
    (fermion_t &ne, fermion_t &pr, double ltemper, o2scl::thermo &locth) {

      return;
    }
    
  public:

    /// \name Basic usage
    //@{
    /// Create a blank Skyrme EOS
    eos_had_lim_holt();

    /// Destructor
    virtual ~eos_had_lim_holt() {};

    /** \brief Equation of state as a function of densities
	at finite temperature
    */
    virtual int calc_temp_e(o2scl::fermion &ne, o2scl::fermion &pr,
			    double temper, o2scl::thermo &th);

    /** \brief Equation of state as a function of densities at 
	zero temperature
    */
    virtual int calc_e(o2scl::fermion &ne, o2scl::fermion &pr,
		       o2scl::thermo &lt);
    
    virtual int calc_p(o2scl::fermion &ne, o2scl::fermion &pr,
		       o2scl::thermo &lt) {
      return 0;
    }

    /// Return string denoting type ("eos_had_lim_holt")
    virtual const char *type() { return "eos_had_lim_holt"; }
    //@}

  };

#endif
