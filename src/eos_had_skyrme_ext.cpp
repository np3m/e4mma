/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2024, Xingfu Du, Zidu Lin, and Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "eos_had_skyrme_ext.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_had_skyrme_ext::eos_had_skyrme_ext() {
}

int eos_had_skyrme_ext::calc_temp_e(fermion &ne, fermion &pr, 
                                    double ltemper, thermo &locth) {

#if !O2SCL_NO_RANGE_CHECK
  check_input(ne,pr,ltemper);
#endif

  double term, term2;
  eff_mass(ne,pr,term,term2);

  if (ne.ms<0.0 || pr.ms<0.0) {
    O2SCL_CONV2_RET("Effective masses negative in ",
		    "eos_had_skyrme_ext::calc_temp_e().",
		    exc_einval,this->err_nonconv);
  }

  // See note in class documentation about zero density
  if (ltemper>0.0 && ne.n==0.0) {
    if (ne.inc_rest_mass) {
      ne.nu=ne.m;
    } else {
      ne.nu=0.0;
    }
    ne.ed=0.0;
    ne.pr=0.0;
    ne.en=0.0;
  } else {
    nrf.calc_density(ne,ltemper);
  }
  if (ltemper>0.0 && pr.n==0.0) {
    if (pr.inc_rest_mass) {
      pr.nu=pr.m;
    } else {
      pr.nu=0.0;
    }
    pr.ed=0.0;
    pr.pr=0.0;
    pr.en=0.0;
  } else {
    nrf.calc_density(pr,ltemper);
  }

  // Compute the coefficients of different powers of density
  // in the hamiltonian
  double ham1, ham2, ham3, ham4, ham5, ham6;
  hamiltonian_coeffs(ham1,ham2,ham3,ham4,ham5,ham6);
  
  // Compute the base thermodynamic properties
  base_thermo(ne,pr,ltemper,locth,term,term2,
	      ham1,ham2,ham3,ham4,ham5,ham6);
  
  return success;
}

int eos_had_skyrme_ext::calc_e(fermion &ne, fermion &pr, thermo &locth) {
  return calc_temp_e(ne,pr,0.0,locth);
}

eos_had_lim_holt::eos_had_lim_holt() {
  fet=&nrf;
}

int eos_had_lim_holt::calc_temp_e(fermion &ne, fermion &pr, 
                                  double temper, thermo &th) {

  double nb=ne.n+pr.n;
  ne.inc_rest_mass=false;
  pr.inc_rest_mass=false;

  double fn=betaL*ne.n+betaU*pr.n+theta*nb*pow(nb,1.0+sigma)+
    thetaL*ne.n*pow(nb,sigma);
  double fp=betaL*pr.n+betaU*ne.n+theta*nb*pow(nb,1.0+sigma)+
    thetaL*pr.n*pow(nb,sigma);
  
  ne.ms=ne.m/(1.0+2.0*ne.m*fn);
  pr.ms=pr.m/(1.0+2.0*pr.m*fp);
  
  if (ne.ms<0.0 || pr.ms<0.0) {
    O2SCL_CONV2_RET("Effective masses negative in ",
		    "eos_had_lim_holt::calc_temp_e().",
		    exc_einval,this->err_nonconv);
  }

  // See note in class documentation about zero density
  if (temper>0.0 && ne.n==0.0) {
    if (ne.inc_rest_mass) {
      ne.nu=ne.m;
    } else {
      ne.nu=0.0;
    }
    ne.ed=0.0;
    ne.pr=0.0;
    ne.en=0.0;
  } else {
    nrf.calc_density(ne,temper);
  }
  if (temper>0.0 && pr.n==0.0) {
    if (pr.inc_rest_mass) {
      pr.nu=pr.m;
    } else {
      pr.nu=0.0;
    }
    pr.ed=0.0;
    pr.pr=0.0;
    pr.en=0.0;
  } else {
    nrf.calc_density(pr,temper);
  }
      
  double np2=ne.n*ne.n+pr.n*pr.n;
      
  double ham1=alphaL*np2;
  double ham2=2.0*alphaU*ne.n*pr.n;
  double ham3=etaL*np2*pow(nb,gamma);
  double ham4=2.0*etaU*ne.n*pr.n*pow(nb,gamma);
  double ham5=zetaL*np2*pow(nb,gamma2);
  double ham6=2.0*zetaU*ne.n*pr.n*pow(nb,gamma2);
      
  double ham=ne.ed+pr.ed;//+ham1+ham2+ham3+ham4+ham5+ham6;
      
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

  // The potential energy contributions to the chemical potentials
  double dhdnn=2.0*alphaL*ne.n+2.0*alphaU*pr.n+
    2.0*etaL*ne.n*pow(nb,gamma)+etaL*gamma*np2*pow(nb,gamma-1.0)+
    2.0*etaU*pr.n*pow(nb,gamma)+
    2.0*gamma*etaU*ne.n*pr.n*pow(nb,gamma-1.0)+
    2.0*zetaL*ne.n*pow(nb,gamma2)+gamma2*zetaL*np2*pow(nb,gamma2-1.0)+
    2.0*zetaU*pr.n*pow(nb,gamma2)+
    2.0*gamma2*zetaU*ne.n*pr.n*pow(nb,gamma2-1.0);
  double dhdnp=2.0*alphaL*pr.n+2.0*alphaU*ne.n+
    2.0*etaL*pr.n*pow(nb,gamma)+etaL*gamma*np2*pow(nb,gamma-1.0)+
    2.0*etaU*ne.n*pow(nb,gamma)+
    2.0*gamma*etaU*pr.n*ne.n*pow(nb,gamma-1.0)+
    2.0*zetaL*pr.n*pow(nb,gamma2)+gamma2*zetaL*np2*pow(nb,gamma2-1.0)+
    2.0*zetaU*ne.n*pow(nb,gamma2)+
    2.0*gamma2*zetaU*pr.n*ne.n*pow(nb,gamma2-1.0);
  dhdnn=0.0;
  dhdnp=0.0;

  // Compute the chemical potentials
  double der=theta*(1.0+sigma)*pow(nb,sigma);

  double dfndnn=betaL+der+thetaL*pow(nb,sigma)+
    thetaL*ne.n*sigma*pow(nb,sigma-1.0);
  double dfndnp=betaU+der+thetaL*ne.n*sigma*pow(nb,sigma-1.0);
  double dfpdnn=betaU+der+thetaL*pr.n*sigma*pow(nb,sigma-1.0);
  double dfpdnp=betaL+der+thetaL*pow(nb,sigma)+
    thetaL*pr.n*sigma*pow(nb,sigma-1.0);

  ne.mu=ne.nu+dhdnn+gn*dfndnn+gp*dfpdnn;
  pr.mu=pr.nu+dhdnp+gp*dfpdnp+gn*dfndnp;
  
  cout << "A: " << endl;
  cout << ne.mu*ne.n << " " << ne.ed+ne.pr << " " << temper << endl;
  cout << pr.mu*pr.n << " " << pr.ed+pr.pr << endl;
  /*
  ne.mu=ne.nu+gn*(betaL+der+thetaL*pow(nb,sigma)+
                  thetaL*ne.n*sigma*pow(nb,sigma-1.0))+
    gp*(betaU+der+sigma*thetaL*ne.n*sigma*pow(nb,sigma-1.0));
  pr.mu=pr.nu+gp*(betaL+der+thetaL*pow(nb,sigma)+
                  thetaL*pr.n*sigma*pow(nb,sigma-1.0))+
    gn*(betaU+der+sigma*thetaL*pr.n*sigma*pow(nb,sigma-1.0));
  */

  // Thermodynamics
  th.ed=ham;
  th.en=ne.en+pr.en;
  th.pr=temper*th.en+ne.mu*ne.n+pr.mu*pr.n-th.ed;

  return success;
}

int eos_had_lim_holt::calc_e(fermion &ne, fermion &pr, thermo &th) {
  ne.n=0.1;
  pr.n=0.06;
  int retx=calc_temp_e(ne,pr,0.0,th);
  double ed1=th.ed;
  //cout << ne.inc_rest_mass<< " " << pr.inc_rest_mass << endl;
  //cout << th.ed*hc_mev_fm/(ne.n+pr.n) << " ";
  cout << ne.mu*hc_mev_fm << " " << pr.mu*hc_mev_fm << endl;
  //cout << ne.ms*hc_mev_fm << " " << pr.ms*hc_mev_fm << endl;
  //cout << ne.nu*hc_mev_fm << " " << pr.nu*hc_mev_fm << endl;
  ne.n+=1.0e-5;
  retx=calc_temp_e(ne,pr,0.0,th);
  cout << (th.ed-ed1)/1.0e-5*hc_mev_fm << " ";
  ne.n-=1.0e-5;
  pr.n+=1.0e-5;  
  retx=calc_temp_e(ne,pr,0.0,th);
  cout << (th.ed-ed1)/1.0e-5*hc_mev_fm << endl;
  exit(-1);
  return 0;
}

