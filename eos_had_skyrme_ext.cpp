/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2022, Xingfu Du, Zidu Lin, and Andrew W. Steiner
  
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

