/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2020, Xingfu Du and Andrew W. Steiner
  
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
#ifndef VIRIAL_SOLVER_DERIV_H
#define VIRIAL_SOLVER_DERIV_H

#include <cmath>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <o2scl/test_mgr.h>
#include <o2scl/mm_funct.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/mroot_cern.h>
#include <o2scl/linear_solver.h>
#include <o2scl/poly.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

/** \brief Virial solver with derivatives
 */
class virial_solver_deriv {
  
 protected:
  
  // Generic polynomial solver
  o2scl::poly_real_coeff_gsl<> quart;
  
  /// Storage for the four roots
  std::complex<double> res[4];

 public:

  /// \name First derivatives of the fugacities
  //@{
  double dzndnn;
  double dzndnp;
  double dzpdnn;
  double dzpdnp;
  double dzndT;
  double dzpdT;
  //@}
  
  /// \name Second derivatives of the fugacities
  //@{
  double d2zndnn2;
  double d2zndnndnp;
  double d2zndnp2;
  double d2zpdnn2;
  double d2zpdnndnp;
  double d2zpdnp2;
  //@}

  /// \name Main functions
  //@{
  /** \brief Solve for the fugacity given the density
   */
  virtual void solve_fugacity(double nn, double np,
			      double lam_n, double lam_p,
			      double b_n, double b_pn,
			      double &zn, double &zp) {

    double npt=pow(lam_n,3)/2.0*np;
    double nnt=pow(lam_p,3)/2.0*nn;

    // At high densities or very low densities, just use the
    // non-interacting result
    //
    // AWS: 9/13/2020: I added the "|| nnt>1.0e5 || npt>1.0e5" option
    // later, and this might not be the best method
    // 
    if (nnt<5.0e-6 || npt<5.0e-6 || nnt>1.0e5 || npt>1.0e5) {
      zn=nnt;
      zp=npt;
      return;
    }

    zn=0.0;
    zp=0.0;
    
    double a=pow(b_n,3)*2.0/b_pn/b_pn-2.0*b_n;
    double b=-1+b_n*b_n*2.0/b_pn/b_pn-b_n/b_pn;
    double c=b_n/(2.0*b_pn*b_pn)-0.5/b_pn-nnt+npt-b_n*b_n*npt*2.0/b_pn/b_pn;
    double d=-b_n*npt/b_pn/b_pn+npt/2.0/b_pn;
    double e=b_n*npt*npt/2.0/b_pn/b_pn;
    
    quart.solve_rc(a,b,c,d,e,res[0],res[1],res[2],res[3]);
    
    std::vector<double> zp_list, zn_list;
    for(size_t k=0;k<4;k++) {
      if (res[k].imag()==0.0 && res[k].real()>0.0 && res[k].real()<500.0) {
	double r0, r1;
	gsl_poly_solve_quadratic(2.0*b_n,2.0*res[k].real()*
				 b_pn+1.0,-nnt,&r0,&r1);
	if (r0>0.0 && r0<500.0) {
	  if (r1>0.0 && r1<500.0) {
	    O2SCL_ERR2("Unexpected pair of roots in ",
		       "virial_solver_deriv::solve_fugacity().",
		       o2scl::exc_einval);
	  }
	  std::cout << res[k] << "," << r0 << " ";
	  zp_list.push_back(res[k].real());
	  zn_list.push_back(r0);
	}
	if (r1>0.0 && r1<500.0) {
	  zp_list.push_back(res[k].real());
	  zn_list.push_back(r1);
	}
      }
    }
    if (zp_list.size()==1) {
      zp=zp_list[0];
      zn=zn_list[0];
    } else if (zp_list.size()==2) {
      double norm_0=zp_list[0]*zp_list[0]+zn_list[0]*zn_list[0];
      double norm_1=zp_list[1]*zp_list[1]+zn_list[1]*zn_list[1];
      if (norm_0<norm_1) {
	zp=zp_list[0];
	zn=zn_list[0];
      } else {
	zp=zp_list[1];
	zn=zn_list[1];
      }
    } else {
      std::cout << "virial_solver_deriv::solve_fugacity "
		<< "multiplicity problem:\n\t"
		<< "res0,res1: " << res[0] << " " << res[1] << "\n\t"
		<< "res2,res3: " << res[2] << " " << res[3] << std::endl;
      std::cout << "\tnn,np,lam_n,lam_p: " << nn << " " << np << " "
		<< lam_n << " " << lam_p << "\n\t"
		<< "nnt,npt,zp_list.size(): " << nnt << " " << npt << " "
		<< zp_list.size() << std::endl;
      O2SCL_ERR2("Unexpected root multiplicity in ",
		 "virial_solver_deriv::solve_fugacity().",o2scl::exc_einval);
    }
    return;
  }

  /** \brief Compute the derivatives given the densities and the
      fugacities (i.e. after a call to \ref solve_fugacity())
  */
  virtual void calc_deriv(double nn, double np,
			  double lam_n, double lam_p,
			  double b_n, double b_pn,
			  double zn, double zp,
			  double dbndT, double dbpndT,
			  double dlamndT, double dlampdT) {
    
    double npt=pow(lam_n,3)/2.0*np;
    double nnt=pow(lam_p,3)/2.0*nn;
    
    // At high densities or very low densities, just use the
    // non-interacting result

    if (nnt<5.0e-6 || npt<5.0e-6) {
      
      dzndnn=nnt/nn;
      dzpdnp=npt/np;
      dzndnp=0.0;
      dzpdnn=0.0;
      
      d2zndnn2=0.0;
      d2zndnndnp=0.0;
      d2zndnp2=0.0;
      d2zpdnn2=0.0;
      d2zpdnndnp=0.0;
      d2zpdnp2=0.0;
      dzndT=1.5*lam_n*lam_n*nn*dlamndT;
      dzpdT=1.5*lam_p*lam_p*np*dlampdT;
      
      return;
    }
    
    dzndnn= -pow(lam_n, 3)*(2*b_n*zp + b_pn*zn + 1.0/2.0)/
      (4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*
       (4*b_n*zp + 2*b_pn*zn + 1)) ;
    dzndnp= b_pn*pow(lam_p, 3)*zn/
      (4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*
       (4*b_n*zp + 2*b_pn*zn + 1)) ;
    dzpdnn= b_pn*pow(lam_n, 3)*zp/
      (4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*
       (4*b_n*zp + 2*b_pn*zn + 1)) ;
    dzpdnp= -pow(lam_p, 3)*(2*b_n*zn + b_pn*zp + 1.0/2.0)/
      (4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*
       (4*b_n*zp + 2*b_pn*zn + 1)) ;
    
    double dzndnn_dzn= -b_pn*pow(lam_n, 3)/(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1)) - pow(lam_n, 3)*(2*b_n*zp + b_pn*zn + 1.0/2.0)*(-4*b_n*(-4*b_n*zp - 2*b_pn*zn - 1) - 4*pow(b_pn, 2)*zp + 2*b_pn*(4*b_n*zn + 2*b_pn*zp + 1))/pow(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1), 2) ;
    double dzndnn_dzp= -2*b_n*pow(lam_n, 3)/(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1)) - pow(lam_n, 3)*(2*b_n*zp + b_pn*zn + 1.0/2.0)*(4*b_n*(4*b_n*zn + 2*b_pn*zp + 1) - 4*pow(b_pn, 2)*zn - 2*b_pn*(-4*b_n*zp - 2*b_pn*zn - 1))/pow(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1), 2) ;
    double dzndnp_dzn= b_pn*pow(lam_p, 3)*zn*(-4*b_n*(-4*b_n*zp - 2*b_pn*zn - 1) - 4*pow(b_pn, 2)*zp + 2*b_pn*(4*b_n*zn + 2*b_pn*zp + 1))/pow(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1), 2) + b_pn*pow(lam_p, 3)/(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1)) ;
    double dzndnp_dzp= b_pn*pow(lam_p, 3)*zn*(4*b_n*(4*b_n*zn + 2*b_pn*zp + 1) - 4*pow(b_pn, 2)*zn - 2*b_pn*(-4*b_n*zp - 2*b_pn*zn - 1))/pow(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1), 2) ;
    double dzpdnn_dzn= b_pn*pow(lam_n, 3)*zp*(-4*b_n*(-4*b_n*zp - 2*b_pn*zn - 1) - 4*pow(b_pn, 2)*zp + 2*b_pn*(4*b_n*zn + 2*b_pn*zp + 1))/pow(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1), 2) ;
    double dzpdnn_dzp= b_pn*pow(lam_n, 3)*zp*(4*b_n*(4*b_n*zn + 2*b_pn*zp + 1) - 4*pow(b_pn, 2)*zn - 2*b_pn*(-4*b_n*zp - 2*b_pn*zn - 1))/pow(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1), 2) + b_pn*pow(lam_n, 3)/(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1)) ;
    double dzpdnp_dzn= -2*b_n*pow(lam_p, 3)/(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1)) - pow(lam_p, 3)*(2*b_n*zn + b_pn*zp + 1.0/2.0)*(-4*b_n*(-4*b_n*zp - 2*b_pn*zn - 1) - 4*pow(b_pn, 2)*zp + 2*b_pn*(4*b_n*zn + 2*b_pn*zp + 1))/pow(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1), 2) ;
    double dzpdnp_dzp= -b_pn*pow(lam_p, 3)/(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1)) - pow(lam_p, 3)*(2*b_n*zn + b_pn*zp + 1.0/2.0)*(4*b_n*(4*b_n*zn + 2*b_pn*zp + 1) - 4*pow(b_pn, 2)*zn - 2*b_pn*(-4*b_n*zp - 2*b_pn*zn - 1))/pow(4*pow(b_pn, 2)*zn*zp - (4*b_n*zn + 2*b_pn*zp + 1)*(4*b_n*zp + 2*b_pn*zn + 1), 2) ;

    d2zndnn2=dzndnn_dzn*dzndnn+dzndnn_dzp*dzpdnn;
    d2zndnndnp=dzndnn_dzn*dzndnp+dzndnn_dzp*dzpdnp;
    d2zndnp2=dzndnp_dzn*dzndnp+dzndnp_dzp*dzpdnp;
    d2zpdnn2=dzpdnn_dzn*dzndnn+dzpdnn_dzp*dzpdnn;
    d2zpdnndnp=dzpdnn_dzn*dzndnp+dzpdnn_dzp*dzpdnp;
    d2zpdnp2=dzpdnp_dzn*dzndnp+dzpdnp_dzp*dzpdnp;

    dzndT= (1.0/2.0)*(-16*b_n*dbndT*pow(zn, 2)*zp - 16*b_n*dbpndT*zn*pow(zp, 2) + 12*b_n*dlamndT*pow(lam_n, 2)*nn*zp - 8*b_pn*dbndT*pow(zn, 3) + 8*b_pn*dbndT*zn*pow(zp, 2) + 6*b_pn*dlamndT*pow(lam_n, 2)*nn*zn - 6*b_pn*dlampdT*pow(lam_p, 2)*np*zn - 4*dbndT*pow(zn, 2) - 4*dbpndT*zn*zp + 3*dlamndT*pow(lam_n, 2)*nn)/(16*pow(b_n, 2)*zn*zp + 8*b_n*b_pn*pow(zn, 2) + 8*b_n*b_pn*pow(zp, 2) + 4*b_n*zn + 4*b_n*zp + 2*b_pn*zn + 2*b_pn*zp + 1) ;
    
    dzpdT= (-b_pn*dzndT*zp + (3.0/4.0)*dlamndT*pow(lam_n, 2)*nn - 1.0/2.0*dzndT - zn*(2*b_n*dzndT + dbndT*zn + dbpndT*zp))/(b_pn*zn) ;
    
    return;
  }
  //@}

};

#endif
