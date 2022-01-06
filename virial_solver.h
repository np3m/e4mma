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
#ifndef VIRIAL_SOLVER_H
#define VIRIAL_SOLVER_H

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

/** \brief Compute the virial EOS
 */
class virial_solver {

 public:

  // Store the number of function and derivative evaluations
  int nf, nd;
  double nn, pn, lambda, T, b_n, b_pn;
  double mfn2_mu_n, mfn2_mu_p, dbndT, zp, zn;
  double dbpndT, dlambdadT, npt, nnt, a, b, c, d, e;
  ubmatrix A;
  ubvector B;

  /// Linear system solver
  o2scl_linalg::linear_solver_LU<ubvector,ubmatrix> lsol;  

  /// Quartic polynomial solver
  o2scl::quartic_real_coeff_cern<> quart;

  // Generic polynomial solver
  o2scl::poly_real_coeff_gsl<> quart2;
  
  /// Storage for the four roots
  std::complex<double> res_zp[4],res_zn[4]; 
  /** \brief Solve for the fugacities given the densities

      This function computes zn and zp from pn and nn
      presuming that lambda, b_n and b_pn have already been specified.
  */
  void solve_fugacity(ubvector &x) {

    npt=pow(lambda,3)/2.0*pn;
    nnt=pow(lambda,3)/2.0*nn;

    if (npt>=nnt) {
      
      // Coefficients for quartic equation of zp in descending order
      
      a=pow(b_n,3)*2.0/b_pn/b_pn-2.0*b_n;
      b=-1.0+b_n*b_n*2.0/b_pn/b_pn-b_n/b_pn;
      c=b_n/(2.0*b_pn*b_pn)-0.5/b_pn-nnt+npt-b_n*b_n*npt*2.0/b_pn/b_pn;
      d=-b_n*npt/b_pn/b_pn+npt/2.0/b_pn;
      e=b_n*npt*npt/2.0/b_pn/b_pn;
      
      quart2.solve_rc(a,b,c,d,e,res_zp[0],res_zp[1],res_zp[2],res_zp[3]);
      //std::cout << "Here1: " << a << " " << b << " " << c << " "
      //<< d << " " << e << std::endl;
      int root_count=0;
      std::complex<double> eval_zp[4];
      ubvector res_zn_real;
      std::complex<double>  eval_zn[4];
      res_zn_real.resize(4);
      /*
        AWS: Changed on 1/9/18 to select largest fugacity instead
        of exiting when multiple roots are found
      */

      for (int i=0;i<4;i++) {

        // Check that the root is positive and that the imaginary
        // part is sufficiently small. Note I use fabs() rather
        // than abs().
        if(res_zp[i].real()>0 && 
                        fabs(res_zp[i].imag()/res_zp[i].real())<1.0e-6) {

	  // Make sure that the specified root is really a solution
	  // of the original polynomial
	  eval_zp[i]=(a*pow(res_zp[i],4.0)+b*pow(res_zp[i],3.0)+
		   c*pow(res_zp[i],2.0)+d*res_zp[i]+e)/e;
	   
	  // Changed from 1e-8 to 1e-6 because at zero temperature no
	  // solutions
	   
	  if (fabs(eval_zp[i].real())<2.0e-6 &&
	                                 fabs(eval_zp[i].imag())<1.0e-8) {
            double r0, r1;
            gsl_poly_solve_quadratic(2.0*b_n,2.0*res_zp[i].real()*b_pn+1.0,-nnt,&r0,&r1);
            std::complex<double> eval_r0=(res_zp[i].real()+2.0
			    *res_zp[i].real()*res_zp[i].real()*b_n
			   +2.0*res_zp[i].real()*r0*b_pn-npt)/npt;
            std::complex<double> eval_r1=(res_zp[i].real()+2.0
			    *res_zp[i].real()*res_zp[i].real()*b_n
			   +2.0*res_zp[i].real()*r1*b_pn-npt)/npt;
            if (fabs(eval_r0.real())<2.0e-6 && r0>0.0) {
	      res_zn_real[i]=r0;
              eval_zn[i]=eval_r0;
	    }
            else if (fabs(eval_r1.real())<2.0e-6 && r1>0.0) {
	      res_zn_real[i]=r1;
              eval_zn[i]=eval_r1;
   	    } else {
              res_zn_real[i]=-1.0;
              //std::cout << "r0,r1:  "<< r0 << " " << r1 << std::endl;
              //std::cout << "eval_r0,r1: " << eval_r0.real() << " " 
	      //			  << eval_r1.real() << std::endl;
            }


            	                       
            if (res_zn_real[i]>0.0) { 
	      root_count++; 
	    }
	  } else {
	    res_zn_real[i]=1.0e10;
	  } 
        } else {
            res_zn_real[i]=1.0e10;
        }
      }
      
      if (root_count==0&&false) {
        std::cout << "Zn/Zp zero roots: " << root_count 
                                                      << std::endl;
        std::cout << "nn: " << nn << " pn: " << pn << " T: " 
                  << T << std::endl;
        std::cout.setf(std::ios::showpos);
        std::cout << "zn: " <<res_zn_real[0] << " " << " ";
        std::cout << eval_zn[0].real() << " " << eval_zn[0].imag() 
                  << std::endl;
        std::cout << "zp: " <<res_zp[0].real() << " " << res_zp[0].imag() 
                  << " ";
        std::cout << eval_zp[0].real() << " " << eval_zp[0].imag() 
                  << std::endl;
        std::cout << "zn: " <<res_zn_real[1] << " ";
        std::cout << eval_zn[1].real() << " " << eval_zn[1].imag() 
                  << std::endl;
        std::cout << "zp: " <<res_zp[1].real() << " " << res_zp[1].imag() 
                  << " ";
        std::cout << eval_zp[1].real() << " " << eval_zp[1].imag() 
                  << std::endl;
        std::cout << "zn: " <<res_zn_real[2] << " ";
        std::cout << eval_zn[2].real() << " " << eval_zn[2].imag() 
                  << std::endl;
        std::cout << "zp: " <<res_zp[2].real() << " " << res_zp[2].imag() 
                  << " ";
        std::cout << eval_zp[2].real() << " " << eval_zp[2].imag() 
                  << std::endl;
        std::cout << "zn: " <<res_zn_real[3] << " ";
        std::cout << eval_zn[3].real() << " " << eval_zn[3].imag() 
                  << std::endl;    
        std::cout << "zp: " <<res_zp[3].real() << " " << res_zp[3].imag() 
                  << " ";
        std::cout << eval_zp[3].real() << " " << eval_zp[3].imag() 
                  << std::endl;          
        std::cout.unsetf(std::ios::showpos);
        O2SCL_ERR("Zero or more than one root in solve_fugacity().",
		  o2scl::exc_efailed);
      }
      int res_index=0;
      double minsq=1.0e100, temp;
      for (int i=0;i<4;i++) {
        if (res_zn_real[i]>0.0 && res_zp[i].real()>0.0) {
          temp = res_zn_real[i]*res_zn_real[i] + res_zp[i].real()
                        *res_zp[i].real();
          if (temp<minsq) {
            minsq=temp;
            res_index=i;
          }
        }
      }
      zn=res_zn_real[res_index];
      zp=res_zp[res_index].real();  
      x[1]=log(zp)*T;
      x[0]=log(zn)*T; 

    } else {
    
     
      // Coefficients for quartic equation of zp in descending order
      
      a=pow(b_n,3)*2.0/b_pn/b_pn-2.0*b_n;
      b=-1.0+b_n*b_n*2.0/b_pn/b_pn-b_n/b_pn;
      c=b_n/(2.0*b_pn*b_pn)-0.5/b_pn-npt+nnt-b_n*b_n*nnt*2.0/b_pn/b_pn;
      d=-b_n*nnt/b_pn/b_pn+nnt/2.0/b_pn;
      e=b_n*nnt*nnt/2.0/b_pn/b_pn;
      
      quart2.solve_rc(a,b,c,d,e,res_zn[0],res_zn[1],res_zn[2],res_zn[3]);
      /*
	std::cout << "Here2: " << a << " " << b << " " << c << " "
	<< d << " " << e << std::endl;
	std::cout << res_zn[0] << " " << res_zn[1] << " " << res_zn[2] << " "
	<< res_zn[3] << std::endl;
      */
      int root_count=0;
      std::complex<double> eval_zn[4];
      ubvector res_zp_real;
      std::complex<double> eval_zp[4];
      res_zp_real.resize(4);
      /*
        AWS: Changed on 1/9/18 to select largest fugacity instead
        of exiting when multiple roots are found
      */

      for (int i=0;i<4;i++) {

        // Check that the root is positive and that the imaginary
        // part is sufficiently small. Note I use fabs() rather
        // than abs().
        if(res_zn[i].real()>0 && 
                        fabs(res_zn[i].imag()/res_zn[i].real())<1.0e-6) {

	  // Make sure that the specified root is really a solution
	  // of the original polynomial
	  eval_zn[i]=(a*pow(res_zn[i],4.0)+b*pow(res_zn[i],3.0)+
		   c*pow(res_zn[i],2.0)+d*res_zn[i]+e)/e;
	   
	  // Changed from 1e-8 to 1e-6 because at zero temperature no
	  // solutions
	   
	  if (fabs(eval_zn[i].real())<2.0e-6 &&
	                                 fabs(eval_zn[i].imag())<1.0e-8) {

            double r0, r1;
            gsl_poly_solve_quadratic(2.0*b_n,2.0*res_zn[i].real()*b_pn+1.0,-npt,&r0,&r1);
            std::complex<double> eval_r0=(res_zn[i].real()+2.0
			   *res_zn[i].real()*res_zn[i].real()*b_n
			   +2.0*res_zn[i].real()*r0*b_pn-nnt)/nnt;
            std::complex<double> eval_r1=(res_zn[i].real()+2.0
			     *res_zn[i].real()*res_zn[i].real()*b_n
			   +2.0*res_zn[i].real()*r1*b_pn-nnt)/nnt;
            if (fabs(eval_r0.real())<2.0e-6 && r0>0.0) {
	      res_zp_real[i]=r0;
              eval_zp[i]=eval_r0;
	    } else if (fabs(eval_r1.real())<2.0e-6 && r1>0.0) {
	      res_zp_real[i]=r1;
	      eval_zp[i]=eval_r1;
	    } else {
	      res_zp_real[i]=-1.0;
	      //std::cout << "r0,r1:  "<< r0 << " " << r1 << std::endl;
              //std::cout << "eval_r0,r1: " << eval_r0.real() << " "
              //                            << eval_r1.real() << std::endl;
            }
  
            if (res_zp_real[i]>0.0) { 
	      root_count++; 
	    }
	  } else {
	    res_zp_real[i]=1.0e10;
	  } 
        } else {
            res_zp_real[i]=1.0e10;
        }
      }
      
      if (root_count==0&&false) {
        std::cout << "Zn/Zp zero roots: " << root_count 
                                                      << std::endl;
        std::cout << "nn: " << nn << " pn: " << pn << " T: " 
                  << T << std::endl;
        std::cout.setf(std::ios::showpos);
        std::cout << "zp: " <<res_zp_real[0] << " ";
        std::cout << eval_zp[0].real() << " " << eval_zp[0].imag() 
                  << std::endl;
        std::cout << "zn: " <<res_zn[0].real() << " " << res_zn[0].imag() 
                  << " ";
        std::cout << eval_zn[0].real() << " " << eval_zn[0].imag() 
                  << std::endl;
        std::cout << "zp: " <<res_zp_real[1] << " ";
        std::cout << eval_zp[1].real() << " " << eval_zp[1].imag() 
                  << std::endl;
        std::cout << "zn: " <<res_zn[1].real() << " " << res_zn[1].imag() 
                  << " ";
        std::cout << eval_zn[1].real() << " " << eval_zn[1].imag() 
                  << std::endl;
        std::cout << "zp: " <<res_zp_real[2] << " ";
        std::cout << eval_zp[2].real() << " " << eval_zp[2].imag() 
                  << std::endl;
        std::cout << "zn: " <<res_zn[2].real() << " " << res_zn[2].imag() 
                  << " ";
        std::cout << eval_zn[2].real() << " " << eval_zn[2].imag() 
                  << std::endl;
        std::cout << "zp: " <<res_zp_real[3] << " ";
        std::cout << eval_zp[3].real() << " " << eval_zp[3].imag() 
                  << std::endl;       
        std::cout << "zn: " <<res_zn[3].real() << " " << res_zn[3].imag() 
                  << " ";
        std::cout << eval_zn[3].real() << " " << eval_zn[3].imag() 
                  << std::endl;          
        std::cout.unsetf(std::ios::showpos);
        O2SCL_ERR("Zero or more than one root in solve_fugacity().",
		  o2scl::exc_efailed);
      }
      int res_index=0;
      double minsq=1.0e100, temp;
      for (int i=0;i<4;i++) {
        if (res_zp_real[i]>0.0 && res_zn[i].real()>0.0) {
          temp = res_zp_real[i]*res_zp_real[i] + res_zn[i].real()
                        *res_zn[i].real();
          if (temp<minsq) {
            minsq=temp;
            res_index=i;
          }
        }
      }
      zp=res_zp_real[res_index];
      zn=res_zn[res_index].real();  
      x[1]=log(zp)*T;
      x[0]=log(zn)*T; 
    
    
    
    
    
    
    }

    //std::cout << "zn,zp: " << zn << " " << zp << std::endl;
    
    return;
  }
   
  void mfn_e(ubvector &x) {
    npt=pow(lambda,3)/2*pn;
    //nnt=pow(lambda,3)/2*nn; pn equals nn;
    //mu_p is x[1] while mu_n is x[0]
    a=2*b_n+2*b_pn;
    b=1;
    c=-npt;
    d=(sqrt(b*b+4*a*c)-b)/2;
    x[0]=log(d)*T;
    zn=d;
    zp=d;
    x[1]=x[0];
    /*std::cout<<"------------------"<<std::endl;
      std::cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<log(d)<<" "<<T<<std::endl;
      std::cout<<x[1]<<std::endl;
      std::cout<<"------------------"<<std::endl;*/
    return;
  }

  /** \brief Here, a brief description of this function
      derivative with respect to nn of mfn (linear solver)
  */
  double mfn21(ubvector &x2) {
    zn=exp(mfn2_mu_n/T);
    zp=exp(mfn2_mu_p/T);
    A.resize(2,2);
    B.resize(2);
    A(0,0)=2/pow(lambda,3)*(zn/T+4*zn*zn*b_n/T+2*zp*zn*b_pn/T);
    A(0,1)=2/pow(lambda,3)*2*zp*zn*b_pn/T;
    A(1,0)=2/pow(lambda,3)*2*zp*zn*b_pn/T;
    A(1,1)=2/pow(lambda,3)*(zp/T+4*zp*zp*b_n/T+2*zp*zn*b_pn/T);
    /*std::cout<<"mfn21 start: "<<std::endl;
      std::cout<<mfn2_mu_n<<" "<<mfn2_mu_p<<std::endl;
      std::cout<<zn<<" "<<zp<<" "<<T<<std::endl;
      std::cout<<"mfn21 end: "<<std::endl;*/
    /*std::cout<<A(0,0)<<" "<<A(0,1)<<" "<<A(1,0)<<" "
      <<A(1,1)<<std::endl;*/
    B(0)=1;
    B(1)=0;
    if (true) {
      double den=A(0,0)*A(1,1)-A(0,1)*A(1,0);
      x2[0]=A(1,1)/den;
      x2[1]=-A(1,0)/den;
    } else {
      lsol.solve(2,A,B,x2);
    }
    return 0;  
     
  }

  /** \brief Here, a brief description of this function
      derivative with respect to pn of mfn (linear solver)
  */
  double mfn31(ubvector &x3) {
    zn=exp(mfn2_mu_n/T);
    zp=exp(mfn2_mu_p/T);
    A.resize(2,2);
    B.resize(2);
    A(0,0)=2/pow(lambda,3)*(zn/T+4*zn*zn*b_n/T+2*zp*zn*b_pn/T);
    A(0,1)=2/pow(lambda,3)*2*zp*zn*b_pn/T;
    A(1,0)=2/pow(lambda,3)*2*zp*zn*b_pn/T;
    A(1,1)=2/pow(lambda,3)*(zp/T+4*zp*zp*b_n/T+2*zp*zn*b_pn/T);
    B(0)=0;
    B(1)=1;
    if (true) {
      double den=A(0,0)*A(1,1)-A(0,1)*A(1,0);
      x3[0]=-A(0,1)/den;
      x3[1]=A(0,0)/den;
    } else {
      lsol.solve(2,A,B,x3);
    }

    return 0;  
  }
   
  /** \brief Here, a brief description of this function
      derivative with respect to T of mfn
  */
  int mfn41(ubvector &x4) {
    // dmundT=x4[0];
    // dmupdT=x4[1];

    A(0,0)=2/pow(lambda,3)*(zn/T+4*zn*zn*b_n/T+2*zp*zn*b_pn/T);
    A(0,1)=2/pow(lambda,3)*2*zp*zn*b_pn/T;
    A(1,0)=2/pow(lambda,3)*2*zp*zn*b_pn/T;
    A(1,1)=2/pow(lambda,3)*(zp/T+4*zp*zp*b_n/T+2*zp*zn*b_pn/T);
    B(0)=-(2/pow(lambda,4)*(-3)*dlambdadT*
	   (zn+2*zn*zn*b_n+
	    2*zp*zn*b_pn)+2/pow(lambda,3)*
	   (zn*(-mfn2_mu_n/T/T)+2*zn*zn*dbndT-
	    4*zn*zn*b_n*mfn2_mu_n/T/T+
	    2*zp*zn*dbpndT-
	    2*zp*zn*b_pn*mfn2_mu_p/T/T-
	    2*zp*zn*b_pn*mfn2_mu_n/T/T));
    B(1)=-(2/pow(lambda,4)*(-3)*dlambdadT*
	   (zp+2*zp*zp*b_n+
	    2*zp*zn*b_pn)+2/pow(lambda,3)*
	   (zp*(-mfn2_mu_p/T/T)+2*zp*zp*dbndT-
	    4*zp*zp*b_n*mfn2_mu_p/T/T+
	    2*zp*zn*dbpndT-
	    2*zp*zn*b_pn*mfn2_mu_p/T/T-
	    2*zp*zn*b_pn*mfn2_mu_n/T/T));
    nf++;
    if (true) {
      double den=A(0,0)*A(1,1)-A(0,1)*A(1,0);
      x4[0]=(A(1,1)*B(0)-A(0,1)*B(1))/den;
      x4[1]=(A(0,0)*B(1)-A(1,0)*B(0))/den;
    } else {
      lsol.solve(2,A,B,x4);
    }
    return 0;
  }

};

#endif
