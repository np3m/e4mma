/// \file OneDimensionalRoot.hpp
/// \authorr lroberts
/// \since Aug 18, 2015
///
/// \brief
///
///

#ifndef EOS_ONEDIMENSIONALROOT_HPP_
#define EOS_ONEDIMENSIONALROOT_HPP_
#ifdef NUOPAC_HAS_GSL

#include <functional> 
#include <stdexcept> 
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h> 
#include <iostream> 
#include <string>
#include <sstream> 

class OneDimensionalRoot { 
public:
  OneDimensionalRoot(double tol=1.e-8, int maxIter=50) :
      mTol(tol), mMaxIter(maxIter) {
    T = gsl_root_fsolver_brent; 
    s = gsl_root_fsolver_alloc(T); 
    gsl_set_error_handler_off();
  } 
  
  ~OneDimensionalRoot() { gsl_root_fsolver_free(s);}

  template<class FUNCTION>
  double operator() (FUNCTION func, const double xlo, const double xhi) {
    
    double flo = func(xlo); 
    double fhi = func(xhi); 

    if (fabs(flo)<mTol) return xlo;
    if (fabs(fhi)<mTol) return xhi;

    if (flo*fhi>0) {
      std::stringstream stout; 
      stout << "Bad interval in one-d root find " << xlo << " " 
          << xhi << " " << flo << " " << fhi << std::endl;
      throw std::runtime_error(stout.str());
    }
    
    
    
    // Similar to stack overflow 13289311 via Jonas 
    gsl_function F; 
    F.function = [] (double x, void * p)->double { 
      return (*static_cast<FUNCTION*>(p))(x);
    };
    F.params = &func;
    
    
    //const gsl_root_fsolver_type *T; 
    //gsl_root_fsolver *s;
    //
    //T = gsl_root_fsolver_brent; 
    //s = gsl_root_fsolver_alloc(T); 
    gsl_root_fsolver_set(s, &F, xlo, xhi); 
    
    int status, status1, status2; 
    int iter = 0; 
    double x_lo, x_hi;
    x_lo = xlo;
    x_hi = xhi;
    do { 
      iter++; 
      status = gsl_root_fsolver_iterate(s);
      x_lo = gsl_root_fsolver_x_lower(s);  
      x_hi = gsl_root_fsolver_x_upper(s); 
      if (fabs((x_lo - x_hi)/x_lo) < 1.e-13) break;
      status = gsl_root_test_interval(x_lo, x_hi, 0, mTol);  
      status1 = gsl_root_test_residual(func(x_lo), mTol);  
      status2 = gsl_root_test_residual(func(x_hi), mTol);  
    } while ((status == GSL_CONTINUE || 
        (status1 == GSL_CONTINUE && status2 == GSL_CONTINUE)) && iter < mMaxIter);
    
    //gsl_root_fsolver_free(s);    
    //gsl_set_error_handler(NULL);
    
    double xcen; 
    if (fabs(func(x_lo))<fabs(func(x_hi))) {
      xcen = x_lo;
    } else {
      xcen = x_hi;
    }
    if (fabs(func(xcen))>fabs(func(0.5*(x_lo+x_hi)))) xcen = 0.5*(x_lo+x_hi);
    
    if ((status != 0 || iter >= mMaxIter) && (fabs((x_hi-x_lo)/x_lo)>1.e-12 
        || fabs(func(xcen)) > 1.e-2)) {
      std::stringstream stout; 
      stout << "Root find did not converge " << x_lo << " " << x_hi << " " 
          << " " << (x_hi - x_lo)/x_lo << " " 
          << func(x_lo) << " " << func(x_hi) << std::endl;
      throw std::runtime_error(stout.str());
    }
    
    //if (fabs(func(xcen))>mTol*1.e2 && ) { 
    //    std::stringstream stout; 
    //    stout << "Root find did not converge 2 " << x_lo << " " << x_hi 
    //        << " " << (x_hi - x_lo)/x_lo
    //        << " " << func(x_lo) << " " << func(x_hi) << std::endl;
    //    throw std::runtime_error(stout.str());
    //}
    return xcen; 
  
  }

protected: 
  double mTol; 
  int mMaxIter; 
  const gsl_root_fsolver_type *T; 
  gsl_root_fsolver *s;
  
};
#endif
#endif // EOS_ONEDIMENSIONALROOT_HPP_

