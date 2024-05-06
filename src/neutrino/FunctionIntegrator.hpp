/// \file FunctionIntegrator.hpp
/// \author jlippuner
/// \since Nov 26, 2014
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_FUNCTIONINTEGRATOR_HPP_
#define SKYNET_UTILITIES_FUNCTIONINTEGRATOR_HPP_

#ifdef NUOPAC_HAS_GSL 

#include <limits>
#include <stdexcept>
#include <string>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

class FunctionIntegrator {
  
public:
  
  FunctionIntegrator(const std::size_t sizeLimit, const double absError,
                     const double relError);

  FunctionIntegrator();

  ~FunctionIntegrator();

  // because this class contains a pointer to allocated memory, we don't allow
  // it to be copied or assigned
  FunctionIntegrator(const FunctionIntegrator&) = delete;
  
  FunctionIntegrator& operator=(const FunctionIntegrator&) = delete;

  template <class FUNCTION>
  
  double Integrate(FUNCTION func, const double lowerLimit,
                   const double upperLimit);

private:
  
  std::size_t mSizeLimit;
  
  double mAbsError;
  
  double mRelError;
  
  gsl_integration_workspace * mpWorkspace;
  
};

class FunctionIntegratorException : public std::runtime_error {
public:
  FunctionIntegratorException(const std::string& what):
    std::runtime_error(what) {}
};

// since this is a template function, the code needs to be in the
// header file, so that the template code is available in the
// compilation unit that calls this function
template<class FUNCTION>
double FunctionIntegrator::Integrate(FUNCTION func, const double lowerLimit,
                                     const double upperLimit) {
  
  if (lowerLimit > upperLimit)
    throw FunctionIntegratorException("Lower limit "
                                      + std::to_string(lowerLimit) +
                                      " is larger than upper limit "
                                      + std::to_string(upperLimit));
  
  if (lowerLimit == upperLimit)
    return 0.0;
  
  // inspired by http://stackoverflow.com/a/13289538/2998298
  gsl_function F;
  F.function = [] (double x, void * p)->double {
                 return (*static_cast<FUNCTION*>(p))(x);
               };
  F.params = &func;

  double result, error;
  int ret = -1;
  
  //gsl_set_error_handler_off(); 

  if (upperLimit == std::numeric_limits<double>::infinity()) {
    if (lowerLimit == -std::numeric_limits<double>::infinity()) {
      ret = gsl_integration_qagi(&F, mAbsError, mRelError, mSizeLimit,
                                 mpWorkspace, &result, &error);
    } else {
      ret = gsl_integration_qagiu(&F, lowerLimit, mAbsError, mRelError,
                                  mSizeLimit, mpWorkspace, &result, &error);
    }
  } else {
    if (lowerLimit == -std::numeric_limits<double>::infinity()) {
      ret = gsl_integration_qagil(&F, upperLimit, mAbsError, mRelError,
                                  mSizeLimit, mpWorkspace, &result, &error);
    } else {
      ret = gsl_integration_qag(&F, lowerLimit, upperLimit, mAbsError,
                                mRelError, mSizeLimit, 6, mpWorkspace,
                                &result, &error);
    }
  }

  if (ret == 0)
    return result;
  else
    throw FunctionIntegratorException("Error occurred in GSL function");
}

#endif

#endif // SKYNET_UTILITIES_FUNCTIONINTEGRATOR_HPP_
