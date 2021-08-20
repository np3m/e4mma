/// \file FunctionIntegrator.cpp
/// \author jlippuner
/// \since Nov 26, 2014
///
/// \brief
///
///

#include "FunctionIntegrator.hpp"

#ifdef NUOPAC_HAS_GSL 

FunctionIntegrator::FunctionIntegrator(const std::size_t sizeLimit,
                                       const double absError,
                                       const double relError) :
  mSizeLimit(sizeLimit),
  mAbsError(absError),
  mRelError(relError) {
  mpWorkspace = gsl_integration_workspace_alloc(mSizeLimit);
}

FunctionIntegrator::FunctionIntegrator():
  FunctionIntegrator(1024, 0.0, 1.0E-12) {}

FunctionIntegrator::~FunctionIntegrator() {
  gsl_integration_workspace_free(mpWorkspace);
}

#endif
