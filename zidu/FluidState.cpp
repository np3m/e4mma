/***********************************************************************
* Copyright (c) 2016 Luke F. Roberts.
* 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to
* the following conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
* LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
* OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
* WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***********************************************************************/
/// \file FluidState.cpp
/// \author lroberts
/// \since Sep 13, 2016
///
/// \brief
///
///

#include <iostream>
#include <array> 
#include <exception> 

#include "FluidState.hpp"
#include "Constants.hpp" 

#ifdef NUOPAC_HAS_GSL 
  #include "FunctionIntegrator.hpp" 
  #include "OneDimensionalRoot.hpp" 
#endif
   
using namespace nuopac;

FluidState::FluidState() : T(0.0), M2(0.0), M3(0.0), M4(0.0), 
        Mu2(0.0), Mu3(0.0), Mu4(0.0), U2(0.0), U4(0.0),  n2(0.0), 
        n3(0.0), n4(0.0) {}

FluidState::FluidState(double T, double M2, double M4, 
    double Mu2, double Mu4, double U2, double U4, double M3, double Mu3) 
    : T(T), M2(M2), M4(M4), Mu2(Mu2), Mu4(Mu4), U2(U2), U4(U4), 
    M3(M3), Mu3(Mu3) {
  n2 = GetDensity(T, M2, Mu2, U2, 2.0);
  n3 = GetDensity(T, M3, Mu3, 0.0, 2.0);
  n3 -= GetDensity(T, M3, -Mu3, 0.0, 2.0);
  n4 = GetDensity(T, M4, Mu4, U4, 2.0);

  // The expression for charged current type rates 
  effectiveDensity = (n2 - n4)/(1.0 - exp((Mu4 - M4 - U4 - Mu2 + M2 + U2)/T));
  // Expression for scattering rates
  if (fabs((n2-n4)/n4) < 1.e-5) 
    effectiveDensity = GetDegenerateElasticReponse(T, M2, Mu2, U2, 2.0); 
}

// Return a state with reversed partine labels and anti-particles for particle 3
FluidState FluidState::ReverseState(FluidState in) {
  return FluidState(in.T, in.M4, in.M2, in.Mu4, in.Mu2, in.U4, in.U2, 
      in.M3, -in.Mu3);
}

// Make a FluidState given two densities  
FluidState FluidState::StateFromDensities(double T, double M2, double M4, 
    double n2, double n4, double U2, double U4, double M3, double n3)  {
  double Mu2 = GetChemicalPotential(T, M2, n2, U2, 2.0);
  double Mu3 = GetChemicalPotential(T, M3, n3, 0.0, 2.0);
  double Mu4 = GetChemicalPotential(T, M4, n4, U4, 2.0); 

  return FluidState(T, M2, M4, Mu2, Mu4, U2, U4, M3, Mu3); 
} 

// Calculate a fluid state in equilibrium with particle species 3   
FluidState  FluidState::BetaEquilibriumConsistentPotential(double T, double M2, 
    double M4, double ntot, double M3, double coup, bool antiParticle) {
  FluidState state = BetaEquilibriumState(T, M2, M4, ntot, 0.0, 0.0, M3, 
      antiParticle, true, coup);
  return state;
}

FluidState FluidState::BetaEquilibriumState(double T, double M2, 
    double M4, double ntot, double U2, double U4, double M3, bool antiParticle,
    bool useCoupling, double coup) {
  
  auto func = [T, M2, M3, M4, U2, U4, ntot, antiParticle, useCoupling, coup] 
      (double y) { 
    double Mu2 = GetChemicalPotential(T, M2, (1-y)*ntot, 0.0, 2.0);  
    double Mu4 = GetChemicalPotential(T, M4, y*ntot, 0.0, 2.0); 
    double deltaU = U2 - U4;
    if (useCoupling) deltaU = coup * ntot * (1.0 - 2.0*y); 
    double Mu3 = Mu2 - Mu4 + deltaU; 
    double n3 = GetDensity(T, M3, Mu3, 0.0, 2.0)
        - GetDensity(T, M3, -Mu3, 0.0, 2.0); 
    if (antiParticle) n3 = -n3; 
    return y - n3/ntot; 
  };
  
#ifdef NUOPAC_HAS_GSL 
  OneDimensionalRoot root1D(tolerance, 100);

  double y = root1D(func, 1.e-10, 1.0-1.e-10); 
  if (useCoupling) {
    U2 = coup * ntot * (0.5 - y);
    U4 =-coup * ntot * (0.5 - y);
  }
  
  double Mu2 = GetChemicalPotential(T, M2, (1-y)*ntot, U2, 2.0);  
  double Mu3 = GetChemicalPotential(T, M3, y*ntot, 0.0, 2.0); 
  double Mu4 = GetChemicalPotential(T, M4, y*ntot, U4, 2.0); 

  return FluidState(T, M2, M4, Mu2, Mu4, U2, U4, M3, Mu3); 
#else 
  std::cerr << std::string("Need GSL to call FluidState::BetaEquilibriumState.") << std::endl;
  throw 0;
  return FluidState(); 
#endif 
}

double FluidState::GetChemicalPotential(double T, double M, double n, 
    double U, double g) {
  
  auto func = [T, M, n, g] (double Mu) { 
    double ng = GetDensity(T, M, Mu, 0.0, g);
    ng -= GetDensity(T, M, -Mu, 0.0, g);
    return 1.0 - ng/n;
  };
  
#ifdef NUOPAC_HAS_GSL 
  OneDimensionalRoot root1D(tolerance, 100);
  double mureldeg = pow(6.0*Constants::Pi*Constants::Pi*fabs(n)/g, 1.0/3.0);
  if (n<0.0) mureldeg = -mureldeg; 
  double muboltz = 
      std::min(M - 1.2*T*log(g/n*pow(M*T/(2.0*Constants::Pi), 1.5)), M);
  muboltz = std::min(muboltz, mureldeg);
  muboltz = 0.0;
  double mudeg = std::max(
      sqrt(pow(3.0*n*2.0*pow(Constants::Pi, 2)/g, 2.0/3.0) + M*M), mureldeg);
  return root1D(func, muboltz, mudeg) + U;  
#else 
  std::cerr << std::string("Need GSL to call FluidState::GetChemicalPotential.") << std::endl;
  throw 0;
  return 0.0; 
#endif
} 

double FluidState::GetDegenerateElasticReponse(double T, double M, 
    double Mu, double U, double g) { 
  
#ifdef NUOPAC_HAS_GSL 
  FunctionIntegrator integrate(512, 0.0, 1.e-1*tolerance); 
  auto integrand = [T, M, Mu, U] (double xx) { 
    double aa = exp(xx + U/T - Mu/T);
    return xx*sqrt(xx*xx-M*M/(T*T))*aa/pow(aa + 1.0, 2) ;
  }; 
  
  double f = 0.0; 
  if (Mu > M + 2.0*T) {
    f = integrate.Integrate(integrand, M/T, Mu/T - 2.0); 
    f += integrate.Integrate(integrand, Mu/T - 2.0, Mu/T + 2.0); 
    f += integrate.Integrate(integrand, Mu/T + 2.0, Mu/T + 30.0); 
  } else { 
    f = integrate.Integrate(integrand, M/T, M/T + 5.0); 
    f += integrate.Integrate(integrand, M/T + 5.0, M/T + 30.0); 
  }

  double nfermi = g/(2.0*Constants::Pi*Constants::Pi) * pow(T, 3) * f;
  
  return nfermi;  
#else 
  std::cerr << std::string("Need GSL to call FluidState::GetDegenerateElasticResponse.") << std::endl;
  throw 0;
  return 0.0; 
#endif
} 
 
double  FluidState::GetDensity(double T, double M, double Mu, double U, 
    double g) { 
  // Make density return fully degenerate density 
  // double nboltz = g*pow(M*T/(2.0*Constants::Pi), 1.5) * exp(Mu/T-M/T) 

#ifdef NUOPAC_HAS_GSL 
  FunctionIntegrator integrate(512, 0.0, 1.e-1*tolerance); 
  auto integrand = [T, M, Mu, U] (double xx) { 
    return xx*sqrt(xx*xx-M*M/(T*T))/(exp(xx + U/T - Mu/T) + 1.0);
  }; 
  
  double f = 0.0; 
  if (Mu/T > M/T + 2.0) {
    f = integrate.Integrate(integrand, M/T, Mu/T - 2.0); 
    f += integrate.Integrate(integrand, Mu/T - 2.0, Mu/T + 2.0); 
    f += integrate.Integrate(integrand, Mu/T + 2.0, Mu/T + 30.0); 
  } else { 
    f = integrate.Integrate(integrand, M/T, M/T + 30.0); 
  }

  double nfermi = g/(2.0*Constants::Pi*Constants::Pi) * pow(T, 3) * f;
  
  return nfermi;  
#else 
  std::cerr << std::string("Need GSL to call FluidState::GetDensity.") << std::endl;
  throw 0; 
  return 0.0; 
#endif
} 

