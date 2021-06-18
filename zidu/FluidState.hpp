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
/// \file FluidState.hpp
/// \author lroberts
/// \since Sep 13, 2016
///
/// \brief
///
///

#ifndef FLUIDSTATE_HPP_
#define FLUIDSTATE_HPP_

#include <limits>
#include <iostream>
#include <array> 

//namespace nuopac {
  ///
  /// Structure that contains the properties of 
  ///
  struct FluidState {
    double T, M2, M3, M4, Mu2, Mu3, Mu4, U2, U4;
    double n2, n3, n4, effectiveDensity; 
    static constexpr const double tolerance = 1.e-6; 
    double SensitiveP;//added by zidu, for testing the sensitivity of pol on residual interactions(or any sensitive parameter) 
    /// Default constructor, which will surely break the Polarization class if 
    /// used 
    FluidState();
    
    /// Construct the fluid state by hand  
    FluidState(double T, double M2, double M4, 
        double Mu2, double Mu4, double U2, double U4, double M3, 
        double Mu3 = -std::numeric_limits<double>::infinity());
    
    /// Return a fluid state with reversed particle labels 
    static FluidState ReverseState(FluidState in); 
      
    /// Return a FluidState given the densities of particles 2 and 4 
   // static FluidState StateFromDensities(double T, double M2, double M4, 
     //   double n2, double n4, double U2, double U4, double M3, double n3 = 1.e-10);//original one
    static FluidState StateFromDensities(double T, double M2, double M4,
        double n2, double n4, double U2, double U4, double M3, double n3 = 1.e-10, double SensitiveP =1.e-10);//added by zidu, for testing the sensitivity of pol on residual interactions

  
    static FluidState BetaEquilibriumConsistentPotential(double T, double M2, 
        double M4, double ntot, double M3, double coup, bool antiParticle = false);
    
    /// Return a fluid state with ntot = n2 + n4 in neutrino free beta-equilibrium 
    /// with particle species 3, i.e. mu_2 = mu_3 + mu_4    
    static FluidState BetaEquilibriumState(double T, double M2, double M4, 
        double ntot, double U2, double U4, double M3, bool antiParticle = false,
        bool useCoupling = false, double coup = 0.0);
      
    /// Get the response function appropriete for zero momentum transfer to the 
    /// the nucleons 
    static double GetDegenerateElasticReponse(double T, double M, double Mu, 
        double U, double g = 2.0);
  
    /// Get the chemical potential for a fermion at given density, temperature, 
    /// mass, etc. 
    static double GetChemicalPotential(double T, double M, double n, double U, 
        double g = 2.0);
    
    /// Calculate the density of a particle given its chemical potential   
    static double GetDensity(double T, double M, double Mu, double U, 
        double g = 2.0);
  };
//}

#endif // FLUIDSTATE_HPP_
