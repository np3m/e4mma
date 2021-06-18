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
/// \file Polarization.hpp
/// \author lroberts
/// \since Apr 02, 2016
///
/// \brief
///
///

#ifndef POLARIZATIONNONREL_HPP_
#define POLARIZATIONNONREL_HPP_

#include <iostream>
#include <array> 
#include "FluidState.hpp"
#include "Tensor.hpp"
#include "Polarization.hpp" 

namespace nuopac {
///
/// Class for calculating the polarization tensor and neutrino interaction 
/// rates assuming a constant weak interaction matrix element. This inherits 
/// from the full Polarization class and mostly uses methods defined there. It 
/// should behave in exactly the same manner, as all of the changes are only 
/// under the hood. 
/// The name of this class is also somewhat a misnomer, since the response 
/// assumes the nucleons are arbitrarily relativistic once the matrix element 
/// has been fixed.
///
class PolarizationNonRel : public Polarization {
public:
  PolarizationNonRel(FluidState st, WeakCouplings wc = WeakCouplings(), 
      bool antiNeutrino = false, bool doReddy = false, bool doBlock = false) 
      : Polarization(st, wc, antiNeutrino, doReddy, doBlock) {} 
   
  std::array<double, 4> CalculateBasePolarizations(double q0, double q) const;
  std::array<double, 4> CalculateBasePolarizationsNeutron(double q0, double q) const;
  std::array<double, 4> CalculateBasePolarizationsProton(double q0, double q) const;
  double GetImPI( double q0, double q) const;
  double GetImPI2(  double q0, double q) const;
  double GetRePIn( double q0, double q) const;
  double GetRePIp(  double q0, double q) const;
  void PrintGetImPI2(  double q0, double q) const;
  double gamma0(double q0, double q) const;
protected:
  void SetPolarizations(double q0, double q, 
      Tensor<double>* piVV, 
      Tensor<double>* piAA, 
      Tensor<double>* piTT, 
      Tensor<double>* piVA, 
      Tensor<double>* piVT, 
      Tensor<double>* piAT) const;
  void SetLeptonTensor(double E1, double q0, double q, Tensor<double>* L) const;
  
};
}
#endif // POLARIZATIONNONREL_HPP_
