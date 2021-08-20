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
/// \file Couplings.cpp
/// \author lroberts
/// \since Apr 02, 2016
///
/// \brief
///
///


#include <iostream>

#include "Couplings.hpp" 
#include "Constants.hpp" 

using namespace nuopac;

WeakCouplings::WeakCouplings() : 
  Cv(-1.0), Ca(-1.23), F2(0.0) {
}

WeakCouplings::WeakCouplings(double cv, double ca, double f2) : 
  Cv(cv), Ca(ca), F2(f2) {
}

WeakCouplings WeakCouplings::NuCapture() {
  return WeakCouplings(1.0, (Constants::DWeak + Constants::FWeak), 
                       3.706);
}

WeakCouplings WeakCouplings::NeutronScattering() {
  return WeakCouplings(-0.5, (-Constants::DWeak - Constants::FWeak)*0.5, 
                       -0.969194*0.5);
}

WeakCouplings WeakCouplings::ProtonScattering() {
  return WeakCouplings((1.0 - 4.0*Constants::sin2ThetaW)*0.5, 
                       (Constants::DWeak + Constants::FWeak)*0.5, 1.019*0.5);
}

WeakCouplings WeakCouplings::NueElectronScattering() {
  return WeakCouplings(0.5 + 2.0*Constants::sin2ThetaW, 0.5, -0.0);
}

WeakCouplings WeakCouplings::NuxElectronScattering() {
  return WeakCouplings(-0.5 + 2.0*Constants::sin2ThetaW, -0.5, -0.0);
}

