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
/// \file Couplings.hpp
/// \author lroberts
/// \since Apr 02, 2016
///
/// \brief
///
///

#ifndef COUPLINGS_HPP_
#define COUPLINGS_HPP_

#include <iostream>

namespace nuopac {

  /// Class defining the coupling constants for various weak interactions, which
  /// is used by the Polarization class. 
  struct WeakCouplings {

    /// Default constructor is charged current couplings without weak magnetism
    WeakCouplings();
  
    /// Constructor for manually specified weak coupling constants 
    WeakCouplings(double cv, double ca, double f2);
  
    /// Return charged current couplings including weak magnetism
    static WeakCouplings NuCapture();

    /// Return neutron scattering couplings including weak magnetism
    static WeakCouplings NeutronScattering();

    /// Return proton scattering couplings including weak magnetism
    static WeakCouplings ProtonScattering();
  
    /// Return nue electron scattering couplings including weak magnetism
    static WeakCouplings NueElectronScattering();
  
    /// Return nux electron scattering couplings including weak magnetism
    static WeakCouplings NuxElectronScattering(); 

    /** \brief Desc
     */
    double Cv;
  
    /** \brief Desc
     */
    double Ca;
  
    /** \brief Desc
     */
    double F2;
  
  };

}
#endif // COUPLINGS_HPP_
