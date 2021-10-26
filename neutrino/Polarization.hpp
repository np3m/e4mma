/*
  -------------------------------------------------------------------
  
  This code is based on nuopac, which was developed originally by Luke
  Roberts. Modifications of nuopac are Copyright (C) 2020-2021, Zidu
  Lin and Andrew W. Steiner.

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

#ifndef POLARIZATION_HPP_
#define POLARIZATION_HPP_

#include <iostream>
#include <array> 

#include "FluidState.hpp"
#include "tensor.h"
#include "Couplings.hpp" 

#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_qng_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>
#include <o2scl/inte_qags_gsl.h>

// Going much lower than 64 seems to degrade accuracy at a few percent level
//#define NPGJ 64
//#define NNPGL 128

extern double tempy;

namespace nuopac {

  /// Class for calculating polarization tensors and interaction rates
  /// with a medium at the mean field level.
  class Polarization {
    
  public:
  
    /// Create an instance of the Polarization class - An instance of
    /// the class is constructed with a fixed fluid state and set of
    /// weak couplings. To calculate the cross section for a new set of
    /// thermodynamic conditions, a new instance of the class must be
    /// created.
    Polarization(FluidState st, WeakCouplings wc=WeakCouplings(), 
                 bool antiNeutrino=false, bool doReddy=false,
                 bool doBlock=false,
                 int NAngularPoints=64, int NQ0Points=256); 
  
    /// Set the fluid state to a new value. Don't want to regenerate
    /// Polarization class every time since it calculates nodes and
    /// weights for integration each time it is instantiated.
    void SetFluidState(FluidState stIn) {st=stIn;}

    /// Return the double differential (with the differential in energy
    /// transfer q0 and momentum transfer magnitude q) interaction rate
    /// with the medium as a function of neutrino energy E1
    double CalculateDifGam(double E1, double q0, double q) const;
  
    /// Return the double differential (with the differential in energy
    /// transfer q0 and electron/neutrino scattering angle mu13)
    /// interaction rate with the medium as a function of neutrino
    /// energy E1
    double CalculateDGamDq0Dmu13(double E1, double q0, double mu13) const;

    /// Return the differential interaction rate (with the differential
    /// in energy transfer q0 and electron/neutrino scattering angle
    /// mu13 integrated over) with the medium as a function of neutrino
    /// energy E1
    double CalculateDGamDq0(double E1, double q0);
  
    /// Return Legendre moment of the the differential interaction rate
    /// (with the differential in energy transfer q0 and
    /// electron/neutrino scattering angle mu13 integrated over) with
    /// the medium as a function of neutrino energy E1
    void CalculateDGamDq0l(double E1, double q0, double* S0,
                           double* S1) const;

    /// Calculate the transport inverse mean free path (opacity) for a
    /// neutrino of energy E1 assuming the neutrino distribution
    /// function is given by the equilibrium distribution function of
    /// particle 3
    double CalculateTransportInverseMFP(double E1) const;
  
    /// Calculate the inverse mean free path for a neutrino of energy E1  
    double CalculateInverseMFP(double E1);
  
    /// Calculate the momentum transfer from a given neutrino/electron 
    /// scattering angle 
    inline double GetqFromMu13(double E1, double q0, double mu13) const {
      double p3=sqrt((E1-q0)*(E1-q0)-st.M3*st.M3);
      return sqrt(E1*E1 + p3*p3 - 2.0*E1*p3*mu13);
    } 
  
    /// Print the value of the response function to screen  
    void PrintResponse(double E1, double q0, double q) const;

    /** \brief Desc 
     */
    double GetResponse(double E1, double q0, double q) const;
    
    /** \brief Desc 
     */
    double GetResponse_mu(double E1, double q0, double x, double delta,
                          double avg, std::vector<double> &xv,
                          std::vector<double> &yv);

    /** \brief Desc 
     */
    double dgamdq0_intl(double E1, double x, double estar,
                                      int sign);
    /** \brief Desc 
     */
    void SetCurrentConservation(bool cons){doCurrentConservation=cons;} 
    
    /** \brief Desc 
     */
    std::array<double, 4> CalculateBasePolarizations(double q0,
                                                     double q) const;

    /** \brief Desc 
     */
    inline double GetCsecPrefactor(double E1, double q0) const;

    /** \brief Desc 
     */
    int NPGJ;
  
    /** \brief Desc 
     */
    int NNPGL;

    /** \brief Desc 
     */
    std::vector<double> xx;
    /** \brief Desc 
     */
    std::vector<double> ww;
    /** \brief Desc 
     */
    std::vector<double> xl;
    /** \brief Desc 
     */
    std::vector<double> wl;
    /** \brief Desc 
     */
    std::vector<double> xgl;
    /** \brief Desc 
     */
    std::vector<double> wgl;

    /// \name Flag for vector and axial part
    //@{
    /// The flag
    int flag;
    /// Value corresponding to axial part
    static const int flag_axial=1;
    /// Value corresponding to vector part
    static const int flag_vector=0;
    //@}

    /// current=0 is neutral current, current=1 is charged current
    int current;
    static const int current_charged=1;
    static const int current_neutral=0;

    /// \name Choose integration method
    //@{
    int integ_method_q0, integ_method_mu;
    static const int integ_base=0;
    static const int integ_o2scl=1;
    static const int integ_compare=2;
    //@}

    /// Adaptive integrator with singularities
    o2scl::inte_qags_gsl<> qags;

    /// Adaptive integrator
    o2scl::inte_qag_gsl<> qag;

    /// Nonadaptive integrator
    o2scl::inte_qng_gsl<> qng;

    /// Adaptive integrator with infinite upper limit
    o2scl::inte_qagiu_gsl<> qagiu;
    
  protected:

    /** \brief Desc 
     */
    virtual void SetPolarizations(double q0, double q, 
                                  Tensor<double>* piVV, 
                                  Tensor<double>* piAA, 
                                  Tensor<double>* piTT, 
                                  Tensor<double>* piVA, 
                                  Tensor<double>* piVT, 
                                  Tensor<double>* piAT) const;
  
    /** \brief Desc 
     */
    virtual void SetLeptonTensor(double E1, double q0, double q,
                                 Tensor<double>* L) const;
  
    // double GetResponse(double E1, double q0, double q) const;
    // inline double GetCsecPrefactor(double E1, double q0) const; 
  
    /** \brief Desc 
     */
    WeakCouplings coup;
  
    /** \brief Desc 
     */
    FluidState st;
  
    /** \brief Desc 
     */
    bool antiNeutrino;
  
    /** \brief Desc 
     */
    bool doReddy;
  
    /** \brief Desc 
     */
    bool doBlock;
  
    /** \brief Desc 
     */
    bool doCurrentConservation;
  
    /** \brief Desc 
     */
    double mG2;
  
  };
}
#endif // POLARIZATION_HPP_
