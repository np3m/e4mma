/*
  -------------------------------------------------------------------
  
  This code is based on nuopac, which was developed originally by Luke
  Roberts. Modifications of nuopac are Copyright (C) 2020-2022, Zidu
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
#include "constants.h"

#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_qng_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>
#include <o2scl/inte_qags_gsl.h>
#include <o2scl/inte_adapt_cern.h>
#include "../inte_custom.h"

// Going much lower than 64 seems to degrade accuracy at a few percent level
//#define NPGJ 64
//#define NNPGL 128

extern double tempy;
extern bool integral_debug;

namespace nuopac {

  // Exact expression
  template<class fp_t> fp_t Fermi0(fp_t eta) {
    if (eta>150.0) return eta;  
    return log(exp(eta) + 1.0);
  }
  
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
    double CalculateDGamDq0Dmu13(double E1, double q0, double mu13,
                                 bool pnm=false) const;

    /// Return the differential interaction rate (with the differential
    /// in energy transfer q0 and electron/neutrino scattering angle
    /// mu13 integrated over) with the medium as a function of neutrino
    /// energy E1
    double CalculateDGamDq0(double E1, double q0, bool pnm=false);
  
    /// Return Legendre moment of the the differential interaction rate
    /// (with the differential in energy transfer q0 and
    /// electron/neutrino scattering angle mu13 integrated over) with
    /// the medium as a function of neutrino energy E1
    void CalculateDGamDq0l(double E1, double q0, double* S0,
                           double* S1, bool pnm=false);

    /// Calculate the transport inverse mean free path (opacity) for a
    /// neutrino of energy E1 assuming the neutrino distribution
    /// function is given by the equilibrium distribution function of
    /// particle 3
    double CalculateTransportInverseMFP(double E1, bool pnm=false);
  
    /// Calculate the inverse mean free path for a neutrino of energy E1  
    double CalculateInverseMFP(double E1, bool pnm=false);
  
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
    double GetResponse(double E1, double q0, double q, bool pnm=false) const;
    
    /** \brief Desc 
     */
    double GetResponse_mu(double E1, double q0, double x, double delta,
                          double avg, bool pnm=false);
    
    /** \brief Desc 
     */
    double GetResponse_mu1(double E1, double q0, double x, double delta,
                          double avg, bool pnm=false);
    /** \brief Desc
     */
    double dgamdq0_intl(double E1, double x, double estar,
                        int sign, bool pnm=false);
    /** \brief Desc 
     */
    void SetCurrentConservation(bool cons) {
      doCurrentConservation=cons;
      return;
    } 
    
    /** \brief Desc 
     */
    std::array<double, 4> CalculateBasePolarizations
    (double q0, double q) const;
    
    /** \brief Desc 
     */
    template<class fp_t=long double>
    std::array<double, 4> base_polarization_new
    (double q0, double q) const {
      
      // Calculate some kinematic factors
      
      fp_t q0t = q0 + st.U2 - st.U4;//orig one
      
      fp_t qa2t = q0t*q0t - q*q;
      
      // [LR]: I don't completely understand this condition, but it seems to be
      // necessary to suppress noise at larger q_0
      
      //orig one
      if (qa2t > pow(st.M2 - st.M4, 2)*0.0) return {0.0, 0.0, 0.0, 0.0};
      
      if (qa2t < 1.e-1 && qa2t >= 0.0) return {0.0, 0.0, 0.0, 0.0}; 
      fp_t beta = 1.0 + (st.M2*st.M2 - st.M4*st.M4)/qa2t;//orig
      
     fp_t arg = beta*beta - 4.0*st.M2*st.M2/qa2t;
      if (arg<0.0) return {0.0, 0.0, 0.0, 0.0};
      fp_t em;
      if (-0.5*beta*q0t + 0.5*q*sqrt(arg)>st.M2) {
        em=-0.5*beta*q0t + 0.5*q*sqrt(arg);
      } else {
        em=st.M2;
      }
      //fp_t em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), st.M2);
      fp_t delta2 = (st.Mu2 - st.U2 - em)/st.T;
      fp_t delta4 = (st.Mu4 - st.U4 - em - q0t)/st.T;//orig one
      
      // following is the new em completely consistent with Redddy's
      // thesis, for non-rel+interacting gas
      fp_t chi=1-st.M4/st.M2; //orig
      
      fp_t c=q0+st.U2-st.U4-q*q/(2*st.M4);
      //the minimum E2 for NC reaction
      //fp_t emNC=std::max((-c*st.M2/q)*(-c*st.M2/q)/(2*st.M2),0.0); 
      fp_t emNC=(-c*st.M2/q)*(-c*st.M2/q)/(2*st.M2);
      if (emNC>0.0) emNC=0.0;
          
      fp_t argCC=1+2*chi*st.M4*c/(q*q);
      if (argCC<0.0) return {0.0, 0.0, 0.0, 0.0};
      //the minimum E2 for CC reaction
      fp_t eminCC=2*q*q/(chi*chi)*(1+chi*st.M4*c/
                                          (q*q)-sqrt(argCC))/(2*st.M2);
      //the maximum E2 for CC reaction
      fp_t emaxCC=2*q*q/(chi*chi)*(1+chi*st.M4*c/
                                          (q*q)+sqrt(argCC))/(2*st.M2); 
      
      fp_t delta2NC=(st.Mu2-st.U2-emNC)/st.T;
      fp_t delta4NC=(st.Mu4-st.U4-emNC-q0t)/st.T;
      
      fp_t delta2minCC=(st.Mu2-st.U2-eminCC)/st.T;
      fp_t delta4minCC=(st.Mu4-st.U4-eminCC-q0t)/st.T;
      
      fp_t delta2maxCC=(st.Mu2-st.U2-emaxCC)/st.T;
      fp_t delta4maxCC=(st.Mu4-st.U4-emaxCC-q0t)/st.T;
      
      // [LR]: Now just need to include some method for calculating these
      // At least at low density, Gamma0 should be the dominant term
      // which looks like the non-relativistic response function 
      // Under non-degenerate conditions (i.e. delta2, delta4 << 0), 
      // Gamma0 = Gamma1 = 0.5*Gamma2 
      // This is exact 
      fp_t Gamma0 = Fermi0(delta2) - Fermi0(delta4);
      
      // reddy nc current
      fp_t xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);
      //Reddy cc current
      fp_t ximinCC=Fermi0(-delta2minCC) - Fermi0(-delta4minCC);
      //Reddy cc current
      fp_t ximaxCC=Fermi0(-delta2maxCC) - Fermi0(-delta4maxCC); 
      
      //if (current==current_neutral) {
      //O2SCL_ERR("Invalid current 1.",o2scl::exc_efailed);
      //}
      
      double PI = o2scl_const::pi;//Constants::Pi; 
      // double piL = st.M2*st.M4*st.T/(PI*q)*Gamma0;//orig one
      
      //neutral current consistent with Reddy's thesis
      // double piL= st.M2*st.M4*st.T/(PI*q)*(xiNC+q0t/st.T);
      
      //charged current consistent with Reddy's thesis
      double piL= ((double)(st.M2*st.M4*st.T/(PI*q)*(ximinCC-ximaxCC)));
      
      tempy=piL;
      
      if (integral_debug) {
        std::cout << "XX: " << q << " " << q0 << " " << st.M2 << " " << st.M4
                  << " "
                  << em << " "
                  << delta2minCC << " " << ximinCC-ximaxCC << " "
                  << tempy << std::endl;
      }
      
      fp_t piQ = 0.0; 
      fp_t piM = 0.0;
      fp_t piT = 0.0;
      
      return {(double)piQ, (double)piL, (double)piM, (double)piT};
    }

    /** \brief Desc 
     */
    template<class fp_t> fp_t GetCsecPrefactor(fp_t E1, fp_t q0) const {
      fp_t p3=sqrt((E1-q0)*(E1-q0)-st.M3*st.M3); 
      const fp_t fac=mG2*pow(Constants::GfMeV, 2)/
        pow(2.0*o2scl_const::pi,3)/16.0;
      fp_t a;
      if (current==1) {
	if (abs(st.Mu4-st.Mu2-q0)<1.0e-15){a=E1*1.0e-15/st.T; std::cout<<"Be Carefule! q0 meets mu4-mu2. it is close to sinqularity!"<<std::endl;}
	else {
        a=E1*(1.0-exp((st.Mu4-st.Mu2-q0)/st.T));
	}
      } else {
	if (abs(q0)<1.0e-15){a=E1*1.0e-15/st.T; std::cout<<"Be Carefule! q0 meets 0.0. it is close to sinqularity!"<<std::endl;}
        else {
        a=E1*(1.0-exp((-q0)/st.T));
	}
      }
      if (doBlock) p3 *= 1.0 - 1.0/(exp((E1-q0-st.Mu3)/st.T) + 1.0);
      return fac*p3/a;
    }

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
    static const int integ_cubature=3;
    static const int integ_mc=4;
    //@}

    /// Adaptive integrator with singularities
    o2scl::inte_qags_gsl<> qags;
    
    // Desc
   // o2scl::inte_adapt_cern<> iac;

    /// Adaptive integrator
    o2scl::inte_qag_gsl<> qag;

    /// Nonadaptive integrator
    o2scl::inte_qng_gsl<> qng;

    /// Adaptive integrator with infinite upper limit
    o2scl::inte_qagiu_gsl<> qagiu;

    /// Desc
    o2scl::inte_custom<> ic;
    
  protected:

    /** \brief Desc 
     */
    virtual void SetPolarizations(double q0, double q, 
                                  Tensor<double>* piVV, 
                                  Tensor<double>* piAA, 
                                  Tensor<double>* piTT, 
                                  Tensor<double>* piVA, 
                                  Tensor<double>* piVT, 
                                  Tensor<double>* piAT, bool pnm=false) const;
  
    /** \brief Desc 
     */
    virtual void SetLeptonTensor(double E1, double q0, double q,
                                 Tensor<double>* L) const;
  
    // double GetResponse(double E1, double q0, double q) const;
    // inline double GetCsecPrefactor(double E1, double q0) const; 
  
    /** \brief Desc 
     */
    WeakCouplings coup;

  public:
    
    /** \brief Desc 
     */
    FluidState st;

  protected:
  
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

    /** \brief Desc 
     */
    double integrand_mc(size_t ndim,
                        const boost::numeric::ublas::vector<double> &xi,
                        void *fdata);

    double integrand_mc1(size_t ndim,
                        const boost::numeric::ublas::vector<double> &xi,
                        void *fdata);

    
  };

  typedef struct integration_params_s {
    Polarization *p;
    double estar;
    int sign;
    double E1;
    std::vector<double> *xv;
    std::vector<double> *yv;
    bool pnm;
  } integration_params;
  
  int integrand_new(unsigned ndim, const double *x, void *fdata,
                    unsigned fdim, double *fval);

}
#endif // POLARIZATION_HPP_
