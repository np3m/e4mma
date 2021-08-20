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
/// \file Polarization.cpp
/// \author lroberts
/// \since Apr 02, 2016
///
/// \brief
///
///
#include <math.h> 
#include <limits> 

//#include "FunctionIntegrator.hpp" 
//#include "OneDimensionalRoot.hpp" 
#include "jacobi_rule.hpp" 
#include "Polarization.hpp"
#include "tensor.h" 
#include "constants.h"

// Use unnamed namespace so these methods are only locally available
namespace {
  
  // Exact expression
  inline double Fermi0(double eta) {
    if (eta>150.0) return eta;  
    return log(exp(eta) + 1.0);
  }

  // Use the approximations of Takahashi, El Eid, and Hillebrandt 1978
  // These functions have maximum errors of ~0.3%, but are much better away from 
  // delta ~ 0. Additionally, these functions are sub-dominant for large masses 
  // so the error induced in the cross sections will be very small  
  inline double Fermi1(double eta) { 
    if (eta>1.e-3) 
        return (pow(eta, 2)*0.5 + 1.6449)/(1.0 + exp(-1.6855*eta));
    return exp(eta)/(1.0 + 0.2159*exp(0.8857*eta));
  }

  inline double Fermi2(double eta) { 
    if (eta>1.e-3) 
        return (pow(eta, 3)/3.0 + 3.2899*eta)/(1.0 - exp(-1.8246*eta)); 
    return 2.0*exp(eta)/(1.0 + 0.1092*exp(0.8908*eta));
  }
  inline std::array<double, 3> FermiAll(double eta) {
    double expeta = exp(eta); 
    if (eta>1.e-3)  
      return {log(expeta + 1.0), 
          (pow(eta, 2)*0.5 + 1.6449)/(1.0 + exp(-1.6855*eta)),
          (pow(eta, 3)/3.0 + 3.2899*eta)/(1.0 - exp(-1.8246*eta))};
    return {log(expeta + 1.0),
           expeta/(1.0 + 0.2159*exp(0.8857*eta)),
           2.0*expeta/(1.0 + 0.1092*exp(0.8908*eta))};
  }
}

using namespace nuopac; 

Polarization::Polarization(FluidState stPol, WeakCouplings wc,
                           bool antiNeutrinoPol, 
                           bool doReddyPol, bool doBlockPol,
                           int NAngularPoints, int NQ0Points) {

  coup=wc;
  st=stPol;
  antiNeutrino=antiNeutrinoPol;
  doReddy=doReddyPol;
  doBlock=doBlockPol;
  doCurrentConservation=false;
  mG2=2.0;
  NPGJ=NAngularPoints;
  NNPGL=NQ0Points;
  xx.resize(NPGJ);
  ww.resize(NPGJ);
  xl.resize(NPGJ);
  wl.resize(NPGJ);
  xgl.resize(NNPGL);
  wgl.resize(NNPGL);
  
  cgqf(NPGJ, 1,  0.0, 0.0, -1.0, 1.0, xl.data(), wl.data()); 
  cgqf(NPGJ, 4, -0.5, 0.0, -1.0, 1.0, xx.data(), ww.data()); 
  cgqf(NNPGL, 5,  0.0, 0.0, 0.0, 1.0, xgl.data(), wgl.data()); 
  for (int i=0; i<NPGJ; ++i) ww[i] = ww[i] * sqrt(1.0 - xx[i]);
  //for (int i=0; i<NPGJ; ++i) wgl[i] = wgl[i] * exp(xgl[i]);
  flag=0;
  current=0;
}

std::array<double, 4> Polarization::CalculateBasePolarizations (
    double q0, double q) const {
//void Polarization::SetPolarizations(double q0, double q) {
  // Calculate some kinematic factors
  double q0t = q0 + st.U2 - st.U4;
  double qa2t = q0t*q0t - q*q; 
  // I don't completely understand this condition, but it seems to be necessary 
  // to suppress noise at larger q_0
  if (qa2t > pow(st.M2 - st.M4, 2)*0.0) return {0.0, 0.0, 0.0, 0.0}; 
  if (qa2t < 1.e-1 && qa2t > 0.0) return {0.0, 0.0, 0.0, 0.0}; 
  double beta = 1.0 + (st.M2*st.M2 - st.M4*st.M4)/qa2t; 
  double arg = beta*beta - 4.0*st.M2*st.M2/qa2t;
  if (arg<0.0) return {0.0, 0.0, 0.0, 0.0}; 
  double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), st.M2);
  double delta2 = (st.Mu2 - st.U2 - em)/st.T;
  double delta4 = (st.Mu4 - st.U4 - em - q0t)/st.T;

  // Now just need to include some method for calculating these
  // At least at low density, Gamma0 should be the dominant term
  // which looks like the non-relativistic response function 
  // Under non-degenerate conditions (i.e. delta2, delta4 << 0), 
  // Gamma0 = Gamma1 = 0.5*Gamma2 
  // This is exact 
  double Gamma0 = Fermi0(delta2) - Fermi0(delta4);
  double Gamma1 = Fermi1(delta2) - Fermi1(delta4);
  double Gamma2 = Fermi2(delta2) - Fermi2(delta4);
  
  // We expect a >> 1 because em ~ M, which means that terms proportional to 
  // Gamma0 are dominant except for under extreme degeneracy 
  if (doReddy) {
    q0t = q0;
    qa2t = q0*q0 - q*q;
  }
  double a = (beta*q0t + 2.0*em)/st.T;
  double PI = Constants::Pi; 
  double piQ = qa2t*st.T/(4.0*PI*q)*Gamma0; 
  double piL = -qa2t*pow(st.T/q, 3)/(4.0*PI) 
      * (a*a*Gamma0 + 4.0*a*Gamma1 + 4.0*Gamma2);
  double piM = -qa2t*pow(st.T/q, 2)/(4.0*PI)*(a*Gamma0 + 2.0*Gamma1);
  double piT = -0.5*piL + (2.0*st.M2*st.M2/qa2t - 0.5*beta*beta)*piQ;
  
  return {piQ, piL, piM, piT};
}

void Polarization::SetPolarizations(double q0, double q,
      Tensor<double>* piVV, 
      Tensor<double>* piAA, 
      Tensor<double>* piTT, 
      Tensor<double>* piVA, 
      Tensor<double>* piVT, 
      Tensor<double>* piAT) const {
  // Calculate the basic parts of the polarization 
  auto pt = CalculateBasePolarizations(q0, q); 
  double piQ = pt[0]; 
  double piL = pt[1]; 
  double piM = pt[2]; 
  double piT = pt[3]; 
  
  // Calculate various kinematic factors
  double q0t = q0 + st.U2 - st.U4; 
  double qa2 = q0t*q0t - q*q; 
  double qa2i = 1.0/qa2;
  double deltam = st.M2 - st.M4;
  double summ = st.M2 + st.M4; 
  double lambda = (st.M2*st.M2 - st.M4*st.M4)*qa2i;
  if (doCurrentConservation || doReddy) {
    q0t = q0;
    qa2 = q0*q0 - q*q; 
    qa2i = 1.0/qa2; 
    summ = 2.0*st.M2; 
    deltam = 0.0; 
    lambda = 0.0;
  }
  double sigmam = 1.0 - deltam*deltam*qa2i; 
  double sigmap = 1.0 - summ*summ*qa2i; 
  double Delta = -deltam/st.M2; 
  double beta = 1.0 + lambda;
  
  // Set the different parts of the polarization 
  piVV->Q  = (lambda*lambda + sigmam - 1.0)*piQ; 
  piVV->L  = piL + sigmam*piQ; 
  piVV->Tp = piT + sigmam*piQ; 
  piVV->Mp = lambda * piM;
  
  piAA->Q  = (lambda*lambda + sigmap - 1.0)*piQ; 
  piAA->L  = piL + sigmap*piQ; 
  piAA->Tp = piT + sigmap*piQ; 
  piAA->Mp = lambda * piM;

  piTT->L  = qa2/(4.0*st.M2*st.M2)*((sigmam - beta*beta 
      + 4.0*st.M2*st.M2*qa2i)*piQ - piL);
  piTT->Tp = qa2/(4.0*st.M2*st.M2)*((sigmam - beta*beta 
      + 4.0*st.M2*st.M2*qa2i)*piQ - piT);
  
  piVT->L  = (2.0 + Delta*beta) * piQ; 
  piVT->Tp = (2.0 + Delta*beta) * piQ;
  piVT->Mp = -0.5*Delta * piM; 
  
  piVA->Tm = 2.0 * piM;
  
  piAT->Tm = (2.0 + Delta) * piM;
}

void Polarization::SetLeptonTensor(double E1, double q0, double q, Tensor<double>* L) const {
  double qa2 = q0*q0 - q*q;
  double E3 = E1 - q0;
  double DeltaU = st.U2 - st.U4; 
  if (doReddy) DeltaU = 0.0;
  double qa2t = (q0+DeltaU)*(q0+DeltaU) - q*q; 

  double p3 = sqrt(E3*E3 - st.M3*st.M3);
  double p1dotp3 = E1*E3 - 0.5*(E1*E1 + p3*p3 - q*q); 
  double p1dotq = -p1dotp3;
  double p1dotnt = -(qa2*E1 + q0*p1dotp3 + (E1*q0+p1dotp3)*DeltaU)/q; 
  double p1dotn = -(qa2*E1 + q0*p1dotp3)/q; 
  double p1dotqt = p1dotq + DeltaU*E1; 
  double qdotqt = q0*(q0 + DeltaU) - q*q; 
  double qdotnt = -q*DeltaU;
  
  L->Tp = 8.0*(p1dotnt*p1dotnt - p1dotnt*qdotnt + p1dotqt*(qdotqt - p1dotqt))/qa2t; 
  L->L = 8.0*(-2.0*p1dotnt*p1dotnt + 2.0*p1dotnt*qdotnt + qa2t*p1dotq)/qa2t;
  L->Q = 8.0*(2.0*p1dotqt*p1dotqt - 2.0*p1dotqt*qdotqt + qa2t*p1dotq)/qa2t; 
  L->Mp = 8.0*(qdotnt*p1dotqt + p1dotnt*(qdotqt - 2.0*p1dotqt))/qa2t; 
  L->Mm = 0.0; // This is actually non-zero, but there is no contribution here 
              // from the polarization tensor in the mean field limit, need to
              // include the actual value here for RPA
  L->Tm = 8.0*(p1dotnt*qdotqt-qdotnt*p1dotqt)/qa2t; 
  
  if (antiNeutrino) {
    L->Tm = -L->Tm;
    L->Mm = -L->Mm;
  } 
}
  
double Polarization::CalculateDGamDq0Dmu13(double E1, double q0, double mu13) const {
  double q = GetqFromMu13(E1, q0, mu13); 
  return GetResponse(E1, q0, q)*2.0*Constants::Pi*GetCsecPrefactor(E1, q0); 
}
  
double Polarization::CalculateDGamDq0(double E1, double q0) const { 
  
  // Only integrate over angles for which |q0| < q
  double p3 = sqrt((E1-q0)*(E1-q0) - st.M3*st.M3);
  double mu13cross = std::max((E1*E1 + p3*p3 - q0*q0)/(2.0*E1*p3), -1.0);
  double delta = (mu13cross + 1.0) / 2.0; 
  double avg = (mu13cross - 1.0) / 2.0;
   
  double integral = 0.0; 
  for (int i=0; i<NPGJ; ++i) { 
      double mu = xx[i]*delta + avg; 
      integral += ww[i] * GetResponse(E1, q0, GetqFromMu13(E1, q0, mu));
  }
  integral *= delta;
  double fac = GetCsecPrefactor(E1, q0);
  double crx=2.0*Constants::Pi*fac*integral;
  if (crx<0.0) crx=0.0;//added by zidu
  return crx; //added by zidu, at high q0, the crx can be small and negative, which might result from calculation accuracy. a negative crx will result in a "nan" of the code output.
 // return 2.0*Constants::Pi*fac*integral; //orig
}

void Polarization::CalculateDGamDq0l(double E1, double q0, double* S0, 
    double* S1) const { 
  
  // Only integrate over angles for which |q0| < q
  double p3 = sqrt((E1-q0)*(E1-q0) - st.M3*st.M3);
  double mu13cross = std::max((E1*E1 + p3*p3 - q0*q0)/(2.0*E1*p3), -1.0);
  double delta = (mu13cross + 1.0) / 2.0; 
  double avg = (mu13cross - 1.0) / 2.0;
   
  double integral0 = 0.0; 
  double integral1 = 0.0; 
  for (int i=0; i<NPGJ; ++i) { 
      double mu = xx[i]*delta + avg;
      double r = GetResponse(E1, q0, GetqFromMu13(E1, q0, mu));
      integral0 += ww[i] * r;
      integral1 += ww[i] * mu * r;
  }
  integral0 *= delta;
  integral1 *= delta;

  double fac = GetCsecPrefactor(E1, q0);
  *S0 = 2.0*Constants::Pi*fac*integral0/2.0;
  *S1 = 2.0*Constants::Pi*fac*integral1*3.0/2.0;
}

double Polarization::CalculateTransportInverseMFP(double E1) const {
  
  double estar = st.M4 + st.U4 - st.M2 - st.U2; 
  estar = std::min(estar, E1 - st.M3);
  double integral = 0.0;  
  for (int sign = -1; sign<2; sign += 2) { // Integrate on both sides of estar
    for (int i=0; i<NNPGL; ++i) {
      double q0 = estar + double(sign)*xgl[i]*st.T;
      if (q0 >= E1 - st.M3) break; 

      double G0, G1; 
      CalculateDGamDq0l(E1, q0, &G0, &G1); 
      double e0 = log(G0) + xgl[i];
      double e1 = log(fabs(G1)) + xgl[i];
      double sgnG1 = (G1 > 0) - (G1 < 0); 
      double emq0 = exp(-q0/st.T); 
      double a1 = exp((E1 - st.Mu3)/st.T); 
      double a3 = exp((E1 - q0 - st.Mu3)/st.T); 
      // fac0 = 1 - f^0_3(1-e^(-q_0/T))
      // fac1 = f_3^1/f_1^1 e^(-q_0/T)(1 - f^0_1(1-e^(q_0/T)))
      // See Pons et al. (1999) eq. 25 or hera++ notes eq. 62 (Feb 12, 2020)
      double fac0 = (1.0+1.0/a1)/(1.0 + 1.0/(a1*emq0));
      double fac1 = 1.0/fac0; // This assumes ratio of f_1^1/f_3^1 = e^(-q_0/T)
                              // where f_i^l is the lth Legendre moment of the 
                              // distribution function at E_i.
                              // At q0 = 0, fac0 = fac1 = 1, and this clearly 
                              // returns the inelastic transport opacity
      integral += wgl[i] * (fac0*exp(e0) - fac1*sgnG1*exp(e1)/3.0); 
    }
  }
  return 2.0*st.T*integral; 
  
}

double Polarization::CalculateInverseMFP(double E1) const {
  
  double estar;
  if (current==1) {
    estar=st.M4 + st.U4 - st.M2 - st.U2;
  } else {
    estar=0.0;
  }
  estar = std::min(estar, E1 - st.M3);
  double integral = 0.0;  
  for (int i=0; i<NNPGL; ++i) {
    double q0 = estar - xgl[i]*st.T;
    if (abs(q0) > 30.0*st.T) break;
    double ee = log(CalculateDGamDq0(E1, q0)) + xgl[i];
    integral += wgl[i] * exp(ee); 
  }

  // This is maybe not the best idea  
  for (int i=0; i<NNPGL; ++i) {
    double q0 = estar + xgl[i]*st.T; 
    if (q0>E1 - st.M3) break; 
    double ee = log(CalculateDGamDq0(E1, q0)) + xgl[i];
    integral += wgl[i] * exp(ee); 
  }
  
  return st.T*integral; 
  
}
 
void Polarization::PrintResponse(double E1, double q0, double q) const {
  
  Tensor<double> piVV, piAA, piTT, piVA, piVT, piAT;   
  SetPolarizations(q0, q, &piVV, &piAA, &piTT, &piVA, &piVT, &piAT);
  Tensor<double> L;
  SetLeptonTensor(E1, q0, q, &L);  
  std::cout << "VV: " << FullContract(L, piVV) << std::endl;
  std::cout << "AA: " << FullContract(L, piAA) << std::endl;
  std::cout << "TT: " << FullContract(L, piTT) << std::endl;
  std::cout << "VA: " << FullContract(L, piVA) << std::endl;
  std::cout << "VT: " << FullContract(L, piVT) << std::endl;
  std::cout << "AT: " << FullContract(L, piAT) << std::endl;
}
  
  

double Polarization::GetResponse(double E1, double q0, double q) const {
  Tensor<double> L;
  SetLeptonTensor(E1, q0, q, &L);
  Tensor<double> piVV, piAA, piTT, piVA, piVT, piAT;   
  SetPolarizations(q0, q, &piVV, &piAA, &piTT, &piVA, &piVT, &piAT);
  
  double pi = 0.0;
  if (flag==0) {
    pi += coup.Cv*coup.Cv*FullContract(L, piVV);
  } else if (flag==1) {
    pi += coup.Ca*coup.Ca*FullContract(L, piAA);
  }
 // delete the last four terms to produce ccOutputT10E10betaRelnoVAVT 
//  pi += coup.F2*coup.F2*FullContract(L, piTT);
//  pi += coup.Cv*coup.Ca*FullContract(L, piVA);
 // pi += coup.Cv*coup.F2*FullContract(L, piVT);
 // pi += coup.Ca*coup.F2*FullContract(L, piAT);
  return pi;
}


inline double Polarization::GetCsecPrefactor(double E1, double q0) const {
  double p3 = sqrt((E1-q0)*(E1-q0) - st.M3*st.M3); 
  const double fac = mG2*pow(Constants::GfMeV, 2)/
    pow(2.0*Constants::Pi, 3)/16.0;
  double a;
  if (current==1) {
    a=E1*(1.0 - exp((st.Mu4 - st.Mu2 - q0)/st.T));
  } else {
    a=E1*(1.0 - exp((0.0 - q0)/st.T));
  }
  if (doBlock) p3 *= 1.0 - 1.0/(exp((E1-q0-st.Mu3)/st.T) + 1.0);
  return fac*p3/a;
}

double Polarization::CalculateDifGam(double E1, double q0, double q) const {
  double fac = GetCsecPrefactor(E1, q0);
  double pi = GetResponse(E1, q0, q);
  return fac*pi;
}

