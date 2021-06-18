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

#include "Tensor.hpp" 
#include "PolarizationNonRel.hpp"
#include "Constants.hpp"
#include "FunctionIntegrator.hpp"

using namespace nuopac; 

// Use unnamed namespace so these methods are only locally available
namespace {
  
  // Exact expression
  inline double Fermi0(double eta) {
    if (eta>40.0) return eta;  
    return log(exp(eta) + 1.0);
  }

  
  

  double f (double x, void * params) {
       double q0=x;
      // double q0fix = ((double *)params)[0];
       double q= ((double *)params)[1];
       double U2=((double *)params)[2];
       double U4=((double *)params)[3];
       double M2=((double *)params)[4];
       double M4=((double *)params)[5];
       double Mu2=((double *)params)[6];
       double Mu4=((double *)params)[7];
       double T=((double *)params)[8];
       
/*       double q0t = q0 + U2 - U4;
       double qa2t = q0t*q0t - q*q;
       if (qa2t > pow(M2 - M4, 2)*0.0) return 0.0;
       if (qa2t < 1.e-1 && qa2t >= 0.0) return 0.0;
       double beta = 1.0 + (M2*M2 - M4*M4)/qa2t;
       double arg = beta*beta - 4.0*M2*M2/qa2t;
       if (arg<0.0) return 0.0;
       double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), M2);
       double delta2 = (Mu2 - U2 - em)/T;
       double delta4 = (Mu4 - U4 - em - q0t)/T;
       
       //following is the new em completely consistent with Redddy's thesis, for non-rel+interacting gas
  double chi=1-M4/M2;
  double c=q0+U2-U4-q*q/(2*M4);
  double emNC=std::max((-c*M2/q)*(-c*M2/q)/(2*M2),0.0); //the minimum E2 for NC reaction

  double argCC=1+2*chi*M4*c/(q*q);
  if (argCC<0.0) return 0.0;
  double eminCC=2*q*q/(chi*chi)*(1+chi*M4*c/(q*q)-sqrt(argCC))/(2*M2); //the minimum E2 for CC reaction
  double emaxCC=2*q*q/(chi*chi)*(1+chi*M4*c/(q*q)+sqrt(argCC))/(2*M2); //the maximum E2 for CC reaction

  double delta2NC=(Mu2 - U2 - emNC)/T;
  double delta4NC=(Mu4 - U4 - emNC-q0t)/T;

  double delta2minCC=(Mu2 - U2 - eminCC)/T;
  double delta4minCC=(Mu4 - U4 - eminCC-q0t)/T;

  double delta2maxCC=(Mu2 - U2 - emaxCC)/T;
  double delta4maxCC=(Mu4 - U4 - emaxCC-q0t)/T;
  
  double xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);  // reddy nc current
  double ximinCC=Fermi0(-delta2minCC) - Fermi0(-delta4minCC); //Reddy cc current
  double ximaxCC=Fermi0(-delta2maxCC) - Fermi0(-delta4maxCC); //Reddy cc current

       double Gamma0 = Fermi0(delta2) - Fermi0(delta4);
       double PI = Constants::Pi;
      // double piL = M2*M4*T/(PI*q)*Gamma0;//orig one

      //  double piL= M2*M4*T/(PI*q)*(xiNC+q0t/T);//neutral current consistent with Reddy's thesis
       double piL= M2*M4*T/(PI*q)*(ximinCC-ximaxCC);//charged current consistent with Reddy's thesis*/
      
       //for neutral current axial part mixed gas p+n
       double q0t = q0 + U2 - U2;
       double qa2t = q0t*q0t - q*q;
       if (qa2t > pow(M2 - M2, 2)*0.0) return 0.0;
       if (qa2t < 1.e-1 && qa2t >= 0.0) return 0.0;
       double beta = 1.0 + (M2*M2 - M2*M2)/qa2t;
       double arg = beta*beta - 4.0*M2*M2/qa2t;
       if (arg<0.0) return 0.0;
       double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), M2);
       double delta2 = (Mu2 - U2 - em)/T;
       double delta4 = (Mu2 - U2 - em - q0t)/T;

       //following is the new em completely consistent with Redddy's thesis, for non-rel+interacting gas
  
       double c=q0+U2-U2-q*q/(2*M2);
       double emNC=std::max((-c*M2/q)*(-c*M2/q)/(2*M2),0.0); //the minimum E2 for NC reaction
       double delta2NC=(Mu2 - U2 - emNC)/T;
       double delta4NC=(Mu2 - U2 - emNC-q0t)/T;
       double xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);  // reddy nc current
       double PI = Constants::Pi;
       double piLN= M2*M2*T/(PI*q)*(xiNC+q0t/T);//neutral current consistent with Reddy's thesis


       double Gamma0 = Fermi0(delta2) - Fermi0(delta4);
      
      // double piLN = M2*M2*T/(PI*q)*Gamma0;

        q0t = q0 + U4 - U4;
        qa2t = q0t*q0t - q*q;
       if (qa2t > pow(M4 - M4, 2)*0.0) return 0.0;
       if (qa2t < 1.e-1 && qa2t >= 0.0) return 0.0;
        beta = 1.0 + (M4*M4 - M4*M4)/qa2t;
        arg = beta*beta - 4.0*M4*M4/qa2t;
       if (arg<0.0) return 0.0;
        em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), M4);
        delta2 = (Mu4 - U4 - em)/T;
        delta4 = (Mu4 - U4 - em - q0t)/T;
	 //following is the new em completely consistent with Redddy's thesis, for non-rel+interacting gas

       c=q0+U4-U4-q*q/(2*M4);
       emNC=std::max((-c*M4/q)*(-c*M4/q)/(2*M4),0.0); //the minimum E2 for NC reaction
       delta2NC=(Mu4 - U4 - emNC)/T;
       delta4NC=(Mu4 - U4 - emNC-q0t)/T;
       xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);  // reddy nc current
      double piLP= M4*M4*T/(PI*q)*(xiNC+q0t/T);//neutral current consistent with Reddy's thesis

        Gamma0 = Fermi0(delta2) - Fermi0(delta4);
       
     //  double piLP = M4*M4*T/(PI*q)*Gamma0;
      
        double piL=piLN+piLP;
       // double piL=piLN;
      
       double f=piL;
      
       return f;
     }
}

std::array<double, 4> PolarizationNonRel::CalculateBasePolarizations(double q0, 
    double q) const {
	
//void Polarization::SetPolarizations(double q0, double q) {
  // Calculate some kinematic factors
 
  double q0t = q0 + st.U2 - st.U4;//orig one
 
  double qa2t = q0t*q0t - q*q; 
  // I don't completely understand this condition, but it seems to be necessary 
  // to suppress noise at larger q_0
  if (qa2t > pow(st.M2 - st.M4, 2)*0.0) return {0.0, 0.0, 0.0, 0.0};//orig one
  
  if (qa2t < 1.e-1 && qa2t >= 0.0) return {0.0, 0.0, 0.0, 0.0}; 
  double beta = 1.0 + (st.M2*st.M2 - st.M4*st.M4)/qa2t;//orig
  
  double arg = beta*beta - 4.0*st.M2*st.M2/qa2t;
  if (arg<0.0) return {0.0, 0.0, 0.0, 0.0}; 
  double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), st.M2);
  double delta2 = (st.Mu2 - st.U2 - em)/st.T;
  double delta4 = (st.Mu4 - st.U4 - em - q0t)/st.T;//orig one
  
  //following is the new em completely consistent with Redddy's thesis, for non-rel+interacting gas
  double chi=1-st.M4/st.M2;
  double c=q0+st.U2-st.U4-q*q/(2*st.M4);
  double emNC=std::max((-c*st.M2/q)*(-c*st.M2/q)/(2*st.M2),0.0); //the minimum E2 for NC reaction

  double argCC=1+2*chi*st.M4*c/(q*q);
  if (argCC<0.0) return {0.0, 0.0, 0.0, 0.0};
  double eminCC=2*q*q/(chi*chi)*(1+chi*st.M4*c/(q*q)-sqrt(argCC))/(2*st.M2); //the minimum E2 for CC reaction
  double emaxCC=2*q*q/(chi*chi)*(1+chi*st.M4*c/(q*q)+sqrt(argCC))/(2*st.M2); //the maximum E2 for CC reaction

  double delta2NC=(st.Mu2 - st.U2 - emNC)/st.T;
  double delta4NC=(st.Mu4 - st.U4 - emNC-q0t)/st.T;

  double delta2minCC=(st.Mu2 - st.U2 - eminCC)/st.T;
  double delta4minCC=(st.Mu4 - st.U4 - eminCC-q0t)/st.T;

  double delta2maxCC=(st.Mu2 - st.U2 - emaxCC)/st.T;
  double delta4maxCC=(st.Mu4 - st.U4 - emaxCC-q0t)/st.T;
 
 
  // Now just need to include some method for calculating these
  // At least at low density, Gamma0 should be the dominant term
  // which looks like the non-relativistic response function 
  // Under non-degenerate conditions (i.e. delta2, delta4 << 0), 
  // Gamma0 = Gamma1 = 0.5*Gamma2 
  // This is exact 
  double Gamma0 = Fermi0(delta2) - Fermi0(delta4);
  
  double xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);  // reddy nc current
  double ximinCC=Fermi0(-delta2minCC) - Fermi0(-delta4minCC); //Reddy cc current
  double ximaxCC=Fermi0(-delta2maxCC) - Fermi0(-delta4maxCC); //Reddy cc current

  double PI = Constants::Pi; 
 // double piL = st.M2*st.M4*st.T/(PI*q)*Gamma0;//orig one

 // double piL= st.M2*st.M4*st.T/(PI*q)*(xiNC+q0t/st.T);//neutral current consistent with Reddy's thesis
  double piL= st.M2*st.M4*st.T/(PI*q)*(ximinCC-ximaxCC);//charged current consistent with Reddy's thesis

  double piQ = 0.0; 
  double piM = 0.0;
  double piT = 0.0;
  
  return {piQ, piL, piM, piT};
}


//only used when calculating NC mixed gas (n+p) axial part,neutron basic polarization use m2 and u2; proton basic polarization use m4 and u4
std::array<double, 4> PolarizationNonRel::CalculateBasePolarizationsNeutron(double q0,
    double q) const {
//void Polarization::SetPolarizations(double q0, double q) {
  // Calculate some kinematic factors
   
  double q0t = q0 + st.U2 - st.U2;
  double qa2t = q0t*q0t - q*q;
  // I don't completely understand this condition, but it seems to be necessary
  // to suppress noise at larger q_0
  if (qa2t > pow(st.M2 - st.M2, 2)*0.0) return {0.0, 0.0, 0.0, 0.0};
  if (qa2t < 1.e-1 && qa2t >= 0.0) return {0.0, 0.0, 0.0, 0.0};
  double beta = 1.0 + (st.M2*st.M2 - st.M2*st.M4)/qa2t;
  double arg = beta*beta - 4.0*st.M2*st.M2/qa2t;
  if (arg<0.0) return {0.0, 0.0, 0.0, 0.0};
  double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), st.M2);
  double delta2 = (st.Mu2 - st.U2 - em)/st.T;
  double delta4 = (st.Mu2 - st.U2 - em - q0t)/st.T;

  //following is the new em completely consistent with Redddy's thesis, for non-rel+interacting gas
  double chi=1-st.M2/st.M2;
  double c=q0+st.U2-st.U2-q*q/(2*st.M2);
  double emNC=std::max((-c*st.M2/q)*(-c*st.M2/q)/(2*st.M2),0.0); //the minimum E2 for NC reaction
 

  double delta2NC=(st.Mu2 - st.U2 - emNC)/st.T;
  double delta4NC=(st.Mu2 - st.U2 - emNC-q0t)/st.T;


  // Now just need to include some method for calculating these
  // At least at low density, Gamma0 should be the dominant term
  // which looks like the non-relativistic response function
  // Under non-degenerate conditions (i.e. delta2, delta4 << 0),
  // Gamma0 = Gamma1 = 0.5*Gamma2
  // This is exact
  double Gamma0 = Fermi0(delta2) - Fermi0(delta4);

  double xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);  // reddy nc current


  double PI = Constants::Pi;
 // double piL = st.M2*st.M2*st.T/(PI*q)*Gamma0;//orig one

  double piL= st.M2*st.M2*st.T/(PI*q)*(xiNC+q0t/st.T);//neutral current consistent with Reddy's thesis
  

  double piQ = 0.0;
  double piM = 0.0;
  double piT = 0.0;


  return {piQ, piL, piM, piT};
}

std::array<double, 4> PolarizationNonRel::CalculateBasePolarizationsProton(double q0,
    double q) const {
//void Polarization::SetPolarizations(double q0, double q) {
  // Calculate some kinematic factors
  
  double q0t = q0 + st.U4 - st.U4;
  double qa2t = q0t*q0t - q*q;
  // I don't completely understand this condition, but it seems to be necessary
  // to suppress noise at larger q_0
  if (qa2t > pow(st.M4 - st.M4, 2)*0.0) return {0.0, 0.0, 0.0, 0.0};
  if (qa2t < 1.e-1 && qa2t >= 0.0) return {0.0, 0.0, 0.0, 0.0};
  double beta = 1.0 + (st.M4*st.M4 - st.M4*st.M4)/qa2t;
  double arg = beta*beta - 4.0*st.M4*st.M4/qa2t;
  if (arg<0.0) return {0.0, 0.0, 0.0, 0.0};
  double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), st.M4);
  double delta2 = (st.Mu4 - st.U4 - em)/st.T;
  double delta4 = (st.Mu4 - st.U4 - em - q0t)/st.T;

   //following is the new em completely consistent with Redddy's thesis, for non-rel+interacting gas
  
  double c=q0+st.U4-st.U4-q*q/(2*st.M4);
  double emNC=std::max((-c*st.M4/q)*(-c*st.M4/q)/(2*st.M4),0.0); //the minimum E2 for NC reaction


  double delta2NC=(st.Mu4 - st.U4 - emNC)/st.T;
  double delta4NC=(st.Mu4 - st.U4 - emNC-q0t)/st.T;


  // Now just need to include some method for calculating these
  // At least at low density, Gamma0 should be the dominant term
  // which looks like the non-relativistic response function
  // Under non-degenerate conditions (i.e. delta2, delta4 << 0),
  // Gamma0 = Gamma1 = 0.5*Gamma2
  // This is exact
  double Gamma0 = Fermi0(delta2) - Fermi0(delta4);

  double xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);  // reddy nc current


  double PI = Constants::Pi;
 // double piL = st.M2*st.M2*st.T/(PI*q)*Gamma0;//orig one

  double piL= st.M4*st.M4*st.T/(PI*q)*(xiNC+q0t/st.T);//neutral current consistent with Reddy's thesis


  double piQ = 0.0;
  double piM = 0.0;
  double piT = 0.0;


  return {piQ, piL, piM, piT};
}


  double PolarizationNonRel::gamma0(double q0, double q) const {
//void Polarization::SetPolarizations(double q0, double q) {
  // Calculate some kinematic factors
  double mu2pol;
  double mu4pol;
  double u2pol;
  double u4pol;
  double m2pol;
  double m4pol;
  mu2pol=st.Mu2;
  mu4pol=st.Mu2;
  u2pol=st.U2;
  u4pol=st.U2;
  m2pol=st.M2;
  m4pol=st.M2;

  double q0t = q0 + u2pol - u4pol;
  double qa2t = q0t*q0t - q*q;
  // I don't completely understand this condition, but it seems to be necessary
  // to suppress noise at larger q_0
 
  
  double beta = 1.0 + (st.M2*st.M2 - st.M4*st.M4)/qa2t;
  double arg = beta*beta - 4.0*st.M2*st.M2/qa2t;
  
  double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), st.M2);
  double delta2 = (mu2pol - u2pol - em)/st.T;
//  double delta4 = (st.Mu4 - st.U4 - em - q0t)/st.T;//orig one
   double delta4 = (mu4pol - u4pol - em - q0t)/st.T;
  // Now just need to include some method for calculating these
  // At least at low density, Gamma0 should be the dominant term
  // which looks like the non-relativistic response function
  // Under non-degenerate conditions (i.e. delta2, delta4 << 0),
  // Gamma0 = Gamma1 = 0.5*Gamma2
  // This is exact
 // double Gamma0 = Fermi0(delta2) - Fermi0(delta4);
    double Gamma0 = Fermi0(delta4);
 // double deltaU=st.U2-st.U4;
 // double deltamu=st.Mu2-st.Mu4;
// double Delta=delta2-delta4;
  return Gamma0;
}


  double PolarizationNonRel::GetImPI( double q0, double q) const {


 // auto pt = CalculateBasePolarizations(q0, q);

 
 
  // double piL = pt[1];
  
   auto ptN=CalculateBasePolarizationsNeutron(q0, q);
   auto ptP=CalculateBasePolarizationsProton(q0, q);
  double piL=ptN[1]+ptP[1];//only for NC mixed gas (N+P) axial channel, in vector channel there is only contribution from neutron, just like pure neutron gas
 // double piL=ptN[1];
  return piL;
}

double PolarizationNonRel::GetImPI2( double q0, double q) const {
  
  double T=10;
  double M2=939; 
  double vel=sqrt((3.0*T)/M2);
  double wmin;
  double wmax;
 // wmin=-3.0*vel*q-0.00000789;
//  wmin=-50.0+1.0E-6;
    wmin=-200.0;
 // wmax=3.0*(vel*q+q*q/(2.0*M2));
 // wmax=50.0;
   wmax=200.0;
  double dw=(wmax-wmin)/1000.0;
  double piLRe=0;
//  FunctionIntegrator integrate(512, 0.0, 1.0E-12);
   
//  auto integrand = [q0, q] (double xx) {
     
//    auto integ=GetImPI(xx,q)/(xx-q0);
 //   return  integ;
 // };
 // piLRe = integrate.Integrate(integrand, wmin, wmax);
 // test Integration**************************************
  // auto integrand = [q0] (double xx) {
    
    //   auto f = log(q0*xx) / sqrt(xx);
    //   return f;
    
  // };

   gsl_integration_workspace * w 
         = gsl_integration_workspace_alloc (1024);
   gsl_function F;
       
       double par[]={q0,q,st.U2,st.U4,st.M2,st.M4,st.Mu2,st.Mu4,st.T};
      // double par[]={q0,q};
       F.function = &f;
       F.params = &par;
      
   double result, error;

  // gsl_integration_qawc (&F, -90, 90, q0, 0, 1e-7, 1000,
   //                          w, &result, &error);
    gsl_integration_qawc (&F, -300, 100.0, q0, 0, 1e-7, 1000,
                             w, &result, &error);
  // gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
    //                         w, &result, &error);


 // ******************************************************
// for (int i=0; i<1000; i++){ 
//   piLRe=piLRe+GetImPI(wmin+dw*i, q)/(wmin+dw*i-q0)*dw; 
// }
 // double piL=GetImPI(q0, q);
 // double vf=(-0.74+2.5)*1.0E-6;
 // double vgt=4.5*1.0E-6;
 // double vrpa=vgt;
 // piL=piL/((1-vrpa*piLRe)*(1-vrpa*piLRe)+vrpa*vrpa*piL*piL);
 // piL=piL/((1-vf*piLRe)*(1-vf*piLRe)+vf*vf*piL*piL);
  piLRe=-result/Constants::Pi;

  gsl_integration_workspace_free (w);

  return piLRe;
}

void PolarizationNonRel::PrintGetImPI2(double q0, double q) const {
 
 /* double RepiL=0;
  double T=10;
  double M2=939;
  double vel=sqrt((3.0*T)/M2);
  double wmin;
  double wmax;
 
  wmin=-3.0*vel*q-0.00000789;
  wmax=3.0*(vel*q+q*q/(2.0*M2));
  double dw=(wmax-wmin)/100.0;
  double piLRe=0;
  for (int i=0; i++; i<100) {
   piLRe=piLRe+GetImPI(wmin+dw, q)/(wmin+dw-q0)*dw;
   std::cout << " piLRe: "<<piLRe<<std::endl;
  }*/

  
}


void PolarizationNonRel::SetPolarizations(double q0, double q,
      Tensor<double>* piVV, 
      Tensor<double>* piAA, 
      Tensor<double>* piTT, 
      Tensor<double>* piVA, 
      Tensor<double>* piVT, 
      Tensor<double>* piAT) const {
 
  // Calculate the basic parts of the polarization 
  auto pt = CalculateBasePolarizations(q0, q);
  double piLRe=GetImPI2(q0, q); 
  double piQ = pt[0]; 
 // double piL = pt[1]; 
  double piM = pt[2]; 
  double piT = pt[3];
  double piL;///only turn this on when the getImPI is used to calculate the NC axial part for mixed gas n+p
  piL=GetImPI(q0, q);//only turn this on when the getImPI is used to calculate the NC axial part for mixed gas n+p 
 // double vf=(-0.74+2.5)*1.0E-5;
 // double vgt=4.5*1.0E-5;
  //vf and vgt for cc comes from arXiv:1205.4066, in unit of MeV^-2, note that in the arxiv paper the unit is wrong, should be fm^2 rather than fm^-2
 // double vf=5.1*1.0E-5;
 
 // double vgt=2.8*1.0E-5;
 
 // double vrpa=vf;
 
 //virial 1st order approximation
  // double vf=-8.45*1.0E-5;
  // double vgt=7.6*1.0E-5;
  // double vrpa=vf;
 
 // ******************insert density dependent Landau Fermi Liquid Parameter********************************************
    double densFac = pow(Constants::HBCFmMeV, 3);

    //**********parameter fitting residual interaction curves*************
    //***nrapr curves*************
    /* double bfsym1=1.9635*1.0E-6;
     double bfsym2=0.00003499;
     double bfsym3=0.00007305;
     double bgtsym1=-2.6304*1.0E-6;
     double bgtsym2=-0.00003932;
     double bgtsym3=-0.00005173;
     double bfpureN1=1.1266645*1.0E-6;
     double bfpureN2=0.00002380;
     double bfpureN3=0.0000262075;
     double bgtpureN1=-3.69989*1.0E-6;
     double bgtpureN2=-0.0000517686;
     double bgtpureN3=-0.00008783369;*/
     //************************************
     //***ft1 curves***********************
    /* double bfsym1=2.461399*1.0E-6;
     double bfsym2=0.00004130;
     double bfsym3=0.00009276;
     double bgtsym1=-8.10023*1.0E-6;
     double bgtsym2=-0.00010199;
     double bgtsym3=-0.000232676;
     double bfpureN1=1.36189*1.0E-6;
     double bfpureN2=0.00002643;
     double bfpureN3=0.00003437;
     double bgtpureN1=-1.20574*1.0E-5;
     double bgtpureN2=-0.000148040;
     double bgtpureN3=-0.000366266;*/
      //***ft2 curves***********************
    /* double bfsym1=8.583149717*1.0E-7;
     double bfsym2=0.000014974664370833527;
     double bfsym3=-0.00004505743409405033;
     double bgtsym1=-9.041128573465587*1.0E-6;
     double bgtsym2=-0.0001052444786;
     double bgtsym3=-0.00017917481;
     double bfpureN1=-7.52635101187*1.0E-7;
     double bfpureN2=-8.708075046558*1.0E-6;
     double bfpureN3=-0.000151316469989;
     double bgtpureN1=-0.0000143483855;
     double bgtpureN2=-0.000167;
     double bgtpureN3=-0.00035834;*/
     //*******************************************
     //***don't fit curves, calulate residual interaction directly************
     double t0,t1,t2,t3,x0,x1,x2,x3,epsilon;
     //****set skyrme parameters********************
     //NRAPR
      /* epsilon= 1.4416*1.0E-01;
       t0 = -2.7197*1.0E+03;
       t1=4.1764*1.0E+02;
       t2= -6.6687*1.0E+01;
       t3= 1.5042*1.0E+04;
       x0= 1.6154*1.0E-01;
       x1= -4.7986*1.0E-02;
       x2= 2.717*1.0E-02;
       x3= 1.3611*1.0E-01;*/

     //ft1
      /* epsilon= 2.2834001699*1.0E-01;
       t0 = -2.2233370185*1.0E+03;
       t1=3.0392867870*1.0E+02;
       t2= 1.2213752541*1.0E+03;
       t3= 1.3992595172*1.0E+04;
       x0= 1.7137205580*1.0E-01;
       x1= -3.6742222012*1.0E+00;
       x2= -1.3703338290*1.0E+00;
       x3= 1.3247942498*1.0E-01;*/


     //ft2
      /* epsilon= 2.1493862350*1.0E-01;
       t0 = -2.3279811323*1.0E+03;
       t1=2.9748565632*1.0E+02;
       t2= 1.8783378926*1.0E+03;
       t3= 1.4661502679*1.0E+04;
       x0= 1.6575377249*1.0E-01;
       x1= -5.3459402501*1.0E+00;
       x2= -1.3322633637*1.0E+00;
       x3= 9.9611039279*1.0E-02;*/
     //ft3
      /* epsilon= 2.5851613108*1.0E-01;
       t0 = -2.0515194890*1.0E+03;
       t1=3.0250225425*1.0E+02;
       t2=1.8419851615*1.0E+03;
       t3=1.3102289407 *1.0E+04;
       x0=2.0105379497 *1.0E-01;
       x1=-5.2329865972 *1.0E+00;
       x2=-1.3266769390 *1.0E+00;
       x3= 1.5916860574*1.0E-01;*/
      //ft4
      /* epsilon= 1.5530331561*1.0E-01;
       t0 = -2.9932037792*1.0E+03;
       t1=2.8339369142*1.0E+02;
       t2=2.3728464374*1.0E+03;
       t3=1.8587502887 *1.0E+04;
       x0=2.0121441704 *1.0E-01;
       x1=-6.7808781839 *1.0E+00;
       x2=-1.3212716203 *1.0E+00;
       x3=1.4007920173 *1.0E-01;*/
      //ft5
      /* epsilon= 1.6949727833*1.0E-01;
       t0 = -2.7182493297*1.0E+03;
       t1=3.1637903698*1.0E+02;
       t2=1.6805765009*1.0E+03;
       t3=1.6565776882 *1.0E+04;
       x0=2.3023295140 *1.0E-01;
       x1=-4.6998516037 *1.0E+00;
       x2=-1.3385260782 *1.0E+00;
       x3=2.4599697984*1.0E-01;*/
       //ft6
       epsilon= 1.4339134197*1.0E-01;
       t0 = -3.1535549664*1.0E+03;
       t1=2.9811820280*1.0E+02;
       t2=2.3615846306*1.0E+03;
       t3=1.9373593217*1.0E+04;
       x0=1.5139360999*1.0E-01;
       x1=-6.5192138295 *1.0E+00;
       x2=-1.3213613298 *1.0E+00;
       x3=5.8876078909*1.0E-01;
       //ft7
      /* epsilon= 1.6112652355*1.0E-01;
       t0=-2.8515584390*1.0E+03;
       t1=2.9926667400*1.0E+02;
       t2=2.1204461670*1.0E+03;
       t3=1.7480390049*1.0E+04;
       x0=9.5708919697*1.0E-02;
       x1=-5.9493754771*1.0E+00;
       x2=-1.3211268948 *1.0E+00;
       x3=3.7117434733*1.0E-02;*/


     //*************************
     double rou=(st.n2+st.n4)/densFac;
     double kf=pow(0.5*rou*3*3.14159*3.14159,1.0/3.0);
     double kfn=pow(rou*3*3.14159*3.14159,1.0/3.0);
     double f0=3.0/4.0*t0+3.0/8.0*t1*kf*kf+5.0/8.0*t2*kf*kf+1.0/2.0*t2*x2*kf*kf+(epsilon+1.0)*(epsilon+2.0)/16.0*t3*pow(rou,epsilon);
     double f0p=-1.0/4.0*t0-1.0/2.0*t0*x0-1.0/8.0*t1*kf*kf-1.0/4.0*t1*x1*kf*kf+1.0/8.0*t2*kf*kf+1.0/4.0*t2*x2*kf*kf-1.0/24.0*t3*pow(rou,epsilon)-1.0/12.0*t3*x3*pow(rou,epsilon);
     double g0=-1.0/4.0*t0+1.0/2.0*x0*t0-1.0/8.0*t1*kf*kf+1.0/4.0*t1*x1*kf*kf+1.0/8.0*t2*kf*kf+1.0/4.0*t2*x2*kf*kf-1.0/24.0*t3*pow(rou,epsilon)+1.0/12.0*t3*x3*pow(rou,epsilon);
     double g0p=-1.0/4.0*t0-1.0/8.0*t1*kf*kf+1.0/8.0*t2*kf*kf-1.0/24.0*t3*pow(rou,epsilon);
     double f0n=1.0/2.0*t0-1.0/2.0*t0*x0+1.0/8.0*t1*kfn*kfn-1.0/4.0*t1*x1*kfn*kfn+3.0/4.0*t2*kfn*kfn+3.0/4.0*t2*x2*kfn*kfn+(epsilon+1.0)*(epsilon+2.0)/24.0*t3*pow(rou,epsilon)-(epsilon+1.0)*(epsilon+2.0)/24.0*t3*x3*pow(rou,epsilon);
     double g0n=-1.0/2.0*t0+1.0/2.0*t0*x0-1.0/4.0*t1*kfn*kfn+1.0/4.0*t1*x1*kfn*kfn+1.0/4.0*t2*kfn*kfn+1.0/4.0*t2*x2*kfn*kfn-1.0/12.0*t3*pow(rou,epsilon)+1.0/12.0*t3*x3*pow(rou,epsilon);
      //***for non-pure matter
      double vfs=f0+f0p;
      double vgts=g0+g0p;
      double vfn=f0n;
      double vgtn=g0n;
      double xx=(st.n2-st.n4)/(st.n2+st.n4);
     // double xx=0.0; //test
      double vf=vfs*(1-xx)+vfn*xx;
      double vgt=vgts*(1-xx)+vgtn*xx;
      double vrpa=vgt/pow(197.3,3);
     
    //*************for non-pure neutron matter************
    //*******symmetric nuclear metter landau fermi liquid***
//    double vfs=bfsym1*log((st.n2+st.n4)/densFac)*log((st.n2+st.n4)/densFac)+bfsym2*log((st.n2+st.n4)/densFac)+bfsym3;
//    double vgts=bgtsym1*log((st.n2+st.n4)/densFac)*log((st.n2+st.n4)/densFac)+bgtsym2*log((st.n2+st.n4)/densFac)+bgtsym3;
    //******pure neutron matter landau fermi liquid***
//    double vfn=bfpureN1*log((st.n2+st.n4)/densFac)*log((st.n2+st.n4)/densFac)+bfpureN2*log((st.n2+st.n4)/densFac)+bfpureN3;
//    double vgtn=bgtpureN1*log((st.n2+st.n4)/densFac)*log((st.n2+st.n4)/densFac)+bgtpureN2*log((st.n2+st.n4)/densFac)+bgtpureN3;
    //linear interpolate any state between symmetric and pure neutron****************
   // double xx=(st.n2-st.n4)/(st.n2+st.n4);
   // double xx=0.0; //test
   // double vf=vfs*(1-xx)+vfn*xx;
   // double vgt=vgts*(1-xx)+vgtn*xx;

    //*************for pure neutron matter**********************
  //  double vf=bfpureN1*log((st.n2)/densFac)*log((st.n2)/densFac)+bfpureN2*log((st.n2)/densFac)+bfpureN3;
  //  double vgt=bgtpureN1*log((st.n2)/densFac)*log((st.n2)/densFac)+bgtpureN2*log((st.n2)/densFac)+bgtpureN3;//for NC pure neutron, st.n2=st.n4=n

    //********at low density, for virial EOS, since free energy only keep density terms up to 2nd order, the virial vf and vgt is not density dependent*********
  //  double vf=-0.000200209;
  //  double vgt=0.000191682;

//    double vrpa=vgt;
 // *************************************************************
  piL=2*(piL/2)/((1-vrpa*(piLRe/2))*(1-vrpa*(piLRe/2))+vrpa*vrpa*(piL/2)*(piL/2));//piL is 2*Im PI  
  // Set the different parts of the polarization 
  piVV->L   = piL; 
  piAA->Tp  = 0.5*piL; 

}

void PolarizationNonRel::SetLeptonTensor(double E1, double q0, double q, Tensor<double>* L) const {
  double E3 = E1 - q0;
  double mu13 = (q*q - E1*E1 - E3*E3)/(2.0*E1*E3);
  L->L  = 8.0*E1*E3*(1 + mu13); 
  L->Tp = 8.0*E1*E3*(3 - mu13); 
}
  
  
