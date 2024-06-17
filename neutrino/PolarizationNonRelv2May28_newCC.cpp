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
/**********************************************************
 * on June 28 2023, zidu totally removed the kinematics constraint that q0^2-q^2>0 is forbidden. It is satisfied
 * when requiring lepton kinematics setting E1=E3+q0 and q=sqrt(E1^2+E3^2-2E1E3costheta). Note this is Kinematic
 * constraint, and this constraint should be respected when picking q and q0 of the ImPi(q0,q). However, this constraint
 * SHOULD NOT be used when calculating the real part of the polarization function. That is to say, all the kinematic constraint
 * in the polarization function in PolarizationNorel should be disabled.
 *
 * ***********************************************************/

/// \file Polarization.cpp
/// \author lroberts
/// \since Apr 02, 2016
///
/// \brief
///
///
#include <math.h> 
#include <limits> 

#include "tensor.h" 
#include "PolarizationNonRel_newCC.hpp"
#include "constants.h"
#include "FunctionIntegrator.hpp"

#include <o2scl/constants.h>

using namespace std; 
using namespace nuopac; 

// Use unnamed namespace so these methods are only locally available
namespace {

  double FermiE (double x, void * params) {
          double epsilon=x;
	  double mu=((double *)params)[0];
	  double T=((double *)params)[1];
	  double integrand=x/(1+exp(x/T-mu/T));
	  return integrand;
  }

  double FermiEwtilde (double x, void * params) {
          double epsilon=x;
          double mu=((double *)params)[0];
          double T=((double *)params)[1];
	  double wtilde=((double *)params)[2];
          double integrand=(x-wtilde)/(1+exp(x/T-mu/T));
          return integrand;
  }

  double FermiE2 (double x, void * params) {
          double epsilon=x;
          double mu=((double *)params)[0];
          double T=((double *)params)[1];
         // double wtilde=((double *)params)[2];
          double integrand=x*x/(1+exp(x/T-mu/T));
          return integrand;
  }

  double FermiE2wtilde (double x, void * params) {
          double epsilon=x;
          double mu=((double *)params)[0];
          double T=((double *)params)[1];
          double wtilde=((double *)params)[2];
          double integrand=(x-wtilde)*(x-wtilde)/(1+exp(x/T-mu/T));
          return integrand;
  }

  double FermiE2kminus (double x, void * params) {
          double epsilon=x;
          double mu=((double *)params)[0];
          double T=((double *)params)[1];
          double wtilde=((double *)params)[2];
	  double kminus=((double *)params)[3];
	  double q=((double *)params)[4];
	  double deltams=((double *)params)[5];
          double integrand=(kminus*q+deltams*x)*(kminus*q+deltams*x)/(1+exp(x/T-mu/T));
          return integrand;
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
       
    double q0t = q0 + U2 - U4;
    double qa2t = q0t*q0t - q*q;
   // if (qa2t > pow(M2 - M4, 2)*0.0) return 0.0;
   // if (qa2t < 1.e-1 && qa2t >= 0.0) return 0.0;
    double beta = 1.0 + (M2*M2 - M4*M4)/qa2t;
    double arg = beta*beta - 4.0*M2*M2/qa2t;
    if (arg<0.0) return 0.0;
    
   // double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), M2);
   // double delta2 = (Mu2 - U2 - em)/T;
   // double delta4 = (Mu4 - U4 - em - q0t)/T;
       
    //following is the new em completely consistent with Redddy's
    //thesis, for non-rel+interacting gas
    double chi=1-M4/M2;
    double c=q0+U2-U4-q*q/(2*M4);
    
    //the minimum E2 for NC reaction
   // double emNC=std::max((-c*M2/q)*(-c*M2/q)/(2*M2),0.0); 

    double argCC=1+2*chi*M4*c/(q*q);
    if (argCC<0.0) return 0.0;
    //the minimum E2 for CC reaction
    double eminCC=2*q*q/(chi*chi)*(1+chi*M4*c/(q*q)-sqrt(argCC))/(2*M2);
    //the maximum E2 for CC reaction
    double emaxCC=2*q*q/(chi*chi)*(1+chi*M4*c/(q*q)+sqrt(argCC))/(2*M2); 

   // double delta2NC=(Mu2 - U2 - emNC)/T;
   // double delta4NC=(Mu4 - U4 - emNC-q0t)/T;

    double delta2minCC=(Mu2 - U2 - eminCC)/T;
    double delta4minCC=(Mu4 - U4 - eminCC-q0t)/T;

    double delta2maxCC=(Mu2 - U2 - emaxCC)/T;
    double delta4maxCC=(Mu4 - U4 - emaxCC-q0t)/T;

    // reddy nc current
   // double xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);
    //Reddy cc current
    double ximinCC=Fermi0(-delta2minCC) - Fermi0(-delta4minCC);
    //Reddy cc current
    double ximaxCC=Fermi0(-delta2maxCC) - Fermi0(-delta4maxCC); 

   // double Gamma0 = Fermi0(delta2) - Fermi0(delta4);
    double PI = o2scl_const::pi;//Constants::Pi;
    // double piL = M2*M4*T/(PI*q)*Gamma0;//orig one

    //if (this->current==this->current_neutral) {
    //O2SCL_ERR("Invalid current 4.",o2scl::exc_efailed);
    //}
    
    //neutral current consistent with Reddy's thesis
    //  double piL= M2*M4*T/(PI*q)*(xiNC+q0t/T);

    //charged current consistent with Reddy's thesis
    double piL= M2*M4*T/(PI*q)*(ximinCC-ximaxCC);
      
    double f=piL;
      
    return f;
  }

  
  double fn(double x, void *params) {
    
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

    double q0t = q0 + U2 - U2;
    double qa2t = q0t*q0t - q*q;
   // if (qa2t > pow(M2 - M2, 2)*0.0) return 0.0;
   // if (qa2t < 1.e-1 && qa2t >= 0.0) return 0.0;
   // double beta = 1.0 + (M2*M2 - M2*M2)/qa2t;
   // double arg = beta*beta - 4.0*M2*M2/qa2t;
   // if (arg<0.0) return 0.0;
   // double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), M2);
   // double delta2 = (Mu2 - U2 - em)/T;
   // double delta4 = (Mu2 - U2 - em - q0t)/T;

    //following is the new em completely consistent with Redddy's
    //thesis, for non-rel+interacting gas

    //if (this->current==this->current_charged) {
    //O2SCL_ERR("Invalid current 4.",o2scl::exc_efailed);
    //}

    double c=q0+U2-U2-q*q/(2*M2);
    //the minimum E2 for NC reaction
    double emNC=std::max((-c*M2/q)*(-c*M2/q)/(2*M2),0.0); 
    double delta2NC=(Mu2 - U2 - emNC)/T;
    double delta4NC=(Mu2 - U2 - emNC-q0t)/T;
    // reddy nc current
    double xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);  
    double PI = o2scl_const::pi;//Constants::Pi;
    //neutral current Neutron PiL consistent with Reddy's thesis
    double piLN= M2*M2*T/(PI*q)*(xiNC+q0t/T);
    return piLN;
  }

  double fp (double x, void * params) {
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

    double q0t = q0 + U4 - U4;
    double qa2t = q0t*q0t - q*q;
   // if (qa2t > pow(M4 - M4, 2)*0.0) return 0.0;
   // if (qa2t < 1.e-1 && qa2t >= 0.0) return 0.0;
   // double  beta = 1.0 + (M4*M4 - M4*M4)/qa2t;
   // double  arg = beta*beta - 4.0*M4*M4/qa2t;
   // if (arg<0.0) return 0.0;
   // double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), M4);
   // double delta2 = (Mu4 - U4 - em)/T;
   // double delta4 = (Mu4 - U4 - em - q0t)/T;
    
    //following is the new em completely consistent with Redddy's
    //thesis, for non-rel+interacting gas

    double c=q0+U4-U4-q*q/(2*M4);
    //the minimum E2 for NC reaction
    double emNC=std::max((-c*M4/q)*(-c*M4/q)/(2*M4),0.0); 
    double delta2NC=(Mu4 - U4 - emNC)/T;
    double delta4NC=(Mu4 - U4 - emNC-q0t)/T;
    // reddy nc current
    double xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);  
    double PI = o2scl_const::pi;//Constants::Pi;
    //neutral current Proton PiL consistent with Reddy's thesis
    double piLP= M4*M4*T/(PI*q)*(xiNC+q0t/T);

    //if (this->current==this->current_charged) {
    //O2SCL_ERR("Invalid current 5.",o2scl::exc_efailed);
    //}
    
    return piLP;
  }

}

std::array<double, 4> PolarizationNonRel::CalculateBasePolarizations
(double q0, double q) const {
  
  // Calculate some kinematic factors
 
  double q0t = q0 + st.U2 - st.U4;//orig one
 
  double qa2t = q0t*q0t - q*q;
  
  // [LR]: I don't completely understand this condition, but it seems to be
  // necessary to suppress noise at larger q_0
  
  //orig one
 // if (qa2t > pow(st.M2 - st.M4, 2)*0.0) return {0.0, 0.0, 0.0, 0.0};
  
 // if (qa2t < 1.e-1 && qa2t >= 0.0) return {0.0, 0.0, 0.0, 0.0}; 
  double beta = 1.0 + (st.M2*st.M2 - st.M4*st.M4)/qa2t;//orig
  
  double arg = beta*beta - 4.0*st.M2*st.M2/qa2t;
  if (arg<0.0) return {0.0, 0.0, 0.0, 0.0}; 
 // double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), st.M2);
 // double delta2 = (st.Mu2 - st.U2 - em)/st.T;
 // double delta4 = (st.Mu4 - st.U4 - em - q0t)/st.T;//orig one
  
  // following is the new em completely consistent with Redddy's
  // thesis, for non-rel+interacting gas
  double chi=1-st.M4/st.M2; //orig
  
  double c=q0+st.U2-st.U4-q*q/(2*st.M4);
  //the minimum E2 for NC reaction
 // double emNC=std::max((-c*st.M2/q)*(-c*st.M2/q)/(2*st.M2),0.0); 

  double argCC=1+2*chi*st.M4*c/(q*q);
  if (argCC<0.0) return {0.0, 0.0, 0.0, 0.0};
  //the minimum E2 for CC reaction
  double eminCC=2*q*q/(chi*chi)*(1+chi*st.M4*c/(q*q)-sqrt(argCC))/(2*st.M2);
  //the maximum E2 for CC reaction
  double emaxCC=2*q*q/(chi*chi)*(1+chi*st.M4*c/(q*q)+sqrt(argCC))/(2*st.M2); 

 // double delta2NC=(st.Mu2 - st.U2 - emNC)/st.T;
 // double delta4NC=(st.Mu4 - st.U4 - emNC-q0t)/st.T;

  double delta2minCC=(st.Mu2 - st.U2 - eminCC)/st.T;
  double delta4minCC=(st.Mu4 - st.U4 - eminCC-q0t)/st.T;

  double delta2maxCC=(st.Mu2 - st.U2 - emaxCC)/st.T;
  double delta4maxCC=(st.Mu4 - st.U4 - emaxCC-q0t)/st.T;
 
  // [LR]: Now just need to include some method for calculating these
  // At least at low density, Gamma0 should be the dominant term
  // which looks like the non-relativistic response function 
  // Under non-degenerate conditions (i.e. delta2, delta4 << 0), 
  // Gamma0 = Gamma1 = 0.5*Gamma2 
  // This is exact 
 // double Gamma0 = Fermi0(delta2) - Fermi0(delta4);

  // reddy nc current
 // double xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);
  //Reddy cc current
  double ximinCC=Fermi0(-delta2minCC) - Fermi0(-delta4minCC);
  //Reddy cc current
  double ximaxCC=Fermi0(-delta2maxCC) - Fermi0(-delta4maxCC); 

  //if (current==current_neutral) {
  //O2SCL_ERR("Invalid current 1.",o2scl::exc_efailed);
  //}
  
  double PI = o2scl_const::pi;//Constants::Pi; 
  // double piL = st.M2*st.M4*st.T/(PI*q)*Gamma0;//orig one

  //neutral current consistent with Reddy's thesis
  // double piL= st.M2*st.M4*st.T/(PI*q)*(xiNC+q0t/st.T);

  //charged current consistent with Reddy's thesis
  double piL= st.M2*st.M4*st.T/(PI*q)*(ximinCC-ximaxCC);

  tempy=piL;
  
  if (integral_debug) {
    cout << "XX: " << q << " " << q0 << " " << st.M2 << " " << st.M4  << " "
         << delta2minCC << " " << ximinCC-ximaxCC << " " << tempy << endl;
  }
  
  double piQ = 0.0; 
  double piM = 0.0;
  double piT = 0.0;
  
  return {piQ, piL, piM, piT};
}

     double PolarizationNonRel::CalculateBaseBeta1
(double q0, double q) const {
     double kminusinit;
     double kminusfin;
     double wtilde;
     double epsiloninitP1;
     double epsiloninitP2;
     double epsilonfinP1;
     double epsilonfinP2;
     double epsilonPartinit;
     wtilde=q0-(st.U4-st.U2);
     kminusinit=st.M2*wtilde/q-q/2.0;
     kminusfin=st.M4*wtilde/q-q/2.0;
     epsilonPartinit=1-2*(st.M4-st.M2)*kminusfin/q/st.M2;
    // epsiloninitP1=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1-sqrt(epsilonPartinit));
     double PI = o2scl_const::pi;
     double beta1;

     if (epsilonPartinit>0) {
     /*
     gsl_integration_workspace * w= gsl_integration_workspace_alloc (1024);
     double parinit[]={st.Mu2,st.T};
     gsl_function Finit;
     Finit.function = &FermiE;
     Finit.params = &parinit;
     epsiloninitP1=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1-sqrt(epsilonPartinit))*(1-sqrt(epsilonPartinit))+st.U2;
     epsiloninitP2=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1+sqrt(epsilonPartinit))*(1+sqrt(epsilonPartinit))+st.U2;
     epsilonfinP1=epsiloninitP1+wtilde+st.U4;
     epsilonfinP2=epsiloninitP2+wtilde+st.U4;
     double InteA;
     double result, error;
     size_t neval;
     gsl_integration_qng (&Finit,epsiloninitP1,epsiloninitP2,0,1e-7,&result,&error,&neval);
     InteA=result;

     gsl_integration_workspace_free (w);

     gsl_integration_workspace * wp= gsl_integration_workspace_alloc (1024);
     double parfin[]={st.Mu4,st.T,wtilde};
     gsl_function Fin;
     Fin.function = &FermiEwtilde;
     Fin.params = &parfin;
     double InteB;
     gsl_integration_qng (&Fin,epsilonfinP1,epsilonfinP2,0,1e-7,&result,&error,&neval);
     InteB=result;
     gsl_integration_workspace_free (wp);*/
     //below is analytical expression
     o2scl::polylog_multip<double , double> poly;
     epsiloninitP1=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1-sqrt(epsilonPartinit))*(1-sqrt(epsilonPartinit))+st.U2;
     epsiloninitP2=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1+sqrt(epsilonPartinit))*(1+sqrt(epsilonPartinit))+st.U2;
     epsilonfinP1=epsiloninitP1+wtilde+st.U4;
     epsilonfinP2=epsiloninitP2+wtilde+st.U4;
     double InteA;
     double InteB;
     double a;
     double b;
     double mu;
     double TT;
     a=epsiloninitP1;
     b=epsiloninitP2;
     mu=st.Mu2;
     TT=st.T;
     InteA=TT*(a*log(1 + exp((-a + mu)/TT)) - b*log(1 + exp((-b + mu)/TT)) - TT*poly.calc(2,-exp((-a + mu)/TT)) + TT*poly.calc(2,-exp((-b + mu)/TT)));
     
     a=epsilonfinP1;
     b=epsilonfinP2;
     mu=st.Mu4;
     InteB=-wtilde*TT*(log(1+exp(-(epsilonfinP1-st.Mu4)/st.T))-log(1+exp(-(epsilonfinP2-st.Mu4)/st.T)));
     InteB=InteB+TT*(a*log(1 + exp((-a + mu)/TT)) - b*log(1 + exp((-b + mu)/TT)) - TT*poly.calc(2,-exp((-a + mu)/TT)) + TT*poly.calc(2,-exp((-b + mu)/TT)));


     // Int_a^b 1/1+exp[(epsilon-mu)/T]depsilon = T*(Log[1+exp[-(a-mu)/T]]-Log[1+exp[-(b-mu)/T]])
     double InteC;
     double InteD;
     InteC= st.T*(log(1+exp(-(epsiloninitP1-st.Mu2)/st.T))-log(1+exp(-(epsiloninitP2-st.Mu2)/st.T)));

     InteD= st.T*(log(1+exp(-(epsilonfinP1-st.Mu4)/st.T))-log(1+exp(-(epsilonfinP2-st.Mu4)/st.T)));

     beta1=kminusfin/q*(InteC-InteD)+(st.M4-st.M2)/q/q*(InteA-InteB);
     beta1=beta1*st.M2*st.M4/(PI*q);
    // cout<<"epsiloninitP1, P2: "<<epsiloninitP1<<" "<<epsiloninitP2<<endl;
    // cout<<"epsilonfinP1, P2: "<<epsilonfinP1<<" "<<epsilonfinP2<<endl;
    // cout<<"Inte A, B, C, D: "<<InteA<<" "<<InteB<<" "<<InteC<<" "<<InteD<<endl;
     } else {beta1=0;}
  // cout<<"epsiloninitP1, P2: "<<epsiloninitP1<<" "<<epsiloninitP2<<endl;
  // cout<<"epsilonfinP1, P2: "<<epsilonfinP1<<" "<<epsilonfinP2<<endl;
  // cout<<"Inte A, B, C, D: "<<InteA<<" "<<InteB<<" "<<InteC<<" "<<InteD<<endl;
   return beta1;

}

 double PolarizationNonRel::CalculateBaseBeta2
(double q0, double q) const {
     double kminusinit;
     double kminusfin;
     double wtilde;
     double epsiloninitP1;
     double epsiloninitP2;
     double epsilonfinP1;
     double epsilonfinP2;
     double epsilonPartinit;
     wtilde=q0-(st.U4-st.U2);
     kminusinit=st.M2*wtilde/q-q/2.0;
     kminusfin=st.M4*wtilde/q-q/2.0;
     epsilonPartinit=1-2*(st.M4-st.M2)*kminusfin/q/st.M2;
    // epsiloninitP1=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1-sqrt(epsilonPartinit));
     double PI = o2scl_const::pi;
     double beta2;

     if (epsilonPartinit>0) {
     o2scl::polylog_multip<double , double> poly;
     epsiloninitP1=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1-sqrt(epsilonPartinit))*(1-sqrt(epsilonPartinit))+st.U2;
     epsiloninitP2=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1+sqrt(epsilonPartinit))*(1+sqrt(epsilonPartinit))+st.U2;
     epsilonfinP1=epsiloninitP1+wtilde+st.U4;
     epsilonfinP2=epsiloninitP2+wtilde+st.U4;
      double InteA;
      double InteB;
      double a;
      double b;
      double mu;
      double TT;
      a=epsiloninitP1;
      b=epsiloninitP2;
      mu=st.Mu2;
      TT=st.T;
      InteA=TT*(a*log(1 + exp((-a + mu)/TT)) - b*log(1 + exp((-b + mu)/TT)) - TT*poly.calc(2,-exp((-a + mu)/TT)) + TT*poly.calc(2,-exp((-b + mu)/TT)));

      a=epsilonfinP1;
      b=epsilonfinP2;
      mu=st.Mu4;
      InteB=-wtilde*TT*(log(1+exp(-(epsilonfinP1-st.Mu4)/st.T))-log(1+exp(-(epsilonfinP2-st.Mu4)/st.T)));
      InteB=InteB+TT*(a*log(1 + exp((-a + mu)/TT)) - b*log(1 + exp((-b + mu)/TT)) - TT*poly.calc(2,-exp((-a + mu)/TT)) + TT*poly.calc(2,-exp((-b + mu)/TT)));
      
      beta2=2*st.M2*(InteA-InteB)/q/q;
      beta2=beta2*st.M2*st.M4/q/PI;
      }else {beta2=0;}

     return beta2;


}

 double PolarizationNonRel::CalculateBaseBeta3
(double q0, double q) const {
     double kminusinit;
     double kminusfin;
     double wtilde;
     double epsiloninitP1;
     double epsiloninitP2;
     double epsilonfinP1;
     double epsilonfinP2;
     double epsilonPartinit;
     wtilde=q0-(st.U4-st.U2);
     kminusinit=st.M2*wtilde/q-q/2.0;
     kminusfin=st.M4*wtilde/q-q/2.0;
     epsilonPartinit=1-2*(st.M4-st.M2)*kminusfin/q/st.M2;
    // epsiloninitP1=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1-sqrt(epsilonPartinit));
     double PI = o2scl_const::pi;
     double beta3;
     if (epsilonPartinit>0) {
     o2scl::polylog_multip<double , double> poly;
     epsiloninitP1=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1-sqrt(epsilonPartinit))*(1-sqrt(epsilonPartinit))+st.U2;
     epsiloninitP2=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1+sqrt(epsilonPartinit))*(1+sqrt(epsilonPartinit))+st.U2;
     epsilonfinP1=epsiloninitP1+wtilde+st.U4;
     epsilonfinP2=epsiloninitP2+wtilde+st.U4;
      double InteA;
      double InteB;
      double a;
      double b;
      double mu;
      double TT;
      double ms1=st.M2;
      double ms2=st.M4;
      a=epsiloninitP1;
      b=epsiloninitP2;
      mu=st.Mu2;
      TT=st.T;
      InteA=TT*(pow(a*(ms1 - ms2) - kminusfin*q,2)*log(1 + exp((-a + mu)/TT)) - pow(b*(ms1 - ms2) - kminusfin*q,2)*log(1 + exp((-b + mu)/TT)) - 2*(ms1 - ms2)*(a*(ms1 - ms2) - kminusfin*q)*TT*poly.calc(2,-exp((-a + mu)/TT)) + 2*(ms1 - ms2)*(b*(ms1 - ms2) - kminusfin*q)*TT*poly.calc(2,-exp((-b + mu)/TT)) - 2*pow(ms1 - ms2,2)*pow(TT,2)*poly.calc(3,-exp((-a + mu)/TT)) + 2*pow(ms1 - ms2,2)*pow(TT,2)*poly.calc(3,-exp((-b + mu)/TT)));
      
      a=epsilonfinP1;
      b=epsilonfinP2;
      mu=st.Mu4;

      InteB=TT*(pow(kminusfin,2)*pow(q,2)*(-log(-2*exp((a - mu)/TT)) + log(-2*exp((b - mu)/TT)) + log(1 + exp((a - mu)/TT)) - log(1 + exp((b - mu)/TT))) + 2*kminusfin*(ms1 - ms2)*q*wtilde*(-log(-2*exp((a - mu)/TT)) + log(-2*exp((b - mu)/TT)) + log(1 + exp((a - mu)/TT)) - log(1 + exp((b - mu)/TT))) + pow(ms1 - ms2,2)*pow(wtilde,2)*(-log(-2*exp((a - mu)/TT)) + log(-2*exp((b - mu)/TT)) + log(1 + exp((a - mu)/TT)) - log(1 + exp((b - mu)/TT))) - 2*kminusfin*(ms1 - ms2)*q*(a*log(1 + exp((-a + mu)/TT)) - b*log(1 + exp((-b + mu)/TT)) - TT*poly.calc(2,-exp((-a + mu)/TT)) + TT*poly.calc(2,-exp((-b + mu)/TT))) - 2*pow(ms1 - ms2,2)*wtilde*(a*log(1 + exp((-a + mu)/TT)) - b*log(1 + exp((-b + mu)/TT)) - TT*poly.calc(2,-exp((-a + mu)/TT)) + TT*poly.calc(2,-exp((-b + mu)/TT))) + pow(ms1 - ms2,2)*(pow(a,2)*log(1 + exp((-a + mu)/TT)) - pow(b,2)*log(1 + exp((-b + mu)/TT)) - 2*a*TT*poly.calc(2,-exp((-a + mu)/TT)) + 2*b*TT*poly.calc(2,-exp((-b + mu)/TT)) - 2*pow(TT,2)*poly.calc(3,-exp((-a + mu)/TT)) + 2*pow(TT,2)*poly.calc(3,-exp((-b + mu)/TT))));

      beta3=(InteA-InteB)/q/q/q/q;
      beta3=st.M2*st.M4/PI/q;

     }else{beta3=0;}

     return beta3;

}  

    double PolarizationNonRel::CalculateBaseBeta4
(double q0, double q) const {
     double kminusinit;
     double kminusfin;
     double wtilde;
     double epsiloninitP1;
     double epsiloninitP2;
     double epsilonfinP1;
     double epsilonfinP2;
     double epsilonPartinit;
     wtilde=q0-(st.U4-st.U2);
     kminusinit=st.M2*wtilde/q-q/2.0;
     kminusfin=st.M4*wtilde/q-q/2.0;
     epsilonPartinit=1-2*(st.M4-st.M2)*kminusfin/q/st.M2;
    // epsiloninitP1=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1-sqrt(epsilonPartinit));
     double PI = o2scl_const::pi;
     double beta4;
     
     if (epsilonPartinit>0) {
     o2scl::polylog_multip<double , double> poly;
     epsiloninitP1=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1-sqrt(epsilonPartinit))*(1-sqrt(epsilonPartinit))+st.U2;
     epsiloninitP2=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1+sqrt(epsilonPartinit))*(1+sqrt(epsilonPartinit))+st.U2;
     epsilonfinP1=epsiloninitP1+wtilde+st.U4;
     epsilonfinP2=epsiloninitP2+wtilde+st.U4;
      double InteA;
      double InteB;
      double a;
      double b;
      double mu;
      double TT;
      double ms1=st.M2;
      double ms2=st.M4;
      a=epsiloninitP1;
      b=epsiloninitP2;
      mu=st.Mu2;
      TT=st.T;
      InteA=TT*(a*log(1 + exp((-a + mu)/TT)) - b*log(1 + exp((-b + mu)/TT)) - TT*poly.calc(2,-exp((-a + mu)/TT)) + TT*poly.calc(2,-exp((-b + mu)/TT)));

      a=epsilonfinP1;
      b=epsilonfinP2;
      mu=st.Mu4;
      InteB=-wtilde*TT*(log(1+exp(-(epsilonfinP1-st.Mu4)/st.T))-log(1+exp(-(epsilonfinP2-st.Mu4)/st.T)));
      InteB=InteB+TT*(a*log(1 + exp((-a + mu)/TT)) - b*log(1 + exp((-b + mu)/TT)) - TT*poly.calc(2,-exp((-a + mu)/TT)) + TT*poly.calc(2,-exp((-b + mu)/TT)));

      double InteC;
      double InteD;
      a=epsiloninitP1;
      b=epsiloninitP2;
      mu=st.Mu2;
      InteC=TT*(pow(a,2)*log(1 + exp((-a + mu)/TT)) - pow(b,2)*log(1 + exp((-b + mu)/TT)) + 2*TT*(-(a*poly.calc(2,-exp((-a + mu)/TT))) + b*poly.calc(2,-exp((-b + mu)/TT)) + TT*(-poly.calc(3,-exp((-a + mu)/TT)) + poly.calc(3,-exp((-b + mu)/TT)))));

      a=epsilonfinP1;
      b=epsilonfinP2;
      mu=st.Mu4;
      InteD=TT*(pow(a - wtilde,2)*log(1 + exp((-a + mu)/TT)) - pow(b - wtilde,2)*log(1 + exp((-b + mu)/TT)) + 2*TT*((-a + wtilde)*poly.calc(2,-exp((-a + mu)/TT)) + (b - wtilde)*poly.calc(2,-exp((-b + mu)/TT)) + TT*(-poly.calc(3,-exp((-a + mu)/TT)) + poly.calc(3,-exp((-b + mu)/TT)))));

      beta4=kminusfin/q*(InteA-InteB)+(st.M4-st.M2)/q/q*(InteC-InteD);
      beta4=beta4*2*st.M2/q/q*st.M2*st.M4/(PI*q);

     }else {beta4=0;}

      return beta4;

     }

     double PolarizationNonRel::CalculateBaseBeta5
(double q0, double q) const {
     double kminusinit;
     double kminusfin;
     double wtilde;
     double epsiloninitP1;
     double epsiloninitP2;
     double epsilonfinP1;
     double epsilonfinP2;
     double epsilonPartinit;
     wtilde=q0-(st.U4-st.U2);
     kminusinit=st.M2*wtilde/q-q/2.0;
     kminusfin=st.M4*wtilde/q-q/2.0;
     epsilonPartinit=1-2*(st.M4-st.M2)*kminusfin/q/st.M2;
    // epsiloninitP1=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1-sqrt(epsilonPartinit));
     double PI = o2scl_const::pi;
     double beta5;

     if (epsilonPartinit>0) {
     o2scl::polylog_multip<double , double> poly;
     epsiloninitP1=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1-sqrt(epsilonPartinit))*(1-sqrt(epsilonPartinit))+st.U2;
     epsiloninitP2=q*q*st.M2/2.0/(st.M4-st.M2)/(st.M4-st.M2)*(1+sqrt(epsilonPartinit))*(1+sqrt(epsilonPartinit))+st.U2;
     epsilonfinP1=epsiloninitP1+wtilde+st.U4;
     epsilonfinP2=epsiloninitP2+wtilde+st.U4;
      double InteA;
      double InteB;
      double a;
      double b;
      double mu;
      double TT;
      double ms1=st.M2;
      double ms2=st.M4;
      a=epsiloninitP1;
      b=epsiloninitP2;
      mu=st.Mu2;
      TT=st.T;
      InteA=TT*(pow(a,2)*log(1 + exp((-a + mu)/TT)) - pow(b,2)*log(1 + exp((-b + mu)/TT)) + 2*TT*(-(a*poly.calc(2,-exp((-a + mu)/TT))) + b*poly.calc(2,-exp((-b + mu)/TT)) + TT*(-poly.calc(3,-exp((-a + mu)/TT)) + poly.calc(3,-exp((-b + mu)/TT)))));

      a=epsilonfinP1;
      b=epsilonfinP2;
      mu=st.Mu4;
      InteB=TT*(pow(a - wtilde,2)*log(1 + exp((-a + mu)/TT)) - pow(b - wtilde,2)*log(1 + exp((-b + mu)/TT)) + 2*TT*((-a + wtilde)*poly.calc(2,-exp((-a + mu)/TT)) + (b - wtilde)*poly.calc(2,-exp((-b + mu)/TT)) + TT*(-poly.calc(3,-exp((-a + mu)/TT)) + poly.calc(3,-exp((-b + mu)/TT)))));

      beta5=4*st.M2*st.M2/q/q/q/q*(InteA-InteB);
      beta5=beta5*st.M2*st.M4/(PI*q);} else {beta5=0;}

      return beta5;



}


// only used when calculating NC mixed gas (n+p) axial part,neutron
// basic polarization use m2 and u2; proton basic polarization use m4
// and u4
std::array<double, 4> PolarizationNonRel::CalculateBasePolarizationsNeutron
(double q0, double q) const {
 
  //void Polarization::SetPolarizations(double q0, double q) {
  // Calculate some kinematic factors
   
  double q0t = q0 + st.U2 - st.U2;
  double qa2t = q0t*q0t - q*q;
  // I don't completely understand this condition, but it seems to be
  // necessary to suppress noise at larger q_0
 // if (qa2t > pow(st.M2 - st.M2, 2)*0.0) return {0.0, 0.0, 0.0, 0.0};
 // if (qa2t < 1.e-1 && qa2t >= 0.0) return {0.0, 0.0, 0.0, 0.0};
  double beta = 1.0 + (st.M2*st.M2 - st.M2*st.M4)/qa2t;
  double arg = beta*beta - 4.0*st.M2*st.M2/qa2t;
  if (arg<0.0) return {0.0, 0.0, 0.0, 0.0};
 // double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), st.M2);
 // double delta2 = (st.Mu2 - st.U2 - em)/st.T;
 // double delta4 = (st.Mu2 - st.U2 - em - q0t)/st.T;

  // following is the new em completely consistent with Redddy's
  // thesis, for non-rel+interacting gas
  double chi=1-st.M2/st.M2;
  double c=q0+st.U2-st.U2-q*q/(2*st.M2);
  //the minimum E2 for NC reaction
  double emNC=std::max((-c*st.M2/q)*(-c*st.M2/q)/(2*st.M2),0.0); 
 

  double delta2NC=(st.Mu2 - st.U2 - emNC)/st.T;
  double delta4NC=(st.Mu2 - st.U2 - emNC-q0t)/st.T;


  // Now just need to include some method for calculating these
  // At least at low density, Gamma0 should be the dominant term
  // which looks like the non-relativistic response function
  // Under non-degenerate conditions (i.e. delta2, delta4 << 0),
  // Gamma0 = Gamma1 = 0.5*Gamma2
  // This is exact
 // double Gamma0 = Fermi0(delta2) - Fermi0(delta4);

  // reddy nc current
  double xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);  


  double PI = o2scl_const::pi;//Constants::Pi;
  // double piL = st.M2*st.M2*st.T/(PI*q)*Gamma0;//orig one

  //neutral current consistent with Reddy's thesis
  double piL= st.M2*st.M2*st.T/(PI*q)*(xiNC+q0t/st.T);

  if (current==current_charged) {
    O2SCL_ERR("Invalid current 2.",o2scl::exc_efailed);
  }
  
  double piQ = 0.0;
  double piM = 0.0;
  double piT = 0.0;


  return {piQ, piL, piM, piT};
}

std::array<double, 4> PolarizationNonRel::CalculateBasePolarizationsProton
(double q0, double q) const {
  
  //void Polarization::SetPolarizations(double q0, double q) {
  // Calculate some kinematic factors
  
  double q0t = q0 + st.U4 - st.U4;
  double qa2t = q0t*q0t - q*q;
  // I don't completely understand this condition, but it seems to be
  // necessary to suppress noise at larger q_0
 // if (qa2t > pow(st.M4 - st.M4, 2)*0.0) return {0.0, 0.0, 0.0, 0.0};
 // if (qa2t < 1.e-1 && qa2t >= 0.0) return {0.0, 0.0, 0.0, 0.0};
  double beta = 1.0 + (st.M4*st.M4 - st.M4*st.M4)/qa2t;
  double arg = beta*beta - 4.0*st.M4*st.M4/qa2t;
  if (arg<0.0) return {0.0, 0.0, 0.0, 0.0};
 // double em = std::max(-0.5*beta*q0t + 0.5*q*sqrt(arg), st.M4);
 // double delta2 = (st.Mu4 - st.U4 - em)/st.T;
 // double delta4 = (st.Mu4 - st.U4 - em - q0t)/st.T;

  //following is the new em completely consistent with Redddy's
  //thesis, for non-rel+interacting gas
  
  double c=q0+st.U4-st.U4-q*q/(2*st.M4);
  //the minimum E2 for NC reaction
  double emNC=std::max((-c*st.M4/q)*(-c*st.M4/q)/(2*st.M4),0.0); 

  double delta2NC=(st.Mu4 - st.U4 - emNC)/st.T;
  double delta4NC=(st.Mu4 - st.U4 - emNC-q0t)/st.T;

  // Now just need to include some method for calculating these
  // At least at low density, Gamma0 should be the dominant term
  // which looks like the non-relativistic response function
  // Under non-degenerate conditions (i.e. delta2, delta4 << 0),
  // Gamma0 = Gamma1 = 0.5*Gamma2
  // This is exact
 // double Gamma0 = Fermi0(delta2) - Fermi0(delta4);

  // reddy nc current
  double xiNC=Fermi0(-delta2NC) - Fermi0(-delta4NC);  

  double PI = o2scl_const::pi;//Constants::Pi;
  // double piL = st.M2*st.M2*st.T/(PI*q)*Gamma0;//orig one

  //neutral current consistent with Reddy's thesis
  double piL= st.M4*st.M4*st.T/(PI*q)*(xiNC+q0t/st.T);

  if (current==current_charged) {
    O2SCL_ERR("Invalid current 3.",o2scl::exc_efailed);
  }
  
  double piQ = 0.0;
  double piM = 0.0;
  double piT = 0.0;

  return {piQ, piL, piM, piT};
}

double PolarizationNonRel::gamma0(double q0, double q) const {
  return 0;
}

double PolarizationNonRel::GetImPI( double q0, double q) const {

  // auto pt = CalculateBasePolarizations(q0, q);

  // double piL = pt[1];
  
  auto ptN=CalculateBasePolarizationsNeutron(q0, q);
  auto ptP=CalculateBasePolarizationsProton(q0, q);
  //only for NC mixed gas (N+P) axial channel, in vector channel there
  // is only contribution from neutron, just like pure neutron gas
  double piL=ptN[1]+ptP[1];
  // double piL=ptN[1];
  return piL;
}

double PolarizationNonRel::GetRePI( double q0, double q) const {

  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1024);
  gsl_function F;
       
  double par[]={q0,q,st.U2,st.U4,st.M2,st.M4,st.Mu2,st.Mu4,st.T};
  // double par[]={q0,q};
  F.function = &f;
  F.params = &par;
      
  double result, error;
  double q0boundaryCCLow, q0boundaryCCHigh;
  //from the kinematics of nucleons, the lowest q0=-(\delta U+\delta m)-q*(2*P_1-q)/(2*M)
  //where P_1 is the momentum of the initial nucleon. Approximately P_1=K_F^1 where K_F^1 is
  //the fermi momentum of the initial nucleon. given k_F, the maximum q*(2*P_1-q)/(2*M), corresponding
  //to the minimum q0, happens when q=k_F. Assuming in simulations, the highest k_F correspond to 10n0,
  //then k_F=713 MeV approximately. Then, the minimum q0=-(\delta U+\delta m)-713(2*713-713)/(2*939),
  //which is q0=-(\delta U+\delta m)-270MeV approximately. Conservatively, we set
  // q0lowThresh=-(\delta U+\delta m)-300MeV, to get RePi at q0=q0lowThresh, the integration range need
  // to be increased. Here I increase it by another 300 MeV

 // q0boundaryCCLow=st.U2-st.U4+st.M2-st.M4+300*2;
   // q0boundaryCCLow=100;
    q0boundaryCCLow=1.5*q+abs(st.U2-st.U4);
  //the upper bound of q0 can not be larger than E1, the largest E1 considered in simulations can be 480 mev, conservatively, 
  //we set the upper bound of q0 to be 500 mev
 // q0boundaryCCHigh=500;
   // q0boundaryCCHigh=100;
   //Large boundary will bring obvious numerical noise of the real piL. set q0bound=1.5q+abs(dU) here

    q0boundaryCCHigh=1.5*q+abs(st.U2-st.U4);
 /* if (10*st.T>100.0) {
          q0boundaryCC=30*st.T;
  } else {
          q0boundaryCC=100.0;
  }*/


  // gsl_integration_qawc (&F, -90, 90, q0, 0, 1e-7, 1000,
  //                          w, &result, &error);
 // gsl_integration_qawc (&F, -300, 100.0, q0, 0, 1e-7, 1000,
   //                     w, &result, &error);
  gsl_integration_qawc (&F, -q0boundaryCCLow, q0boundaryCCHigh, q0, 0, 1e-7, 1000,
                        w, &result, &error);
  // gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
  //                         w, &result, &error);


  // ******************************************************
  double piLRe=-result/o2scl_const::pi;//Constants::Pi;

  gsl_integration_workspace_free (w);

  return piLRe;
}

double PolarizationNonRel::GetRePIn( double q0, double q) const {


  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1024);
  gsl_function F;

  double par[]={q0,q,st.U2,st.U4,st.M2,st.M4,st.Mu2,st.Mu4,st.T};
  // double par[]={q0,q};
  F.function = &fn;
  F.params = &par;

  double result, error;
 /* double q0boundaryNC;
  if (10*st.T>100.0) {
          q0boundaryNC=10*st.T;
  } else {
          q0boundaryNC=100.0;
  }*/

   double q0boundaryNCLow, q0boundaryNCHigh;
  //from the kinematics of nucleons, the lowest q0=-(\delta U+\delta m)-q*(2*P_1-q)/(2*M)
  //where P_1 is the momentum of the initial nucleon. Approximately P_1=K_F^1 where K_F^1 is
  //the fermi momentum of the initial nucleon. given k_F, the maximum q*(2*P_1-q)/(2*M), corresponding
  //to the minimum q0, happens when q=k_F. Assuming in simulations, the highest k_F correspond to 10n0,
  //then k_F=713 MeV approximately. Then, the minimum q0=-(\delta U+\delta m)-713(2*713-713)/(2*939),
  //which is q0=-(\delta U+\delta m)-270MeV approximately. Conservatively, we set
  // q0lowThresh=-(\delta U+\delta m)-300MeV, to get RePi at q0=q0lowThresh, the integration range need
  // to be increased. Here I increase it by another 300 MeV

  //q0boundaryNCLow=st.U2-st.U2+st.M2-st.M2+300*2;
   q0boundaryNCLow=1.5*q;
  //the upper bound of q0 can not be larger than E1, the largest E1 considered in simulations can be 480 mev, conservatively,
  //we set the upper bound of q0 to be 500 mev
  //q0boundaryNCHigh=500;
  ////Large boundary will bring obvious numerical noise of the real piL. set q0bound=1.5q here

   q0boundaryNCHigh=1.5*q;


  // gsl_integration_qawc (&F, -90, 90, q0, 0, 1e-7, 1000,
  //                          w, &result, &error);
 // gsl_integration_qawc (&F, -100, 100.0, q0, 0, 1e-7, 1000,
  //                      w, &result, &error);//usually use this one
     gsl_integration_qawc (&F, -q0boundaryNCLow, q0boundaryNCHigh, q0, 0, 1e-7, 1000,
                        w, &result, &error);
  // gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
  //                         w, &result, &error);


  // ******************************************************
  double  piLRe=-result/o2scl_const::pi;//Constants::Pi;

  gsl_integration_workspace_free (w);

  return piLRe;
}

double PolarizationNonRel::GetRePIp( double q0, double q) const {


  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1024);
  gsl_function F;

  double par[]={q0,q,st.U2,st.U4,st.M2,st.M4,st.Mu2,st.Mu4,st.T};
  // double par[]={q0,q};
  
  F.function = &fp;
  F.params = &par;

  double result, error;
  
  // gsl_integration_qawc (&F, -90, 90, q0, 0, 1e-7, 1000,
  //                          w,&result,&error);
 /* double q0boundaryNC;
  if (10*st.T>100.0) {
	  q0boundaryNC=10*st.T;
  } else {
	  q0boundaryNC=100.0;
  }*/
   double q0boundaryNCLow, q0boundaryNCHigh;
  //from the kinematics of nucleons, the lowest q0=-(\delta U+\delta m)-q*(2*P_1-q)/(2*M)
  //where P_1 is the momentum of the initial nucleon. Approximately P_1=K_F^1 where K_F^1 is
  //the fermi momentum of the initial nucleon. given k_F, the maximum q*(2*P_1-q)/(2*M), corresponding
  //to the minimum q0, happens when q=k_F. Assuming in simulations, the highest k_F correspond to 10n0,
  //then k_F=713 MeV approximately. Then, the minimum q0=-(\delta U+\delta m)-713(2*713-713)/(2*939),
  //which is q0=-(\delta U+\delta m)-270MeV approximately. Conservatively, we set
  // q0lowThresh=-(\delta U+\delta m)-300MeV, to get RePi at q0=q0lowThresh, the integration range need
  // to be increased. Here I increase it by another 300 MeV

 // q0boundaryNCLow=st.U4-st.U4+st.M4-st.M4+300*2;
  q0boundaryNCLow=1.5*q;
  //the upper bound of q0 can not be larger than E1, the largest E1 considered in simulations can be 480 mev, conservatively,
  //we set the upper bound of q0 to be 500 mev
  //q0boundaryNCHigh=500;
  //Large boundary will bring obvious numerical noise of the real piL. set q0bound=1.5q here
  q0boundaryNCHigh=1.5*q;



//  gsl_integration_qawc(&F,-100,100.0,q0,0,1e-7,1000,
  //                     w,&result,&error);// usually use this one
  gsl_integration_qawc(&F,-q0boundaryNCLow,q0boundaryNCHigh,q0,0,1e-7,1000,
                       w,&result,&error);

  // gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
  //                         w, &result, &error);

  double piLRe=-result/o2scl_const::pi;//Constants::Pi;

  gsl_integration_workspace_free (w);

  return piLRe;
}



void PolarizationNonRel::PrintGetImPI2(double q0, double q) const {
  
  /* 
     double RepiL=0;
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
     }
  */
  
  return;
}


void PolarizationNonRel::SetPolarizations(double q0, double q,
                                          Tensor<double>* piVV, 
                                          Tensor<double>* piAA, 
                                          Tensor<double>* piTT, 
                                          Tensor<double>* piVA, 
                                          Tensor<double>* piVT, 
                                          Tensor<double>* piAT,
                                          bool pnm) const {
  
  if (pnm==true || current==current_neutral) {
    double piLn;
    double piLp;
    double piLnRe;
    double piLpRe;
    double piRPAvec;
    double piRPAax;
    double piL;
    SetPolarizations_neutral(q0,q,piVV,piAA,piTT,piVA,piVT,piAT,
                             piLn,piLp,piLnRe,piLpRe,
                             piRPAvec,piRPAax,piL,pnm);
  } else if (pnm==false) {
    double piLIm;
    double piLRe;
    double piRPAvec;
    double piRPAax;
    double piL;
    double Imbeta1;
    SetPolarizations_charged(q0,q,piVV,piAA,piTT,piVA,piVT,piAT,piLIm,
                             piLRe,piRPAvec,piRPAax,piL,Imbeta1);
  }

  return; 
}
 
void PolarizationNonRel::SetPolarizations_neutral
(double q0, double q,
 Tensor<double>* piVV, Tensor<double>* piAA, Tensor<double>* piTT, 
 Tensor<double>* piVA, Tensor<double>* piVT, Tensor<double>* piAT,
 double &piLn, double &piLp, double &piLnRe, double &piLpRe,
 double &piRPAvec, double &piRPAax, double &piL, bool pnm) const {

  // Calculate the basic parts of the polarization
  std::array<double,4> ptN=CalculateBasePolarizationsNeutron(q0,q);
  std::array<double,4> ptP;
  if (!pnm) {
    ptP=CalculateBasePolarizationsProton(q0,q);
  }
  
  piLn=ptN[1];

  if (!pnm) {
    piLp=ptP[1];
  }
  
  piLnRe=GetRePIn(q0,q);
  if (!pnm) {
    piLpRe=GetRePIp(q0,q);
  }
  
  double impin, impip, repin, repip;
  double pirpaVec, pirpaAx;

  if (pnm) {
    impin=piLn/2.0;
    impip=0.0;
    repin=piLnRe/2.0;
    repip=0.0;
  } else {
    impin=piLn/2.0;
    impip=piLp/2.0;
    repin=piLnRe/2.0;
    repip=piLpRe/2.0;
  }
  
  // Coulomb correction for fpp
  double e2=1.0/137.0*4.0*o2scl_const::pi;
  double qtf2=4.0*e2*cbrt(o2scl_const::pi)*
    pow(3.0*xn_proton*pow(o2scl_const::hc_mev_fm,3),2.0/3.0);
  double coulombf=e2*4.0*o2scl_const::pi/(q*q+qtf2);

  double xfppCoul=xfpp+coulombf;
  
  if (pnm) {
    
    // Vector polarization in pure neutron matter
    piRPAvec=impin/(impin*impin*xfnn*xfnn+pow(1.0-repin*xfnn,2.0));
    // Axial polarization in pure neutron matter
    piRPAax=impin/(impin*impin*xgnn*xgnn+pow(1.0-repin*xgnn,2.0));
    
  } else {
    
    // Vector polarization
    piRPAvec=(impin+
              xfnp*xfnp*impin*impin*impip+
              xfppCoul*xfppCoul*impin*impip*impip+
              xfnp*xfnp*impip*repin*repin-
              2.0*xfppCoul*impin*repip+
              xfppCoul*xfppCoul*impin*repip*repip)/
      ((-xfnn*impin-
        xfppCoul*impip-
        xfnp*xfnp*impip*repin+
        xfnn*xfppCoul*impip*repin-
        xfnp*xfnp*impin*repip+
        xfnn*xfppCoul*impin*repip)*
       (-xfnn*impin-
        xfppCoul*impip-
        xfnp*xfnp*impip*repin+
        xfnn*xfppCoul*impip*repin-
        xfnp*xfnp*impin*repip+
        xfnn*xfppCoul*impin*repip)+
       (1+
        xfnp*xfnp*impin*impip-
        xfnn*xfppCoul*impin*impip-
        xfnn*repin-
        xfppCoul*repip-
        xfnp*xfnp*repin*repip+
        xfnn*xfppCoul*repin*repip)*
       (1+
        xfnp*xfnp*impin*impip-
        xfnn*xfppCoul*impin*impip-
        xfnn*repin-
        xfppCoul*repip-
        xfnp*xfnp*repin*repip+
        xfnn*xfppCoul*repin*repip));
  
    // Axial polarization
    piRPAax=((xgnn+xgnp)*(xgnn+xgnp)*impin*impin*impip+impip*
             (-1.0+xgnn*repin+xgnp*repin)*
             (-1.0+xgnn*repin+xgnp*repin)+impin*
             (1.0-2.0*xgnp*repip-
              2.0*xgpp*repip+xgnp*xgnp*(impip*impip+repip*repip)+
              2.0*xgnp*xgpp*(impip*impip+repip*repip)+
              xgpp*xgpp*(impip*impip+repip*repip)))/
      ((-xgnn*impin-xgpp*impip-xgnp*xgnp*impip*repin+
        xgnn*xgpp*impip*repin-xgnp*xgnp*impin*repip+
        xgnn*xgpp*impin*repip)*
       (-xgnn*impin-xgpp*impip-xgnp*xgnp*impip*repin+
        xgnn*xgpp*impip*repin-xgnp*xgnp*impin*repip+
        xgnn*xgpp*impin*repip)+
       (1.0+xgnp*xgnp*impin*impip-xgnn*xgpp*impin*impip-
        xgnn*repin-xgpp*repip-xgnp*xgnp*repin*repip+
        xgnn*xgpp*repin*repip)*
       (1.0+xgnp*xgnp*impin*impip-xgnn*xgpp*impin*impip-
        xgnn*repin-xgpp*repip-xgnp*xgnp*repin*repip+
        xgnn*xgpp*repin*repip));
  }

  if (flag==flag_vector) {
    piL=2.0*piRPAvec;
  } else {
    piL=2.0*piRPAax;
  }

  // Set the different parts of the polarization 
  piVV->L=piL; 
  piAA->Tp=0.5*piL; 

  return;
}

void PolarizationNonRel::SetPolarizations_charged
(double q0, double q,
 Tensor<double>* piVV, Tensor<double>* piAA, Tensor<double>* piTT, 
 Tensor<double>* piVA, Tensor<double>* piVT, Tensor<double>* piAT, double &piLIm,
 double &piLRe, double &piRPAvec, double &piRPAax, double &piL, double &Imbeta1) const {

  // Calculate the basic parts of the polarization
  std::array<double, 4> pt;
  if (false) {
    pt=CalculateBasePolarizations(q0,q);
  } else {
    pt=base_polarization_new<>(q0,q);
  }
  
  piL=pt[1];
  piLIm=piL;//record the MF imaginary PIL
  Imbeta1=CalculateBaseBeta1(q0,q);
  
  piLRe=GetRePI(q0,q);

  piRPAvec=(piL/2)/((1-xvf*(piLRe/2))*(1-xvf*(piLRe/2))+
                      xvf*xvf*(piL/2)*(piL/2));
  
  piRPAax=(piL/2)/((1-xvgt*(piLRe/2))*(1-xvgt*(piLRe/2))+
                     xvgt*xvgt*(piL/2)*(piL/2));
  
  //piL is 2*Im PI, and it is now RPA piL
  if (flag==flag_vector) {
    piL=2*piRPAvec;
  } else {
    piL=2*piRPAax;
  }
  //test
  //cout<<"q0,q,piL: "<<q0<<" "<<q<<" "<<piL<<endl;  
  // Set the different parts of the polarization 
  piVV->L=piL; 
  piAA->Tp=0.5*piL;

  return;
}

void PolarizationNonRel::SetLeptonTensor(double E1, double q0, double q,
                                         Tensor<double>* L) const {
  double E3 = E1 - q0;
 // double mu13 = (q*q - E1*E1 - E3*E3)/(2.0*E1*E3);//original one
  double mu13 = -(q*q - E1*E1 - E3*E3)/(2.0*E1*E3);
  L->L  = 8.0*E1*E3*(1 + mu13); 
  L->Tp = 8.0*E1*E3*(3 - mu13);
  return;
}
  
  
