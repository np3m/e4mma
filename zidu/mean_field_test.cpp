#include <iostream> 
#include <math.h> 
#include <vector>
#include<fstream>

#include "Tensor.hpp"
#include "Polarization.hpp"
#include "PolarizationNonRel.hpp"
#include "Constants.hpp"


using namespace nuopac; 

int main() {
//*********************condition set*******************************************	
  double T = 10;
  double M2 = 939.0;
  double M4 = 937.7;
  double U2 = 0.0;
  double U4 = 0.0;
  double M3 = 0.0;
   
  
  const double densFac = pow(Constants::HBCFmMeV, 3);
 //****************Input File******************************************** 
  char outputnc[200];
  std::ofstream datanc;
  sprintf(outputnc, "%s", "ncOutputT10E12NonRelWholeDe_NRAPREoS_Reddy_RPAvectorCompleteVersion+CoulombF_P+N_test_v2test_dynamicAx+Vec_q3T_null.txt");
//  sprintf(outputnc, "%s", "ncOutputT10E12NonRel_test.txt");
  datanc.open(outputnc, std::ios::out);
  char outputcc[200];
  std::ofstream datacc;
  sprintf(outputcc, "%s", "ccOutputT10E12NonRelWholeDen_NRAPREoS_Reddy_RPAAxialResIntYpDependent_v2test_dynamicAx+Vec_q10T.txt");
 // sprintf(outputcc, "%s", "ccOutputT10E12NonRelAxial_wholeDen_test.txt");
  datacc.open(outputcc, std::ios::out);
  char trajfilename[]="skyrme_NRAPR_eos_Read.txt";
 // char trajfilename[]="ft_7_Read.txt";
 // char trajfilename[]="VirialEoS.txt";
  FILE*fpconfig;
  fpconfig=fopen(trajfilename,"r");
  //********************variable set***********************************************
  double n2;
  double n4;
  FluidState st;// for specified Y_e state
  FluidState beta;// for beta equillibrium state
  double E1;
  double sig0;
  double elastic;
  double elastic1;
  double piL;
  double piLRe;
    // directly read EoS data table**************************

double *ed;
ed=(double *) malloc (sizeof(double)*1000);
double *pr;
pr=(double *) malloc (sizeof(double)*1000);
double *nb;
nb=(double *) malloc (sizeof(double)*1000);

double *mu2eos;
mu2eos=(double *) malloc (sizeof(double)*1000);
double *mu4eos;
mu4eos=(double *) malloc (sizeof(double)*1000);
double *mu3eos;
mu3eos=(double *) malloc (sizeof(double)*1000);
double *u2eos;
u2eos=(double *) malloc (sizeof(double)*1000);
double *u4eos;
u4eos=(double *) malloc (sizeof(double)*1000);
double *n2eos;
n2eos=(double *) malloc (sizeof(double)*1000);
double *n4eos;
n4eos=(double *) malloc (sizeof(double)*1000);
double *n3eos;
n3eos=(double *) malloc (sizeof(double)*1000);
double *kfn;
kfn=(double *) malloc (sizeof(double)*1000);
double *kfp;
kfp=(double *) malloc (sizeof(double)*1000);
double *kfe;
kfe=(double *) malloc (sizeof(double)*1000);
double *fcs2;
fcs2=(double *) malloc (sizeof(double)*1000);
double *dednb_Ye, *dPdnb_Ye, *mumu, *nmu, *kfmu, *cs2, *logp, *loge, *s, *urca, *ad_index, *msn, *msp, *dnndnun, *dnpdnup, *dmundnn, *dmudn_mixed, *dmupdnp;
dednb_Ye=(double *) malloc (sizeof(double)*1000);
dPdnb_Ye=(double *) malloc (sizeof(double)*1000);
mumu=(double *) malloc (sizeof(double)*1000);
nmu=(double *) malloc (sizeof(double)*1000);
kfmu=(double *) malloc (sizeof(double)*1000);
cs2=(double *) malloc (sizeof(double)*1000);
logp=(double *) malloc (sizeof(double)*1000);
loge=(double *) malloc (sizeof(double)*1000);
s=(double *) malloc (sizeof(double)*1000);
urca=(double *) malloc (sizeof(double)*1000);
ad_index=(double *) malloc (sizeof(double)*1000);
msn=(double *) malloc (sizeof(double)*1000);
msp=(double *) malloc (sizeof(double)*1000);
dnndnun=(double *) malloc (sizeof(double)*1000);
dnpdnup=(double *) malloc (sizeof(double)*1000);
dmundnn=(double *) malloc (sizeof(double)*1000);
dmudn_mixed=(double *) malloc (sizeof(double)*1000);
dmupdnp=(double *) malloc (sizeof(double)*1000);

double mevFac=197.326;
//***reading nrapr EoS************
for (int i = 0; i < 100; i++)
                {
                        fscanf(fpconfig, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &ed[i], &pr[i], &nb[i], &mu2eos[i], &mu4eos[i], &mu3eos[i],  &n2eos[i],  &n4eos[i],  &n3eos[i], &kfn[i], &kfp[i], &kfe[i], &fcs2[i], &dednb_Ye[i], &dPdnb_Ye[i], &mumu[i], &nmu[i], &kfmu[i], &cs2[i], &logp[i], &loge[i], &s[i], &urca[i], &ad_index[i], &msn[i], &msp[i], &u2eos[i], &u4eos[i], &dnndnun[i], &dnpdnup[i], &dmundnn[i], &dmudn_mixed[i], &dmupdnp[i]);

                }
//****reading ft_x EoS*******************
/*for (int i = 0; i < 100; i++)
                {
                        fscanf(fpconfig, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &ed[i], &pr[i], &nb[i], &mu2eos[i], &mu4eos[i], &mu3eos[i],  &n2eos[i],  &n4eos[i],  &n3eos[i], &kfn[i], &kfp[i], &kfe[i], &cs2[i], &logp[i], &loge[i], &s[i], &urca[i], &ad_index[i], &msn[i], &msp[i], &u2eos[i], &u4eos[i], &dnndnun[i], &dnpdnup[i], &dmundnn[i], &dmudn_mixed[i], &dmupdnp[i]);

                }*/

//****reading virial EoS, the table reads: n2, n4[fm^-3], u2[MeV], u4, mu2[MeV], mu4,*******************
/*for (int i = 0; i < 100; i++)
                {
                        fscanf(fpconfig, "%lf %lf %lf %lf %lf %lf", &n2eos[i],  &n4eos[i], &u2eos[i], &u4eos[i],&mu2eos[i], &mu4eos[i]);
                }*/



/*for (int i=0; i<10; i++)
{
   std::cout <<" n2: "<<n2eos[i]<<" n4: "<<n4eos[i]<<" U2: "<<u2eos[i]<<" u4: "<<u4eos[i]<<" U2-U4: "<<(u2eos[i]-u4eos[i])*mevFac<<" mu2: "<<mu2eos[i]<<" mu4: "<<mu4eos[i]<<" mu2-mu4: "<<(mu2eos[i]-mu4eos[i])*mevFac<<" m2eff: "<<msn[i]<<" m4eff: "<<msp[i]<<std::endl;
}*/
//virial output below
for (int i=0; i<10; i++)
{
   std::cout <<" n2: "<<n2eos[i]<<" n4: "<<n4eos[i]<<" U2: "<<u2eos[i]<<" u4: "<<u4eos[i]<<" U2-U4: "<<(u2eos[i]-u4eos[i])<<" mu2: "<<mu2eos[i]<<" mu4: "<<mu4eos[i]<<" mu2-mu4: "<<(mu2eos[i]-mu4eos[i])<<std::endl;
}


//***********************************************************************************************************************************


  //**********************neutral current*************************************************
  for (int t=1;t<1;t++)
  {
   T = 10;
   M2 = 939.0;
  // M4 = 939.0;
   
   
   M4 = 937.7;//if use nrapr fluidstate with n+p mixture, or calculate pure proton
   U2 = 0.0;
   U4 = 0.0;
   M3 = 0.0;
  
  // n2 = 1.e-3*densFac*t;
  // n4 = 1.e-3*densFac*t;
  //****set nrapr eos*****************
    n2=n2eos[t]*densFac;
    n4=n4eos[t]*densFac;//if use nrapr
    u2eos[t]=mu2eos[t]*mevFac-(kfn[t]*mevFac)*(kfn[t]*mevFac)/(2.0*msn[t]*mevFac)-M2;//u2eos is with rest mass. now remove the rest mass and Ui=mui-mu_free, here mu_free is the fermi energy
    u4eos[t]=mu4eos[t]*mevFac-(kfp[t]*mevFac)*(kfp[t]*mevFac)/(2.0*msp[t]*mevFac)-M4;
 // *******************************************************************************
    FluidState betaEoS;
   // betaEoS=FluidState::StateFromDensities(T, M2, M4, n2, n4, U2, U4, M3, n4);
  // betaEoS=FluidState::StateFromDensities(T, M2, M4, n2, n4, u2eos[t], u4eos[t], M3, n4);
  // betaEoS=FluidState::StateFromDensities(T, msn[t]*mevFac, msp[t]*mevFac, n2, n4, u2eos[t], u4eos[t], M3, n4);
 
 // betaEoS=FluidState::StateFromDensities(T, msn[t]*mevFac, msn[t]*mevFac, n2, n2, 0.0, 0.0, M3, n4);//calculate the neutrino-neutron response in nrapr, forget proton, simply set n2=n4,u2=u4=0, and m2=m4
  betaEoS=FluidState::StateFromDensities(T, msn[t]*mevFac, msp[t]*mevFac, n2, n4, 1.0*u2eos[t], 1.0*u4eos[t], M3, n4);//for mixture gas n+p
 // betaEoS=FluidState::StateFromDensities(T, 0.9*M2, 0.9*M4, n2, n4, 1.0*u2eos[t], 1.0*u4eos[t], M3, n4);//for mixture gas n+p,test Effmass influence 
 // betaEoS=FluidState::StateFromDensities(T, msp[t]*mevFac, msp[t]*mevFac, n4, n4, 1.0*u4eos[t], 1.0*u4eos[t], M3, n2);//pure proton 
 // betaEoS=FluidState::StateFromDensities(T, msn[t]*mevFac, msn[t]*mevFac, n2, n2, 1.0*u2eos[t], 1.0*u2eos[t], M3, n2);//calculate the neutrino-neutron response in nrapr, forget proton, simply set n2=n4,u2=u4=0, and m2=m4, note that in Reddy NC, u2 and u4 has to be the one in the EoS, while if for nuopac NC calculation, u2 is only required to be equal to u4.
 // betaEoS=FluidState::StateFromDensities(T, M2, M2, n2, n2, 20.0*n2/(0.16*densFac), 20.0*n2/(0.16*densFac), M3, n2);//calculate the neutrino-neutron response in nuopac, forget proton, simply set n2=n4, and m2=m4
 //  betaEoS.Mu2=betaEoS.Mu2-msn[t]*mevFac;
 
 //  betaEoS.Mu4=betaEoS.Mu2-msn[t]*mevFac;// if use nuopac to determine mu, and reddy's non-rel to find response, the chemical potential should be tranformed to non-rel version, by substract ing effective mass from it
   
   betaEoS.Mu2=mu2eos[t]*mevFac-939.0;
  // betaEoS.Mu2=mu4eos[t]*mevFac-937.7;//for pure proton part
  // betaEoS.Mu4=mu2eos[t]*mevFac-939.0;// if use mu from nrapr, and reddy's non-rel to find response, the chemical potential should be tranformed to non-rel version, by substracting bare  mass from it
   betaEoS.Mu4=mu4eos[t]*mevFac-937.7;//for mixture gas n+p, and for proton part 

  // betaEoS.Mu2=betaEoS.Mu2-939.0;
  // betaEoS.Mu4=betaEoS.Mu4-939.0;// if use nuopac with bare mass EoS and its mu, and reddy's non-rel to find response, the chemical potential should be tranformed to non-rel version, by substracting bare mass from it
 
  //*****Virial**************************
  /* betaEoS=FluidState::StateFromDensities(T, M2, M4, n2eos[t]*densFac, n4eos[t]*densFac, u2eos[t], u4eos[t], M3);//for mixture gas n+p
   betaEoS.Mu2=mu2eos[t];
   betaEoS.Mu4=mu4eos[t];*/
  //*****************************************

 // FluidState st = FluidState::StateFromDensities(T, M2, M4, n2, n4, 20.0*n2/(0.16*densFac), -20.0*n2/(0.16*densFac), M3, n2);
// st = FluidState::StateFromDensities(T, M2, M4, n2, n4, 20.0*n2/(0.16*densFac), 20.0*n2/(0.16*densFac), M3, n2);
 
 // FluidState stInverse(T, M2, M4, st.Mu2, st.Mu4, 20.0*n2/(0.16*densFac),  -20.0*n2/(0.16*densFac), M3, st.Mu3);
// FluidState stInverse(T, M2, M4, st.Mu2, st.Mu4, 20.0*n2/(0.16*densFac), 20.0*n2/(0.16*densFac), M3, st.Mu3);
// FluidState st = FluidState::StateFromDensities(T, M2, M4, n2, n4, 0.0*n2/(0.16*densFac), -0.0*n2/(0.16*densFac), M3, n2);

/* FluidState st = FluidState::StateFromDensities(T,  M2, M4, n2, n4, 0.0*n2/(0.16*densFac), -0.0*n2/(0.16*densFac), M3, n2);
 FluidState stInverse(T, M2, M4, st.Mu2, st.Mu4, 0.0*n2/(0.16*densFac),  -0.0*n2/(0.16*densFac), M3, st.Mu3);
  if (fabs(1.0 - stInverse.n2/n2)>1.e-3) {
    std::cerr << "Bad density find in n2: " << n2
    << " " << stInverse.n2 << std::endl;
    return 1;
  }
  if (fabs(1.0 - stInverse.n4/n4)>1.e-3) {
    std::cerr << "Bad density find in n4: " << n4
    << " " << stInverse.n4 << std::endl;
    return 1;
  }*/

 
  WeakCouplings nscat = WeakCouplings::NeutronScattering();
  nscat.F2 = 0.0;
 // PolarizationNonRel pol(beta, nscat, false);
 // PolarizationNonRel pol(st, nscat, false);
 // PolarizationNonRel pol(st, nscat, false);
  PolarizationNonRel pol(betaEoS, nscat, false, false, false);//turn off pauli blocking for NC
  E1 = 12.0;
  double full = pol.CalculateInverseMFP(E1)/Constants::HBCFmMeV*1.e13;


  sig0 =
      4.0*pow(Constants::GfMeV*E1*Constants::HBCFmMeV, 2)/Constants::Pi;
   // sig0 = sig0*(pow(nscat.Cv, 2) + 3.0*pow(nscat.Ca, 2))/4.0;
  sig0 = sig0*(1.0*pow(nscat.Cv, 2) + 0.0*pow(nscat.Ca, 2))/4.0;
   // double elastic = sig0*beta.effectiveDensity/densFac*1.e13;
   // double elastic1 = sig0*beta.n2/densFac*1.e13;
 // betaEoS.n4=betaEoS.n2;//to make elastic calculation work
  elastic = sig0*betaEoS.effectiveDensity/densFac*1.e13;//for mixture gas state with n2 not equal to n2, it doesn't work
  elastic1 = sig0*betaEoS.n2/densFac*1.e13;//for mixture gas state with n2 not equal to n2, it doesn't work
  
 // piL=pol.GetImPI(-30,30);
  
  
 // piLRe=pol.GetImPI2(0, 30);
 

//**************test rpa ph interactions*************************************
   double t0,t1,t2,t3,x0,x1,x2,x3,epsilon;
   //ft6
     /*  epsilon= 1.4339134197*1.0E-01;
       t0 = -3.1535549664*1.0E+03;
       t1=2.9811820280*1.0E+02;
       t2=2.3615846306*1.0E+03;
       t3=1.9373593217*1.0E+04;
       x0=1.5139360999*1.0E-01;
       x1=-6.5192138295 *1.0E+00;
       x2=-1.3213613298 *1.0E+00;
       x3=5.8876078909*1.0E-01;*/
    //NRAPR
       epsilon= 1.4416*1.0E-01;
       t0 = -2.7197*1.0E+03;
       t1=4.1764*1.0E+02;
       t2= -6.6687*1.0E+01;
       t3= 1.5042*1.0E+04;
       x0= 1.6154*1.0E-01;
       x1= -4.7986*1.0E-02;
       x2= 2.717*1.0E-02;
       x3= 1.3611*1.0E-01;

     double rou=(betaEoS.n2+betaEoS.n4)/densFac;
     double roun=betaEoS.n2/densFac;
     double roup=betaEoS.n4/densFac;
      double kf=pow(0.5*rou*3*3.14159*3.14159,1.0/3.0);
     double kfnTest=pow(roun*3*3.14159*3.14159,1.0/3.0);
     double kfpTest=pow(roup*3*3.14159*3.14159,1.0/3.0);
      double fnn=0.5*(t0*(1.0-x0)+1.0/6.0*t3*pow(rou,epsilon)*(1.0-x3)+2.0/3.0*epsilon*t3*pow(rou,epsilon-1)*((1+x3/2.0)*rou-(1.0/2.0+x3)*roun)+1.0/6.0*epsilon*(epsilon-1.0)*t3*pow(rou,epsilon-2.0)*((1+x3/2.0)*pow(rou,2.0)-(0.5+x3)*(roun*roun+roup*roup)))+0.25*(t1*(1-x1)+3*t2*(1+x2))*kfnTest*kfnTest;
     double fpp=0.5*(t0*(1.0-x0)+1.0/6.0*t3*pow(rou,epsilon)*(1.0-x3)+2.0/3.0*epsilon*t3*pow(rou,epsilon-1)*((1+x3/2.0)*rou-(1.0/2.0+x3)*roup)+1.0/6.0*epsilon*(epsilon-1.0)*t3*pow(rou,epsilon-2.0)*((1+x3/2.0)*pow(rou,2.0)-(0.5+x3)*(roun*roun+roup*roup)))+0.25*(t1*(1-x1)+3*t2*(1+x2))*kfpTest*kfpTest;
     double gnn=0.5*(t0*(x0-1)+1.0/6.0*t3*pow(rou,epsilon)*(x3-1.0))+0.25*(t1*(x1-1)+t2*(1+x2))*kfnTest*kfnTest;
     double gpp=0.5*(t0*(x0-1)+1.0/6.0*t3*pow(rou,epsilon)*(x3-1.0))+0.25*(t1*(x1-1)+t2*(1+x2))*kfpTest*kfpTest;

     fnn=fnn/pow(197.3,3);
     fpp=fpp/pow(197.3,3);
     gnn=gnn/pow(197.3,3);
     gpp=gpp/pow(197.3,3);

  /*  double fnn;
    double fpp;
    double gnn;
    double gpp;
    fnn=fpp=-0.0000844949;
    gnn=gpp=0.0000759672;
  double f0p=0.000118591;
  double g0p=0.0000730909;
  double vf=2.0*f0p;
  double vgt=2.0*g0p;
  double vrpa=vf;//already in unit of MeV-2*/



//***********conclusion*******************************

 // piLRe=pol.GetImPI2(wmin+1.0E-6, 30);
 
 // std::cout << "Neutral current: "  <<st.n2/densFac<<" full: "<< full << " elastic: " << elastic << " elastic1: "
   //   <<elastic1<<" full/elastic:"<< full/elastic<< " full/elastic1: "<< full/elastic1<<" M3 " <<st.M3<<" piL "<< piL<<" piLRPA "<<piLRe<<" w "<<wmin<<std::endl;
 // datanc << "Neutral current: "  <<st.n2/densFac<<" full: "<< full << " elastic: " << elastic << " elastic1: "
 //     <<elastic1<<" full/elastic: "<< full/elastic<< " full/elastic1: "<< full/elastic1 << std::endl;
 // std::cout <<"neutral current: "<<(st.n2)/densFac<<" full: "<< full  << " elastic: "<<elastic<<" elastic1: "<<elastic1<<" full/elastc: " << full/elastic << " full/elastic1: "<< full/elastic1<<" M3: "<<st.M3<<" deltaU: "<<st.U2-st.U4<<" mu3 "<<st.Mu3<< " mu2: "<<st.Mu2<<" mu4: "<<st.Mu4<<std::endl;
 // datanc <<"neutral current: "<<(st.n2)/densFac<<" fullAxialPart: "<< full  << " elastic: "<<elastic<<" elastic1: "<<elastic1<<" full/elastc: " << full/elastic << " full/elastic1: "<< full/elastic1<<" M3: "<<st.M3<<" deltaU: "<<st.U2-st.U4<<" mu3 "<<st.Mu3<<" mu2: "<<st.Mu2<<" Mu4: "<<st.Mu4<< std::endl;
  std::cout <<"neutral current: "<<(betaEoS.n2)/densFac<<" "<<(betaEoS.n3)/densFac<<" "<<(betaEoS.n4)/densFac<<" full: "<< full  << " elastic: "<<elastic<<" elastic1: "<<elastic1<<" full/elastc: " << full/elastic << " full/elastic1: "<< full/elastic1<<" M3: "<<betaEoS.M3<<" deltaU: "<<betaEoS.U2-betaEoS.U4<<" U2: "<<betaEoS.U2<<" U4: "<<betaEoS.U4<<" mu3 "<<betaEoS.Mu3<< " mu2: "<<betaEoS.Mu2<<" mu4: "<<betaEoS.Mu4<<" fnn: "<<fnn<<" fpp: "<<fpp<<" gnn: "<<gnn<<" gpp: "<<gpp<<std::endl;
  datanc <<"neutral current: "<<(betaEoS.n2)/densFac<<" "<<(betaEoS.n3)/densFac<<" "<<(betaEoS.n4)/densFac<<" full: "<< full  << " elastic: "<<elastic<<" elastic1: "<<elastic1<<" full/elastc: " << full/elastic << " full/elastic1: "<< full/elastic1<<" M3: "<<betaEoS.M3<<" deltaU: "<<betaEoS.U2-betaEoS.U4<<" mu3 "<<betaEoS.Mu3<<" mu2: "<<betaEoS.Mu2<<" Mu4: "<<betaEoS.Mu4<<" fnn: "<<fnn<<" fpp: "<<fpp<<" gnn: "<<gnn<<" gpp: "<<gpp<< std::endl;
 //********************detailed q0 or cos dependence*****************************************
 //**********Set integration range***********
  double vel=sqrt((3.0*T)/betaEoS.M2);
  double wmin;
  double wmax;
  wmin=-3.0*vel*3*T-0.00000789;
 // wmin=-50.0+1.0e-3;
  wmax=3.0*(vel*3*T+3*T*3*T/(2.0*betaEoS.M2));
 // wmax=50.0+1.0e-3;
  double dw=(wmax-wmin)/100;

 
  for (int k=1;k<100;k++)
  {
 
  double w=wmin+dw*k;
 // double w=30.0;
 // piL=pol.GetImPI(w,3*T);//piL is 2*ImPI0
 // piLRe=pol.GetImPI2(w, 3*T);//piLRe is 2*RePI0
 //**********new version**************
    auto ptN= pol.CalculateBasePolarizationsNeutron(w, 3*T);//neutral current
    auto ptP= pol.CalculateBasePolarizationsProton(w, 3*T); //neutral current
    double piLn=ptN[1];
    double piLp=ptP[1];
 //***********PiLRe**********************
    double piLnRe=pol.GetRePIn(w,3*T);
    double piLpRe=pol.GetRePIp(w,3*T);   

 
 //************************************
 
 // double vf=(-0.74+2.5)*1.0E-5;
 // double vgt=4.5*1.0E-5;
 // double vrpa=vf;
 
 
   
 //
// double t0,t1,t2,t3,x0,x1,x2,x3,epsilon;
 //ft6
     /*  epsilon= 1.4339134197*1.0E-01;
       t0 = -3.1535549664*1.0E+03;
       t1=2.9811820280*1.0E+02;
       t2=2.3615846306*1.0E+03;
       t3=1.9373593217*1.0E+04;
       x0=1.5139360999*1.0E-01;
       x1=-6.5192138295 *1.0E+00;
       x2=-1.3213613298 *1.0E+00;
       x3=5.8876078909*1.0E-01;*/
  //NRAPR
       epsilon= 1.4416*1.0E-01;
       t0 = -2.7197*1.0E+03;
       t1=4.1764*1.0E+02;
       t2= -6.6687*1.0E+01;
       t3= 1.5042*1.0E+04;
       x0= 1.6154*1.0E-01;
       x1= -4.7986*1.0E-02;
       x2= 2.717*1.0E-02;
       x3= 1.3611*1.0E-01;

      rou=(betaEoS.n2+betaEoS.n4)/densFac;
      kf=pow(0.5*rou*3*3.14159*3.14159,1.0/3.0);
    // double kfn=pow(rou*3*3.14159*3.14159,1.0/3.0);


 //for high density skyrme
     
     double f0=3.0/4.0*t0+3.0/8.0*t1*kf*kf+5.0/8.0*t2*kf*kf+1.0/2.0*t2*x2*kf*kf+(epsilon+1.0)*(epsilon+2.0)/16.0*t3*pow(rou,epsilon);
     double f0p=-1.0/4.0*t0-1.0/2.0*t0*x0-1.0/8.0*t1*kf*kf-1.0/4.0*t1*x1*kf*kf+1.0/8.0*t2*kf*kf+1.0/4.0*t2*x2*kf*kf-1.0/24.0*t3*pow(rou,epsilon)-1.0/12.0*t3*x3*pow(rou,epsilon);
     double g0=-1.0/4.0*t0+1.0/2.0*x0*t0-1.0/8.0*t1*kf*kf+1.0/4.0*t1*x1*kf*kf+1.0/8.0*t2*kf*kf+1.0/4.0*t2*x2*kf*kf-1.0/24.0*t3*pow(rou,epsilon)+1.0/12.0*t3*x3*pow(rou,epsilon);
     double g0p=-1.0/4.0*t0-1.0/8.0*t1*kf*kf+1.0/8.0*t2*kf*kf-1.0/24.0*t3*pow(rou,epsilon);
    // double f0n=1.0/2.0*t0-1.0/2.0*t0*x0+1.0/8.0*t1*kfn*kfn-1.0/4.0*t1*x1*kfn*kfn+3.0/4.0*t2*kfn*kfn+3.0/4.0*t2*x2*kfn*kfn+(epsilon+1.0)*(epsilon+2.0)/24.0*t3*pow(rou,epsilon)-(epsilon+1.0)*(epsilon+2.0)/24.0*t3*x3*pow(rou,epsilon);
    // double g0n=-1.0/2.0*t0+1.0/2.0*t0*x0-1.0/4.0*t1*kfn*kfn+1.0/4.0*t1*x1*kfn*kfn+1.0/4.0*t2*kfn*kfn+1.0/4.0*t2*x2*kfn*kfn-1.0/12.0*t3*pow(rou,epsilon)+1.0/12.0*t3*x3*pow(rou,epsilon);
      fnn=0.5*(t0*(1.0-x0)+1.0/6.0*t3*pow(rou,epsilon)*(1.0-x3)+2.0/3.0*epsilon*t3*pow(rou,epsilon-1)*((1+x3/2.0)*rou-(1.0/2.0+x3)*roun)+1.0/6.0*epsilon*(epsilon-1.0)*t3*pow(rou,epsilon-2.0)*((1+x3/2.0)*pow(rou,2.0)-(0.5+x3)*(roun*roun+roup*roup)))+0.25*(t1*(1-x1)+3*t2*(1+x2))*kfnTest*kfnTest;
      fpp=0.5*(t0*(1.0-x0)+1.0/6.0*t3*pow(rou,epsilon)*(1.0-x3)+2.0/3.0*epsilon*t3*pow(rou,epsilon-1)*((1+x3/2.0)*rou-(1.0/2.0+x3)*roup)+1.0/6.0*epsilon*(epsilon-1.0)*t3*pow(rou,epsilon-2.0)*((1+x3/2.0)*pow(rou,2.0)-(0.5+x3)*(roun*roun+roup*roup)))+0.25*(t1*(1-x1)+3*t2*(1+x2))*kfpTest*kfpTest;
      gnn=0.5*(t0*(x0-1)+1.0/6.0*t3*pow(rou,epsilon)*(x3-1.0))+0.25*(t1*(x1-1)+t2*(1+x2))*kfnTest*kfnTest;
      gpp=0.5*(t0*(x0-1)+1.0/6.0*t3*pow(rou,epsilon)*(x3-1.0))+0.25*(t1*(x1-1)+t2*(1+x2))*kfpTest*kfpTest;
     double fnp=0.5*(t0*(2.0+x0)+1.0/6.0*t3*pow(rou,epsilon)*(2.0+x3)+1.0/2.0*epsilon*t3*pow(rou,epsilon)+1.0/6.0*epsilon*(epsilon-1.0)*t3*pow(rou,epsilon-2.0)*((1+x3/2.0)*pow(rou,2.0)-(0.5+x3)*(roun*roun+roup*roup)))+0.5*0.25*(t1*(2.0+x1)+t2*(2.0+x2))*(kfnTest*kfnTest+kfpTest*kfpTest);
     double gnp=0.5*(t0*x0+1.0/6.0*t3*pow(rou,epsilon)*x3)+0.5*0.25*(t1*x1+t2*x2)*(kfnTest*kfnTest+kfpTest*kfpTest); 
     fnn=fnn/pow(197.3,3);
     fpp=fpp/pow(197.3,3);
     gnn=gnn/pow(197.3,3);
     gpp=gpp/pow(197.3,3);
     fnp=fnp/pow(197.3,3);
     gnp=gnp/pow(197.3,3);


 //for low density virial 
 /*   double fnn;
    double fpp;
    double gnn;
    double gpp;
    double fnp;
    double gnp;
    fnn=fpp=-0.0000844949;
    gnn=gpp=0.0000759672;
    fnp=-0.000321676;
    gnp=-0.0000702146;
  double f0p=0.000118591;
  double g0p=0.0000730909;
   
  double vf=2.0*f0p;
  double vgt=2.0*g0p;
  double vrpa=vf;//already in unit of MeV-2*/

   



  double piLRPA;
  double piLRPAvec;
  double response;
  double FermiF;
  double zz;
  /***********test************************
  double q0t = w + st.U2 - st.U4;
  double q=3*T;
  double qa2t = q0t*q0t - q*q;
  // I don't completely understand this condition, but it seems to be necessary
  // to suppress noise at larger q_0
 // if (qa2t > pow(st.M2 - st.M4, 2)*0.0) return {0.0, 0.0, 0.0, 0.0};
 // if (qa2t < 1.e-1 && qa2t > 0.0) return {0.0, 0.0, 0.0, 0.0};
  double bet = 1.0 + (st.M2*st.M2 - st.M4*st.M4)/qa2t;
  double arg = bet*bet - 4.0*st.M2*st.M2/qa2t;
 // if (arg<0.0) return {0.0, 0.0, 0.0, 0.0};
  double em = std::max(-0.5*bet*q0t + 0.5*q*sqrt(arg), st.M2);
  double delta2 = (st.Mu2 - st.U2 - em)/st.T;
  double delta4 = (st.Mu4 - st.U4 - em - q0t)/st.T;

  // Now just need to include some method for calculating these
  // At least at low density, Gamma0 should be the dominant term
  // which looks like the non-relativistic response function
  // Under non-degenerate conditions (i.e. delta2, delta4 << 0),
  // Gamma0 = Gamma1 = 0.5*Gamma2
  // This is exact
  double Fermi0delta2;
    if (delta2>40.0){Fermi0delta2=delta2;}
    else {Fermi0delta2=log(exp(delta2) + 1.0);}

  double Fermi0delta4;
    if (delta4>40.0){Fermi0delta4=delta4;}
    else {Fermi0delta4=log(exp(delta4) + 1.0);}
    
  

  double Gamma0 =Fermi0delta2-Fermi0delta4;*/

  //*************************************
  zz=(w+0.0)/T;//neutral current
  FermiF=1/(1-exp(-zz));
 // piLRPA=2*(piL/2)/((1-vrpa*(piLRe/2))*(1-vrpa*(piLRe/2))+vrpa*vrpa*(piL/2)*(piL/2))*FermiF;
 // piLRPA=2*(piL/2)*FermiF;
 // piLRPA=2*((piLn+piLp)/2+(gnn-gpp)*(piLn/2.0*piLpRe/2.0-piLp/2.0*piLnRe/2.0))/((1-gnn*(piLnRe/2)-gpp*(piLpRe/2))*(1-gnn*(piLnRe/2)-gpp*(piLpRe/2))+(piLn/2.0*gnn+piLp/2.0*gpp)*(piLn/2.0*gnn+piLp/2.0*gpp))*FermiF;
 // piLRPAvec=2*(piLn/2)/((1-fnn*(piLnRe/2))*(1-fnn*(piLnRe/2))+(piLn/2.0*fnn)*(piLn/2.0*fnn))*FermiF;
  //piLRPA is s(q0,q)
  //
  //
  //*********************complete version of neutral current n+p gas polarization functions**********************
  double impin,impip,repin,repip;
  double pirpaVec,pirpaAx;
  impin=piLn/2.0;
  impip=piLp/2.0;
  repin=piLnRe/2.0;
  repip=piLpRe/2.0;
  //vector polarization complete version
  //********adding coulomb force in fpp only for NC vector part****************
  double e2,qtf2;
  double piconst;
  double coulombf;
  e2=1.0/137.0*4.0*piconst;
  piconst=3.1415926;
  qtf2=4.0*e2*pow(piconst,0.333333)*pow(3.0*rou*densFac,2.0/3.0);
  coulombf=e2*4.0*piconst/((3.0*T)*(3.0*T)+qtf2);

  fpp=fpp+coulombf;


  //**************************
  pirpaVec=(impin+fnp*fnp*impin*impin*impip+fpp*fpp*impin*impip*impip+fnp*fnp*impip*repin*repin-2.0*fpp*impin*repip+fpp*fpp*impin*repip*repip)/((-fnn*impin-fpp*impip-fnp*fnp*impip*repin+fnn*fpp*impip*repin-fnp*fnp*impin*repip+fnn*fpp*impin*repip)*(-fnn*impin-fpp*impip-fnp*fnp*impip*repin+fnn*fpp*impip*repin-fnp*fnp*impin*repip+fnn*fpp*impin*repip)+(1+fnp*fnp*impin*impip-fnn*fpp*impin*impip-fnn*repin-fpp*repip-fnp*fnp*repin*repip+fnn*fpp*repin*repip)*(1+fnp*fnp*impin*impip-fnn*fpp*impin*impip-fnn*repin-fpp*repip-fnp*fnp*repin*repip+fnn*fpp*repin*repip));
//
//
//
//axial polarization complete version
  pirpaAx=((gnn+gnp)*(gnn+gnp)*impin*impin*impip+impip*(-1.0+gnn*repin+gnp*repin)*(-1.0+gnn*repin+gnp*repin)+impin*(1.0-2.0*gnp*repip-2.0*gpp*repip+gnp*gnp*(impip*impip+repip*repip)+2.0*gnp*gpp*(impip*impip+repip*repip)+gpp*gpp*(impip*impip+repip*repip)))/((-gnn*impin-gpp*impip-gnp*gnp*impip*repin+gnn*gpp*impip*repin-gnp*gnp*impin*repip+gnn*gpp*impin*repip)*(-gnn*impin-gpp*impip-gnp*gnp*impip*repin+gnn*gpp*impip*repin-gnp*gnp*impin*repip+gnn*gpp*impin*repip)+(1.0+gnp*gnp*impin*impip-gnn*gpp*impin*impip-gnn*repin-gpp*repip-gnp*gnp*repin*repip+gnn*gpp*repin*repip)*(1.0+gnp*gnp*impin*impip-gnn*gpp*impin*impip-gnn*repin-gpp*repip-gnp*gnp*repin*repip+gnn*gpp*repin*repip));

 piLRPA=2.0*pirpaAx*FermiF;
 piLRPAvec=2.0*pirpaVec*FermiF;


 // std::cout << "neutral current density: " <<betaEoS.n2/densFac<<" vrpa: "<<vrpa<<" full "<< full <<" w "<<w<< " piL(w) " <<piL<<" piLRe(w) "<<piLRe<<" piLRPA "<<piLRPA<<" response "<<response<<" FermiF: "<<FermiF<< std::endl;
 // datanc << "neutral current density: " <<betaEoS.n2/densFac<<" full "<< full <<" vrpa: "<<vrpa<< " w "<<w<< " piL(w) " <<piL<<" piLRe(w) "<<piLRe<<" piLRPA "<<piLRPA<<" response "<<response<<" FermiF "<<FermiF<< std::endl;
 //
 // new version******
 //   std::cout << "neutral current density: " <<betaEoS.n2/densFac<<" fnn: "<<fnn<<" fpp: "<<fpp<<" gnn: "<<gnn<<" gpp: "<<gpp<<" full "<< full <<" w "<<w<< " piLn(w) " <<piLn<<" piLp(w) " <<piLp<<" piLRen(w) "<<piLnRe<<" piLRep(w) "<<piLpRe<<" piLRPAaxial "<<piLRPA<<" piLRPAvector "<<piLRPAvec<<" FermiF: "<<FermiF<< std::endl;
  datanc<< "neutral current density: " <<betaEoS.n2/densFac<<" fnn: "<<fnn<<" fpp: "<<fpp<<" gnn: "<<gnn<<" gpp: "<<gpp<<" full "<< full <<" w "<<w<< " piLn(w) " <<piLn<<" piLp(w) " <<piLp<<" piLRen(w) "<<piLnRe<<" piLRep(w) "<<piLpRe<<" piLRPA "<<piLRPA<<" piLRPAvector "<<piLRPAvec<<" FermiF: "<<FermiF<< std::endl;
  }

 // output dgammadq0*************************
  
 // double estar = betaEoS.M4 + betaEoS.U2 - betaEoS.M2 - betaEoS.U2;
  double estar=0.0;//neutral current
  estar = std::min(estar, E1 - betaEoS.M3);
  double integral = 0.0;
/*  for (int i=0; i<pol.NNPGL; ++i) {
    double q0 = estar - pol.xgl[i]*betaEoS.T;
   // double q0t = q0 + betaEoS.U2 - betaEoS.U4;
   // double q= pol.GetqFromMu13(E1, q0, 1.0);//minimum q given q0 and E1
   // double qa2t = q0t*q0t - q*q;//max qa2t given q0 and E1
    double ee;
    if (abs(q0) > 30.0*T) break;
    ee = log(pol.CalculateDGamDq0(E1, q0)) + pol.xgl[i];
    integral += pol.wgl[i] * exp(ee);
    std::cout<<" estar: "<<estar<<" q0: "<<q0<<" crx integrand: "<< pol.wgl[i] * exp(ee) <<" crx: "<<integral<<std::endl;
    datanc <<" estar: "<<estar<<" q0: "<<q0<<" crx integrand: "<< pol.wgl[i] * exp(ee) <<" crx: "<<integral<<std::endl;
  }

  for (int i=0; i<pol.NNPGL; ++i) {
    double q0 = estar + pol.xgl[i]*betaEoS.T;
    if (q0>E1 - betaEoS.M3) break;
    double ee = log(pol.CalculateDGamDq0(E1, q0)) + pol.xgl[i];
    integral += pol.wgl[i] * exp(ee);
    std::cout<<" estar: "<<estar<<" q0: "<<q0<<" crx integrand: "<< pol.wgl[i] * exp(ee) <<" crx: "<<integral<<std::endl;
    datanc <<" estar: "<<estar<<" q0: "<<q0<<" crx integrand: "<< pol.wgl[i] * exp(ee) <<" crx: "<<integral<<std::endl;
  }*/

 //***************output dgammadq0dcos for various q0********************************************************************
 /* for (int i=0;i<1;++i) {
  double q0 = estar - pol.xgl[i]*betaEoS.T;
  double p3 = sqrt((E1-q0)*(E1-q0) - st.M3*betaEoS.M3);
  double mu13cross = std::max((E1*E1 + p3*p3 - q0*q0)/(2.0*E1*p3), -1.0);
  double delta = (mu13cross + 1.0) / 2.0;
  double avg = (mu13cross - 1.0) / 2.0;

  double fac = pol.GetCsecPrefactor(E1, q0);
  double integral = 0.0;
  for (int i=0; i<pol.NPGJ; ++i) {
      double mu = pol.xx[i]*delta + avg;
      integral += pol.ww[i] * pol.GetResponse(E1, q0, pol.GetqFromMu13(E1, q0, mu));
      
      std::cout<<" q0: "<<q0<<" cos: "<<mu<<" q: "<< pol.GetqFromMu13(E1, q0, mu)<<" ImPI: "<<pol.GetImPI(q0, pol.GetqFromMu13(E1, q0, mu))<<" crx integrand: "<<pol.ww[i] * pol.GetResponse(E1, q0, pol.GetqFromMu13(E1, q0, mu))*delta <<" fac: "<<fac<<" crx: "<<2.0*Constants::Pi*fac*integral<<std::endl;
      datanc <<" q0: "<<q0<<" cos: "<<mu<<" q: "<< pol.GetqFromMu13(E1, q0, mu)<<" crx integrand: "<<pol.ww[i] * pol.GetResponse(E1, q0, pol.GetqFromMu13(E1, q0, mu))*delta <<" fac: "<<fac<<" crx: "<<2.0*Constants::Pi*fac*integral<<std::endl;

  }
  }*/
  
 
 /* for (int c=0;c<100;c++)
  {
          cos=cos+0.02;
          w=-50.0;
          for (int ww=0;ww<100;ww++)
          {
                  w=w+0.7;
                  qq = pol.GetqFromMu13(E1, w, cos);
		  response=pol.GetResponse( E1, w, qq);
		  prefac=pol.GetCsecPrefactor( E1, w);
		  gamma=pol.gamma0(w, qq);
                  crx=pol.CalculateDGamDq0Dmu13( E1, w, cos);
                  crxtot=crxtot+crx*0.7*0.02;
                  datanc << "cos: " <<cos<<" q: "<<qq<<" w: "<<w<<" crx: "<< crx <<" response: "<<response<<" prefac: "<<prefac<<" gamma0 "<<gamma<< std::endl;
          }
  }
  datanc << " crxtot: " <<crxtot << std::endl;*/
  }

  //****************************Charged current*********************************************************
  double cos=-1.0;
  double w=-50.0;
  double crx=0.0;
  double qq;
  double crxtot=0;
  double response;
  double prefac;
  double gamma;
  
  

  E1 = 12;
  double deltaM = 1.293;
  bool antiNu = false;
  M4 = M2 - deltaM;
  M3 = 0.511;


//***********************************************************************************************************************************

  for (int t=1;t<20;t++)
  {
  
 // n2 = 1*1.e-3*densFac*0.6*t;
 // n4 = 1*1.e-3*densFac*0.4*t;
 n2=n2eos[t]*densFac;
 n4=n4eos[t]*densFac;

 // U2=40.0/(0.16*densFac)*(n2+n4)*(0.5-0.4);
 // U4=-40.0/(0.16*densFac)*(n2+n4)*(0.5-0.4);

 // U2=40.0/(0.16*densFac)*(n2+n4)*(0.5-n4eos[t]/(n2eos[t]+n4eos[t]));
 // U4=-40.0/(0.16*densFac)*(n2+n4)*(0.5-n4eos[t]/(n2eos[t]+n4eos[t]));

  //adding isoscaler part of the potential to see if the final outcome will be changed
  // U2=-50.0/(0.16*densFac)*(n2+n4)+40.0/(0.16*densFac)*(n2+n4)*(0.5-n4eos[t]/(n2eos[t]+n4eos[t]));
  // U4=-50.0/(0.16*densFac)*(n2+n4)-40.0/(0.16*densFac)*(n2+n4)*(0.5-n4eos[t]/(n2eos[t]+n4eos[t]));

 // *****************Plug in EoS U1, U2, effective mass*****************************

 
 FluidState betaEoS;
  
   u2eos[t]=mu2eos[t]*mevFac-(kfn[t]*mevFac)*(kfn[t]*mevFac)/(2.0*msn[t]*mevFac)-M2;//u2eos is with rest mass. now remove the rest mass and Ui=mui-mu_free, here mu_free is the fermi energy
   u4eos[t]=mu4eos[t]*mevFac-(kfp[t]*mevFac)*(kfp[t]*mevFac)/(2.0*msp[t]*mevFac)-M4;

  //  betaEoS=FluidState::StateFromDensities(T, M2, M4, n2, n4, U2, U4, M3, n4);
 // betaEoS=FluidState::StateFromDensities(T, M2, M4, n2, n4, u2eos[t], u4eos[t], M3);
  betaEoS=FluidState::StateFromDensities(T, msn[t]*mevFac, msp[t]*mevFac, n2, n4, u2eos[t], u4eos[t], M3, n4);
 // betaEoS=FluidState::StateFromDensities(T, msn[t]*mevFac, msp[t]*mevFac, n2, n4, U2, U4, M3);
 // betaEoS=FluidState::StateFromDensities(T, msn[t]*mevFac, M4*msn[t]*mevFac/M2, n2, n4, u2eos[t], u4eos[t], M3);
 // betaEoS=FluidState::StateFromDensities(T, msn[t]*mevFac, M4*msn[t]*mevFac/M2, n2, n4, u2eos[t], u4eos[t], M3);// betaEoS=FluidState::BetaEquilibriumConsistentPotential( T,  M2, M4, n2+n4, M3, 40.0/(0.16*densFac));
//
  // betaEoS.Mu2=betaEoS.Mu2-msn[t]*mevFac;
  // betaEoS.Mu4=betaEoS.Mu4-msp[t]*mevFac;// if use nuopac to determine mu, and reddy's non-rel to find response, the chemical potential should be tranformed to non-rel version, by substract ing effective mass from it
  // betaEoS.Mu3=betaEoS.Mu3-0.511;

   betaEoS.Mu2=mu2eos[t]*mevFac-939.0;
   betaEoS.Mu4=mu4eos[t]*mevFac-937.7;// if use mu from nrapr, and reddy's non-rel to find response, the chemical potential should be tranformed to non-rel version, by substracting bare  mass from it
   betaEoS.Mu3=mu3eos[t]*mevFac-0.511;

  // betaEoS.Mu2=betaEoS.Mu2-939.0;
  // betaEoS.Mu4=betaEoS.Mu4-937.7;// if use mu from nuopac with bare mass, and reddy's non-rel to find response, the chemical potential should be tranformed to non-rel version, by substracting bare  mass from it
  // betaEoS.Mu3=betaEoS.Mu3-0.511;
  //
   //*****Virial**************************
  /* betaEoS=FluidState::StateFromDensities(T, M2, M4, n2eos[t]*densFac, n4eos[t]*densFac, u2eos[t], u4eos[t], M3);
   betaEoS.Mu2=mu2eos[t];
   betaEoS.Mu4=mu4eos[t];
   betaEoS.Mu3=mu2eos[t]-mu4eos[t];*/
  //*****************************************

  //******test residual interaction*****************
      double t0,t1,t2,t3,x0,x1,x2,x3,epsilon;
     // ****set skyrme parameters********************
     //ft2
     /*  epsilon= 2.1493862350*1.0E-01;
       t0 = -2.3279811323*1.0E+03;
       t1=2.9748565632*1.0E+02;
       t2= 1.8783378926*1.0E+03;
       t3= 1.4661502679*1.0E+04;
       x0= 1.6575377249*1.0E-01;
       x1= -5.3459402501*1.0E+00;
       x2= -1.3322633637*1.0E+00;
       x3= 9.9611039279*1.0E-02;*/

     //NRAPR
    /*   epsilon= 1.4416*1.0E-01;
       t0 = -2.7197*1.0E+03;
       t1=4.1764*1.0E+02;
       t2= -6.6687*1.0E+01;
       t3= 1.5042*1.0E+04;
       x0= 1.6154*1.0E-01;
       x1= -4.7986*1.0E-02;
       x2= 2.717*1.0E-02;
       x3= 1.3611*1.0E-01;*/

     // *************************
  /*   double rou=(betaEoS.n2+betaEoS.n4)/densFac;
     double kf=pow(0.5*rou*3*3.14159*3.14159,1.0/3.0);
    // double kfn=pow(rou*3*3.14159*3.14159,1.0/3.0);
     double f0=3.0/4.0*t0+3.0/8.0*t1*kf*kf+5.0/8.0*t2*kf*kf+1.0/2.0*t2*x2*kf*kf+(epsilon+1.0)*(epsilon+2.0)/16.0*t3*pow(rou,epsilon);
     double f0p=-1.0/4.0*t0-1.0/2.0*t0*x0-1.0/8.0*t1*kf*kf-1.0/4.0*t1*x1*kf*kf+1.0/8.0*t2*kf*kf+1.0/4.0*t2*x2*kf*kf-1.0/24.0*t3*pow(rou,epsilon)-1.0/12.0*t3*x3*pow(rou,epsilon);
     double g0=-1.0/4.0*t0+1.0/2.0*x0*t0-1.0/8.0*t1*kf*kf+1.0/4.0*t1*x1*kf*kf+1.0/8.0*t2*kf*kf+1.0/4.0*t2*x2*kf*kf-1.0/24.0*t3*pow(rou,epsilon)+1.0/12.0*t3*x3*pow(rou,epsilon);
     double g0p=-1.0/4.0*t0-1.0/8.0*t1*kf*kf+1.0/8.0*t2*kf*kf-1.0/24.0*t3*pow(rou,epsilon);
    // double f0n=1.0/2.0*t0-1.0/2.0*t0*x0+1.0/8.0*t1*kfn*kfn-1.0/4.0*t1*x1*kfn*kfn+3.0/4.0*t2*kfn*kfn+3.0/4.0*t2*x2*kfn*kfn+(epsilon+1.0)*(epsilon+2.0)/24.0*t3*pow(rou,epsilon)-(epsilon+1.0)*(epsilon+2.0)/24.0*t3*x3*pow(rou,epsilon);
    // double g0n=-1.0/2.0*t0+1.0/2.0*t0*x0-1.0/4.0*t1*kfn*kfn+1.0/4.0*t1*x1*kfn*kfn+1.0/4.0*t2*kfn*kfn+1.0/4.0*t2*x2*kfn*kfn-1.0/12.0*t3*pow(rou,epsilon)+1.0/12.0*t3*x3*pow(rou,epsilon);
  
      double vf=2*f0p;
      double vgt=2*g0p;
      double vrpa=vgt/pow(197.3,3);*/

    double fnn;
    double fpp;
    double gnn;
    double gpp;
    fnn=fpp=-0.0000844949;
    gnn=gpp=0.0000759672;
  double f0p=0.000118591;
  double g0p=0.0000730909;
  double vf=2.0*f0p;
  double vgt=2.0*g0p;
  double vrpa=vgt;//already in unit of MeV-2


  //************************
std::cout <<" nrparEos: "<<" baryon den: "<<(betaEoS.n2+betaEoS.n4)/densFac<<" y_N: "<<betaEoS.n2/(betaEoS.n2+betaEoS.n4)<< " n2: "<<betaEoS.n2/densFac<<" n4: "<<betaEoS.n4/densFac<<" M3: "<<betaEoS.M3<<" m2eff: "<<betaEoS.M2<<" m4eff: "<<betaEoS.M4<<" U2: "<<betaEoS.U2<<" U4: "<<betaEoS.U4<<" deltaU: "<<betaEoS.U2-betaEoS.U4<<" mu2: "<<betaEoS.Mu2<<" Mu4: "<<betaEoS.Mu4<<" mu2-mu4 "<<betaEoS.Mu2-betaEoS.Mu4<<" Mu3: "<<betaEoS.Mu3<<" n3: "<<betaEoS.n3/densFac<<" vf: "<<vf<<" vgt: "<<vgt<< std::endl;
datacc <<" nrparEos: "<<" baryon den: "<<(betaEoS.n2+betaEoS.n4)/densFac<<" y_N: "<<betaEoS.n2/(betaEoS.n2+betaEoS.n4)<< " n2: "<<betaEoS.n2/densFac<<" n4: "<<betaEoS.n4/densFac<<" M3: "<<betaEoS.M3<<" m2eff: "<<betaEoS.M2<<" m4eff: "<<betaEoS.M4<<" U2: "<<betaEoS.U2<<" U4: "<<betaEoS.U4<<" deltaU: "<<betaEoS.U2-betaEoS.U4<<" mu2: "<<betaEoS.Mu2<<" Mu4: "<<betaEoS.Mu4<<" mu2-mu4 "<<betaEoS.Mu2-betaEoS.Mu4<<" Mu3: "<<betaEoS.Mu3<<" n3: "<<betaEoS.n3/densFac<<" vf symmetric: "<<vf<<" vgt symmetric: "<<vgt<< std::endl;
 
 
 
 
 
 // ********************************************************************************
 
 // beta=FluidState::BetaEquilibriumConsistentPotential( T,  M2,  M4, 1.e-3*densFac*t,  M3, 40.0/(0.16*densFac));
 // U2=40.0/(0.16*densFac)*(n2+n4)*(0.5-0.4)-50;
 // U4=-40.0/(0.16*densFac)*(n2+n4)*(0.5-0.4)-50;
 // U2=0.0;
 // U4=0.0;
 // st = FluidState::StateFromDensities(T, M2, M4, n2, n4, U2, U4, M3);

  WeakCouplings ncap = WeakCouplings::NuCapture();
  ncap.F2 = 0.0;

 // PolarizationNonRel pol = PolarizationNonRel(st, ncap, antiNu, false, true);
 // PolarizationNonRel pol = PolarizationNonRel(beta, ncap, antiNu, false, true);
  PolarizationNonRel pol = PolarizationNonRel(betaEoS, ncap, antiNu, false, true);
// PolarizationNonRel pol = PolarizationNonRel(betaEoS, ncap, antiNu, false, false);//turn off pauli blocking
  double noWm = pol.CalculateInverseMFP(E1)/Constants::HBCFmMeV*1.e13;

  ncap = WeakCouplings::NuCapture();
 // pol = PolarizationNonRel(st, ncap, antiNu, false, true);
 // pol = PolarizationNonRel(beta, ncap, antiNu, false, true);
  pol = PolarizationNonRel(betaEoS, ncap, antiNu, false, true);
 
  double Wm = pol.CalculateInverseMFP(E1)/Constants::HBCFmMeV*1.e13;
 // double Ee=E1+deltaM+beta.U2-beta.U4;
  double Ee=E1+deltaM+betaEoS.U2-betaEoS.U4;
 // double Ee=E1+deltaM+st.U2-st.U4;
 // double pe=sqrt(Ee*Ee-M3*M3)*(1-1/(exp((Ee-beta.Mu3)/T)+1));
  double pe=sqrt(Ee*Ee-M3*M3)*(1-1/(exp((Ee-betaEoS.Mu3)/T)+1));
 // double pe=sqrt(Ee*Ee-M3*M3)*(1-1/(exp((Ee-st.Mu3)/T)+1));
  // Elastic approximation to the cross section
 // sig0 = 4.0*pow(Constants::GfMeV*(E1+deltaM)*Constants::HBCFmMeV, 2)
   //   / Constants::Pi;
  sig0 = 4.0*pow(Constants::GfMeV*(0.511)*Constants::HBCFmMeV, 2)
      / Constants::Pi;
  // We divide by 16.0 instead of 4.0 because of our convention for the coupling
  // constants
  sig0 *= (0.0*pow(ncap.Cv, 2) + 3.0*pow(ncap.Ca, 2))/4.0;
 // sig0 *= (1.0*pow(ncap.Cv, 2) + 3.0*pow(ncap.Ca, 2))/4.0;
 // sig0 *= sqrt(1.0 - pow(0.511/(E1+deltaM),2));
 // elastic1 = sig0*beta.n2/densFac/0.511/0.511*Ee*Ee*1.e13;
  elastic1 = sig0*betaEoS.n2/densFac/0.511/0.511*Ee*Ee*1.e13;
 // elastic1 = sig0*st.n2/densFac/0.511/0.511*Ee*Ee*1.e13; 
 // sig0 *= sqrt(1.0 - pow(0.511/(E1+deltaM+beta.U2-beta.U4),2));
  sig0 *= sqrt(1.0 - pow(0.511/(E1+deltaM+betaEoS.U2-betaEoS.U4),2));
 // sig0 *= sqrt(1.0 - pow(0.511/(E1+deltaM+st.U2-st.U4),2));
  
 
 // elastic = sig0*beta.effectiveDensity/densFac*pe*Ee/0.511/0.511*1.e13;
  elastic = sig0*betaEoS.effectiveDensity/densFac*pe*Ee/0.511/0.511*1.e13;
 // elastic = sig0*st.effectiveDensity/densFac*pe*Ee/0.511/0.511*1.e13;
 // WeakCouplings coup = WeakCouplings::NuCapture();
// std::cout <<"charged current: "<<(beta.n2+beta.n4)/densFac<<" "<<beta.n2/(beta.n2+beta.n4)<<" noWm: "<< noWm << " Wm: " << Wm << " elastic: "<<elastic<<" elastic1: "<<elastic1<<" noWm/elastc: " << noWm/elastic << " noWm/elastic1: "<< noWm/elastic1<<" M3: "<<beta.M3<<" neff/n: "<<st.effectiveDensity/beta.n2<<" pe/Ee: "<<pe/Ee<<" deltaU: "<<beta.U2-beta.U4<<" mu3 "<<beta.Mu3<< std::endl;
  
//std::cout <<"charged current: "<<(st.n2+st.n4)/densFac<<" "<<st.n2/(st.n2+st.n4)<<" noWm: "<< noWm << " Wm: " << Wm << " elastic: "<<elastic<<" elastic1: "<<elastic1<<" noWm/elastc: " << noWm/elastic << " noWm/elastic1: "<< noWm/elastic1<<" M3: "<<st.M3<<" neff/n: "<<st.effectiveDensity/st.n2<<" pe/Ee: "<<pe/Ee<<" deltaU: "<<st.U2-st.U4<<" mu2-mu4 "<<st.Mu2-st.Mu4<< std::endl;

//datacc <<"charged current: "<<(st.n2+st.n4)/densFac<<" "<<st.n2/(st.n2+st.n4)<<" noWm: "<< noWm << " Wm: " << Wm << " elastic: "<<elastic<<" elastic1: "<<elastic1<<" noWm/elastc: " << noWm/elastic << " noWm/elastic1: "<< noWm/elastic1<<" M3: "<<st.M3<<" neff/n: "<<st.effectiveDensity/st.n2<<" pe/Ee: "<<pe/Ee<<" deltaU: "<<st.U2-st.U4<<" mu3 "<<st.Mu3<< std::endl;
//std::cout <<"charged current: "<<(beta.n2+beta.n4)/densFac<<" "<<beta.n2/(beta.n2+beta.n4)<<" noWm: "<< noWm << " Wm: " << Wm << " elastic: "<<elastic<<" elastic1: "<<elastic1<<" noWm/elastc: " << noWm/elastic << " noWm/elastic1: "<< noWm/elastic1<<" M3: "<<beta.M3<<" neff/n: "<<beta.effectiveDensity/beta.n2<<" pe/Ee: "<<pe/Ee<<" deltaU: "<<beta.U2-beta.U4<<" mu2-mu4 "<<beta.Mu2-beta.Mu4<< std::endl;
//datacc <<"charged current: "<<(beta.n2+beta.n4)/densFac<<" "<<beta.n2/(beta.n2+beta.n4)<<" noWm: "<< noWm << " Wm: " << Wm << " elastic: "<<elastic<<" elastic1: "<<elastic1<<" noWm/elastc: " << noWm/elastic << " noWm/elastic1: "<< noWm/elastic1<<" M3: "<<beta.M3<<" neff/n: "<<beta.effectiveDensity/beta.n2<<" pe/Ee: "<<pe/Ee<<" deltaU: "<<beta.U2-beta.U4<<" mu2-mu4 "<<beta.Mu2-beta.Mu4<< std::endl;
std::cout <<"charged current: "<<(betaEoS.n2+betaEoS.n4)/densFac<<" "<<betaEoS.n2/(betaEoS.n2+betaEoS.n4)<<" noWm: "<< noWm << " Wm: " << Wm << " elastic: "<<elastic<<" elastic1: "<<elastic1<<" noWm/elastc: " << noWm/elastic << " noWm/elastic1: "<< noWm/elastic1<<" M3: "<<betaEoS.M3<<" neff/n: "<<betaEoS.effectiveDensity/betaEoS.n2<<" pe/Ee: "<<pe/Ee<<" deltaU: "<<betaEoS.U2-betaEoS.U4<<" mu2-mu4 "<<betaEoS.Mu2-betaEoS.Mu4<<" ImPi0: "<<pol.GetImPI(-20.9,30)<< std::endl;
datacc <<"charged current: "<<(betaEoS.n2+betaEoS.n4)/densFac<<" "<<betaEoS.n2/(betaEoS.n2+betaEoS.n4)<<" noWm: "<< noWm << " Wm: " << Wm << " elastic: "<<elastic<<" elastic1: "<<elastic1<<" noWm/elastc: " << noWm/elastic << " noWm/elastic1: "<< noWm/elastic1<<" M3: "<<betaEoS.M3<<" neff/n: "<<betaEoS.effectiveDensity/betaEoS.n2<<" pe/Ee: "<<pe/Ee<<" deltaU: "<<betaEoS.U2-betaEoS.U4<<" U2: "<<betaEoS.U2<<" U4: "<<betaEoS.U4<<" mu2-mu4 "<<betaEoS.Mu2-betaEoS.Mu4<<" Mu2: "<<betaEoS.Mu2<<" Mu4: "<<betaEoS.Mu4<< std::endl;
//**********Set integration range***********
  double vel=sqrt((3.0*T)/M2);
  double wmin;
  double wmax;
 // wmin=-3.0*vel*3*T-0.00000789;
  wmin=-100.0;
 // wmax=3.0*(vel*3*T+3*T*3*T/(2.0*M2));
  wmax=50.0;
  double dw=(wmax-wmin)/100;
//******************************************

 for (int k=1;k<100;k++)
  {

    w=wmin+dw*k;
//   **********new version**************
    auto pt= pol.CalculateBasePolarizations(w, 10*T);//charged current
   
    piL=pt[1];
   
//  ***********PiLRe**********************
    piLRe=pol.GetImPI2(w,10*T);
 

 // piL=pol.GetImPI(w,5*T);//piL is 2*ImPI0
 // piLRe=pol.GetImPI2(w, 5*T);//piLRe is 2*RePI0
  //vf and vgt for cc comes from arXiv:1205.4066, in unit of MeV^-2, note that in the arxiv paper the unit is wrong, should be fm^2 rather than fm^-2
 /* double vf=5.1*1.0E-5;
  double vgt=2.8*1.0E-5;
  double vrpa=vf;*/

 // double t0,t1,t2,t3,x0,x1,x2,x3,epsilon;
 /* //ft7
       epsilon= 1.6112652355*1.0E-01;
       t0=-2.8515584390*1.0E+03;
       t1=2.9926667400*1.0E+02;
       t2=2.1204461670*1.0E+03;
       t3=1.7480390049*1.0E+04;
       x0=9.5708919697*1.0E-02;
       x1=-5.9493754771*1.0E+00;
       x2=-1.3211268948 *1.0E+00;
       x3=3.7117434733*1.0E-02;*/
  //NRAPR
       epsilon= 1.4416*1.0E-01;
       t0 = -2.7197*1.0E+03;
       t1=4.1764*1.0E+02;
       t2= -6.6687*1.0E+01;
       t3= 1.5042*1.0E+04;
       x0= 1.6154*1.0E-01;
       x1= -4.7986*1.0E-02;
       x2= 2.717*1.0E-02;
       x3= 1.3611*1.0E-01;


     double rou=(betaEoS.n2+betaEoS.n4)/densFac;
     double roun=betaEoS.n2/densFac;
     double roup=betaEoS.n4/densFac;
    
     double kfnTest=pow(roun*3*3.14159*3.14159,1.0/3.0);
     double kfpTest=pow(roup*3*3.14159*3.14159,1.0/3.0);
     double kf=pow(0.5*rou*3*3.14159*3.14159,1.0/3.0);

    // double f0p=-1.0/4.0*t0-1.0/2.0*t0*x0-1.0/8.0*t1*kf*kf-1.0/4.0*t1*x1*kf*kf+1.0/8.0*t2*kf*kf+1.0/4.0*t2*x2*kf*kf-1.0/24.0*t3*pow(rou,epsilon)-1.0/12.0*t3*x3*pow(rou,epsilon);
    // double g0p=-1.0/4.0*t0-1.0/8.0*t1*kf*kf+1.0/8.0*t2*kf*kf-1.0/24.0*t3*pow(rou,epsilon);
    // double vf=2.0*f0p;// for symmetric nuclear matter
    // double vgt=2.0*g0p;// for symmetric nuclear matter
    // double fnn=0.5*(t0*(1.0-x0)+1.0/6.0*t3*pow(rou,epsilon)*(1.0-x3)+2.0/3.0*epsilon*t3*pow(rou,epsilon-1)*((1+x3/2.0)*rou-(1.0/2.0+x3)*roun)+1.0/6.0*epsilon*(epsilon-1.0)*t3*pow(rou,epsilon-2.0)*((1+x3/2.0)*pow(rou,2.0)-(0.5+x3)*(roun*roun+roup*roup)))+0.25*(t1*(1-x1)+3*t2*(1+x2))*kfn*kfn;
    // double fnp=0.5*(t0*(2.0+x0)+1.0/6.0*t3*pow(rou,epsilon)*(2.0+x3)+1.0/2.0*epsilon*t3*pow(rou,epsilon)+1.0/6.0*epsilon*(epsilon-1.0)*t3*pow(rou,epsilon-2.0)*((1+x3/2.0)*pow(rou,2.0)-(0.5+x3)*(roun*roun+roup*roup)))+0.5*0.25*(t1*(2.0+x1)+t2*(2.0+x2))*(kfn*kfn+kfp*kfp);
    // double gnn=0.5*(t0*(x0-1)+1.0/6.0*t3*pow(rou,epsilon)*(x3-1.0))+0.25*(t1*(x1-1)+t2*(1+x2))*kfn*kfn;
    // double gnp=0.5*(t0*x0+1.0/6.0*t3*pow(rou,epsilon)*x3)+0.5*0.25*(t1*x1+t2*x2)*(kfn*kfn+kfp*kfp);
   //  double vf=fnn-fnp;
   //  double vgt=gnn-gnp;
       double w1nnVec,w1npVec,w1nnAx,w1npAx,w2nnVec,w2npVec,w2nnAx,w2npAx;
      w1nnVec=t0*(1.0-x0)+1.0/6.0*t3*pow(rou,epsilon)*(1.0-x3)+2.0/3.0*epsilon*t3*pow(rou,epsilon-1)*((1+x3/2.0)*rou-(1.0/2.0+x3)*roun)+1.0/6.0*epsilon*(epsilon-1.0)*t3*pow(rou,epsilon-2.0)*((1+x3/2.0)*pow(rou,2.0)-(0.5+x3)*(roun*roun+roup*roup));
      w1npVec=t0*(2.0+x0)+1.0/6.0*t3*pow(rou,epsilon)*(2.0+x3)+1.0/2.0*epsilon*t3*pow(rou,epsilon)+1.0/6.0*epsilon*(epsilon-1.0)*t3*pow(rou,epsilon-2.0)*((1+x3/2.0)*pow(rou,2.0)-(0.5+x3)*(roun*roun+roup*roup));
      w1nnAx=t0*(x0-1)+1.0/6.0*t3*pow(rou,epsilon)*(x3-1.0);
      w1npAx=t0*x0+1.0/6.0*t3*pow(rou,epsilon)*x3;
      w2nnVec=0.25*(t1*(1-x1)+3*t2*(1+x2));
      w2npVec=0.25*(t1*(2.0+x1)+t2*(2.0+x2));
      w2nnAx=0.25*(t1*(x1-1)+t2*(1+x2));
      w2npAx=0.25*(t1*x1+t2*x2);

       vf=0.5*(w1nnVec-w1npVec)+(w2nnVec-w2npVec)*kfnTest*kfnTest;//kf should be the hole momenta at fermi see surface, here the transition is (pn^-1,pn^-1), the hole is neutron hole
       vgt=0.5*(w1nnAx-w1npAx)+(w2nnAx-w2npAx)*kfnTest*kfnTest;
       vrpa=vgt/pow(197.3,3);//change unit to MeV-2
      double  vrpaVec=vf/pow(197.3,3);//change unit to MeV-2 

     //virial
  /*   double fnn;
    double fpp;
    double gnn;
    double gpp;
    fnn=fpp=-0.0000844949;
    gnn=gpp=0.0000759672;
  double f0p=0.000118591;
  double g0p=0.0000730909;
  double vf=2.0*f0p;
  double vgt=2.0*g0p;
  double vrpa=vgt;//already in unit of MeV-2
  double vrpaVec=vf; */
  
  double piLRPA;
 
  double FermiF;
  double zz;
  zz=(w+betaEoS.Mu2-betaEoS.Mu4)/T;
 // zz=(w+beta.Mu2-beta.Mu4)/T;
  FermiF=1/(1-exp(-zz));
  response=pol.GetResponse( E1, w, 10*T);
  piLRPA=2*(piL/2)/((1-vrpa*(piLRe/2))*(1-vrpa*(piLRe/2))+vrpa*vrpa*(piL/2)*(piL/2))*FermiF;//piLRPA is s(q0,q)
  double piLRPAVec=2*(piL/2)/((1-vrpaVec*(piLRe/2))*(1-vrpaVec*(piLRe/2))+vrpaVec*vrpaVec*(piL/2)*(piL/2))*FermiF;
 // piLRPA=2*(piL/2)*FermiF;//piLRPA is s(q0,q)
 // std::cout << "charged current density: " <<n2+n4<<" full "<< Wm <<" w "<<w<< " piL(w) " <<piL<<" piLRe(w) "<<piLRe<<" piLRPA "<<piLRPA<<" piLRPA vector "<<piLRPAVec<<" response "<<response<<" FermiF: "<<FermiF<< std::endl;
  datacc << "charged current density: " <<betaEoS.n2<<" vgt: "<<vgt <<" vf: "<<vf<< " w "<<w<< " piL(w) " <<piL<<" piLRe(w) "<<piLRe<<" piLRPA "<<piLRPA<<" piLRPA vector "<<piLRPAVec<<" FermiF "<<FermiF<< std::endl;
 
  }

//print S(q0,q)
  /* cos=-1.0;
   w=-150.0;
   double piL;
   double piLRe;
   double piLRPA;
   double FermiF;
   double vf=1.9635*1.0E-6*log((n2+n4)/densFac)*log((n2+n4)/densFac)+0.00003499*log((n2+n4)/densFac)+0.00007305;
   double vgt=-2.6304*1.0E-6*log((n2+n4)/densFac)*log((n2+n4)/densFac)-0.00003932*log((n2+n4)/densFac)-0.00005173;
   double vrpa=vf;*/
 /* for (int c=0;c<100;c++)
  {
          cos=cos+0.02;
          w=-150.0;
          for (int ww=0;ww<200;ww++)
          {
                  w=w+1.0;
                  qq = pol.GetqFromMu13(E1, w, cos);
                  piL=pol.GetImPI(w,qq);//piL is 2*ImPI0
                 // piLRe=pol.GetImPI2(w, qq);//piLRe is 2*RePI0
		  double zz;
                  zz=(w+betaEoS.Mu2-betaEoS.Mu4)/T;
                  FermiF=1/(1-exp(-zz));
                 // piLRPA=2*(piL/2)/((1-vrpa*(piLRe/2))*(1-vrpa*(piLRe/2))+vrpa*vrpa*(piL/2)*(piL/2))*FermiF;//piLRPA is s(q0,q)
                  piLRPA=2*(piL/2)*FermiF;//piLRPA is s(q0,q)
                   
		   double em;
		   double q0t = w + betaEoS.U2 - betaEoS.U4;//orig one
                   // double q0t = q0 + u2pol - u4pol;
                   double qa2t = q0t*q0t - qq*qq;
                   // I don't completely understand this condition, but it seems to be necessary
                   // to suppress noise at larger q_0
                   if (qa2t > pow(betaEoS.M2 - betaEoS.M4, 2)*0.0) em=0;//orig one
                   // if (qa2t > pow(m2pol - m4pol, 2)*0.0) return {0.0, 0.0, 0.0, 0.0};
                   if (qa2t < 1.e-1 && qa2t >= 0.0) em=0.0;
                   double beta = 1.0 + ( betaEoS.M2* betaEoS.M2 -  betaEoS.M4* betaEoS.M4)/qa2t;//orig
 // double beta = 1.0 + (m2pol*m2pol - m4pol*m4pol)/qa2t;
                   double arg = beta*beta - 4.0* betaEoS.M2* betaEoS.M2/qa2t;
                   if (arg<0.0) em=0.0;
                   em = std::max(-0.5*beta*q0t + 0.5*qq*sqrt(arg), betaEoS.M2);

                  datacc << "cos: " <<cos<<" q: "<<qq<<" w: "<<w<<" piL: "<< piL <<" piLRe: "<<piLRe<<" piLRPA: "<<piLRPA<<" em: "<<em<< std::endl;
          }
  }*/
 // datacc << " crxtot: " <<crxtot << std::endl;
 // double resp=0;


//print response
  
  /* cos=-1.0;
   w=-300.0;
   crx=0.0;
   qq;
   crxtot=0;
  for (int c=0;c<100;c++)
  {
          cos=cos+0.02;
	  w=-300.0;
	  for (int ww=0;ww<400;ww++)
	  {
		  w=w+1.0;
		  qq = pol.GetqFromMu13(E1, w, cos);
		  response=pol.GetResponse( E1, w, qq);
		  prefac=pol.GetCsecPrefactor( E1, w);
                  crx=pol.CalculateDGamDq0Dmu13( E1, w, cos);
		  crxtot=crxtot+crx*1.0*0.02;
		  datacc << "cos: " <<cos<<" q: "<<qq<<" w: "<<w<<" crx: "<< crx <<" response: "<<response<<" prefac: "<<prefac<< std::endl;
	  }
  }
  datacc << " crxtot: " <<crxtot << std::endl;
  double resp=0;*/
  // output dgammadq0*************************
 // double estar = st.M4 + st.U4 - st.M2 - st.U2;
 // estar = std::min(estar, E1 - st.M3);
  double estar = betaEoS.M4 + betaEoS.U4 - betaEoS.M2 - betaEoS.U2;
  estar = std::min(estar, E1 - betaEoS.M3);
 /* double integral = 0.0;
  for (int i=0; i<pol.NNPGL; ++i) {
    double q0 = estar - pol.xgl[i]*betaEoS.T;
    if (abs(q0)>30*betaEoS.T) break;
    double ee = log(pol.CalculateDGamDq0(E1, q0)) + pol.xgl[i];
    integral += pol.wgl[i] * exp(ee);
    std::cout<<" estar: "<<estar<<" q0: "<<q0<<" crx integrand: "<< pol.wgl[i] * exp(ee) <<" crx: "<<integral<<std::endl;
    datacc <<" estar: "<<estar<<" q0: "<<q0<<" crx integrand: "<< pol.wgl[i] * exp(ee) <<" crx: "<<integral<<std::endl;
   // std::cout <<" estar: "<<estar<<" q0: "<<q0<<" crx integrand: "<< pol.wgl[i] * exp(ee) <<" crx: "<<integral<<std::endl;
  }

  for (int i=0; i<pol.NNPGL; ++i) {
    double q0 = estar + pol.xgl[i]*betaEoS.T;
    if (q0>E1 - betaEoS.M3) break;
    double ee = log(pol.CalculateDGamDq0(E1, q0)) + pol.xgl[i];
    integral += pol.wgl[i] * exp(ee);
    datacc <<" estar: "<<estar<<" q0: "<<q0<<" crx integrand: "<< pol.wgl[i] * exp(ee) <<" crx: "<<integral<<std::endl;
    std::cout <<" estar: "<<estar<<" q0: "<<q0<<" crx integrand: "<< pol.wgl[i] * exp(ee) <<" crx: "<<integral<<std::endl;
  }*/
//***************output dgammadq0dcos for various q0********************************************************************
 /* for (int i=29;i<30;++i) {
  double q0 = estar + pol.xgl[i]*betaEoS.T;
  double p3 = sqrt((E1-q0)*(E1-q0) - betaEoS.M3*betaEoS.M3);
  double mu13cross = std::max((E1*E1 + p3*p3 - q0*q0)/(2.0*E1*p3), -1.0);
  double delta = (mu13cross + 1.0) / 2.0;
  double avg = (mu13cross - 1.0) / 2.0;
  
  double fac = pol.GetCsecPrefactor(E1, q0);
  double integral = 0.0;
  for (int i=0; i<pol.NPGJ; ++i) {
      double mu = pol.xx[i]*delta + avg;
      integral += pol.ww[i] * pol.GetResponse(E1, q0, pol.GetqFromMu13(E1, q0, mu));
      std::cout <<" q0: "<<q0<<" cos: "<<mu<<" crx integrand: "<<pol.ww[i] * pol.GetResponse(E1, q0, pol.GetqFromMu13(E1, q0, mu))*delta <<" fac: "<<fac<<" crx: "<<2.0*Constants::Pi*fac*integral<<std::endl;
      datacc <<" q0: "<<q0<<" cos: "<<mu<<" crx integrand: "<<pol.ww[i] * pol.GetResponse(E1, q0, pol.GetqFromMu13(E1, q0, mu))*delta <<" fac: "<<fac<<" crx: "<<2.0*Constants::Pi*fac*integral<<std::endl;

  }
  }*/
 /* for (int i=0;i<16;++i) {
  double q0 = estar + pol.xgl[i]*st.T;
  double p3 = sqrt((E1-q0)*(E1-q0) - st.M3*st.M3);
  double mu13cross = std::max((E1*E1 + p3*p3 - q0*q0)/(2.0*E1*p3), -1.0);
  double delta = (mu13cross + 1.0) / 2.0;
  double avg = (mu13cross - 1.0) / 2.0;

  double fac = pol.GetCsecPrefactor(E1, q0);
  double integral = 0.0;
  for (int i=0; i<pol.NPGJ; ++i) {
      double mu = pol.xx[i]*delta + avg;
      integral += pol.ww[i] * pol.GetResponse(E1, q0, pol.GetqFromMu13(E1, q0, mu));

   //  double vf=5.1*1.0E-5;
   //  double vgt=2.8*1.0E-5;
    double vf=1.9635*1.0E-6*log((n2+n4)/densFac)*log((n2+n4)/densFac)+0.00003499*log((n2+n4)/densFac)+0.00007305;
    double vgt=-2.6304*1.0E-6*log((n2+n4)/densFac)*log((n2+n4)/densFac)-0.00003932*log((n2+n4)/densFac)-0.00005173;
     double vrpa=vgt;

     double piL=pol.GetImPI(q0,pol.GetqFromMu13(E1, q0, mu));//piL is 2*ImPI0
     double piLRe=pol.GetImPI2(q0, pol.GetqFromMu13(E1, q0, mu));//piLRe is 2*RePI0
     double piLmod=2*(piL/2)/((1-vrpa*(piLRe/2))*(1-vrpa*(piLRe/2))+vrpa*vrpa*(piL/2)*(piL/2));

      datacc <<" q: "<<pol.GetqFromMu13(E1, q0, mu)<<" q0: "<<q0<<" cos: "<<mu<<" Response: "<< pol.GetResponse(E1, q0, pol.GetqFromMu13(E1, q0, mu)) <<" piL: "<<piL<<" piLRe: "<<piLRe<<" PiLmod: "<<piLmod<<std::endl;

  }

  
  }*/


 /* for (int k=1;k<21;k++)
  {
 //charged current  
  pol = PolarizationNonRel(st, ncap, antiNu, false, true);
  double q=10;  
  double q0=-10+k*1.0;
  resp=pol.GetResponse(E1, q0, q);
  M3=0.511;
  double E3 =E1-q0;
  double p3 = sqrt((E1-q0)*(E1-q0)-st.M3*st.M3);
  double theta13=(q*q-E1*E1-p3*p3)/(-2*E1*p3);
  double Term= 8.0*E1*E3*(ncap.Cv*ncap.Cv*(1+theta13)+ncap.Ca*ncap.Ca*(3-theta13)); 
  std::cout << "charged current: " <<1.e-3*densFac*t<<" "<<st.n2/(1.e-4*densFac*t) << " LI: " << resp <<" piL=LI/Term: "<<resp/Term<<" theta13: "<<theta13<< std::endl;
  //neutral current
  M3=0.0;
  st = FluidState::StateFromDensities(T, M2, M4, n2, n4, 20.0*n2/(0.16*densFac), -20.0*n2/(0.16*densFac), M3, n2);
  pol= PolarizationNonRel(st, nscat, false);  
  p3 = sqrt((E1-q0)*(E1-q0)-st.M3*st.M3);
  theta13=sqrt((q*q-E1*E1-p3*p3)/(-2*E1*p3));
  Term= 8.0*E1*E3*(nscat.Cv*nscat.Cv*(1+theta13)+nscat.Ca*nscat.Ca*(3-theta13));
  std::cout << "neutral current: " <<1.e-3*densFac*t<<" "<<st.n2/(1.e-4*densFac*t) << " LI: " << resp <<" piL=LI/Term: "<<resp/Term<<" theta13: "<<theta13<< std::endl;

  }*/
  }
 
  // Test that the simple EoS stuff is working correctly
/* for (int i=1; i<100;i++)
 {  
  double n2 = 1.e-4*densFac*i;
  double n4 = 1.e-4*densFac*i;
  FluidState st = FluidState::StateFromDensities(T, M2, M4, n2, n4, U2, U4, M3, n2);
  FluidState stInverse(T, M2, M4, st.Mu2, st.Mu4, U2, U4, M3, st.Mu3);   
  if (fabs(1.0 - stInverse.n2/n2)>1.e-3) { 
    std::cerr << "Bad density find in n2: " << n2 
    << " " << stInverse.n2 << std::endl;
    return 1; 
  } 
  if (fabs(1.0 - stInverse.n4/n4)>1.e-3) {
    std::cerr << "Bad density find in n4: " << n4
    << " " << stInverse.n4 << std::endl; 
    return 1;
  }   
 
  // Now test a neutral current cross section in the low density, low energy limit//

 
   
  WeakCouplings nscat = WeakCouplings::NeutronScattering(); 
  nscat.F2 = 0.0; 
  Polarization pol(st, nscat, false);
  
  double E1 = 10.0; 
  double full = pol.CalculateInverseMFP(E1)/Constants::HBCFmMeV*1.e13; 
  
  sig0 = 
      4.0*pow(Constants::GfMeV*E1*Constants::HBCFmMeV, 2)/Constants::Pi; 
  sig0 = sig0*(pow(nscat.Cv, 2) + 3.0*pow(nscat.Ca, 2))/4.0;
  elastic = sig0*st.effectiveDensity/densFac*1.e13; 
//  std::cout << "Neutral current: " << n2 <<" "<< full << " " << elastic << " " 
//      << 1.0-full/elastic << std::endl;
//  datanc << "Neutral current: " << n2 <<" "<< full << " " << elastic << " "
//      << 1.0-full/elastic << std::endl;
  if (fabs((full-elastic)/full)>1.e-2) return 1;
  //datanc.close();
  // Now test a charged current cross section


  E1 = 1.0;
  deltaM = 1.293; 
  antiNu = false; 
  M4 = M2 - deltaM; 
  M3 = 0.511;
  n2 = 1.e-4*densFac*0.9*i;
  n4 = 1.e-4*densFac*0.1*i;
  st = FluidState::StateFromDensities(T, M2, M4, n2, n4, U2, U4, M3); 
  
  WeakCouplings ncap = WeakCouplings::NuCapture(); 
  ncap.F2 = 0.0;
  pol = Polarization(st, ncap, antiNu, false); 
  double noWm = pol.CalculateInverseMFP(E1)/Constants::HBCFmMeV*1.e13; 
  
  ncap = WeakCouplings::NuCapture(); 
  pol = Polarization(st, ncap, antiNu); 
  double Wm = pol.CalculateInverseMFP(E1)/Constants::HBCFmMeV*1.e13; 
  
  // Elastic approximation to the cross section  
  sig0 = 4.0*pow(Constants::GfMeV*(E1+deltaM)*Constants::HBCFmMeV, 2)
      / Constants::Pi; 
  // We divide by 16.0 instead of 4.0 because of our convention for the coupling 
  // constants
  sig0 *= (1.0*pow(ncap.Cv, 2) + 3.0*pow(ncap.Ca, 2))/4.0;
  sig0 *= sqrt(1.0 - pow(0.511/(E1+deltaM),2));
  elastic = sig0*st.n2/densFac*1.e13;

  
//  std::cout <<"charged current test: "<<n2<<" "<< noWm << " " << Wm << " " << (Wm/elastic-1.0)*M2/E1 << " " 
//      << elastic << " " << elastic/noWm - 1.0 << std::endl;
//  datacc << "charged current: " <<n2<<" "<< noWm << " " << Wm << " " << (Wm/elastic-1.0)*M2/E1 << " "
//      << elastic << " " << elastic/noWm - 1.0 << std::endl;
  
  if (fabs((noWm-elastic)/noWm)>1.e-2) return 1;
   
  return 0;
 }*/
  datacc.close();
  datanc.close();
  fclose(fpconfig);

}

