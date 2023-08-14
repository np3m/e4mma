/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2023, Xingfu Du, Zidu Lin, and Andrew W. Steiner
  
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
#include "eos.h"

#include <time.h>

#include <gsl/gsl_sf_hyperg.h>

#include <o2scl/xml.h>
#include <o2scl/format_float.h>
#include <o2scl/cloud_file.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/deriv_cern.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/fit_nonlin.h>
#include <o2scl/nstar_cold.h>
#include <o2scl/nucmass_densmat.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

double eos_crust_virial_v2::bn0_free() {
  return 0.0;
}

double eos_crust_virial_v2::bpn0_free() {
  return 0.0;
}

double eos_crust_virial_v2::bn1_free() {
  return -1.0/4.0/sqrt(2.0);
}

double eos_crust_virial_v2::bpn1_free() {
  return 0.0;
}

/// The temperature must be specified in MeV
double eos_crust_virial_v2::bna(double T) {
  if (T<1.0 || T>20.0) {
    cout << "Problem 1: " << T << endl;
    O2SCL_ERR("Cannot extrapolate spin virial.",o2scl::exc_efailed);
  }
  return iv_ba.eval(T);
}

/// The temperature must be specified in MeV
double eos_crust_virial_v2::bpna(double T) {
  if (T<1.0 || T>20.0) {
    cout << "Problem 2: " << T << endl;
    O2SCL_ERR("Cannot extrapolate spin virial.",o2scl::exc_efailed);
  }
  return iv_bpna.eval(T);
}

double eos_crust_virial_v2::f0(double lambda, double T) {
  return -lambda*lambda*lambda/2.0*T*
    (bn0(T)-bn0_free()+bn1(T)-bn1_free()+
     bpn0(T)-bpn0_free()+bpn1(T)-bpn1_free());
}

double eos_crust_virial_v2::f0p(double lambda, double T) {
  return -lambda*lambda*lambda/2.0*T*
    (bn0(T)-bn0_free()+bn1(T)-bn1_free()-
     bpn0(T)+bpn0_free()-bpn1(T)+bpn1_free());
}

double eos_crust_virial_v2::g0(double lambda, double T) {
  return lambda*lambda*lambda/2.0*T*
    (bn0(T)-bn0_free()-bn1(T)+bn1_free()+
     bpn0(T)-bpn0_free()-bpn1(T)+bpn1_free());
}

double eos_crust_virial_v2::g0p(double lambda, double T) {
  return lambda*lambda*lambda/2.0*T*
    (bn0(T)-bn0_free()-bn1(T)+bn1_free()-
     bpn0(T)+bpn0_free()+bpn1(T)-bpn1_free());
}

eos_crust_virial_v2::eos_crust_virial_v2() {
  bn_params.resize(10);
  bn_params[0]=2.874487202922e-01;
  bn_params[1]=2.200575070883e-03;
  bn_params[2]=-2.621025627694e-05;
  bn_params[3]=-6.061665959200e-08;
  bn_params[4]=1.059451872186e-02;
  bn_params[5]=5.673374476876e-02;
  bn_params[6]=3.492489364849e+00;
  bn_params[7]=-2.710552654167e-03;
  bn_params[8]=3.140521199464e+00;
  bn_params[9]=1.200987113605e+00;
  bpn_params.resize(6);
  bpn_params[0]=1.527316309589e+00;
  bpn_params[1]=1.748834077357e-04;
  bpn_params[2]=1.754991542102e+01;
  bpn_params[3]=4.510380054238e-01;
  bpn_params[4]=2.751333759925e-01;
  bpn_params[5]=-1.125035495140e+00;

  ba_T={1,2,3,4,5,6,7,8,9,10,12,14,16,18,20};
  vba={-0.638,-0.653,-0.651,-0.648,-0.643,-0.637,
      -0.631,-0.625,-0.620,-0.615,-0.6065,-0.597,
      -0.589,-0.583,-0.577};
  vbpna={6.18,1.74,1.05,0.785,0.640,0.561,
        0.504,0.463,0.432,0.408,0.374,0.352,
        0.336,0.324,0.315};
  iv_ba.set(ba_T.size(),ba_T,vba,itp_steffen);
  iv_bpna.set(ba_T.size(),ba_T,vbpna,itp_steffen);

  include_deuteron=false;
}

double eos_crust_virial_v2::bn_func(size_t np, const vector<double> &par,
				    double T) {
  return par[0]+par[1]*T+par[2]*T*T+par[3]*T*T*T+par[4]*
    exp(-par[5]*pow(T-par[6],2.0))+par[7]*exp(-par[8]*(T-par[9]));
}

double eos_crust_virial_v2::bpn_func(size_t np, const vector<double> &par,
				     double T) {
  double ret=par[0]*exp(-par[1]*(T+par[2])*(T+par[2]))+
    par[3]*exp(-par[4]*(T+par[5]));
  if (include_deuteron) {
    ret+=3.0/sqrt(2.0)*(exp(2.22/T)-1.0);
  }
  return ret;
}

double eos_crust_virial_v2::bn_f(double T) {
  return bn_func(bn_params.size(),bn_params,T);
}

double eos_crust_virial_v2::bpn_f(double T) {
  return bpn_func(bpn_params.size(),bpn_params,T);
}

double eos_crust_virial_v2::dbndT_f(double T) {
  return bn_params[1]+2.0*bn_params[2]*T+3.0*bn_params[3]*T*T-
    2.0*bn_params[4]*bn_params[5]*(T-bn_params[6])*
    exp(-bn_params[5]*pow(T-bn_params[6],2.0))-
    bn_params[7]*bn_params[8]*exp(-bn_params[8]*(T-bn_params[9]));
}

double eos_crust_virial_v2::dbpndT_f(double T) {
  return -bpn_params[0]*bpn_params[1]*2.0*(bpn_params[2]+T)*
    exp(-bpn_params[1]*(T+bpn_params[2])*(T+bpn_params[2]))-
    bpn_params[3]*bpn_params[4]*exp(-bpn_params[4]*(T+bpn_params[5]));
}

void eos_crust_virial_v2::fit(bool show_fit) {

  // Chi squared
  double chi2=0.0;

  // Fitter class
  fit_nonlin<chi_fit_funct<vector<double>,ubmatrix,std::function<
    double(size_t,const std::vector<double> &, double)> >,
	     vector<double>,ubmatrix> fitter;

  // --------------------------------------------
  // Fit neutron virial coefficient

  // Use the data from Horowitz and Schwenk (2006), PLB
  
  static const size_t neut_data=26;
  std::vector<double> Tv_neut(neut_data);
  Tv_neut[0]=0.1;
  Tv_neut[1]=0.5;
  Tv_neut[2]=1.0;
  Tv_neut[3]=2.0;
  Tv_neut[4]=3.0;
  Tv_neut[5]=4.0;
  Tv_neut[6]=5.0;
  Tv_neut[7]=6.0;
  Tv_neut[8]=7.0;
  Tv_neut[9]=8.0;
  Tv_neut[10]=9.0;
  Tv_neut[11]=10.0;
  Tv_neut[12]=12.0;
  Tv_neut[13]=14.0;
  Tv_neut[14]=16.0;
  Tv_neut[15]=18.0;
  Tv_neut[16]=20.0;
  Tv_neut[17]=22.0;
  Tv_neut[18]=24.0;
  Tv_neut[19]=25.0;
  Tv_neut[20]=30.0;
  Tv_neut[21]=35.0;
  Tv_neut[22]=40.0;
  Tv_neut[23]=45.0;
  Tv_neut[24]=50.0;
  Tv_neut[25]=150.0;

  bnv.resize(neut_data);
  bnv[0]=0.207;
  bnv[1]=0.272;
  bnv[2]=0.288;
  bnv[3]=0.303;
  bnv[4]=0.306;
  bnv[5]=0.306;
  bnv[6]=0.306;
  bnv[7]=0.306;
  bnv[8]=0.307;
  bnv[9]=0.307;
  bnv[10]=0.308;
  bnv[11]=0.309;
  bnv[12]=0.310;
  bnv[13]=0.313;
  bnv[14]=0.315;
  bnv[15]=0.318;
  bnv[16]=0.320;
  bnv[17]=0.322;
  bnv[18]=0.324;
  bnv[19]=0.325;
  bnv[20]=0.329;
  bnv[21]=0.330;
  bnv[22]=0.330;
  bnv[23]=0.328;
  bnv[24]=0.324;
  bnv[25]=-pow(2.0,-5.0/2.0);
  
  vector<double> bn_err(neut_data);
  for(size_t i=0;i<neut_data;i++) {
    bn_err[i]=fabs(bnv[i])/1.0e2;
  }

  // Fitting function
  std::function<
    double(size_t,const std::vector<double> &, double)> ff_neutron=
    std::bind(std::mem_fn<double(size_t,const vector<double> &,double)>
	      (&eos_crust_virial_v2::bn_func),
	      this,std::placeholders::_1,std::placeholders::_2,
	      std::placeholders::_3);
  
  // Chi-squared and fitting data
  chi_fit_funct<vector<double>,ubmatrix,std::function<
    double(size_t,const std::vector<double> &, double)> > 
    cff(neut_data,Tv_neut,bnv,bn_err,ff_neutron);
  
  cout << "Neutron virial coefficient:\n" << endl;
  
  ubmatrix covar(bn_np,bn_np);

  cout << "Initial chi-squared: " << cff.chi2(bn_np,bn_params) << endl;
  
  fitter.fit(bn_np,bn_params,covar,chi2,cff);

  if (show_fit) {
    cout << "Final chi-squared: " << chi2 << endl;
    cout << "params: " << endl;
    cout.precision(12);
    for(size_t j=0;j<bn_np;j++) cout << "bn_params[" << j << "]="
				     << bn_params[j] << ";" << endl;
    cout.precision(6);
    cout << endl;
    
    table<> t;
    t.line_of_names("T bn bn_err bn_fit");
    for(size_t j=0;j<neut_data;j++) {
      cout << Tv_neut[j] << " " << bnv[j] << " " << bn_err[j] << " "
	   << ff_neutron(bn_np,bn_params,Tv_neut[j]) << endl;
      double line[4]={Tv_neut[j],bnv[j],bn_err[j],
		      ff_neutron(bn_np,bn_params,Tv_neut[j])};
      t.line_of_data(4,line);
    }
    cout << endl;

    hdf_file hf;
    hf.open_or_create("fit_neut.o2");
    hdf_output(hf,t,"fit_neut");
    hf.close();
  }
  
  // --------------------------------------------
  // Fit neutron-proton virial coefficient
  
  cout << "Neutron-proton virial coefficient:\n" << endl;

  static const size_t nuc_data=17;

  std::vector<double> Tv_nuc(nuc_data);
  Tv_nuc[0]=0.1;
  Tv_nuc[1]=1.0;
  Tv_nuc[2]=2.0;
  Tv_nuc[3]=3.0;
  Tv_nuc[4]=4.0;
  Tv_nuc[5]=5.0;
  Tv_nuc[6]=6.0;
  Tv_nuc[7]=7.0;
  Tv_nuc[8]=8.0;
  Tv_nuc[9]=9.0;
  Tv_nuc[10]=10.0;
  Tv_nuc[11]=12.0;
  Tv_nuc[12]=14.0;
  Tv_nuc[13]=16.0;
  Tv_nuc[14]=18.0;
  Tv_nuc[15]=20.0;
  Tv_nuc[16]=150.0;
  
  ubmatrix covar2(bpn_np,bpn_np);

  // T=0.1 MeV point from effective range expansion
  bpnv[0]=2.046;

  // Data from Horowitz paper
  bpnv[1]=19.4;
  bpnv[2]=6.1;
  bpnv[3]=4.018;
  bpnv[4]=3.19;
  bpnv[5]=2.74;
  bpnv[6]=2.46;
  bpnv[7]=2.26;
  bpnv[8]=2.11;
  bpnv[9]=2.00;
  bpnv[10]=1.91;
  bpnv[11]=1.76;
  bpnv[12]=1.66;
  bpnv[13]=1.57;
  bpnv[14]=1.51;
  bpnv[15]=1.45;

  // T=150 MeV point setting virial coefficient to zero
  bpnv[16]=0.0;

  // Subtract off deuteron contribution
  for(size_t i=1;i<16;i++) {
    bpnv[i]-=3.0/sqrt(2.0)*(exp(2.224/Tv_nuc[i])-1.0);
  }

  // Determine uncertainties
  std::vector<double> bpn_err(nuc_data);
  for(size_t i=0;i<16;i++) {
    bpn_err[i]=fabs(bpnv[i])/1.0e2;
  }
  bpn_err[16]=0.04;

  // Fitting function
  std::function<
    double(size_t,const std::vector<double> &, double)> ff_nuc=
    std::bind(std::mem_fn<double(size_t,const vector<double> &,double)>
	      (&eos_crust_virial_v2::bpn_func),
	      this,std::placeholders::_1,std::placeholders::_2,
	      std::placeholders::_3);
  
  // Chi-squared and fitting data
  chi_fit_funct<vector<double>,ubmatrix,std::function<
    double(size_t,const std::vector<double> &, double)> > 
    cff_nuc(nuc_data,Tv_nuc,bpnv,bpn_err,ff_nuc);
  
  cout << "Initial chi-squared: " << cff_nuc.chi2(bpn_np,bpn_params) << endl;

  fitter.fit(bpn_np,bpn_params,covar2,chi2,cff_nuc);

  if (show_fit) {
    cout << "Final chi-squared: " << chi2 << endl;
    cout << "params: " << endl;
    cout.precision(12);
    for(size_t j=0;j<bpn_np;j++) cout << "bpn_params[" << j << "]="
				      << bpn_params[j] << ";" << endl;
    cout.precision(6);
    cout << endl;
    
    table<> t;
    t.line_of_names("T bpn bpn_err bpn_fit");
    for(size_t j=0;j<nuc_data;j++) {
      cout << Tv_nuc[j] << " " << bpnv[j] << " " << bpn_err[j] << " "
	   << ff_nuc(bpn_np,bpn_params,Tv_nuc[j]) << endl;
      double line[4]={Tv_nuc[j],bpnv[j],bpn_err[j],
		      ff_nuc(bpn_np,bpn_params,Tv_nuc[j])};
      t.line_of_data(4,line);
    }

    hdf_file hf;
    hf.open_or_create("fit_nuc.o2");
    hdf_output(hf,t,"fit_nuc");
    hf.close();
  }
  
  // End of eos_crust_virial_v2::fit()
  return;
}

int eos::hrg_load(std::vector<std::string> &sv, bool itive_com) {

  std::string fn=sv[1];

  ifstream fin;
  fin.open(fn.c_str());
  std::string stmp;

  for(size_t k=0;k<12;k++) fin >> stmp;

  int dtemp;
  while (fin >> dtemp) { 
    part_pdg_db::pdg_entry p;
    p.id=dtemp;
    fin >> p.name >> p.mass >> p.mass_errp >> p.spin_deg >>
      p.baryon >> p.strangeness >> stmp >> stmp >> stmp >> p.charge;
    //cout << "Read " << p.name << " " << p.mass << " " << p.charge << endl;
    getline(fin,stmp);
    p.mass_errm=p.mass_errp;
    p.width=0.0;
    p.width_errp=0.0;
    p.width_errm=0.0;
    part_db.push_back(p);
  }
  
  fin.close();

  cout << "Read " << part_db.size() << " resonances." << endl;

  int nferm=0;
  int nbos=0;
  for(size_t j=0;j<part_db.size();j++) {
    // Skip the photon because it's a boson with two degrees of freedom
    if (part_db[j].spin_deg%2==0 && part_db[j].id!=22) {
      fermion f;
      f.m=part_db[j].mass*1.0e3/hc_mev_fm;
      f.g=part_db[j].spin_deg;
      f.n=0.0;
      f.ed=0.0;
      f.en=0.0;
      f.pr=0.0;
      f.mu=0.0;
      f.nu=0.0;
      f.ms=0.0;
      res_f.push_back(f);
      nferm++;
      //cout << j << " " << part_db[j].id << " "
      //<< nferm << " " << f.m*hc_mev_fm << " " << f.g << endl;
    } else {
      boson b;
      b.m=part_db[j].mass*1.0e3/hc_mev_fm;
      b.g=part_db[j].spin_deg;
      b.n=0.0;
      b.ed=0.0;
      b.en=0.0;
      b.pr=0.0;
      b.mu=0.0;
      b.nu=0.0;
      b.ms=0.0;
      res_b.push_back(b);
      nbos++;
    }
  }

  cout << "Read " << nferm << " fermions and " << nbos << " bosons."
       << endl;
  
  return 0;
}

double eos::fit_fun(size_t np, const std::vector<double> &parms,
		    double nb) {
  if (old_ns_fit) {
    return (sqrt(nb)*parms[0]+nb*parms[1]+
	    nb*sqrt(nb)*parms[2]+nb*nb*parms[3]+
	    nb*nb*nb*parms[4]);
  }
  return (nb*parms[0]+nb*nb*parms[1]+nb*nb*nb*parms[2]+
	  nb*nb*nb*nb*parms[3]+nb*nb*nb*nb*nb*parms[4]);
}

double eos::ed_fit(double nb) {
  return fit_fun(5,ns_fit_parms,nb)*nb/hc_mev_fm;
}

double eos::mu_fit(double nb) {
  if (old_ns_fit) {
    return (1.5*sqrt(nb)*ns_fit_parms[0]+2.0*nb*ns_fit_parms[1]+
	    2.5*nb*sqrt(nb)*ns_fit_parms[2]+3.0*nb*nb*ns_fit_parms[3]+
	    4.0*nb*nb*nb*ns_fit_parms[4])/hc_mev_fm;
  }
  return (2.0*nb*ns_fit_parms[0]+3.0*nb*nb*ns_fit_parms[1]+
	  4.0*nb*nb*nb*ns_fit_parms[2]+
	  5.0*nb*nb*nb*nb*ns_fit_parms[3]+
	  6.0*nb*nb*nb*nb*nb*ns_fit_parms[4])
    /hc_mev_fm;
}

double eos::dmudn_fit(double nb) {
  if (old_ns_fit) {
    return (0.75/sqrt(nb)*ns_fit_parms[0]+2.0*ns_fit_parms[1]+
	    3.75*sqrt(nb)*ns_fit_parms[2]+6.0*nb*ns_fit_parms[3]+
	    12.0*nb*nb*ns_fit_parms[4])/hc_mev_fm;
  }
  return (2.0*ns_fit_parms[0]+6.0*nb*ns_fit_parms[1]+
	  12.0*nb*nb*ns_fit_parms[2]+20.0*nb*nb*nb*
	  ns_fit_parms[3]+30.0*nb*nb*nb*nb*ns_fit_parms[4])/hc_mev_fm;
}

double eos::cs2_fit(double nb) {
  return dmudn_fit(nb)*nb/(mu_fit(nb)+939.565/hc_mev_fm);
}

void eos::min_max_cs2(double &cs2_min, double &cs2_max) {
  double nb=0.08;
  double cs2=cs2_fit(nb);
  cs2_min=cs2;
  cs2_max=cs2;
  for(nb=0.04;nb<ns_nb_max;nb+=0.02) {
    cs2=cs2_fit(nb);
    if (cs2<cs2_min) cs2_min=cs2;
    if (cs2>cs2_max) cs2_max=cs2;
  }
  return;
}

void eos::ns_fit(int row) {

  if (row>=((int)(nstar_tab.get_nlines()))) {
    O2SCL_ERR("Row not allowed in ns_fit().",
	      exc_efailed);
  }

  // Set the class data member to the appropriate row
  i_ns=row;
  
  Lambda_bar_14=nstar_tab.get("lambda_bar14",row);
  
  // Initial values of the parameters
  ns_fit_parms.resize(5);
  ns_fit_parms[0]=-6.102748e3;
  ns_fit_parms[1]=3.053497e3;
  ns_fit_parms[2]=4.662834e3;
  ns_fit_parms[3]=-8.371958e2;
  ns_fit_parms[4]=-5.228209e2;

  // This table stores the neutron star EOS in a form that
  // can be used by the fitting object 
  table_units<> nstar_high;  
  nstar_high.line_of_names("nb EoA");
  ns_nb_max=nstar_tab.get((string)"nb_max",i_ns);
  for(size_t i=0;i<100 && (0.04+((double)i)*0.012)<(ns_nb_max+0.000001);i++) {
    double EoA=nstar_tab.get(((string)"EoA_")+szttos(i),i_ns)*hc_mev_fm;
    if (fabs(nstar_tab.get(((string)"EoA_")+szttos(i),i_ns))>1.0e-2) {
      double line[2]={0.04+((double)i)*0.012,
		      nstar_tab.get(((string)"EoA_")+szttos(i),i_ns)*
		      hc_mev_fm};
      nstar_high.line_of_data(2,line);
    }
  }
  
  // Fit the energy per baryon 
  nstar_high.function_column("abs(EoA)/100","Eerr");
  typedef std::function<
    double(size_t,const std::vector<double> &, double)> fit_funct2;
  fit_funct2 ff=std::bind
    (std::mem_fn<double(size_t,const std::vector<double> &,double)>
     (&eos::fit_fun),this,std::placeholders::_1,
     std::placeholders::_2,std::placeholders::_3);

  // Collect vector references to specify the data as
  // a fit function
  size_t ndat=nstar_high.get_nlines();

  if (false) {
    //plot nstarhigh 
    hdf_file nshf;
    nshf.open_or_create("nbcs2.o2");
    hdf_output(nshf,nstar_high,"nbcs2");
    nshf.close();
    //int sret=system("o2graph -read nbcs2.o2 -plot nb EoA -show");
  }

  const std::vector<double> &xdat=nstar_high["nb"];
  const std::vector<double> &ydat=nstar_high["EoA"];
  const std::vector<double> &yerr=nstar_high["Eerr"];

  // Set up the fitting function (just a simple chi^2 fit)
  size_t nparms=5;
  chi_fit_funct<std::vector<double>,ubmatrix,fit_funct2> 
    cff(ndat,xdat,ydat,yerr,ff);
  
  // The fit covariance matrix
  ubmatrix covar(nparms,nparms);

  // The fit chi-squared
  if (verbose>0) {
    cout << "ns_fit_parms[0]=" << ns_fit_parms[0] << ";" << endl;
    cout << "ns_fit_parms[1]=" << ns_fit_parms[1] << ";" << endl;
    cout << "ns_fit_parms[2]=" << ns_fit_parms[2] << ";" << endl;
    cout << "ns_fit_parms[3]=" << ns_fit_parms[3] << ";" << endl;
    cout << "ns_fit_parms[4]=" << ns_fit_parms[4] << ";" << endl;
  }

  // Perform the fit
  fit_nonlin<chi_fit_funct<std::vector<double>,ubmatrix,fit_funct2>,
	     std::vector<double>,ubmatrix> fn;
  fn.fit(nparms,ns_fit_parms,covar,chi2_ns,cff);
  
  // Store the results of the fit in column named "EoA_fit", then compute
  // energy density and chemical potential//
  nstar_high.new_column("EoA_fit");
  nstar_high.new_column("ed_fit");
  nstar_high.new_column("mu_fit");
  nstar_high.new_column("cs2_fit");
  for(size_t i=0;i<ndat;i++) {
    double nb=nstar_high.get("nb",i);
    double EoA=ff(nparms,ns_fit_parms,nb);
    nstar_high.set("EoA_fit",i,EoA);
    nstar_high.set("ed_fit",i,(EoA+939.0)/hc_mev_fm*nstar_high.get("nb",i));
    double mu=mu_fit(nb)+939.0/hc_mev_fm;
    nstar_high.set("mu_fit",i,mu);
    double cs2=cs2_fit(nb);
    nstar_high.set("cs2_fit",i,cs2);
  }
  nstar_high.function_column("(EoA+939)/197.33*nb","ed");
  nstar_high.deriv("nb","ed","mu");
  nstar_high.deriv2("nb","ed","dmudn");
  nstar_high.function_column("nb/mu*dmudn","cs2");

  // Readjust ns_nb_max to ensure it's lower than the point
  // at which c_s^2 becomes superluminal
  double nb_new=0.0;
  for(size_t j=0;j<nstar_high.get_nlines()-1;j++) {
    if (nstar_high.get("cs2_fit",j)<1.0 &&
	nstar_high.get("cs2_fit",j+1)>1.0) {
      nb_new=nstar_high.get("nb",j)+
	(nstar_high.get("nb",j+1)-nstar_high.get("nb",j))*
	(1.0-nstar_high.get("cs2_fit",j))/
	(nstar_high.get("cs2_fit",j+1)-nstar_high.get("cs2_fit",j));
    }
  }
  if (nb_new>0.01) ns_nb_max=nb_new;
  
  // Output the fit results to the screen
  if (verbose>0) {
    cout << "Parameters: " << endl;
    cout << "ns_fit_parms[0]=" << ns_fit_parms[0] << ";" << endl;
    cout << "ns_fit_parms[1]=" << ns_fit_parms[1] << ";" << endl;
    cout << "ns_fit_parms[2]=" << ns_fit_parms[2] << ";" << endl;
    cout << "ns_fit_parms[3]=" << ns_fit_parms[3] << ";" << endl;
    cout << "ns_fit_parms[4]=" << ns_fit_parms[4] << ";" << endl;
    cout << "chi2: " << chi2_ns << endl;
  }

  // If true, record the results of the fit
  if (ns_record) {
    table_units<> tab2=nstar_high;
    // If ns_record is true, this file stores the original neutron star
    // EOS and the results of the fit
    o2scl_hdf::hdf_file hf;
    hf.open_or_create("ns_fit.o2");
    hdf_output(hf,tab2,"ns_fit");
    hf.close();
    exit(-1);
  }
  
  return;
}

eos::eos() {

  cs2_verbose=0;
  
  // Nucleon init
  neutron.init(o2scl_settings.get_convert_units().convert
	       ("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  proton.init(o2scl_settings.get_convert_units().convert
	      ("kg","1/fm",o2scl_mks::mass_proton),2.0);
  neutron.non_interacting=false;
  proton.non_interacting=false;
  neutron.inc_rest_mass=false;
  proton.inc_rest_mass=false; 

  // Nucleon init (take 2)
  n_chiral.init(o2scl_settings.get_convert_units().convert
		("kg","1/fm",o2scl_mks::mass_neutron),2.0);
  p_chiral.init(o2scl_settings.get_convert_units().convert
		("kg","1/fm",o2scl_mks::mass_proton),2.0);
  n_chiral.non_interacting=false;
  p_chiral.non_interacting=false;
  n_chiral.inc_rest_mass=false;
  p_chiral.inc_rest_mass=false;

  // Lepton and photon inits
  electron.init(o2scl_settings.get_convert_units().convert
		("kg","1/fm",o2scl_mks::mass_electron),2.0);
  muon.init(o2scl_settings.get_convert_units().convert
	    ("kg","1/fm",o2scl_mks::mass_muon),2.0);
  neutrino.init(0.0,1.0);
  photon.init(0.0,2.0);

  // Default settings
  verbose=0;
  test_ns_cs2=false;
  include_muons=false;
  output_files=true; 
  old_ns_fit=true;
  ns_record=false;
  model_selected=false;
  select_cs2_test=true;
  a_virial=10.0;
  b_virial=10.0;

  // Initial parameter values
  i_ns=-1;
  i_skyrme=-1;
  qmc_alpha=0.48;
  qmc_beta=3.45;
  qmc_a=12.7;
  qmc_b=2.12;
  qmc_n0=0.16;

  // Temporary string for object names
  string name;

  int mpi_rank=0, mpi_size=1;

#ifdef O2SCL_MPI
  // Get MPI rank, etc.
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  
  // Ensure that multiple MPI ranks aren't reading from the
  // filesystem at the same time
  int tag=0, buffer=0;
  if (mpi_size>1 && mpi_rank>=1) {
    MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,
	     tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
#endif

  // Open the neutron star data file
  std::string ns_file="data/qmc_twop_10_0_out";
  o2scl_hdf::hdf_file hf;
  hf.open(ns_file);
  o2scl_hdf::hdf_input(hf,nstar_tab,name);
  hf.close();

  // Open the Skyrme data file
  std::string UNEDF_file="data/thetaANL-1002x12.o2";
  hf.open(UNEDF_file);
  o2scl_hdf::hdf_input(hf,UNEDF_tab,name);
  hf.close();

  use_alt_eos=false;
  o2scl_hdf::skyrme_load(sk_alt,"NRAPR");

#ifdef O2SCL_MPI
  // Send a message to the next MPI rank
  if (mpi_size>1 && mpi_rank<mpi_size-1) {
    MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
	     tag,MPI_COMM_WORLD);
  }
#endif

  nB_grid_spec="301,10^(i*0.04-12)*2.0";
  Ye_grid_spec="70,0.01*(i+1)";
  T_grid_spec="160,0.1*1.046^i";
  S_grid_spec="5,0.1*i";

  // Skyrme couplings
  sk_Tcorr.t0=5.067286719233e+03;
  sk_Tcorr.t1=1.749251370992e+00;
  sk_Tcorr.t2=-4.721193938990e-01;
  sk_Tcorr.t3=-1.945964529505e+05;
  sk_Tcorr.x0=4.197555064408e+01;
  sk_Tcorr.x1=-6.947915483747e-02;
  sk_Tcorr.x2=4.192016722695e-01;
  sk_Tcorr.x3=-2.877974634128e+01;
  sk_Tcorr.alpha=0.144165;

  elep.include_photons=true;
  elep.include_muons=false;

  // Seed the random number generator with the clock time
  rng.clock_seed();

  old_version=false;

  eos_Tcorr=&sk_Tcorr;
  eosp_alt=&sk_alt;
  strange_axis=false;
  alt_name="";
  rmf_fields=false;
}

double eos::energy_density_qmc(double nn, double np) {
  
  double e_qmc=(qmc_a*pow((nn+np)/qmc_n0,qmc_alpha)+
		qmc_b*pow((nn+np)/qmc_n0,qmc_beta))*
    (nn+np)/hc_mev_fm;
  
  return e_qmc;
}
  
double eos::energy_density_ns(double nn) {
  return ed_fit(nn);
}    

double eos::free_energy_density_virial
(fermion &n, fermion &p, double T, thermo &th, double &dmundnn,
 double &dmundnp, double &dmupdnn, double &dmupdnp, double &dmundT,
 double &dmupdT) {
 
  double nn=n.n;
  double np=p.n;
    
  double T_MeV=T*hc_mev_fm;
  
  double b_n=ecv.bn_f(T_MeV);
  double dbndT=ecv.dbndT_f(T_MeV)*hc_mev_fm;
  double b_pn=ecv.bpn_f(T_MeV);
  double dbpndT=ecv.dbpndT_f(T_MeV)*hc_mev_fm;

  double lambda=sqrt(4.0*o2scl_const::pi/(n.m+p.m)/T);
  double dlambdadT=-sqrt(o2scl_const::pi/(n.m+p.m))/pow(sqrt(T),3.0);

  double zn, zp;

  if (nn*pow(lambda,3.0)>1.0e-5 || np*pow(lambda,3.0)>1.0e-5) {

    // If the densities are large enough, then compute the
    // virial result
    vsd.solve_fugacity(nn,np,lambda,lambda,b_n,b_pn,zn,zp);
    
    vsd.calc_deriv(nn,np,lambda,lambda,
		   b_n,b_pn,zn,zp,
		   dbndT,dbpndT,dlambdadT,dlambdadT);
    
    dmundnn=T/zn*vsd.dzndnn;
    dmundnp=T/zn*vsd.dzndnp;
    dmupdnn=T/zp*vsd.dzpdnn;
    dmupdnp=T/zp*vsd.dzpdnp;
    dmundT=log(zn)+T/zn*vsd.dzndT;
    dmupdT=log(zp)+T/zp*vsd.dzpdT;
    
    double zn_old=zn, zp_old=zp;
    double dmundnn_old=dmundnn;
    double dmundnp_old=dmundnp;
    double dmupdnn_old=dmupdnn;
    double dmupdnp_old=dmupdnp;
    double dmundT_old=dmundT;
    double dmupdT_old=dmupdT;

    if (false) {
    
      if (fabs(zn-zn_old)/fabs(zn_old)>1.0e-5 ||
	  fabs(zp-zp_old)/fabs(zp_old)>1.0e-5) {
	cout << "Disagreement." << endl;
	cout << "nn,np,TMeV: " << n.n << " " << p.n << " "
	     << T*hc_mev_fm << endl;
	cout << "zn,zp(old): " << zn_old << " " << zp_old << endl;
	cout << "zn,zp(new): " << zn << " " << zp << endl;
	exit(-1);
      }
    
      if (fabs(dmundnn-dmundnn_old)/fabs(dmundnn_old)>1.0e-5 ||
	  fabs(dmundnp-dmundnp_old)/fabs(dmundnp_old)>1.0e-5) {
	cout << "Disagreement 2." << endl;
	cout << "nn,np,TMeV: " << n.n << " " << p.n << " "
	     << T*hc_mev_fm << endl;
	cout << "dmundnn,dmundnp(old): " << dmundnn_old << " "
	     << dmundnp_old << endl;
	cout << "dmundnn,dmundnp(new): " << dmundnn << " "
	     << dmundnp << endl;
	exit(-1);
      }
    
      if (fabs(dmupdnn-dmupdnn_old)/fabs(dmupdnn_old)>1.0e-5 ||
	  fabs(dmupdnp-dmupdnp_old)/fabs(dmupdnp_old)>1.0e-5) {
	cout << "Disagreement 3." << endl;
	cout << "nn,np,TMeV: " << n.n << " " << p.n << " "
	     << T*hc_mev_fm << endl;
	cout << "zn,zp(old): " << zn_old << " " << zp_old << endl;
	cout << "zn,zp(new): " << zn << " " << zp << endl;
	exit(-1);
      }
    
      if (fabs(dmupdT-dmupdT_old)/fabs(dmupdT_old)>1.0e-5 ||
	  fabs(dmupdT-dmupdT_old)/fabs(dmupdT_old)>1.0e-5) {
	cout << "Disagreement 4." << endl;
	cout << "nn,np,TMeV: " << n.n << " " << p.n << " "
	     << T*hc_mev_fm << endl;
	cout << "zn,zp(old): " << zn_old << " " << zp_old << endl;
	cout << "zn,zp(new): " << zn << " " << zp << endl;
	exit(-1);
      }

    }
    
  } else {
    
    // Otherwise, the virial correction is negligable, so
    // just use the classical result
    zn=nn*pow(lambda,3.0)/2.0;
    zp=np*pow(lambda,3.0)/2.0;
    
    dmundnn=T*pow(lambda,3.0)/2.0/zn;
    dmundnp=0.0;
    dmupdnn=0.0;
    dmupdnp=T*pow(lambda,3.0)/2.0/zp;
    dmundT=n.mu/T-1.5;
    dmupdT=p.mu/T-1.5;
    
  }

  n.mu=log(zn)*T;
  p.mu=log(zp)*T;
  
  if (nn*pow(lambda,3.0)>1.0e-5 || np*pow(lambda,3.0)>1.0e-5) {
  
    th.pr=2.0*T/pow(lambda,3.0)*(zn+zp+(zn*zn+zp*zp)*b_n+
				 2.0*zp*zn*b_pn);
    th.en=5.0*th.pr/2.0/T-nn*log(zn)-np*log(zp)+
      2.0*T/pow(lambda,3.0)*((zn*zn+zp*zp)*dbndT
			     +2.0*zp*zn*dbpndT);
  } else {

    th.pr=2.0*T/pow(lambda,3.0)*(zn+zp);
    th.en=5.0*th.pr/2.0/T-nn*log(zn)-np*log(zp);

  }

  double f_vir=n.mu*nn+p.mu*np-th.pr;
  
  th.ed=f_vir+T*th.en;

  return f_vir;
}

int eos::solve_coeff_big(size_t nv, const ubvector &x, ubvector &y, 
			 double ns_nb_max_l, double cs_ns_2,
			 double cs_ns_last) {
  double a1l=x[0];
  double a2l=x[1];
  y[0]=1.0-a1l+(a1l*a2l*pow(ns_nb_max_l,a1l))/
    (1.0+a2l*pow(ns_nb_max_l,a1l))-cs_ns_last;
  y[1]=1.0-a1l+(a1l*a2l*pow(2.0,a1l))/(1.0+a2l*pow(2.0,a1l))-cs_ns_2;
  
  return 0;
}

int eos::solve_coeff_small(size_t nv, const ubvector &x,
			   ubvector &y, double ns_nb_max_l,
			   double cs_ns_2, double cs_ns_last) {
                                    
  double a1l=x[0];
  double a2l=x[1];
  y[0]=(a1l-a1l*a2l*pow(ns_nb_max_l,a1l)/
        (1.0+a2l*pow(ns_nb_max_l,a1l)))-cs_ns_last;
  y[1]=(a1l-a1l*a2l*pow(2.0,a1l)/(1.0+a2l*pow(2.0,a1l)))-cs_ns_2;
  
  return 0;
}

int eos::new_ns_eos(double nb, fermion &n,
		    double &e_ns, double &densdnn) {
  
  double cs_ns_last;
  double cs_ns_2;
  double a1l, a2l;
  double c1l, c2l;

  double e_ns_last=ed_fit(ns_nb_max);
  double p_ns_last=mu_fit(ns_nb_max)*ns_nb_max-e_ns_last;

  // -----------------------------------------------------
  // Solve for a1l and a2l

  bool success=true;
  mroot_hybrids<> mh;
  mh.err_nonconv=false;

  // Initial guess for a1l and a2l
  ubvector mx(2), my(2);
  
  if (nb<(ns_nb_max-1.0e-6)) {
    
    // If we're in the region that the neutron star
    // EOS is causal, just use that result
    e_ns=energy_density_ns(nb);
    
    densdnn=mu_fit(nb);
    
  } else {

    // If the speed of sound is increasing
    // at high densities
    
    cs_ns_last=cs2_fit(ns_nb_max);
    cs_ns_2=phi;
    
    if (cs_ns_2>cs_ns_last) {
      mx[0]=1.0;
      mx[1]=1.0;
      mm_funct mfbig=std::bind
	(std::mem_fn<int(size_t,const ubvector &,
			 ubvector &, double, double, double)>
	 (&eos::solve_coeff_big),
	 this,std::placeholders::_1,
	 std::placeholders::_2,
	 std::placeholders::_3,ns_nb_max,cs_ns_2,cs_ns_last);
      int mret=mh.msolve(2,mx,mfbig);
      if (mret!=0) success=false;
      
      a1l=mx[0];
      a2l=mx[1];
      
      // solve for c1, c2
      c1l=(e_ns_last+n.m*ns_nb_max+p_ns_last)/((ns_nb_max*ns_nb_max)*
					       (a2l+pow(ns_nb_max,-a1l)));
      c2l=0.5*(e_ns_last+n.m*ns_nb_max-p_ns_last+
	       a1l*(e_ns_last+n.m*ns_nb_max+p_ns_last)
	       /((a1l-2.0)*(1.0+a2l*pow(ns_nb_max,a1l))));
      e_ns=-n.m*nb+(a2l*nb*nb/2.0+pow(nb,2.0-a1l)/(2.0-a1l))*c1l+c2l;
      densdnn=-n.m+c1l*(a2l*nb+pow(nb,1.0-a1l));

    } else if (cs_ns_2<cs_ns_last) {

      // If the speed of sound is decreasing
      // at high densities

      mx[0]=2.5;
      mx[1]=1.0;
      mm_funct mfsmall=std::bind
	(std::mem_fn<int(size_t,const ubvector &,
			 ubvector &, double, double ,double)>
	 (&eos::solve_coeff_small),
	 this,std::placeholders::_1,
	 std::placeholders::_2,
	 std::placeholders::_3,ns_nb_max,cs_ns_2,cs_ns_last);
      int mret=mh.msolve(2,mx,mfsmall);
      if (mret!=0) success=false;
      a1l=mx[0];
      a2l=mx[1];
      
      double hyperg=gsl_sf_hyperg_2F1
	(1.0,1.0,1.0-1.0/a1l,1.0/a2l*pow(nb,-a1l)/(1.0/a2l*pow(nb,-a1l)+1.0));
      double hyperg_max=gsl_sf_hyperg_2F1
	(1.0,1.0,1.0-1.0/a1l,1.0/a2l*pow(ns_nb_max,-a1l)/
	 (1.0/a2l*pow(ns_nb_max,-a1l)+1.0));
      
      // Transform hyperg to hyperg_new see notes based on Pfaff
      // transformation
      double hyperg_new=hyperg*(1.0/(1.0+1.0/a2l*pow(nb,-a1l)));
      double hyperg_max_new=hyperg_max*(1.0/(1.0+1.0/a2l*pow(ns_nb_max,-a1l)));

      // solve for c1l, c2l
      c1l=pow(ns_nb_max,-a1l-1.0)*(a2l*pow(ns_nb_max,a1l)+1.0)*
	(e_ns_last+n.m*ns_nb_max+p_ns_last);
      c2l=pow(ns_nb_max,-a1l)*
	(a2l*pow(ns_nb_max,a1l)*(e_ns_last+n.m*ns_nb_max)-
	 (a2l*pow(ns_nb_max,a1l)+1.0)
	 *hyperg_max_new*(e_ns_last+n.m*ns_nb_max+p_ns_last))/a2l;
      e_ns=(c1l*nb*hyperg_new)/a2l+c2l-n.m*nb;
      densdnn=-(a2l*n.m*pow(nb,a1l)-c1l*pow(nb,a1l)+n.m)/(a2l*pow(nb,a1l)+1.0);

    } else if (cs_ns_2==cs_ns_last) {

      // If the speed of sound is independent of density at
      // high densities
      
      e_ns=-n.m*nb+(e_ns_last+n.m*ns_nb_max+p_ns_last)
	/(1.0+cs_ns_last)*pow((nb/ns_nb_max),cs_ns_last+1.0)
	+(cs_ns_last*(e_ns_last+n.m*ns_nb_max)-p_ns_last)
	/(1.0+cs_ns_last);
      densdnn=-n.m+(e_ns_last+n.m*ns_nb_max+p_ns_last)
	*pow((nb/ns_nb_max),cs_ns_last)/ns_nb_max;
    }
  }

  if (success==false) {
    cout << "Warning NS EOS solver failing." << endl;
    return 1;
  }
  return 0;
}

double eos::free_energy_density
(fermion &n, fermion &p, double T, thermo &th) {
  std::map<std::string,double> vdet;
  return free_energy_density_detail(n,p,T,th,vdet);
}

double eos::free_energy_density_s
(fermion &n, fermion &p, double Y_s, double T, thermo &th) {
  std::map<std::string,double> vdet;
  return free_energy_density_detail_s(n,p,Y_s,T,th,vdet);
}

double eos::free_energy_density_detail_s
(o2scl::fermion &n, o2scl::fermion &p, double Y_s, double T, o2scl::thermo &th,
 std::map<std::string,double> &vdet) {

  if (!strange_axis) {
    return free_energy_density_detail(n,p,T,th,vdet);
  }
  
  if (use_alt_eos) {
    rmf_hyp.calc_hyp_e_nobeta_np(Y_s,n,p,rmf_hyp.def_lambda,
                                 rmf_hyp.def_sigma_p,
                                 rmf_hyp.def_sigma_z,
                                 rmf_hyp.def_sigma_m,
                                 rmf_hyp.def_cascade_z,
                                 rmf_hyp.def_cascade_m,th);
  } else {
    O2SCL_ERR("Strangness only supported for alternate EOSs.",
              o2scl::exc_eunimpl);
  }
  
  return 0.0;
}

double eos::free_energy_density_detail
(o2scl::fermion &n, o2scl::fermion &p, double T, o2scl::thermo &th,
 std::map<std::string,double> &vdet) {
  
  double zn;
  double zp;
  double f1;
  double f2;
  double f3;
  double f4;
  double g_virial;
  double dgvirialdT;
  double dgvirialdnn;
  double dgvirialdnp;
  
  if (use_alt_eos) {

    if (eosp_alt==&rmf_hyp) {
      rmf_hyp.calc_temp_hyp_e(n.n+p.n,p.n,n,p,rmf_hyp.def_lambda,
                              rmf_hyp.def_sigma_p,
                              rmf_hyp.def_sigma_z,
                              rmf_hyp.def_sigma_m,
                              rmf_hyp.def_cascade_z,
                              rmf_hyp.def_cascade_m,
                              T,th);
    } else {
      eosp_alt->calc_temp_e(n,p,T,th);
    }
    zn=0.0;
    zp=0.0;
    f1=0.0;
    f2=0.0;
    f3=0.0;
    f4=0.0;
    g_virial=0.0;
    dgvirialdT=0.0;
    double fr=th.ed-T*th.en;
    vdet["zn"]=0.0;
    vdet["zp"]=0.0;
    vdet["F1"]=0.0;
    vdet["F2"]=0.0;
    vdet["F3"]=0.0;
    vdet["F4"]=0.0;
    vdet["g"]=0.0;
    vdet["dgdT"]=0.0;
    vdet["dgdnn"]=0.0;
    vdet["dgdnp"]=0.0;
    vdet["msn"]=neutron.ms*hc_mev_fm;
    vdet["msp"]=proton.ms*hc_mev_fm;
    vdet["Un"]=neutron.nu*hc_mev_fm;
    vdet["Up"]=proton.nu*hc_mev_fm;
    if (rmf_fields) {
      double sigma, omega, rho;
      rmf.get_fields(sigma,omega,rho);
      vdet["sigma"]=sigma;
      vdet["omega"]=omega;
      vdet["rho"]=rho;
    }
    return fr;
  }
  
  if (model_selected==false) {
    O2SCL_ERR("No model selected in free_energy_density().",
	      o2scl::exc_einval);
  }

  double nn=n.n;
  double pn=p.n;
  double nb=nn+pn;
  double ye=pn/nb;

  double n0=0.16;

  // ----------------------------------------------------------------
  // Compute the virial EOS

  double dmundT, dmupdT, dmundnn, dmupdnn, dmundnp, dmupdnp;
  if (T==0.0) {
    f_virial=0.0;
    s_virial=0.0;
    dmundT=0.0;
    dmupdT=0.0;
    dmundnn=0.0;
    dmupdnn=0.0;
    dmundnp=0.0;
    dmupdnp=0.0;
  } else {
    f_virial=free_energy_density_virial
      (n,p,T,th,dmundnn,dmundnp,dmupdnn,dmupdnp,dmundT,dmupdT);
    s_virial=th.en;
  }

  double dfvirialdT=-th.en;
  double P_virial=th.pr;
  double mu_n_virial=n.mu;
  double mu_p_virial=p.mu;
  zn=exp(mu_n_virial/T);
  zp=exp(mu_p_virial/T);
  g_virial=1.0/(a_virial*zn*zn+a_virial*zp*zp+b_virial*zn*zp+1.0);
  if (T==0.0) {
    dfvirialdT=0.0;
    P_virial=0.0;
    mu_n_virial=0.0;
    mu_p_virial=0.0;
    zn=0.0;
    zp=0.0;
    g_virial=0.0;
  }

  // Compute the T=0 effective masses by computing the asymmetric
  // Skyrme EOS at T=0:
  
  n.n=nn;
  p.n=pn;
  sk.calc_e(n,p,th);
  
  double msn_T0=n.ms, msp_T0=p.ms;
  
  // ----------------------------------------------------------------
  // Compute the Skyrme EOS in nuclear matter at T=0

  n.n=(nn+pn)/2.0;
  p.n=(nn+pn)/2.0;
  n.mu=n.m;
  p.mu=p.m;
  
  sk.calc_e(n,p,th);
  
  double mu_n_skyrme_eqdenT0=n.mu;
  double mu_p_skyrme_eqdenT0=p.mu;
  double P_skyrme_eqdenT0=-th.ed+mu_n_skyrme_eqdenT0*n.n+
    mu_p_skyrme_eqdenT0*p.n;
  double f_skyrme_eqdenT0=th.ed;

  double f_skyrme_T=0.0, f_skyrme_T0=0.0;
  double f_skyrme_neut_T, f_skyrme_neut_T0;
  double f_skyrme_eqden_T, f_skyrme_eqden_T0;
  double mu_n_skyrme_T=0.0, mu_n_skyrme_T0=0.0; 
  double mu_p_skyrme_T=0.0, mu_p_skyrme_T0=0.0; 
  double mu_n_neut_T, mu_n_neut_T0; 
  double mu_p_neut_T, mu_p_neut_T0; 
  double mu_n_eqden_T, mu_n_eqden_T0; 
  double mu_p_eqden_T, mu_p_eqden_T0; 
  double s_skyrme_T=0.0, s_eqden_T, s_neut_T;

  double msn_Tcorr_T0=0.0, msp_Tcorr_T0=0.0;
  double msn_Tcorr_T=0.0, msp_Tcorr_T=0.0;
  
  if (old_version==false) {
    
    // ----------------------------------------------------------------
    // Compute the Skyrme EOS in nuclear matter at T=0 for chiral fit
    
    n.n=nn;
    p.n=pn;
    n.mu=n.m;
    p.mu=p.m;
  
    eos_Tcorr->calc_e(n,p,th);
    
    mu_n_skyrme_T0=n.mu;
    mu_p_skyrme_T0=p.mu;
    double P_skyrme_T0=-th.ed+mu_n_skyrme_T0*n.n+mu_p_skyrme_T0*p.n;
    f_skyrme_T0=th.ed;
    msn_Tcorr_T0=n.ms;
    msp_Tcorr_T0=p.ms;
    
    // ----------------------------------------------------------------
    // Compute the Skyrme EOS in nuclear matter at finite T for chiral fit
    
    n.n=nn;
    p.n=pn;
    n.mu=n.m;
    p.mu=p.m;

    eos_Tcorr->calc_temp_e(n,p,T,th);
    
    mu_n_skyrme_T=n.mu;
    mu_p_skyrme_T=p.mu;
    double P_skyrme_T=-th.ed+mu_n_skyrme_T0*n.n+mu_p_skyrme_T0*p.n;
    f_skyrme_T=th.ed-T*th.en;
    s_skyrme_T=th.en;
    
    msn_Tcorr_T=n.ms;
    msp_Tcorr_T=p.ms;
    
    eos_Tcorr->calc_temp_e(n,p,T,th);
    
    // ----------------------------------------------------------------
    // Next, compute the Skyrme EOS at the specified density, proton
    // fraction, and temperature
    
    n_chiral.n=(nn+pn)/2.0;
    p_chiral.n=(nn+pn)/2.0;
    
    n_chiral.mu=n_chiral.m;
    p_chiral.mu=p_chiral.m;
    
    eos_Tcorr->calc_temp_e(n_chiral,p_chiral,T,th_chiral);
    
    f_skyrme_eqden_T=th_chiral.ed-T*th_chiral.en; 
    mu_p_eqden_T=p_chiral.mu;
    mu_n_eqden_T=n_chiral.mu;
    s_eqden_T=th_chiral.en;
    
    eos_Tcorr->calc_e(n_chiral,p_chiral,th_chiral);
    
    f_skyrme_eqden_T0=th_chiral.ed;
    mu_p_eqden_T0=p_chiral.mu;
    mu_n_eqden_T0=n_chiral.mu;
    
    n_chiral.n=nn+pn;
    p_chiral.n=0.0;
    p_chiral.n=1.0e-10;
    
    eos_Tcorr->calc_temp_e(n_chiral,p_chiral,T,th_chiral);
    
    f_skyrme_neut_T=th_chiral.ed-T*th_chiral.en; 
    mu_p_neut_T=p_chiral.mu;
    mu_n_neut_T=n_chiral.mu;
    s_neut_T=th_chiral.en;
    
    eos_Tcorr->calc_e(n_chiral,p_chiral,th_chiral);
    
    f_skyrme_neut_T0=th_chiral.ed;
    mu_p_neut_T0=p_chiral.mu;
    mu_n_neut_T0=n_chiral.mu;
    
  } else {

    // ----------------------------------------------------------------
    // Next, compute the Skyrme EOS at the specified density, proton
    // fraction, and temperature
    
    n_chiral.n=(nn+pn)/2.0;
    p_chiral.n=(nn+pn)/2.0;
    
    eos_Tcorr->calc_temp_e(n_chiral,p_chiral,T,th_chiral);
    
    f_skyrme_eqden_T=th_chiral.ed-T*th_chiral.en; 
    mu_p_eqden_T=p_chiral.mu;
    mu_n_eqden_T=n_chiral.mu;
    s_eqden_T=th_chiral.en;
    
    eos_Tcorr->calc_e(n_chiral,p_chiral,th_chiral);
    
    f_skyrme_eqden_T0=th_chiral.ed;
    mu_p_eqden_T0=p_chiral.mu;
    mu_n_eqden_T0=n_chiral.mu;
    
    n_chiral.n=nn+pn;
    p_chiral.n=0.0;
    
    eos_Tcorr->calc_temp_e(n_chiral,p_chiral,T,th_chiral);
    
    f_skyrme_neut_T=th_chiral.ed-T*th_chiral.en; 
    mu_p_neut_T=p_chiral.mu;
    mu_n_neut_T=n_chiral.mu;
    s_neut_T=th_chiral.en;
    
    eos_Tcorr->calc_e(n_chiral,p_chiral,th_chiral);
    
    f_skyrme_neut_T0=th_chiral.ed;
    mu_p_neut_T0=p_chiral.mu;
    mu_n_neut_T0=n_chiral.mu;

  }
  
  // ----------------------------------------------------------------
  // QMC EOS
  
  double e_qmc=energy_density_qmc(nn,pn);   

  // ----------------------------------------------------------------
  // Neutron star EOS
  
  double e_ns;
  double densdnn;

  new_ns_eos(nb,n,e_ns,densdnn);

  if (test_ns_cs2) {
    cout << "ns_nb_max: " << ns_nb_max << endl;
    cout << "phi: " << phi << endl;
    table_units<> tx;
    tx.line_of_names("nb ed mu");
    for(nb=0.08;nb<2.0;nb+=0.01) {
      new_ns_eos(nb,n,e_ns,densdnn);
      vector<double> line={nb,e_ns,densdnn};
      tx.line_of_data(3,line);
    }

    tx.function_column("ed+939.0/197.33*nb","edf");
    tx.function_column("mu+939.0/197.33","muf");
    tx.deriv("nb","ed","mu2");
    tx.deriv("nb","edf","muf2");
    tx.deriv("nb","muf2","dmufdn");
    tx.function_column("nb*dmufdn/muf","cs2");

    hdf_file hf;
    hf.open_or_create("ns2test.o2");
    hdf_output(hf,tx,"ns2test");
    hf.close();

    int sret=system("o2graph -read ns2test.o2 -plot nb cs2 -show");
  }
  
  // ----------------------------------------------------------------
  // Combine all the results to get the full free energy density
  // and put it in f_total

  double T_MeV=T*hc_mev_fm;
  double gamma=20.0;
  double h=1.0/(1.0+exp(gamma*(nn+pn-n0*1.5)));
  double e_combine=e_qmc*h+e_ns*(1.0-h);
  double e_sym=e_combine-f_skyrme_eqdenT0;
  double dyednn=-pn/nb/nb;
  double dyednp=nn/nb/nb;
  double delta2=(1.0-2.0*ye)*(1.0-2.0*ye);
  double ddelta2dnn=2.0*(1.0-2.0*ye)*(-2.0*dyednn);
  double ddelta2dnp=2.0*(1.0-2.0*ye)*(-2.0*dyednp);
  if (old_version==false) {
    f_deg=f_skyrme_eqdenT0+delta2*e_sym+
      f_skyrme_T-f_skyrme_T0;
  } else {
    f_deg=f_skyrme_eqdenT0+delta2*e_sym+
      delta2*(f_skyrme_neut_T-f_skyrme_neut_T0)+
      (1.0-delta2)*(f_skyrme_eqden_T-f_skyrme_eqden_T0);    
  }
  double f_total=f_virial*g_virial+f_deg*(1.0-g_virial);
  f1=f_virial*g_virial*(1.0-g_virial);
  f2=f_skyrme_eqdenT0*(1.0-g_virial);
  f3=delta2*e_sym*(1.0-g_virial);
  f4=(f_skyrme_T-f_skyrme_T0)*(1.0-g_virial);
  
  // -------------------------------------------------------------
  // Compute derivatives for chemical potentials
  
  dgvirialdnn=-(1.0/pow(a_virial*zn*zn+a_virial*zp*zp
			       +b_virial*zn*zp+1.0,2.0))*
    (2.0*a_virial*zn*zn/T*dmundnn
     +2.0*a_virial*zp*zp/T*dmupdnn+b_virial*zn*zp/T*dmundnn
     +b_virial*zn*zp/T*dmupdnn);
  dgvirialdnp=-(1.0/pow(a_virial*zn*zn+a_virial*zp*zp
			       +b_virial*zn*zp+1.0,2.0))*
    (2.0*a_virial*zn*zn/T*dmundnp
     +2.0*a_virial*zp*zp/T*dmupdnp+b_virial*zn*zp/T*dmundnp
     +b_virial*zn*zp/T*dmupdnp);
  
  double dfvirialdnn=mu_n_virial;
  double dfvirialdnp=mu_p_virial;
  if (T==0.0) {
    dgvirialdnn=0.0;
    dgvirialdnp=0.0;
    dfvirialdnn=0.0;
    dfvirialdnp=0.0;
  }
  
  double dfskyrme_eqden_T0dnn=(mu_n_skyrme_eqdenT0+mu_p_skyrme_eqdenT0)/2.0;
  double dfskyrme_eqden_T0dnp=dfskyrme_eqden_T0dnn;
	
  double dhdnn=-gamma*exp(gamma*(nn+pn-1.5*n0))/
    (1.0+exp(gamma*(nn+pn-1.5*n0)))/
    (1.0+exp(gamma*(nn+pn-1.5*n0)));
  
  double desymdnn=((qmc_a*pow((nn+pn)/qmc_n0,qmc_alpha)*(qmc_alpha+1.0)+
		    qmc_b*pow((nn+pn)/qmc_n0,qmc_beta)*(qmc_beta+1.0))/
		   hc_mev_fm)*h+e_qmc*dhdnn+
    densdnn*(1.0-h)-e_ns*dhdnn-dfskyrme_eqden_T0dnp/2.0-
    dfskyrme_eqden_T0dnn/2.0;
  double desymdnp=desymdnn;

  double dfdegdnn, dfdegdnp;
  if (old_version==false) {
    dfdegdnn=dfskyrme_eqden_T0dnn+(1.0-2.0*ye)
      *(1.0-2.0*ye)*desymdnn+ddelta2dnn*e_sym
      +mu_n_skyrme_T-mu_n_skyrme_T0;
    dfdegdnp=dfskyrme_eqden_T0dnp+(1.0-2.0*ye)
      *(1.0-2.0*ye)*desymdnp+ddelta2dnp*e_sym
      +mu_p_skyrme_T-mu_p_skyrme_T0;
  } else {
    dfdegdnn=dfskyrme_eqden_T0dnn+(1.0-2.0*ye)
      *(1.0-2.0*ye)*desymdnn+ddelta2dnn*e_sym
      +delta2*(mu_n_neut_T-mu_n_neut_T0)
      +ddelta2dnn*(f_skyrme_neut_T-f_skyrme_neut_T0)
      +(1.0-delta2)*(mu_n_eqden_T/2.0+mu_p_eqden_T/2.0-
		     mu_n_eqden_T0/2.0-mu_p_eqden_T0/2.0)-
      ddelta2dnn*(f_skyrme_eqden_T-f_skyrme_eqden_T0);
    dfdegdnp=dfskyrme_eqden_T0dnp+
      (1.0-2.0*ye)*(1.0-2.0*ye)*desymdnp+
      ddelta2dnp*e_sym
      +delta2*(mu_n_neut_T-mu_n_neut_T0)
      +ddelta2dnp*(f_skyrme_neut_T-f_skyrme_neut_T0)
      +(1.0-delta2)*(mu_p_eqden_T/2.0+mu_n_eqden_T/2.0-
		     mu_p_eqden_T0/2.0-mu_n_eqden_T0/2.0)-
      ddelta2dnp*(f_skyrme_eqden_T-f_skyrme_eqden_T0);    
  }

  // Final chemical potentials
  n.mu=dfvirialdnn*g_virial+f_virial*dgvirialdnn+dfdegdnn*(1.0-g_virial)
    +f_deg*(-dgvirialdnn);
  p.mu=dfvirialdnp*g_virial+f_virial*dgvirialdnp+dfdegdnp*(1.0-g_virial)
    +f_deg*(-dgvirialdnp);

  // Final effective masses
  n.ms=neutron.m*g_virial+(1.0-g_virial)*(msn_T0+msn_Tcorr_T-msn_Tcorr_T0);
  p.ms=proton.m*g_virial+(1.0-g_virial)*(msp_T0+msp_Tcorr_T-msp_Tcorr_T0);
  
  // -------------------------------------------------------------
  // Compute derivatives for entropy

  dgvirialdT=-(1.0/pow(1.0+a_virial*zn*zn+a_virial*zp*zp+
                       b_virial*zn*zp,2.0))*
    (2.0*a_virial*zn*zn*dmundT/T-2.0*a_virial*zn*zn*mu_n_virial/T/T
     +2.0*a_virial*zp*zp*dmupdT/T-2.0*a_virial*zp*zp*mu_p_virial/T/T
     +b_virial*zn*zp*dmundT/T+b_virial*zn*zp*dmupdT/T-
     b_virial*zn*zp*mu_n_virial/T/T-b_virial*zn*zp*mu_p_virial/T/T);
  if (T==0.0) dgvirialdT=0.0;
  
  // Restore p.n and n.n 
  n.n=nn;
  p.n=pn;
 
  double dfdegdT;
  if (old_version==false) {
    dfdegdT=-s_skyrme_T;
  } else {
    dfdegdT=delta2*(-s_neut_T)+(1.0-delta2)*(-s_eqden_T);
  }
  
  th.en=-(dfvirialdT*g_virial+f_virial*dgvirialdT+dfdegdT*(1-g_virial)
          +f_deg*(-dgvirialdT));
  th.pr=-f_total+n.n*n.mu+p.n*p.mu;
  th.ed=f_total+T*th.en;

  // Ensure the entropy is exactly zero at T=0 
  if (T==0.0) th.en=0.0;
  
  if (verbose>=1) {
    cout << endl;
    cout << "i_ns,i_skyrme= " << i_ns << " " << i_skyrme << endl;
    cout << endl;
    
    cout << "zn,zp= " << zn << " " << zp << endl;
    cout << "g_virial= " << g_virial << " (g=1 means full virial EOS) dgdT= "
	 << dgvirialdT << endl;
    cout << "h=        " << h
	 << " (h=1 means full QMC, h=0 means full NS)" << endl;
    cout << endl;

    cout.setf(ios::showpos);
    
    // Three contributions to the symmetry energy
    cout << "d2*e_qmc*h,d2*E_qmc*h     = " << delta2*e_qmc*h << " 1/fm^4 "
	 << delta2*e_qmc/nb*h*hc_mev_fm << " MeV" << endl;
    cout << "d2*e_ns(1-h),d2*E_ns(1-h) = "
	 << delta2*e_ns*(1.0-h) << " 1/fm^4 "
	 << delta2*e_ns*(1.0-h)/nb*hc_mev_fm << " MeV" << endl;
    cout << "d2*(-f_sk_nuc_T0)         = "
	 << -delta2*f_skyrme_eqdenT0 << " 1/fm^4 "
	 << -delta2*f_skyrme_eqdenT0/nb*hc_mev_fm << " MeV" << endl;
    cout << endl;

    // The four contributions to the free energy density of degenerate
    // matter
    cout << "f_nucT0, F_nucT0       = " << f_skyrme_eqdenT0 << " 1/fm^4 "
	 << f_skyrme_eqdenT0/nb*hc_mev_fm << " MeV" << endl;
    cout << "d2*f_symT0, d2*F_symT0 = " << delta2*e_sym << " 1/fm^4 "
	 << delta2*e_sym/nb*hc_mev_fm << " MeV" << endl;
    cout << "f_nucT, F_nucT         = "
	 << f_skyrme_eqden_T-f_skyrme_eqden_T0
	 << " 1/fm^4 "
	 << (f_skyrme_eqden_T-f_skyrme_eqden_T0)/nb*hc_mev_fm
	 << " MeV" << endl;
    cout << "d2*f_symT, d2*F_symT   = "
	 << delta2*(f_skyrme_neut_T-f_skyrme_neut_T0-
		    f_skyrme_eqden_T+f_skyrme_eqden_T0) << " 1/fm^4 "
	 << delta2*(f_skyrme_neut_T-f_skyrme_neut_T0-
		    f_skyrme_eqden_T+f_skyrme_eqden_T0)/nb*hc_mev_fm
	 << " MeV" << endl;
    cout << "f_deg, F_deg           = " << f_deg << " 1/fm^4 "
	 << f_deg/nb*hc_mev_fm << " MeV" << endl;
    cout << endl;

    // Virial and degenerate contributions and total free energy
    cout << "(1-g)*f_deg,(1-g)*F_deg = " << (1.0-g_virial)*f_deg << " 1/fm^4 "
	 << (1.0-g_virial)*f_deg/nb*hc_mev_fm << " MeV" << endl;
    cout << "g*f_virial,g*F_virial   = " << g_virial*f_virial << " 1/fm^4 "
	 << g_virial*f_virial/nb*hc_mev_fm << " MeV" << endl;
    cout << "f_total,F_total         = " << f_total << " 1/fm^4 "
	 << f_total/nb*hc_mev_fm << " MeV" << endl;
    cout << endl;

    cout << "ed (w/rm), pr= "
	 << th.ed+neutron.n*neutron.m+proton.n*proton.m
	 << " 1/fm^4 " << th.pr << " 1/fm^4" << endl;
    cout << "entropy, s per baryon= " << th.en << " 1/fm^3 "
	 << th.en/nb << endl;
    //cout << "entropy from sk_Tcorr= " << s_neut_T << " "
    //<< s_eqden_T << endl;
    //cout << "s_virial= " << s_virial << endl;
    //cout << "dg_virial_dT= " << dgvirialdT << endl;
    cout << "mu_n, mu_p = " << n.mu*hc_mev_fm << " MeV "
	 << p.mu*hc_mev_fm << " MeV" << endl;
    cout << endl;
    if (false) {
      cout << -dfvirialdT*g_virial-f_virial*dgvirialdT << " "
	   << -dfdegdT*(1-g_virial)+f_deg*(dgvirialdT) << endl;
      cout << -dfdegdT << " " << (1-g_virial) << endl;
      cout << f_skyrme_eqdenT0+delta2*e_sym << " "
	   << delta2 << " " << (f_skyrme_neut_T-f_skyrme_neut_T0) << " "
	   << (1.0-delta2) << " "
	   << (f_skyrme_eqden_T-f_skyrme_eqden_T0) << endl;
      cout << endl;
    }
    cout.unsetf(ios::showpos);
    
  }

  vdet["zn"]=zn;
  vdet["zp"]=zp;
  vdet["F1"]=f1/nb*hc_mev_fm;
  vdet["F2"]=f2/nb*hc_mev_fm;
  vdet["F3"]=f3/nb*hc_mev_fm;
  vdet["F4"]=f4/nb*hc_mev_fm;
  vdet["g"]=g_virial;
  vdet["dgdT"]=dgvirialdT/hc_mev_fm;
  vdet["dgdnn"]=dgvirialdnn;
  vdet["dgdnp"]=dgvirialdnp;
  vdet["msn"]=neutron.ms*hc_mev_fm;
  vdet["msp"]=proton.ms*hc_mev_fm;
  vdet["Un"]=neutron.nu*hc_mev_fm;
  vdet["Up"]=proton.nu*hc_mev_fm;
  
  return f_total;
}
  
double eos::free_energy_density_alt
(fermion &n, fermion &p, double nn, double np, double T, thermo &th) {
  n.n=nn;
  p.n=np;
  return free_energy_density(n,p,T,th);
}

double eos::free_energy_density_ep
(double nn, double np, double T) {
  neutron.n=nn;
  proton.n=np;
  
  elep.include_muons=include_muons;
  elep.pair_density_eq(proton.n,T);
  //cout << "three: " << elep.include_muons << " "
  //<< proton.n << " " << elep.e.n << endl;
  
  double frnp=free_energy_density(neutron,proton,T,th2);
  
  th2.ed+=elep.e.ed+elep.ph.ed;
  th2.pr+=elep.e.pr+elep.ph.pr;
  th2.en+=elep.e.en+elep.ph.en;
  
  return frnp+elep.e.ed-elep.e.en*T+elep.ph.ed-T*elep.ph.en;
}

double eos::entropy(fermion &n, fermion &p, double nn,
		    double pn, double T, thermo &th) {

  n.n=nn;
  p.n=pn;
  n.mu=n.m;
  p.mu=p.m;
  
  free_energy_density(n,p,T,th);
  
  elep.include_muons=include_muons;
  elep.pair_density_eq(p.n,T);
  
  return th.en+elep.e.en+elep.ph.en;
}

double eos::ed(fermion &n, fermion &p, double nn,
	       double pn, double T, thermo &th) {
  n.n=nn;
  p.n=pn;
  n.mu=n.m;
  p.mu=p.m;
  
  free_energy_density(n,p,T,th);
  
  elep.include_muons=include_muons;
  elep.pair_density_eq(p.n,T);

  return th.ed+elep.e.ed+photon.ed+n.m*nn+p.m*pn;
}

double eos::dfdnn_total(fermion &n, fermion &p, double nn, 
			double pn, double T, thermo &th) {
  
  n.n=nn;
  p.n=pn;
  n.mu=n.m;
  p.mu=p.m;
  
  free_energy_density(n,p,T,th);
  
  elep.include_muons=include_muons;
  elep.pair_density_eq(p.n,T);

  return n.mu+n.m;
}

double eos::dfdnp_total(fermion &n, fermion &p, double nn, 
			double pn, double T, thermo &th) {

  n.n=nn;
  p.n=pn;
  n.mu=n.m;
  p.mu=p.m;
  
  free_energy_density(n,p,T,th);
  
  elep.include_muons=include_muons;
  elep.pair_density_eq(p.n,T);
  
  return p.mu+elep.e.mu+p.m;
}

double eos::cs2_func(fermion &n, fermion &p, double T, thermo &th) {
  
  deriv_gsl<> gd;
  double nn=n.n;
  double np=p.n;
  double nb=n.n+p.n;

  double frx=free_energy_density_ep(nn,np,T);
  thermo thx;
  double enx=entropy(n,p,nn,np,T,thx);
  
  // Include the nucleon rest masses in the terms involving the
  // chemical potentials. We have to remember that the chemical
  // potentials are currently stored in neutron.mu and proton.mu
  // not in n.mu and p.mu (this is a bit confusing).
  double den=th2.en*T+(neutron.mu+n.m)*n.n+(proton.mu+p.m)*p.n+
    elep.e.mu*elep.e.n;
  double en=th2.en;
  if (cs2_verbose>0) {
    cout << "fr (hom): " << frx << endl;
    cout << "en (hom): " << enx << endl;
    cout << "mun,mup,mue (hom): " << neutron.mu << " " << proton.mu
         << " " << elep.e.mu << endl;
    cout << "den1,den2,den3,den4 (hom):\n  " << th2.en*T << " "
         << (n.mu+n.m)*n.n << " "
         << (p.mu+p.m)*p.n << " "
         << elep.e.mu << " " << elep.e.n << endl;
  }

  // Numerically compute required second derivatives
  double fac=1.0e3;

  // Note that the nucleon rest mass contributions are included in the
  // functions dfdnn_total and dfdnp_total.
  
  // d^2f/dnn^2
  std::function<double(double)> f_nnnn_func=
    std::bind(std::mem_fn<double(fermion &, fermion &, double,
				 double, double, thermo &)>
	      (&eos::dfdnn_total),
	      this,std::ref(n),std::ref(p),std::placeholders::_1,
	      np,T,std::ref(th));
  gd.h=fabs(nn)/fac;
  double f_nnnn=gd.deriv(nn,f_nnnn_func);
  
  // d^2f/dnn/dT
  std::function<double(double)> f_nnT_func=
    std::bind(std::mem_fn<double(fermion &, fermion &, double,
				 double, double, thermo &)>
	      (&eos::dfdnn_total),
	      this,std::ref(n),std::ref(p),nn,np,std::placeholders::_1,
	      std::ref(th));
  gd.h=fabs(T)/fac;
  double f_nnT=gd.deriv(T,f_nnT_func);

  // d^2f/dnp^2
  std::function<double(double)> f_npnp_func=
    std::bind(std::mem_fn<double(fermion &, fermion &, double,
				 double, double, thermo &)>
	      (&eos::dfdnp_total),
	      this,std::ref(n),std::ref(p),nn,std::placeholders::_1,
	      T,std::ref(th));
  gd.h=fabs(np)/fac;
  double f_npnp=gd.deriv(np,f_npnp_func);
  
  // d^2f/dnnnp
  std::function<double(double)> f_nnnp_func=
    std::bind(std::mem_fn<double(fermion &, fermion &, double,
				 double, double, thermo &)>
	      (&eos::dfdnn_total),
	      this,std::ref(n),std::ref(p),nn,std::placeholders::_1,
	      T,std::ref(th));
  gd.h=fabs(np)/fac;
  double f_nnnp=gd.deriv(np,f_nnnp_func);
  
  // d^2f/dnp/dT
  std::function<double(double)> f_npT_func=
    std::bind(std::mem_fn<double(fermion &, fermion &, double,
				 double, double, thermo &)>
	      (&eos::dfdnp_total),
	      this,std::ref(n),std::ref(p),nn,np,std::placeholders::_1,
	      std::ref(th));
  gd.h=fabs(T)/fac;
  double f_npT=gd.deriv(T,f_npT_func);

  // d^2f/dT^2
  std::function<double(double)> f_TT_func=
    std::bind(std::mem_fn<double(fermion &, fermion &, double,
				 double, double, thermo &)>
	      (&eos::entropy),
	      this,std::ref(n),std::ref(p),nn,np,std::placeholders::_1,
	      std::ref(th));
  gd.h=fabs(T)/fac;
  double f_TT=-gd.deriv(T,f_TT_func);

  if (cs2_verbose>0) {
    cout << "nn,np,en (hom): " << nn << " " << np << " " << en << endl;
    cout << "f_nnnn, f_nnnp, f_npnp, f_nnT, f_npT, f_TT, den (hom):\n  "
         << f_nnnn << " " << f_nnnp << " " << f_npnp << " "
         << f_nnT << "\n  " << f_npT << " " << f_TT << " "
         << den << endl;
  }
  
  double cs_sq=(nn*nn*(f_nnnn-f_nnT*f_nnT/f_TT)+
		2.0*nn*np*(f_nnnp-f_nnT*f_npT/f_TT)+
		np*np*(f_npnp-f_npT*f_npT/f_TT)-
		2.0*en*(nn*f_nnT/f_TT+np*f_npT/f_TT)-en*en/f_TT)/den;
  
  return cs_sq;
}

int eos::table_Ye(std::vector<std::string> &sv, bool itive_com) {

  std::string fname=sv[1];
  double Ye=o2scl::stod(sv[2]);

  size_t n_nB=301;
  size_t n_T=160;
    
  std::string nB_grid_spec2="10^(i*0.04-12)*2.0";
  std::string T_grid_spec2="0.2+0.81*i";

  vector<double> nB_grid, T_grid;
  
  calc_utf8<> calc;
  std::map<std::string,double> vars;
  
  calc.compile(nB_grid_spec2.c_str());
  for(size_t i=0;i<n_nB;i++) {
    vars["i"]=((double)i);
    nB_grid.push_back(calc.eval(&vars));
  }
  
  calc.compile(T_grid_spec2.c_str());
  for(size_t i=0;i<n_T;i++) {
    vars["i"]=((double)i);
    T_grid.push_back(calc.eval(&vars));
  }

  table3d t;
  t.set_xy("nB",n_nB,nB_grid,"T",n_T,T_grid);
  t.new_slice("Fint");
  t.new_slice("Pint");
  t.new_slice("Sint");
  t.new_slice("Eint");
  t.new_slice("F");
  t.new_slice("P");
  t.new_slice("S");
  t.new_slice("E");
  t.new_slice("g");
  t.new_slice("msn");
  t.new_slice("msp");
  t.new_slice("mun");
  t.new_slice("mup");
  t.new_slice("mue");
  t.new_slice("cs2");

  elep.e.mu=elep.e.m;
  
  for(int i=n_nB-1;i>=0;i--) {
    cout << i << "/" << n_nB << endl;
    for(size_t j=0;j<n_T;j++) {
      
      neutron.n=nB_grid[i]*(1.0-Ye);
      proton.n=nB_grid[i]*Ye;
      
      free_energy_density(neutron,proton,T_grid[j]/hc_mev_fm,th2);
      double foa_hc=(hc_mev_fm*th2.ed-T_grid[j]*th2.en)/nB_grid[i];
      
      elep.include_muons=include_muons;
      elep.e.n=nB_grid[i]*Ye/2.0;
      elep.pair_density_eq(nB_grid[i]*Ye,T_grid[j]/hc_mev_fm);

      t.set(i,j,"Fint",foa_hc);
      t.set(i,j,"Sint",th2.en/nB_grid[i]);
      t.set(i,j,"Eint",th2.ed/nB_grid[i]*hc_mev_fm);
      t.set(i,j,"Pint",th2.pr*hc_mev_fm);
      t.set(i,j,"g",0.0);
      
      t.set(i,j,"msn",neutron.ms/neutron.m);
      t.set(i,j,"msp",proton.ms/proton.m);

      t.set(i,j,"mun",neutron.mu*hc_mev_fm);
      t.set(i,j,"mup",proton.mu*hc_mev_fm);
      t.set(i,j,"mue",elep.e.mu*hc_mev_fm);

      t.set(i,j,"F",foa_hc+(elep.th.ed*hc_mev_fm-
                            elep.th.en*T_grid[j])/nB_grid[i]);
      t.set(i,j,"E",th2.ed/nB_grid[i]*hc_mev_fm+
            (elep.th.ed*hc_mev_fm)/nB_grid[i]);
      t.set(i,j,"P",th2.pr+elep.th.pr*hc_mev_fm);
      t.set(i,j,"S",th2.en/nB_grid[i]+elep.th.en/nB_grid[i]);

      double cs2_val=cs2_func(neutron,proton,T_grid[j]/hc_mev_fm,th2);
      t.set(i,j,"cs2",cs2_val);
      
    }
  }

  hdf_file hf;
  hf.open_or_create(fname);
  hdf_output(hf,(const table3d &)t,"table_Ye");
  hf.close();
  
  return 0;
}

int eos::table_nB(std::vector<std::string> &sv, bool itive_com) {

  std::string fname=sv[1];
  double nB=o2scl::stod(sv[2]);

  size_t n_Ye=99;
  size_t n_T=160;
    
  std::string Ye_grid_spec2="0.01*(i+1)";
  std::string T_grid_spec2="0.2+0.81*i";

  vector<double> Ye_grid, T_grid;
  
  calc_utf8<> calc;
  std::map<std::string,double> vars;
  
  calc.compile(Ye_grid_spec2.c_str());
  for(size_t i=0;i<n_Ye;i++) {
    vars["i"]=((double)i);
    Ye_grid.push_back(calc.eval(&vars));
  }
  
  calc.compile(T_grid_spec2.c_str());
  for(size_t i=0;i<n_T;i++) {
    vars["i"]=((double)i);
    T_grid.push_back(calc.eval(&vars));
  }

  table3d t;
  t.set_xy("Ye",n_Ye,Ye_grid,"T",n_T,T_grid);
  t.line_of_names("Fint Sint Pint g msn msp cs2");

  for(size_t i=0;i<n_Ye;i++) {
    cout << i << "/" << n_Ye << endl;
    for(size_t j=0;j<n_T;j++) {
      neutron.n=nB*(1.0-Ye_grid[i]);
      proton.n=nB*Ye_grid[i];
      double t1, t2, t3, t4, t5;
      free_energy_density(neutron,proton,T_grid[j]/hc_mev_fm,th2);
      double foa_hc=hc_mev_fm*(th2.ed-T_grid[j]/hc_mev_fm*th2.en)/
	(neutron.n+proton.n);
      t.set(i,j,"Fint",foa_hc);
      t.set(i,j,"Sint",th2.en/nB);
      t.set(i,j,"Pint",th2.pr*hc_mev_fm);
      t.set(i,j,"g",0.0);
      t.set(i,j,"msn",neutron.ms/neutron.m);
      t.set(i,j,"msp",proton.ms/proton.m);
      double cs2_val=cs2_func(neutron,proton,T_grid[j]/hc_mev_fm,th2);
      t.set(i,j,"cs2",cs2_val);
    }
  }

  hdf_file hf;
  hf.open_or_create(fname);
  hdf_output(hf,(const table3d &)t,"table_nB");
  hf.close();
  
  return 0;
}

int eos::process_grid_spec() {
  
  size_t st[3]={n_nB2,n_Ye2,n_T2};
  cout << "Processing: " << nB_grid_spec << " " << Ye_grid_spec << " "
       << T_grid_spec << endl;
  
  calc_utf8<> calc;
  std::map<std::string,double> vars;
  
  vector<std::string> split_res;

  split_string_delim(nB_grid_spec,split_res,',');
  n_nB2=stoszt(split_res[0]);
  nB_grid2.resize(0);
  cout << "n_nB: " << n_nB2 << endl;
  
  calc.compile(split_res[1].c_str());
  for(size_t i=0;i<n_nB2;i++) {
    vars["i"]=((double)i);
    nB_grid2.push_back(calc.eval(&vars));
  }
  
  split_string_delim(Ye_grid_spec,split_res,',');
  n_Ye2=stoszt(split_res[0]);
  Ye_grid2.resize(0);
  cout << "n_Ye: " << n_Ye2 << endl;
  
  calc.compile(split_res[1].c_str());
  for(size_t i=0;i<n_Ye2;i++) {
    vars["i"]=((double)i);
    Ye_grid2.push_back(calc.eval(&vars));
  }
  
  split_string_delim(T_grid_spec,split_res,',');
  n_T2=stoszt(split_res[0]);
  T_grid2.resize(0);
  cout << "n_T: " << n_T2 << endl;
  
  calc.compile(split_res[1].c_str());
  for(size_t i=0;i<n_T2;i++) {
    vars["i"]=((double)i);
    T_grid2.push_back(calc.eval(&vars));
  }
  
  split_string_delim(S_grid_spec,split_res,',');
  n_S2=stoszt(split_res[0]);
  S_grid2.resize(0);
  cout << "n_S: " << n_S2 << endl;
  
  calc.compile(split_res[1].c_str());
  for(size_t i=0;i<n_S2;i++) {
    vars["i"]=((double)i);
    S_grid2.push_back(calc.eval(&vars));
  }
  
  return 0;
}

int eos::table_full(std::vector<std::string> &sv, bool itive_com) {

  std::string fname=sv[1];

  process_grid_spec();
    
  size_t size_arr[3]={n_nB2,n_Ye2,n_T2};
  vector<vector<double> > grid_arr={nB_grid2,Ye_grid2,T_grid2};

  tensor_grid3<> t_F(n_nB2,n_Ye2,n_T2);
  t_F.set_grid(grid_arr);
  tensor_grid3<> t_Fint(n_nB2,n_Ye2,n_T2);
  t_Fint.set_grid(grid_arr);
  tensor_grid3<> t_E(n_nB2,n_Ye2,n_T2);
  t_E.set_grid(grid_arr);
  tensor_grid3<> t_Eint(n_nB2,n_Ye2,n_T2);
  t_Eint.set_grid(grid_arr);
  tensor_grid3<> t_P(n_nB2,n_Ye2,n_T2);
  t_P.set_grid(grid_arr);
  tensor_grid3<> t_Pint(n_nB2,n_Ye2,n_T2);
  t_Pint.set_grid(grid_arr);
  tensor_grid3<> t_S(n_nB2,n_Ye2,n_T2);
  t_S.set_grid(grid_arr);
  tensor_grid3<> t_Sint(n_nB2,n_Ye2,n_T2);
  t_Sint.set_grid(grid_arr);
  tensor_grid3<> t_mun(n_nB2,n_Ye2,n_T2);
  t_mun.set_grid(grid_arr);
  tensor_grid3<> t_mup(n_nB2,n_Ye2,n_T2);
  t_mup.set_grid(grid_arr);
  tensor_grid3<> t_cs2(n_nB2,n_Ye2,n_T2);
  t_cs2.set_grid(grid_arr);
  tensor_grid3<> t_mue(n_nB2,n_Ye2,n_T2);
  t_mue.set_grid(grid_arr);

  for(int i=n_nB2-1;i>=0;i--) {
    cout << "i_nB,n_nB,nB[i]: " << n_nB2-1-i << " " << n_nB2 << " "
	 << nB_grid2[i] << endl;
    for(size_t j=0;j<n_Ye2;j++) {
      for(size_t k=0;k<n_T2;k++) {

	// Hadronic part
	neutron.n=nB_grid2[i]*(1.0-Ye_grid2[j]);
	proton.n=nB_grid2[i]*Ye_grid2[j];
        free_energy_density(neutron,proton,T_grid2[k]/hc_mev_fm,th2);

        elep.include_muons=include_muons;
        elep.include_photons=true;
        
        elep.pair_density_eq(nB_grid2[i]*Ye_grid2[j],T_grid2[k]/hc_mev_fm);

	t_Fint.set(i,j,k,hc_mev_fm*(th2.ed-T_grid2[k]/hc_mev_fm*th2.en)/
		   (neutron.n+proton.n));
	t_F.set(i,j,k,hc_mev_fm*(th2.ed+elep.th.ed-T_grid2[k]/hc_mev_fm*th2.en)/
		(neutron.n+proton.n));
	t_Eint.set(i,j,k,hc_mev_fm*(th2.ed)/(neutron.n+proton.n));
	t_E.set(i,j,k,hc_mev_fm*(th2.ed+elep.th.ed)/(neutron.n+proton.n));
	t_Pint.set(i,j,k,hc_mev_fm*(th2.pr));
	t_P.set(i,j,k,hc_mev_fm*(th2.pr+elep.th.pr));
	t_Sint.set(i,j,k,hc_mev_fm*(th2.en)/(neutron.n+proton.n));
	t_S.set(i,j,k,hc_mev_fm*(th2.en+elep.th.en)/(neutron.n+proton.n));
	t_mun.set(i,j,k,hc_mev_fm*neutron.mu);
	t_mup.set(i,j,k,hc_mev_fm*proton.mu);

	double cs2_val=cs2_func(neutron,proton,T_grid2[k]/hc_mev_fm,th2);
	t_cs2.set(i,j,k,cs2_val);
	t_mue.set(i,j,k,elep.e.mu);

	if (!std::isfinite(th2.ed)) {
	  cout << "Hadronic energy density not finite." << endl;
	  cout << "n_B: " << nB_grid2[i] << " Y_e: " << Ye_grid2[j] << " T: "
	       << T_grid2[k] << endl;
	  cout << "hadrons ed: " << th2.ed << " pr: " << th2.pr
	       << " en: " << th2.en << endl;
	  exit(-1);
	}
	if (!std::isfinite(th2.pr)) {
	  cout << "Hadronic pressure not finite." << endl;
	  cout << "n_B: " << nB_grid2[i] << " Y_e: " << Ye_grid2[j] << " T: "
	       << T_grid2[k] << endl;
	  cout << "hadrons ed: " << th2.ed << " pr: " << th2.pr
	       << " en: " << th2.en << endl;
	  exit(-1);
	}
	if (!std::isfinite(th2.en)) {
	  cout << "Hadronic entropy density not finite." << endl;
	  cout << "n_B: " << nB_grid2[i] << " Y_e: " << Ye_grid2[j] << " T: "
	       << T_grid2[k] << endl;
	  cout << "hadrons ed: " << th2.ed << " pr: " << th2.pr
	       << " en: " << th2.en << endl;
	  exit(-1);
	}
	if (!std::isfinite(elep.th.ed)) {
	  cout << "Leptonic energy density not finite." << endl;
	  cout << "n_B: " << nB_grid2[i] << " Y_e: " << Ye_grid2[j] << " T: "
	       << T_grid2[k] << endl;
	  cout << "leptons ed: " << elep.th.ed << " pr: " << elep.th.pr
	       << " en: " << elep.th.en << endl;
	  exit(-1);
	}
	if (!std::isfinite(elep.th.pr)) {
	  cout << "Leptonic pressure not finite." << endl;
	  cout << "n_B: " << nB_grid2[i] << " Y_e: " << Ye_grid2[j] << " T: "
	       << T_grid2[k] << endl;
	  cout << "leptons ed: " << elep.th.ed << " pr: " << elep.th.pr
	       << " en: " << elep.th.en << endl;
	  exit(-1);
	}
	if (!std::isfinite(elep.th.en)) {
	  cout << "Leptonic entropy density not finite." << endl;
	  cout << "n_B: " << nB_grid2[i] << " Y_e: " << Ye_grid2[j] << " T: "
	       << T_grid2[k] << endl;
	  cout << "leptons ed: " << elep.th.ed << " pr: " << elep.th.pr
	       << " en: " << elep.th.en << endl;
	  exit(-1);
	}
	if (th2.en+elep.th.en<0.0 && th2.pr>0.0) {
	  cout << "Entropy negative where pressure is positive." << endl;
	  cout << "n_B: " << nB_grid2[i] << " Y_e: " << Ye_grid2[j] << " T: "
	       << T_grid2[k] << endl;
	  cout << "hadrons ed: " << th2.ed << " pr: " << th2.pr
	       << " en: " << th2.en << endl;
	  cout << "leptons ed: " << elep.th.ed << " pr: " << elep.th.pr
	       << " en: " << elep.th.en << endl;
	  exit(-1);
	}

      }
    }
  }

  hdf_file hf;
  hf.open_or_create(fname);
  hf.set_szt("n_nB",n_nB2);
  hf.set_szt("n_Ye",n_Ye2);
  hf.set_szt("n_T",n_T2);
  hf.setd_vec("nB_grid",nB_grid2);
  hf.setd_vec("Ye_grid",Ye_grid2);
  hf.setd_vec("T_grid",T_grid2);
  hdf_output(hf,t_Fint,"Fint");
  hdf_output(hf,t_F,"F");
  hdf_output(hf,t_Eint,"Eint");
  hdf_output(hf,t_E,"E");
  hdf_output(hf,t_Pint,"Pint");
  hdf_output(hf,t_P,"P");
  hdf_output(hf,t_Sint,"Sint");
  hdf_output(hf,t_S,"S");
  hdf_output(hf,t_mun,"mun");
  hdf_output(hf,t_mup,"mup");
  hdf_output(hf,t_cs2,"cs2");
  hdf_output(hf,t_mue,"mue");
  hf.seti("with_leptons",1);
  hf.seti("baryons_only",1);
  hf.seti("derivs_computed",1);
  
  vector<string> oth_names={"cs2"};
  vector<string> oth_units={""};
  size_t n_oth=oth_names.size();
  hf.set_szt("n_oth",n_oth);
  hf.sets_vec_copy("oth_names",oth_names);
  hf.sets_vec_copy("oth_units",oth_units);
  
  hf.close();
  
  return 0;
}

int eos::test_deriv(std::vector<std::string> &sv, bool itive_com) {

  if (false) {
    // Make the Skyrme EOSs more accurate
    sk.nrf.def_density_root.tol_rel/=1.0e2;
    sk.nrf.def_density_root.tol_abs/=1.0e2;
    sk_Tcorr.nrf.def_density_root.tol_rel/=1.0e2;
    sk_Tcorr.nrf.def_density_root.tol_abs/=1.0e2;
  }

  cout.precision(5);
  
  test_mgr t;
  t.set_output_level(1);

  if (model_selected==false && use_alt_eos==false) {
    cerr << "No model selected." << endl;
    return 1;
  }

  double err_fac;
    
  deriv_gsl<> gd;
  double mun_num, mup_num, en_num=0.0, mun_err, mup_err, en_err=0.0;
  double T, nn, np, nb, avg1, avg2, avg3, avg_sum=0.0;
  size_t count;

  // ---------------------------------------------------------------------

  cout << "T=0.1 MeV, Ye=0.01" << endl;
  cout << "n_B         mu_n        mup_p       s" << endl;
  
  err_fac=1.0e6;
  avg1=0.0;
  avg2=0.0;
  avg3=0.0;
  for(count=0,nb=1.0e-10;nb<1.6;nb*=1.3,count++) {

    // Neutron chemical potential
    neutron.n=nb*0.99;
    proton.n=nb*0.01;
    T=0.1/hc_mev_fm;
    std::function<double(double)> fef_nn=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),
		std::placeholders::_1,proton.n,T,
		std::ref(th2));
    gd.h=neutron.n/1.0e2;
    gd.deriv_err(neutron.n,fef_nn,mun_num,mun_err);
    
    // Proton chemical potential
    neutron.n=nb*0.99;
    proton.n=nb*0.01;
    T=0.1/hc_mev_fm;
    std::function<double(double)> fef_np=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),neutron.n,
		std::placeholders::_1,T,std::ref(th2));
    gd.h=proton.n/1.0e2;
    gd.deriv_err(proton.n,fef_np,mup_num,mup_err);
    
    // Entropy
    neutron.n=nb*0.99;
    proton.n=nb*0.01;
    T=0.1/hc_mev_fm;
    std::function<double(double)> fef_en=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),neutron.n,proton.n,
		std::placeholders::_1,std::ref(th2));
    gd.h=T/1.0e2;
    gd.deriv_err(T,fef_en,en_num,en_err);
    en_num*=-1.0;

    // Now compute analytical results from function
    neutron.n=nb*0.99;
    proton.n=nb*0.01;
    T=0.1/hc_mev_fm;
    free_energy_density(neutron,proton,T,th2);

    double rat1=fabs((mun_num-neutron.mu)/mun_err);
    double rat2=fabs((mup_num-proton.mu)/mup_err);
    double rat3=fabs((en_num-th2.en)/en_err);
    if (count%10==0 || rat1>10.0 || rat2>10.0 || rat3>10.0) {
      cout << nb << " ";
      cout << rat1 << " ";
      cout << rat2 << " ";
      cout << rat3 << endl;
    }
    avg1+=rat1;
    avg2+=rat2;
    avg3+=rat3;
    avg_sum+=avg1+avg2+avg3;
    
    t.test_abs(neutron.mu,mun_num,mun_err*err_fac,"mun, T=0.1, Ye=0.01");
    t.test_abs(proton.mu,mup_num,mup_err*err_fac,"mup, T=0.1, Ye=0.01");
    t.test_abs(th2.en,en_num,en_err*err_fac,"en, T=0.1, Ye=0.01");
  }
  cout << "  averages: " << avg1/((double)count) << " ";
  cout << avg2/((double)count) << " ";
  cout << avg3/((double)count) << endl;
  cout << endl;

  // ---------------------------------------------------------------------
  
  cout << "T=0.1 MeV, Ye=0.49" << endl;
  cout << "n_B         mu_n        mup_p       s" << endl;

  err_fac=1.0e6;
  avg1=0.0;
  avg2=0.0;
  avg3=0.0;
  for(count=0,nb=1.0e-10;nb<1.6;nb*=1.3,count++) {
    neutron.n=nb*0.51;
    proton.n=nb*0.49;
    T=0.1/hc_mev_fm;
    
    // Neutron chemical potential
    std::function<double(double)> fef_nn=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),
		std::placeholders::_1,proton.n,T,
		std::ref(th2));
    gd.h=neutron.n/1.0e2;
    gd.deriv_err(neutron.n,fef_nn,mun_num,mun_err);
    
    // Proton chemical potential
    neutron.n=nb*0.51;
    proton.n=nb*0.49;
    T=0.1/hc_mev_fm;
    std::function<double(double)> fef_np=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),neutron.n,
		std::placeholders::_1,T,std::ref(th2));
    gd.h=proton.n/1.0e2;
    gd.deriv_err(proton.n,fef_np,mup_num,mup_err);
    
    // Entropy
    neutron.n=nb*0.51;
    proton.n=nb*0.49;
    T=0.1/hc_mev_fm;
    std::function<double(double)> fef_en=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),neutron.n,proton.n,
		std::placeholders::_1,std::ref(th2));
    gd.h=T/1.0e2;
    gd.deriv_err(T,fef_en,en_num,en_err);
    en_num*=-1.0;
    
    // Now compute analytical results from function
    neutron.n=nb*0.51;
    proton.n=nb*0.49;
    T=0.1/hc_mev_fm;
    free_energy_density(neutron,proton,T,th2);
    
    double rat1=fabs((mun_num-neutron.mu)/mun_err);
    double rat2=fabs((mup_num-proton.mu)/mup_err);
    double rat3=fabs((en_num-th2.en)/en_err);
    if (count%10==0 || rat1>10.0 || rat2>10.0 || rat3>10.0) {
      cout << nb << " ";
      cout << rat1 << " ";
      cout << rat2 << " ";
      cout << rat3 << endl;
    }
    avg1+=rat1;
    avg2+=rat2;
    avg3+=rat3;
    avg_sum+=avg1+avg2+avg3;
    
    t.test_abs(neutron.mu,mun_num,mun_err*err_fac,"mun, T=0.1, Ye=0.49");
    t.test_abs(proton.mu,mup_num,mup_err*err_fac,"mup, T=0.1, Ye=0.49");
    t.test_abs(th2.en,en_num,en_err*err_fac,"en, T=0.1, Ye=0.49");
  }
  cout << "  averages: " << avg1/((double)count) << " ";
  cout << avg2/((double)count) << " ";
  cout << avg3/((double)count) << endl;
  cout << endl;
  
  // ---------------------------------------------------------------------

  cout << "T=1 MeV, Ye=0.01" << endl;
  cout << "n_B         mu_n        mup_p       s" << endl;

  err_fac=1.0e6;
  avg1=0.0;
  avg2=0.0;
  avg3=0.0;
  for(count=0,nb=1.0e-10;nb<1.6;nb*=1.3,count++) {

    // Neutron chemical potential
    neutron.n=nb*0.99;
    proton.n=nb*0.01;
    T=1.0/hc_mev_fm;
    std::function<double(double)> fef_nn=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),
		std::placeholders::_1,proton.n,T,
		std::ref(th2));
    gd.h=neutron.n/1.0e2;
    gd.deriv_err(neutron.n,fef_nn,mun_num,mun_err);
          
    // Proton chemical potential
    neutron.n=nb*0.99;
    proton.n=nb*0.01;
    T=1.0/hc_mev_fm;
    std::function<double(double)> fef_np=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),neutron.n,
		std::placeholders::_1,T,std::ref(th2));
    gd.h=proton.n/1.0e2;
    gd.deriv_err(proton.n,fef_np,mup_num,mup_err);
          
    // Entropy
    neutron.n=nb*0.99;
    proton.n=nb*0.01;
    T=1.0/hc_mev_fm;
    std::function<double(double)> fef_en=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),neutron.n,proton.n,
		std::placeholders::_1,std::ref(th2));
    gd.h=T/1.0e2;
    gd.deriv_err(T,fef_en,en_num,en_err);
    en_num*=-1.0;

    // Now compute analytical results from function
    neutron.n=nb*0.99;
    proton.n=nb*0.01;
    T=1.0/hc_mev_fm;
    free_energy_density(neutron,proton,T,th2);
          
    double rat1=fabs((mun_num-neutron.mu)/mun_err);
    double rat2=fabs((mup_num-proton.mu)/mup_err);
    double rat3=fabs((en_num-th2.en)/en_err);
    if (count%10==0 || rat1>10.0 || rat2>10.0 || rat3>10.0) {
      cout << nb << " ";
      cout << rat1 << " ";
      cout << rat2 << " ";
      cout << rat3 << endl;
    }
    avg1+=rat1;
    avg2+=rat2;
    avg3+=rat3;
    avg_sum+=avg1+avg2+avg3;

    t.test_abs(neutron.mu,mun_num,mun_err*err_fac,"mun, T=1, Ye=0.01");
    t.test_abs(proton.mu,mup_num,mup_err*err_fac,"mup, T=1, Ye=0.01");
    t.test_abs(th2.en,en_num,en_err*err_fac,"en, T=1, Ye=0.01");
  }
  cout << "  averages: " << avg1/((double)count) << " ";
  cout << avg2/((double)count) << " ";
  cout << avg3/((double)count) << endl;
  cout << endl;

  // ---------------------------------------------------------------------

  cout << "T=1 MeV, Ye=0.49" << endl;
  cout << "n_B         mu_n        mup_p       s" << endl;

  err_fac=1.0e6;
  avg1=0.0;
  avg2=0.0;
  avg3=0.0;
  for(count=0,nb=1.0e-10;nb<1.6;nb*=1.3,count++) {

    // Neutron chemical potential
    neutron.n=nb*0.51;
    proton.n=nb*0.49;
    T=1.0/hc_mev_fm;
    std::function<double(double)> fef_nn=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),
		std::placeholders::_1,proton.n,T,
		std::ref(th2));
    gd.h=neutron.n/1.0e2;
    gd.deriv_err(neutron.n,fef_nn,mun_num,mun_err);
          
    // Proton chemical potential
    neutron.n=nb*0.51;
    proton.n=nb*0.49;
    T=1.0/hc_mev_fm;
    std::function<double(double)> fef_np=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),neutron.n,
		std::placeholders::_1,T,std::ref(th2));
    gd.h=proton.n/1.0e2;
    gd.deriv_err(proton.n,fef_np,mup_num,mup_err);
          
    // Entropy
    neutron.n=nb*0.51;
    proton.n=nb*0.49;
    T=1.0/hc_mev_fm;
    std::function<double(double)> fef_en=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),neutron.n,proton.n,
		std::placeholders::_1,std::ref(th2));
    gd.h=T/1.0e2;
    gd.deriv_err(T,fef_en,en_num,en_err);
    en_num*=-1.0;

    // Now compute analytical results from function
    neutron.n=nb*0.51;
    proton.n=nb*0.49;
    T=1.0/hc_mev_fm;
    free_energy_density(neutron,proton,T,th2);
          
    double rat1=fabs((mun_num-neutron.mu)/mun_err);
    double rat2=fabs((mup_num-proton.mu)/mup_err);
    double rat3=fabs((en_num-th2.en)/en_err);
    if (count%10==0 || rat1>10.0 || rat2>10.0 || rat3>10.0) {
      cout << nb << " ";
      cout << rat1 << " ";
      cout << rat2 << " ";
      cout << rat3 << endl;
    }
    avg1+=rat1;
    avg2+=rat2;
    avg3+=rat3;
    avg_sum+=avg1+avg2+avg3;

    t.test_abs(neutron.mu,mun_num,mun_err*err_fac,"mun, T=1, Ye=0.49");
    t.test_abs(proton.mu,mup_num,mup_err*err_fac,"mup, T=1, Ye=0.49");
    t.test_abs(th2.en,en_num,en_err*err_fac,"en, T=1, Ye=0.49");
  }
  cout << "  averages: " << avg1/((double)count) << " ";
  cout << avg2/((double)count) << " ";
  cout << avg3/((double)count) << endl;
  cout << endl;

  // ---------------------------------------------------------------------

  cout << "T=30 MeV, Ye=0.01" << endl;
  cout << "n_B         mu_n        mup_p       s" << endl;

  err_fac=1.0e6;
  avg1=0.0;
  avg2=0.0;
  avg3=0.0;
  for(count=0,nb=1.0e-10;nb<1.6;nb*=1.3,count++) {
          
    // Neutron chemical potential
    neutron.n=nb*0.99;
    proton.n=nb*0.01;
    T=30.0/hc_mev_fm;
    std::function<double(double)> fef_nn=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),
		std::placeholders::_1,proton.n,T,
		std::ref(th2));
    gd.h=neutron.n/1.0e2;
    gd.deriv_err(neutron.n,fef_nn,mun_num,mun_err);
          
    // Proton chemical potential
    neutron.n=nb*0.99;
    proton.n=nb*0.01;
    T=30.0/hc_mev_fm;
    std::function<double(double)> fef_np=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),neutron.n,
		std::placeholders::_1,T,std::ref(th2));
    gd.h=proton.n/1.0e2;
    gd.deriv_err(proton.n,fef_np,mup_num,mup_err);
          
    // Entropy
    neutron.n=nb*0.99;
    proton.n=nb*0.01;
    T=30.0/hc_mev_fm;
    std::function<double(double)> fef_en=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),neutron.n,proton.n,
		std::placeholders::_1,std::ref(th2));
    gd.h=T/1.0e2;
    gd.deriv_err(T,fef_en,en_num,en_err);
    en_num*=-1.0;

    // Now compute analytical results from function
    neutron.n=nb*0.99;
    proton.n=nb*0.01;
    T=30.0/hc_mev_fm;
    free_energy_density(neutron,proton,T,th2);
          
    double rat1=fabs((mun_num-neutron.mu)/mun_err);
    double rat2=fabs((mup_num-proton.mu)/mup_err);
    double rat3=fabs((en_num-th2.en)/en_err);
    if (count%10==0 || rat1>10.0 || rat2>10.0 || rat3>10.0) {
      cout << nb << " ";
      cout << rat1 << " ";
      cout << rat2 << " ";
      cout << rat3 << endl;
    }
    avg1+=rat1;
    avg2+=rat2;
    avg3+=rat3;
    avg_sum+=avg1+avg2+avg3;

    t.test_abs(neutron.mu,mun_num,mun_err*err_fac,"mun, T=30, Ye=0.01");
    t.test_abs(proton.mu,mup_num,mup_err*err_fac,"mup, T=30, Ye=0.01");
    t.test_abs(th2.en,en_num,en_err*err_fac,"en, T=30, Ye=0.01");
  }
  cout << "  averages: " << avg1/((double)count) << " ";
  cout << avg2/((double)count) << " ";
  cout << avg3/((double)count) << endl;
  cout << endl;
  
  // ---------------------------------------------------------------------

  cout << "T=30 MeV, Ye=0.49" << endl;
  cout << "n_B         mu_n        mup_p       s" << endl;
  
  err_fac=1.0e6;
  avg1=0.0;
  avg2=0.0;
  avg3=0.0;
  for(count=0,nb=1.0e-10;nb<1.6;nb*=1.3,count++) {

    // Neutron chemical potential
    neutron.n=nb*0.51;
    proton.n=nb*0.49;
    T=30.0/hc_mev_fm;
    std::function<double(double)> fef_nn=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),
		std::placeholders::_1,proton.n,T,
		std::ref(th2));
    gd.h=neutron.n/1.0e2;
    gd.deriv_err(neutron.n,fef_nn,mun_num,mun_err);
          
    // Proton chemical potential
    neutron.n=nb*0.51;
    proton.n=nb*0.49;
    T=30.0/hc_mev_fm;
    std::function<double(double)> fef_np=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),neutron.n,
		std::placeholders::_1,T,std::ref(th2));
    gd.h=proton.n/1.0e2;
    gd.deriv_err(proton.n,fef_np,mup_num,mup_err);
          
    // Entropy
    neutron.n=nb*0.51;
    proton.n=nb*0.49;
    T=30.0/hc_mev_fm;
    std::function<double(double)> fef_en=
      std::bind(std::mem_fn<double(fermion &,fermion &,double,
				   double, double, thermo &)>
		(&eos::free_energy_density_alt),this,
		std::ref(neutron),std::ref(proton),neutron.n,proton.n,
		std::placeholders::_1,std::ref(th2));
    gd.h=T/1.0e2;
    gd.deriv_err(T,fef_en,en_num,en_err);
    en_num*=-1.0;

    // Now compute analytical results from function
    neutron.n=nb*0.51;
    proton.n=nb*0.49;
    T=30.0/hc_mev_fm;
    free_energy_density(neutron,proton,T,th2);
          
    double rat1=fabs((mun_num-neutron.mu)/mun_err);
    double rat2=fabs((mup_num-proton.mu)/mup_err);
    double rat3=fabs((en_num-th2.en)/en_err);
    if (count%10==0 || rat1>10.0 || rat2>10.0 || rat3>10.0) {
      cout << nb << " ";
      cout << rat1 << " ";
      cout << rat2 << " ";
      cout << rat3 << endl;
    }
    avg1+=rat1;
    avg2+=rat2;
    avg3+=rat3;
    avg_sum+=avg1+avg2+avg3;

    t.test_abs(neutron.mu,mun_num,mun_err*err_fac,"mun, T=30, Ye=0.49");
    t.test_abs(proton.mu,mup_num,mup_err*err_fac,"mup, T=30, Ye=0.49");
    t.test_abs(th2.en,en_num,en_err*err_fac,"en, T=30, Ye=0.49");
  }
  cout << "  averages: " << avg1/((double)count) << " ";
  cout << avg2/((double)count) << " ";
  cout << avg3/((double)count) << endl;
  cout << endl;
  cout << "avg_sum=" << avg_sum << endl;
  cout << endl;

  // ---------------------------------------------------------------------
  // Summarize

  t.report();
  
  return 0;
}

int eos::eos_sn(std::vector<std::string> &sv, bool itive_com) {

  eos_sn_oo eso;
  eso.verbose=0;
  hdf_file hf;

  /* Attempt to download the EOS file
   */
  cloud_file cf;
  cf.verbose=2;
  std::string sha=((std::string)"d8c4d4f1315942a663e96fc6452f66d90fc")+
    "87f283e0ed552c8141d1ddba34c19";
  cf.hash_type=cloud_file::sha256;
  cf.hdf5_open_hash(hf,"LS220_234r_136t_50y_analmu_20091212_SVNr26.h5",
		    ((string)"https://isospin.roam.utk.edu/")+
		    "public/eos_tables/scollapse/LS220_234r_136t_50y_"+
		    "analmu_20091212_SVNr26.h5",sha,"data");
  
  vector<vector<double> > pts={{0.16,0.01,0.1},
			       {0.16,0.01,10.0},
			       {0.48,0.5,0.1},
			       {0.48,0.01,0.1},
			       {0.01,0.5,10.0},
			       {0.01,0.5,0.1},
			       {0.004,0.5,10.0},
			       {0.004,0.5,0.1},
			       {0.001,0.5,10.0},
			       {0.001,0.5,0.1},
			       {0.0004,0.5,10.0},
			       {0.0004,0.5,0.1},
			       {0.0001,0.5,10.0},
			       {0.0001,0.5,0.1}};
  
  // -----------------------------------------------------------------
  // LS220
  
  double f, F_eg, fint;

  eso.load("data/LS220_234r_136t_50y_analmu_20091212_SVNr26.h5",
	   eos_sn_oo::ls_mode);
  
  cout << endl;
  cout << "LS220: " << endl;
  cout << endl;

  for(size_t i=0;i<pts.size();i++) {
    double lnB=pts[i][0];
    double lYe=pts[i][1];
    double lT=pts[i][2];
    
    cout << "nB: " << lnB << " Ye: " << lYe
	 << " T: " << lT << endl;
    f=eso.F.interp_linear(lnB,lYe,lT);
    thermo lep;
    double mue2;
    eso.compute_eg_point(lnB,lYe,lT,lep,mue2);
    F_eg=(lep.ed-lep.en*lT)/lnB;
    cout << "F_full,F_eg,Xn,Xa: " 
	 << f << " " << F_eg << " "
	 << eso.Xn.interp_linear(lnB,lYe,lT) << " "
	 << eso.Xalpha.interp_linear(lnB,lYe,lT) << endl;
    cout << "F_int: " << f-F_eg << endl;
    cout << endl;
  }

  // -----------------------------------------------------------------
  // SFHo (from RMF EOS object)

  cout << endl;
  cout << "SFHo (eos_had_rmf): " << endl;
  cout << endl;

  rmf_load(rmf,"SFHo");
  
  for(size_t i=0;i<pts.size();i++) {
    if (i!=4 && i!=6 && i!=8) {
      double lnB=pts[i][0];
      double lYe=pts[i][1];
      double lT=pts[i][2];
      
      cout << i << " nB: " << lnB << " Ye: " << lYe
	   << " T: " << lT << endl;
      neutron.n=lnB*(1.0-lYe);
      proton.n=lnB*lYe;
      rmf.calc_temp_e(neutron,proton,lT/hc_mev_fm,th2);
      f=(th2.ed-lT/hc_mev_fm*th2.en)*hc_mev_fm/lnB;
      thermo lep;
      double mue2;
      eso.compute_eg_point(lnB,lYe,lT,lep,mue2);
      F_eg=(lep.ed-lep.en*lT)/lnB;
      cout << "F_full,F_eg,Xn,Xa: " 
	   << f+F_eg << " " << F_eg << " " << 0.0 << " "
	   << 0.0 << endl;
      cout << "F_int: " << f << endl;
      cout << endl;
    }
  }

  // -----------------------------------------------------------------
  // SFHo (O'Connor's table)

  eso.load("data/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5",
	   eos_sn_oo::hfsl_mode);
  
  cout << endl;
  cout << "SFHo (O'Connor): " << endl;
  cout << endl;

  for(size_t i=0;i<pts.size();i++) {
    double lnB=pts[i][0];
    double lYe=pts[i][1];
    double lT=pts[i][2];
    
    cout << "nB: " << lnB << " Ye: " << lYe
	 << " T: " << lT << endl;
    f=eso.F.interp_linear(lnB,lYe,lT);
    thermo lep;
    double mue2;
    eso.compute_eg_point(lnB,lYe,lT,lep,mue2);
    F_eg=(lep.ed-lep.en*lT)/lnB;
    cout << "F_full,F_eg,Xn,Xa: " 
	 << f << " " << F_eg << " "
	 << eso.Xn.interp_linear(lnB,lYe,lT) << " "
	 << eso.Xalpha.interp_linear(lnB,lYe,lT) << endl;
    cout << "F_int: " << f-F_eg << endl;
    cout << endl;
  }

  // --------------------------------------------------------------------
  // SFHo (Hempel's table)
  
  eos_sn_hfsl esh;
  esh.verbose=0;
  esh.load("data/sfho_frdm_shen98_v1.03.tab");
  
  cout << endl;
  cout << "SFHo (Hempel): " << endl;
  cout << endl;

  for(size_t i=0;i<pts.size();i++) {
    double lnB=pts[i][0];
    double lYe=pts[i][1];
    double lT=pts[i][2];
    
    cout << "nB: " << lnB << " Ye: " << lYe
	 << " T: " << lT << endl;
    fint=eso.Fint.interp_linear(lnB,lYe,lT);
    thermo lep;
    double mue2;
    eso.compute_eg_point(lnB,lYe,lT,lep,mue2);
    F_eg=(lep.ed-lep.en*lT)/lnB;
    cout << "F_full,F_eg,Xn,Xa: " 
	 << fint+F_eg << " " << F_eg << " "
	 << eso.Xn.interp_linear(lnB,lYe,lT) << " "
	 << eso.Xalpha.interp_linear(lnB,lYe,lT) << endl;
    cout << "F_int: " << fint << endl;
    cout << endl;
  }

  // -----------------------------------------------------------------
  // SFHx (from RMF EOS object)

  cout << endl;
  cout << "SFHx (eos_had_rmf): " << endl;
  cout << endl;

  rmf_load(rmf,"SFHx");
  
  for(size_t i=0;i<pts.size();i++) {
    if (i!=12) {
      double lnB=pts[i][0];
      double lYe=pts[i][1];
      double lT=pts[i][2];
      
      cout << i << " nB: " << lnB << " Ye: " << lYe
	   << " T: " << lT << endl;
      neutron.n=lnB*(1.0-lYe);
      proton.n=lnB*lYe;
      rmf.calc_temp_e(neutron,proton,lT/hc_mev_fm,th2);
      f=(th2.ed-lT/hc_mev_fm*th2.en)*hc_mev_fm/lnB;
      thermo lep;
      double mue2;
      eso.compute_eg_point(lnB,lYe,lT,lep,mue2);
      F_eg=(lep.ed-lep.en*lT)/lnB;
      cout << "F_full,F_eg,Xn,Xa: " 
	   << f+F_eg << " " << F_eg << " " << 0.0 << " "
	   << 0.0 << endl;
      cout << "F_int: " << f << endl;
      cout << endl;
    }
  }

  // -----------------------------------------------------------------
  // SFHx

  eso.load("data/Hempel_SFHxEOS_rho234_temp180_ye60_version_1.1_20120817.h5",
	   eos_sn_oo::hfsl_mode);

  cout << endl;
  cout << "SFHx (O'Connor): " << endl;
  cout << endl;

  for(size_t i=0;i<pts.size();i++) {
    if (i!=12) {
      double lnB=pts[i][0];
      double lYe=pts[i][1];
      double lT=pts[i][2];
      
      cout << "nB: " << lnB << " Ye: " << lYe
	   << " T: " << lT << endl;
      fint=eso.Fint.interp_linear(lnB,lYe,lT);
      thermo lep;
      double mue2;
      eso.compute_eg_point(lnB,lYe,lT,lep,mue2);
      F_eg=(lep.ed-lep.en*lT)/lnB;
      cout << "F_full,F_eg,Xn,Xa: " 
	   << fint+F_eg << " " << F_eg << " "
	   << eso.Xn.interp_linear(lnB,lYe,lT) << " "
	   << eso.Xalpha.interp_linear(lnB,lYe,lT) << endl;
      cout << "F_int: " << fint << endl;
      cout << endl;
    }
  }

  // -----------------------------------------------------------------
  // IUFSU (from RMF EOS object)

  cout << endl;
  cout << "IUFSU (eos_had_rmf): " << endl;
  cout << endl;

  rmf_load(rmf,"IUFSU");
  
  for(size_t i=0;i<pts.size();i++) {
    if (i!=6) {
      double lnB=pts[i][0];
      double lYe=pts[i][1];
      double lT=pts[i][2];
      
      cout << i << " nB: " << lnB << " Ye: " << lYe
	   << " T: " << lT << endl;
      neutron.n=lnB*(1.0-lYe);
      proton.n=lnB*lYe;
      rmf.calc_temp_e(neutron,proton,lT/hc_mev_fm,th2);
      f=(th2.ed-lT/hc_mev_fm*th2.en)*hc_mev_fm/lnB;
      thermo lep;
      double mue2;
      eso.compute_eg_point(lnB,lYe,lT,lep,mue2);
      F_eg=(lep.ed-lep.en*lT)/lnB;
      cout << "F_full,F_eg,Xn,Xa: " 
	   << f+F_eg << " " << F_eg << " " << 0.0 << " "
	   << 0.0 << endl;
      cout << "F_int: " << f << endl;
      cout << endl;
    }
  }

  // -----------------------------------------------------------------
  // IUFSU

  eso.load("data/Hempel_IUFEOS_rho234_temp180_ye60_version_1.1_20140129.h5",
	   eos_sn_oo::hfsl_mode);

  cout << endl;
  cout << "IUFSU (O'Connor): " << endl;
  cout << endl;

  for(size_t i=0;i<pts.size();i++) {
    double lnB=pts[i][0];
    double lYe=pts[i][1];
    double lT=pts[i][2];
    
    cout << "nB: " << lnB << " Ye: " << lYe
	 << " T: " << lT << endl;
    fint=eso.Fint.interp_linear(lnB,lYe,lT);
    thermo lep;
    double mue2;
    eso.compute_eg_point(lnB,lYe,lT,lep,mue2);
    F_eg=(lep.ed-lep.en*lT)/lnB;
    cout << "F_full,F_eg,Xn,Xa: " 
	 << fint+F_eg << " " << F_eg << " "
	 << eso.Xn.interp_linear(lnB,lYe,lT) << " "
	 << eso.Xalpha.interp_linear(lnB,lYe,lT) << endl;
    cout << "F_int: " << fint << endl;
    cout << endl;
  }

  return 0;
}

int eos::solve_Ye(size_t nv, const ubvector &x, ubvector &y,
		  double nb, double T, double muL) {

  // The temperature T should be in 1/fm

  double Ye=x[0];
  neutron.n=nb*(1.0-Ye);
  proton.n=nb*Ye;

  if (use_alt_eos==false) {
    double t1, t2;
    sk.eff_mass(neutron,proton,t1,t2);
    if (neutron.ms<0.0 || proton.ms<0.0) return 1;
    if (neutron.n<0.0 || proton.n<0.0) return 2;
  }

  free_energy_density(neutron,proton,T,th2);
  
  photon.massless_calc(T);

  electron.n=proton.n;
  electron.mu=electron.m; 
  relf.pair_density(electron,T);
  
  y[0]=neutron.mu-proton.mu-electron.mu+muL+neutron.m-proton.m;
  return 0;
}

int eos::solve_fixed_sonb_YL(size_t nv, const ubvector &x, ubvector &y,
			     double nB, double sonb, double YL) {

  // The temperature T should be in 1/fm

  double Ye=x[0];
  double T=x[1];
  
  if (x[0]<1.0e-5 || x[0]>0.6) return 1;
  if (sonb>0.0 && (x[1]<1.0e-5 || x[1]>1.0)) return 2;

  neutron.n=nB*(1.0-Ye);
  proton.n=nB*Ye;

  // This is an additional temporary temperature variable we use
  // to make sure the results are exact for s=0
  double T2=T;
  if (sonb==0.0) T2=0.0;

  if (use_alt_eos==false) {
    double t1, t2;
    sk.eff_mass(neutron,proton,t1,t2);
    if (neutron.ms<0.0 || proton.ms<0.0) return 1;
  }
    
  free_energy_density(neutron,proton,T2,th2);
  
  photon.massless_calc(T2);

  electron.n=proton.n;
  electron.mu=electron.m; 
  relf.pair_density(electron,T2);

  double en=th2.en+electron.en+photon.en;
  
  if (YL<1.0e-4) {
    
    y[0]=proton.mu+proton.m+electron.mu-neutron.mu-neutron.m;
    
  } else {
    
    neutrino.mu=proton.mu+proton.m+electron.mu-neutron.mu-neutron.m;
    relf.massless_pair_mu(neutrino,T2);
    double YL2=(neutrino.n+electron.n)/nB;

    y[0]=YL2-YL;
  }

  if (sonb==0.0) {
    y[1]=T;
  } else {
    y[1]=en/nB-sonb;
  }

  return 0;
}

int eos::solve_T(size_t nv, const ubvector &x, ubvector &y,
		 double nb, double Ye, double sonb) {
  
  // The temperature T should be in 1/fm
  double T=x[0];
  
  neutron.n=nb*(1.0-Ye);
  proton.n=nb*Ye;

  if (use_alt_eos==false) {
    double t1, t2;
    sk.eff_mass(neutron,proton,t1,t2);
    if (neutron.ms<0.0 || proton.ms<0.0) return 1;
  }
  
  free_energy_density(neutron,proton,T,th2);
    
  photon.massless_calc(T);

  electron.n=proton.n;
  electron.mu=electron.m; 
  relf.pair_density(electron,T);

  y[0]=(th2.en+electron.en+photon.en)/nb-sonb;
  return 0;
}

int eos::mcarlo_data(std::vector<std::string> &sv, bool itive_com) {

  table<> t;
  t.line_of_names(((string)"index S L qmc_a qmc_b qmc_alpha ")+
		  "qmc_beta i_ns i_skyrme phi eos_n0 eos_EoA eos_K chi2_ns "+
		  "ns_fit0 ns_fit1 ns_fit2 ns_fit3 ns_fit4 "+
		  "F_0004_50_10 F_016_01_01 F_016_01_10 "+
		  "F_048_01_01 F_048_50_01 F_100_50_10 "+
		  "ns_min_cs2 ns_max_cs2 Lambda_bar_14");

  vector<double> nB_arr={0.004,0.16,0.16,0.48,0.48,1.0};
  vector<double> Ye_arr={0.5,0.01,0.01,0.01,0.5,0.5};
  vector<double> T_arr={10.0,0.1,10.0,0.1,0.1,10.0};

  static const int N=10000;
  for(int j=0;j<N;j++){
    
    std::vector<std::string> obj;
    random(obj,false);

    vector<double> line={((double)j),eos_S,eos_L,qmc_a,qmc_b,qmc_alpha,
			 qmc_beta,((double)i_ns),((double)i_skyrme),phi,
			 eos_n0,eos_EoA,eos_K,chi2_ns,ns_fit_parms[0],
			 ns_fit_parms[1],ns_fit_parms[2],ns_fit_parms[3],
			 ns_fit_parms[4]};
    
    for(size_t k=0;k<6;k++) {
      neutron.n=nB_arr[k]*(1.0-Ye_arr[k]);
      proton.n=nB_arr[k]*Ye_arr[k];
      double T=T_arr[k]/hc_mev_fm;
      line.push_back(free_energy_density(neutron,proton,T,th2)/
                     nB_arr[k]*hc_mev_fm);
    }

    double ns_min_cs2, ns_max_cs2;
    min_max_cs2(ns_min_cs2,ns_max_cs2);
    line.push_back(ns_min_cs2);
    line.push_back(ns_max_cs2);

    line.push_back(Lambda_bar_14);
    
    cout << "Line: ";
    for(size_t i=0;i<line.size();i++) {
      cout << line[i] << " ";
    }
    cout << endl;
    
    t.line_of_data(line.size(),line);
    if (line.size()!=t.get_ncolumns()) {
      O2SCL_ERR("Table sync error in mcarlo_data().",exc_esanity);
    }

    if (j%10==0 || j==N-1) {
      hdf_file hf1;
      std::string fname="mcarlo_data";
      if (sv.size()>1) fname+="_"+sv[1];
      fname+=".o2";
      hf1.open_or_create(fname);
      o2scl_hdf::hdf_output(hf1,t,"mcarlo");
      hf1.close();
    }    
  }
    
  return 0;
}

int eos::vir_comp(std::vector<std::string> &sv, bool itive_com) {

  table<> t, t2;
  t.line_of_names("log_nB F");
  t2.line_of_names("log_nB zn F_vir");

  for(int j=0;j<1000;j++){
    std::vector<std::string> obj;
    random(obj,false);

    for(double nb=1.0e-4;nb<4.001e-1;nb*=pow(4000.0,1.0/99.0)) {
      neutron.n=nb/2.0;
      proton.n=nb/2.0;
      double T=5.0/hc_mev_fm;
      double line[2]={log10(nb),free_energy_density(neutron,proton,T,th2)/
		      nb*hc_mev_fm};
      t.line_of_data(2,line);
      if (j==0) {
	double F_vir=free_energy_density_virial(neutron,proton,T,th2)/
	  nb*hc_mev_fm;
	double zn=exp(neutron.mu/T);
	vector<double> line2={log10(nb),zn,F_vir};
	t2.line_of_data(3,line2);
      }
    }
  }
  
  hdf_file hf;
  hf.open_or_create("vir_comp.o2");
  hdf_output(hf,t,"vir_comp");
  hdf_output(hf,t2,"vir_comp2");
  hf.close();
  
  return 0;
}

int eos::pns_eos(std::vector<std::string> &sv, bool itive_com) {

  table_units<> eost;
  eost.line_of_names("nB Ye T ed pr nn np mun mup ne mue");
  eost.line_of_units(((std::string)"1/fm^3 . 1/fm 1/fm^4 1/fm^4 1/fm^3 ")+
		     "1/fm^3 1/fm 1/fm 1/fm^3 1/fm");

  double sonb=o2scl::stod(sv[1]);
  double YL=o2scl::stod(sv[2]);
  
  ubvector x(2);
  x[0]=0.1;
  x[1]=10.0/hc_mev_fm;

  double Ye0=0.0, T0=0.0;
  
  for(double nB=0.08;nB<1.5;nB+=0.01) {
    
    mroot_hybrids<> mh;
    mm_funct pns=std::bind
      (std::mem_fn<int(size_t,const ubvector &,
		       ubvector &, double, double, double)>
       (&eos::solve_fixed_sonb_YL),
       this,std::placeholders::_1,
       std::placeholders::_2,
       std::placeholders::_3,nB,sonb,YL);
    mh.verbose=0;
    mh.msolve(2,x,pns);

    cout << nB << " " << x[0] << " " << x[1] << endl;
    if (nB<0.08*1.001) {
      Ye0=x[0];
      T0=x[1];
    }
    
    vector<double> line={nB,x[0],x[1],th2.ed+electron.ed+photon.ed+
			 neutron.m*neutron.n+proton.m*proton.n,
			 th2.pr+electron.pr+photon.pr,neutron.n,
			 proton.n,neutron.mu,proton.mu,electron.n,
			 electron.mu};
    eost.line_of_data(line);
    
  }

  x[0]=Ye0;
  x[1]=T0;
  
  eos_tov_interp eti;

  bool cold_crust=true;
  
  if (cold_crust) {
    
    eti.default_low_dens_eos();
    
  } else {

    /* 
       This section is currently unused
    */
    for(double nB=0.07;nB>=1.0e-6;nB/=1.3) {
      
      mroot_hybrids<> mh;
      mm_funct pns=std::bind
	(std::mem_fn<int(size_t,const ubvector &,
			 ubvector &, double, double, double)>
	 (&eos::solve_fixed_sonb_YL),
	 this,std::placeholders::_1,
	 std::placeholders::_2,
	 std::placeholders::_3,nB,sonb,YL);
      mh.verbose=0;
      mh.msolve(2,x,pns);
      
      cout << nB << " " << x[0] << " " << x[1] << endl;
      
      vector<double> line={nB,x[0],x[1],th2.ed+electron.ed+photon.ed+
			   neutron.m*neutron.n+proton.m*proton.n,
			   th2.pr+electron.pr+photon.pr,neutron.n,
			   proton.n,neutron.mu,proton.mu,electron.n,
			   electron.mu};
      eost.line_of_data(line);
    }
    
    eost.sort_table("nB");
    
    exit(-1);
    
  }

  tov_solve ts;
  eti.read_table(eost,"ed","pr","nB");
  ts.set_eos(eti);
  ts.mvsr();

  shared_ptr<table_units<> > mvsrt=ts.get_results();
  cout << mvsrt->get_unit("ed") << endl;

  //for(size_t j=0;j<mvsrt->get_nlines();j++) {
  //cout << mvsrt->get("r",j) << " " << mvsrt->get("gm",j) << endl;
  //}

  hdf_file hf;
  hf.open_or_create(sv[3]);
  hdf_output(hf,eost,"eos");
  hdf_output(hf,*mvsrt,"mvsr");
  hf.close();

  return 0;
}

int eos::select_model(std::vector<std::string> &sv, bool itive_com) {

  i_ns=o2scl::stod(sv[1]);
  i_skyrme=o2scl::stod(sv[2]);
  qmc_alpha=o2scl::stod(sv[3]);
  qmc_a=o2scl::stod(sv[4]);
  eos_L=o2scl::stod(sv[5]);
  eos_S=o2scl::stod(sv[6]);
  phi=o2scl::stod(sv[7]);

  int iret=select_internal(i_ns,i_skyrme,qmc_alpha,qmc_a,eos_L,eos_S,phi);
  if (iret!=0) {
    cerr << "Model is unphysical (iret=" << iret << ")." << endl;
    return 1;
  }
  
  return 0;
}

int eos::select_internal(int i_ns_loc, int i_skyrme_loc,
			 double qmc_alpha_loc, double qmc_a_loc,
			 double eos_L_loc, double eos_S_loc,
			 double phi_loc) {
  
  i_ns=i_ns_loc;
  i_skyrme=i_skyrme_loc;
  qmc_alpha=qmc_alpha_loc;
  qmc_a=qmc_a_loc;
  eos_L=eos_L_loc;
  eos_S=eos_S_loc;
  phi=phi_loc;
  
  model_selected=true;

  ns_fit(i_ns);
  
  double ns_min_cs2, ns_max_cs2;
  min_max_cs2(ns_min_cs2,ns_max_cs2);
  if (ns_min_cs2<0.0) {
    model_selected=false;
    return 1;
  }

  if (9.17*eos_S-266.0>eos_L || 14.3*eos_S-379.0<eos_L) {
    model_selected=false;
    return 2;
  }

  double rho0=UNEDF_tab.get(((string)"rho0"),i_skyrme);
  double K=UNEDF_tab.get(((string)"K"),i_skyrme);
  double Ms_inv=UNEDF_tab.get(((string)"Ms_inv"),i_skyrme);
  double EoA=UNEDF_tab.get(((string)"EoA"),i_skyrme);
  double unedf_a=UNEDF_tab.get(((string)"a"),i_skyrme);
  double unedf_L=UNEDF_tab.get(((string)"L"),i_skyrme);
  
  double Crdr0=UNEDF_tab.get(((string)"Crdr0"),i_skyrme);
  double Crdr1=UNEDF_tab.get(((string)"Crdr1"),i_skyrme);
  double CrdJ0=UNEDF_tab.get(((string)"CrdJ0"),i_skyrme);
  double CrdJ1=UNEDF_tab.get(((string)"CrdJ1"),i_skyrme);
  double Vp=UNEDF_tab.get(((string)"Vp"),i_skyrme);
  double Vn=UNEDF_tab.get(((string)"Vn"),i_skyrme);

  // Store some of the nuclear matter parameters
  eos_n0=rho0;
  eos_EoA=EoA;
  eos_K=K;
    
  // Determine QMC coefficients
  qmc_b=eos_S+EoA-qmc_a;
  qmc_beta=(eos_L/3.0-qmc_a*qmc_alpha)/qmc_b;
    
  if (qmc_b<0.0 || qmc_beta>5.0) {
    model_selected=false;
    return 3;
  }
  
  double Ms_star=1/Ms_inv;

  sk.alt_params_saturation(rho0,EoA/hc_mev_fm,K/hc_mev_fm,Ms_star,
			   eos_S/hc_mev_fm,eos_L/hc_mev_fm,1.0/1.249,
			   Crdr0/hc_mev_fm,Crdr1/hc_mev_fm,
			   CrdJ0/hc_mev_fm,CrdJ1/hc_mev_fm);

  /*
    double f0, g0, f0p, g0p, f1, g1, f1p, g1p;
    sk.landau_nuclear(rho0,939.0/197.33,f0,g0,f0p,g0p,f1,g1,f1p,g1p);
    cout << f0 << " " << g0 << " " << f0p << " " << g0p << " " << f1 << " "
    << g1 << " " << f1p << " " << g1p << endl;
    sk.landau_neutron(rho0,939.0/197.33,f0,g0,f1,g1);
    cout << f0 << " " << g0 << " " << f1 << " " << g1 << endl;
    exit(-1);
  */

  /*
    sk.saturation();
    cout << "Here " << 0.15 << " " << rho0 << " "
    << sk.f_effm_neut(0.15,(0.1441255-0.00758352)/0.15) << " "
    << sk.f_effm_prot(0.15,(0.1441255-0.00758352)/0.15) << endl;
  */
  
  // Test to make sure dineutrons are not bound
  for(double nb=0.01;nb<0.16;nb+=0.001) {
    neutron.n=nb;
    proton.n=0.0;
    sk.calc_e(neutron,proton,th2);
    if (th2.ed/nb<0.0) {
      model_selected=false;
      return 4;
    }
  }

  // Ensure effective masses are positive
  
  fermion &n=neutron;
  fermion &p=proton;
  thermo &th=th2;
  double term, term2;
  
  n.n=1.0;
  p.n=1.0;
  
  term=0.25*(sk.t1*(1.0+sk.x1/2.0)+sk.t2*(1.0+sk.x2/2.0));
  term2=0.25*(sk.t2*(0.5+sk.x2)-sk.t1*(0.5+sk.x1));
  n.ms=n.m/(1.0+2.0*((n.n+p.n)*term+n.n*term2)*n.m);
  p.ms=p.m/(1.0+2.0*((n.n+p.n)*term+p.n*term2)*p.m);
  
  if (n.ms<0.0 || p.ms<0.0) {
    model_selected=false;
    return 5;
  }
  
  n.n=2.0;
  p.n=0.0;
  
  term=0.25*(sk.t1*(1.0+sk.x1/2.0)+sk.t2*(1.0+sk.x2/2.0));
  term2=0.25*(sk.t2*(0.5+sk.x2)-sk.t1*(0.5+sk.x1));
  n.ms=n.m/(1.0+2.0*((n.n+p.n)*term+n.n*term2)*n.m);
  p.ms=p.m/(1.0+2.0*((n.n+p.n)*term+p.n*term2)*p.m);
  
  if (n.ms<0.0 || p.ms<0.0) {
    model_selected=false;
    return 6;
  }
  
  n.n=0.0;
  p.n=2.0;
  
  term=0.25*(sk.t1*(1.0+sk.x1/2.0)+sk.t2*(1.0+sk.x2/2.0));
  term2=0.25*(sk.t2*(0.5+sk.x2)-sk.t1*(0.5+sk.x1));
  n.ms=n.m/(1.0+2.0*((n.n+p.n)*term+n.n*term2)*n.m);
  p.ms=p.m/(1.0+2.0*((n.n+p.n)*term+p.n*term2)*p.m);
  
  if (n.ms<0.0 || p.ms<0.0) {
    model_selected=false;
    return 7;
  }

  // --------------------------------------------------------
  // Test beta equilibrium
  
  // Loop over baryon densities
  cout << "Going to beta-eq test: " << endl;
  for(double nbx=0.1;nbx<2.00001;nbx+=0.05) {
    
    // Beta equilibrium at T=1 MeV
    mm_funct mf=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &,
		       double, double,double)>
       (&eos::solve_Ye),this,std::placeholders::_1, 
       std::placeholders::_2,std::placeholders::_3,nbx,1.0/hc_mev_fm,0.0);
    
    int verbose_store=verbose;
    verbose=0;

    ubvector Ye_trial(1);
    Ye_trial[0]=0.05;
    mroot_hybrids<> mh;
    mh.err_nonconv=false;
    mh.def_jac.err_nonconv=false;
    int ret=mh.msolve(1,Ye_trial,mf);
    if (ret!=0) {
      model_selected=false;
      return 8;
    }
    
    double Ye=Ye_trial[0];
    
    verbose=verbose_store;
    
    if (Ye<0.0 || Ye>1.0) {
      model_selected=false;
      return 9;
    }
    
  }

  // --------------------------------------------------------
  // Test cs2
  
  model_selected=true;

  if (select_cs2_test) {
    cout << "Going to cs2 test: " << endl;
    for(double nbx=0.1;nbx<2.00001;nbx+=0.05) {
      for(double yex=0.05;yex<0.4501;yex+=0.1) {
	for(double Tx=1.0/hc_mev_fm;Tx<10.01/hc_mev_fm;Tx+=9.0/hc_mev_fm) {
	  neutron.n=nbx*(1.0-yex);
	  proton.n=nbx*yex;
	  double cs2x=cs2_func(neutron,proton,Tx,th2);
	  if (cs2x<0.0) {
	    model_selected=false;
	    if (true) {
	      cout << "Negative speed of sound." << endl;
	      cout << nbx << " " << yex << " " << Tx*hc_mev_fm << " "
		   << cs2x << endl;
	      exit(-1);
	    }	    
	    return 10;
	  }
	}
      }
    }
  }

  return 0;
}

int eos::random(std::vector<std::string> &sv, bool itive_com) {

  // This function never fails, and it requires a call to
  // free_energy_density(), so we set this to true
  model_selected=true;

  bool done=false;
  while (done==false) {

    done=true;
    
    if (verbose>0) {
      cout << "Selecting random model." << endl;
    }

    // Select a random value for phi
    phi=rng.random();

    // Random neutron star EOS
    i_ns=rng.random_int(nstar_tab.get_nlines());

    // Select a random QMC two-body interaction
    qmc_alpha=(rng.random()*0.06)+0.47;
    qmc_a=(rng.random()*1.0)+12.5;
   
    // Select a random value of S and L according
    // to perscription in PRC 91, 015804 (2015)
    eos_L=rng.random()*21.0+44.0;
    eos_S=rng.random()*6.6+29.5;

    // Select a random Skyrme model
    i_skyrme=rng.random_int(UNEDF_tab.get_nlines());
    
    if (true || verbose>1) {
      cout << "Trying random model: " << endl;
      cout << "i_ns= " << i_ns << endl;
      cout << "i_skyrme= " << i_skyrme << endl;
      cout << "alpha= " << qmc_alpha << endl;
      cout << "a= " << qmc_a << endl;
      cout << "eos_L= " << eos_L << endl;
      cout << "eos_S= " << eos_S << endl;
      cout << "phi= " << phi << endl;
    }

    int ret=select_internal(i_ns,i_skyrme,qmc_alpha,qmc_a,
			    eos_L,eos_S,phi);
    if (true || verbose>1) {
      if (ret==0) {
	cout << "Success." << endl;
      } else {
	cout << "Failed (" << ret << "). Selecting new random model." << endl;
      }
    }

    if (ret!=0) done=false;
  }
  
  if (true) {
    cout << "Function eos::random() selected parameters: " << endl;
    cout << i_ns << " " << i_skyrme << " " << qmc_alpha
	 << " " << qmc_a << " " << eos_L << " " << eos_S << " "
	 << phi << endl;
  }
  
  return 0;
}

int eos::point(std::vector<std::string> &sv, bool itive_com) {

  double nB=o2scl::stod(sv[1]);
  double Ye=o2scl::stod(sv[2]);
  double T=o2scl::stod(sv[3])/hc_mev_fm;

  neutron.n=nB*(1.0-Ye);
  proton.n=nB*Ye;
  free_energy_density(neutron,proton,T,th2);
  if (verbose>=1) {
    
    cout.setf(ios::showpos);
    double f_total=th2.ed-T*th2.en;
    
    cout << "f_total,F_total         = " << f_total << " 1/fm^4 "
         << f_total/nB*hc_mev_fm << " MeV" << endl;
    cout << endl;
    
    cout << "ed (w/rm), pr= "
         << th2.ed+neutron.n*neutron.m+proton.n*proton.m
         << " 1/fm^4 " << th2.pr << " 1/fm^4" << endl;
    cout << "entropy, s per baryon= " << th2.en << " 1/fm^3 "
         << th2.en/nB << endl;
    cout << "mu_n, mu_p = " << neutron.mu*hc_mev_fm << " MeV "
         << proton.mu*hc_mev_fm << " MeV" << endl;
    cout << endl;
    cout.unsetf(ios::showpos);
    
  }
  
  return 0;
}

int eos::test_eg(std::vector<std::string> &sv,
		 bool itive_com) {

  string fname;
  if (sv.size()>1) {
    fname=sv[1];
    cout << "Will output to file named " << fname << " ." << endl;
  }
  
  // Choose a larger grid for a more complete test
  size_t n_nB=326;
  size_t n_Ye=101;
  size_t n_T=161;
    
  std::string nB_grid_spec2="10^(i*0.04-12)";
  std::string Ye_grid_spec2="0.01*i";
  // This temperature grid includes zero, but is otherwise
  // similar to Hempel's logarithmic temperature grid
  std::string T_grid_spec2="i*5^(i*0.02-3.2)";

  vector<double> nB_grid, T_grid, Ye_grid;
  
  calc_utf8<> calc;
  std::map<std::string,double> vars;
  
  calc.compile(nB_grid_spec2.c_str());
  for(size_t i=0;i<n_nB;i++) {
    vars["i"]=((double)i);
    nB_grid.push_back(calc.eval(&vars));
  }
  
  calc.compile(Ye_grid_spec2.c_str());
  for(size_t i=0;i<n_Ye;i++) {
    vars["i"]=((double)i);
    Ye_grid.push_back(calc.eval(&vars));
  }
  
  calc.compile(T_grid_spec2.c_str());
  for(size_t i=0;i<n_T;i++) {
    vars["i"]=((double)i);
    T_grid.push_back(calc.eval(&vars));
  }

  size_t size_arr[3]={n_nB,n_Ye,n_T};
  vector<vector<double> > grid_arr={nB_grid,Ye_grid,T_grid};

  tensor_grid3<> t_F(n_nB,n_Ye,n_T);
  t_F.set_grid(grid_arr);
  tensor_grid3<> t_E(n_nB,n_Ye,n_T);
  t_E.set_grid(grid_arr);
  tensor_grid3<> t_P(n_nB,n_Ye,n_T);
  t_P.set_grid(grid_arr);
  tensor_grid3<> t_S(n_nB,n_Ye,n_T);
  t_S.set_grid(grid_arr);
  tensor_grid3<> t_mue(n_nB,n_Ye,n_T);
  t_mue.set_grid(grid_arr);
  
  eos_sn_base eso;
  eso.include_muons=include_muons;
  elep.include_muons=include_muons;

  for(size_t i=0;i<n_nB;i++) {
    double nB=nB_grid[i];
    if (i%10==0) {
      cout << "i,nB: ";
      cout.width(3);
      cout << i << " " << nB << endl;
    }
    for(size_t j=1;j<n_Ye-1;j++) {
      double Ye=Ye_grid[j];

      // Construct a new guess for the electron chemical
      elep.e.mu=elep.e.m;
      elep.e.n=nB*Ye/2.0;
  
      for(size_t k=0;k<n_T;k++) {
        
	double T_MeV=T_grid[k];
	thermo lep;

        //elep.verbose=2;
        elep.pair_density_eq(nB*Ye,T_MeV/hc_mev_fm);
        if (verbose>1) {
          cout << nB << " " << Ye << " " << T_MeV << " " << elep.e.n << " "
               << elep.mu.n << " " << elep.e.n+elep.mu.n << endl;
        }
        t_F.set(i,j,k,(hc_mev_fm*elep.th.ed-T_grid[k]*elep.th.en)/nB);
        t_E.set(i,j,k,hc_mev_fm*elep.th.ed/nB);
        t_P.set(i,j,k,hc_mev_fm*elep.th.pr);
        t_S.set(i,j,k,hc_mev_fm*elep.th.en/nB);
        t_mue.set(i,j,k,hc_mev_fm*elep.e.mu);

      }
    }
  }

  if (fname.length()>0) {
    cout << "Writing file " << fname << " ." << endl;
    hdf_file hf;
    hf.open_or_create(fname);
    hf.seti("include_muons",include_muons);
    hf.set_szt("n_nB",n_nB);
    hf.set_szt("n_Ye",n_Ye);
    hf.set_szt("n_T",n_T);
    hf.setd_vec("nB_grid",nB_grid);
    hf.setd_vec("Ye_grid",Ye_grid);
    hf.setd_vec("T_grid",T_grid);
    hdf_output(hf,t_F,"F_eg");
    hdf_output(hf,t_E,"E_eg");
    hdf_output(hf,t_P,"P_eg");
    hdf_output(hf,t_S,"S_eg");
    hdf_output(hf,t_mue,"mue");
    hf.close();
  }
  
  return 0;
}

int eos::vir_fit(std::vector<std::string> &sv,
		 bool itive_com) {
  ecv.fit(true);
  return 0;
}

int eos::alt_model(std::vector<std::string> &sv,
                     bool itive_com) {

  if (sv.size()<2) {
    cout << "Not enough parameters for alt-model." << endl;
  }
  
  if (sv[1]=="Skyrme" || sv[1]=="skyrme") {
    if (sv.size()<3) {
      cout << "Not enough parameters for alt-model of type Skyrme."
           << endl;
    }
    o2scl_hdf::skyrme_load(sk_alt,sv[2]);
    eosp_alt=&sk_alt;
    use_alt_eos=true;
  } else if (sv[1]=="RMF" || sv[1]=="rmf") {
    o2scl_hdf::rmf_load(rmf,sv[2]);
    eosp_alt=&rmf;
    use_alt_eos=true;
  } else if (sv[1]=="RMFH" || sv[1]=="rmfh") {
    o2scl_hdf::rmf_load(rmf_hyp,sv[2]);
    eosp_alt=&rmf;
    use_alt_eos=true;
  } else {
    cerr << "Model " << sv[1] << " not understood." << endl;
    return 2;
  }
  return 0;
}

void eos::setup_cli(o2scl::cli &cl, bool read_docs) {

  static const int nopt=15;
  o2scl::comm_option_s options[nopt]=
    {{0,"test-deriv","",0,0,"","",
       new o2scl::comm_option_mfptr<eos>
       (this,&eos::test_deriv),o2scl::cli::comm_option_both,
       1,"","eos","test_deriv","doc/xml/classeos.xml"},
     {0,"table-Ye","",2,2,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::table_Ye),o2scl::cli::comm_option_both,
      1,"","eos","table_Ye","doc/xml/classeos.xml"},
     {0,"table-nB","",2,2,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::table_nB),o2scl::cli::comm_option_both,
      1,"","eos","table_nB","doc/xml/classeos.xml"},      
     {0,"table-full","",1,1,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::table_full),o2scl::cli::comm_option_both,
      1,"","eos","table_full","doc/xml/classeos.xml"},      
     {0,"vir-fit","",0,0,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::vir_fit),o2scl::cli::comm_option_both,
      1,"","eos","vir_fit","doc/xml/classeos.xml"},      
     {0,"eos-sn","",0,0,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::eos_sn),o2scl::cli::comm_option_both,
      1,"","eos","eos_sn","doc/xml/classeos.xml"},      
     {0,"mcarlo-data","",0,1,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::mcarlo_data),o2scl::cli::comm_option_both,
      1,"","eos","mcarlo_data","doc/xml/classeos.xml"},      
     {0,"point","",0,3,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::point),o2scl::cli::comm_option_both,
      1,"","eos","point","doc/xml/classeos.xml"},      
     {0,"random","",0,0,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::random),o2scl::cli::comm_option_both,
      1,"","eos","random","doc/xml/classeos.xml"},      
     {0,"pns-eos","",3,3,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::pns_eos),o2scl::cli::comm_option_both,
      1,"","eos","pns_eos","doc/xml/classeos.xml"},      
     {0,"select-model","",7,7,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::select_model),o2scl::cli::comm_option_both,
      1,"","eos","select_model","doc/xml/classeos.xml"},
     {0,"test-eg","",0,1,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::test_eg),o2scl::cli::comm_option_both,
      1,"","eos","test_eg","doc/xml/classeos.xml"},      
     {0,"vir-comp","",0,0,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::vir_comp),o2scl::cli::comm_option_both,
      1,"","eos","vir_comp","doc/xml/classeos.xml"},      
     {0,"alt-model","",1,2,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::alt_model),o2scl::cli::comm_option_both,
      1,"","eos","alt_model","doc/xml/classeos.xml"},
     {0,"hrg-load","",1,1,"","",
      new o2scl::comm_option_mfptr<eos>
      (this,&eos::hrg_load),o2scl::cli::comm_option_both,
      1,"","eos","hrg_load","doc/xml/classeos.xml"}
    };

  cl.doc_o2_file="data/eos_nuclei_docs.o2";
  
  p_verbose.i=&verbose;
  p_verbose.help="";
  p_verbose.doc_class="eos";
  p_verbose.doc_name="verbose";
  p_verbose.doc_xml_file="doc/xml/classeos.xml";
  cl.par_list.insert(make_pair("verbose",&p_verbose));

  p_old_ns_fit.b=&old_ns_fit;
  p_old_ns_fit.help="";
  p_old_ns_fit.doc_class="eos";
  p_old_ns_fit.doc_name="old_ns_fit";
  p_old_ns_fit.doc_xml_file="doc/xml/classeos.xml";
  cl.par_list.insert(make_pair("old_ns_fit",&p_old_ns_fit));

  p_ns_record.b=&ns_record;
  p_ns_record.help="";
  p_ns_record.doc_class="eos";
  p_ns_record.doc_name="ns_record";
  p_ns_record.doc_xml_file="doc/xml/classeos.xml";
  cl.par_list.insert(make_pair("ns_record",&p_ns_record));

  p_include_muons.b=&include_muons;
  p_include_muons.help="";
  p_include_muons.doc_class="eos";
  p_include_muons.doc_name="include_muons";
  p_include_muons.doc_xml_file="doc/xml/classeos.xml";
  cl.par_list.insert(make_pair("include_muons",&p_include_muons));

  p_select_cs2_test.b=&select_cs2_test;
  p_select_cs2_test.help="";
  p_select_cs2_test.doc_class="eos";
  p_select_cs2_test.doc_name="select_cs2_test";
  p_select_cs2_test.doc_xml_file="doc/xml/classeos.xml";
  cl.par_list.insert(make_pair("select_cs2_test",&p_select_cs2_test));

  p_test_ns_cs2.b=&test_ns_cs2;
  p_test_ns_cs2.help="";
  p_test_ns_cs2.doc_class="eos";
  p_test_ns_cs2.doc_name="test_ns_cs2";
  p_test_ns_cs2.doc_xml_file="doc/xml/classeos.xml";
  cl.par_list.insert(make_pair("test_ns_cs2",&p_test_ns_cs2));

  p_use_alt_eos.b=&use_alt_eos;
  p_use_alt_eos.help="";
  p_use_alt_eos.doc_class="eos";
  p_use_alt_eos.doc_name="use_alt_eos";
  p_use_alt_eos.doc_xml_file="doc/xml/classeos.xml";
  cl.par_list.insert(make_pair("use_alt_eos",&p_use_alt_eos));

  p_a_virial.d=&a_virial;
  p_a_virial.help="";
  p_a_virial.doc_class="eos";
  p_a_virial.doc_name="a_virial";
  p_a_virial.doc_xml_file="doc/xml/classeos.xml";
  cl.par_list.insert(make_pair("a_virial",&p_a_virial));

  p_b_virial.d=&b_virial;
  p_b_virial.help="";
  p_b_virial.doc_class="eos";
  p_b_virial.doc_name="b_virial";
  p_b_virial.doc_xml_file="doc/xml/classeos.xml";
  cl.par_list.insert(make_pair("b_virial",&p_b_virial));

  cl.set_comm_option_vec(nopt,options);
  
  if (read_docs) {
    
    // The default xml-to-o2 function overwrites all the eos_nuclei
    // documentation so we disable it here if we're not calling
    // this function from there.
    cl.remove_comm_option("xml-to-o2");
    
    cl.read_docs();
  }

  return;
}
