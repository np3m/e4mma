/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2024, Zidu Lin and Andrew W. Steiner
  
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
#include "eos_nuclei.h"

#include "neutrino/tensor.h"
#include "neutrino/Polarization.hpp"
#include "neutrino/PolarizationNonRel.hpp"
#include "neutrino/constants.h"

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;
using namespace nuopac;

int eos_nuclei::mcarlo_nuclei(std::vector<std::string> &sv, bool itive_com) {

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
  for(int j=0;j<N;j++) {

    cout << "j: " << j << endl;
    
    if (j==0) {
      select_internal(470,738,0.5,13.0,62.4,32.8,0.9);
    } else {
      std::vector<std::string> obj;
      random(obj,false);
    }

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

    // This value is set in EOS
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

int eos_nuclei::mcarlo_nuclei2(std::vector<std::string> &sv, 
			       bool itive_com) {
  
  /*
    t.line_of_names(((string)"index S L qmc_a qmc_b qmc_alpha ")+
    "qmc_beta i_ns i_skyrme phi eos_n0 eos_EoA eos_K "+
    "chi2_ns ns_fit0 ns_fit1 ns_fit2 ns_fit3 ns_fit4 "+
    "Fint Eint Sint Pint munfull mupfull log_xn log_xp Z A "+
    "Amin Amax NmZ_max NmZ_min Xn Xp Xalpha Xd Xt XHe3 XLi4 "+
    "Xnuclei ns_min_cs2 ns_max_cs2 Lambda_bar_14");
  */

  double nB=function_to_double(sv[1]);
  double Ye=function_to_double(sv[2]);
  double T=function_to_double(sv[3])/hc_mev_fm;
  
  if (loaded==false) {
    cerr << "Requires EOS guess." << endl;
    return 2;
  }

  size_t inB=vector_lookup(n_nB2,nB_grid2,nB);
  nB=nB_grid2[inB];
  size_t iYe=vector_lookup(n_Ye2,Ye_grid2,Ye);
  Ye=Ye_grid2[iYe];
  size_t iT=vector_lookup(n_T2,T_grid2,T*hc_mev_fm);
  T=T_grid2[iT]/hc_mev_fm;
  vector<size_t> ix={inB,iYe,iT};

  if (tg_flag.get(ix)<9.9) {
    cerr << "Point not converged." << endl;
    return 3;
  }

  table<> tab;
  tab.line_of_names("log_xn log_xp Z N A ZoA Xnuclei");
  
  static const int N=o2scl::stoi(sv[4]);
  for(int j=0;j<N;j++) {
    std::cout << "j: " << j << endl;
    
    std::vector<std::string> obj;
    random(obj,false);
      
    double log_xn, log_xp;
    int A_min=((int)(tg_A_min.get(ix)));
    int A_max=((int)(tg_A_max.get(ix)));
    int NmZ_min=((int)(tg_NmZ_min.get(ix)));
    int NmZ_max=((int)(tg_NmZ_max.get(ix)));
    double mun_full, mup_full;
    thermo thx;
    bool dist_changed=true;
    bool no_nuclei=false;

    double log_xn_best=0.0, log_xp_best=0.0;
    double fr_best=1.0e10;
    double Zbar, Nbar;
    
    for(size_t k=0;k<10;k++) {
      
      log_xn=tg_log_xn.get(ix)+(rng.random()*6.0-3.0);
      log_xp=tg_log_xp.get(ix)+(rng.random()*6.0-3.0);
      
      map<string,double> vdet;
      int ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
			    thx,mun_full,mup_full,
			    A_min,A_max,NmZ_min,NmZ_max,vdet,
			    dist_changed,no_nuclei);
      cout << "ret,log_xn,log_xp,Z,N,A,Z/A:\n  "
	   << ret << " " << log_xn << " " << log_xp << " "
	   << Zbar << " " << Nbar << " " << Zbar+Nbar << " " 
	   << Zbar/(Zbar+Nbar) << endl;
      
      if (ret==0 && thx.ed-T*thx.en<fr_best) {
	fr_best=thx.ed-T*thx.en;
	log_xn_best=log_xn;
	log_xp_best=log_xp;        
      }
    }

    log_xn=log_xn_best;
    log_xp=log_xp_best;
    
    map<string,double> vdet;
    int ret2=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
			   thx,mun_full,mup_full,
			   A_min,A_max,NmZ_min,NmZ_max,vdet,
			   dist_changed,no_nuclei);
    ubvector X;
    compute_X(nB,X);
    
    cout << "AWS: " << ret2 << " " << log_xn << " " << log_xp << " "
	 << Zbar << " " << Nbar << " " << Zbar+Nbar << " " 
	 << Zbar/(Zbar+Nbar) << " " << X[5] << endl;
    
    if (ret2==0) {
      vector<double> line={log_xn,log_xp,Zbar,Nbar,Zbar+Nbar,
        Zbar/(Zbar+Nbar),X[5]};
      tab.line_of_data(7,line);
    }

    if (j%100==0) {
      hdf_file hf;
      hf.open_or_create(sv[5]);
      hdf_output(hf,tab,"mn2");
      hf.close();
    }
    
  }

  return 0;
}

int eos_nuclei::mcarlo_beta(std::vector<std::string> &sv, 
                            bool itive_com) {
  
  if (loaded==false) {
    cerr << "Requires EOS guess." << endl;
    return 2;
  }

  size_t n_point=3;
  if (sv.size()>=3) {
    n_point=stoszt(sv[2]);
  }

  /*
    Questions for Zidu:
    1) is there really a vec and axvec mfp? does it make sense to
    combine them?
    2) in the NC response, I do piRPAvec=2*piRPAvec*fermiF, should I
    do that in CC as well?
    3) Check virial result
    4) Make sure we're tabulating the right response functions

    Todos:
    0) Fix the nucleon effective masses at low density
  */
  
  map<string,double> vdet;
  
  table_units<> tab;
  tab.line_of_names("i_ns i_skyrme qmc_alpha qmc_a phi ");
  tab.line_of_names("t0 t1 t2 t3 x0 x1 x2 x3 epsilon ");
  tab.line_of_units(((string)". . MeV MeV . 1/fm^2 1/fm^4 1/fm^4 ")+
                    "1/fm^(3*a+2) . . . . .");
  
  vector<string> col_list={"nB","nn","np",
    "g","dgdnn","dgdnp","msn","msp","mun","mup","mu_n_nonint",
    "mu_p_nonint","mue",
    "U2","U4","log_xn","log_xp","Z","N",
    "A","ZoA","Xnuclei","Ye_best",
    "fnn_sk","fpp_sk","fnp_sk",
    "gnn_sk","gpp_sk","gnp_sk",
    "fnn_virial","fpp_virial","fnp_virial",
    "gnn_virial","gpp_virial","gnp_virial",
    "fnn","fpp","fnp","fnn_dg0","fpp_dg0","fnp_dg0",
    "gnn","gpp","gnp",
    "vf","vf_dg0","vgt",
    "cc_vec_imfp","cc_vec_imfp_dg0","cc_axvec_imfp",
    "cc_vec_imfp_norpa","cc_axvec_imfp_norpa",
    "nc_vec_imfp","nc_vec_imfp_dg0","nc_axvec_imfp",
    "nc_vec_imfp_norpa","nc_axvec_imfp_norpa"};
  
  vector<string> unit_list={"1/fm^3","1/fm^3","1/fm^3",
    "","1/MeV^3","1/MeV^3","MeV","MeV","MeV","MeV","MeV",
    "MeV","MeV",
    "MeV","MeV","","","","",
    "","","","",
    "1/MeV^2","1/MeV^2","1/MeV^2",
    "1/MeV^2","1/MeV^2","1/MeV^2",
    "1/MeV^2","1/MeV^2","1/MeV^2",
    "1/MeV^2","1/MeV^2","1/MeV^2",
    "1/MeV^2","1/MeV^2","1/MeV^2","1/MeV^2","1/MeV^2","1/MeV^2",
    "1/MeV^2","1/MeV^2","1/MeV^2",
    "1/MeV^2","1/MeV^2","1/MeV^2",
    "1/cm","1/cm","1/cm","1/cm","1/cm",
    "1/cm","1/cm","1/cm","1/cm","1/cm"};

  if (unit_list.size()!=col_list.size()) {
    cout << col_list.size() << " " << unit_list.size() << endl;
    O2SCL_ERR("Table sync 1.",o2scl::exc_einval);
  }

  for(size_t ipoint=0;ipoint<n_point;ipoint++) {

    for(size_t ik=0;ik<col_list.size();ik++) {
    std:string temp=col_list[ik]+"_"+o2scl::szttos(ipoint);
      tab.new_column(temp);
      tab.set_unit(temp,unit_list[ik]);
    }
    if (n_point<5) {
      for(size_t ik=0;ik<100;ik++) {
        tab.new_column(((string)"nc_piRPAvec_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
        tab.new_column(((string)"nc_piRPAax_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
        tab.new_column(((string)"nc_resp_RPAvec_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
        tab.new_column(((string)"nc_resp_RPAax_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
      }
      for(size_t ik=0;ik<100;ik++) {
        tab.new_column(((string)"cc_piRPAvec_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
        tab.new_column(((string)"cc_piRPAax_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
        tab.new_column(((string)"cc_resp_RPAvec_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
        tab.new_column(((string)"cc_resp_RPAax_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
      }
    }
  }
  
  // 1.0e-4 is well into the virial region, 5.0e-3 gives g \approx 0.6,
  // and 0.15 is near saturation density and far from the virial region
  
  vector<double> nB_list={1.0e-4,5.0e-3,0.15};
  vector<double> TMeV_list={10,10,10};
  include_detail=true;

  if (n_point>5) {
    nB_list.clear();
    TMeV_list.clear();
    for(size_t j=0;j<100;j++) {
      nB_list.push_back(1.0e-4*pow(0.15/1.0e-4,((double)j)/99.0));
      TMeV_list.push_back(10.0);
    }
  }
  
  static const int N=10000;
  for(int j=0;j<N;j++) {

    std::cout << "j: " << j << endl;

    if (j==0) {
      use_alt_eos=true;
      vector<string> args={"alt-model","Skyrme","NRAPR"};
      alt_model(args,true);
    } else if (j==1) {
      use_alt_eos=true;
      vector<string> args={"alt-model","Skyrme","SGII"};
      alt_model(args,true);
    } else if (j==2) {
      use_alt_eos=true;
      vector<string> args={"alt-model","Skyrme","UNEDF0"};
      alt_model(args,true);
    } else if (j==3) {
      use_alt_eos=true;
      vector<string> args={"alt-model","Skyrme","UNEDF2"};
      alt_model(args,true);
    } else if (j==4) {
      use_alt_eos=true;
      vector<string> args={"alt-model","Skyrme","SV-min"};
      alt_model(args,true);
    } else {
      use_alt_eos=false;
      // Create a random EOS
      std::vector<std::string> obj;
      random(obj,false);
    }
    if (use_alt_eos) {
      // Copy the couplings to the 'sk' object so we can use
      // those for the Fermi Liquid parameters
      sk.t0=sk_alt.t0;
      sk.t1=sk_alt.t1;
      sk.t2=sk_alt.t2;
      sk.t3=sk_alt.t3;
      sk.x0=sk_alt.x0;
      sk.x1=sk_alt.x1;
      sk.x2=sk_alt.x2;
      sk.x3=sk_alt.x3;
      sk.alpha=sk_alt.alpha;
      cout << "t0,t1: " << sk.t0*hc_mev_fm << " " << sk.t1*hc_mev_fm << endl;
      cout << "t2,t3: " << sk.t2*hc_mev_fm << " " << sk.t3*hc_mev_fm << endl;
      cout << "x0,x1: " << sk.x0 << " " << sk.x1 << endl;
      cout << "x2,x3: " << sk.x2 << " " << sk.x3 << endl;
      cout << "alpha: " << sk.alpha << endl;
    }

    vector<double> line={((double)i_ns),((double)i_skyrme),
      qmc_alpha,qmc_a,phi,
      sk.t0*hc_mev_fm,sk.t1*hc_mev_fm,
      sk.t2*hc_mev_fm,sk.t3*hc_mev_fm,
      sk.x0,sk.x1,sk.x2,sk.x3,sk.alpha};

    if (true) {
      hdf_file hf;
      hf.open_or_create(sv[1]);
      hf.setd_vec("nB_list",nB_list);
      hf.close();
    }
    
    for(size_t ipoint=0;ipoint<n_point;ipoint++) {      

      double nB=nB_list[ipoint];
      double T=TMeV_list[ipoint]/hc_mev_fm;
      double T_MeV=TMeV_list[ipoint];

      size_t inB=vector_lookup(n_nB2,nB_grid2,nB);
      //nB=nB_grid2[inB];
      size_t iT=vector_lookup(n_T2,T_grid2,T*hc_mev_fm);
      //T=T_grid2[iT]/hc_mev_fm;
      
      cout << "Using nB = " << nB << " 1/fm^3 and T = " << T*hc_mev_fm
           << " MeV for\n  ipoint = " << ipoint << " out of total "
           << n_point << endl;
      
      double log_xn_best=0.0, log_xp_best=0.0;
      double fr_best=1.0e10;
      size_t iYe_best=0;
      double log_xn, log_xp;
      double Zbar, Nbar;
      thermo thx;
      double mun_full, mup_full;
      bool dist_changed=true;
      bool no_nuclei=false;
    
      eos_sn_base eso;
      eso.include_muons=false;
      thermo lep;
      
      size_t Ye_min, Ye_max;
      if (n_point>5) {
        double Ye_min_d=29.0*exp(-((double)(ipoint))/24.0);
        double Ye_max_d=10.0+29.0*exp(-((double)(ipoint))/24.0);
        if (Ye_min_d<0.0) Ye_min_d=0.0;
        if (Ye_max_d>40.0) Ye_max_d=40.0;
        Ye_min=((size_t)Ye_min_d);
        Ye_max=((size_t)Ye_max_d);
        cout << ipoint << " " << Ye_min_d << " " << Ye_max_d << " "
             << Ye_min << " " << Ye_max << endl;
      } else {
        if (ipoint==0) {
          Ye_min=30;
          Ye_max=40;
        } else if (ipoint==1) {
          Ye_min=2;
          Ye_max=12;
        } else {
          Ye_min=0;
          Ye_max=10;
        }
      }

      double Ye_best;
      
      if (false) {
        
        for(size_t iYe=Ye_min;iYe<Ye_max;iYe++) {
          vector<size_t> ix={inB,iYe,iT};
          
          //for(size_t iYe=0;iYe<10;iYe++) {
          double Ye=Ye_grid2[iYe];
          
          if (tg_flag.get(ix)<9.9) {
            cerr << "Point not converged." << endl;
            return 3;
          }
          
          int A_min=((int)(tg_A_min.get(ix)));
          int A_max=((int)(tg_A_max.get(ix)));
          int NmZ_min=((int)(tg_NmZ_min.get(ix)));
          int NmZ_max=((int)(tg_NmZ_max.get(ix)));
          
          double fr_Ye=1.0e100;
          
          for(size_t k=0;k<10;k++) {
            
            // Ensure a random initial guess
            
            //log_xn=tg_log_xn.get(ix)+(rng.random()*6.0-3.0);
            //log_xp=tg_log_xp.get(ix)+(rng.random()*6.0-3.0);
            log_xn=tg_log_xn.get(ix)+(rng.random()*0.5-0.25);
            log_xp=tg_log_xp.get(ix)+(rng.random()*0.5-0.25);
            
            // Compute the EOS

            int ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
                                  thx,mun_full,mup_full,
                                  A_min,A_max,NmZ_min,NmZ_max,vdet,
                                  dist_changed,no_nuclei);
            
            
            if (false) {
              cout << "ret,log_xn,log_xp,Z,N,A,Z/A,fr:\n  "
                   << ret << " " << log_xn << " " << log_xp << " "
                   << Zbar << " " << Nbar << " " << Zbar+Nbar << " " 
                   << Zbar/(Zbar+Nbar) << " "
                   << thx.ed-T*thx.en << endl;
            }
            
            // Compute the number density of free neutrons and protons
            double n_fraction, p_fraction;
            if (nB<0.16) {
              double xn=0.0;
              if (log_xn>-300.0) {
                xn=pow(10.0,log_xn);
              }
              double xp=0.0;
              if (log_xp>-300.0) {
                xp=pow(10.0,log_xp);
              }
              double n0_loc=0.16;
              neutron.n=xn*(1.0-nB/n0_loc)/(1.0-nB*xn/n0_loc-nB*xp/n0_loc)*nB;
              proton.n=xp*(1.0-nB/n0_loc)/(1.0-nB*xn/n0_loc-nB*xp/n0_loc)*nB;    
            } else {
              neutron.n=(1.0-Ye)*nB;
              proton.n=Ye*nB;
            }
            
            double mun_gas=vdet["mun_gas"];
            double mup_gas=vdet["mup_gas"];
            
            // Make sure to compute kf, which is not always computed at
            // finite temperature
            sk.def_fet.kf_from_density(neutron);
            sk.def_fet.kf_from_density(proton);
            
            if (false) {
              cout << "mu2,mu4,mu2-mu4: " << mun_gas*hc_mev_fm << " "
                   << mup_gas*hc_mev_fm << " "
                   << (mun_gas-mup_gas)*hc_mev_fm << endl;
            }
            
            // Add the electrons
            double mue=electron.m;
            eso.compute_eg_point(nB_grid2[inB],Ye_grid2[iYe],
                                 T_grid2[iT],lep,mue);
            
            // Copy the electron results to the local electron object
            electron.n=nB_grid2[inB]*Ye_grid2[iYe];
            electron.mu=mue;
            
            if (false) {
              cout << "mue: " << mue*hc_mev_fm << endl;
            }
            
            // Compute the total free energy
            double fr=thx.ed-T*thx.en+lep.ed-T*lep.en;
            
            // Compute the lowest free energy for this value of Ye
            if (ret==0 && fr<fr_Ye) {
              fr_Ye=fr;
            }
            
            // Compute the overall lowest free energy
            if (ret==0 && fr<fr_best) {
              fr_best=fr;
              log_xn_best=log_xn;
              log_xp_best=log_xp;
              iYe_best=iYe;
            }
            //cout << "iYe_best: " << iYe_best << endl;
            
            // End of for(size_t k=0;k<10;k++) {
          }
          
          cout << "iYe,Ye,iYe_best,fr_Ye: " << iYe << " " << Ye << " "
               << iYe_best << " " << fr_Ye << endl;
          cout << ipoint << " " << Ye_min << " " << Ye_max << " "
               << iYe_best << endl;
          if (false) {
            cout << endl;
          }
          
          // End of loop for(size_t iYe=0;iYe<n_Ye2;iYe++) {
        }
        
        Ye_best=Ye_grid2[iYe_best];
        
      } else {
        
        double ii=log10(nB)*31.485+126;
        Ye_best=0.05+0.28*exp(-ii/24.0);
        iYe_best=((size_t)(Ye_best*100.0-1.0));
        cout << "iYe_best,Ye_best: " << iYe_best << " " << Ye_best << endl;
      }
      
      // Now compute the EOS at the optimal Ye
      
      vector<size_t> ix_best={inB,iYe_best,iT};
      int A_min_best=((int)(tg_A_min.get(ix_best)));
      int A_max_best=((int)(tg_A_max.get(ix_best)));
      int NmZ_min_best=((int)(tg_NmZ_min.get(ix_best)));
      int NmZ_max_best=((int)(tg_NmZ_max.get(ix_best)));

      int ret2=10;
      for(size_t k=0;k<10 && ret2!=0;k++) {
        log_xn=tg_log_xn.get(ix_best)+(rng.random()*0.5-0.25);
        log_xp=tg_log_xp.get(ix_best)+(rng.random()*0.5-0.25);
        ret2=eos_vary_dist(nB,Ye_best,T,log_xn,log_xp,Zbar,Nbar,
                           thx,mun_full,mup_full,
                           A_min_best,A_max_best,NmZ_min_best,
                           NmZ_max_best,vdet,
                           dist_changed,no_nuclei);
      }
      if (ret2!=0) {
        cout << "Point failed." << endl;
        exit(-1);
      }

      // Compute the number density of free neutrons and protons
      double n_fraction, p_fraction;
      if (nB<0.16) {
        double xn=0.0;
        if (log_xn>-300.0) {
          xn=pow(10.0,log_xn);
        }
        double xp=0.0;
        if (log_xp>-300.0) {
          xp=pow(10.0,log_xp);
        }
        double n0_loc=0.16;
        neutron.n=xn*(1.0-nB/n0_loc)/(1.0-nB*xn/n0_loc-nB*xp/n0_loc)*nB;
        proton.n=xp*(1.0-nB/n0_loc)/(1.0-nB*xn/n0_loc-nB*xp/n0_loc)*nB;    
      } else {
        neutron.n=(1.0-Ye_best)*nB;
        proton.n=Ye_best*nB;
      }

      double mun_gas=vdet["mun_gas"];
      double mup_gas=vdet["mup_gas"];
      cout << "mun_gas [MeV], mup_gas [MeV]: " << mun_gas*hc_mev_fm << " "
           << mup_gas*hc_mev_fm << endl;
      
      // Make sure to compute kf, which is not always computed at
      // finite temperature
      sk.def_fet.kf_from_density(neutron);
      sk.def_fet.kf_from_density(proton);
      
      // Set the electron chemical potential (with the electron rest mass)
      double mue=mun_gas+neutron.m-mup_gas-proton.m;
      
      // Copy the electron results to the local electron object
      electron.n=nB*Ye_best;
      electron.mu=mue;
      cout << "Ye_best, mue (with rest mass) [MeV]: " << Ye_best << " "
           << mue*hc_mev_fm << endl;
      
      ubvector X;
      compute_X(nB,X);
    
      cout << "Beta-eq point (ret2,Ye_best): " << ret2 << " "
           << Ye_best << " " << "\n  log_xn,log_xp: " 
           << log_xn << " " << log_xp
           << "\n  Zbar,Nbar: " << Zbar << " " << Nbar
           << "\n  Abar,Yebar,Xnuclei:"
           << Zbar+Nbar << " " << Zbar/(Zbar+Nbar) << " " << X[5] << endl;

      double mu_n_nonint, mu_p_nonint;
      if (true) {
        // Noninteracting fermions, but with the same mass as the
        // effective mass of the original neutron and proton
        fermion n2(vdet["msn"]/hc_mev_fm,2.0),
          p2(vdet["msp"]/hc_mev_fm,2.0);
        n2.n=neutron.n;
        p2.n=proton.n;
        n2.mu=mun_gas;
        p2.mu=mup_gas;
        n2.inc_rest_mass=false;
        p2.inc_rest_mass=false;
        fermion_nonrel fnr;
        fnr.calc_density(n2,T);
        fnr.calc_density(p2,T);
        mu_n_nonint=n2.mu;
        mu_p_nonint=p2.mu;
      }
      
      if (ret2==0) {

        cout << "mun: " << neutron.mu*hc_mev_fm << endl;
        cout << "mup: " << proton.mu*hc_mev_fm << endl;
        cout << "msn: " << vdet["msn"] << " "
             << vdet_units.find("msn")->second << endl;
        cout << "msp: " << vdet["msp"] << " "
             << vdet_units.find("msp")->second << endl;
        cout << "nn: " << neutron.n << endl;
        cout << "np: " << proton.n << endl;
        cout << "g,dgdnn [fm^3],dgdnp [fm^3]: " << vdet["g"] << " "
             << vdet["dgdnn"] << " " << vdet["dgdnp"] << endl;

        double u2eos=neutron.mu*hc_mev_fm-mu_n_nonint*hc_mev_fm;
        cout << "U2 [MeV]: " << u2eos << endl;
        double u4eos=proton.mu*hc_mev_fm-mu_p_nonint*hc_mev_fm;
        cout << "U4 [MeV]: " << u4eos << endl;
        cout << "T [MeV]: " << T*hc_mev_fm << endl;
        
        if (false) {
          vdet["msn"]=neutron.m*hc_mev_fm;
          vdet["msp"]=proton.m*hc_mev_fm;
          neutron.n=0.721726*0.0002;
          proton.n=0.0002-neutron.n;
          cout << "neutron.n proton.n: ";
          cout << neutron.n << " " << proton.n << endl;
          u2eos=-0.230804;
          u4eos=-0.392108;
          neutron.mu=-46.6625/hc_mev_fm;
          proton.mu=-56.375/hc_mev_fm;
          electron.mu=9.7124/hc_mev_fm;
          electron.n=proton.n;
        }
        
        FluidState betaEoS;
        betaEoS=FluidState::StateFromDensities
          (T*hc_mev_fm,vdet["msn"],vdet["msp"],
           neutron.n*pow(hc_mev_fm,3.0),proton.n*pow(hc_mev_fm,3.0),
           u2eos,u4eos,electron.m*hc_mev_fm,electron.n*pow(hc_mev_fm,3.0));
      
        WeakCouplings nscat=WeakCouplings::NeutronScattering();
        nscat.F2=0.0;
      
        WeakCouplings ncap=WeakCouplings::NuCapture();
        ncap.F2=0.0;
      
        // Incoming neutrino energy
        double E1=30.0;
      
        betaEoS.Mu2=neutron.mu*hc_mev_fm;
        betaEoS.Mu4=proton.mu*hc_mev_fm;
        betaEoS.Mu3=(electron.mu-electron.m)*hc_mev_fm;
        cout << "mu2 [MeV], mu4 [MeV], mu3 [MeV] (without rest mass): "
             << betaEoS.Mu2 << " "
             << betaEoS.Mu4 << " "
             << betaEoS.Mu3 << endl;

        //PolarizationNonRel pol_cc(betaEoS, ncap, false, false, true);
        //pol_cc.current=Polarization::current_charged;
        
        PolarizationNonRel pol_cc(betaEoS, ncap, false, false, true);
        pol_cc.current=Polarization::current_charged;
        
        PolarizationNonRel pol_nc(betaEoS, nscat, false, false, false);
        pol_nc.current=Polarization::current_neutral;

        if (n_point>5 && ipoint>50) {
          pol_nc.qagiu.tol_abs=4.0e-19;
          pol_cc.qagiu.tol_abs=4.0e-19;
        } else {
          pol_nc.qagiu.tol_abs=1.0e-10;
          pol_cc.qagiu.tol_abs=1.0e-10;
        }
          
        // [fm^2]
        double fnn_sk=0.5*(sk.t0*(1.0-sk.x0)+
                           1.0/6.0*sk.t3*pow((neutron.n+proton.n),sk.alpha)*
                           (1.0-sk.x3)+2.0/3.0*sk.alpha*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha-1)*
                           ((1+sk.x3/2.0)*(neutron.n+proton.n)-
                            (1.0/2.0+sk.x3)*neutron.n)+1.0/6.0*
                           sk.alpha*(sk.alpha-1.0)*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha-2.0)*
                           ((1+sk.x3/2.0)*pow((neutron.n+proton.n),2.0)-
                            (0.5+sk.x3)*
                            (neutron.n*neutron.n+proton.n*proton.n)))+
          0.25*(sk.t1*(1-sk.x1)+3*sk.t2*(1+sk.x2))*neutron.kf*neutron.kf;
        
        // [fm^2]
        double w1nn_vec_sk=(sk.t0*(1.0-sk.x0)+
                            1.0/6.0*sk.t3*pow((neutron.n+proton.n),
                                              sk.alpha)*    
                            (1.0-sk.x3)+2.0/3.0*sk.alpha*sk.t3*
                            pow((neutron.n+proton.n),sk.alpha-1)* 
                            ((1+sk.x3/2.0)*(neutron.n+proton.n)-           
                             (1.0/2.0+sk.x3)*neutron.n)+1.0/6.0*            
                            sk.alpha*(sk.alpha-1.0)*sk.t3*                   
                            pow((neutron.n+proton.n),sk.alpha-2.0)*          
                            ((1+sk.x3/2.0)*pow((neutron.n+proton.n),2.0)-    
                             (0.5+sk.x3)*                                     
                             (neutron.n*neutron.n+proton.n*proton.n)));
        
        // [fm^4]
        double w2nn_vec_sk=0.25*(sk.t1*(1-sk.x1)+3*sk.t2*(1+sk.x2));
        
        // [fm^2]
        double fpp_sk=0.5*(sk.t0*(1.0-sk.x0)+
                           1.0/6.0*sk.t3*pow((neutron.n+proton.n),sk.alpha)*
                           (1.0-sk.x3)+2.0/3.0*sk.alpha*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha-1)*
                           ((1+sk.x3/2.0)*(neutron.n+proton.n)-
                            (1.0/2.0+sk.x3)*proton.n)+1.0/6.0*sk.alpha*
                           (sk.alpha-1.0)*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha-2.0)*
                           ((1+sk.x3/2.0)*pow((neutron.n+proton.n),2.0)-
                            (0.5+sk.x3)*
                            (neutron.n*neutron.n+proton.n*proton.n)))+
          0.25*(sk.t1*(1-sk.x1)+3*sk.t2*(1+sk.x2))*proton.kf*proton.kf;
        
        // [fm^2]
        double gnn_sk=0.5*(sk.t0*(sk.x0-1)+
                           1.0/6.0*sk.t3*pow((neutron.n+proton.n),sk.alpha)*
                           (sk.x3-1.0))+
          0.25*(sk.t1*(sk.x1-1)+sk.t2*(1+sk.x2))*neutron.kf*neutron.kf;

        // [fm^2]
        double w1nn_ax_sk=(sk.t0*(sk.x0-1)+
                           1.0/6.0*sk.t3*pow((neutron.n+proton.n),sk.alpha)*
                           (sk.x3-1.0));
        
        // [fm^4]
        double w2nn_ax_sk=0.25*(sk.t1*(sk.x1-1)+sk.t2*(1+sk.x2));       
        
        // [fm^2]
        double gpp_sk=0.5*(sk.t0*(sk.x0-1)+
                           1.0/6.0*sk.t3*pow((neutron.n+proton.n),sk.alpha)*
                           (sk.x3-1.0))+
          0.25*(sk.t1*(sk.x1-1)+sk.t2*(1+sk.x2))*proton.kf*proton.kf;
        
        // [fm^2]
        double fnp_sk=0.5*(sk.t0*(2.0+sk.x0)+1.0/6.0*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha)*(2.0+sk.x3)+
                           1.0/2.0*sk.alpha*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha)+
                           1.0/6.0*sk.alpha*(sk.alpha-1.0)*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha-2.0)*
                           ((1+sk.x3/2.0)*pow((neutron.n+proton.n),2.0)-
                            (0.5+sk.x3)*
                            (neutron.n*neutron.n+proton.n*proton.n)))+
          0.5*0.25*(sk.t1*(2.0+sk.x1)+sk.t2*(2.0+sk.x2))*
          (neutron.kf*neutron.kf+proton.kf*proton.kf);

        // [fm^2]
        double w1np_vec_sk=(sk.t0*(2.0+sk.x0)+1.0/6.0*sk.t3*
                            pow((neutron.n+proton.n),sk.alpha)*(2.0+sk.x3)+
                            1.0/2.0*sk.alpha*sk.t3*
                            pow((neutron.n+proton.n),sk.alpha)+
                            1.0/6.0*sk.alpha*(sk.alpha-1.0)*sk.t3*
                            pow((neutron.n+proton.n),sk.alpha-2.0)*
                            ((1+sk.x3/2.0)*pow((neutron.n+proton.n),2.0)-
                             (0.5+sk.x3)*
                             (neutron.n*neutron.n+proton.n*proton.n)));
        
        // [fm^4]
        double w2np_vec_sk=0.25*(sk.t1*(2.0+sk.x1)+sk.t2*(2.0+sk.x2));
        
        // [fm^2]
        double gnp_sk=0.5*(sk.t0*sk.x0+1.0/6.0*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha)*sk.x3)+
          0.5*0.25*(sk.t1*sk.x1+sk.t2*sk.x2)*
          (neutron.kf*neutron.kf+proton.kf*proton.kf);
        
        // [fm^2]
        double w1np_ax_sk=(sk.t0*sk.x0+1.0/6.0*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha)*sk.x3);
        
        // [fm^4]
        double w2np_ax_sk=0.25*(sk.t1*sk.x1+sk.t2*sk.x2);
      
        // Convert these to 1/MeV^2 by dividing by (hbar*c)^2
        fnn_sk/=pow(hc_mev_fm,2);
        fnp_sk/=pow(hc_mev_fm,2);
        fpp_sk/=pow(hc_mev_fm,2);
        gnn_sk/=pow(hc_mev_fm,2);
        gnp_sk/=pow(hc_mev_fm,2);
        gpp_sk/=pow(hc_mev_fm,2);
        w1nn_vec_sk/=pow(hc_mev_fm,2);
        w1nn_ax_sk/=pow(hc_mev_fm,2);
        w1np_vec_sk/=pow(hc_mev_fm,2);
        w1np_ax_sk/=pow(hc_mev_fm,2);
          
        ecv.include_deuteron=true;
        double b_n=ecv.bn_f(T*hc_mev_fm);
        double b_pn=ecv.bpn_f(T*hc_mev_fm);
      
        // [1/MeV]
        double lambda=sqrt(4.0*o2scl_const::pi/(neutron.m+proton.m)/T/
                           hc_mev_fm/hc_mev_fm);
        
        // [1/MeV^3]
        double lambda3=lambda*lambda*lambda;
      
        // [1/MeV^2]
        double f0=ecv.f0(lambda,T*hc_mev_fm);
        double f0p=ecv.f0p(lambda,T*hc_mev_fm);
        double g0=ecv.g0(lambda,T*hc_mev_fm);
        double g0p=ecv.g0p(lambda,T*hc_mev_fm);
      
        cout << "lambda: " << lambda << endl;
        cout << "bn0,bn0free,bn1,bn1free: "
             << ecv.bn0(T*hc_mev_fm) << " " << ecv.bn0_free() << " "
             << ecv.bn1(T*hc_mev_fm) << " " << ecv.bn1_free() << endl;
        cout << "bpn0,bpnfree,bpn1,bpn1free: "
             << ecv.bpn0(T*hc_mev_fm) << " " << ecv.bpn0_free() << " "
             << ecv.bpn1(T*hc_mev_fm) << " " << ecv.bpn1_free() << endl;
        cout << "f0,f0p,g0,g0p: " << f0 << " " << f0p << " " << g0 << " "
             << g0p << endl;
        ecv.include_deuteron=false;

        // [1/MeV^2]
        double fnn_virial=f0+f0p;
        double fnp_virial=f0-f0p;
        double fpp_virial=fnn_virial;
        double gnn_virial=g0+g0p;
        double gnp_virial=g0-g0p;
        double gpp_virial=gnn_virial;
      
        double g_virial=vdet["g"];

        double bn0_hat=ecv.bn0(T_MeV)-ecv.bn0_free();
        double bn1_hat=ecv.bn1(T_MeV)-ecv.bn1_free();
        double bpn0_hat=ecv.bpn0(T_MeV)-ecv.bpn0_free();
        double bpn1_hat=ecv.bpn1(T_MeV)-ecv.bpn1_free();

        // Both in [MeV]
        double dUdnn_vir=(-bpn0_hat*T_MeV*lambda3*proton.n-
                          bpn1_hat*T_MeV*lambda3*proton.n-
                          bn0_hat*T_MeV*lambda3*neutron.n-
                          bn1_hat*T_MeV*lambda3*neutron.n)*pow(hc_mev_fm,3.0);
        double dUdnp_vir=(-bpn0_hat*T_MeV*lambda3*neutron.n-
                          bpn1_hat*T_MeV*lambda3*neutron.n-
                          bn0_hat*T_MeV*lambda3*proton.n-
                          bn1_hat*T_MeV*lambda3*proton.n)*pow(hc_mev_fm,3.0);
        
        // Both in [1/MeV]
        double dtau_dtaun_vir=1.0/neutron.m/hc_mev_fm;
        double dtau_dtaup_vir=1.0/proton.m/hc_mev_fm;
        
        // Both in [MeV]
        double dUdnn_sk=((neutron.n+proton.n)*sk.t0*(1.0+sk.t0/2.0)-
                         neutron.n*sk.t0*(0.5+sk.t0)-
                         neutron.n*pow(neutron.n+proton.n,sk.alpha)/6.0*
                         sk.t3*(0.5+sk.x3)-
                         pow(neutron.n+proton.n,-1.0+sk.alpha)/12.0*
                         sk.t3*(0.5+sk.x3)*
                         sk.alpha*(neutron.n*neutron.n+proton.n*proton.n)+
                         pow(neutron.n+proton.n,1.0+sk.alpha)/12.0*
                         sk.t3*(0.5+sk.x3)*
                         (2.0+sk.alpha))*hc_mev_fm;
        double dUdnp_sk=((neutron.n+proton.n)*sk.t0*(1.0+sk.t0/2.0)-
                         proton.n*sk.t0*(0.5+sk.t0)-
                         proton.n*pow(neutron.n+proton.n,sk.alpha)/
                         6.0*sk.t3*(0.5+sk.x3)-
                         pow(neutron.n+proton.n,-1.0+sk.alpha)/
                         12.0*sk.t3*(0.5+sk.x3)*
                         sk.alpha*(neutron.n*neutron.n+proton.n*proton.n)+
                         pow(neutron.n+proton.n,1.0+sk.alpha)/
                         12.0*sk.t3*(0.5+sk.x3)*
                         (2.0+sk.alpha))*hc_mev_fm;
        
        // Both in [1/MeV]
        double dtau_dtaun_sk=(1.0/neutron.m+
                              2.0*(0.25*(neutron.n+proton.n)*
                                   (sk.t1*(1.0+sk.x1/2.0)+
                                    sk.t2*(1.0+sk.x2/2.0))+
                                   0.25*neutron.n*
                                   (-sk.t1*(0.5+sk.x1)+
                                    sk.t2*(0.5+sk.x2))))/hc_mev_fm;
        double dtau_dtaup_sk=(1.0/proton.m+
                              2.0*(0.25*(neutron.n+proton.n)*
                                   (sk.t1*(1.0+sk.x1/2.0)+
                                    sk.t2*(1.0+sk.x2/2.0))+
                                   0.25*proton.n*
                                   (-sk.t1*(0.5+sk.x1)+
                                    sk.t2*(0.5+sk.x2))))/hc_mev_fm;
        
        // Coulomb correction for fpp
        double e2=1.0/137.0*4.0*o2scl_const::pi;
        
        // qtf2 has units of MeV^2
        double qtf2=4.0*e2*cbrt(o2scl_const::pi)*
          pow(3.0*proton.n*pow(hc_mev_fm,3),2.0/3.0);
        double q=3.0*T_MeV;
        
        // Variable coulombf has units of 1/MeV^2
        double coulombf=e2*4.0*o2scl_const::pi/(q*q+qtf2);

        // Convert dgdnn and dgdnp to [1/MeV^3]
        vdet["dgdnn"]/=pow(hc_mev_fm,3.0);
        vdet["dgdnp"]/=pow(hc_mev_fm,3.0);

        // [1/MeV^2]
        double fnn=fnn_virial*g_virial+fnn_sk*(1.0-g_virial)+
          2.0*vdet["dgdnn"]*dUdnn_vir-2.0*vdet["dgdnn"]*dUdnn_sk+
          (-vdet["dgdnn"]*dtau_dtaun_sk+vdet["dgdnn"]*dtau_dtaun_vir)*
          neutron.kf*neutron.kf*pow(hc_mev_fm,2.0);

        // "dg0" means, terms without dgdnn and dgdnp terms
        
        // [1/MeV^2]
        double fnn_dg0=fnn_virial*g_virial+fnn_sk*(1.0-g_virial);
        
        // [1/MeV^2]
        double w1nn_vec_general=2.0*fnn_virial*g_virial+
          (1.0-g_virial)*w1nn_vec_sk/pow(hc_mev_fm,2.0)+
          2.0*(2.0*vdet["dgdnn"]*dUdnn_vir-2.0*vdet["dgdnn"]*dUdnn_sk);
        
        // [fm^2/MeV^2]
        double w2nn_vec_general=(1.0-g_virial)*
          w2nn_vec_sk/pow(hc_mev_fm,2.0)+
          (-vdet["dgdnn"]*dtau_dtaun_sk+vdet["dgdnn"]*
           dtau_dtaun_vir)*pow(hc_mev_fm,2.0);

        // [1/MeV^2]
        double fnp=fnp_virial*g_virial+fnp_sk*(1.0-g_virial)+
          vdet["dgdnn"]*(dUdnp_vir-dUdnp_sk)+
          vdet["dgdnp"]*(dUdnn_vir-dUdnn_sk)+
          0.5*(vdet["dgdnn"]*(dtau_dtaup_vir-dtau_dtaup_sk)*
               proton.kf*proton.kf+
               vdet["dgdnp"]*(dtau_dtaun_vir-dtau_dtaun_sk)*
               neutron.kf*neutron.kf)*pow(hc_mev_fm,2.0);

        // [1/MeV^2]
        double fnp_dg0=fnp_virial*g_virial+fnp_sk*(1.0-g_virial);
        
        // [1/MeV^2]
        double w1np_vec_general=2.0*fnp_virial*g_virial+
          (1.0-g_virial)*w1np_vec_sk/pow(hc_mev_fm,2.0)+
          2.0*(vdet["dgdnn"]*(dUdnp_vir-dUdnp_sk)+
               vdet["dgdnp"]*(dUdnn_vir-dUdnn_sk));
        
        // [fm^2/MeV^2]
        double w2np_vec_general=(1.0-g_virial)*w2np_vec_sk/pow(hc_mev_fm,2.0)+
          vdet["dgdnn"]*(dtau_dtaup_vir-dtau_dtaup_sk)*pow(hc_mev_fm,2.0)*
          proton.kf*proton.kf/(proton.kf*proton.kf+neutron.kf*neutron.kf)+
          vdet["dgdnp"]*(dtau_dtaun_vir-dtau_dtaun_sk)*pow(hc_mev_fm,2.0)*
          neutron.kf*neutron.kf/(proton.kf*proton.kf+neutron.kf*neutron.kf);
 
        // [1/MeV^2]
        double fpp=fpp_virial*g_virial+fpp_sk*(1.0-g_virial)+
          2.0*vdet["dgdnp"]*dUdnp_vir-2.0*vdet["dgdnp"]*dUdnp_sk+
          (-vdet["dgdnp"]*dtau_dtaun_sk+vdet["dgdnp"]*dtau_dtaun_vir)*
          proton.kf*proton.kf*pow(hc_mev_fm,2.0);
        
        // [1/MeV^2]
        double fpp_dg0=fpp_virial*g_virial+fpp_sk*(1.0-g_virial);
        
        // [1/MeV^2]
        double gnn=gnn_virial*g_virial+gnn_sk*(1.0-g_virial);

        // [1/MeV^2]
        double w1nn_ax_general=2.0*gnn_virial*g_virial+
          (1.0-g_virial)*w1nn_ax_sk;

        // [fm^2/MeV^2]
        double w2nn_ax_general=(1.0-g_virial)*w2nn_ax_sk/pow(hc_mev_fm,2.0);
 
        // [1/MeV^2]
        double gnp=gnp_virial*g_virial+gnp_sk*(1.0-g_virial);

        // [1/MeV^2]
        double w1np_ax_general=2.0*gnp_virial*g_virial+
          (1.0-g_virial)*w1np_ax_sk;
        
        // [fm^2/MeV^2]
        double w2np_ax_general=(1.0-g_virial)*w2np_ax_sk/pow(hc_mev_fm,2.0);
        
        // [1/MeV^2]
        double gpp=gpp_virial*g_virial+gpp_sk*(1.0-g_virial);

        cout << "fnn [1/MeV^2], fnn_dg0 [1/MeV^2], fnp [1/MeV^2], "
             << "fnp_dg0 [1/MeV^2], fpp [1/MeV^2], fpp_dg0 [1/MeV^2]: "
             << fnn << " " << fnn_dg0 << " " << fnp << " "
             << fnp_dg0 << " " << fpp << " " << fpp_dg0 << endl;

        cout << "gnn [1/MeV^2], gnp [1/MeV^2], gpp [1/MeV^2]: "
             << gnn << " " << gnp << " " << gpp << endl;
      
        // kf should be the hole momenta at fermi see surface, here the
        // transition is (pn^-1,pn^-1), the hole is neutron hole

        // Rearrangement terms

        // Units of 1/MeV^2
        double rea=(sk.t3*sk.alpha/3.0*pow(neutron.n+proton.n,sk.alpha-1.0)*
                    ((neutron.n+proton.n)*(1.0+sk.x3/2.0)-
                     neutron.n*(sk.x3+0.5)))/pow(hc_mev_fm,2.0);
        
        // Units of MeV
        double reb=(sk.t3*sk.alpha/12.0*pow(neutron.n+proton.n,sk.alpha-1.0)*
                    (pow(neutron.n+proton.n,2.0)*(1.0+sk.x3/2.0)-
                     (neutron.n*neutron.n+proton.n*proton.n)*(sk.x3+0.5)))*
          hc_mev_fm;
        
        // Units of 1/MeV^2
        double rec=(0.25*sk.alpha*sk.t3*pow(neutron.n+proton.n,sk.alpha))/
          pow(hc_mev_fm,2.0);
        
        // Units of 1/MeV^2
        double fnn_tilde=fnn-(1.0-g_virial)*rea+vdet["dgdnn"]*reb*2.0;
        
        // Units of 1/MeV^2
        double w1nn_vec_sktilde=w1nn_vec_sk-2.0*rea;
        
        // Units of 1/MeV^2
        double w1nn_vec_general_tilde=2.0*fnn_virial*g_virial+
          (1.0-g_virial)*w1nn_vec_sktilde+
          2.0*(2.0*vdet["dgdnn"]*dUdnn_vir-2.0*vdet["dgdnn"]*dUdnn_sk);

        // Units of 1/MeV^2
        double fnp_tilde=fnp-(1.0-g_virial)*rec+
          (vdet["dgdnn"]+vdet["dgdnp"])*reb;

        // Units of 1/MeV^2
        double w1np_vec_sktilde=w1np_vec_sk-2.0*rec;
        
        // Units of 1/MeV^2
        double w1np_vec_general_tilde=2.0*fnp_virial*g_virial+
          (1.0-g_virial)*w1np_vec_sktilde+
          2.0*vdet["dgdnn"]*(dUdnp_vir-dUdnp_sk)+
          vdet["dgdnp"]*(dUdnn_vir-dUdnn_sk);
        
        // Units of 1/MeV^2
        double vf_old=fnn_tilde-fnp_tilde;

        // Units of 1/MeV^2
        double vf=0.5*(w1nn_vec_general_tilde-w1np_vec_general_tilde)+
          (vdet["dgdnn"]*reb*2.0-(vdet["dgdnn"]+vdet["dgdnp"])*reb)+
          (w2nn_vec_general-w2np_vec_general)*neutron.kf*neutron.kf;
        
        // Units of 1/MeV^2
        double w1nn_vec_general_tilde_dg0=2.0*fnn_virial*g_virial+
          (1.0-g_virial)*w1nn_vec_sktilde;
        
        // Units of 1/MeV^2
        double w1np_vec_general_tilde_dg0=2.0*fnp_virial*g_virial+
          (1.0-g_virial)*w1np_vec_sktilde;
        
        // Units of fm^2/MeV^2
        double w2nn_vec_general_dg0=(1.0-g_virial)*
          w2nn_vec_sk/pow(hc_mev_fm,2.0);

        // Units of fm^2/MeV^2
        double w2np_vec_general_dg0=(1.0-g_virial)*
          w2np_vec_sk/pow(hc_mev_fm,2.0);

        // Units of 1/MeV^2
        double vf_dg0=0.5*(w1nn_vec_general_tilde_dg0-
                           w1np_vec_general_tilde_dg0)+
          (w2nn_vec_general_dg0-w2np_vec_general_dg0)*neutron.kf*neutron.kf;
        
        // Units of 1/MeV^2
        double vgt_old=gnn-gnp;
        
        // Units of 1/MeV^2
        double vgt=0.5*(w1nn_ax_general-w1np_ax_general)+
          (w2nn_ax_general-w2np_ax_general)*neutron.kf*neutron.kf;

        cout << "vf_old [1/MeV^2], vf [1/MeV^2], "
             << "vf_dg0 [1/MeV^2], vgt_old [1/MeV^2], "
             << "vgt [1/MeV^2]: " << vf_old << " " << vf << " "
             << vf_dg0 << " " << vgt_old << " " << vgt << endl;
 
        // -----------------------------------------------------------------
        // Charged current mean free path

        pol_cc.integ_method_mu=Polarization::integ_mc;
        pol_cc.integ_method_q0=Polarization::integ_mc;
        //pol_cc.integ_method_mu=Polarization::integ_cubature;
        //pol_cc.integ_method_q0=Polarization::integ_cubature;
        //pol_cc.integ_method_mu=Polarization::integ_o2scl;
        //pol_cc.integ_method_q0=Polarization::integ_o2scl;
        //pol_cc.integ_method_mu=Polarization::integ_base;
        //pol_cc.integ_method_q0=Polarization::integ_base;
        
        pol_cc.set_residual(fnn,fnp,fpp,gnn,gnp,gpp,
                            vf_dg0,vgt,proton.n);
        
        pol_cc.flag=Polarization::flag_vector;
        double cc_vec_mfp_dg0=pol_cc.CalculateInverseMFP(E1)/hc_mev_fm*1.e13;
        cout << "charged current, vector part, no dgdn terms: "
             << cc_vec_mfp_dg0 << endl;
      
        pol_cc.set_residual(fnn,fnp,fpp,gnn,gnp,gpp,vf,vgt,proton.n);
        
        pol_cc.flag=Polarization::flag_vector;
        double cc_vec_mfp=pol_cc.CalculateInverseMFP(E1)/hc_mev_fm*1.e13;
        cout << "charged current, vector part: " << cc_vec_mfp << endl;
      
        pol_cc.flag=Polarization::flag_axial;
        double cc_axvec_mfp=pol_cc.CalculateInverseMFP(E1)/hc_mev_fm*1.e13;
        cout << "charged current, axial part: " << cc_axvec_mfp << endl;
      
        // -----------------------------------------------------------------
        // Charged current mean free path without RPA
        
        pol_cc.integ_method_mu=Polarization::integ_mc;
        pol_cc.integ_method_q0=Polarization::integ_mc;
        //pol_cc.integ_method_mu=Polarization::integ_cubature;
        //pol_cc.integ_method_q0=Polarization::integ_cubature;
        //pol_cc.integ_method_mu=Polarization::integ_o2scl;
        //pol_cc.integ_method_q0=Polarization::integ_o2scl;
        //pol_cc.integ_method_mu=Polarization::integ_base;
        //pol_cc.integ_method_q0=Polarization::integ_base;
        
        pol_cc.set_residual(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,proton.n);
        
        pol_cc.flag=Polarization::flag_vector;
        double cc_vec_mfp_norpa=pol_cc.CalculateInverseMFP(E1)/hc_mev_fm*1.e13;
        cout << "charged current, vector part, no RPA: " << cc_vec_mfp << endl;
      
        pol_cc.flag=Polarization::flag_axial;
        double cc_axvec_mfp_norpa=pol_cc.CalculateInverseMFP(E1)/
          hc_mev_fm*1.e13;
        cout << "charged current, axial part, no RPA: " << cc_axvec_mfp << endl;
      
        // -----------------------------------------------------------------
        // Neutral current mean free path
      
        pol_nc.integ_method_mu=Polarization::integ_mc;
        pol_nc.integ_method_q0=Polarization::integ_mc;
        
        pol_nc.set_residual(fnn_dg0,fnp_dg0,fpp_dg0,gnn,gnp,gpp,
                            vf,vgt,proton.n);
      
        pol_nc.flag=Polarization::flag_vector;
        double nc_vec_mfp_dg0=pol_nc.CalculateInverseMFP(E1)/hc_mev_fm*1.e13;
        cout << "neutral current, vector part, no dgdn terms: "
             << nc_vec_mfp_dg0 << endl;
      
        pol_nc.set_residual(fnn,fnp,fpp,gnn,gnp,gpp,vf,vgt,proton.n);
      
        pol_nc.flag=Polarization::flag_vector;
        double nc_vec_mfp=pol_nc.CalculateInverseMFP(E1)/hc_mev_fm*1.e13;
        cout << "neutral current, vector part: " << nc_vec_mfp << endl;
      
        pol_nc.flag=Polarization::flag_axial;
        double nc_axvec_mfp=pol_nc.CalculateInverseMFP(E1)/hc_mev_fm*1.e13;
        cout << "neutral current, axial part: " << nc_axvec_mfp << endl;

        // -----------------------------------------------------------------
        // Neutral current mean free path no RPA
      
        pol_nc.integ_method_mu=Polarization::integ_mc;
        pol_nc.integ_method_q0=Polarization::integ_mc;
        
        pol_nc.set_residual(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,proton.n);
      
        pol_nc.flag=Polarization::flag_vector;
        double nc_vec_mfp_norpa=pol_nc.CalculateInverseMFP(E1)/hc_mev_fm*1.e13;
        cout << "neutral current, vector part, no RPA: " << nc_vec_mfp << endl;
      
        pol_nc.flag=Polarization::flag_axial;
        double nc_axvec_mfp_norpa=pol_nc.CalculateInverseMFP(E1)/
          hc_mev_fm*1.e13;
        cout << "neutral current, axial part, no RPA: " << nc_axvec_mfp << endl;

        line.push_back(nB);
        line.push_back(neutron.n);
        line.push_back(proton.n);
        line.push_back(vdet["g"]);
        line.push_back(vdet["dgdnn"]);
        line.push_back(vdet["dgdnp"]);
        line.push_back(neutron.ms*hc_mev_fm);
        line.push_back(proton.ms*hc_mev_fm);
        line.push_back(neutron.mu*hc_mev_fm);
        line.push_back(proton.mu*hc_mev_fm);
        line.push_back(mu_n_nonint*hc_mev_fm);
        line.push_back(mu_p_nonint*hc_mev_fm);
        line.push_back(electron.mu*hc_mev_fm);
        line.push_back(u2eos);
        line.push_back(u4eos);
        line.push_back(log_xn);
        line.push_back(log_xp);
        line.push_back(Zbar);
        line.push_back(Nbar);
        line.push_back(Zbar+Nbar);
        line.push_back(Zbar/(Zbar+Nbar));

        line.push_back(X[5]);
        line.push_back(Ye_best);
        line.push_back(fnn_sk);
        line.push_back(fpp_sk+coulombf);
        line.push_back(fnp_sk);
        line.push_back(gnn_sk);
        line.push_back(gpp_sk);
        line.push_back(gnp_sk);
        line.push_back(fnn_virial);
        line.push_back(fpp_virial+coulombf);
        line.push_back(fnp_virial);
        line.push_back(gnn_virial);
        line.push_back(gpp_virial);
        line.push_back(gnp_virial);
        line.push_back(fnn);
        line.push_back(fpp+coulombf);
        line.push_back(fnp);
        line.push_back(fnn_dg0);
        line.push_back(fpp_dg0+coulombf);
        line.push_back(fnp_dg0);
        line.push_back(gnn);
        line.push_back(gpp);
        line.push_back(gnp);
        //line.push_back(vf_sk);
        //line.push_back(vgt_sk);
        //line.push_back(vf_virial);
        //line.push_back(vgt_virial);
        line.push_back(vf);
        line.push_back(vf_dg0);
        line.push_back(vgt);
        line.push_back(cc_vec_mfp);
        line.push_back(cc_vec_mfp_dg0);
        line.push_back(cc_axvec_mfp);
        line.push_back(cc_vec_mfp_norpa);
        line.push_back(cc_axvec_mfp_norpa);
        line.push_back(nc_vec_mfp);
        line.push_back(nc_vec_mfp_dg0);
        line.push_back(nc_axvec_mfp);
        line.push_back(nc_vec_mfp_norpa);
        line.push_back(nc_axvec_mfp_norpa);

        // This counting is wrong
        //if (line.size()!=col_list.size()) {
        //cout << line.size() << " " << col_list.size() << endl;
        //O2SCL_ERR("table sync 2.",o2scl::exc_einval);
        //}
        
        // -----------------------------------------------------------------
        // Neutral current dynamic response at q0=w, q=3*T
        
        if (n_point<5) {
          
          //double T_MeV=T*hc_mev_fm;
          
          double vel=sqrt((3.0*T_MeV)/betaEoS.M2);
          double wmin;
          double wmax;
          wmin=-3.0*vel*3*T_MeV-0.00000789;
          // wmin=-50.0+1.0e-3;
          wmax=3.0*(vel*3*T_MeV+3*T_MeV*3*T_MeV/(2.0*betaEoS.M2));
          // wmax=50.0+1.0e-3;
          double dw=(wmax-wmin)/99;
          
          vector<double> w_nc;
          
          for (int k=0;k<100;k++) {
            
            double w=wmin+dw*k;
            w_nc.push_back(w);
            
            Tensor<double> piVV, piAA, piTT, piVA, piVT, piAT;
            double piL, piLn, piLp, piLnRe, piLpRe, piRPAax, piRPAvec;
            
            pol_nc.SetPolarizations_neutral(w,3*T_MeV,&piVV,&piAA,
                                            &piTT,&piVA,&piVT,&piAT,
                                            piLn,piLp,piLnRe,piLpRe,
                                            piRPAvec,piRPAax,piL);
            
            //neutral current
            double zz=(w+0.0)/T_MeV;
            double FermiF=1/(1-exp(-zz));
            
            double resp_RPAvec=2.0*piRPAvec*FermiF;
            double resp_RPAax=2.0*piRPAax*FermiF;
            
            line.push_back(piRPAvec);
            line.push_back(piRPAax);
            line.push_back(resp_RPAvec);
            line.push_back(resp_RPAax);
          }

          if (j==0) {
            hdf_file hf;
            hf.open_or_create(sv[1]);
            hf.setd_vec("w_nc",w_nc);
            hf.close();
          }
          
        }
        
        // -----------------------------------------------------------------
        // Charged current dynamic response, at q0=w, q=3*T
        
        if (n_point<5) {
          
          //double T_MeV=T*hc_mev_fm;
          static const double densFac=pow(hc_mev_fm,3.0);
          
          // Set integration range
          
          double wmin;
          double wmax;
          wmin=-100.0;
          wmax=50.0;
          double dw=(wmax-wmin)/99;
          
          vector<double> w_cc;

          if (false) {
            double w=-90.8929;
            double dcos=2.0/30.0;
            for(double ii=0.0;ii<30.1;ii+=1.0) {
              double mu=-1.0+dcos*ii;
              pol_cc.flag=Polarization::flag_axial;
              double rtue=pol_cc.GetResponse(E1,w,
                                             pol_cc.GetqFromMu13(E1,w,mu));
              cout << ii << " " << mu << " " << rtue << endl;
            }
            exit(-1);
          }
          
          for (int k=0;k<100;k++) {
            
            double w=wmin+dw*k;
            w_cc.push_back(w);
            
            Tensor<double> piVV, piAA, piTT, piVA, piVT, piAT;
            double piL, piLRe, piRPAax, piRPAvec;
            
            pol_cc.SetPolarizations_charged(w,3*T_MeV,&piVV,&piAA,&piTT,
                                            &piVA,&piVT,&piAT,
                                            piLRe,piRPAvec,piRPAax,piL);
            
            double zz=(w+betaEoS.Mu2-betaEoS.Mu4)/T_MeV;
            double FermiF=1/(1-exp(-zz));
            
            double response=pol_cc.GetResponse(E1,w,10*T_MeV);   
            
            double resp_RPAvec=2.0*piRPAvec*FermiF;
            double resp_RPAax=2.0*piRPAax*FermiF;
            
            line.push_back(piRPAvec);
            line.push_back(piRPAax);
            line.push_back(resp_RPAvec);
            line.push_back(resp_RPAax);
          }

          if (j==0) {
            hdf_file hf;
            hf.open_or_create(sv[1]);
            hf.setd_vec("w_cc",w_cc);
            hf.close();
          }
          
        }
      }
    }
    
    if (tab.get_ncolumns()!=line.size()) {
      O2SCL_ERR("Mismatch of columns.",o2scl::exc_einval);
    }
    if (false) {
      for(size_t jj=0;jj<line.size();jj++) {
        cout << tab.get_column_name(jj) << " ";
        cout << line[jj] << endl;
      }
      exit(-1);
    }
    tab.line_of_data(line.size(),line);
    
    if (true || j%100==0) {
      hdf_file hf;
      hf.open_or_create(sv[1]);
      hdf_output(hf,tab,"mb");
      hf.close();
    }

  }

  return 0;
}

int eos_nuclei::mcarlo_neutron(std::vector<std::string> &sv, 
                               bool itive_com) {
  
  if (loaded==false) {
    cerr << "Requires EOS guess." << endl;
    return 2;
  }

  size_t n_point=15;
  if (sv.size()>=3) {
    n_point=stoszt(sv[2]);
  }

  map<string,double> vdet;
  
  table_units<> tab;
  tab.line_of_names("i_ns i_skyrme qmc_alpha qmc_a phi ");
  tab.line_of_names("t0 t1 t2 t3 x0 x1 x2 x3 epsilon ");
  tab.line_of_units(((string)". . MeV MeV . 1/fm^2 1/fm^4 1/fm^4 ")+
                    "1/fm^(3*a+2) . . . . .");
  
  vector<string> col_list={"nB","nn","np",
    "g","dgdnn","dgdnp","msn","msp","mun","mup","mu_n_nonint",
    "mu_p_nonint",
    "U2","U4","log_xn","log_xp",
    "fnn_sk","fpp_sk","fnp_sk",
    "gnn_sk","gpp_sk","gnp_sk",
    "fnn_virial","fpp_virial","fnp_virial",
    "gnn_virial","gpp_virial","gnp_virial",
    "fnn","fpp","fnp","fnn_dg0","fpp_dg0","fnp_dg0",
    "gnn","gpp","gnp",
    "vf","vf_dg0","vgt",
    "nc_vec_imfp","nc_vec_imfp_dg0","nc_axvec_imfp",
    "nc_vec_imfp_norpa","nc_axvec_imfp_norpa"};
  
  vector<string> unit_list={"1/fm^3","1/fm^3","1/fm^3",
    "","1/MeV^3","1/MeV^3","MeV","MeV","MeV","MeV","MeV",
    "MeV",
    "MeV","MeV","","",
    "1/MeV^2","1/MeV^2","1/MeV^2",
    "1/MeV^2","1/MeV^2","1/MeV^2",
    "1/MeV^2","1/MeV^2","1/MeV^2",
    "1/MeV^2","1/MeV^2","1/MeV^2",
    "1/MeV^2","1/MeV^2","1/MeV^2","1/MeV^2","1/MeV^2","1/MeV^2",
    "1/MeV^2","1/MeV^2","1/MeV^2",
    "1/MeV^2","1/MeV^2","1/MeV^2",
    "1/cm","1/cm","1/cm","1/cm","1/cm"};

  if (unit_list.size()!=col_list.size()) {
    cout << col_list.size() << " " << unit_list.size() << endl;
    O2SCL_ERR("Table sync 1 mcarlo_neutron.",o2scl::exc_einval);
  }

  for(size_t ipoint=0;ipoint<n_point;ipoint++) {

    for(size_t ik=0;ik<col_list.size();ik++) {
    std:string temp=col_list[ik]+"_"+o2scl::szttos(ipoint);
      tab.new_column(temp);
      tab.set_unit(temp,unit_list[ik]);
    }
    if (n_point<20) {
      for(size_t ik=0;ik<100;ik++) {
        tab.new_column(((string)"nc_piRPAvec_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
        tab.new_column(((string)"nc_piRPAax_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
        tab.new_column(((string)"nc_piLn_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
        tab.new_column(((string)"nc_piLnRe_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
        tab.new_column(((string)"nc_resp_RPAvec_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
        tab.new_column(((string)"nc_resp_RPAax_")+o2scl::szttos(ik)+"_"+
                       o2scl::szttos(ipoint));
      }
      tab.new_column(((string)"sum_vec_")+o2scl::szttos(ipoint));
      tab.new_column(((string)"sum_ax_")+o2scl::szttos(ipoint));
    }
  }
  
  // 1.0e-4 is well into the virial region, 5.0e-3 gives g \approx 0.6,
  // and 0.15 is near saturation density and far from the virial region
  
  vector<double> nB_list={1.0e-4,5.0e-3,0.016,0.16,
    0.01364,0.01608,0.01860,0.02160,0.02549,0.02947,0.03347,0.03754,
    0.04151,0.04549,0.04952};
  vector<double> TMeV_list={10,10,10,10,
    20,20,20,20,20,20,20,20,
    20,20,20};
  include_detail=true;

  if (n_point>20) {
    nB_list.clear();
    TMeV_list.clear();
    for(size_t j=0;j<100;j++) {
      nB_list.push_back(1.0e-4*pow(0.15/1.0e-4,((double)j)/99.0));
      TMeV_list.push_back(10.0);
    }
  }
  
  static const int N=10000;
  for(int j=0;j<N;j++) {

    std::cout << "j: " << j << endl;

    if (j==0) {
      use_alt_eos=true;
      vector<string> args={"alt-model","Skyrme","NRAPR"};
      alt_model(args,true);
    } else if (j==1) {
      use_alt_eos=true;
      vector<string> args={"alt-model","Skyrme","SGII"};
      alt_model(args,true);
    } else if (j==2) {
      use_alt_eos=true;
      vector<string> args={"alt-model","Skyrme","UNEDF0"};
      alt_model(args,true);
    } else if (j==3) {
      use_alt_eos=true;
      vector<string> args={"alt-model","Skyrme","UNEDF2"};
      alt_model(args,true);
    } else if (j==4) {
      use_alt_eos=true;
      vector<string> args={"alt-model","Skyrme","SV-min"};
      alt_model(args,true);
    } else {
      use_alt_eos=false;
      // Create a random EOS
      std::vector<std::string> obj;
      random(obj,false);
    }
    if (use_alt_eos) {
      // Copy the couplings to the 'sk' object so we can use
      // those for the Fermi Liquid parameters
      sk.t0=sk_alt.t0;
      sk.t1=sk_alt.t1;
      sk.t2=sk_alt.t2;
      sk.t3=sk_alt.t3;
      sk.x0=sk_alt.x0;
      sk.x1=sk_alt.x1;
      sk.x2=sk_alt.x2;
      sk.x3=sk_alt.x3;
      sk.alpha=sk_alt.alpha;
      cout << "t0,t1: " << sk.t0*hc_mev_fm << " " << sk.t1*hc_mev_fm << endl;
      cout << "t2,t3: " << sk.t2*hc_mev_fm << " " << sk.t3*hc_mev_fm << endl;
      cout << "x0,x1: " << sk.x0 << " " << sk.x1 << endl;
      cout << "x2,x3: " << sk.x2 << " " << sk.x3 << endl;
      cout << "alpha: " << sk.alpha << endl;
    }

    vector<double> line={((double)i_ns),((double)i_skyrme),
      qmc_alpha,qmc_a,phi,
      sk.t0*hc_mev_fm,sk.t1*hc_mev_fm,
      sk.t2*hc_mev_fm,sk.t3*hc_mev_fm,
      sk.x0,sk.x1,sk.x2,sk.x3,sk.alpha};

    if (true) {
      hdf_file hf;
      hf.open_or_create(sv[1]);
      hf.setd_vec("nB_list",nB_list);
      hf.close();
    }
    
    for(size_t ipoint=0;ipoint<n_point;ipoint++) {      

      double nB=nB_list[ipoint];
      double T=TMeV_list[ipoint]/hc_mev_fm;
      double T_MeV=TMeV_list[ipoint];

      size_t inB=vector_lookup(n_nB2,nB_grid2,nB);
      //nB=nB_grid2[inB];
      size_t iT=vector_lookup(n_T2,T_grid2,T*hc_mev_fm);
      //T=T_grid2[iT]/hc_mev_fm;
      
      cout << "Using nB = " << nB << " 1/fm^3 and T = " << T*hc_mev_fm
           << " MeV for\n  ipoint = " << ipoint << " out of total "
           << n_point << endl;
      
      double log_xn_best=0.0, log_xp_best=0.0;
      double fr_best=1.0e10;
      size_t iYe_best=0;
      double log_xn, log_xp;
      double Zbar, Nbar;
      thermo thx;
      double mun_full, mup_full;
      bool dist_changed=true;
      bool no_nuclei=false;
    
      eos_sn_base eso;
      eso.include_muons=false;
      thermo lep;
      
      // Now compute the EOS at the optimal Ye
      
      vector<size_t> ix_best={inB,0,iT};

      thermo th_gas;
      neutron.n=nB;
      proton.n=nB/1.0e5;
      double fr=free_energy_density_detail(neutron,proton,T,th_gas,vdet);
      //vdet["g"]=1.0;
      //vdet["dgdnn"]=0.0;
      //vdet["dgdnp"]=0.0;
        /*
        ret2=eos_vary_dist(nB,0.0,T,log_xn,log_xp,Zbar,Nbar,
                           thx,mun_full,mup_full,
                           A_min_best,A_max_best,NmZ_min_best,
                           NmZ_max_best,vdet,
                           dist_changed,no_nuclei);
        */
      
      double mun_gas=neutron.mu;
      double mup_gas=proton.mu;
      cout << "mun_gas [MeV], mup_gas [MeV]: " << mun_gas*hc_mev_fm << " "
           << mup_gas*hc_mev_fm << endl;
      
      // Make sure to compute kf, which is not always computed at
      // finite temperature
      sk.def_fet.kf_from_density(neutron);

      /*
      cout << "Beta-eq point (ret2,Ye_best): " << ret2 << " "
           << Ye_best << " " << "\n  log_xn,log_xp: " 
           << log_xn << " " << log_xp
           << "\n  Zbar,Nbar: " << Zbar << " " << Nbar
           << "\n  Abar,Yebar,Xnuclei:"
           << Zbar+Nbar << " " << Zbar/(Zbar+Nbar) << " " << X[5] << endl;
      */

      double mu_n_nonint, mu_p_nonint;
      if (true) {
        // Noninteracting fermions, but with the same mass as the
        // effective mass of the original neutron and proton
        fermion n2(vdet["msn"]/hc_mev_fm,2.0);
        fermion p2(vdet["msp"]/hc_mev_fm,2.0);
        n2.n=neutron.n;
        p2.n=proton.n;
        n2.mu=mun_gas;
        p2.mu=mup_gas;
        n2.inc_rest_mass=false;
        p2.inc_rest_mass=false;
        fermion_nonrel fnr;
        fnr.calc_density(n2,T);
        fnr.calc_density(p2,T);
        mu_n_nonint=n2.mu;
        mu_p_nonint=p2.mu;
      }
      
      if (true) {

        cout << "nn,np: " << neutron.n << " " << proton.n << endl;
        cout << "mun: " << neutron.mu*hc_mev_fm << endl;
        cout << "mup: " << proton.mu*hc_mev_fm << endl;
        cout << "msn: " << vdet["msn"] << " "
             << vdet_units.find("msn")->second << endl;
        cout << "msp: " << vdet["msp"] << " "
             << vdet_units.find("msp")->second << endl;
        cout << "nn: " << neutron.n << endl;
        cout << "np: " << proton.n << endl;
        cout << "g,dgdnn [fm^3],dgdnp [fm^3]: " << vdet["g"] << " "
             << vdet["dgdnn"] << " " << vdet["dgdnp"] << endl;

        double u2eos=neutron.mu*hc_mev_fm-mu_n_nonint*hc_mev_fm;
        cout << "U2 [MeV]: " << u2eos << endl;
        double u4eos=proton.mu*hc_mev_fm-mu_p_nonint*hc_mev_fm;
        cout << "U4 [MeV]: " << u4eos << endl;
        cout << "T [MeV]: " << T*hc_mev_fm << endl;
        
        if (false) {
          vdet["msn"]=neutron.m*hc_mev_fm;
          vdet["msp"]=proton.m*hc_mev_fm;
          neutron.n=0.721726*0.0002;
          proton.n=0.0002-neutron.n;
          cout << "neutron.n proton.n: ";
          cout << neutron.n << " " << proton.n << endl;
          u2eos=-0.230804;
          u4eos=-0.392108;
          neutron.mu=-46.6625/hc_mev_fm;
          proton.mu=-56.375/hc_mev_fm;
          electron.mu=9.7124/hc_mev_fm;
          electron.n=proton.n;
        }
        
        FluidState betaEoS;
        betaEoS=FluidState::StateFromDensities
          (T*hc_mev_fm,vdet["msn"],vdet["msp"],
           neutron.n*pow(hc_mev_fm,3.0),proton.n*pow(hc_mev_fm,3.0),
           u2eos,u4eos,electron.m*hc_mev_fm,electron.n*pow(hc_mev_fm,3.0));
      
        WeakCouplings nscat=WeakCouplings::NeutronScattering();
        nscat.F2=0.0;
      
        WeakCouplings ncap=WeakCouplings::NuCapture();
        ncap.F2=0.0;
      
        // Incoming neutrino energy
        double E1=30.0;
      
        betaEoS.Mu2=neutron.mu*hc_mev_fm;
        betaEoS.Mu4=proton.mu*hc_mev_fm;
        betaEoS.Mu3=0.0;
        cout << "mu2 [MeV], mu4 [MeV], mu3 [MeV] (without rest mass): "
             << betaEoS.Mu2 << " "
             << betaEoS.Mu4 << " "
             << betaEoS.Mu3 << endl;

        PolarizationNonRel pol_nc(betaEoS,nscat,false,false,false);
        pol_nc.current=Polarization::current_neutral;

        if (n_point>20 && ipoint>50) {
          pol_nc.qagiu.tol_abs=4.0e-19;
        } else {
          pol_nc.qagiu.tol_abs=1.0e-10;
        }
          
        // [fm^2]
        double fnn_sk=0.5*(sk.t0*(1.0-sk.x0)+
                           1.0/6.0*sk.t3*pow((neutron.n+proton.n),sk.alpha)*
                           (1.0-sk.x3)+2.0/3.0*sk.alpha*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha-1)*
                           ((1+sk.x3/2.0)*(neutron.n+proton.n)-
                            (1.0/2.0+sk.x3)*neutron.n)+1.0/6.0*
                           sk.alpha*(sk.alpha-1.0)*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha-2.0)*
                           ((1+sk.x3/2.0)*pow((neutron.n+proton.n),2.0)-
                            (0.5+sk.x3)*
                            (neutron.n*neutron.n+proton.n*proton.n)))+
          0.25*(sk.t1*(1-sk.x1)+3*sk.t2*(1+sk.x2))*neutron.kf*neutron.kf;
        
        // [fm^2]
        double w1nn_vec_sk=(sk.t0*(1.0-sk.x0)+
                            1.0/6.0*sk.t3*pow((neutron.n+proton.n),
                                              sk.alpha)*    
                            (1.0-sk.x3)+2.0/3.0*sk.alpha*sk.t3*
                            pow((neutron.n+proton.n),sk.alpha-1)* 
                            ((1+sk.x3/2.0)*(neutron.n+proton.n)-           
                             (1.0/2.0+sk.x3)*neutron.n)+1.0/6.0*            
                            sk.alpha*(sk.alpha-1.0)*sk.t3*                   
                            pow((neutron.n+proton.n),sk.alpha-2.0)*          
                            ((1+sk.x3/2.0)*pow((neutron.n+proton.n),2.0)-    
                             (0.5+sk.x3)*                                     
                             (neutron.n*neutron.n+proton.n*proton.n)));
        
        // [fm^4]
        double w2nn_vec_sk=0.25*(sk.t1*(1-sk.x1)+3*sk.t2*(1+sk.x2));
        
        // [fm^2]
        double fpp_sk=0.5*(sk.t0*(1.0-sk.x0)+
                           1.0/6.0*sk.t3*pow((neutron.n+proton.n),sk.alpha)*
                           (1.0-sk.x3)+2.0/3.0*sk.alpha*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha-1)*
                           ((1+sk.x3/2.0)*(neutron.n+proton.n)-
                            (1.0/2.0+sk.x3)*proton.n)+1.0/6.0*sk.alpha*
                           (sk.alpha-1.0)*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha-2.0)*
                           ((1+sk.x3/2.0)*pow((neutron.n+proton.n),2.0)-
                            (0.5+sk.x3)*
                            (neutron.n*neutron.n+proton.n*proton.n)))+
          0.25*(sk.t1*(1-sk.x1)+3*sk.t2*(1+sk.x2))*proton.kf*proton.kf;
        
        // [fm^2]
        double gnn_sk=0.5*(sk.t0*(sk.x0-1)+
                           1.0/6.0*sk.t3*pow((neutron.n+proton.n),sk.alpha)*
                           (sk.x3-1.0))+
          0.25*(sk.t1*(sk.x1-1)+sk.t2*(1+sk.x2))*neutron.kf*neutron.kf;

        // [fm^2]
        double w1nn_ax_sk=(sk.t0*(sk.x0-1)+
                           1.0/6.0*sk.t3*pow((neutron.n+proton.n),sk.alpha)*
                           (sk.x3-1.0));
        
        // [fm^4]
        double w2nn_ax_sk=0.25*(sk.t1*(sk.x1-1)+sk.t2*(1+sk.x2));       
        
        // [fm^2]
        double gpp_sk=0.5*(sk.t0*(sk.x0-1)+
                           1.0/6.0*sk.t3*pow((neutron.n+proton.n),sk.alpha)*
                           (sk.x3-1.0))+
          0.25*(sk.t1*(sk.x1-1)+sk.t2*(1+sk.x2))*proton.kf*proton.kf;
        
        // [fm^2]
        double fnp_sk=0.5*(sk.t0*(2.0+sk.x0)+1.0/6.0*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha)*(2.0+sk.x3)+
                           1.0/2.0*sk.alpha*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha)+
                           1.0/6.0*sk.alpha*(sk.alpha-1.0)*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha-2.0)*
                           ((1+sk.x3/2.0)*pow((neutron.n+proton.n),2.0)-
                            (0.5+sk.x3)*
                            (neutron.n*neutron.n+proton.n*proton.n)))+
          0.5*0.25*(sk.t1*(2.0+sk.x1)+sk.t2*(2.0+sk.x2))*
          (neutron.kf*neutron.kf+proton.kf*proton.kf);

        // [fm^2]
        double w1np_vec_sk=(sk.t0*(2.0+sk.x0)+1.0/6.0*sk.t3*
                            pow((neutron.n+proton.n),sk.alpha)*(2.0+sk.x3)+
                            1.0/2.0*sk.alpha*sk.t3*
                            pow((neutron.n+proton.n),sk.alpha)+
                            1.0/6.0*sk.alpha*(sk.alpha-1.0)*sk.t3*
                            pow((neutron.n+proton.n),sk.alpha-2.0)*
                            ((1+sk.x3/2.0)*pow((neutron.n+proton.n),2.0)-
                             (0.5+sk.x3)*
                             (neutron.n*neutron.n+proton.n*proton.n)));
        
        // [fm^4]
        double w2np_vec_sk=0.25*(sk.t1*(2.0+sk.x1)+sk.t2*(2.0+sk.x2));
        
        // [fm^2]
        double gnp_sk=0.5*(sk.t0*sk.x0+1.0/6.0*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha)*sk.x3)+
          0.5*0.25*(sk.t1*sk.x1+sk.t2*sk.x2)*
          (neutron.kf*neutron.kf+proton.kf*proton.kf);
        
        // [fm^2]
        double w1np_ax_sk=(sk.t0*sk.x0+1.0/6.0*sk.t3*
                           pow((neutron.n+proton.n),sk.alpha)*sk.x3);
        
        // [fm^4]
        double w2np_ax_sk=0.25*(sk.t1*sk.x1+sk.t2*sk.x2);
      
        // Convert these to 1/MeV^2 by dividing by (hbar*c)^2
        fnn_sk/=pow(hc_mev_fm,2);
        fnp_sk/=pow(hc_mev_fm,2);
        fpp_sk/=pow(hc_mev_fm,2);
        gnn_sk/=pow(hc_mev_fm,2);
        gnp_sk/=pow(hc_mev_fm,2);
        gpp_sk/=pow(hc_mev_fm,2);
        w1nn_vec_sk/=pow(hc_mev_fm,2);
        w1nn_ax_sk/=pow(hc_mev_fm,2);
        w1np_vec_sk/=pow(hc_mev_fm,2);
        w1np_ax_sk/=pow(hc_mev_fm,2);
          
        ecv.include_deuteron=true;
        double b_n=ecv.bn_f(T*hc_mev_fm);
        double b_pn=ecv.bpn_f(T*hc_mev_fm);
      
        // [1/MeV]
        double lambda=sqrt(4.0*o2scl_const::pi/(neutron.m+proton.m)/T/
                           hc_mev_fm/hc_mev_fm);
        
        // [1/MeV^3]
        double lambda3=lambda*lambda*lambda;
      
        // [1/MeV^2]
        double f0=ecv.f0(lambda,T*hc_mev_fm);
        double f0p=ecv.f0p(lambda,T*hc_mev_fm);
        double g0=ecv.g0(lambda,T*hc_mev_fm);
        double g0p=ecv.g0p(lambda,T*hc_mev_fm);
      
        cout << "lambda: " << lambda << endl;
        cout << "bn0,bn0free,bn1,bn1free: "
             << ecv.bn0(T*hc_mev_fm) << " " << ecv.bn0_free() << " "
             << ecv.bn1(T*hc_mev_fm) << " " << ecv.bn1_free() << endl;
        cout << "bpn0,bpnfree,bpn1,bpn1free: "
             << ecv.bpn0(T*hc_mev_fm) << " " << ecv.bpn0_free() << " "
             << ecv.bpn1(T*hc_mev_fm) << " " << ecv.bpn1_free() << endl;
        cout << "f0,f0p,g0,g0p: " << f0 << " " << f0p << " " << g0 << " "
             << g0p << endl;
        ecv.include_deuteron=false;

        // [1/MeV^2]
        double fnn_virial=f0+f0p;
        double fnp_virial=f0-f0p;
        double fpp_virial=fnn_virial;
        double gnn_virial=g0+g0p;
        double gnp_virial=g0-g0p;
        double gpp_virial=gnn_virial;
      
        double g_virial=vdet["g"];

        double bn0_hat=ecv.bn0(T_MeV)-ecv.bn0_free();
        double bn1_hat=ecv.bn1(T_MeV)-ecv.bn1_free();
        double bpn0_hat=ecv.bpn0(T_MeV)-ecv.bpn0_free();
        double bpn1_hat=ecv.bpn1(T_MeV)-ecv.bpn1_free();

        // Both in [MeV]
        double dUdnn_vir=(-bpn0_hat*T_MeV*lambda3*proton.n-
                          bpn1_hat*T_MeV*lambda3*proton.n-
                          bn0_hat*T_MeV*lambda3*neutron.n-
                          bn1_hat*T_MeV*lambda3*neutron.n)*pow(hc_mev_fm,3.0);
        double dUdnp_vir=(-bpn0_hat*T_MeV*lambda3*neutron.n-
                          bpn1_hat*T_MeV*lambda3*neutron.n-
                          bn0_hat*T_MeV*lambda3*proton.n-
                          bn1_hat*T_MeV*lambda3*proton.n)*pow(hc_mev_fm,3.0);
        
        // Both in [1/MeV]
        double dtau_dtaun_vir=1.0/neutron.m/hc_mev_fm;
        double dtau_dtaup_vir=1.0/proton.m/hc_mev_fm;
        
        // Both in [MeV]
        double dUdnn_sk=((neutron.n+proton.n)*sk.t0*(1.0+sk.t0/2.0)-
                         neutron.n*sk.t0*(0.5+sk.t0)-
                         neutron.n*pow(neutron.n+proton.n,sk.alpha)/6.0*
                         sk.t3*(0.5+sk.x3)-
                         pow(neutron.n+proton.n,-1.0+sk.alpha)/12.0*
                         sk.t3*(0.5+sk.x3)*
                         sk.alpha*(neutron.n*neutron.n+proton.n*proton.n)+
                         pow(neutron.n+proton.n,1.0+sk.alpha)/12.0*
                         sk.t3*(0.5+sk.x3)*
                         (2.0+sk.alpha))*hc_mev_fm;
        double dUdnp_sk=((neutron.n+proton.n)*sk.t0*(1.0+sk.t0/2.0)-
                         proton.n*sk.t0*(0.5+sk.t0)-
                         proton.n*pow(neutron.n+proton.n,sk.alpha)/
                         6.0*sk.t3*(0.5+sk.x3)-
                         pow(neutron.n+proton.n,-1.0+sk.alpha)/
                         12.0*sk.t3*(0.5+sk.x3)*
                         sk.alpha*(neutron.n*neutron.n+proton.n*proton.n)+
                         pow(neutron.n+proton.n,1.0+sk.alpha)/
                         12.0*sk.t3*(0.5+sk.x3)*
                         (2.0+sk.alpha))*hc_mev_fm;
        
        // Both in [1/MeV]
        double dtau_dtaun_sk=(1.0/neutron.m+
                              2.0*(0.25*(neutron.n+proton.n)*
                                   (sk.t1*(1.0+sk.x1/2.0)+
                                    sk.t2*(1.0+sk.x2/2.0))+
                                   0.25*neutron.n*
                                   (-sk.t1*(0.5+sk.x1)+
                                    sk.t2*(0.5+sk.x2))))/hc_mev_fm;
        double dtau_dtaup_sk=(1.0/proton.m+
                              2.0*(0.25*(neutron.n+proton.n)*
                                   (sk.t1*(1.0+sk.x1/2.0)+
                                    sk.t2*(1.0+sk.x2/2.0))+
                                   0.25*proton.n*
                                   (-sk.t1*(0.5+sk.x1)+
                                    sk.t2*(0.5+sk.x2))))/hc_mev_fm;
        
        // Coulomb correction for fpp
        double e2=1.0/137.0*4.0*o2scl_const::pi;
        
        // qtf2 has units of MeV^2
        double qtf2=4.0*e2*cbrt(o2scl_const::pi)*
          pow(3.0*proton.n*pow(hc_mev_fm,3),2.0/3.0);
        double q=3.0*T_MeV;
        
        // Variable coulombf has units of 1/MeV^2
        double coulombf=e2*4.0*o2scl_const::pi/(q*q+qtf2);

        // Convert dgdnn and dgdnp to [1/MeV^3]
        vdet["dgdnn"]/=pow(hc_mev_fm,3.0);
        vdet["dgdnp"]/=pow(hc_mev_fm,3.0);

        // [1/MeV^2]
        double fnn=fnn_virial*g_virial+fnn_sk*(1.0-g_virial)+
          2.0*vdet["dgdnn"]*dUdnn_vir-2.0*vdet["dgdnn"]*dUdnn_sk+
          (-vdet["dgdnn"]*dtau_dtaun_sk+vdet["dgdnn"]*dtau_dtaun_vir)*
          neutron.kf*neutron.kf*pow(hc_mev_fm,2.0);

        // "dg0" means, terms without dgdnn and dgdnp terms
        
        // [1/MeV^2]
        double fnn_dg0=fnn_virial*g_virial+fnn_sk*(1.0-g_virial);
        
        // [1/MeV^2]
        double w1nn_vec_general=2.0*fnn_virial*g_virial+
          (1.0-g_virial)*w1nn_vec_sk/pow(hc_mev_fm,2.0)+
          2.0*(2.0*vdet["dgdnn"]*dUdnn_vir-2.0*vdet["dgdnn"]*dUdnn_sk);
        
        // [fm^2/MeV^2]
        double w2nn_vec_general=(1.0-g_virial)*
          w2nn_vec_sk/pow(hc_mev_fm,2.0)+
          (-vdet["dgdnn"]*dtau_dtaun_sk+vdet["dgdnn"]*
           dtau_dtaun_vir)*pow(hc_mev_fm,2.0);

        // [1/MeV^2]
        double fnp=fnp_virial*g_virial+fnp_sk*(1.0-g_virial)+
          vdet["dgdnn"]*(dUdnp_vir-dUdnp_sk)+
          vdet["dgdnp"]*(dUdnn_vir-dUdnn_sk)+
          0.5*(vdet["dgdnn"]*(dtau_dtaup_vir-dtau_dtaup_sk)*
               proton.kf*proton.kf+
               vdet["dgdnp"]*(dtau_dtaun_vir-dtau_dtaun_sk)*
               neutron.kf*neutron.kf)*pow(hc_mev_fm,2.0);

        // [1/MeV^2]
        double fnp_dg0=fnp_virial*g_virial+fnp_sk*(1.0-g_virial);
        
        // [1/MeV^2]
        double w1np_vec_general=2.0*fnp_virial*g_virial+
          (1.0-g_virial)*w1np_vec_sk/pow(hc_mev_fm,2.0)+
          2.0*(vdet["dgdnn"]*(dUdnp_vir-dUdnp_sk)+
               vdet["dgdnp"]*(dUdnn_vir-dUdnn_sk));
        
        // [fm^2/MeV^2]
        double w2np_vec_general=(1.0-g_virial)*w2np_vec_sk/pow(hc_mev_fm,2.0)+
          vdet["dgdnn"]*(dtau_dtaup_vir-dtau_dtaup_sk)*pow(hc_mev_fm,2.0)*
          proton.kf*proton.kf/(proton.kf*proton.kf+neutron.kf*neutron.kf)+
          vdet["dgdnp"]*(dtau_dtaun_vir-dtau_dtaun_sk)*pow(hc_mev_fm,2.0)*
          neutron.kf*neutron.kf/(proton.kf*proton.kf+neutron.kf*neutron.kf);
 
        // [1/MeV^2]
        double fpp=fpp_virial*g_virial+fpp_sk*(1.0-g_virial)+
          2.0*vdet["dgdnp"]*dUdnp_vir-2.0*vdet["dgdnp"]*dUdnp_sk+
          (-vdet["dgdnp"]*dtau_dtaun_sk+vdet["dgdnp"]*dtau_dtaun_vir)*
          proton.kf*proton.kf*pow(hc_mev_fm,2.0);
        
        // [1/MeV^2]
        double fpp_dg0=fpp_virial*g_virial+fpp_sk*(1.0-g_virial);
        
        // [1/MeV^2]
        double gnn=gnn_virial*g_virial+gnn_sk*(1.0-g_virial);

        // [1/MeV^2]
        double w1nn_ax_general=2.0*gnn_virial*g_virial+
          (1.0-g_virial)*w1nn_ax_sk;

        // [fm^2/MeV^2]
        double w2nn_ax_general=(1.0-g_virial)*w2nn_ax_sk/pow(hc_mev_fm,2.0);
 
        // [1/MeV^2]
        double gnp=gnp_virial*g_virial+gnp_sk*(1.0-g_virial);

        // [1/MeV^2]
        double w1np_ax_general=2.0*gnp_virial*g_virial+
          (1.0-g_virial)*w1np_ax_sk;
        
        // [fm^2/MeV^2]
        double w2np_ax_general=(1.0-g_virial)*w2np_ax_sk/pow(hc_mev_fm,2.0);
        
        // [1/MeV^2]
        double gpp=gpp_virial*g_virial+gpp_sk*(1.0-g_virial);

        cout << "fnn [1/MeV^2], fnn_dg0 [1/MeV^2], fnp [1/MeV^2], "
             << "fnp_dg0 [1/MeV^2], fpp [1/MeV^2], fpp_dg0 [1/MeV^2]: "
             << fnn << " " << fnn_dg0 << " " << fnp << " "
             << fnp_dg0 << " " << fpp << " " << fpp_dg0 << endl;

        cout << "gnn [1/MeV^2], gnp [1/MeV^2], gpp [1/MeV^2]: "
             << gnn << " " << gnp << " " << gpp << endl;
      
        // kf should be the hole momenta at fermi see surface, here the
        // transition is (pn^-1,pn^-1), the hole is neutron hole

        // Rearrangement terms

        // Units of 1/MeV^2
        double rea=(sk.t3*sk.alpha/3.0*pow(neutron.n+proton.n,sk.alpha-1.0)*
                    ((neutron.n+proton.n)*(1.0+sk.x3/2.0)-
                     neutron.n*(sk.x3+0.5)))/pow(hc_mev_fm,2.0);
        
        // Units of MeV
        double reb=(sk.t3*sk.alpha/12.0*pow(neutron.n+proton.n,sk.alpha-1.0)*
                    (pow(neutron.n+proton.n,2.0)*(1.0+sk.x3/2.0)-
                     (neutron.n*neutron.n+proton.n*proton.n)*(sk.x3+0.5)))*
          hc_mev_fm;
        
        // Units of 1/MeV^2
        double rec=(0.25*sk.alpha*sk.t3*pow(neutron.n+proton.n,sk.alpha))/
          pow(hc_mev_fm,2.0);
        
        // Units of 1/MeV^2
        double fnn_tilde=fnn-(1.0-g_virial)*rea+vdet["dgdnn"]*reb*2.0;
        
        // Units of 1/MeV^2
        double w1nn_vec_sktilde=w1nn_vec_sk-2.0*rea;
        
        // Units of 1/MeV^2
        double w1nn_vec_general_tilde=2.0*fnn_virial*g_virial+
          (1.0-g_virial)*w1nn_vec_sktilde+
          2.0*(2.0*vdet["dgdnn"]*dUdnn_vir-2.0*vdet["dgdnn"]*dUdnn_sk);

        // Units of 1/MeV^2
        double fnp_tilde=fnp-(1.0-g_virial)*rec+
          (vdet["dgdnn"]+vdet["dgdnp"])*reb;

        // Units of 1/MeV^2
        double w1np_vec_sktilde=w1np_vec_sk-2.0*rec;
        
        // Units of 1/MeV^2
        double w1np_vec_general_tilde=2.0*fnp_virial*g_virial+
          (1.0-g_virial)*w1np_vec_sktilde+
          2.0*vdet["dgdnn"]*(dUdnp_vir-dUdnp_sk)+
          vdet["dgdnp"]*(dUdnn_vir-dUdnn_sk);
        
        // Units of 1/MeV^2
        double vf_old=fnn_tilde-fnp_tilde;

        // Units of 1/MeV^2
        double vf=0.5*(w1nn_vec_general_tilde-w1np_vec_general_tilde)+
          (vdet["dgdnn"]*reb*2.0-(vdet["dgdnn"]+vdet["dgdnp"])*reb)+
          (w2nn_vec_general-w2np_vec_general)*neutron.kf*neutron.kf;
        
        // Units of 1/MeV^2
        double w1nn_vec_general_tilde_dg0=2.0*fnn_virial*g_virial+
          (1.0-g_virial)*w1nn_vec_sktilde;
        
        // Units of 1/MeV^2
        double w1np_vec_general_tilde_dg0=2.0*fnp_virial*g_virial+
          (1.0-g_virial)*w1np_vec_sktilde;
        
        // Units of fm^2/MeV^2
        double w2nn_vec_general_dg0=(1.0-g_virial)*
          w2nn_vec_sk/pow(hc_mev_fm,2.0);

        // Units of fm^2/MeV^2
        double w2np_vec_general_dg0=(1.0-g_virial)*
          w2np_vec_sk/pow(hc_mev_fm,2.0);

        // Units of 1/MeV^2
        double vf_dg0=0.5*(w1nn_vec_general_tilde_dg0-
                           w1np_vec_general_tilde_dg0)+
          (w2nn_vec_general_dg0-w2np_vec_general_dg0)*neutron.kf*neutron.kf;
        
        // Units of 1/MeV^2
        double vgt_old=gnn-gnp;
        
        // Units of 1/MeV^2
        double vgt=0.5*(w1nn_ax_general-w1np_ax_general)+
          (w2nn_ax_general-w2np_ax_general)*neutron.kf*neutron.kf;

        cout << "vf_old [1/MeV^2], vf [1/MeV^2], "
             << "vf_dg0 [1/MeV^2], vgt_old [1/MeV^2], "
             << "vgt [1/MeV^2]: " << vf_old << " " << vf << " "
             << vf_dg0 << " " << vgt_old << " " << vgt << endl;
 
        // -----------------------------------------------------------------
        // Neutral current mean free path
      
        pol_nc.integ_method_mu=Polarization::integ_mc;
        pol_nc.integ_method_q0=Polarization::integ_mc;
        
        pol_nc.set_residual(fnn_dg0,fnp_dg0,fpp_dg0,gnn,gnp,gpp,
                            vf,vgt,proton.n);
      
        pol_nc.flag=Polarization::flag_vector;
        double nc_vec_mfp_dg0=pol_nc.CalculateInverseMFP(E1,true)/
          hc_mev_fm*1.e13;
        cout << "neutral current, vector part, no dgdn terms: "
             << nc_vec_mfp_dg0 << endl;
      
        pol_nc.set_residual(fnn,fnp,fpp,gnn,gnp,gpp,vf,vgt,proton.n);
        
        pol_nc.flag=Polarization::flag_vector;
        double nc_vec_mfp=pol_nc.CalculateInverseMFP(E1,true)/
          hc_mev_fm*1.e13;
        cout << "neutral current, vector part: " << nc_vec_mfp << endl;
      
        pol_nc.flag=Polarization::flag_axial;
        double nc_axvec_mfp=pol_nc.CalculateInverseMFP(E1,true)/
          hc_mev_fm*1.e13;
        cout << "neutral current, axial part: " << nc_axvec_mfp << endl;

        // -----------------------------------------------------------------
        // Neutral current mean free path no RPA
      
        pol_nc.integ_method_mu=Polarization::integ_mc;
        pol_nc.integ_method_q0=Polarization::integ_mc;
        
        pol_nc.set_residual(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,proton.n);
      
        pol_nc.flag=Polarization::flag_vector;
        double nc_vec_mfp_norpa=pol_nc.CalculateInverseMFP(E1,true)/
          hc_mev_fm*1.e13;
        cout << "neutral current, vector part, no RPA: " << nc_vec_mfp << endl;
      
        pol_nc.flag=Polarization::flag_axial;
        double nc_axvec_mfp_norpa=pol_nc.CalculateInverseMFP(E1,true)/
          hc_mev_fm*1.e13;
        cout << "neutral current, axial part, no RPA: " << nc_axvec_mfp << endl;

        // Go back to the "with RPA" calculations to compute the
        // response below
        pol_nc.set_residual(fnn_dg0,fnp_dg0,fpp_dg0,gnn,gnp,gpp,
                            vf,vgt,proton.n);
        
        if (true) {
          line.push_back(nB);
          line.push_back(neutron.n);
          line.push_back(proton.n);
          line.push_back(vdet["g"]);
          line.push_back(vdet["dgdnn"]);
          line.push_back(vdet["dgdnp"]);
          line.push_back(neutron.ms*hc_mev_fm);
          line.push_back(proton.ms*hc_mev_fm);
          line.push_back(neutron.mu*hc_mev_fm);
          line.push_back(proton.mu*hc_mev_fm);
          line.push_back(mu_n_nonint*hc_mev_fm);
          line.push_back(mu_p_nonint*hc_mev_fm);
          line.push_back(u2eos);
          line.push_back(u4eos);
          line.push_back(log_xn);
          line.push_back(log_xp);
          
          //line.push_back(X[5]);
          //line.push_back(Ye_best);
          line.push_back(fnn_sk);
          line.push_back(fpp_sk+coulombf);
          line.push_back(fnp_sk);
          line.push_back(gnn_sk);
          line.push_back(gpp_sk);
          line.push_back(gnp_sk);
          line.push_back(fnn_virial);
          line.push_back(fpp_virial+coulombf);
          line.push_back(fnp_virial);
          line.push_back(gnn_virial);
          line.push_back(gpp_virial);
          line.push_back(gnp_virial);
          line.push_back(fnn);
          line.push_back(fpp+coulombf);
          line.push_back(fnp);
          line.push_back(fnn_dg0);
          line.push_back(fpp_dg0+coulombf);
          line.push_back(fnp_dg0);
          line.push_back(gnn);
          line.push_back(gpp);
          line.push_back(gnp);
          //line.push_back(vf_sk);
          //line.push_back(vgt_sk);
          //line.push_back(vf_virial);
          //line.push_back(vgt_virial);
          line.push_back(vf);
          line.push_back(vf_dg0);
          line.push_back(vgt);
          
          line.push_back(nc_vec_mfp);
          line.push_back(nc_vec_mfp_dg0);
          line.push_back(nc_axvec_mfp);
          line.push_back(nc_vec_mfp_norpa);
          line.push_back(nc_axvec_mfp_norpa);
        }
        
        // This counting is wrong
        //if (line.size()!=col_list.size()) {
        //cout << line.size() << " " << col_list.size() << endl;
        //O2SCL_ERR("table sync 2.",o2scl::exc_einval);
        //}
        
        // -----------------------------------------------------------------
        // Neutral current dynamic response at q0=w, q=3*T
        
        if (n_point<20) {
          
          //double T_MeV=T*hc_mev_fm;
          
          double vel=sqrt((3.0*T_MeV)/betaEoS.M2);
          double wmin;
          double wmax;
          wmin=-3.0*vel*3*T_MeV-0.00000789;
          // wmin=-50.0+1.0e-3;
          wmax=3.0*(vel*3*T_MeV+3*T_MeV*3*T_MeV/(2.0*betaEoS.M2));
          // wmax=50.0+1.0e-3;
          double dw=(wmax-wmin)/99;
          
          vector<double> w_nc;

          double sum_vec=0.0;
          double sum_ax=0.0;
          
          for (int k=0;k<100;k++) {
            
            double w=wmin+dw*k;
            w_nc.push_back(w);
            
            Tensor<double> piVV, piAA, piTT, piVA, piVT, piAT;
            double piL, piLn, piLp, piLnRe, piLpRe, piRPAax, piRPAvec;
            
            pol_nc.SetPolarizations_neutral(w,3*T_MeV,&piVV,&piAA,
                                            &piTT,&piVA,&piVT,&piAT,
                                            piLn,piLp,piLnRe,piLpRe,
                                            piRPAvec,piRPAax,piL,true);
            
            //neutral current
            double zz=(w+0.0)/T_MeV;
            double FermiF=1/(1-exp(-zz));
            
            double resp_RPAvec=2.0*piRPAvec*FermiF;
            double resp_RPAax=2.0*piRPAax*FermiF;

            line.push_back(piRPAvec);
            line.push_back(piRPAax);
            line.push_back(piLn);
            line.push_back(piLnRe);
            line.push_back(resp_RPAvec);
            line.push_back(resp_RPAax);

            sum_vec+=resp_RPAvec*dw;
            sum_ax+=resp_RPAax*dw;
          }

          sum_vec/=2*pi*nB*pow(hc_mev_fm,3.0);
          sum_ax/=2*pi*nB*pow(hc_mev_fm,3.0);
          
          line.push_back(sum_vec);
          line.push_back(sum_ax);
          
          cout << "sum_vec,sum_ax: " << sum_vec << " " << sum_ax << endl;

          if (j==0) {
            hdf_file hf;
            hf.open_or_create(sv[1]);
            hf.setd_vec("w_nc",w_nc);
            hf.close();
          }

        }
      }
    }

    if (true) {
      for(size_t jj=0;jj<line.size() || jj<tab.get_ncolumns();jj++) {
        if (jj<tab.get_ncolumns()) {
          cout << tab.get_column_name(jj) << " ";
        } else {
          cout << "<not present>" << " ";
        }
        if (jj<line.size()) {
          cout << line[jj] << endl;
        } else {
          cout << "<not present>" << endl;
        }
      }
    }
    
    if (tab.get_ncolumns()!=line.size()) {
      O2SCL_ERR("Mismatch of columns.",o2scl::exc_einval);
    }
    
    tab.line_of_data(line.size(),line);
    
    if (true || j%100==0) {
      hdf_file hf;
      hf.open_or_create(sv[1]);
      hdf_output(hf,tab,"mb");
      hf.close();
    }

  }

  return 0;
}

