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
#include "eos_nuclei.h"
#include <yaml-cpp/yaml.h>

#ifdef O2SCL_EIGEN
#include <eigen3/Eigen/Dense>
#endif

#include <o2scl/root_brent_gsl.h>
#include <o2scl/classical.h>
#include <o2scl/interp.h>
#include <o2scl/svd.h>
#include <o2scl/vector.h>
#include <o2scl/mmin_bfgs2.h>
#include <o2scl/diff_evo_adapt.h>
#include <o2scl/rng.h>
#include <o2scl/xml.h>
#include <o2scl/auto_format.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;
using namespace o2scl_auto_format;

eos_nuclei::eos_nuclei() {

  fd_A_max=600;
  
  nuclei.resize(6);
  nuc_alpha=&nuclei[0];
  nuc_deut=&nuclei[1];
  nuc_trit=&nuclei[2];
  nuc_he3=&nuclei[3];
  nuc_li4=&nuclei[4];
  nuc_heavy=&nuclei[5];

  nB_grid_spec="301,10^(i*0.04-12)*2.0";
  Ye_grid_spec="70,0.01*(i+1)";
  T_grid_spec="160,0.1*1.046^i";
  S_grid_spec="5,0.1*i";

  extend_frdm=false;
  show_all_nuclei=false;
  recompute=false;
  verify_only=false;
  edge_list="";
  six_neighbors=0;
  propagate_points=true;
  survey_eqs=false;
  
  mh.ntrial=10000;
  mh.err_nonconv=false;
  mh.def_jac.err_nonconv=false;
  // Note these two lines are different!
  mh.tol_rel=1.0e-6;
  mh_tol_rel=1.0e-6;
  rbg.err_nonconv=false;

  baryons_only=true;
  derivs_computed=false;
  with_leptons=false;

  //nucleon_func="(nb<0.04)*((i<100)*10+(i>=100)*sqrt(i))+(nb>=0.04)*100";
  nucleon_func="(i<100)*10+(i>=100)*sqrt(i)";
  
  vi.append({"msg","index","flag","nB","Ye","T","log_xn","log_xp",
      "Z","N","fr","ed","pr","en","mun","mup","no_nuclei",
      "inB","iYe","iT","A_min","A_max","NmZ_min","NmZ_max",
      "Xalpha","Xd","Xt","XHe3","XLi4","Xnuclei",
      "zn","zp","F1","F2","F3","F4","Un","Up","msn",
      "msp","g","dgdT"});

  select_high_T_internal(6);
  
  ubvector fit_params(10);
  fit_params[0]=7.331653147364e-01;
  fit_params[1]=1.165181857191e+00;
  fit_params[2]=2.633495904779e+01;
  fit_params[3]=1.613488012262e+01;
  fit_params[4]=2.236420424613e+01;
  fit_params[5]=3.265413010129e+01;
  fit_params[6]=2.896993239440e+01;
  fit_params[7]=4.221684396819e-01;
  fit_params[8]=1.624162891236e+00;
  fit_params[9]=7.161943683319e-01;
  frdm.fit_fun(10,fit_params);

  max_time=0.0;
  alg_mode=4;
  fixed_dist_alg=1111;
  rnuc_less_rws=true;

  function_verbose=11111;
  max_ratio=7.0;
  loaded=false;
  file_update_time=1800;
  file_update_iters=100000;

  // These function calls do nothing if these environment variables
  // are not defined
  slack.set_channel_from_env("O2SCL_SLACK_CHANNEL");
  slack.set_username_from_env("O2SCL_SLACK_USERNAME");

  baryons_only=true;

  n_nB2=0;
  n_Ye2=0;
  n_T2=0;
  n_S2=0;

  ext_guess="";
  include_detail=false;

  vdet_units.insert(make_pair("zn",""));
  vdet_units.insert(make_pair("zp",""));
  // Note that vdet has different units than tg_mue!
  vdet_units.insert(make_pair("mue","1/fm"));
  vdet_units.insert(make_pair("Ymu",""));
  vdet_units.insert(make_pair("F1","MeV"));
  vdet_units.insert(make_pair("F2","MeV"));
  vdet_units.insert(make_pair("F3","MeV"));
  vdet_units.insert(make_pair("F4","MeV"));
  vdet_units.insert(make_pair("msn","MeV"));
  vdet_units.insert(make_pair("msp","MeV"));
  vdet_units.insert(make_pair("Un","MeV"));
  vdet_units.insert(make_pair("Up","MeV"));
  vdet_units.insert(make_pair("g",""));
  vdet_units.insert(make_pair("dgdT","1/MeV"));
  vdet_units.insert(make_pair("dgdnn","fm^3"));
  vdet_units.insert(make_pair("dgdnp","fm^3"));
  vdet_units.insert(make_pair("sigma","1/fm"));
  vdet_units.insert(make_pair("omega","1/fm"));
  vdet_units.insert(make_pair("rho","1/fm"));
  vdet_units.insert(make_pair("mun_gas","1/fm"));
  vdet_units.insert(make_pair("mup_gas","1/fm"));
  
  inc_hrg=false;

  pfuncs.spin_deg_mode=1;
}

eos_nuclei::~eos_nuclei() {
}

void eos_nuclei::compute_X(double nB, ubvector &X) {
  X.resize(6);
  
  if (nuclei.size()<6) {
    O2SCL_ERR("Nuclei array not properly sized.",o2scl::exc_esanity);
  }
  
  nuc_alpha=&nuclei[0];
  nuc_deut=&nuclei[1];
  nuc_trit=&nuclei[2];
  nuc_he3=&nuclei[3];
  nuc_li4=&nuclei[4];
  nuc_heavy=&nuclei[5];
  
  X[0]=nuc_alpha->n*4.0/nB;
  X[1]=nuc_deut->n*2.0/nB;
  X[2]=nuc_trit->n*3.0/nB;
  X[3]=nuc_he3->n*3.0/nB;
  X[4]=nuc_li4->n*4.0/nB;
  
  if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
    X[5]=0.0;
    for(size_t i=5;i<nuclei.size();i++) {
      X[5]+=nuclei[i].n*(nuclei[i].Z+nuclei[i].N)/nB;
    }
  } else {
    X[5]=nuc_heavy->n*(nuc_heavy->Z+nuc_heavy->N)/nB;
  }

  return;
}

void eos_nuclei::load_nuclei() {
  int mpi_rank, mpi_size;
#ifndef NO_MPI
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
  
  // Load the nuclear masses
  o2scl_hdf::ame_load(ame);
  o2scl_hdf::mnmsk_load(m95);
  o2scl_hdf::hfb_sp_load(hfb,27);
  pfuncs.load();
  
#ifndef NO_MPI
  // Send a message to the next MPI rank
  if (mpi_size>1 && mpi_rank<mpi_size-1) {
    MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
	     tag,MPI_COMM_WORLD);
  }
#endif

  ame.get_nucleus(2,2,*nuc_alpha);
  ame.get_nucleus(1,1,*nuc_deut);
  ame.get_nucleus(1,2,*nuc_trit);
  ame.get_nucleus(2,1,*nuc_he3);
  ame.get_nucleus(3,1,*nuc_li4);

  nuc_alpha->g=1.0;
  nuc_deut->g=3.0;
  nuc_trit->g=2.0;
  nuc_li4->g=5.0;
  nuc_he3->g=2.0;
 
  // Compute neutron and proton separation energies in MeV. Note that
  // Shen et al. (2010) defines Sn and Sp to be positive for stable
  // nuclei. We only compute the light nuclei here, the heavy nuclei
  // will be added later.
  o2scl::nucleus nuc_temp;
  Sneut.resize(6);
  Sprot.resize(6);
  for(size_t i=0;i<5;i++) {
    if (ame.is_included(nuclei[i].Z,nuclei[i].N-1) &&
        ame.is_included(nuclei[i].Z-1,nuclei[i].N)) {
      ame.get_nucleus(nuclei[i].Z,nuclei[i].N-1,nuc_temp);
      Sneut[i]=-(nuclei[i].be-nuc_temp.be)*hc_mev_fm;
      ame.get_nucleus(nuclei[i].Z-1,nuclei[i].N,nuc_temp);
      Sprot[i]=-(nuclei[i].be-nuc_temp.be)*hc_mev_fm;
    } else {
      Sneut[i]=-1.0;
      Sprot[i]=-1.0;
    }
  }
  
  //for(size_t i=0;i<5;i++) {
  //ame.get_nucleus(nuclei[i].Z,nuclei[i].N-1,nuc_temp);
  //Sneut[i]=-(nuclei[i].be-nuc_temp.be)*hc_mev_fm;
  //ame.get_nucleus(nuclei[i].Z-1,nuclei[i].N,nuc_temp);
  //Sprot[i]=-(nuclei[i].be-nuc_temp.be)*hc_mev_fm;  
  //}
  
  vomega.resize(6);
  vomega_prime.resize(6);
  Ec.resize(6);

  return;
}

double partition_func::delta_small_iand(double E) {
  if (E<1.0e-200) return 0.0;
  // This integrand is in units of 1/MeV? Check this.
  double ret=sqrt(pi)/12.0*exp(2.0*sqrt(a*(E-delta)))/
    pow(a,1.0/4.0)/pow((E-delta),5.0/4.0)*exp(-E/T_MeV);
  if (!std::isfinite(ret)) {
    cout << "a,delta,T_MeV,E: "
	 << a << " " << delta << " " << T_MeV << " " << E << endl;
    cout << exp(2.0*sqrt(a*(E-delta))) << " " << exp(-E/T_MeV) << " " 
	 << pow((E-delta),5.0/4.0) << " " << ret << endl;
    O2SCL_ERR2("Value of delta_small_iand is not finite ",
	       "in partition_func::delta_small_iand().",o2scl::exc_efailed);
  }
  return ret;
}

double partition_func::delta_small_iand_prime(double E) {
  if (E<1.0e-200) return 0.0;
  double ret=E/T_MeV*sqrt(pi)/12.0*exp(2.0*sqrt(a*(E-delta)))/
    pow(a,1.0/4.0)/pow((E-delta),5.0/4.0)*exp(-E/T_MeV);
  if (!std::isfinite(ret)) {
    cout << "a,delta,T_MeV,E: "
	 << a << " " << delta << " " << T_MeV << " " << E << endl;
    O2SCL_ERR2("Value of delta_small_iand_prime is not finite ",
	       "in partition_func::delta_small_iand_prime().",
	       o2scl::exc_efailed);
  }
  return ret;
}
  
double partition_func::delta_large_iand(double E) {
  if (E<1.0e-200) return 0.0;
  double ret=C*exp((E-delta)/Tc)*exp(-E/T_MeV); 
  if (!std::isfinite(ret)) {
    cout << "a,delta,T_MeV,E: "
	 << a << " " << delta << " " << T_MeV << " " << E << endl;
    O2SCL_ERR2("Value of delta_large_iand is not finite ",
	       "in partition_func::delta_large_iand().",o2scl::exc_efailed);
  }
  return ret;
} 

double partition_func::delta_large_iand_prime(double E) {
  if (E<1.0e-200) return 0.0;
  double ret=E/T_MeV*C*exp((E-delta)/Tc)*exp(-E/T_MeV); 
  if (!std::isfinite(ret)) {
    cout << "a,delta,T_MeV,E: "
	 << a << " " << delta << " " << T_MeV << " " << E << endl;
    O2SCL_ERR2("Value of delta_large_iand_prime is not finite ",
	       "in partition_func::delta_large_iand_prime().",
	       o2scl::exc_efailed);
  }
  return ret;
} 

double eos_nuclei::f_min_search(size_t nvar,const ubvector &x,
                                double nb, double ye, double T) {
  
  //we need here nb, ye, T, nuclei[i], vomega[i];
  double n0=0.16;

  /*if(x[0]>1.0 || x[0]<0.0) {
    return 1.0e100;
    }
    if(x[1]>1.0 || x[1]<0.0) {
    return 1.0e100;
    }*/
  if(x[0]>0.0 || x[0]<-150.0) {
    return 1.0e100;
  }
  if(x[1]>0.0 || x[1]<-150.0) {
    return 1.0e100;
  }
  double x0=pow(10.0,x[0]),x1=pow(10.0,x[1]);
  //double x0=x[0],x1=x[1];
  double nnprime1=nb*x0;
  double npprime1=nb*x1;
  neutron.n=nnprime1;
  proton.n=npprime1;
  free_energy_density(neutron,proton,T,th2);
  double mun_old1=neutron.mu,mup_old1=proton.mu;
  double p0_nuc1=th2.pr;
  double fnp=th2.ed-T*th2.en;
  double kappa1=1.0-nb/n0;
  double xi1=kappa1/(1.0-nb*x0/n0-nb*x1/n0);
	    
  double f_total1=0.0,f_c=0.0,sum_nuc=0.0;

  //compute deuteron triton he3 li4
  ubvector Ec1(6),f(6);
  for(size_t i=1;i<5;i++) {
    double A=nuclei[i].A;
    double rA=cbrt(3.0*A/4.0/pi/n0);
    double xx=cbrt(nb*ye/n0*A/nuclei[i].Z);
    Ec1[i]=-0.6*nuclei[i].Z*nuclei[i].Z*fine_structure_f<double>()/rA*
      (3.0/2.0*xx-1.0/2.0*pow(xx,3.0));
    double lambda=sqrt(2.0*pi/nuclei[i].m/T);
    double vv=(nuclei[i].N+nuclei[i].Z)/n0;
    nuclei[i].n=kappa1*vomega[i]/pow(lambda,3.0)*
      exp((((double)(nuclei[i].N))*mun_old1
	   +((double)(nuclei[i].Z))*mup_old1
	   -nuclei[i].be-Ec1[i]-p0_nuc1*vv)/T);

    f[i]=-T*(log(vomega[i]/nuclei[i].n/pow(lambda,3.0))+1.0)
      *nuclei[i].n;
    sum_nuc+=nuclei[i].n;
    f[i]+=nuclei[i].n*nuclei[i].be;
    f_c+=nuclei[i].n*Ec1[i];
    f_total1+=f[i];
  }
            
  double A=nuclei[5].A,Z=nuclei[5].Z;
  double nn1=nnprime1*xi1,np1=npprime1*xi1;
  double nalpha=((nb*ye-nuclei[1].n*nuclei[1].Z-nuclei[2].n
		  *nuclei[2].Z-nuclei[3].n*nuclei[3].Z
		  -nuclei[4].n*nuclei[4].Z)*A
		 -(nb-nuclei[1].n*nuclei[1].A-nuclei[2].n*nuclei[2].A
		   -nuclei[3].n*nuclei[3].A-nuclei[4].n*nuclei[4].A)
		 *Z-np1*A+nn1*Z+np1*Z)/(nuclei[0].Z*A
					-nuclei[0].A*Z);
  double rA_alpha=cbrt(3.0*nuclei[0].A/4.0/pi/n0);
  double xx_alpha=cbrt(nb*ye/n0*nuclei[0].A/nuclei[0].Z);
  double Ecalpha=-0.6*nuclei[0].Z*nuclei[0].Z*fine_structure_f<double>()
    /rA_alpha*(3.0/2.0*xx_alpha
	       -1.0/2.0*pow(xx_alpha,3.0));
  double lambda_alpha=sqrt(2.0*pi/nuclei[0].m/T);

  double nheavy=((nb-nuclei[1].n*nuclei[1].A-nuclei[2].n*nuclei[2].A
		  -nuclei[3].n*nuclei[3].A-nuclei[4].n*nuclei[4].A)
		 -nnprime1*xi1-npprime1*xi1-nalpha*nuclei[0].A)
    /nuclei[5].A;
  if(nheavy<0.0 || nalpha<0.0) {
    return 1.0e100;
  }
  double rA_heavy=cbrt(3.0*nuclei[5].A/4.0/pi/n0);
  double xx_heavy=cbrt(nb*ye/n0*nuclei[5].A/nuclei[5].Z);
  double Echeavy=-0.6*nuclei[5].Z*nuclei[5].Z*fine_structure_f<double>()
    /rA_heavy*(3.0/2.0*xx_heavy
               -1.0/2.0*pow(xx_heavy,3.0));
  double lambda_heavy=sqrt(2.0*pi/nuclei[5].m/T);
	  
  f_total1+=xi1*fnp-T*(log(vomega[5]/nheavy
			   /pow(lambda_heavy,3.0))+1.0)*
    nheavy-T*log(kappa1)*nheavy
    +Echeavy*nheavy+nuclei[5].be*nheavy
    -T*(log(vomega[0]/nalpha
            /pow(lambda_alpha,3.0))+1.0)*nalpha-
    T*log(kappa1)*nalpha
    +Ecalpha*nalpha+nuclei[0].be*nalpha
    +f_c-T*sum_nuc*log(kappa1);
  return f_total1;

}

int eos_nuclei::init_function(size_t dim, const ubvector &x, ubvector &y) {
  for (size_t i=0;i<dim;i++) {
    y[i]=0.1*rng.random();
  }
  return 0;
}

int eos_nuclei::load(std::vector<std::string> &sv,
		     bool itive_com) {
  if (sv.size()<2) {
    cerr << "No filename in load." << endl;
    return 1;
  }
  cout << "Loading: " << sv[1] << endl;
  read_results(sv[1]);
  return 0;
}

int eos_nuclei::output(std::vector<std::string> &sv,
		       bool itive_com) {
  if (sv.size()<2) {
    cerr << "No filename in output." << endl;
    return 1;
  }
  write_results(sv[1]);	
  return 0;
}

int eos_nuclei::max_fun(size_t nv, const ubvector &x, ubvector &y,
                        interp_vec<vector<double>,ubvector> &itp_P, 
                        interp_vec<vector<double>,ubvector> &itp_mu) {
  double nb_low=x[0];
  double nb_high=x[1];

  cout << "I: " << x[0] << " " << x[1] << endl;
  
  if (x[0]<0.01 || x[0]>0.16 || x[1]<0.01 || x[1]>0.16) return 1;
  if (nb_high-nb_low<0.002) return 2;

  double P_low=itp_P.eval(nb_low);
  double P_high=itp_P.eval(nb_high);
  double mu_low=itp_mu.eval(nb_low);
  double mu_high=itp_mu.eval(nb_high);

  y[0]=(P_low-P_high)/fabs(P_high);
  y[1]=(mu_low-mu_high)/fabs(mu_high);
  
  return 0;
}

int eos_nuclei::maxwell(std::vector<std::string> &sv,
                        bool itive_com) {
  
  int method=4;
  
  if (method==4) {
    
    if (with_leptons==false || derivs_computed==false) {
      cout << "maxwell_test only works with leptons and derivatives." << endl;
      return 0;
    }
  
    for (size_t iYe=0;iYe<n_Ye2;iYe++) {
      for (size_t iT=0;iT<n_T2;iT++) {

        if (verbose>1) {
          for(size_t inB=0;inB<n_nB2;inB++) {
            vector<size_t> ix={inB,iYe,iT};
            cout << nB_grid2[inB] << " " << tg_P.get(ix) << endl;
          }
        }

        double Ye=Ye_grid2[iYe];
        double T_MeV=T_grid2[iT];
      
        ubvector negative(n_nB2);
      
        // attempt 2
        for(size_t inB=0;inB<n_nB2;inB++) {
          negative[inB]=0.0;
        }
        
        for(size_t inB=0;inB<n_nB2-1;inB++) {
          vector<size_t> ix={inB,iYe,iT};
          vector<size_t> ixp1={inB+1,iYe,iT};
          if (tg_P.get(ix)>tg_P.get(ixp1)) {
            negative[inB]=1.0;
            negative[inB+1]=1.0;
          }
        }
        // end of attempt 2
      
        double ll=nB_grid2[n_nB2-1];
        double ul=nB_grid2[0];
        size_t li=n_nB2-1;
        size_t ui=0;
        for(size_t inB=0;inB<n_nB2;inB++) {
          if (negative[inB]>0.5) {
            if (nB_grid2[inB]<ll) {
              ll=nB_grid2[inB];
              li=inB;
            }
            if (nB_grid2[inB]>ul) {
              ul=nB_grid2[inB];
              ui=inB;
            }
          }
        }
      
        if (((int)ui)-((int)li)>=10) {
          cout << "Problem." << endl;
          exit(-1);
        }
        if (ui>1 && li<n_nB2-2) {
          //li--;
          //ui++;
          
          size_t left=li-1;
          size_t right=ui+1;
          vector<size_t> ileft={left,iYe,iT};
          vector<size_t> iright={right,iYe,iT};
          
          double nB_left=nB_grid2[left];
          double nB_right=nB_grid2[right];
          double F_left=tg_F.get(ileft);
          double F_right=tg_F.get(iright);
          double P_left=tg_P.get(ileft);
          double P_right=tg_P.get(iright);
          double E_left=tg_E.get(ileft);
          double E_right=tg_E.get(iright);

          while (P_left>P_right) {
            li--;
            ui++;
            
            left=li-1;
            right=ui+1;
            ileft[0]=left;
            iright[0]=right;
          
            nB_left=nB_grid2[left];
            nB_right=nB_grid2[right];
            F_left=tg_F.get(ileft);
            F_right=tg_F.get(iright);
            P_left=tg_P.get(ileft);
            P_right=tg_P.get(iright);
            E_left=tg_E.get(ileft);
            E_right=tg_E.get(iright);
          }
          
          if (true || verbose>1) {
            std::cout << "P_left, P_right: "
                      << P_left << " " << P_right << endl;
          }
          
          cout.precision(4);
          cout << "Adjusting from nB = " 
               << nB_grid2[li] << " to " << nB_grid2[ui]
               << " at Ye = " << Ye_grid2[iYe]
               << " and T = " << T_grid2[iT] << endl;
          
          cout.precision(6);
          
          for(size_t inB=left;inB<=right;inB++) {
            double nb_frac=(nB_grid2[inB]-nB_grid2[left])/
              (nB_grid2[right]-nB_grid2[left]);
            
            vector<size_t> ix={inB,iYe,iT};

            if (verbose>1) {
              cout << inB << " " << nb_frac << " ";
              cout << tg_Pint.get(ix) << " " << tg_P.get(ix) << " ";
            }
             
            double P_new=P_left+(P_right-P_left)*nb_frac;
            tg_Pint.get(ix)+=(P_new-tg_P.get(ix));
            tg_P.get(ix)=P_new;
            double E_new=E_left+(E_right-E_left)*nb_frac;
            tg_Eint.get(ix)+=(E_new-tg_E.get(ix));
            tg_E.get(ix)=E_new;

            if (verbose>1) {
              cout << tg_Pint.get(ix) << " " << tg_P.get(ix) << endl;
            }

            // F = - P + mun * nn + mup * np
            tg_Fint.get(ix)=-tg_Pint.get(ix)/nB_grid2[inB]+
              (1.0-Ye_grid2[iYe])*tg_mun.get(ix)+
              Ye_grid2[iYe]*tg_mup.get(ix);
            
            // E = F + T S 
            tg_Eint.get(ix)=tg_Fint.get(ix)+T_grid2[iT]*tg_Sint.get(ix);
            
            // F = - P + mun * nn + mup * np
            tg_F.get(ix)=-tg_P.get(ix)/nB_grid2[inB]+
              (1.0-Ye_grid2[iYe])*tg_mun.get(ix)+
              Ye_grid2[iYe]*tg_mup.get(ix)+
              Ye_grid2[iYe]*tg_mue.get(ix);
            
            // E = F + T S 
            tg_E.get(ix)=tg_F.get(ix)+T_grid2[iT]*tg_S.get(ix);
            
          }

          bool fail=false;
          for(size_t inB=0;inB<n_nB2-1;inB++) {
            vector<size_t> ix={inB,iYe,iT};
            vector<size_t> ixp1={inB+1,iYe,iT};
            double P1=tg_P.get(ix);
            double P2=tg_P.get(ixp1);
            if (P2<P1) {
              fail=true;
              cout << "Failed between " << nB_grid2[inB] << " and "
                   << nB_grid2[inB+1] << endl;
            }
          }
          
          
          if (verbose>1 || fail) {
            for(size_t inB=0;inB<n_nB2;inB++) {
              vector<size_t> ix={inB,iYe,iT};
              cout << nB_grid2[inB] << " " << tg_P.get(ix) << endl;
            }
            char ch;
            cin >> ch;
          }
          
        } else if (ui>1 || li<n_nB2-2) {
          vector_out(cout,negative,true);
          cout << "iYe,iT: " << iYe << " " << iT << endl;
          cout << "Error (maxwell4): " << li << " " << ui << endl;
          exit(-1);
        }
      
      }
    }
  
  } else if (method==0) {
    
    if (with_leptons==false || derivs_computed==false) {
      cout << "maxwell_test only works with leptons and derivatives." << endl;
      return 0;
    }
  
    for (size_t iYe=0;iYe<n_Ye2;iYe++) {
      for (size_t iT=0;iT<n_T2;iT++) {
      
        double Ye=Ye_grid2[iYe];
        double T_MeV=T_grid2[iT];
      
        ubvector dPdnb(n_nB2);
      
        for(size_t inB=0;inB<n_nB2;inB++) {
          if (inB>0 && inB<n_nB2-1) {
            vector<size_t> ixm1={inB-1,iYe,iT};
            vector<size_t> ixp1={inB+1,iYe,iT};
            dPdnb[inB]=(tg_P.get(ixp1)-tg_P.get(ixm1))/
              (nB_grid2[inB+1]-nB_grid2[inB-1]);
          } else if (inB==0) {
            vector<size_t> ix={inB,iYe,iT};
            vector<size_t> ixp1={inB+1,iYe,iT};
            dPdnb[0]=(tg_P.get(ixp1)-tg_P.get(ix))/
              (nB_grid2[inB+1]-nB_grid2[inB-1]);
          } else {
            dPdnb[n_nB2-1]=dPdnb[n_nB2-2];
          }
        }
      
        ubvector negative(n_nB2);
      
        for(size_t inB=0;inB<n_nB2;inB++) {
          if (dPdnb[inB]<0.0) negative[inB]=1.0;
          else negative[inB]=0.0;
        }
      
        double ll=nB_grid2[n_nB2-1];
        double ul=nB_grid2[0];
        size_t li=n_nB2-1;
        size_t ui=0;
        for(size_t inB=0;inB<n_nB2;inB++) {
          if (negative[inB]>0.5) {
            if (nB_grid2[inB]<ll) {
              ll=nB_grid2[inB];
              li=inB;
            }
            if (nB_grid2[inB]>ul) {
              ul=nB_grid2[inB];
              ui=inB;
            }
          }
        }
      
        if (((int)ui)-((int)li)>=10) {
          cout << "Problem." << endl;
          exit(-1);
        }
        if (ui>1 && li<n_nB2-2) {
          //ui++;
          //li--;
          size_t left=li-1;
          size_t right=ui+1;
          vector<size_t> ileft={left,iYe,iT};
          vector<size_t> iright={right,iYe,iT};
          double nB_left=nB_grid2[left];
          double nB_right=nB_grid2[right];
          double F_left=tg_F.get(ileft)*nB_grid2[left];
          double F_right=tg_F.get(iright)*nB_grid2[right];
          double gamma=log(F_left/F_right)/log(nB_left/nB_right);
          double coeff=F_left/pow(nB_left,gamma);
          cout << Ye << " " << T_MeV << " " << nB_left << " "
               << nB_right << endl;
          if (false) {
            cout << left << " " << right << endl;
            cout << nB_left << " "
                 << (nB_left+nB_right)/2.0 << " " << nB_right << endl;
            cout << F_left << " " << coeff*pow((nB_left+nB_right)/2.0,gamma)
                 << " " << F_right << endl;
          }
          for(size_t inB=li;inB<=ui;inB++) {
            vector<size_t> ix={inB,iYe,iT};
            double shift=coeff*pow(nB_grid2[inB],gamma)-
              tg_F.get(ix)*nB_grid2[inB];
            tg_Fint.get(ix)+=shift/nB_grid2[inB];
            //cout << "inB,shift: " << inB << " " << shift << endl;
          }
          //char ch;
          //cin >> ch;
        } else if (ui>1 || li<n_nB2-2) {
          vector_out(cout,negative,true);
          cout << "Error maxwell 0: " << li << " " << ui << endl;
          exit(-1);
        }
      
      }
    }
  
    with_leptons==false;
    derivs_computed==false;

  } else if (method==1) {
  
    if (with_leptons==false || derivs_computed==false) {
      cout << "maxwell_test only works with leptons and derivatives." << endl;
      return 0;
    }
  
    for (size_t iYe=0;iYe<n_Ye2;iYe++) {
      for (size_t iT=0;iT<n_T2;iT++) {
      
        double Ye=Ye_grid2[iYe];
        double T_MeV=T_grid2[iT];

        cout << "Ye,Y_MeV: " << Ye << " " << T_MeV << endl;

        ubvector P(n_nB2);
        ubvector mu(n_nB2);
        ubvector mun(n_nB2);
        ubvector mup(n_nB2);
        ubvector mue(n_nB2);
        ubvector dPdnb(n_nB2);
      
        for(size_t inB=0;inB<n_nB2;inB++) {
          vector<size_t> ix={inB,iYe,iT};
        
          P[inB]=tg_P.get(ix);
          mun[inB]=tg_mun.get(ix);
          mup[inB]=tg_mup.get(ix);
          mue[inB]=tg_mue.get(ix);
          mu[inB]=(1.0-Ye_grid2[iYe])*tg_mun.get(ix)+
            Ye_grid2[iYe]*(tg_mup.get(ix)+tg_mue.get(ix));
        
          if (inB>0 && inB<n_nB2-1) {
            vector<size_t> ixm1={inB-1,iYe,iT};
            vector<size_t> ixp1={inB+1,iYe,iT};
            dPdnb[inB]=(tg_P.get(ixp1)-tg_P.get(ixm1))/
              (nB_grid2[inB+1]-nB_grid2[inB-1]);
          } else if (inB==0) {
            vector<size_t> ixp1={inB+1,iYe,iT};
            dPdnb[0]=(tg_P.get(ixp1)-tg_P.get(ix))/
              (nB_grid2[inB+1]-nB_grid2[inB-1]);
          } else {
            dPdnb[n_nB2-1]=dPdnb[n_nB2-2];
          }
        }
      
        ubvector negative(n_nB2);
      
        for(size_t inB=0;inB<n_nB2;inB++) {
          if (dPdnb[inB]<0.0) negative[inB]=1.0;
          else negative[inB]=0.0;
        }
      
        double ll=nB_grid2[n_nB2-1];
        double ul=nB_grid2[0];
        size_t li=n_nB2-1;
        size_t ui=0;
        for(size_t inB=0;inB<n_nB2;inB++) {
          if (negative[inB]>0.5) {
            if (nB_grid2[inB]<ll) {
              ll=nB_grid2[inB];
              li=inB;
            }
            if (nB_grid2[inB]>ul) {
              ul=nB_grid2[inB];
              ui=inB;
            }
          }
        }
      
        if (((int)ui)-((int)li)>=10) {
          cout << "Problem." << endl;
          exit(-1);
        }
        if (ui>1 && li<n_nB2-2) {

          double nb_low=nB_grid2[li];
          double nb_high=nB_grid2[ui];
        
          interp_vec<vector<double>,ubvector>
            itp_P(n_nB2,nB_grid2,P,itp_linear);
          interp_vec<vector<double>,ubvector>
            itp_mu(n_nB2,nB_grid2,mu,itp_linear);
          interp_vec<vector<double>,ubvector>
            itp_mun(n_nB2,nB_grid2,mun,itp_linear);
          interp_vec<vector<double>,ubvector>
            itp_mup(n_nB2,nB_grid2,mup,itp_linear);
          interp_vec<vector<double>,ubvector>
            itp_mue(n_nB2,nB_grid2,mue,itp_linear);

          cout << li << " " << ui << " " << nb_low << " " << nb_high << endl;
          ofstream fout;
          fout.open("max2.out");
          fout << "nb P mu" << endl;
          cout << "nb P mu" << endl;
          for(double nbx=0.01;nbx<0.16;nbx+=0.001) {
            fout << nbx << " " << itp_P.eval(nbx) << " "
                 << itp_mu.eval(nbx) << endl;
            cout << nbx << " " << itp_P.eval(nbx) << " "
                 << itp_mu.eval(nbx) << " "
                 << itp_mun.eval(nbx) << " "
                 << itp_mup.eval(nbx) << " "
                 << itp_mue.eval(nbx) << endl;
          }
          fout.close();

          mh.verbose=2;

          mm_funct func=std::bind
            (std::mem_fn<int(size_t,const ubvector &,ubvector&,
                             interp_vec<vector<double>,ubvector> &,
                             interp_vec<vector<double>,ubvector> &)>
             (&eos_nuclei::max_fun),this,std::placeholders::_1,
             std::placeholders::_2,std::placeholders::_3,
             std::ref(itp_P),std::ref(itp_mun));

          ubvector x(2), y(2);
          x[0]=nb_low;
          x[1]=nb_high;

          mh.msolve(2,x,func);

          cout << "H: " << x[0] << " " << x[1] << endl;
        
          cout << "Done." << endl;
          exit(-1);
        
        } else if (ui>1 || li<n_nB2-2) {
          vector_out(cout,negative,true);
          cout << "Error: maxwell 1 " << li << " " << ui << endl;
          exit(-1);
        }
      
      }
    }
  
    with_leptons==false;
    derivs_computed==false;

  } else if (method==2) {
  
    if (with_leptons==false || derivs_computed==false) {
      cout << "maxwell_test only works with leptons and derivatives." << endl;
      return 0;
    }
  
    for (size_t iYe=0;iYe<n_Ye2;iYe++) {
      for (size_t iT=0;iT<n_T2;iT++) {
      
        double Ye=Ye_grid2[iYe];
        double T_MeV=T_grid2[iT];

        cout << "Ye,Y_MeV: " << Ye << " " << T_MeV << endl;

        ubvector P(n_nB2);
        ubvector mu(n_nB2);
        ubvector oonB(n_nB2);
      
        for(size_t inB=1;inB<n_nB2;inB++) {
          vector<size_t> ixm1={inB-1,iYe,iT};
          vector<size_t> ix={inB,iYe,iT};
          double Pint1=tg_Pint.get(ixm1);
          double Pint2=tg_Pint.get(ix);

          oonB[inB]=1.0/nB_grid2[inB];
          mu[inB]=(1.0-Ye_grid2[iYe])*tg_mun.get(ix)+
            Ye_grid2[iYe]*(tg_mup.get(ix)+tg_mue.get(ix));
        }

        vector<double> signc;
        vector<size_t> signix;
        vector<double> signP;
        ofstream fout;
        fout.open("max3.out");
        for(size_t inB=1;inB<n_nB2-1;inB++) {
          cout << inB << " " << P[inB]-P[inB-1] << " " << P[inB] << " "
               << mu[inB] << endl;
          fout << nB_grid2[inB] << " " << P[inB] << endl;
          if ((P[inB]-P[inB-1])*(P[inB+1]-P[inB])<0.0) {
            signc.push_back(nB_grid2[inB]);
            signix.push_back(inB);
            signP.push_back(P[inB]);
          }
        }
        fout.close();

        if (signc.size()>0) {
          vector_out(cout,signix,true);
          vector_out(cout,signc,true);

          size_t inB_min=vector_min_value<vector<size_t>,size_t>(signix);
          size_t inB_max=vector_max_value<vector<size_t>,size_t>(signix);
          signP.push_back(P[inB_min-1]);
          signP.push_back(P[inB_max+1]);
          vector_out(cout,signP,true);

          double P_min=vector_min_value<vector<double>,double>(signP);
          double P_max=vector_max_value<vector<double>,double>(signP);
          cout << P_min << " " << P_max << endl;

          interp_vec<ubvector,ubvector> itp_P(n_nB2,oonB,P,itp_linear);

          for(double Pt=P_min;Pt<P_max*1.0000001;Pt+=(P_max-P_min)/60.0) {

            vector<double> locs;
            vector_find_level(Pt,n_nB2,oonB,P,locs);
          
            cout << Pt << " " << locs.size() << endl;

            if (locs.size()==3) {
              cout << itp_P.integ(locs[0],locs[1]) << " "
                   << itp_P.integ(locs[1],locs[2]) << " "
                   << itp_P.integ(locs[0],locs[2]) << endl;
           
            }
          
          }
        
          exit(-1);
        }
      
      }
    }
  
    with_leptons==false;
    derivs_computed==false;

  }
  
  return 0;
}
  
int eos_nuclei::new_muons(size_t nv, const ubvector &x, ubvector &y,
                          double nB, double Ye, double T,
                          map<string,double> &vdet, eos_sn_base &eso) {
  
  eso.electron.n=x[0];
  
  thermo lep;
  eso.compute_eg_point(nB,eso.electron.n/nB,T*hc_mev_fm,lep,vdet["mue"]);
  
  y[0]=(eso.electron.n+eso.muon.n-nB*Ye)/(nB*Ye);
  //cout << "\t" << eso.electron.n << " " << eso.muon.n << " " << nB*Ye << endl;

  vdet["Ymu"]=eso.muon.n/nB;
  
  return 0;
}


int eos_nuclei::add_eg(std::vector<std::string> &sv,
		       bool itive_com) {

  if (loaded==false || n_nB2==0 || n_Ye2==0 || n_T2==0) {
    cerr << "No EOS loaded." << endl;
    return 2;
  }
  
  if (derivs_computed==false) {
    cerr << "add_eg requires derivs_computed==true." << endl;
    return 3;
  }

  size_t st[3]={n_nB2,n_Ye2,n_T2};
  vector<double> packed;
  for(size_t i=0;i<n_nB2;i++) {
    packed.push_back(nB_grid2[i]);
  }
  for(size_t i=0;i<n_Ye2;i++) {
    packed.push_back(Ye_grid2[i]);
  }
  for(size_t i=0;i<n_T2;i++) {
    packed.push_back(T_grid2[i]);
  }

  if (tg_E.total_size()==0) {
    tg_E.resize(3,st);
    tg_E.set_grid_packed(packed);
    tg_E.set_all(0.0);
  }
  if (tg_P.total_size()==0) {
    tg_P.resize(3,st);
    tg_P.set_grid_packed(packed);
    tg_P.set_all(0.0);
  }
  if (tg_S.total_size()==0) {
    tg_S.resize(3,st);
    tg_S.set_grid_packed(packed);
    tg_S.set_all(0.0);
  }
  if (tg_F.total_size()==0) {
    tg_F.resize(3,st);
    tg_F.set_grid_packed(packed);
    tg_F.set_all(0.0);
  }
  if (tg_mue.total_size()==0) {
    tg_mue.resize(3,st);
    tg_mue.set_grid_packed(packed);
    tg_mue.set_all(0.0);
  }
  if (include_muons && tg_Ymu.total_size()==0) {
    tg_Ymu.resize(3,st);
    tg_Ymu.set_grid_packed(packed);
    tg_Ymu.set_all(0.0);
  }
  
  elep.include_muons=include_muons;

  map<std::string,double> vdet;
  
  for (size_t i=0;i<n_nB2;i++) {
    for (size_t j=0;j<n_Ye2;j++) {
      double mue_last;
      for (size_t k=0;k<n_T2;k++) {
	
	thermo lep;
        double nB=nB_grid2[i];
        double Ye=Ye_grid2[j];
        double T_MeV=T_grid2[k];
        double T=T_grid2[k]/hc_mev_fm;
        
        // Note that this function accepts the temperature in MeV
        // and the electron chemical potential is returned in 1/fm.
        if (k==0) {
          // For the lowest temperature grid point, give an initial
          // guess for the electron chemical potential
          vdet["mue"]=0.511/197.33;
        } else {
          vdet["mue"]=mue_last;
        }

        if (include_muons) {
          
          cout << nB << " " << Ye << " " << T_MeV << endl;

          elep.pair_density_eq(nB*Ye,T_MeV/hc_mev_fm);
          if (verbose>1) {
            cout << nB << " " << Ye << " " << T_MeV << " "
                 << elep.e.n << " " << elep.mu.n << " "
                 << elep.e.n+elep.mu.n << endl;
          }

        } else {
        
          elep.pair_density_eq(nB*Ye,T_MeV/hc_mev_fm);
          if (verbose>1) {
            cout << "ed+pr: " << elep.th.ed+elep.th.pr 
                 << "en*T+n*mu: "
                 << elep.th.en*T_MeV/hc_mev_fm+elep.e.mu*elep.e.n
                 << endl;
          }
          
        }

        vdet["mue"]=elep.e.mu;
        mue_last=vdet["mue"];
        
        vector<size_t> ix={i,j,k};
	tg_F.set(ix,tg_Fint.get(ix)+
                 (hc_mev_fm*elep.th.ed-T_MeV*elep.th.en)/nB);
	tg_E.set(ix,tg_Eint.get(ix)+hc_mev_fm*elep.th.ed/nB);
	tg_P.set(ix,tg_Pint.get(ix)+hc_mev_fm*elep.th.pr);
	tg_S.set(ix,tg_Sint.get(ix)+elep.th.en/nB);
	tg_mue.set(ix,hc_mev_fm*vdet["mue"]);
        
        double np=nB*Ye;
        double nn=nB*(1.0-Ye);
        
        if (include_muons) {
          // Set muon density
          tg_Ymu.set(ix,elep.mu.n/nB);
        }

        if (true) {
          
          double scale=fabs(tg_E.get(ix));
          if (scale<10.0) scale=10.0;
          
          double ti_check=tg_E.get(ix)*nB+tg_P.get(ix)-
            T_MeV*tg_S.get(ix)*nB-nn*tg_mun.get(ix)-np*tg_mup.get(ix);
          if (this->include_muons) {
            ti_check-=(elep.e.n+elep.mu.n)*tg_mue.get(ix);
          } else {
            ti_check-=Ye*nB*tg_mue.get(ix);
          }            
          ti_check/=scale*nB;
          
          if (fabs(ti_check)>1.0e-4) {
            
            double ti_check2;
            ti_check2=(elep.th.ed*hc_mev_fm+elep.th.pr*hc_mev_fm-
                       T_MeV*elep.th.en-
                       Ye*nB*vdet["mue"]*hc_mev_fm)/scale/nB;
            double ti_int_check=(tg_Eint.get(ix)*nB+
                                 tg_Pint.get(ix)-
                                 T_MeV*tg_Sint.get(ix)*nB-
                                 nn*tg_mun.get(ix)-
                                 np*tg_mup.get(ix))/scale/nB;
            
            cout << "ti check failed." << endl;
            cout << "nB,Ye,T_MeV,ti_check,scale: " << nB << " " << Ye << " "
                 << T_MeV << " " << ti_check << " " << scale << endl;
            cout << "ti_int_check, ti_check2, Ye*nB: "
                 << ti_int_check << " " << ti_check2 << " " << Ye*nB << endl;
            cout << "terms in ti_check: " << tg_E.get(ix)*nB << " "
                 << tg_P.get(ix) << " "
                 << T_MeV*tg_S.get(ix)*nB << " "
                 << nn*tg_mun.get(ix) << " "
                 << np*tg_mup.get(ix) << " "
                 << np*tg_mue.get(ix) << endl;
            cout << "terms in ti_check2: "
                 << elep.th.ed << " " << elep.th.pr << " "
                 << T_MeV/hc_mev_fm*elep.th.en << " "
                 << np*vdet["mue"] << endl;
            cout << "  " << elep.th.ed+elep.th.pr << " "
                 << T_MeV/hc_mev_fm*elep.th.en+np*vdet["mue"] << endl;
            cout << "  " << elep.th.ed+elep.th.pr-
              T_MeV/hc_mev_fm*elep.th.en-np*vdet["mue"] << endl;
            cout << "  "
                 << (elep.th.ed+elep.th.pr-T_MeV/hc_mev_fm*
                     elep.th.en-np*vdet["mue"])/scale/
              nB*hc_mev_fm << endl;
            cout << "mue [MeV]: " << tg_mue.get(ix) << endl;
            exit(-1);
          }
          //char ch;
          //cin >> ch;
        }          

      }
      
    }
    cout << "add_eg(): " << i+1 << "/" << n_nB2 << endl;
  }

  with_leptons=true;

  return 0;
}

int eos_nuclei::eg_table(std::vector<std::string> &sv,
			 bool itive_com) {

  if (sv.size()<2) {
    cerr << "Need filename." << endl;
    return 2;
  }
  
  // This is similar to eos::eg_table(), so I need
  // to figure out how to combine them sensibly.

  if (loaded==false || n_nB2==0 || n_Ye2==0 || n_T2==0) {
    new_table();
  }
  
  o2scl::tensor_grid<> mu_e;
  o2scl::tensor_grid<> E;
  o2scl::tensor_grid<> Ymu;
  o2scl::tensor_grid<> P;
  o2scl::tensor_grid<> S;
  o2scl::tensor_grid<> F;
  
  size_t st[3]={n_nB2,n_Ye2,n_T2};
  vector<vector<double> > grid={nB_grid2,Ye_grid2,T_grid2};

  mu_e.resize(3,st);
  Ymu.resize(3,st);
  E.resize(3,st);
  P.resize(3,st);
  S.resize(3,st);
  F.resize(3,st);

  mu_e.set_grid(grid);
  Ymu.set_grid(grid);
  E.set_grid(grid);
  P.set_grid(grid);
  S.set_grid(grid);
  F.set_grid(grid);

  eos_sn_base esb;
  esb.include_muons=include_muons;

  for (size_t i=0;i<n_nB2;i++) {
    double nB=nB_grid2[i];
    for (size_t j=0;j<n_Ye2;j++) {
      for (size_t k=0;k<n_T2;k++) {
        vector<size_t> ix={i,j,k};
        
	double T_MeV=T_grid2[k];
	
	thermo lep;
	double mue;
	esb.compute_eg_point(nB_grid2[i],Ye_grid2[j],T_grid2[k],lep,mue);

	mu_e.set(ix,mue*hc_mev_fm);
	E.set(ix,lep.ed/nB*hc_mev_fm);
	P.set(ix,lep.pr*hc_mev_fm);
	S.set(ix,lep.en/nB);
	F.set(ix,(lep.ed-T_MeV*lep.en)/nB*hc_mev_fm);
	if (include_muons) {
	  Ymu.set(ix,muon.n*hc_mev_fm/nB);
	}
      }
    }
    cout << "eg_table(): " << i+1 << "/" << n_nB2 << endl;
  }

  hdf_file hf;
  hf.open_or_create(sv[1]);
  hdf_output(hf,mu_e,"mue");
  if (include_muons) {
    hdf_output(hf,Ymu,"Ymu");
  }
  hdf_output(hf,F,"F");
  hdf_output(hf,E,"E");
  hdf_output(hf,P,"P");
  hdf_output(hf,S,"S");
  hf.close();
  
  return 0;
}

int eos_nuclei::eos_deriv(std::vector<std::string> &sv,
                          bool itive_com) {

  std::cout << "Computing derivatives." << endl;
  
  // -----------------------------------------------------
  // Read table
    
  if (derivs_computed==false) {
    
    derivs_computed=true;
    
    size_t st[3]={n_nB2,n_Ye2,n_T2};
    
    calc_utf8<> calc;
    std::map<std::string,double> vars;
    
    vector<double> packed;
    vector<std::string> split_res;

    split_string_delim(nB_grid_spec,split_res,',');
    n_nB2=stoszt(split_res[0]);
    
    calc.compile(split_res[1].c_str());
    for(size_t i=0;i<n_nB2;i++) {
      vars["i"]=((double)i);
      packed.push_back(nB_grid2[i]);
    }
    
    split_string_delim(Ye_grid_spec,split_res,',');
    n_Ye2=stoszt(split_res[0]);
    
    calc.compile(split_res[1].c_str());
    for(size_t i=0;i<n_Ye2;i++) {
      vars["i"]=((double)i);
      packed.push_back(Ye_grid2[i]);
    }
    
    split_string_delim(T_grid_spec,split_res,',');
    n_T2=stoszt(split_res[0]);
    
    calc.compile(split_res[1].c_str());
    for(size_t i=0;i<n_T2;i++) {
      vars["i"]=((double)i);
      packed.push_back(T_grid2[i]);
    }
    
    tg_Eint.resize(3,st);
    tg_Pint.resize(3,st);
    tg_Sint.resize(3,st);
    tg_mun.resize(3,st);
    tg_mup.resize(3,st);
    
    tg_Eint.set_grid_packed(packed);
    tg_Pint.set_grid_packed(packed);
    tg_Sint.set_grid_packed(packed);
    tg_mun.set_grid_packed(packed);
    tg_mup.set_grid_packed(packed);
    
    tg_Eint.set_all(0.0);
    tg_Pint.set_all(0.0);
    tg_Sint.set_all(0.0);
    tg_mun.set_all(0.0);
    tg_mup.set_all(0.0);
  }

  // Temporary vectors which store the free energy density in
  // MeV/fm^3
  ubvector fint_vec_nB(n_nB2), fint_vec_Ye(n_Ye2), fint_vec_T(n_T2);

  // The interpolator
  interp_steffen<vector<double>,ubvector> itp_stf;

  // Temporary vectors for derivatives in the nB, Ye, and T directions
  fint_vec_nB.resize(n_nB2);
  fint_vec_Ye.resize(n_Ye2);
  fint_vec_T.resize(n_T2);

  // First, compute the temperature derivative, which is just the
  // entropy times minus one
  for (size_t i=0;i<n_nB2;i++) {
    for (size_t j=0;j<n_Ye2;j++) {
      for(size_t k=0;k<n_T2;k++) {
        vector<size_t> ix={i,j,k};
	fint_vec_T[k]=tg_Fint.get(ix)*nB_grid2[i];
      }
      itp_stf.set(n_T2,T_grid2,fint_vec_T);
      for(size_t k=0;k<n_T2;k++) {
        vector<size_t> ix={i,j,k};
	// Note that fint_vec_T above is stored in MeV/fm^3, so
	// when we take a derivative wrt to temperature (stored
	// in MeV), we get the correct units of 1/fm^3, and then
	// we divide by the baryon density in units of 1/fm^3
	tg_Sint.set(ix,-itp_stf.deriv(T_grid2[k])/nB_grid2[i]);
      }
    }
  }

  // Second, compute the electron fraction derivative, which we
  // temporarily store in tg_mup
  for (size_t i=0;i<n_nB2;i++) {
    for(size_t k=0;k<n_T2;k++) {
      for (size_t j=0;j<n_Ye2;j++) {
        vector<size_t> ix={i,j,k};
	fint_vec_Ye[j]=tg_Fint.get(ix)*nB_grid2[i];
      }
      itp_stf.set(n_Ye2,Ye_grid2,fint_vec_Ye);
      for (size_t j=0;j<n_Ye2;j++) {
        vector<size_t> ix={i,j,k};
	tg_mup.set(ix,itp_stf.deriv(Ye_grid2[j]));
      }
    }
  }

  // Third, compute the baryon density derivative, which we
  // temporarily store in tg_mun
  for (size_t j=0;j<n_Ye2;j++) {
    for(size_t k=0;k<n_T2;k++) {
      for (size_t i=0;i<n_nB2;i++) {
        vector<size_t> ix={i,j,k};
	fint_vec_nB[i]=tg_Fint.get(ix)*nB_grid2[i];
      }
      itp_stf.set(n_nB2,nB_grid2,fint_vec_nB);
      for (size_t i=0;i<n_nB2;i++) {
        vector<size_t> ix={i,j,k};
	tg_mun.set(ix,itp_stf.deriv(nB_grid2[i]));
      }
    }
  }

  // Now go through every point and compute the remaining
  // quantities
  for (size_t i=0;i<n_nB2;i++) {
    for (size_t j=0;j<n_Ye2;j++) {
      for (size_t k=0;k<n_T2;k++) {
        vector<size_t> ix={i,j,k};
	
	double en=tg_Sint.get(ix)*nB_grid2[i];
	double dfdnB=tg_mun.get(ix);
	double dfdYe=tg_mup.get(ix);

	if (false && nB_grid2[i]>0.3) {
	  cout << tg_mun.get(ix) << endl;
	  cout << tg_mup.get(ix) << endl;
	}
	
	tg_mun.get(ix)=dfdnB-dfdYe*Ye_grid2[j]/nB_grid2[i];
	tg_mup.get(ix)=dfdnB-dfdYe*(Ye_grid2[j]-1.0)/nB_grid2[i];

	if (false && nB_grid2[i]>0.3) {
	  cout << tg_mun.get(ix) << endl;
	  cout << tg_mup.get(ix) << endl;
	  exit(-1);
	}
	
	// E = F + T S
	tg_Eint.get(ix)=tg_Fint.get(ix)+
	  T_grid2[k]*tg_Sint.get(ix);
	
	// P = - F + mun * nn + mup * np
	tg_Pint.get(ix)=-tg_Fint.get(ix)*nB_grid2[i]+
	  nB_grid2[i]*(1.0-Ye_grid2[j])*tg_mun.get(ix)+
	  nB_grid2[i]*Ye_grid2[j]*tg_mup.get(ix);
      }
    }
  }
  
  std::cout << "Finished computing derivatives." << endl;

  return 0;
    
}

int eos_nuclei::eos_deriv_v2(std::vector<std::string> &sv,
                             bool itive_com) {
  
  std::cout << "Computing derivatives." << endl;
  
  // -----------------------------------------------------
  // Read table
    
  if (derivs_computed==false) {
    
    derivs_computed=true;
    
    size_t st[3]={n_nB2,n_Ye2,n_T2};
    
    calc_utf8<> calc;
    std::map<std::string,double> vars;
    
    vector<double> packed;
    vector<std::string> split_res;

    split_string_delim(nB_grid_spec,split_res,',');
    n_nB2=stoszt(split_res[0]);
    
    calc.compile(split_res[1].c_str());
    for(size_t i=0;i<n_nB2;i++) {
      vars["i"]=((double)i);
      packed.push_back(nB_grid2[i]);
    }
    
    split_string_delim(Ye_grid_spec,split_res,',');
    n_Ye2=stoszt(split_res[0]);
    
    calc.compile(split_res[1].c_str());
    for(size_t i=0;i<n_Ye2;i++) {
      vars["i"]=((double)i);
      packed.push_back(Ye_grid2[i]);
    }
    
    split_string_delim(T_grid_spec,split_res,',');
    n_T2=stoszt(split_res[0]);
    
    calc.compile(split_res[1].c_str());
    for(size_t i=0;i<n_T2;i++) {
      vars["i"]=((double)i);
      packed.push_back(T_grid2[i]);
    }
    
    tg_Eint.resize(3,st);
    tg_Pint.resize(3,st);
    tg_Sint.resize(3,st);
    tg_mun.resize(3,st);
    tg_mup.resize(3,st);
    
    tg_Eint.set_grid_packed(packed);
    tg_Pint.set_grid_packed(packed);
    tg_Sint.set_grid_packed(packed);
    tg_mun.set_grid_packed(packed);
    tg_mup.set_grid_packed(packed);
    
    tg_Eint.set_all(0.0);
    tg_Pint.set_all(0.0);
    tg_Sint.set_all(0.0);
    tg_mun.set_all(0.0);
    tg_mup.set_all(0.0);
  }

  o2scl::tensor_grid<> tg_dmundnB;
  o2scl::tensor_grid<> tg_dmundYe;
  o2scl::tensor_grid<> tg_dmundT;
  o2scl::tensor_grid<> tg_dmupdYe;
  o2scl::tensor_grid<> tg_dmupdT;
  o2scl::tensor_grid<> tg_dsdT;
  
  size_t st[3]={n_nB2,n_Ye2,n_T2};
  vector<vector<double> > grid={nB_grid2,Ye_grid2,T_grid2};
  
  tg_dmundnB.resize(3,st);
  tg_dmundYe.resize(3,st);
  tg_dmundT.resize(3,st);
  tg_dmupdYe.resize(3,st);
  tg_dmupdT.resize(3,st);
  tg_dsdT.resize(3,st);
  
  tg_dmundnB.set_grid(grid);
  tg_dmundYe.set_grid(grid);
  tg_dmundT.set_grid(grid);
  tg_dmupdYe.set_grid(grid);
  tg_dmupdT.set_grid(grid);
  tg_dsdT.set_grid(grid);
  
  // Temporary vectors which store the free energy density in
  // MeV/fm^3
  ubvector fint_vec_nB(n_nB2), fint_vec_Ye(n_Ye2), fint_vec_T(n_T2);

  // The interpolator
  interp_steffen<vector<double>,ubvector> itp_stf;

  // Temporary vectors for derivatives in the nB, Ye, and T directions
  fint_vec_nB.resize(n_nB2);
  fint_vec_Ye.resize(n_Ye2);
  fint_vec_T.resize(n_T2);

  // First, compute the temperature derivative, which is just the
  // entropy times minus one
  for (size_t i=0;i<n_nB2;i++) {
    for (size_t j=0;j<n_Ye2;j++) {
      for(size_t k=0;k<n_T2;k++) {
        vector<size_t> ix={i,j,k};
	fint_vec_T[k]=tg_Fint.get(ix)*nB_grid2[i];
      }
      itp_stf.set(n_T2,T_grid2,fint_vec_T);
      for(size_t k=0;k<n_T2;k++) {
        vector<size_t> ix={i,j,k};
	// Note that fint_vec_T above is stored in MeV/fm^3, so
	// when we take a derivative wrt to temperature (stored
	// in MeV), we get the correct units of 1/fm^3, and then
	// we divide by the baryon density in units of 1/fm^3
	tg_Sint.set(ix,-itp_stf.deriv(T_grid2[k])/nB_grid2[i]);
      }
    }
  }

  // Second, compute the electron fraction derivative, which we
  // temporarily store in tg_mup
  for (size_t i=0;i<n_nB2;i++) {
    for(size_t k=0;k<n_T2;k++) {
      for (size_t j=0;j<n_Ye2;j++) {
        vector<size_t> ix={i,j,k};
	fint_vec_Ye[j]=tg_Fint.get(ix)*nB_grid2[i];
      }
      itp_stf.set(n_Ye2,Ye_grid2,fint_vec_Ye);
      for (size_t j=0;j<n_Ye2;j++) {
        vector<size_t> ix={i,j,k};
	tg_mup.set(ix,itp_stf.deriv(Ye_grid2[j]));
      }
    }
  }

  // Third, compute the baryon density derivative, which we
  // temporarily store in tg_mun
  for (size_t j=0;j<n_Ye2;j++) {
    for(size_t k=0;k<n_T2;k++) {
      for (size_t i=0;i<n_nB2;i++) {
        vector<size_t> ix={i,j,k};
	fint_vec_nB[i]=tg_Fint.get(ix)*nB_grid2[i];
      }
      itp_stf.set(n_nB2,nB_grid2,fint_vec_nB);
      for (size_t i=0;i<n_nB2;i++) {
        vector<size_t> ix={i,j,k};
	tg_mun.set(ix,itp_stf.deriv(nB_grid2[i]));
      }
    }
  }

  // Now go through every point and compute the remaining
  // quantities
  for (size_t i=0;i<n_nB2;i++) {
    for (size_t j=0;j<n_Ye2;j++) {
      for (size_t k=0;k<n_T2;k++) {
	
        vector<size_t> ix={i,j,k};
        
	double en=tg_Sint.get(ix)*nB_grid2[i];
	double dfdnB=tg_mun.get(ix);
	double dfdYe=tg_mup.get(ix);

	tg_mun.get(ix)=dfdnB-dfdYe*Ye_grid2[j]/nB_grid2[i];
	tg_mup.get(ix)=dfdnB-dfdYe*(Ye_grid2[j]-1.0)/nB_grid2[i];

	// E = F + T S
	tg_Eint.get(ix)=tg_Fint.get(ix)+
	  T_grid2[k]*tg_Sint.get(ix);
	
	// P = - F + mun * nn + mup * np
	tg_Pint.get(ix)=-tg_Fint.get(ix)*nB_grid2[i]+
	  nB_grid2[i]*(1.0-Ye_grid2[j])*tg_mun.get(ix)+
	  nB_grid2[i]*Ye_grid2[j]*tg_mup.get(ix);
      }
    }
  }
  
  std::cout << "Finished computing derivatives." << endl;

  return 0;
    
}

int eos_nuclei::eos_second_deriv(std::vector<std::string> &sv,
				 bool itive_com) {

  std::cout << "Computing second derivatives." << endl;
  
  // -----------------------------------------------------
  // Read table
    
  if (derivs_computed==false) {
    cerr << "Fail." << endl;
    return 3;
  }

  o2scl::tensor_grid<> tg_dmundnB;
  o2scl::tensor_grid<> tg_dmundYe;
  o2scl::tensor_grid<> tg_dmundT;
  o2scl::tensor_grid<> tg_dmupdYe;
  o2scl::tensor_grid<> tg_dmupdT;
  o2scl::tensor_grid<> tg_dsdT;
  
  derivs_computed=true;
  
  size_t st[3]={n_nB2,n_Ye2,n_T2};
  vector<vector<double> > grid={nB_grid2,Ye_grid2,T_grid2};
  
  tg_dmundnB.resize(3,st);
  tg_dmundYe.resize(3,st);
  tg_dmundT.resize(3,st);
  tg_dmupdYe.resize(3,st);
  tg_dmupdT.resize(3,st);
  tg_dsdT.resize(3,st);
  
  tg_dmundnB.set_grid(grid);
  tg_dmundYe.set_grid(grid);
  tg_dmundT.set_grid(grid);
  tg_dmupdYe.set_grid(grid);
  tg_dmupdT.set_grid(grid);
  tg_dsdT.set_grid(grid);
  
  // Temporary vectors
  ubvector mun_vec_nB(n_nB2), mun_vec_Ye(n_Ye2), mun_vec_T(n_T2);
  ubvector mup_vec_Ye(n_Ye2), mup_vec_T(n_T2);
  ubvector s_vec_T(n_T2);

  // The interpolator
  interp_steffen<vector<double>,ubvector> itp_stf;

  // First, compute the temperature derivatives
  for (size_t i=0;i<n_nB2;i++) {
    double nB=nB_grid2[i];
    for (size_t j=0;j<n_Ye2;j++) {
      for(size_t k=0;k<n_T2;k++) {
        vector<size_t> ix={i,j,k};
	mun_vec_T[k]=tg_mun.get(ix);
	mup_vec_T[k]=tg_mup.get(ix);
	s_vec_T[k]=tg_S.get(ix)/nB;
      }
      itp_stf.set(n_T2,T_grid2,mun_vec_T);
      for(size_t k=0;k<n_T2;k++) {
        vector<size_t> ix={i,j,k};
	tg_dmundT.set(ix,itp_stf.deriv(T_grid2[k]));
      }
      itp_stf.set(n_T2,T_grid2,mup_vec_T);
      for(size_t k=0;k<n_T2;k++) {
        vector<size_t> ix={i,j,k};
	tg_dmupdT.set(ix,itp_stf.deriv(T_grid2[k]));
      }
      itp_stf.set(n_T2,T_grid2,s_vec_T);
      for(size_t k=0;k<n_T2;k++) {
        vector<size_t> ix={i,j,k};
	tg_dsdT.set(ix,itp_stf.deriv(T_grid2[k]));
      }
    }
  }

  // Second, compute the Ye derivatives
  for (size_t i=0;i<n_nB2;i++) {
    for(size_t k=0;k<n_T2;k++) {
      for (size_t j=0;j<n_Ye2;j++) {
        vector<size_t> ix={i,j,k};
	mun_vec_Ye[j]=tg_mun.get(ix);
	mup_vec_Ye[j]=tg_mup.get(ix);
      }
      itp_stf.set(n_Ye2,Ye_grid2,mun_vec_Ye);
      for (size_t j=0;j<n_Ye2;j++) {
        vector<size_t> ix={i,j,k};
	tg_dmundYe.set(ix,itp_stf.deriv(Ye_grid2[j]));
      }
      itp_stf.set(n_Ye2,Ye_grid2,mup_vec_Ye);
      for (size_t j=0;j<n_Ye2;j++) {
        vector<size_t> ix={i,j,k};
	tg_dmupdYe.set(ix,itp_stf.deriv(Ye_grid2[j]));
      }
    }
  }

  // Third, compute the baryon density derivative
  for (size_t j=0;j<n_Ye2;j++) {
    for(size_t k=0;k<n_T2;k++) {
      for (size_t i=0;i<n_nB2;i++) {
        vector<size_t> ix={i,j,k};
	mun_vec_nB[i]=tg_mun.get(ix);
      }
      itp_stf.set(n_nB2,nB_grid2,mun_vec_nB);
      for (size_t i=0;i<n_nB2;i++) {
        vector<size_t> ix={i,j,k};
	tg_dmundnB.set(ix,itp_stf.deriv(nB_grid2[i]));
      }
    }
  }

  // Now go through every point and compute the remaining
  // quantities
  
  for (size_t i=0;i<n_nB2;i++) {
    double nB=nB_grid2[i];
    for (size_t j=0;j<n_Ye2;j++) {
      double Ye=Ye_grid2[j];
      for (size_t k=0;k<n_T2;k++) {
	double T_MeV=T_grid2[k];
        vector<size_t> ix={i,j,k};

	double dmupdnB=tg_dmundnB.get(ix)+
	  tg_dmundYe.get(ix)*(1.0-Ye)/nB+
	  tg_dmupdYe.get(ix)*Ye/nB;
	  
	double f_nnnn=tg_dmundnB.get(ix)-
	  Ye*tg_dmundYe.get(ix)/nB;
	double f_nnnp=tg_dmundnB.get(ix)+
	  (1.0-Ye)*tg_dmundYe.get(ix)/nB;
	double f_npnp=dmupdnB+
	  (1.0-Ye)*tg_dmupdYe.get(ix)/nB;
	double f_nnT=tg_dmundT.get(ix);
	double f_npT=tg_dmupdT.get(ix);
	double f_TT=-tg_dsdT.get(ix);
	double nn=nB*(1.0-Ye);
	double np=nB*Ye;
	double en=tg_Sint.get(ix)*nB;

	double mue=0.0;
	double den=en*T_MeV+(tg_mun.get(ix)+neutron.m)*nn+
	  (tg_mup.get(ix)+proton.m+mue)*np;
	
	double cs_sq=(nn*nn*(f_nnnn-f_nnT*f_nnT/f_TT)+
		      2.0*nn*np*(f_nnnp-f_nnT*f_npT/f_TT)+
		      np*np*(f_npnp-f_npT*f_npT/f_TT)-
		      2.0*en*(nn*f_nnT/f_TT+np*f_npT/f_TT)-en*en/f_TT)/den;
	
      }
    }
  }
  
  std::cout << "Finished computing second derivatives." << endl;

  return 0;
    
}

int eos_nuclei::stability(std::vector<std::string> &sv,
                          bool itive_com) {


  if (with_leptons==false) {
    cerr << "The 'stability' command requires an EOS table with "
         << "leptons and photons." << endl;
    return 2;
  }
  
  // The six second derivatives we need to compute
  o2scl::tensor_grid<> dmundYe, dmundnB, dmupdYe, dsdT, dsdnB, dsdYe;
  o2scl::tensor_grid<> egv[4], cs2, cs2_hom;

  interp_vec<vector<double> > itp_sta_a, itp_sta_b, itp_sta_c;

  vector<double> packed;
  for(size_t i=0;i<n_nB2;i++) {
    packed.push_back(nB_grid2[i]);
  }
  for(size_t i=0;i<n_Ye2;i++) {
    packed.push_back(Ye_grid2[i]);
  }
  for(size_t i=0;i<n_T2;i++) {
    packed.push_back(T_grid2[i]);
  }
  size_t st[3]={n_nB2,n_Ye2,n_T2};
  
  dmundYe.resize(3,st);
  dmundnB.resize(3,st);
  dmupdYe.resize(3,st);
  dsdT.resize(3,st);
  dsdnB.resize(3,st);
  dsdYe.resize(3,st);
  for(size_t i=0;i<4;i++) {
    egv[i].resize(3,st);
  }
  cs2.resize(3,st);
  cs2_hom.resize(3,st);

  dmundYe.set_grid_packed(packed);
  dmundnB.set_grid_packed(packed);
  dmupdYe.set_grid_packed(packed);
  dsdT.set_grid_packed(packed);
  dsdnB.set_grid_packed(packed);
  dsdYe.set_grid_packed(packed);
  for(size_t i=0;i<4;i++) {
    egv[i].set_grid_packed(packed);
  }
  cs2.set_grid_packed(packed);
  cs2_hom.set_grid_packed(packed);
  
  // The baryon density derivatives
  for (size_t j=0;j<n_Ye2;j++) {
    for (size_t k=0;k<n_T2;k++) {
      vector<double> mun_of_nB, s_of_nB;
      for (size_t i=0;i<n_nB2;i++) {
        vector<size_t> ix={i,j,k};
        mun_of_nB.push_back(tg_mun.get(ix)/hc_mev_fm+
                            neutron.m);
        s_of_nB.push_back(tg_S.get(ix)*nB_grid2[i]);
      }
      itp_sta_a.set(n_nB2,nB_grid2,mun_of_nB,itp_steffen);
      itp_sta_b.set(n_nB2,nB_grid2,s_of_nB,itp_steffen);
      for (size_t i=0;i<n_nB2;i++) {
        vector<size_t> ix={i,j,k};
        dmundnB.get(ix)=itp_sta_a.deriv(nB_grid2[i]);
        dsdnB.get(ix)=itp_sta_b.deriv(nB_grid2[i]);
      }
    }
  }
  
  // The electron fraction derivatives
  for (size_t i=0;i<n_nB2;i++) {
    double nB=nB_grid2[i];
    for (size_t k=0;k<n_T2;k++) {
      vector<double> mun_of_Ye, mup_of_Ye, s_of_Ye;
      for (size_t j=0;j<n_Ye2;j++) {
        vector<size_t> ix={i,j,k};
        mun_of_Ye.push_back(tg_mun.get(ix)/hc_mev_fm+
                            neutron.m);
        mup_of_Ye.push_back(tg_mup.get(ix)/hc_mev_fm+
                            proton.m+tg_mue.get(ix)/hc_mev_fm);
        s_of_Ye.push_back(tg_S.get(ix)*nB);
      }
      itp_sta_a.set(n_Ye2,Ye_grid2,mun_of_Ye,itp_steffen);
      itp_sta_b.set(n_Ye2,Ye_grid2,mup_of_Ye,itp_steffen);
      itp_sta_c.set(n_Ye2,Ye_grid2,s_of_Ye,itp_steffen);
      for (size_t j=0;j<n_Ye2;j++) {
        vector<size_t> ix={i,j,k};
        dmundYe.get(ix)=itp_sta_a.deriv(Ye_grid2[j]);
        dmupdYe.get(ix)=itp_sta_b.deriv(Ye_grid2[j]);
        dsdYe.get(ix)=itp_sta_b.deriv(Ye_grid2[j]);
      }
    }
  }
  
  // The temperature derivative
  for (size_t i=0;i<n_nB2;i++) {
    double nB=nB_grid2[i];
    for (size_t j=0;j<n_Ye2;j++) {
      vector<double> s_of_T;
      for (size_t k=0;k<n_T2;k++) {
        vector<size_t> ix={i,j,k};
        s_of_T.push_back(tg_S.get(ix)*nB);
      }
      itp_sta_a.set(n_T2,T_grid2,s_of_T,itp_steffen);
      for (size_t k=0;k<n_T2;k++) {
        vector<size_t> ix={i,j,k};
        dsdT.get(ix)=itp_sta_a.deriv(T_grid2[k])*hc_mev_fm;
      }
    }
  }
  
  // Storage for the matrix and SVD
  ubmatrix mat(4,4), V(4,4);
  ubvector sing(4), work(4);
  
  int unstable_count=0;
  int superlum_count=0;
  int dPdnB_negative_count=0;
  int PS_negative_count=0;
  int total_count=0;
  
  vector<size_t> i_nB_fix;
  vector<size_t> i_Ye_fix;
  vector<size_t> i_T_fix;
  vector<size_t> type_fix;
  
  string outfile;
  size_t ilo=0, ihi=n_nB2;
  size_t jlo=0, jhi=n_Ye2;
  size_t klo=0, khi=n_T2;
  ilo=260;
  ihi=261;
  if (sv.size()>=4) {
    double nBx=o2scl::function_to_double(sv[1]);
    double Yex=o2scl::function_to_double(sv[2]);
    double Tx=o2scl::function_to_double(sv[3]);
    ilo=vector_lookup(n_nB2,nB_grid2,nBx);
    ihi=ilo+1;
    jlo=vector_lookup(n_Ye2,Ye_grid2,Yex);
    jhi=jlo+1;
    klo=vector_lookup(n_T2,T_grid2,Tx);
    khi=klo+1;
    cout << "Density " << nBx << " index " << ilo << endl;
    cout << "Electron fraction " << Yex << " index " << jlo << endl;
    cout << "Temperature " << Tx << " index " << klo << endl;
  } else {
    outfile=sv[1];
  }
  
  /// Compute the stability matrix and its eigenvalues at each point
  for (size_t i=ilo;i<ihi;i++) {
    double nB=nB_grid2[i];
    cout << "nB at " << i << " is " << nB << endl;
    for (size_t j=jlo;j<jhi;j++) {
      double Ye=Ye_grid2[j];
      for (size_t k=klo;k<khi;k++) {
        
        total_count++;
        
        double T_MeV=T_grid2[k];
        vector<size_t> ix={i,j,k};
        
        // Check entropy and pressure are positive
        if (tg_P.get(ix)<0.0 || tg_S.get(ix)<0.0) {
          cout << "P or S <0: nB,Ye,T[MeV]:\n  "
               << nB << " " << Ye << " " << T_MeV << " "
               << tg_P.get(ix) << " " << tg_S.get(ix) << endl;
          PS_negative_count++;
          i_nB_fix.push_back(i);
          i_Ye_fix.push_back(j);
          i_T_fix.push_back(k);
          type_fix.push_back(1);
        }
        
        // Check dPdnB
        if (i<n_nB2) {
          vector<size_t> ixp1={i+1,j,k};
          double dP=tg_P.get(ixp1)-tg_P.get(ix);
          if (dP<0.0) {
            cout << "dPdnB<0: nB,Ye,T[MeV]:\n  "
                 << nB << " " << Ye << " " << T_MeV << endl;
            if (true) {
              cout << "ix  ,P: ";
              vector_out(cout,ix,false);
              cout << " " << tg_P.get(ix) << endl;
              cout << "ix+1,P: ";
              vector_out(cout,ixp1,false);
              cout << " " << tg_P.get(ixp1) << endl;
              exit(-1);
            }
            dPdnB_negative_count++;
            i_nB_fix.push_back(i);
            i_Ye_fix.push_back(j);
            i_T_fix.push_back(k);
            type_fix.push_back(2);
            i_nB_fix.push_back(i+1);
            i_Ye_fix.push_back(j);
            i_T_fix.push_back(k);
            type_fix.push_back(2);
          }
        }
        
        // Entropy and densities
        double en=tg_S.get(ix)*nB;
        double nn2=(1.0-Ye)*nB;
        double np2=Ye*nB;
	
        // Temporarily store the six derivatives
        double dmundnBv=dmundnB.get(ix);
        double dmundYev=dmundYe.get(ix);
        double dmupdYev=dmupdYe.get(ix);
        double dsdnBv=dsdnB.get(ix);
        double dsdYev=dsdYe.get(ix);
        double dsdTv=dsdT.get(ix);
        
        // Compute dmupdnB, which is related to the other three
        // by a Maxwell relation
        double dmupdnBv=Ye/nB*dmupdYev+dmundnBv+(1.0-Ye)/nB*dmundYev;
        
        // Transform from (Ye,nB) to (nn,np)
        double dmundnn=dmundnBv-Ye/nB*dmundYev;
        double dmundnp=dmundnBv+(1.0-Ye)/nB*dmundYev;
        double dmupdnn=dmundnp;
        double dmupdnp=dmupdnBv+(1.0-Ye)/nB*dmupdYev;
        
        // Use dmundT = -dsdnn and dmupdT = -dsdnp
        double dmundT=-dsdnBv+Ye/nB*dsdYev;
        double dmupdT=-dsdnBv-(1.0-Ye)/nB*dsdYev;
        
        // For the matrix, we need some additional derivatives
        // at fixed entropy
        double dTdnn_s=dmundT/dsdTv;
        double dTdnp_s=dmupdT/dsdTv;
        double dmundnn_s=dmundT*dTdnn_s+dmundnn;
        double dmundnp_s=dmundT*dTdnp_s+dmundnp;
        double dmupdnn_s=dmupdT*dTdnn_s+dmupdnn;
        double dmupdnp_s=dmupdT*dTdnp_s+dmupdnp;
        double dPdnn_s=en*dTdnn_s+nn2*dmundnn_s+np2*dmupdnn_s;
        double dPdnp_s=en*dTdnp_s+nn2*dmundnp_s+np2*dmupdnp_s;
        
        // And one derivative of the pressure with respect
        // to entropy
        double dPds=(en+nn2*dmundT+np2*dmupdT)/dsdTv;
	
        // Now, construct the matrix
        mat(0,0)=nn2*dPdnn_s+np2*dPdnp_s+en*dPds;
        mat(0,1)=-en*dPds;
        mat(0,2)=-nn2*dPdnn_s;
        mat(0,3)=-np2*dPdnp_s;
        mat(1,0)=mat(0,1);
        mat(1,1)=en*en/dsdTv;
        mat(1,2)=nn2*en*dTdnn_s;
        mat(1,3)=np2*en*dTdnp_s;
        mat(2,0)=mat(0,2);
        mat(2,1)=mat(1,2);
        mat(2,2)=nn2*nn2*dmundnn_s;
        mat(2,3)=nn2*np2*dmundnp_s;
        mat(3,0)=mat(0,3);
        mat(3,1)=mat(1,3);
        mat(3,2)=mat(2,3);
        mat(3,3)=np2*np2*dmupdnp_s;
        
#ifdef NEVER_DEFINED
        Eigen::MatrixXd mat2=Eigen::MatrixXd::Ones(4,4);
        for(size_t jj=0;jj<4;jj++) {
          for(size_t kk=0;kk<4;kk++) {
            mat2(jj,kk)=mat(jj,kk);
          }
        }
        Eigen::VectorXcd eivals=mat2.eigenvalues();
#endif
        
        // Compute eigenvalues using SVD
        o2scl_linalg::SV_decomp(4,4,mat,V,sing,work);
        
        for(size_t ij=0;ij<4;ij++) {
          egv[ij].get(ix)=sing[ij];
        }
        
        if (sing[0]<0.0 || sing[1]<0.0 || sing[2]<0.0 ||
            sing[3]<0.0) {
          cout << "Unstable: " << nB << " " << Ye << " " << T_MeV << " "
               << sing[0] << " " << sing[1] << " " << sing[2]
               << " " << sing[3] << endl;
          unstable_count++;
          i_nB_fix.push_back(i);
          i_Ye_fix.push_back(j);
          i_T_fix.push_back(k);
          type_fix.push_back(3);
          //char ch;
          //cin >> ch;
        }
        
        // Compute squared speed of sound
        double f_nnnn=dmundnn;
        double f_nnnp=dmundnp;
        double f_npnp=dmupdnp;
        double f_nnT=dmundT;
        double f_npT=dmupdT;
        double f_TT=-dsdTv;
        double den=en*T_MeV/hc_mev_fm+
          (tg_mun.get(ix)/hc_mev_fm+neutron.m)*nn2+
          (tg_mup.get(ix)/hc_mev_fm+proton.m)*np2+
          tg_mue.get(ix)*np2/hc_mev_fm;
        if (cs2_verbose>0) {
          cout << endl;
          cout << "fr: " << tg_F.get(ix)/hc_mev_fm*nB << endl;
          cout << "en: " << tg_S.get(ix)*nB << endl;
          cout << "mun[1/fm],mup[1/fm],mue[1/fm]: "
               << tg_mun.get(ix)/hc_mev_fm << " "
               << tg_mup.get(ix)/hc_mev_fm << " "
               << tg_mue.get(ix)/hc_mev_fm
               << endl;
          cout << "den1,den2,den3,den4 (all [1/fm^3]):\n  "
               << en*T_MeV/hc_mev_fm << " " 
               << (tg_mun.get(ix)/hc_mev_fm+neutron.m)*nn2 << " " 
               << (tg_mup.get(ix)/hc_mev_fm+proton.m)*np2 << " "
               << tg_mue.get(ix)/hc_mev_fm << " " << np2 << endl;
          cout << "nn,np,en: " << nn2 << " " << np2 << " " << en << endl;
          cout << "f_nnnn, f_nnnp, f_npnp, f_nnT, f_npT, f_TT, den:\n  "
               << f_nnnn << " " << f_nnnp << " " << f_npnp << " "
               << f_nnT << "\n  " << f_npT << " " << f_TT << " "
               << den << endl;
        }
        double cs_sq=(nn2*nn2*(f_nnnn-f_nnT*f_nnT/f_TT)+
                      2.0*nn2*np2*(f_nnnp-f_nnT*f_npT/f_TT)+
                      np2*np2*(f_npnp-f_npT*f_npT/f_TT)-
                      2.0*en*(nn2*f_nnT/f_TT+np2*f_npT/f_TT)-en*en/f_TT)/den;
        
        cs2.get(ix)=cs_sq;
        
        // This code requires a model to compute the homogeneous cs2
        if (true) {
          neutron.n=nB*(1.0-Ye);
          proton.n=nB*Ye;
          thermo th;
          cs2_hom.get(ix)=cs2_func(neutron,proton,T_MeV/hc_mev_fm,th);
          
          if (cs2_verbose>0 || sv.size()>=4) {
            cout << "cs2 (het,hom): " << cs_sq << " "
                 << cs2_hom.get(ix) << endl;
          }
        }
        
        if (cs_sq>1.0 || !std::isfinite(cs_sq)) {
          //cout << tg_mun.get(ix) << " " << neutron.m << " "
          //<< electron.mu << " " << electron.n << endl;
          cout << "Superluminal: nB,Ye,T[MeV],cs2,cs2_hom:\n  "
               << nB << " " << Ye << " " << T_MeV << " "
               << cs_sq << " " << cs2_hom.get(ix) << endl;
          superlum_count++;
          i_nB_fix.push_back(i);
          i_Ye_fix.push_back(j);
          i_T_fix.push_back(k);
          type_fix.push_back(4);
          //exit(-1);
        }
        
      }
    }
  }
  
  cout << "Unstable count (type 3): " << unstable_count << endl;
  cout << "Superluminal count (type 4): " << superlum_count << endl;
  cout << "dPdnB<0 count (type 2): " << dPdnB_negative_count << endl;
  cout << "P|S<0 count (type 1): " << PS_negative_count << endl;
  cout << "Total count: " << total_count << endl;
  
  recompute=true;
  
  if (true) {

    derivs_computed=false;
    with_leptons=false;
    
    for(size_t i=0;i<i_nB_fix.size();i++) {
      
      if (type_fix[i]==1 || type_fix[i]==2) {
        
        double nB=nB_grid2[i_nB_fix[i]];
        double Ye=Ye_grid2[i_Ye_fix[i]];
        double T_MeV=T_grid2[i_T_fix[i]];
        
        cout << "nB, Ye, T [MeV]: " << nB << " " << Ye << " "
             << T_MeV << endl;
        
        vector<string> sv2={"",o2scl::dtos(nB),
          o2scl::dtos(Ye),o2scl::dtos(T_MeV)};
        
        vector<size_t> ix={i_nB_fix[i],i_Ye_fix[i],i_T_fix[i]};
        tg_flag.get(ix)=5;
        
        mh_tol_rel=1.0e-8;
        point_nuclei(sv2,false);
        
        if (tg_flag.get(ix)<10) {
          
          mh_tol_rel=1.0e-7;
          point_nuclei(sv2,false);
          
          if (tg_flag.get(ix)<10) {
            
            mh_tol_rel=1.0e-6;
            point_nuclei(sv2,false);
            
            if (tg_flag.get(ix)<10) {
              
              mh_tol_rel=1.0e-4;
              point_nuclei(sv2,false);
            }
            
          }  
        }
        
        cout << endl;
      }
    }
    
  }
  
  if (false) {
    
    derivs_computed=false;
    with_leptons=false;
    
    cout << "Here: " << i_nB_fix.size() << endl;
    
    vector<bool> in_a_group(i_nB_fix.size());
    vector_set_all(in_a_group,false);
    
    vector<vector<size_t>> i_nB_groups, i_Ye_groups, i_T_groups;
    bool group_done=false;
    
    int countx=0;
    
    while (group_done==false) {
      
      bool found_one=false;
      
      for(size_t i=0;i<i_nB_fix.size();i++) {
        
        if (in_a_group[i]==false) {
          found_one=true;
          
          for(size_t j=0;j<i_nB_groups.size() &&
                in_a_group[i]==false;j++) {
            for(size_t k=0;k<i_nB_groups[j].size() &&
                  in_a_group[i]==false;k++) {
              
              int i1=i_nB_fix[i];
              int i2=i_Ye_fix[i];
              int i3=i_T_fix[i];
              int i4=i_nB_groups[j][k];
              int i5=i_Ye_groups[j][k];
              int i6=i_T_groups[j][k];
              if (abs(i1-i4)+abs(i2-i5)+abs(i3-i6)<=2) {
                /*
                  cout << "Add to group: " << i << " " << i_nB_fix[i] << " "
                  << i_Ye_fix[i] << " " << i_T_fix[i] << " with "
                  << i4 << " " << i5 << " " << i6 << endl;
                */
                i_nB_groups[j].push_back(i1);
                i_Ye_groups[j].push_back(i2);
                i_T_groups[j].push_back(i3);
                in_a_group[i]=true;
              }
            }
          }
          
          if (in_a_group[i]==false) {
            //cout << "Start group: " << i << " " << i_nB_fix[i] << " "
            //<< i_Ye_fix[i] << " " << i_T_fix[i] << endl;
            i_nB_groups.push_back({i_nB_fix[i]});
            i_Ye_groups.push_back({i_Ye_fix[i]});
            i_T_groups.push_back({i_T_fix[i]});
            in_a_group[i]=true;
            //char ch;
            //std::cin >> ch;
          }
          
        }
        
      }
      
      if (found_one==false) {
        group_done=true;
      }
      
      int tot=0;
      for(size_t i=0;i<i_nB_groups.size();i++) {
        cout << "Group: ";
        cout.width(3);
        cout << i << " ";
        cout.width(3);
        cout << i_nB_groups[i].size() << " ";
        cout.width(3);
        cout << i_Ye_groups[i].size() << " ";
        cout.width(3);
        cout << i_T_groups[i].size() << " ";
        tot+=i_nB_groups[i].size();
        if (i_nB_groups[i].size()>2) {
          cout << vector_mean<vector<size_t>>(i_nB_groups[i]) << " "
               << vector_stddev<vector<size_t>>(i_nB_groups[i]) << " ";
          cout << vector_mean<vector<size_t>>(i_Ye_groups[i]) << " "
               << vector_stddev<vector<size_t>>(i_Ye_groups[i]) << " ";
          cout << vector_mean<vector<size_t>>(i_T_groups[i]) << " "
               << vector_stddev<vector<size_t>>(i_T_groups[i]) << endl;
        } else {
          cout << 0.0 << " " << 0.0 << " ";
          cout << 0.0 << " " << 0.0 << " ";
          cout << 0.0 << " " << 0.0 << endl;
        }
      }
      cout << "tot: " << tot << endl;
      
    }
  }
  
  if (outfile.length()>0) {
    
    hdf_file hf;
    hf.open_or_create(outfile);
    hdf_output(hf,dmundnB,"dmundnB");
    hdf_output(hf,dmundYe,"dmundYe");
    hdf_output(hf,dmupdYe,"dmupdYe");
    hdf_output(hf,dsdnB,"dsdnB");
    hdf_output(hf,dsdYe,"dsdYe");
    hdf_output(hf,dsdT,"dsdT");
    hdf_output(hf,egv[0],"egv0");
    hdf_output(hf,egv[1],"egv1");
    hdf_output(hf,egv[2],"egv2");
    hdf_output(hf,egv[3],"egv3");
    hdf_output(hf,cs2,"cs2");
    hdf_output(hf,cs2_hom,"cs2_hom");
    hf.close();
  }
  
  return 0;
}

double eos_nuclei::solve_nuclei_ld
(double x2, size_t nv, const ubvector &x, double nb, double ye, double T,
 int ix, double &mun_gas, double &mup_gas, thermo &th_gas) {
				       
  ubvector xp(nv), yp(nv);
  xp[0]=x[0];
  xp[1]=x[1];
  xp[ix]=x2;
  map<std::string,double> vdet;
  int ret=solve_nuclei(nv,xp,yp,nb,ye,T,0,mun_gas,mup_gas,th_gas,vdet);
  if (ret!=0) return pow(10.0,80.0+((double(ret))));
  return yp[ix];
}

double eos_nuclei::solve_nuclei_min
(size_t nv, const ubvector &x, double nb, double ye, double T,
 double &mun_gas, double &mup_gas, thermo &th_gas) {

  double retval;
  ubvector yp(2);
  map<string,double> vdet;
  int ret=solve_nuclei(nv,x,yp,nb,ye,T,0,mun_gas,mup_gas,th_gas,vdet);
  retval=yp[0]*yp[0]+yp[1]*yp[1];
  if (ret!=0) retval=pow(10.0,80.0+((double(ret))));
  if (!std::isfinite(yp[0]) || !std::isfinite(yp[1]) ||
      !std::isfinite(retval)) {
    return 1.0e99;
  }
  return retval;
}

void eos_nuclei::store_hrg(double mun, double mup,
                           double nn, double np, double T,
                           table_units<> &tab) {
  
  tab.clear();
  
  tab.line_of_names("id mass spin_deg mu n ed pr en baryon charge");
  tab.line_of_units(". MeV . MeV 1/fm^3 MeV/fm^4 MeV/fm^4 1/fm^3 . .");
  
  int iferm=0;
  int ibos=0;
  for(size_t j=0;j<part_db.size();j++) {
    if (part_db[j].id==22) {
      photon.massless_calc(T);
      vector<double> line={((double)22),0.0,2,
        photon.mu,photon.n,photon.ed,photon.pr,photon.en,
        ((double)part_db[j].baryon),((double)part_db[j].charge)};
      tab.line_of_data(line.size(),line);
      ibos++;
    } else if (part_db[j].id==2212 && part_db[j].charge==1) {
      vector<double> line={((double)2212),neutron.m*hc_mev_fm,2,
        mun*hc_mev_fm,nn,0.0,0.0,0.0,
        ((double)part_db[j].baryon),((double)part_db[j].charge)};
      tab.line_of_data(line.size(),line);
      iferm++;
    } else if (part_db[j].id==2212 && part_db[j].charge==0) {
      vector<double> line={((double)2212),neutron.m*hc_mev_fm,2,
        mun*hc_mev_fm,nn,0.0,0.0,0.0,
        ((double)part_db[j].baryon),((double)part_db[j].charge)};
      tab.line_of_data(line.size(),line);
      iferm++;
    } else if (part_db[j].spin_deg%2==0) {
      vector<double> line={((double)part_db[j].id),
        res_f[iferm].m*hc_mev_fm,
        res_f[iferm].g,res_f[iferm].mu*hc_mev_fm,res_f[iferm].n,
        res_f[iferm].ed*hc_mev_fm,res_f[iferm].pr*hc_mev_fm,
        res_f[iferm].en,
        ((double)part_db[j].baryon),((double)part_db[j].charge)};
      tab.line_of_data(line.size(),line);
      iferm++;
    } else {
      vector<double> line={((double)part_db[j].id),
        res_b[ibos].m*hc_mev_fm,
        res_b[ibos].g,res_b[ibos].mu*hc_mev_fm,res_b[ibos].n,
        res_b[ibos].ed*hc_mev_fm,res_b[ibos].pr*hc_mev_fm,
        res_b[ibos].en,
        ((double)part_db[j].baryon),((double)part_db[j].charge)};
      tab.line_of_data(line.size(),line);
      ibos++;
    }
  }
  
  return;
}

int eos_nuclei::solve_nuclei_mu
(size_t nv, const ubvector &x, ubvector &y, double mun, double mup, double T,
 double &mun_gas, double &mup_gas, thermo &th_gas) {
  
  ubvector xt(2), yt(2);
  int ret;
  double fnb1, fnb2, fye1, fye2;
  
  double nb=x[0];
  double ye=x[1];
  double log_xn_0=x[2];
  double log_xp_0=x[3];
  double log_xn_1=x[4];
  double log_xp_1=x[5];
  double log_xn_2=x[6];
  double log_xp_2=x[7];
  double log_xn_3=x[8];
  double log_xp_3=x[9];

  if (nb<1.0e-12 || nb>2.0) {
    return 10;
  }
  if (ye<0.0 || ye>0.7) {
    return 11;
  }

  map<std::string,double> vdet;
  thermo th;
  
  xt[0]=log_xn_0;
  xt[1]=log_xp_0;
  ret=solve_nuclei(2,xt,yt,nb*(1.0-1.0e-4),ye,T,0,mun_gas,mup_gas,
                   th_gas,vdet);
  if (ret!=0) {
    cerr << "First solve_nuclei function failed." << endl;
    return ret;
  }
  fnb1=compute_fr_nuclei(nb,ye,T,xt[0],xt[1],th,th_gas);
  y[2]=yt[0];
  y[3]=yt[1];

  xt[0]=log_xn_1;
  xt[1]=log_xp_1;
  ret=solve_nuclei(2,xt,yt,nb*(1.0+1.0e-4),ye,T,0,mun_gas,mup_gas,
                   th_gas,vdet);
  if (ret!=0) {
    cerr << "Second solve_nuclei function failed." << endl;
    return ret;
  }
  fnb2=compute_fr_nuclei(nb,ye,T,xt[0],xt[1],th,th_gas);
  y[4]=yt[0];
  y[5]=yt[1];

  xt[0]=log_xn_2;
  xt[1]=log_xp_2;
  ret=solve_nuclei(2,xt,yt,nb,ye*(1.0-1.0e-4),T,0,mun_gas,mup_gas,
                   th_gas,vdet);
  if (ret!=0) {
    cerr << "Third solve_nuclei function failed." << endl;
    return ret;
  }
  fye1=compute_fr_nuclei(nb,ye,T,xt[0],xt[1],th,th_gas);
  y[6]=yt[0];
  y[7]=yt[1];

  xt[0]=log_xn_3;
  xt[1]=log_xp_3;
  ret=solve_nuclei(2,xt,yt,nb,ye*(1.0+1.0e-4),T,0,mun_gas,mup_gas,
                   th_gas,vdet);
  if (ret!=0) {
    cerr << "Fourth solve_nuclei function failed." << endl;
    return ret;
  }
  fye2=compute_fr_nuclei(nb,ye,T,xt[0],xt[1],th,th_gas);
  y[8]=yt[0];
  y[9]=yt[1];
  
  double dfdnB=(fnb2-fnb1)/(nb*(1.0+2.0e-4));
  double dfdYe=(fye2-fye1)/(ye*(1.0+2.0e-4));
  
  double mun2=dfdnB-dfdYe*ye/nb;
  double mup2=dfdnB-dfdYe*(ye-1.0)/nb;

  y[0]=mun2-mun;
  y[1]=mup2-mup;

  if (true) {
    // For debugging
    std::cout << "mu: " << mun2 << " " << mun << " " << mup2 << " " << mup
              << std::endl;
    std::cout << "x: " << std::flush;
    vector_out(cout,x,true);
    std::cout << "y: " << std::flush;
    vector_out(cout,y,true);
  }

  return 0;
}

double eos_nuclei::compute_fr_nuclei
(double nB, double Ye,
 double T, double log_xn, double log_xp, thermo &thx, thermo &th_gas) {

  // quantities left
  size_t n_nuclei=nuclei.size();
  
  static const double n0=0.16;
  
  // -------------------------------------------------------------
  // Compute free energy density and entropy density
  
  double xn=pow(10.0,log_xn);
  double xp=pow(10.0,log_xp);

  double kappa=1.0-nB/n0;
  double xi=kappa/(1.0-nB*xn/n0-nB*xp/n0);

  // The total of the number density over all nuclei
  double sum_nuc=0.0;

  // Begin with zero and then add up contributions. Free energy
  // density in fm^{-4} and entropy density in fm^{-3}
  double f=0.0;

  thx.en=0.0;

  for (size_t i=0;i<n_nuclei;i++) {
    
    double en_nuc, fr_nuc;
    
    if (nuclei[i].n>1.0e-300) {
      double lambda=sqrt(2.0*pi/nuclei[i].m/T);
      fr_nuc=-T*(log(vomega[i]/nuclei[i].n/pow(lambda,3.0))+1.0)*
       nuclei[i].n;
      // Note that everything here, including vomega_prime, is in
      // units of powers of femtometers
      en_nuc=nuclei[i].n*(log(vomega[i]/nuclei[i].n/pow(lambda,3.0))+
                         5.0/2.0+vomega_prime[i]*T/vomega[i]);
      if (!std::isfinite(fr_nuc)) {
       if (nuclei[i].n<1.0e-200) {
         nuclei[i].n=0.0;
         fr_nuc=0.0;
         en_nuc=0.0;
       } else {
         cout << "Nuclear free energy not finite in eos_fixed_dist()."
               << endl;
         cout << nuclei[i].n << " " << nuclei[i].be << " " << lambda
              << " " << vomega[i] << endl;
         exit(-1);
       }
      }
    } else {
      nuclei[i].n=0.0;
      fr_nuc=0.0;
      en_nuc=0.0;
    }
    fr_nuc+=nuclei[i].n*nuclei[i].be+1.433e-05*
      pow(nuclei[i].Z,2.39)/hc_mev_fm*nuclei[i].n;
    sum_nuc+=nuclei[i].n;
    
    f+=fr_nuc+nuclei[i].n*Ec[i];
    thx.en+=en_nuc;
  }

  // Final calculations of free energy density, entropy
  // density, and energy density
  
  f+=xi*(th_gas.ed-T*th_gas.en)-T*sum_nuc*log(kappa);
  thx.en+=xi*th_gas.en+sum_nuc*log(kappa);
  thx.ed=f+T*thx.en;

  return f;
}

int eos_nuclei::solve_nuclei(size_t nv, const ubvector &x, ubvector &y,
			     double nB, double Ye, double T,
			     int loc_verbose, double &mun_gas,
			     double &mup_gas, thermo &th_gas,
			     std::map<string,double> &vdet) {
  
  double log_xn=x[0];
  double log_xp=x[1];
  loc_verbose=0;

  // Set two flags help us handle low-density regimes
  bool nn_small=false;
  bool np_small=false;

  double xn;
  if (log_xn<-300.0) {
    if (alg_mode==1 || alg_mode==3) return 10;
    nn_small=true;
    xn=1.0e-100;
  } else {
    xn=pow(10.0,log_xn);
  }
  double xp;
  if (log_xp<-300.0) {
    if (alg_mode==1 || alg_mode==3) return 11;
    np_small=true;
    xp=1.0e-100;
  } else {
    xp=pow(10.0,log_xp);
  }
  
  // Compute kappa and xi
  double n0=0.16;
  double T_MeV=T*hc_mev_fm;
  double kappa=1.0-nB/n0;
  double xi=kappa/(1.0-nB*xn/n0-nB*xp/n0);

  // Disallow unphysical solutions
  if (xn<=0.0 || xn>1.0) {
    if (loc_verbose>0) {
      cout << "Failure 1." << endl;
    }
    return 1;
  }
  if (xp<=0.0 || xp>1.0) {
    if (loc_verbose>0) {
      cout << "Failure 2." << endl;
    }
    return 2;
  }
  if (kappa<0.0 || kappa>1.0) {
    if (loc_verbose>0) {
      cout << "Failure 3." << endl;
    }
    return 3;
  }
  if (xi<0.0 || xi-1.0>1.0e-15) {
    if (loc_verbose>0) {
      cout << "Failure 4." << endl;
    }
    return 4;
  }
  
  size_t n_nuclei=nuclei.size();

  // Compute Coulomb energy for each nucleus
  ubvector xx(n_nuclei);
  for(size_t i=0;i<n_nuclei;i++) {
    double A=nuclei[i].A;
    // Nuclear radius in fm
    double rA=cbrt(3.0*A/4.0/pi/n0);
    xx[i]=cbrt(nB*Ye/n0*A/nuclei[i].Z);
    // Coulomb energy in 1/fm
    Ec[i]=-0.6*nuclei[i].Z*nuclei[i].Z*fine_structure_f<double>()/rA*
      (1.5*xx[i]-0.5*pow(xx[i],3.0));
  }

  // Compute mun_old and mup_old from DSH
  double nn_prime=nB*xn, np_prime=nB*xp;

  neutron.n=nn_prime;
  proton.n=np_prime;

  // Ensure that this function is deterministic
  neutron.mu=neutron.m;
  proton.mu=proton.m;
  
  if (np_small && nn_small) {

    // Use shift to compute correct proton chemical potential
    proton.mu=T*log(1.0/proton.g*pow(2.0*pi/proton.ms/T,1.5))+
      log_xp*T*log(10.0);
    neutron.mu=T*log(1.0/neutron.g*pow(2.0*pi/neutron.ms/T,1.5))+
      log_xp*T*log(10.0);
    
  } else if (np_small) {
    
    // Compute chemical potential shift at a fiducial proton density
    proton.n=1.0e-100;

    free_energy_density(neutron,proton,T,th2);
    double mup_shift=proton.mu-
      T*log(1.0/proton.g*pow(2.0*pi/proton.ms/T,1.5))+100.0*T*log(10.0);

    // Use shift to compute correct proton chemical potential
    proton.mu=mup_shift+T*log(1.0/proton.g*pow(2.0*pi/proton.ms/T,1.5))+
      log_xp*T*log(10.0);
    
  } else if (nn_small) {

    // Compute chemical potential shift at a fiducial neutron density
    neutron.n=1.0e-100;

    free_energy_density(neutron,proton,T,th2);
    double mun_shift=neutron.mu-
      T*log(1.0/neutron.g*pow(2.0*pi/neutron.ms/T,1.5))+100.0*T*log(10.0);

    // Use shift to compute correct neutron chemical potential
    neutron.mu=mun_shift+T*log(1.0/neutron.g*pow(2.0*pi/neutron.ms/T,1.5))+
      log_xp*T*log(10.0);
    
  } else {
    free_energy_density_detail(neutron,proton,T,th_gas,vdet);
  }

  mun_gas=neutron.mu;
  mup_gas=proton.mu;
  double nn=nn_prime*xi;
  double np=np_prime*xi;

  // Set densities and accumulate results for nn_tilde and np_tilde
  double nn_tilde=nn;
  double np_tilde=np;

  for(size_t i=0;i<n_nuclei;i++) {

    // AWS 4/2/19: If x>1, then the nuclear radius is larger than
    // the WS cell size. Updated on 8/30
    if ((rnuc_less_rws && xx[i]>1.0) || nuclei[i].N>max_ratio*nuclei[i].Z ||
	nuclei[i].Z>max_ratio*nuclei[i].N) {

      nuclei[i].n=0.0;
      
    } else {
      
      // The nucleon de Broglie wavelength in 1/fm
      double lambda=sqrt(2.0*pi/nuclei[i].m/T);
      // Volume in fm^3
      double vv=(nuclei[i].A)/n0;

      // Section for testing removal of electron binding energy
      if (false) {
	cout.precision(10);
	frdm.get_nucleus(nuclei[i].Z,nuclei[i].N,nuclei[i]);
	cout << "Binding energy with electron binding: " 
	     << nuclei[i].be << endl;
	frdm.ael=0.0;
	frdm.get_nucleus(nuclei[i].Z,nuclei[i].N,nuclei[i]);
	cout << "Binding energy without electron binding: " 
	     << nuclei[i].be << endl;
	frdm.ael=1.433e-5;
	frdm.get_nucleus(nuclei[i].Z,nuclei[i].N,nuclei[i]);
	cout << "Binding energy expression: " << 
	  nuclei[i].be+1.433e-5*pow(nuclei[i].Z,2.39)/hc_mev_fm << endl;
	exit(-1);
      }
      
      double arg=(((double)nuclei[i].N)*mun_gas+
		  ((double)nuclei[i].Z)*mup_gas
		  -nuclei[i].be-1.433e-05*pow(nuclei[i].Z,2.39)/hc_mev_fm
		  -Ec[i]-th_gas.pr*vv)/T;
      if (arg>700.0) {
	return 7;
      } else if (arg<-700.0) {
	nuclei[i].n=0.0;
      } else {
	nuclei[i].n=kappa*vomega[i]/pow(lambda,3.0)*exp(arg);
      }

      /*
	AWS, 8/30/19: This is unnecessary because configurations with
	n*A>nB are impossible if the equations are solved. Removing
	this will probably help the solver.
      */
      if ((alg_mode==1 || alg_mode==3) &&
	  nuclei[i].n*nuclei[i].A>1.0e4*nB) {
        if (loc_verbose>0) {
	  cout << "Failure 5." << endl;
        }
 	return 5;
      }
      
      double A=nuclei[i].A;
      // Nuclear radius in fm
      double rA=cbrt(3.0*A/4.0/pi/n0);

      // Add nuclear contribution to total number of neutrons and protons
      nn_tilde+=nuclei[i].N*nuclei[i].n;
      np_tilde+=nuclei[i].Z*nuclei[i].n;
      
      if (loc_verbose>0) {
	cout << "i,N,Z,mun_gas,mup_gas,be: "
	     << i << " " << nuclei[i].N << " " << nuclei[i].Z << " "
	     << mun_gas << " " << mup_gas << " " << nuclei[i].be << endl;
	cout << "Ec,P*V,mun cont,mup cont,n: " << Ec[i] << " "
	     << 1.433e-05*pow(nuclei[i].Z,2.39) << " "
	     << th_gas.pr*vv << " "
	     << (((double)nuclei[i].N)*mun_gas+
		 ((double)nuclei[i].Z)*mup_gas
		 -nuclei[i].be-Ec[i]-1.433e-05*pow(nuclei[i].Z,2.39)
		 -th_gas.pr*vv) << " " << nuclei[i].n << endl;
	cout << "mun , mup " << ((double)nuclei[i].N)*mun_gas << " "
	     << ((double)nuclei[i].Z)*mup_gas << endl;
	cout << "binding " << nuclei[i].be << endl;
	cout << "th_gas.pr, vv" << th_gas.pr << " " << vv << endl;
	cout << "xn, xp " << xn << " " << xp << endl;
	cout << "f_gas " << neutron.mu*nn_prime+proton.mu*np_prime
	  -th_gas.pr << endl;
        cout << mun_gas*n0-th_gas.pr << " " << mup_gas*n0-th_gas.pr << endl;
      }
    }

    // End of loop over nuclei
  }
  
  //inc_hrg loop
  if (false) {

    double nn_tilde_old=nn_tilde;
    double np_tilde_old=np_tilde;
    
    double nB2=0.0, Ye2=0.0;
    
    int iferm=0;
    int ibos=0;
    for(size_t j=0;j<part_db.size();j++) {
      if (part_db[j].id==2212 && part_db[j].charge==1) {
        if (iferm>=((int)res_f.size())) {
          O2SCL_ERR("Indexing problem with fermions.",o2scl::exc_efailed);
        }
        res_f[iferm]=proton;
        iferm++;
      } else if (part_db[j].id==2212 && part_db[j].charge==0) {
        if (iferm>=((int)res_f.size())) {
          O2SCL_ERR("Indexing problem with fermions.",o2scl::exc_efailed);
        }
        res_f[iferm]=neutron;
        iferm++;
      } else if (part_db[j].id==22 && part_db[j].charge==0) {
        if (ibos>=((int)res_b.size())) {
          O2SCL_ERR("Indexing problem with bosons.",o2scl::exc_efailed);
        }
        photon.massless_calc(T);
        res_b[ibos]=photon;
        ibos++;
      } else if (part_db[j].spin_deg%2==0) {
        if (iferm>=((int)res_f.size())) {
          O2SCL_ERR("Indexing problem with fermions.",o2scl::exc_efailed);
        }
        res_f[iferm].mu=part_db[j].baryon*(neutron.mu+neutron.m)+
          part_db[j].charge*(proton.mu+proton.m-neutron.mu-neutron.m);
        relf.calc_mu(res_f[iferm],T);
        if (j>700) {
          //cout << j << " " << part_db[j].baryon << " "
            //   << part_db[j].charge << " "
            //   << res_f[iferm].mu << " " << res_f[iferm].n << endl;
        }
        nB2+=part_db[j].baryon*res_f[iferm].n;
        Ye2+=part_db[j].charge*res_f[iferm].n;
        iferm++;
      } else if (part_db[j].id==-211 && part_db[j].charge==-1){
        cout <<"solve_nuclei: " << j << part_db[j].id << " " << part_db[j].name << endl;
        if (ibos>=((int)res_b.size())) {
          O2SCL_ERR("Indexing problem with bosons.",o2scl::exc_efailed);
        }
        cout << "pion mass: " << res_b[ibos].m*hc_mev_fm << endl;
        res_b[ibos].mu=part_db[j].baryon*(neutron.mu+neutron.m)+
          part_db[j].charge*(proton.mu+proton.m-neutron.mu-neutron.m);
        //cout <<"Pion chem pot: " << res_b[ibos].mu*hc_mev_fm << endl;
        fr.calc_mu(res_b[ibos], neutron, proton, T, nB);
        //cout <<"n_e: " << electron.n << ", n_mu: " << muon.n << ", n_pi: " << res_b[ibos].n << endl;
        //cout << "n_p: " << proton.n << ",n_tot: " << electron.n+muon.n+res_b[ibos].n << endl;
        nB2+=part_db[j].baryon*res_b[ibos].n;
        Ye2+=part_db[j].charge*res_b[ibos].n;
        ibos++;
      } else {
      res_b[ibos].mu=part_db[j].baryon*(neutron.mu+neutron.m)+
        part_db[j].charge*(proton.mu+proton.m-neutron.mu-neutron.m);
      effb.calc_mu(res_b[ibos],T);
      nB2+=part_db[j].baryon*res_b[ibos].n;
      Ye2+=part_db[j].charge*res_b[ibos].n;
      ibos++;
      }
      if (!std::isfinite(nB2) ||
          !std::isfinite(Ye2)) {
        cout << j << " " << part_db[j].id << " "
             << nB2 << " " << Ye2 << endl;
        O2SCL_ERR("HRG problem.",o2scl::exc_einval);
      }
    }
    Ye2/=nB2;

    // I'm not 100% sure this is right
    nn_tilde+=nB2*(1.0-Ye2);
    np_tilde+=nB2*Ye2;

    //nB2=neutron.n+proton.n+delta_pp.n;
    //Ye2=(proton.n-res_b[ibos].n)/nB2;
    cout << "HRG: mu_n: " << neutron.mu << ", mu_p: " << proton.mu << ", nn~_old: "
        << nn_tilde_old << endl;
    cout << "HRG: np~_old: " << np_tilde_old << ", nn~: "
        << nn_tilde << ", np~: " << np_tilde << endl;
    if (!std::isfinite(nn_tilde) ||
        !std::isfinite(np_tilde)) {
      O2SCL_ERR("HRG problem.",o2scl::exc_einval);
    }
    //char ch;
    //cin >> ch;
    
  }
  
  double Ymu=0.0;
  
  if (include_muons) {
    
    eos_sn_base eso;
    eso.include_muons=true;
    thermo lep;
    cout << "vdet[mume]: " << vdet["mue"] << endl;
    eso.compute_eg_point(nB,Ye,T*hc_mev_fm,lep,vdet["mue"]);
    Ymu=eso.muon.n/nB;
    vdet["Ymu"]=Ymu;

    /*
      cout << "e,mu: " << eso.electron.mu << " " << eso.electron.m << " "
      << eso.electron.n << endl;
      cout << "e,mu: " << eso.muon.mu << " " << eso.muon.m << " "
      << eso.muon.n << endl;
    */
    
    // Ensure that we have the correct values of nB and Ye
    y[0]=nn_tilde/nB/(1.0-Ye-Ymu)-1.0;
    y[1]=np_tilde/nB/(Ye+Ymu)-1.0;

  } else {
  
    // Ensure that we have the correct values of nB and Ye
    y[0]=nn_tilde/nB/(1.0-Ye)-1.0;
    y[1]=np_tilde/nB/Ye-1.0;

  }

  /*
    cout << "x,y: " << x[0] << " " << x[1] << " "
    << y[0] << " " << y[1] << " " << nB << " " << Ye << " " << Ymu << endl;
    cout << nn_tilde << " " << nB*(1.0-Ye-Ymu) << endl;
    cout << np_tilde << " " << nB*(Ye+Ymu) << endl;
  */
    
  if (loc_verbose>0) {
    cout << "nn,np: " << nn << " " << np << endl;
    cout << "x: " << x[0] << " " << x[1] << endl;
    cout << "y: " << y[0] << " " << y[1] << endl;
  }

  //cout << "Z: " << x[0] << " " << x[1] << " " << y[0] << " "
  //<< y[1] << " " << fabs(y[0])+fabs(y[1]) << endl;
  
  // Don't allow infinite values
  if (!std::isfinite(y[0]) || !std::isfinite(y[1])) {
    if (loc_verbose>0) {
      cout << "Failure 6." << endl;
    }
    return 6;
  }
    
  if (loc_verbose>2) {
    exit(-1);
  }

  return 0;
}

int eos_nuclei::eos_fixed_ZN(double nB, double Ye, double T,
			     double &log_xn, double &log_xp,
			     size_t nucZ1, size_t nucN1,
			     thermo &thx, double &mun_full,
			     double &mup_full) {

  bool debug=false;
  if (function_verbose%10>1) {
    debug=true;
  }

  // The free energy density in 1/fm^4
  double f;

  double T_MeV=T*hc_mev_fm;
  double n0=0.16;
  
  // Chemical potentials and thermodynamic quantities for homogeneous
  // matter
  double mun_gas, mup_gas;
  thermo th_gas;
  
  if (nucZ1<8 || nucN1<8) {
    return 1;
  }
  
  // Obtain heavy nucleus and neutron and proton separation energies.
  // Ensuring that all three nuclei are in the same mass formula
  // avoids problems with the neutron and separation energies at the
  // boundary between the mass formulae.

  o2scl::nucleus nuc_temp;
  
  if (ame.is_included(nucZ1,nucN1) &&
      ame.is_included(nucZ1-1,nucN1) &&
      ame.is_included(nucZ1,nucN1-1)) {
    
    ame.get_nucleus(nucZ1,nucN1,*nuc_heavy);
    ame.get_nucleus(nucZ1,nucN1-1,nuc_temp);
    Sneut[5]=-(nuclei[5].be-nuc_temp.be)*hc_mev_fm;
    ame.get_nucleus(nucZ1-1,nucN1,nuc_temp);
    Sprot[5]=-(nuclei[5].be-nuc_temp.be)*hc_mev_fm;
    
  } else if (m95.is_included(nucZ1,nucN1) &&
	     m95.is_included(nucZ1-1,nucN1) &&
	     m95.is_included(nucZ1,nucN1-1)) {

    m95.get_nucleus(nucZ1,nucN1,*nuc_heavy);
    m95.get_nucleus(nucZ1,nucN1-1,nuc_temp);
    Sneut[5]=-(nuclei[5].be-nuc_temp.be)*hc_mev_fm;
    m95.get_nucleus(nucZ1-1,nucN1,nuc_temp);
    Sprot[5]=-(nuclei[5].be-nuc_temp.be)*hc_mev_fm;
    
  } else if (extend_frdm) {
    
    frdm.get_nucleus(nucZ1,nucN1,*nuc_heavy);
    frdm.get_nucleus(nucZ1,nucN1-1,nuc_temp);
    Sneut[5]=-(nuclei[5].be-nuc_temp.be)*hc_mev_fm;
    frdm.get_nucleus(nucZ1-1,nucN1,nuc_temp);
    Sprot[5]=-(nuclei[5].be-nuc_temp.be)*hc_mev_fm;

  } else {

    O2SCL_ERR("No nucleus mass information available.",
              o2scl::exc_einval);
    
  }
  
  // Set spin degeneracy of heavy nucleus
  if (hfb.is_included(nucZ1,nucN1)) {
    if (hfb.get_ZN(nucZ1,nucN1).Jexp<99) {
      nuc_heavy->g=2.0*hfb.get_ZN(nucZ1,nucN1).Jexp+1.0;
    } else {
      nuc_heavy->g=2.0*hfb.get_ZN(nucZ1,nucN1).Jth+1.0;
    }
  } else {
    if (nucZ1%2==0 && nucN1%2==0) {
      nuc_heavy->g=1.0;
      //} else if (nucZ1%2==1 && nucN1%2==1) {
      //nuc_heavy->g=3.0;
    } else {
      nuc_heavy->g=2.0;
    }
  }

  // Set up integrand for partition function
  part_func.T_MeV=T_MeV;
  funct f1=std::bind(std::mem_fn<double(double)>
		     (&partition_func::delta_small_iand),&part_func,
		     std::placeholders::_1);
  funct f1_prime=std::bind(std::mem_fn<double(double)>
			   (&partition_func::delta_small_iand_prime),&part_func,
			   std::placeholders::_1);
  funct f2=std::bind(std::mem_fn<double(double)>
		     (&partition_func::delta_large_iand),&part_func,
		     std::placeholders::_1);
  funct f2_prime=std::bind(std::mem_fn<double(double)>
			   (&partition_func::delta_large_iand_prime),
			   &part_func,std::placeholders::_1);
  
  double res, err, ret, ret_prime, res_prime, err_prime;
  
  for (int i=0;i<6;i++) {

    if (i==1) {
      
      vomega[1]=3.0;
      vomega_prime[i]=0.0;
      
    } else {

      // MeV
      double zEd=min(Sneut[i],Sprot[i])/2.0;
      // fm
      double zR=1.25*cbrt(nuclei[i].Z+nuclei[i].N-1.0);
      // MeV since nuclei[i].m is in fm^{-1}
      double zER=0.5/nuclei[i].m/zR/zR*hc_mev_fm;
      // MeV
      double zEc=(nuclei[i].Z-1.0)*fine_structure_f<double>()/zR*hc_mev_fm;
      // MeV
      double zEt=min(Sneut[i]+zER,Sprot[i]+zER+zEc/2.0);
      // fm
      double zrA=cbrt(3.0*(nuclei[i].N+nuclei[i].Z)/4.0/pi/n0);
      // MeV
      double delta_p=11.0/sqrt(nuclei[i].Z+nuclei[i].N)*
	(1.0+pow(-1.0,nuclei[i].Z)/2.0+pow(-1.0,nuclei[i].N)/2.0);
      if (nuclei[i].Z<=30) {
	// 1/MeV
	part_func.a=0.052*pow(nuclei[i].N+nuclei[i].Z,1.2);
	// MeV
	part_func.delta=delta_p-80.0/(nuclei[i].Z+nuclei[i].N);
      } else {
	// 1/MeV
	part_func.a=0.125*(nuclei[i].N+nuclei[i].Z);
	// MeV
	part_func.delta=delta_p-80.0/(nuclei[i].Z+nuclei[i].N)-0.5;
      }
      if (!std::isfinite(part_func.a)) {
	cout << nuclei[i].N << " " << nuclei[i].Z << " " << delta_p << " "
	     << part_func.a << endl;
	cout << "a not finite." << endl;
	exit(-1);
      }
      if (!std::isfinite(part_func.delta)) {
	cout << nuclei[i].N << " " << nuclei[i].Z << " "
	     << part_func.delta << endl;
	cout << "delta not finite." << endl;
	exit(-1);
      }
      
      if (Sneut[i]<0.0 || Sprot[i]<0.0) {
	
	vomega[i]=nuclei[i].g;
	vomega_prime[i]=0.0;
	
      } else {
	
	if (part_func.delta>zEd) {
	  // MeV
	  part_func.delta=zEd;
	  // MeV
	  part_func.Tc=1.0/(-1.25/part_func.delta+sqrt(part_func.a)/
			    sqrt(part_func.delta));
	  // 1/MeV
	  part_func.C=sqrt(pi)/12.0*pow(part_func.a,-0.25)*
	    pow(part_func.delta,-1.25)*
	    exp(1.25+sqrt(part_func.a*part_func.delta));
	  
	  if(2.0*zEd<zEt) {
	    if (!std::isfinite(zEd)) {
	      cout << "zEd not finite (1)." << endl;
	      exit(-1);
	    }
	    if (!std::isfinite(zEt)) {
	      cout << "zEt not finite (1)." << endl;
	      exit(-1);
	    }
	    ret=iqg.integ_err(f1,2.0*zEd,zEt,res,err)+
	      iqg.integ_err(f2,zEd,2.0*zEd,res,err);
	    ret_prime=iqg.integ_err(f1_prime,2.0*zEd,zEt,res_prime,err_prime)+
	      iqg.integ_err(f2_prime,zEd,2.0*zEd,res_prime,err_prime);
	    
	  } else {
	    
	    ret=iqg.integ_err(f2,zEd,zEt,res,err);
	    ret_prime=iqg.integ_err(f2_prime,zEd,zEt,res_prime,err_prime);
	    
	  }
	} else {
	  if (!std::isfinite(zEd)) {
	    cout << "zEd not finite (2)." << endl;
	    exit(-1);
	  }
	  if (!std::isfinite(zEt)) {
	    cout << "zEt not finite (2)." << endl;
	    exit(-1);
	  }
	  ret=iqg.integ_err(f1,zEd,zEt,res,err);
	  ret_prime=iqg.integ_err(f1_prime,zEd,zEt,res_prime,err_prime);
	}
	vomega[i]=nuclei[i].g+res;
	vomega_prime[i]=res_prime/T;
      }
      
    }
  }

  ubvector x1(2), y1(2);

  map<string,double> vdet;
  mm_funct sn_func=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector&,double,double,
		     double,int,double &,double &,thermo &,
		     std::map<string,double> &)>
     (&eos_nuclei::solve_nuclei),this,std::placeholders::_1,
     std::placeholders::_2,std::placeholders::_3,nB,Ye,T,0,
     std::ref(mun_gas),std::ref(mup_gas),std::ref(th_gas),std::ref(vdet));
  
  x1[0]=log_xn;
  x1[1]=log_xp;
  
  mh.tol_abs=mh.tol_rel/1.0e4;
  if (debug) mh.verbose=1;
  if (function_verbose%10>2) mh.verbose=2;
  
  int mh_ret=mh.msolve(2,x1,sn_func);
  
  if (mh_ret!=0) {
    mh_ret=mh.msolve(2,x1,sn_func);
  }

  if (mh_ret!=0) {
    mh_ret=mh.msolve(2,x1,sn_func);
  }
  
  if (mh_ret!=0) {
    
    bool rbg_done=false;
    double qual_last=0.0;

    size_t k;
    for(k=0;k<10 && rbg_done==false;k++) {

      for(size_t j=0;j<2;j++) {

	funct ld_func=std::bind
	  (std::mem_fn<double(double,size_t,const ubvector &,
			      double,double,double,int,double &,double &,
			      thermo &)>
	   (&eos_nuclei::solve_nuclei_ld),this,
	   std::placeholders::_1,2,std::ref(x1),nB,Ye,T,j,
	   std::ref(mun_gas),std::ref(mup_gas),std::ref(th_gas));
	
	double blow=x1[j];
	double bhigh=x1[j]+0.001;
	double ylow=ld_func(blow);
	double yhigh=ld_func(bhigh);
	for(size_t i=0;i<40 && ylow*yhigh>0.0;i++) {
	  double diff=bhigh-blow;
	  blow-=diff*1.1;
	  bhigh+=diff*1.1;
	  ylow=ld_func(blow);
	  yhigh=ld_func(bhigh);
	  if (debug) {
	    cout << blow << " " << ylow << " "
		 << bhigh << " " << yhigh << endl;
	  }
	}
	
	if ((yhigh<0.0 && ylow>0.0) || (yhigh>0.0 && ylow<0.0)) {
	  if (debug) rbg.verbose=1;
	  if (function_verbose%10>2) rbg.verbose=2;
	  int rbg_ret=rbg.solve_bkt(blow,bhigh,ld_func);
	  if (rbg_ret!=0) {
	    return 2;
	  }
	} else {
	  return 3;
	}
	
	x1[j]=blow;
      }

      sn_func(2,x1,y1);
      
      double qual=sqrt(y1[0]*y1[0]+y1[1]*y1[1]);
      if (debug) {
	cout << "qual: " << qual << endl;
      }
      if (qual<1.0e-6) {
	rbg_done=true;
      } else {
	if (k==0) {
	  qual_last=qual;
	} else if (qual_last==qual) {
	  return 4;
	} else {
	  qual_last=qual;
	}
      }
    }

    if (k==10 && rbg_done==false) {
      return 5;
    }
    
  }

  // Perform a final call to solve_nuclei() to ensure
  // the chemical potentials are updated
  sn_func(2,x1,y1);
  
  if (debug) {
    sn_func(2,x1,y1);
    cout << "Success in eos_fixed_ZN():\n\t" << nucZ1 << " "
	 << nucN1 << " " << x1[0] << " " << x1[1] << " " << y1[0] << " "
	 << y1[1] << endl;
  }

  log_xn=x1[0];
  log_xp=x1[1];
  double xn=pow(10.0,x1[0]);
  double xp=pow(10.0,x1[1]);

  cout << "Fixme: nn and np not defined." << endl;
  exit(-1);
  double nn=1.0e10, np=1.0e10;
  
  double kappa=1.0-nB/n0;
  double xi=kappa+nn/n0+np/n0;
  
  double sum_nuc=0.0;
  double f_c=0.0, p_c2=0.0;
  ubvector fnuc(6), en(6), xx(6);
  f=0.0;
  thx.en=0.0;
  
  for (size_t i=0;i<6;i++) {
    if (nuclei[i].n>1.0e-300) {
      double lambda=sqrt(2.0*pi/nuclei[i].m/T);
      fnuc[i]=-T*(log(vomega[i]/nuclei[i].n/pow(lambda,3.0))+1.0)*
	nuclei[i].n;
      en[i]=nuclei[i].n*(log(vomega[i]/nuclei[i].n/pow(lambda,3.0))+
			 5.0/2.0+vomega_prime[i]*T/vomega[i]);
      if (!std::isfinite(fnuc[i])) {
	if (nuclei[i].n<1.0e-200) {
	  nuclei[i].n=0.0;
	  fnuc[i]=0.0;
	  en[i]=0.0;
	} else {
	  cout << "Nuclear free energy not finite in eos_fixed_ZN()." << endl;
	  cout << nuclei[i].n << " " << nuclei[i].be << " " << lambda
	       << " " << vomega[i] << endl;
	  exit(-1);
	}
      }
    } else {
      nuclei[i].n=0.0;
      fnuc[i]=0.0;
      en[i]=0.0;
    }
    fnuc[i]+=nuclei[i].n*nuclei[i].be;
    sum_nuc+=nuclei[i].n;
    f_c+=nuclei[i].n*Ec[i];
    f+=fnuc[i];
    thx.en+=en[i];
    // Coulomb contribution to pressure for nuclei
    xx[i]=cbrt(nB*Ye/n0*nuclei[i].A/nuclei[i].Z);
    double rA=cbrt(3.0*nuclei[i].A/4.0/pi/n0);
    p_c2+=-0.6*nuclei[i].n*nuclei[i].Z*nuclei[i].Z*
      fine_structure_f<double>()/rA*
      (0.5*xx[i]-0.5*pow(xx[i],3.0));
  }
  f+=f_c+xi*(th_gas.ed-T*th_gas.en)-T*sum_nuc*log(kappa);
  thx.en+=xi*th_gas.en+sum_nuc*log(kappa);
  thx.ed=f+T*thx.en;
  double p0_nuc=0.0;
  thx.pr=p0_nuc+sum_nuc*T/kappa+p_c2;

  if (!std::isfinite(f)) {
    cout << "Free energy not finite." << endl;
    vector_out(cout,fnuc,true);
    cout << log_xn << " " << log_xp << endl;
    cout << f_c << " " << th_gas.ed-T*th_gas.en << " " << kappa << endl;
    cout << f << endl;
    exit(-1);
    return 6;
  }
  
  mun_full=neutron.m;
  mup_full=proton.m;

  return 0;
}

int eos_nuclei::solve_hrg(size_t nv, const ubvector &x,
                          ubvector &y, double nB, double Ye, double T) {

  //exit(-1);
  if(false) {

  double nB2, Ye2;

  neutron.n=x[0];
  proton.n=x[1];

  thermo th;
  if (neutron.n<0.0 || proton.n<0.0) {
    return 1;
  }
  eosp_alt->calc_temp_e(neutron,proton,T,th);

  double nn_total=neutron.n, np_total=proton.n;

  double nn_tilde=neutron.n;
  double np_tilde=proton.n;

  double nn_tilde_old=nn_tilde;
  double np_tilde_old=np_tilde;
    
  nB2=0.0, Ye2=0.0;
    
  int iferm=0;
  int ibos=0;
  for(size_t j=0;j<part_db.size();j++) {
    if (part_db[j].id==2212 && part_db[j].charge==1) {
      if (iferm>=((int)res_f.size())) {
        O2SCL_ERR("Indexing problem with fermions.",o2scl::exc_efailed);
      }
      res_f[iferm]=proton;
      iferm++;
    } else if (part_db[j].id==2212 && part_db[j].charge==0) {
      if (iferm>=((int)res_f.size())) {
        O2SCL_ERR("Indexing problem with fermions.",o2scl::exc_efailed);
      }
      res_f[iferm]=neutron;
      iferm++;
    } else if (part_db[j].id==22 && part_db[j].charge==0) {
      if (ibos>=((int)res_b.size())) {
        O2SCL_ERR("Indexing problem with bosons.",o2scl::exc_efailed);
      }
      photon.massless_calc(T);
      res_b[ibos]=photon;
      ibos++;
    } else if (part_db[j].spin_deg%2==0) {
      if (iferm>=((int)res_f.size())) {
        O2SCL_ERR("Indexing problem with fermions.",o2scl::exc_efailed);
      }
      res_f[iferm].mu=part_db[j].baryon*(neutron.mu+neutron.m)+
        part_db[j].charge*(proton.mu+proton.m-neutron.mu-neutron.m);
      relf.calc_mu(res_f[iferm],T);
      if (j>700) {
        //cout << j << " " << part_db[j].baryon << " "
          //   << part_db[j].charge << " "
          //   << res_f[iferm].mu << " " << res_f[iferm].n << endl;
      }
      nB2+=part_db[j].baryon*res_f[iferm].n;
      Ye2+=part_db[j].charge*res_f[iferm].n;
      iferm++;
    } else if (part_db[j].id==-211 && part_db[j].charge==-1){
      cout <<"solve_nuclei: " << j << part_db[j].id << " " << part_db[j].name << endl;
      if (ibos>=((int)res_b.size())) {
        O2SCL_ERR("Indexing problem with bosons.",o2scl::exc_efailed);
      }
      cout << "pion mass: " << res_b[ibos].m*hc_mev_fm << endl;
      res_b[ibos].mu=part_db[j].baryon*(neutron.mu+neutron.m)+
        part_db[j].charge*(proton.mu+proton.m-neutron.mu-neutron.m);
      cout <<"Pion chem pot: " << res_b[ibos].mu*hc_mev_fm << endl;
      //cout << "Un: " << vdet["Un"] << " "
      //       << vdet_units.find("Un")->second << endl;
	    //cout << "Up: " << vdet["Up"] << " "
      //       << vdet_units.find("Up")->second << endl;
      //fr.calc_mu(res_b[ibos], neutron, proton, T, nB, vdet["Un"], vdet["Up"]);
      nB2+=part_db[j].baryon*res_b[ibos].n;
      Ye2+=part_db[j].charge*res_b[ibos].n;
      ibos++;
    } else {
      res_b[ibos].mu=part_db[j].baryon*(neutron.mu+neutron.m)+
        part_db[j].charge*(proton.mu+proton.m-neutron.mu-neutron.m);
      effb.calc_mu(res_b[ibos],T);
      nB2+=part_db[j].baryon*res_b[ibos].n;
      Ye2+=part_db[j].charge*res_b[ibos].n;
      ibos++;
    }
    if (!std::isfinite(nB2) ||
          !std::isfinite(Ye2)) {
        cout << j << " " << part_db[j].id << " "
             << nB2 << " " << Ye2 << endl;
        O2SCL_ERR("HRG problem.",o2scl::exc_einval);
    }
  }
  Ye2/=nB2;
  
  // I'm not 100% sure this is right
  nn_tilde+=nB2*(1.0-Ye2);
  np_tilde+=nB2*Ye2;
  
  //nB2=0.0;
  //Ye2=0.0;
  //nB2=neutron.n+proton.n+delta_pp.n;
  //Ye2=(proton.n-pi_minus.n+pi_plus.n)/nB2;
  cout << "HRG: mu_n: " << neutron.mu << ", mu_p: " << proton.mu << ", nn~_old: "
        << nn_tilde_old << endl;
  cout << "HRG: np~_old: " << np_tilde_old << ", nn~: "
        << nn_tilde << ", np~: " << np_tilde << endl;
  if (!std::isfinite(nn_tilde) ||
      !std::isfinite(np_tilde)) {
     O2SCL_ERR("HRG problem.",o2scl::exc_einval);
  }

  y[0]=(nB2-nB)/nB;
  y[1]=(Ye2-Ye)/Ye;

  /*
    cout << "HRG: " << nB2 << " " << Ye2 << " "
    << nB << " " << Ye << " "
    << neutron.n << " " << proton.n << " "
    << pi_minus.n << " " << pi_plus.n << " " << delta_pp.n << endl;
  */
  }
  return 0;
}

int eos_nuclei::nuc_matter(double nB, double Ye, double T,
			   double &log_xn, double &log_xp,
			   double &Zbar, double &Nbar, thermo &thx,
			   double &mun_full, double &mup_full,
			   int &A_min, int &A_max,
			   int &NmZ_min, int &NmZ_max,
			   map<string,double> &vdet) {

  if (include_muons) {
    
    mm_funct nm_func=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector&,double,double,
                       double,map<string,double> &)>
       (&eos_nuclei::nuc_matter_muons),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,nB,Ye,T,
       std::ref(vdet));
    
    ubvector x(2), y(2);
    
    x[0]=pow(10.0,log_xn)*nB;
    x[1]=pow(10.0,log_xp)*nB;

    int mret=mh.msolve(2,x,nm_func);
    if (mret!=0) {
      O2SCL_ERR("nuc matter muons failed.",o2scl::exc_einval);
    }

    // Ensure the last function evaluation stores the correct values
    // in 'vdet':
    nm_func(2,x,y);
    
    neutron.n=x[0];
    proton.n=x[1];
    
  } else {
    
    neutron.n=nB*(1.0-Ye);
    proton.n=nB*Ye;
    
  }    
  
  log_xn=log10(neutron.n/nB);
  log_xp=log10(proton.n/nB);
  
  Zbar=0.0;
  Nbar=0.0;
  
  if (alg_mode==4) {
    // If alg_mode is 4, then just set these to the fixed
    // distribution values
    A_min=5;
    A_max=fd_A_max;
    NmZ_min=-200;
    NmZ_max=200;
  } else {
    // Otherwise, these are the default guesses for
    // the "vary_dist" mode
    A_min=8;
    A_max=9;
    NmZ_min=-1;
    NmZ_max=1;
  }

  // Set all the nuclear densities to zero (this does appear
  // to be necessary not to confuse the other functions).
  for(size_t i=0;i<nuclei.size();i++) {
    nuclei[i].n=0.0;
  }

  if (inc_hrg) {

    mm_funct hrg_func=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector&,double,double,
                       double)>
       (&eos_nuclei::solve_hrg),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,nB,Ye,T);

    ubvector x(2), y(2);
    
    x[0]=nB*(1.0-Ye);
    x[1]=nB*Ye;

    mh.verbose=2;
    mh.msolve(2,x,hrg_func);
    
  } else {
    
    // Now compute the full EOS
    double f1, f2, f3, f4;
    //free_energy_density_detail(neutron,proton,T,thx,vdet);

    double nB2=0.0, Ye2=0.0;
    int iferm=0;
    int ibos=0;

    eos_sn_base eso;
    eso.verbose=2;
    thermo lep;
    elep.pair_density_eq(nB*Ye,T);
    //vdet["mue"]=electron.m;
    //cout << "nuc_matter::nB,Ye,T: " << nB << " " << Ye << " " << T << endl;
    eso.compute_eg_point(nB,Ye,T*hc_mev_fm,lep,electron.mu);
    cout << "try: " << elep.e.mu << endl;
    cout << "mue: " << electron.mu << endl;
    cout << "nuc_matter::mu_e, m_e: " << electron.mu << " " << electron.m
         << endl;

    for(size_t j=0;j<part_db.size();j++) {
      if (part_db[j].id==2212 && part_db[j].charge==1) {
        if (iferm>=((int)res_f.size())) {
          O2SCL_ERR("Indexing problem with fermions.",o2scl::exc_efailed);
        }
        res_f[iferm]=proton;
        iferm++;
      } else if (part_db[j].id==2212 && part_db[j].charge==0) {
        if (iferm>=((int)res_f.size())) {
          O2SCL_ERR("Indexing problem with fermions.",o2scl::exc_efailed);
        }
        res_f[iferm]=neutron;
        iferm++;
      } else if (part_db[j].id==22 && part_db[j].charge==0) {
        if (ibos>=((int)res_b.size())) {
          O2SCL_ERR("Indexing problem with bosons.",o2scl::exc_efailed);
        }
        photon.massless_calc(T);
        res_b[ibos]=photon;
        ibos++;
      } else if (part_db[j].spin_deg%2==0) {
        if (iferm>=((int)res_f.size())) {
          O2SCL_ERR("Indexing problem with fermions.",o2scl::exc_efailed);
        }
        res_f[iferm].mu=part_db[j].baryon*(neutron.mu+neutron.m)+
          part_db[j].charge*(proton.mu+proton.m-neutron.mu-neutron.m);
        relf.calc_mu(res_f[iferm],T);
        if (j>700) {
          //cout << j << " " << part_db[j].baryon << " "
            //   << part_db[j].charge << " "
            //   << res_f[iferm].mu << " " << res_f[iferm].n << endl;
        }
        nB2+=part_db[j].baryon*res_f[iferm].n;
        Ye2+=part_db[j].charge*res_f[iferm].n;
        iferm++;
      } else if (part_db[j].id==-211 && part_db[j].charge==-1){
        // Including pion calculations
        cout <<"nuc_matter: " << j << part_db[j].id << " " << part_db[j].name << endl;
        if (ibos>=((int)res_b.size())) {
          O2SCL_ERR("Indexing problem with bosons.",o2scl::exc_efailed);
        }
        cout << "pion mass: " << res_b[ibos].m*hc_mev_fm << endl;
        //res_b[ibos].mu=part_db[j].baryon*(neutron.mu+neutron.m)+
        //  part_db[j].charge*(proton.mu+proton.m-neutron.mu-neutron.m);
        res_b[ibos].mu=elep.e.mu;
        //cout <<"Pion chem pot: " << res_b[ibos].mu*hc_mev_fm << endl;
        funct function=std::bind
          (std::mem_fn<double(double, boson &, fermion &, 
           fermion &, fermion &, double, double, thermo &,
           map<string,double> &)>(&eos_nuclei::solve_pion),
           this,std::placeholders::_1,std::ref(res_b[ibos]),
           std::ref(neutron),std::ref(proton), 
           std::ref(elep.e),T,nB,std::ref(thx),std::ref(vdet));
        
        rbg.verbose=1;
        int ret=rbg.solve(proton.n, function);
        //cout <<"n_e: " << electron.n << ", n_mu: " << muon.n << ", n_pi: " << res_b[ibos].n << endl;
        //cout << "n_p: " << proton.n << ",n_tot: " << electron.n+muon.n+res_b[ibos].n << endl;
        nB2+=part_db[j].baryon*res_b[ibos].n;
        Ye2+=part_db[j].charge*res_b[ibos].n;
        ibos++;
      } else {
      res_b[ibos].mu=part_db[j].baryon*(neutron.mu+neutron.m)+
        part_db[j].charge*(proton.mu+proton.m-neutron.mu-neutron.m);
      effb.calc_mu(res_b[ibos],T);
      nB2+=part_db[j].baryon*res_b[ibos].n;
      Ye2+=part_db[j].charge*res_b[ibos].n;
      ibos++;
      }
      if (!std::isfinite(nB2) ||
          !std::isfinite(Ye2)) {
        cout << j << " " << part_db[j].id << " "
             << nB2 << " " << Ye2 << endl;
        O2SCL_ERR("HRG problem.",o2scl::exc_einval);
      }
    }

    mun_full=neutron.mu;
    mup_full=proton.mu;
    cout <<"nB2, Ye2: " << nB2 << " " << Ye2 << endl;
    
  }
  
  return 0;
}

double eos_nuclei::solve_pion(double x, boson &b, fermion &n, fermion &p, 
                              fermion &e, double T, double nB, thermo &thx, 
                              map<string,double> &vdet){
  p.n=x;
  n.n=nB-x;
  fr.verbose=2;
  free_energy_density_detail(n,p,T,thx,vdet);
  fr.calc_mu(b, n, p, T, nB);
  cout << "solve pion: n_p, n_n, n_e, n_pi: " << p.n  << " " << n.n << " " << e.n << " " << b.n << endl;
  return x-e.n-b.n;
}

int eos_nuclei::nuc_matter_muons(size_t nv, const ubvector &x, ubvector &y,
                                 double nB, double Ye, double T,
                                 map<string,double> &vdet) {

  neutron.n=x[0];
  proton.n=x[1];
  
  eos_sn_base eso;
  eso.include_muons=true;
  thermo lep;
  eso.compute_eg_point(nB,Ye,T*hc_mev_fm,lep,vdet["mue"]);

  y[0]=(neutron.n+proton.n)/nB-1.0;
  y[1]=(eso.electron.n+eso.muon.n-proton.n)/proton.n;
  vdet["Ymu"]=eso.muon.n/nB;
  
  return 0;
}

int eos_nuclei::eos_vary_dist
(double nB, double Ye, double T, double &log_xn, double &log_xp,
 double &Zbar, double &Nbar, thermo &thx, double &mun_full,
 double &mup_full, int &A_min, int &A_max, int &NmZ_min, int &NmZ_max,
 map<string,double> &vdet, bool dist_changed, bool no_nuclei) {
  
  int loc_verbose=function_verbose/1000%10;
  
  // We want homogeneous matter if we're above the saturation density
  // or if the 'no_nuclei' flag is true. 

  if (no_nuclei==true || nB>0.16) {
    nuc_matter(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,thx,mun_full,
	       mup_full,A_min,A_max,NmZ_min,NmZ_max,vdet);
    return 0;
  }
  
  int mpi_rank=0;
#ifndef NO_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
#endif
  
  bool done=false;
  bool expanded_A_min=false, expanded_A_max=false;
  bool expanded_NmZ_min=false, expanded_NmZ_max=false;

  if (alg_mode==4) {
    // If alg_mode is 4, then don't vary the distibution at all,
    // just choose these values
    A_min=5;
    A_max=fd_A_max;
    NmZ_min=-200;
    NmZ_max=200;
  } else if (A_min==0 || A_max==0 ||
	     NmZ_min==0 || NmZ_max==0) {
    A_min=8;
    A_max=9;
    NmZ_min=-1;
    NmZ_max=1;
  }    
  
  do {

    //cout << "X: " << A_min << " " << A_max << " "
    //<< NmZ_min << " " << NmZ_max << endl;

    if (A_min<4.5) {
      O2SCL_ERR2("Function eos_vary_dist() does not support ",
		 "A_min less than 4.5.",o2scl::exc_esanity);
    }

    // Compute the EOS with the current nuclear distribution
    int ret=eos_fixed_dist
      (nB,Ye,T,log_xn,log_xp,thx,mun_full,mup_full,
       A_min,A_max,NmZ_min,NmZ_max,vdet,dist_changed,no_nuclei);
    if (ret!=0) {
      return ret;
    }

    // Compute the contribution to the baryon density from the 
    // edges of the distribution
    
    double nB_A_min=0.0, nB_A_max=0.0;
    double nB_NmZ_min=0.0, nB_NmZ_max=0.0, nB_total=0.0;

    for(size_t index=5;index<nuclei.size();index++) {
	
      int A=nuclei[index].A;
      int N=nuclei[index].N;
      int Z=nuclei[index].Z;
      int NmZ=N-Z;
	
      double nB_temp=(N+Z)*nuclei[index].n;
      if (nB_temp<0.0) {
	cout << "A1: " << Z << " " << N << " " << index << " "
	     << nuclei[index].n << " " << nuclei.size() << endl;
	O2SCL_ERR2("Baryon density less than zero in ",
		   "vary_dist.",o2scl::exc_efailed);
      }
      nB_total+=nB_temp;
      if (A==A_min) {
	nB_A_min+=nB_temp;
      }
      if (A==A_max) {
	nB_A_max+=nB_temp;
      }
      if (NmZ==NmZ_min) {
	nB_NmZ_min+=nB_temp;
      }
      if (NmZ==NmZ_max) {
	nB_NmZ_max+=nB_temp;
      }
    }
    
    done=true;
    
    if (alg_mode==2) {

      // If the small A contribution is too large, then decrease A_min
      if (nB_A_min>nB*1.0e-6 && A_min>8) {
	A_min--;
	expanded_A_min=true;
	done=false;
      }
      // If the small A contribution is too small, then increase A_min
      if (expanded_A_min==false && nB_A_min<nB*1.0e-12 &&
	  A_min<400 && A_min<A_max) {
	A_min++;
	done=false;
      }
      
      // If the large A contribution is too large, then increase A_max
      if (nB_A_max>nB*1.0e-6 && A_max<399) {
	A_max++;
	expanded_A_max=true;
	done=false;
      }
      // If the large A contribution is too small, then decrease A_max
      if (expanded_A_max==false && nB_A_max<nB*1.0e-12 &&
	  A_max>8 && A_max>A_min) {
	A_max--;
	done=false;
      }

      // If the small N-Z contribution is too large, then decrease NmZ_min
      if (nB_NmZ_min>nB*1.0e-6) {
	NmZ_min--;
	expanded_NmZ_min=true;
	done=false;
      }
      // If the small N-Z contribution is too small, then increase NmZ_min
      if (expanded_NmZ_min==false && nB_NmZ_min<nB*1.0e-12 &&
	  NmZ_max-NmZ_min>1) {
	NmZ_min++;
	done=false;
      }
      
      // If the large N-Z contribution is too large, then increase NmZ_max
      if (nB_NmZ_max>nB*1.0e-6) {
	NmZ_max++;
	expanded_NmZ_max=true;
	done=false;
      }
      // If the large N-Z contribution is too small, then decrease NmZ_max
      if (expanded_NmZ_max==false && nB_NmZ_max<nB*1.0e-12 &&
	  NmZ_max-NmZ_min>1) {
	NmZ_max--;
	done=false;
      }
    }

    if (loc_verbose>1 && alg_mode==2) {
      cout << "Rank " << mpi_rank << " nB_A_min,max: " << nB_A_min
	   << " " << nB_A_max << " nB_NmZ_min,max: "
	   << nB_NmZ_min << " " << nB_NmZ_max << endl;
      cout << "Rank " << mpi_rank << "A_min,max: "
	   << A_min << " " << A_max << " NmZ_min,NmZ_max: "
	   << NmZ_min << " " << NmZ_max << " done: " << done << endl;
    }
    //char ch;
    //cin >> ch;
    
    if (done==false) dist_changed=true;

  } while (done==false);

  // Determine average N and Z
  
  double nt=0.0;
  Nbar=0.0;
  Zbar=0.0;
  for(size_t i=0;i<nuclei.size();i++) {
    nt+=nuclei[i].n;
    Zbar+=nuclei[i].Z*nuclei[i].n;
    Nbar+=nuclei[i].N*nuclei[i].n;
    if (loc_verbose>2) {
      cout << "vary_dist nucleus index,Z,N,n: " << i << " "
	   << nuclei[i].Z << " " << nuclei[i].N << " "
	   << nuclei[i].n << endl;
    }
  }
  Zbar/=nt;
  Nbar/=nt;
  
  // Compare with the homogeneous free energy and if it's favored,
  // then set the values accordingly

  if (nB>0.01) {
    
    thermo thx2;
    double log_xn2, log_xp2, Zbar2, Nbar2, mun_full2, mup_full2;
    int A_min2, A_max2, NmZ_min2, NmZ_max2;
    std::map<std::string,double> vdet2;

    nuc_matter(nB,Ye,T,log_xn2,log_xp2,Zbar2,Nbar2,thx2,mun_full2,
	       mup_full2,A_min2,A_max2,NmZ_min2,NmZ_max2,vdet2);

    if (loc_verbose>1) {
      cout << "Comparing: " << thx2.ed-T*thx2.en << " "
           << thx.ed-T*thx.en << endl;
    }
    if (thx2.ed-T*thx2.en<thx.ed-T*thx.en) {
      if (loc_verbose>1) {
        cout << "Nuclear matter preferred." << endl;
      }
      log_xn=log_xn2;
      log_xp=log_xp2;
      Zbar=Zbar2;
      Nbar=Nbar2;
      thx=thx2;
      mun_full=mun_full2;
      mup_full=mup_full2;
      A_min=A_min2;
      A_max=A_max2;
      NmZ_min=NmZ_min2;
      NmZ_max=NmZ_max2;
      vdet=vdet2;
      return 0;
    }

    // If nuclear matter was not preferred, then just go back
    // to the nuclei solution
    int ret=eos_fixed_dist
      (nB,Ye,T,log_xn,log_xp,thx,mun_full,mup_full,
       A_min,A_max,NmZ_min,NmZ_max,vdet,false,no_nuclei);
    if (ret!=0) {
      cout << "Fail X." << endl;
      exit(-1);
    }
    
  }

  if (loc_verbose>1) {
    cout << "Rank,size,Zbar,Nbar " << mpi_rank << " "
	 << nuclei.size() << " " << Zbar << " " << Nbar << endl;
  }

  if (loc_verbose==8) {
    char ch;
    cin >> ch;
  } else if (loc_verbose>=9) {
    exit(-1);
  }
  
  return 0;
}

int eos_nuclei::eos_fixed_dist
(double nB, double Ye, double T, double &log_xn, double &log_xp,
 thermo &thx, double &mun_full, double &mup_full, int &A_min,
 int &A_max, int &NmZ_min, int &NmZ_max,
 map<string,double> &vdet, bool dist_changed,
 bool no_nuclei) {

  int mpi_rank=0, mpi_size=1;
#ifndef NO_MPI
  // Get MPI rank, etc.
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
#endif
  
  int loc_verbose=function_verbose/100%10;
  double T_MeV=T*hc_mev_fm;
  double n0=0.16;
  // Chemical potentials for homogeneous matter
  double mun_gas, mup_gas;
  thermo th_gas;

  int n_solves=(fixed_dist_alg%10)+1;
  int n_brackets=((fixed_dist_alg/10)%10)*10;
  int n_minimizes=(fixed_dist_alg/100)%10*4;
  int n_randoms=((fixed_dist_alg/1000)%10)*1000;
  if (mpi_size==1 && loc_verbose>1) {
    cout << "n_solves,n_brackets,n_minimizes,n_randoms: "
	 << n_solves << " " << n_brackets << " " << n_minimizes << " "
	 << n_randoms << endl;
  }

  // ----------------------------------------------------------------
  // If the distribution has changed, set up the nuclei array and the
  // neutron separation energies

  if (dist_changed) {

    // Count nuclei
    size_t count=5;
    for(int A=A_min;A<=A_max;A++) {
      for(int NmZ=NmZ_min;NmZ<=NmZ_max;NmZ+=2) {
	int N=(A+NmZ)/2;
	int Z=A-N;
	if (N>2 && Z>2 && (rnuc_less_rws==false || Z>=nB*Ye*A/n0) &&
	    N<=max_ratio*Z && Z<=max_ratio*N) {
	  if (ame.is_included(Z,N) &&
	      ame.is_included(Z-1,N) &&
	      ame.is_included(Z,N-1)) {
            count++;
	  } else if (m95.is_included(Z,N) &&
		     m95.is_included(Z-1,N) &&
		     m95.is_included(Z,N-1)) {
            count++;
	  } else if (extend_frdm) {
            count++;
          }
	}
      }
    }

    if (count<5) {
      O2SCL_ERR2("Size of nuclei vector less than 5 in ",
		 "eos_nuclei::eos_fixed_dist().",o2scl::exc_esanity);
    }
    
    // Resize nuclei array
    nuclei.resize(count);

    // Set up the light nuclei
    nuc_alpha=&nuclei[0];
    nuc_deut=&nuclei[1];
    nuc_trit=&nuclei[2];
    nuc_he3=&nuclei[3];
    nuc_li4=&nuclei[4];
    
    ame.get_nucleus(2,2,*nuc_alpha);
    ame.get_nucleus(1,1,*nuc_deut);
    ame.get_nucleus(1,2,*nuc_trit);
    ame.get_nucleus(2,1,*nuc_he3);
    ame.get_nucleus(3,1,*nuc_li4);
    
    nuc_alpha->g=1.0;
    nuc_deut->g=3.0;
    nuc_trit->g=2.0;
    nuc_li4->g=5.0;
    nuc_he3->g=2.0;
    
    // Resize all of the vectors
    Sneut.resize(count);
    Sprot.resize(count);
    vomega.resize(count);
    vomega_prime.resize(count);
    Ec.resize(count);
    
    // Set the neutron and proton separation energies for light nuclei
    nucleus nuc_temp;
    for(size_t i=0;i<5;i++) {
      
      if (ame.is_included(nuclei[i].Z,nuclei[i].N-1) &&
          ame.is_included(nuclei[i].Z-1,nuclei[i].N)) {
        ame.get_nucleus(nuclei[i].Z,nuclei[i].N-1,nuc_temp);
        Sneut[i]=-(nuclei[i].be-nuc_temp.be)*hc_mev_fm;
        ame.get_nucleus(nuclei[i].Z-1,nuclei[i].N,nuc_temp);
        Sprot[i]=-(nuclei[i].be-nuc_temp.be)*hc_mev_fm;
      } else {
        Sneut[i]=-1.0;
        Sprot[i]=-1.0;
      }
      
    }
    
    // Set the neutron and proton separation energies and
    // spin degeneracies for remaining nuclei
    size_t index=5;
    for(int A=A_min;A<=A_max;A++) {
      for(int NmZ=NmZ_min;NmZ<=NmZ_max;NmZ+=2) {
	int N=(A+NmZ)/2;
	int Z=A-N;
        
	if (N>2 && Z>2 && (rnuc_less_rws==false || Z>=nB*Ye*A/n0) &&
	    N<=max_ratio*Z && Z<=max_ratio*N) {
          
          bool skip_nucleus=false;
          
	  // Obtain heavy nucleus and neutron and proton separation
	  // energies. Ensuring that all three nuclei are in the same
	  // mass formula avoids problems with the neutron and
	  // separation energies at the boundary between the mass
	  // formulae.

	  if (ame.is_included(Z,N) &&
	      ame.is_included(Z-1,N) &&
	      ame.is_included(Z,N-1)) {
	    
	    ame.get_nucleus(Z,N,nuclei[index]);
	    ame.get_nucleus(Z,N-1,nuc_temp);
	    Sneut[index]=-(nuclei[index].be-nuc_temp.be)*hc_mev_fm;
	    ame.get_nucleus(Z-1,N,nuc_temp);
	    Sprot[index]=-(nuclei[index].be-nuc_temp.be)*hc_mev_fm;
	    
	  } else if (m95.is_included(Z,N) &&
		     m95.is_included(Z-1,N) &&
		     m95.is_included(Z,N-1)) {
	    
	    m95.get_nucleus(Z,N,nuclei[index]);
	    m95.get_nucleus(Z,N-1,nuc_temp);
	    Sneut[index]=-(nuclei[index].be-nuc_temp.be)*hc_mev_fm;
	    m95.get_nucleus(Z-1,N,nuc_temp);
	    Sprot[index]=-(nuclei[index].be-nuc_temp.be)*hc_mev_fm;

	  } else if (extend_frdm) {
	    
	    frdm.get_nucleus(Z,N,nuclei[index]);
	    frdm.get_nucleus(Z,N-1,nuc_temp);
	    Sneut[index]=-(nuclei[index].be-nuc_temp.be)*hc_mev_fm;
	    frdm.get_nucleus(Z-1,N,nuc_temp);
	    Sprot[index]=-(nuclei[index].be-nuc_temp.be)*hc_mev_fm;
	    
	  } else {
            
            skip_nucleus=true;
            
          }            

          if (skip_nucleus==false) {

            if (index>=nuclei.size()) {
              std::cout << "Index: " << index << " nuclei.size(): "
                        << nuclei.size() << std::endl;
              O2SCL_ERR2("Indexing problem in ",
                         "eos_nuclei::eos_fixed_dist().",o2scl::exc_esanity);
            }
            
            // Set spin degeneracy of heavy nucleus
            if (hfb.is_included(Z,N)) {
              if (hfb.get_ZN(Z,N).Jexp<99) {
                nuclei[index].g=2.0*hfb.get_ZN(Z,N).Jexp+1.0;
              } else {
                nuclei[index].g=2.0*hfb.get_ZN(Z,N).Jth+1.0;
              }
            } else {
              if (Z%2==0 && N%2==0) {
                nuclei[index].g=1.0;
                //} else if (Z%2==1 && N%2==1) {
                //nuclei[index].g=3.0;
              } else {
                nuclei[index].g=2.0;
              }
            }
            
            index++;
          }
	}
      }
    }
    
    // End of if dist_changed
  }

  // ----------------------------------------------------------------
  // Compute the partition functions for this distribution and this
  // temperature
  
  size_t n_nuclei=nuclei.size();

  // Set up integrand for partition function
  part_func.T_MeV=T_MeV;
  funct f1=std::bind(std::mem_fn<double(double)>
		     (&partition_func::delta_small_iand),&part_func,
		     std::placeholders::_1);
  funct f1_prime=std::bind(std::mem_fn<double(double)>
			   (&partition_func::delta_small_iand_prime),
			   &part_func,std::placeholders::_1);
  funct f2=std::bind(std::mem_fn<double(double)>
		     (&partition_func::delta_large_iand),&part_func,
		     std::placeholders::_1);
  funct f2_prime=std::bind(std::mem_fn<double(double)>
			   (&partition_func::delta_large_iand_prime),
			   &part_func,std::placeholders::_1);

  
  double res, err, ret, ret_prime, res_prime, err_prime;
  
  for (size_t i=0;i<n_nuclei;i++) {

    if (i==1) {
      
      vomega[i]=3.0;
      vomega_prime[i]=0.0;
      
    } else {
      
      double zEd=min(Sneut[i],Sprot[i])/2.0;
      double zR=1.25*cbrt(nuclei[i].Z+nuclei[i].N-1.0);
      double zER=0.5/nuclei[i].m/zR/zR*hc_mev_fm;
      double zEc=(nuclei[i].Z-1.0)*fine_structure_f<double>()/zR*hc_mev_fm;
      double zEt=min(Sneut[i]+zER,Sprot[i]+zER+zEc/2.0);
      double zrA=cbrt(3.0*(nuclei[i].N+nuclei[i].Z)/4.0/pi/n0);
      double delta_p=11.0/sqrt(nuclei[i].Z+nuclei[i].N)*
	(1.0+pow(-1.0,nuclei[i].Z)/2.0+pow(-1.0,nuclei[i].N)/2.0);
      if (nuclei[i].Z<=30) {
	part_func.a=0.052*pow(nuclei[i].N+nuclei[i].Z,1.2);
	part_func.delta=delta_p-80.0/(nuclei[i].Z+nuclei[i].N);
      } else {
	part_func.a=0.125*(nuclei[i].N+nuclei[i].Z);
	part_func.delta=delta_p-80.0/(nuclei[i].Z+nuclei[i].N)-0.5;
      }
      if (!std::isfinite(part_func.a)) {
	cout << nuclei[i].N << " " << nuclei[i].Z << " " << delta_p << " "
	     << part_func.a << endl;
	O2SCL_ERR2("Variable a not finite in ",
		   " eos_nuclei::eos_fixed_dist()",
		   o2scl::exc_esanity);
      }
      if (!std::isfinite(part_func.delta)) {
	cout << nuclei[i].N << " " << nuclei[i].Z << " "
	     << part_func.delta << endl;
        cout << endl;
        cout << delta_p << endl;
        cout << endl;
	O2SCL_ERR2("Variable delta not finite",
		   " eos_nuclei::eos_fixed_dist()",
		   o2scl::exc_esanity);
      }
      
      if (Sneut[i]<0.0 || Sprot[i]<0.0) {
	
	vomega[i]=nuclei[i].g;
	vomega_prime[i]=0.0;
	
      } else {
	
	if (part_func.delta>zEd) {
	  part_func.delta=zEd;
	  part_func.Tc=1.0/(-1.25/part_func.delta+sqrt(part_func.a)/
			    sqrt(part_func.delta));
	  part_func.C=sqrt(pi)/12.0*pow(part_func.a,-0.25)*
	    pow(part_func.delta,-1.25)*
	    exp(1.25+sqrt(part_func.a*part_func.delta));
	  if(2.0*zEd<zEt) {
	    if (!std::isfinite(zEd)) {
	      O2SCL_ERR2("Variable zEd not finite",
			 "eos_nuclei::eos_fixed_dist()",
			 o2scl::exc_esanity);
	    }
	    if (!std::isfinite(zEt)) {
	      O2SCL_ERR2("Variable zEt not finite",
			 "eos_nuclei::eos_fixed_dist()",
			 o2scl::exc_esanity);
	    }
	    ret=iqg.integ_err(f1,2.0*zEd,zEt,res,err)+
	      iqg.integ_err(f2,zEd,2.0*zEd,res,err);
	    ret_prime=iqg.integ_err(f1_prime,2.0*zEd,zEt,res_prime,err_prime)+
	      iqg.integ_err(f2_prime,zEd,2.0*zEd,res_prime,err_prime);
	  } else {
	    ret=iqg.integ_err(f2,zEd,zEt,res,err);
	    ret_prime=iqg.integ_err(f2_prime,zEd,zEt,res_prime,err_prime);
	  }
	} else {
	  if (!std::isfinite(zEd)) {
	    O2SCL_ERR2("Variable zEd not finite",
		       "eos_nuclei::eos_fixed_dist()",
		       o2scl::exc_esanity);
	  }
	  if (!std::isfinite(zEt)) {
	    O2SCL_ERR2("Variable zEt not finite",
		       "eos_nuclei::eos_fixed_dist()",
		       o2scl::exc_esanity);
	  }
	  ret=iqg.integ_err(f1,zEd,zEt,res,err);
	  ret_prime=iqg.integ_err(f1_prime,zEd,zEt,res_prime,err_prime);
	}
	vomega[i]=nuclei[i].g+res;
	vomega_prime[i]=res_prime/T;
      }
      
    }

    if (false) {
      double T_K=o2scl_settings.get_convert_units().convert
        ("MeV","K",T_MeV);
      
      double v, vop;
      pfuncs.few78(nuclei[i].Z,nuclei[i].N,T_K,v,vop);
      vop/=T;
      
      if (fabs(vomega[i]-v)/fabs(v)>1.0e-5 ||
          fabs(vomega_prime[i]-vop)/fabs(vop)>1.0e-5) {
        cout << nuclei[i].Z << " " << nuclei[i].N << " "
             << vomega[i] << " " << vomega_prime[i] << " ";
        cout << v << " " << vop << endl;
        exit(-1);
      }
    }
    
  }

  // ---------------------------------------------------------------
  // Set up for calling the solver
  
  ubvector x1(2), y1(2);

  mm_funct sn_func=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector&,double,double,
		     double,int,double &,double &,thermo &,
		     map<string,double> &)>
     (&eos_nuclei::solve_nuclei),this,std::placeholders::_1,
     std::placeholders::_2,std::placeholders::_3,nB,Ye,T,0,
     std::ref(mun_gas),std::ref(mup_gas),std::ref(th_gas),
     std::ref(vdet));

  x1[0]=log_xn;
  x1[1]=log_xp;

  if (mpi_size==1 && loc_verbose>=2) {
    cout << "Initial guess (log_xn,log_xp): "
         << x1[0] << " " << x1[1] << endl;
  }

  if ((alg_mode==2 || alg_mode==4) && nB<1.0e-11) {
    mh.tol_rel=mh_tol_rel/1.0e2;
  } else {
    mh.tol_rel=mh_tol_rel;
  }

  // 8/27/19: updated
  // 8/13/2020 AWS: I don't think this does anything because
  // the mroot_hybrids solver ignores tol_abs. 
  if (alg_mode==2 || alg_mode==4 || nB<1.0e-09) {
    mh.tol_abs=mh.tol_rel/1.0e4;
  }

  int mh_ret=-1;
  
  if (verify_only) {
    sn_func(2,x1,y1);
    if (fabs(y1[0])+fabs(y1[1])<mh.tol_rel) {
      mh_ret=0;
    } else {
      cout << "Verify failed." << endl;
      cout << "  nB,Ye,T[MeV]: " << nB << " " << Ye << " "
	   << T*hc_mev_fm << endl;
      cout << "  y[0],y[1],tol: " << y1[0] << " " << y1[1] << " "
	   << mh.tol_rel << endl;
      return 1;
    }
  }
  
  if (loc_verbose>2) mh.verbose=1;

  double qual_best=1.0e100;

  if (mh_ret!=0) {
    cout << "x1 size: " << x1.size() << endl;
    mh_ret=mh.msolve(2,x1,sn_func);
    if (mh_ret==0 && mpi_size==1) {
      if (loc_verbose>1) {
        cout << "Rank " << mpi_rank
             << " finished after initial solve." << endl;
      }
    }
  }

  if (alg_mode==2 || alg_mode==4) {

    for(int k=0;k<n_solves && mh_ret!=0;k++) {
      mh_ret=mh.msolve(2,x1,sn_func);
      int iret=sn_func(2,x1,y1);
      if (iret==0 && fabs(y1[0])+fabs(y1[1])<qual_best) {
	qual_best=fabs(y1[0])+fabs(y1[1]);
      }
      if (loc_verbose>=2 && mpi_size==1) {
	cout << "Rank " << mpi_rank << " solve " << k+1
	     << "/" << n_solves << " x1[0],x1[1],qual: ";
        cout.precision(5);
        cout << x1[0] << " " << x1[1] << " "
	     << qual_best << endl;
        cout.precision(6);
      }
      if (mh_ret==0 && mpi_size==1 && loc_verbose>0) {
	cout << "Rank " << mpi_rank << " finished after solve " << k
	     << "." << endl;
      }
    }
    
  } else {
    double count=0;
    double count_total=5000;
    if(nB<2.2e-12) {
      count_total=200;
    }
    if(nB<0.08 && nB>0.01 && T*197.33<10.0) {
      while((mh_ret!=0||((x1[0]>log_xn*0.9||x1[1]>log_xp*0.9)&& 
			 pow(10.0,x1[0])+pow(10.0,x1[1])>0.50))
	    && count<count_total) {
        mh.tol_abs=mh.tol_rel/1.0e4;
        x1[0]=log_xn*(1.0+1.0*rng.random_int(100)*0.01);
        x1[1]=log_xp*(1.0+1.0*rng.random_int(100)*0.01);
        mh_ret=mh.msolve(2,x1,sn_func);
        count++;
      }
    } else {
      while(mh_ret!=0 && count<count_total) {

        mh.tol_abs=mh.tol_rel/1.0e4;
        x1[0]=log_xn*(1.0+0.5*rng.random_int(100)*0.01);
        x1[1]=log_xp*(1.0+0.5*rng.random_int(100)*0.01);
        mh_ret=mh.msolve(2,x1,sn_func);
        count++;
      }
    }
  }

  // ---------------------------------------------------------------
  // If the solver failed, call the bracketing solver

  if (mh_ret!=0 && (nB<2.2e-12 || alg_mode==2 || alg_mode==4)) {
    
    bool bracketing_done=false;
    double qual_last=0.0;
    vector<double> ranges={1.0e99,1.0e99,1.0e99,1.0e99};

    ubvector x1_old(2);
    x1_old[0]=x1[0];
    x1_old[1]=x1[1];
    
    // Try 10 iterations solving one variable at a time
    int k;
    for(k=0;k<n_brackets && bracketing_done==false;k++) {

      // j=0 (1) for neutrons (protons)
      for(size_t j=0;j<2;j++) {

	// First, find a suitable bracket
	funct ld_func=std::bind
	  (std::mem_fn<double(double,size_t,const ubvector &,
			      double,double,double,int,double &,double &,
			      thermo &)>
	   (&eos_nuclei::solve_nuclei_ld),this,
	   std::placeholders::_1,2,std::ref(x1),nB,Ye,T,j,
	   std::ref(mun_gas),std::ref(mup_gas),std::ref(th_gas));
	
	double blow=x1[j];
	double bhigh=x1[j]+0.001;
	double ylow=ld_func(blow);
	double yhigh=ld_func(bhigh);
	for(size_t i=0;i<40 && ylow*yhigh>0.0;i++) {
	  double diff=bhigh-blow;
	  blow-=diff*1.1;
	  bhigh+=diff*1.1;
	  ylow=ld_func(blow);
	  yhigh=ld_func(bhigh);
	  
	  if (fabs(ylow)<1.0e70 && fabs(yhigh)<1.0e70) {
	    // Store the bracketing ranges for later use
	    ranges[j*2]=blow;
	    ranges[j*2+1]=bhigh;
	  }
	  if (loc_verbose>2) {
	    cout << blow << " " << ylow << " "
		 << bhigh << " " << yhigh << endl;
	  }
	}
	
	if ((yhigh<0.0 && ylow>0.0) || (yhigh>0.0 && ylow<0.0)) {
	  // If we've found a good bracket, then use
	  // Brent's method.
	  if (loc_verbose>2) rbg.verbose=1;
	  int rbg_ret=rbg.solve_bkt(blow,bhigh,ld_func);
	  if (rbg_ret!=0) {
	    // If the bracketing solver failed, then stop
	    j=2;
	    bracketing_done=true;
	  } else {
	    // Update the current root
	    x1[j]=blow;
	  }
	} else {
	  // We couldn't find a bracket, so just stop
	  j=2;
	  bracketing_done=true;
	}

	// End of for (j=0;j<2;j++) {
      }

      if (loc_verbose>=2 && mpi_size==1) {
	cout << "Rank " << mpi_rank << " bracket " << k+1
	     << "/" << n_brackets << " x1[0],x1[1]: ";
        cout.precision(5);
        cout << x1[0] << " " << x1[1] << endl;
        cout.precision(6);
      }
      
      if (bracketing_done==false) {

	// Check the quality of the current root
	int iret=sn_func(2,x1,y1);

	if (iret==0 && fabs(y1[0])+fabs(y1[1])<qual_best) {
	  qual_best=fabs(y1[0])+fabs(y1[1]);
	}
	  
	double qual=fabs(y1[0])+fabs(y1[1]);
	if (loc_verbose>2) {
	  cout << "qual: " << qual << " mh.tol_rel: "
	       << mh.tol_rel << endl;
	}
	
	if (iret==0 && qual<mh.tol_rel) {
	  // Stop, and set mh_ret to zero to indicate that we're done
	  k=n_brackets;
	  mh_ret=0;
	  if (mpi_size==1 && loc_verbose>0) {
	    cout << "Rank " << mpi_rank << " finished after "
		 << k << " brackets." << endl;
	  }
	} else if (k>0 && qual_last==qual) {
	  // If our quality hasn't improved, stop
	  k=10;
          bracketing_done=true;
	} else {
	  // Update qual_last
	  qual_last=qual;
	}
	
      }
      
      // End of for (k=0;k<10;k++) {
    }

    {
      // Sometimes, the bracketing fails and results in a point
      // which returns an error, so in that case we restore the
      // best guess before bracketing
      int yret=sn_func(2,x1,y1);
      if (yret!=0) {
	x1[0]=x1_old[0];
	x1[1]=x1_old[1];
      }
    }

    // ---------------------------------------------------------------
    // Adaptive differential evolution
    
    if (false) {

      multi_funct min_func=std::bind
	(std::mem_fn<double(size_t,const ubvector &,
			    double,double,double,double &,double &,
			    thermo &)>
	 (&eos_nuclei::solve_nuclei_min),this,
	 std::placeholders::_1,std::placeholders::_2,nB,Ye,T,
	 std::ref(mun_gas),std::ref(mup_gas),std::ref(th_gas));
      
      if (fd_rand_ranges.size()>=4) {
	ranges[0]=fd_rand_ranges[0];
	ranges[1]=fd_rand_ranges[1];
	ranges[2]=fd_rand_ranges[2];
	ranges[3]=fd_rand_ranges[3];
      }
      
      // If the ranges don't include the best point so far, expand them
      if (x1[0]<ranges[0]) ranges[0]=x1[0]-(ranges[1]-ranges[0]);
      if (x1[0]>ranges[1]) ranges[1]=x1[0]+(ranges[1]-ranges[0]);
      if (x1[1]<ranges[2]) ranges[2]=x1[1]-(ranges[3]-ranges[2]);
      if (x1[1]>ranges[3]) ranges[3]=x1[1]+(ranges[3]-ranges[2]);

      if (ranges[1]>0.0) ranges[1]=0.0;
      if (ranges[3]>0.0) ranges[3]=0.0;

      diff_evo_adapt<> de;
      de.err_nonconv=false;
      ubvector step(2);
      step[0]=ranges[1]-ranges[0];
      step[1]=ranges[3]-ranges[2];
      double ymin;
      //de.verbose=1;
      int de_ret=de.mmin(2,x1,ymin,min_func);

      if (de_ret==0) {
	int iret=sn_func(2,x1,y1);
	
	if (iret==0 && fabs(y1[0])+fabs(y1[1])<qual_best) {
	  qual_best=fabs(y1[0])+fabs(y1[1]);
	}
	
	if (qual_best<mh.tol_rel) {
	  mh_ret=0;
	}
      }
      
    }
    
    // ---------------------------------------------------------------
    // Try the minimizer

    if (n_minimizes>0 && mh_ret!=0) {
      
      multi_funct min_func=std::bind
	(std::mem_fn<double(size_t,const ubvector &,
			    double,double,double,double &,double &,
			    thermo &)>
	 (&eos_nuclei::solve_nuclei_min),this,
	 std::placeholders::_1,std::placeholders::_2,nB,Ye,T,
	 std::ref(mun_gas),std::ref(mup_gas),std::ref(th_gas));
      
      mmin_simp2<> ms;
      ms.ntrial*=10;
      double ymin;
      ms.err_nonconv=false;
      //ms.verbose=1;
      ms.tol_abs/=1.0e4;
      ubvector step(1);
      step[0]=0.01;
      ms.set_step(1,step);
      int mret=1;
      int jk;
      for(jk=0;jk<n_minimizes && qual_best>mh.tol_rel;jk++) {
	
	mret=ms.mmin(2,x1,ymin,min_func);

	// AWS 9/18/2020: This doesn't seem to help
	if (false) {
	  mh_ret=mh.msolve(2,x1,sn_func);
	}
	
	int iret=sn_func(2,x1,y1);

	if (iret==0 && fabs(y1[0])+fabs(y1[1])<qual_best) {
	  qual_best=fabs(y1[0])+fabs(y1[1]);
	}

	if (loc_verbose>=2 && mpi_size==1) {
	  cout << "Rank " << mpi_rank << " min " << jk+1
	       << "/" << n_minimizes << " ret,qual,x1[0],x1[1]: ";
          cout.precision(5);
          cout << mret << " " << qual_best
	       << " " << x1[0] << " " << x1[1] << endl;
          cout.precision(6);
	}
	
      }

      // Evaluate the equations and see if they have been
      // solved, independent of the value of mret
      int iret2=sn_func(2,x1,y1);
      
      double qual=fabs(y1[0])+fabs(y1[1]);
      if (loc_verbose>2) {
	cout << "qual: " << qual << " " << y1[0] << " " << y1[1] << " "
	     << mret << endl;
      }
      
      if (iret2==0 && qual<mh.tol_rel) {
	// Set mh_ret to zero to indicate that we're done
	mh_ret=0;
	if (mpi_size==1 && loc_verbose>0) {
	  cout << "Rank " << mpi_rank << " finished after "
	       << jk << " minimizer calls." << endl;
	}
      }
      
    }
    
    // ---------------------------------------------------------------
    // Set up the solution ranges
    
    if (mh_ret!=0) {
      
      if (fd_rand_ranges.size()>=4) {
	ranges[0]=fd_rand_ranges[0];
	ranges[1]=fd_rand_ranges[1];
	ranges[2]=fd_rand_ranges[2];
	ranges[3]=fd_rand_ranges[3];
      }
      
      // If the ranges don't include the best point so far, expand them
      if (x1[0]<ranges[0]) ranges[0]=x1[0]-(ranges[1]-ranges[0]);
      if (x1[0]>ranges[1]) ranges[1]=x1[0]+(ranges[1]-ranges[0]);
      if (x1[1]<ranges[2]) ranges[2]=x1[1]-(ranges[3]-ranges[2]);
      if (x1[1]>ranges[3]) ranges[3]=x1[1]+(ranges[3]-ranges[2]);

      if (ranges[1]>0.0) ranges[1]=0.0;
      if (ranges[3]>0.0) ranges[3]=0.0;
      
      if (mpi_size==1 && loc_verbose>1) {
	cout << "x1,ranges: ";
        cout << x1[0] << " " << x1[1] << " "
	     << ranges[0] << " " << ranges[1] << "\n  "
	     << ranges[2] << " " << ranges[3] << endl;
      }
      
    }
    
    // ---------------------------------------------------------------
    // If the minimizer didn't work, try random initial guesses

    if (mh_ret!=0 && ranges[0]<1.0e50 && ranges[1]<1.0e50 &&
	ranges[2]<1.0e50 && ranges[3]<1.0e50) {
      for(int kk=0;kk<n_randoms && mh_ret!=0;kk++) {
	x1[0]=ranges[0]+rng.random()*(ranges[1]-ranges[0]);
	x1[1]=ranges[2]+rng.random()*(ranges[3]-ranges[2]);
	mh_ret=mh.msolve(2,x1,sn_func);
	if (loc_verbose>2) {
	  cout << kk << " " << x1[0] << " " << x1[1] << " "
	       << mh_ret << endl;
	}
	int iret=sn_func(2,x1,y1);
	
	if (iret==0 && fabs(y1[0])+fabs(y1[1])<qual_best) {
	  qual_best=fabs(y1[0])+fabs(y1[1]);
	}
	
	if (loc_verbose>=2 && mpi_size==1 && kk%100==99) {
	  cout << "Rank " << mpi_rank << " solve " << kk+1
	       << "/" << n_randoms << " x1[0],x1[1],qual: ";
          cout.precision(4);
          cout << x1[0] << " " << x1[1] << " "
	       << qual_best << endl;
          cout.precision(6);
	}
	if (mh_ret==0 && mpi_size==1 && loc_verbose>0) {
	  cout << "Rank " << mpi_rank << " finished after "
	       << kk << " random solves." << endl;
	}
	//cout << kk << " " << mh_ret << endl;

	// AWS 8/31/2020: I'm not sure if this helps or not
	if (false && kk%1000==999) {
	  double d1=ranges[1]-ranges[0];
	  double d2=ranges[3]-ranges[2];
	  ranges[0]-=d1;
	  ranges[1]+=d1;
	  ranges[2]-=d2;
	  ranges[3]+=d2;
	  if (ranges[1]>0.0) ranges[1]=0.0;
	  if (ranges[3]>0.0) ranges[3]=0.0;
	  if (mpi_size==1 && loc_verbose>1) {
	    cout << "ranges: "
		 << ranges[0] << " " << ranges[1] << " "
		 << ranges[2] << " " << ranges[3] << endl;
	  }
	}

	//cout << "mh_ret: " << mh_ret << endl;
	//char ch;
	//cin >> ch;
	
      }
    }

    // ---------------------------------------------------------------
    // Survey the solution space

    if (survey_eqs) {

      static const size_t surv_N=100;
      
      vector<double> g_forward, g_backward;
      for(double xt=1.0e-8;xt<0.01001;xt*=pow(1.0e6,1.0/((double)surv_N-1))) {
        g_forward.push_back(xt);
        g_backward.push_back(xt);
      }
      vector_reverse<vector<double>,double>(surv_N,g_backward);
      
      // Divide the space into four quadrants, ++, +-, -+, and --
      table3d surv_pp;
      surv_pp.set_xy("dx",surv_N,g_forward,"dy",surv_N,g_forward);
      surv_pp.line_of_names("e1 e2 s ret");
      table3d surv_pm;
      surv_pm.set_xy("dx",surv_N,g_forward,"dy",surv_N,g_backward);
      surv_pm.line_of_names("e1 e2 s ret");
      table3d surv_mm;
      surv_mm.set_xy("dx",surv_N,g_backward,"dy",surv_N,g_backward);
      surv_mm.line_of_names("e1 e2 s ret");
      table3d surv_mp;
      surv_mp.set_xy("dx",surv_N,g_backward,"dy",surv_N,g_forward);
      surv_mp.line_of_names("e1 e2 s ret");
      
      ubvector x2(2);
      
      for(size_t sj=0;sj<surv_N;sj++) {
        cout << sj << "/" << surv_N << endl;
	for(size_t suk=0;suk<surv_N;suk++) {

          // Handle the ++ quadrant
	  x2[0]=x1[0]+surv_pp.get_grid_x(sj);
	  x2[1]=x1[1]+surv_pp.get_grid_y(suk);
	  int siret=sn_func(2,x2,y1);
          
	  //cout << "Survey: " << sj << " " << suk << " ";
	  //cout << siret << " " << y1[0] << " " << y1[1] << " ";
          
	  if (siret==0) {
	    surv_pp.set(sj,suk,"e1",y1[0]);
	    surv_pp.set(sj,suk,"e2",y1[1]);
	    surv_pp.set(sj,suk,"s",fabs(y1[0])+fabs(y1[1]));
	  } else {
	    surv_pp.set(sj,suk,"e1",0.0);
	    surv_pp.set(sj,suk,"e2",0.0);
	    surv_pp.set(sj,suk,"s",0.0);
	  }
	  surv_pp.set(sj,suk,"ret",siret);

          // Handle the +- quadrant
	  x2[0]=x1[0]+surv_pm.get_grid_x(sj);
	  x2[1]=x1[1]-surv_pm.get_grid_y(suk);
	  siret=sn_func(2,x2,y1);
          
	  //cout << siret << " " << y1[0] << " " << y1[1] << " ";
          
	  if (siret==0) {
	    surv_pm.set(sj,suk,"e1",y1[0]);
	    surv_pm.set(sj,suk,"e2",y1[1]);
	    surv_pm.set(sj,suk,"s",fabs(y1[0])+fabs(y1[1]));
	  } else {
	    surv_pm.set(sj,suk,"e1",0.0);
	    surv_pm.set(sj,suk,"e2",0.0);
	    surv_pm.set(sj,suk,"s",0.0);
	  }
	  surv_pm.set(sj,suk,"ret",siret);

          // Handle the -+ quadrant
          
	  x2[0]=x1[0]-surv_mp.get_grid_x(sj);
	  x2[1]=x1[1]+surv_mp.get_grid_y(suk);
	  siret=sn_func(2,x2,y1);
          
	  //cout << siret << " " << y1[0] << " " << y1[1] << " ";
          
	  if (siret==0) {
	    surv_mp.set(sj,suk,"e1",y1[0]);
	    surv_mp.set(sj,suk,"e2",y1[1]);
	    surv_mp.set(sj,suk,"s",fabs(y1[0])+fabs(y1[1]));
	  } else {
	    surv_mp.set(sj,suk,"e1",0.0);
	    surv_mp.set(sj,suk,"e2",0.0);
	    surv_mp.set(sj,suk,"s",0.0);
	  }
	  surv_mp.set(sj,suk,"ret",siret);

          // Handle the -- quadrant
	  x2[0]=x1[0]-surv_mm.get_grid_x(sj);
	  x2[1]=x1[1]-surv_mm.get_grid_y(suk);
	  siret=sn_func(2,x2,y1);
	  //cout << siret << " " << y1[0] << " " << y1[1] << endl;
	  if (siret==0) {
	    surv_mm.set(sj,suk,"e1",y1[0]);
	    surv_mm.set(sj,suk,"e2",y1[1]);
	    surv_mm.set(sj,suk,"s",fabs(y1[0])+fabs(y1[1]));
	  } else {
	    surv_mm.set(sj,suk,"e1",0.0);
	    surv_mm.set(sj,suk,"e2",0.0);
	    surv_mm.set(sj,suk,"s",0.0);
	  }
	  surv_mm.set(sj,suk,"ret",siret);
	}
      }
      cout << endl;

      // Compute the maximum value of the solution quality in each
      // quadrant. If the return value is non-zero, Set the solution
      // quality to the maximum value. Also, if any value is larger
      // than 1, then set it equal to 1. 
      
      double pps_max=matrix_max_value<ubmatrix,double>
	(surv_pp.get_slice("s"));
      double pps_min=matrix_min_value<ubmatrix,double>
	(surv_pp.get_slice("s"));
      for(size_t sj=0;sj<surv_N;sj++) {
	for(size_t suk=0;suk<surv_N;suk++) {
	  if (fabs(surv_pp.get(sj,suk,"ret"))>1.0e-4) {
	    surv_pp.set(sj,suk,"s",pps_max);
	  }
	  if (surv_pp.get(sj,suk,"s")>1.0) {
	    surv_pp.set(sj,suk,"s",1.0);
	  }
	}
      }
      double pms_max=matrix_max_value<ubmatrix,double>
	(surv_pm.get_slice("s"));
      double pms_min=matrix_min_value<ubmatrix,double>
	(surv_pm.get_slice("s"));
      for(size_t sj=0;sj<surv_N;sj++) {
	for(size_t suk=0;suk<surv_N;suk++) {
	  if (fabs(surv_pm.get(sj,suk,"ret"))>1.0e-4) {
	    surv_pm.set(sj,suk,"s",pms_max);
	  }
	  if (surv_pm.get(sj,suk,"s")>1.0) {
	    surv_pm.set(sj,suk,"s",1.0);
	  }
	}
      }
      double mps_max=matrix_max_value<ubmatrix,double>
	(surv_mp.get_slice("s"));
      double mps_min=matrix_min_value<ubmatrix,double>
	(surv_mp.get_slice("s"));
      for(size_t sj=0;sj<surv_N;sj++) {
	for(size_t suk=0;suk<surv_N;suk++) {
	  if (fabs(surv_mp.get(sj,suk,"ret"))>1.0e-4) {
	    surv_mp.set(sj,suk,"s",mps_max);
	  }
	  if (surv_mp.get(sj,suk,"s")>1.0) {
	    surv_mp.set(sj,suk,"s",1.0);
	  }
	}
      }
      double mms_max=matrix_max_value<ubmatrix,double>
	(surv_mm.get_slice("s"));
      double mms_min=matrix_min_value<ubmatrix,double>
	(surv_mm.get_slice("s"));
      for(size_t sj=0;sj<surv_N;sj++) {
	for(size_t suk=0;suk<surv_N;suk++) {
	  if (fabs(surv_mm.get(sj,suk,"ret"))>1.0e-4) {
	    surv_mm.set(sj,suk,"s",mms_max);
	  }
	  if (surv_mm.get(sj,suk,"s")>1.0) {
	    surv_mm.set(sj,suk,"s",1.0);
	  }
	}
      }

      double min_min=1.0, min_x=0.0, min_y=0.0;
      
      // Compute the optimal point in each quadrant
      cout.setf(ios::showpos);
      cout.precision(10);
      double val;
      size_t suoi, suoj;
      matrix_min_index<ubmatrix,double>(surv_pp.get_slice("s"),
					suoi,suoj,val);
      cout << "quad  optimal           delta x           delta y       "
           << "    x_opt            y_opt" << endl;
      cout << "++   " << val << " " << surv_pp.get_grid_x(suoi) << " "
	   << surv_pp.get_grid_y(suoj) << " " 
	   << x1[0]+surv_pp.get_grid_x(suoi) << " "
	   << x1[1]+surv_pp.get_grid_y(suoj) << endl;
      if (val<min_min) {
        min_min=val;
        min_x=x1[0]+surv_pp.get_grid_x(suoi);
        min_y=x1[1]+surv_pp.get_grid_y(suoj);
      }
      matrix_min_index<ubmatrix,double>(surv_pm.get_slice("s"),
					suoi,suoj,val);
      cout << "+-   " << val << " " << surv_pm.get_grid_x(suoi) << " "
	   << -surv_pm.get_grid_y(suoj) << " " 
	   << x1[0]+surv_pm.get_grid_x(suoi) << " "
	   << x1[1]-surv_pm.get_grid_y(suoj) << endl;
      if (val<min_min) {
        min_min=val;
        min_x=x1[0]+surv_pm.get_grid_x(suoi);
        min_y=x1[1]-surv_pm.get_grid_y(suoj);
      }
      matrix_min_index<ubmatrix,double>(surv_mp.get_slice("s"),
					suoi,suoj,val);
      cout << "-+   " << val << " " << -surv_mp.get_grid_x(suoi) << " "
	   << surv_mp.get_grid_y(suoj) << " " 
	   << x1[0]-surv_mp.get_grid_x(suoi) << " "
	   << x1[1]+surv_mp.get_grid_y(suoj) << endl;
      if (val<min_min) {
        min_min=val;
        min_x=x1[0]-surv_mp.get_grid_x(suoi);
        min_y=x1[1]+surv_mp.get_grid_y(suoj);
      }
      matrix_min_index<ubmatrix,double>(surv_mm.get_slice("s"),
					suoi,suoj,val);
      cout << "--   " << val << " " << -surv_mm.get_grid_x(suoi) << " "
	   << -surv_mm.get_grid_y(suoj) << " " 
	   << x1[0]-surv_mm.get_grid_x(suoi) << " "
	   << x1[1]-surv_mm.get_grid_y(suoj) << endl;
      if (val<min_min) {
        min_min=val;
        min_x=x1[0]-surv_mm.get_grid_x(suoi);
        min_y=x1[1]-surv_mm.get_grid_y(suoj);
      }
      cout << endl;
      cout.unsetf(ios::showpos);

      cout << "min_min: " << min_min << " " << min_x << " " << min_y << endl;
      
      cout.precision(6);

      hdf_file hf;
      hf.open_or_create("survey.o2");
      hdf_output(hf,(const table3d &)surv_pp,"surv_pp");
      hdf_output(hf,(const table3d &)surv_pm,"surv_pm");
      hdf_output(hf,(const table3d &)surv_mp,"surv_mp");
      hdf_output(hf,(const table3d &)surv_mm,"surv_mm");
      hf.close();
    }
    
    if (mh_ret!=0) {
      if (loc_verbose>1) {
	cout << "Failed in eos_fixed_dist(), "
	     << "x1[0], x1[1], y1[0], y1[1]:\n  " 
	     << x1[0] << " " << x1[1] << " " << y1[0] << " "
	     << y1[1] << endl;
	cout << "  qual_best: " << qual_best << endl;
	cout << "  nB,Ye,T[MeV]: " << nB << " " << Ye << " "
             << T*hc_mev_fm << endl;
      }
      if (loc_verbose==8) {
	char ch;
	cin >> ch;
      } else if (loc_verbose>=9) {
	exit(-1);
      }
      return 5;
    }

    // End of if (mh_ret!=0 && (nB<2.2e-12 || alg_mode==2 ||
    // alg_mode==4)) {
  }

  // Perform a final call to solve_nuclei() to ensure
  // the chemical potentials are updated
  sn_func(2,x1,y1);

  // Store the chemical potentials for homogeneous matter
  vdet["mun_gas"]=mun_gas;
  vdet["mup_gas"]=mup_gas;
  
  if (mpi_size==1 && loc_verbose>1) {
    sn_func(2,x1,y1);
    cout << "Success in eos_fixed_dist(), "
	 << "x1[0], x1[1], y1[0], y1[1]:\n  " 
	 << x1[0] << " " << x1[1] << " " << y1[0] << " "
	 << y1[1] << endl;
    cout << "nB,Ye,T[MeV]: " << nB << " " << Ye << " "
         << T*hc_mev_fm << endl;
  }

  // 8/27: log_xn and log_xp are used below, so it's important that
  // they are set to the correct values, independent of the value of
  // mh_ret.
  if (alg_mode==0 || alg_mode==2 || alg_mode==4) {
    log_xn=x1[0];
    log_xp=x1[1];
  } else {
    if(mh_ret==0) {
      log_xn=x1[0];
      log_xp=x1[1];
    }
  }

  // -------------------------------------------------------------
  // Compute free energy density and entropy density
  
  double xn=pow(10.0,log_xn);
  double xp=pow(10.0,log_xp);

  double kappa=1.0-nB/n0;
  double xi=kappa/(1.0-nB*xn/n0-nB*xp/n0);

  // The total of the number density over all nuclei
  double sum_nuc=0.0;

  // Begin with zero and then add up contributions. Free energy
  // density in fm^{-4} and entropy density in fm^{-3}
  double f=0.0;

  thx.en=0.0;

  for (size_t i=0;i<n_nuclei;i++) {
    
    double en_nuc, fr_nuc;
    
    if (nuclei[i].n>1.0e-300) {
      double lambda=sqrt(2.0*pi/nuclei[i].m/T);
      fr_nuc=-T*(log(vomega[i]/nuclei[i].n/pow(lambda,3.0))+1.0)*
	nuclei[i].n;
      // Note that everything here, including vomega_prime, is in
      // units of powers of femtometers
      en_nuc=nuclei[i].n*(log(vomega[i]/nuclei[i].n/pow(lambda,3.0))+
			  5.0/2.0+vomega_prime[i]*T/vomega[i]);
      if (!std::isfinite(fr_nuc)) {
	if (nuclei[i].n<1.0e-200) {
	  nuclei[i].n=0.0;
	  fr_nuc=0.0;
	  en_nuc=0.0;
	} else {
	  cout << "Nuclear free energy not finite in eos_fixed_dist()."
               << endl;
	  cout << nuclei[i].n << " " << nuclei[i].be << " " << lambda
	       << " " << vomega[i] << endl;
	  exit(-1);
	}
      }
    } else {
      nuclei[i].n=0.0;
      fr_nuc=0.0;
      en_nuc=0.0;
    }
    fr_nuc+=nuclei[i].n*nuclei[i].be+1.433e-05*
      pow(nuclei[i].Z,2.39)/hc_mev_fm*nuclei[i].n;
    sum_nuc+=nuclei[i].n;
    
    f+=fr_nuc+nuclei[i].n*Ec[i];
    thx.en+=en_nuc;
  }

  // Final calculations of free energy density, entropy
  // density, and energy density
  
  f+=xi*(th_gas.ed-T*th_gas.en)-T*sum_nuc*log(kappa);
  thx.en+=xi*th_gas.en+sum_nuc*log(kappa);
  thx.ed=f+T*thx.en;

  // Add the HRG contribution, except for neutrons and protons
  // which are already included
  
  if (inc_hrg) {
    for(size_t j=0;j<res_f.size();j++) {
      if (fabs(res_f[j].m-939.0/hc_mev_fm)>1.0) {
        f+=(res_f[j].ed-T*res_f[j].en);
        thx.ed+=res_f[j].ed;
        thx.en+=res_f[j].en;
      }
    }
    for(size_t j=0;j<res_b.size();j++) {
      f+=(res_b[j].ed-T*res_b[j].en);
      thx.ed+=res_b[j].ed;
      thx.en+=res_b[j].en;
    }
  }
  
  // -------------------------------------------------------------

  // Nucleon chemical potentials and pressure
  // are not currently computed 
  mun_full=neutron.m;
  mup_full=proton.m;
  thx.pr=0.0;

  if (loc_verbose==8) {
    char ch;
    cin >> ch;
  } else if (loc_verbose>=9) {
    exit(-1);
  }
  
  return 0;
}


int eos_nuclei::eos_vary_ZN
(double nB, double Ye, double T, double &log_xn, double &log_xp,
 size_t &nuc_Z1, size_t &nuc_N1, 
 thermo &thx, double &mun_full, double &mup_full, bool no_nuclei) {

  bool debug=false;
  if (function_verbose/10%10>1) debug=true;
  
  if (debug) {
    cout << "vary_ZN() computing nB, Ye, T(1/fm): " << nB << " "
	 << Ye << " " << T << endl;
    cout << "  log_xn, log_xp, Z, N: "
	 << log_xn << " " << log_xp << " " << nuc_Z1 << " "
	 << nuc_N1 << endl;
  }
  
  // We want homogeneous matter if we're above the saturation
  // density or if the 'no_nuclei' flag is true
  if (no_nuclei==true || nB>0.16) {
    int A_min, A_max, NmZ_min, NmZ_max;
    double Zbar, Nbar;
    map<string,double> vdet;
    nuc_matter(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,thx,mun_full,
	       mup_full,A_min,A_max,NmZ_min,NmZ_max,vdet);
    nuc_Z1=0;
    nuc_N1=0;
    return 0;
  }
  
  if (nuc_Z1==0 && nuc_N1==0) {
    // If no nucleus is given as an initial guess, just choose Nickel 56
    nuc_Z1=28;
    nuc_N1=28;
  }

  double log_xn_min=0.0, log_xp_min=0.0, f_min=1.0e100;
  size_t Z_min=0, N_min=0;
  bool min_set=false;

  calc_utf8<> calc;
  std::map<std::string,double> vars;
  calc.compile(nucleon_func.c_str());

  vars["nB"]=nB;
  vars["Ye"]=Ye;
  vars["T"]=T*hc_mev_fm;
  vars["i"]=((double)nuc_Z1);
  int delta_Z=((double)calc.eval(&vars));
  vars["i"]=((double)nuc_N1);
  int delta_N=((double)calc.eval(&vars));

  int iZ0=((int)(nuc_Z1+1.0e-6));
  int iN0=((int)(nuc_N1+1.0e-6));

  double log_xn_guess=log_xn;
  double log_xp_guess=log_xp;
  
  for(int iZ=iZ0-delta_Z;iZ<iZ0+delta_Z+0.01;iZ++) {
    for(int iN=iN0-delta_N;iN<iN0+delta_N+0.01;iN++) {
      
      if (iZ>7 && iN>7) {
	
	// Ensure every nucleus starts with the same initial guess for
	// log_xn and log_xp. I'm not sure if this is the right
	// choice.
	log_xn=log_xn_guess;
	log_xp=log_xp_guess;
	int eret=eos_fixed_ZN(nB,Ye,T,log_xn,log_xp,iZ,iN,thx,
			      mun_full,mup_full);
	double f=thx.ed-T*thx.en;
	if (show_all_nuclei) {
	  if (eret==0) {
	    cout.precision(4);
	    cout.setf(ios::showpos);
	    cout << f << " " << f/nB*hc_mev_fm << " ";
	    cout.unsetf(ios::showpos);
	    cout << iZ << " " << iN << " " << log_xn << " " << log_xp << " ";
	    for(size_t j=0;j<6;j++) {
	      cout << nuclei[j].n*nuclei[j].A/nB << " ";
	    }
	    cout << nuclei[5].n*nuclei[5].A/nB;
	    cout << endl;
	  } else {
	    cout.precision(4);
	    cout << eret << " ";
	    cout << iZ << " " << iN << " " << log_xn << " " << log_xp << " ";
	    for(size_t j=0;j<6;j++) {
	      cout << nuclei[j].n*nuclei[j].A/nB << " ";
	    }
	    cout << nuclei[5].n*nuclei[5].A/nB;
	    cout << endl;
	  }
	}
	cout.precision(6);
	if (std::isfinite(f_min) && eret==0 &&
	    (min_set==false || f<f_min)) {
	  min_set=true;
	  f_min=f;
	  Z_min=iZ;
	  N_min=iN;
	  log_xn_min=log_xn;
	  log_xp_min=log_xp;
	}
      }
    }
  }

  if (min_set==false) {
    if (debug) {
      cout << "All nuclei failed: " << nB << " " << Ye << " " << T << endl;
    }
    return 1;
  }
  
  // After the loop, do a final call of eos_fixed_ZN() to make
  // sure the heavy nucleus, thx, and chemical potentials are set
  // properly
  log_xn=log_xn_min;
  log_xp=log_xp_min;
  nuc_Z1=Z_min;
  nuc_N1=N_min;
  
  int lret=eos_fixed_ZN(nB,Ye,T,log_xn,log_xp,nuc_Z1,nuc_N1,thx,
			mun_full,mup_full);
  
  if (lret!=0) {
    return 2;
  }
  if (!std::isfinite(thx.ed) || !std::isfinite(thx.en)) {
    return 3;
  }

  // The values might have changed, so we update them
  log_xn=log_xn_min;
  log_xp=log_xp_min;

  if (debug) {
    cout << "vary_ZN() success, nB, Ye, T(1/fm), f_min: "
	 << nB << " " << Ye << " " << T << " "
	 << f_min << endl;
  }
  
  return 0;
}

int eos_nuclei::store_point
(size_t i_nB, size_t i_Ye, size_t i_T, double nB, double Ye, double T,
 thermo &th, double log_xn, double log_xp, double Zbar, double Nbar,
 double mun_full, double mup_full, ubvector &X,
 double A_min, double A_max, double NmZ_min, double NmZ_max,
 double loc_flag, map<string,double> &vdet) {

  int loc_verbose=function_verbose/10000%10;

  vector<size_t> ix={i_nB,i_Ye,i_T};
  
  double fr=th.ed-T*th.en;
  if (!std::isfinite(fr)) {
    cout << "Free energy not finite in store_point(), (nB,Ye,T)=("
	 << nB << "," << Ye << "," << T*hc_mev_fm << "). Skipping."
	 << endl;
    return 0;
  }

  if (!std::isfinite(Nbar) || !std::isfinite(Zbar)) {
    cout << "Z or N not finite in store_point(), (nB,Ye,T)=("
	 << nB << "," << Ye << "," << T*hc_mev_fm << "). Skipping."
	 << endl;
    return 0;
  }

  int iflag=((int)(tg_flag.get(ix)*(1.0+1.0e-12)));

  if (iflag==iflag_done) {
    double fr_old=tg_Fint.get(ix)/hc_mev_fm*nB;
    if (fr>=fr_old) {
      cout << "Old point has smaller free energy. Old: " << fr_old
	   << " New: " << fr << endl;
      
      // If a good point is already stored, and the current free
      // energy is smaller than the new free energy, then just return
      // without storing anything.
      
      // AWS: 08/05/20: This return statement was commented out,
      // I'm not sure why, but I'm pretty sure it needs to be here,
      // so I'm uncommenting it for now.
      return 0;
      
    } else {
      cout << "New point has smaller free energy. Old: " << fr_old
	   << " New: " << fr << endl;
    }      
  }
  
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
    double n0=0.16;
    n_fraction=xn*(1.0-nB/n0)/(1.0-nB*xn/n0-nB*xp/n0);
    p_fraction=xp*(1.0-nB/n0)/(1.0-nB*xn/n0-nB*xp/n0);    
  } else {
    n_fraction=1.0-Ye;
    p_fraction=Ye;
  }

  if (X.size()<6) {
    O2SCL_ERR("X vector incorrectly sized in store_point().",
	      o2scl::exc_esanity);
  }
  
  // Test to make sure we're storing a good point
  if (alg_mode!=3) {
    double nb_frac_sum=n_fraction+p_fraction+X[0]+X[1]+
      X[2]+X[3]+X[4]+X[5];
    if (fabs(nb_frac_sum-1.0)>mh_tol_rel) {
      cout << "nB fraction mismatch." << endl;
      cout << "  " << n_fraction << " " << p_fraction << " ";
      vector_out(cout,X,true);
      cout << "  nB_frac_sum,mh_tol_rel: " << nb_frac_sum << " "
	   << mh_tol_rel << endl;
      cout << "  nB,Ye,T[1/fm]: " << nB << " " << Ye << " " << T << endl;
      cout << "  Not storing point." << endl;
      exit(-1);
      return 1;
    }
    
    if (alg_mode!=2 && alg_mode!=4) {
      
      // X[0]*nB/4.0 is the number density of alpha particles, so
      // X[0]*nB/4.0*2.0 is the number of protons contributed from
      // alpha particles, then we just divide by nB.
      
      double Ye_sum=p_fraction+X[0]/4.0*2.0+X[1]/2.0+X[2]/3.0+
	X[3]/3.0*2.0+X[4]/4.0*3.0+X[5]/(Nbar+Zbar)*Zbar;
      if (fabs(Ye_sum-Ye)>mh_tol_rel) {
	cout << "Ye fraction mismatch." << endl;
	cout << "  Ye sum,Ye,mh_tol_rel: " << Ye_sum << " " << Ye << " "
	     << mh_tol_rel << endl;
	cout << "  Not storing point." << endl;
	return 2;
      }
    }
  }
    
  tg_log_xn.set(ix,log_xn);
  tg_log_xp.set(ix,log_xp);
  tg_Z.set(ix,Zbar);
  tg_A.set(ix,Zbar+Nbar);
  tg_flag.set(ix,loc_flag);
  tg_Fint.set(ix,(th.ed-T*th.en)/nB*hc_mev_fm);
  tg_Eint.set(ix,th.ed/nB*hc_mev_fm);
  tg_Sint.set(ix,th.en/nB);

  if (loc_verbose>1) {
    cout << "store_point() output:\n  nB, Ye, T [MeV]: "
	 << nB << " " << Ye << " " << T*hc_mev_fm << endl;
    cout << "  log_xn, log_xp, Z, A: "
	 << log_xn << " " << log_xp << " " << Zbar << " "
	 << Zbar+Nbar << endl;
    cout << "  flag Fint Xn Xp: " << loc_flag << " "
	 << (th.ed-T*th.en)/nB*hc_mev_fm << " "
	 << n_fraction << " " << p_fraction << endl;
  }

  tg_Xn.set(ix,n_fraction);
  tg_Xp.set(ix,p_fraction);

  tg_Xalpha.set(ix,X[0]);
  tg_Xd.set(ix,X[1]);
  tg_Xt.set(ix,X[2]);
  tg_XHe3.set(ix,X[3]);
  tg_XLi4.set(ix,X[4]);
  tg_Xnuclei.set(ix,X[5]);

  if (rmf_fields) {
    tg_sigma.set(ix,vdet["sigma"]*hc_mev_fm);
    tg_omega.set(ix,vdet["omega"]*hc_mev_fm);
    tg_rho.set(ix,vdet["rho"]*hc_mev_fm);
  }
  
  if (loc_verbose>1) {
    cout << "  Xalpha, Xd, Xt: " << X[0] << " " << X[1] << " "
	 << X[2] << endl;
    cout << "  XHe3, XLi4, Xnuclei: " << X[3] << " " << X[4] << " "
	 << X[5] << endl;
  }
  
  if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
    tg_A_min.set(ix,A_min);
    tg_A_max.set(ix,A_max);
    tg_NmZ_min.set(ix,NmZ_min);
    tg_NmZ_max.set(ix,NmZ_max);

    if (loc_verbose>1) {
      cout << "  A_min, A_max, NmZ_min, NmZ_max: ";
      cout.precision(4);
      cout << A_min << " " << A_max << " " << NmZ_min << " "
	   << NmZ_max << endl;
      cout.precision(6);
    }
  }

  tg_Sint.set(ix,th.en/nB);

  if (include_muons || with_leptons) {
    tg_mue.set(ix,vdet["mue"]*hc_mev_fm);
  }
  if (include_muons) {
    tg_Ymu.set(ix,vdet["Ymu"]);
  }
  
  if (include_detail) {
    tg_zn.set(ix,vdet["zn"]);
    tg_zp.set(ix,vdet["zp"]);
    tg_F1.set(ix,vdet["F1"]);
    tg_F2.set(ix,vdet["F2"]);
    tg_F3.set(ix,vdet["F3"]);
    tg_F4.set(ix,vdet["F4"]);
    tg_Un.set(ix,vdet["Un"]);
    tg_Up.set(ix,vdet["Up"]);
    tg_msn.set(ix,vdet["msn"]);
    tg_msp.set(ix,vdet["msp"]);
    tg_g.set(ix,vdet["g"]);
    tg_dgdT.set(ix,vdet["dgdT"]);
  }
  
  // AWS 8/4/2020: This section is commented out because the code does
  // not correctly analytically compute mun, mup, Eint and Pint
  if (false && derivs_computed) {
    
    tg_Pint.set(ix,th.pr*hc_mev_fm);

    if (loc_verbose>1) {
      cout << "  Eint, Pint, Sint: " << th.ed/nB*hc_mev_fm << " "
	   << th.pr*hc_mev_fm << " " << th.en/nB << endl;
    }
    tg_mun.set(ix,mun_full*hc_mev_fm);
    tg_mup.set(ix,mup_full*hc_mev_fm);

    if (loc_verbose>1) {
      cout << "  mun, mup: " << mun_full*hc_mev_fm << " "
	   << mup_full*hc_mev_fm << endl;
    }
    
    if (with_leptons) {
      
      electron.n=nB*Ye;
      electron.mu=electron.m;
      relf.pair_density(electron,T);
      photon.massless_calc(T);

      tg_F.set(ix,(fr+electron.ed+photon.ed-
                   T*(electron.en+photon.en))/nB*hc_mev_fm);
      tg_E.set(ix,(th.ed+electron.ed+photon.ed)/nB*hc_mev_fm);
      tg_P.set(ix,(th.pr+electron.pr+photon.pr)*hc_mev_fm);
      tg_S.set(ix,(th.en+electron.en+photon.en)/nB);
      tg_mue.set(ix,electron.mu*hc_mev_fm);
      
    }
    
  }

  if (loc_verbose>8) {
    exit(-1);
  }

  return 0;
}

int eos_nuclei::select_high_T(std::vector<std::string> &sv,
                              bool itive_com) {
  if (sv.size()<2) {
    cerr << "Command \"select-high-T\" needs an integer argument."
	 << endl;
    return 1;
  }
  return select_high_T_internal(o2scl::stoi(sv[1]));
}

int eos_nuclei::select_high_T_internal(int option) {
  if (option==0) {
    // The original DSH fit
    sk_Tcorr.t0=5.067286719233e+03;
    sk_Tcorr.t1=1.749251370992e+00;
    sk_Tcorr.t2=-4.721193938990e-01;
    sk_Tcorr.t3=-1.945964529505e+05;
    sk_Tcorr.x0=4.197555064408e+01;
    sk_Tcorr.x1=-6.947915483747e-02;
    sk_Tcorr.x2=4.192016722695e-01;
    sk_Tcorr.x3=-2.877974634128e+01;
    sk_Tcorr.alpha=0.144165;
    sk_Tcorr.a=0.0;
    sk_Tcorr.b=1.0;
    sk_Tcorr.W0=0.0;
    sk_Tcorr.b4=0.0;
    sk_Tcorr.b4p=0.0;
    eos_Tcorr=&sk_Tcorr;
  } else if (option==1) {
    // NRAPR from Steiner et al. (2005)
    sk_Tcorr.t0=sk_alt.t0;
    sk_Tcorr.t1=sk_alt.t1;
    sk_Tcorr.t2=sk_alt.t2;
    sk_Tcorr.t3=sk_alt.t3;
    sk_Tcorr.x0=sk_alt.x0;
    sk_Tcorr.x1=sk_alt.x1;
    sk_Tcorr.x2=sk_alt.x2;
    sk_Tcorr.x3=sk_alt.x3;
    sk_Tcorr.alpha=sk_alt.alpha;
    sk_Tcorr.a=sk_alt.a;
    sk_Tcorr.b=sk_alt.b;
    sk_Tcorr.W0=sk_alt.W0;
    sk_Tcorr.b4=sk_alt.b4;
    sk_Tcorr.b4p=sk_alt.b4p;
    eos_Tcorr=&sk_Tcorr;
  } else if (option>=2 && option<=5) {
    
    skyrme_ext.a=0.0;
    skyrme_ext.b=1.0;

    if (option==2) {
      // Skchi414 from Lim and Holt (2017)
      skyrme_ext.t0=-1734.0261/hc_mev_fm;
      skyrme_ext.t1=255.6550/hc_mev_fm;
      skyrme_ext.t2=-264.0678/hc_mev_fm;
      skyrme_ext.t3=12219.5884/hc_mev_fm;
      skyrme_ext.t4=556.1320/hc_mev_fm;
      skyrme_ext.t5=0.0;
      skyrme_ext.x0=0.4679;
      skyrme_ext.x1=-0.5756;
      skyrme_ext.x2=-0.3955;
      skyrme_ext.x3=0.7687;
      skyrme_ext.x4=-15.8761;
      skyrme_ext.x5=0.0;
      skyrme_ext.alpha=1.0/3.0;
      skyrme_ext.alpha2=1.0;
      skyrme_ext.alpha3=1.0;
      skyrme_ext.W0=93.7236/hc_mev_fm;
    } else if (option==3) {
      // Skchi450 from Lim and Holt (2017)
      skyrme_ext.t0=-1803.2928/hc_mev_fm;
      skyrme_ext.t1=301.8208/hc_mev_fm;
      skyrme_ext.t2=-273.2827/hc_mev_fm;
      skyrme_ext.t3=12783.8619/hc_mev_fm;
      skyrme_ext.t4=564.1049/hc_mev_fm;
      skyrme_ext.t5=0.0;
      skyrme_ext.x0=0.4430;
      skyrme_ext.x1=-0.3622;
      skyrme_ext.x2=-0.4105;
      skyrme_ext.x3=0.6545;
      skyrme_ext.x4=-11.3160;
      skyrme_ext.x5=0.0;
      skyrme_ext.alpha=1.0/3.0;
      skyrme_ext.alpha2=1.0;
      skyrme_ext.alpha3=1.0;
      skyrme_ext.W0=106.4288/hc_mev_fm;
    } else if (option==4) {
      // Skchi500 from Lim and Holt (2017)
      skyrme_ext.t0=-1747.48258/hc_mev_fm;
      skyrme_ext.t1=241.31968/hc_mev_fm;
      skyrme_ext.t2=-331.04118/hc_mev_fm;
      skyrme_ext.t3=12491.50533/hc_mev_fm;
      skyrme_ext.t4=405.03174/hc_mev_fm;
      skyrme_ext.t5=0.0;
      skyrme_ext.x0=0.59530;
      skyrme_ext.x1=-1.15893;
      skyrme_ext.x2=-0.58432;
      skyrme_ext.x3=1.20050;
      skyrme_ext.x4=-25.49381;
      skyrme_ext.x5=0.0;
      skyrme_ext.alpha=1.0/3.0;
      skyrme_ext.alpha2=1.0;
      skyrme_ext.alpha3=1.0;
      skyrme_ext.W0=98.08897/hc_mev_fm;
    } else if (option==5) {
      skyrme_ext.t0=-1652.9500/hc_mev_fm;
      skyrme_ext.t1=405.3203/hc_mev_fm;
      skyrme_ext.t2=61.7621/hc_mev_fm;
      skyrme_ext.t3=9769.4593/hc_mev_fm;
      skyrme_ext.t4=188.5656/hc_mev_fm;
      skyrme_ext.t5=758.1091/hc_mev_fm;
      skyrme_ext.x0=0.3956;
      skyrme_ext.x1=-0.7307;
      skyrme_ext.x2=-3.8054;
      skyrme_ext.x3=0.7012;
      skyrme_ext.x4=2.3804;
      skyrme_ext.x5=-4.3589;
      skyrme_ext.alpha=1.0/3.0;
      skyrme_ext.alpha2=2.0/3.0;
      skyrme_ext.alpha3=1.0;
      skyrme_ext.W0=115.5764/hc_mev_fm;
    }
    
    skyrme_ext.parent_method=true;
    skyrme_ext.saturation();
    eos_Tcorr=&skyrme_ext;

    test_mgr tm;
    tm.set_output_level(2);
    cout << skyrme_ext.n0 << endl;
    tm.test_abs(skyrme_ext.n0,0.16,0.02,"n0");
    cout << skyrme_ext.eoa*hc_mev_fm << endl;
    tm.test_abs(skyrme_ext.eoa*hc_mev_fm,-16.0,2.0,"eoa");
    cout << skyrme_ext.comp*hc_mev_fm << endl;
    tm.test_abs(skyrme_ext.comp*hc_mev_fm,240.0,40.0,"comp");
    cout << skyrme_ext.esym*hc_mev_fm << endl;
    tm.test_abs(skyrme_ext.esym*hc_mev_fm,30.0,5.0,"esym");
    cout << skyrme_ext.msom << endl;
    tm.test_abs(skyrme_ext.msom,0.8,0.3,"msom");
    tm.report();
    
  } else if (option==6) {
    
    // Sk chi m*
    sk_Tcorr.t0=-2260.7/hc_mev_fm;
    sk_Tcorr.t1=433.189/hc_mev_fm;
    sk_Tcorr.t2=274.553/hc_mev_fm;
    sk_Tcorr.t3=12984.4/hc_mev_fm;
    sk_Tcorr.x0=0.327488;
    sk_Tcorr.x1=-1.088968;
    sk_Tcorr.x2=-1.822404;
    sk_Tcorr.x3=0.4429;
    sk_Tcorr.alpha=0.198029;
    sk_Tcorr.a=0.0;
    sk_Tcorr.b=1.0;
    sk_Tcorr.W0=0.0;
    sk_Tcorr.b4=0.0;
    sk_Tcorr.b4p=0.0;
    eos_Tcorr=&sk_Tcorr;

    sk_Tcorr.saturation();

    test_mgr tm;
    tm.set_output_level(0);
    tm.test_abs(sk_Tcorr.n0,0.16,0.02,"n0");
    tm.test_abs(sk_Tcorr.eoa*hc_mev_fm,-16.0,2.0,"eoa");
    tm.test_abs(sk_Tcorr.comp*hc_mev_fm,240.0,40.0,"comp");
    tm.test_abs(sk_Tcorr.esym*hc_mev_fm,30.0,5.0,"esym");
    tm.test_abs(sk_Tcorr.msom,0.8,0.2,"msom");
    if (false) {
      cout << "n0: " << sk_Tcorr.n0 << endl;
      cout << "EoA: " << sk_Tcorr.eoa*hc_mev_fm << endl;
      cout << "K: " << sk_Tcorr.comp*hc_mev_fm << endl;
      cout << "Esym: " << sk_Tcorr.esym*hc_mev_fm << endl;
      cout << "msom: " << sk_Tcorr.msom << endl;
    }
    if (tm.report()==false) {
      O2SCL_ERR("sk_Tcorr failed.",o2scl::exc_esanity);
    }
    
  } else if (option==7) {
    
    lim_holt.alphaL=-328.6611/hc_mev_fm;
    lim_holt.alphaU=-938.9511/hc_mev_fm;
    lim_holt.betaL=5.4410/hc_mev_fm;
    lim_holt.betaU=82.0118/hc_mev_fm;
    lim_holt.etaL=366.8049/hc_mev_fm;
    lim_holt.etaU=1033.4339/hc_mev_fm;
    lim_holt.zetaL=88.1020/hc_mev_fm;
    lim_holt.zetaU=-286.3289/hc_mev_fm;
    lim_holt.theta=2.5554/hc_mev_fm;
    lim_holt.thetaL=38.9387/hc_mev_fm;
    lim_holt.gamma=0.3333;
    lim_holt.gamma2=1.0000;
    lim_holt.sigma=1.0000;
    
    eos_Tcorr=&lim_holt;

    lim_holt.def_sat_mroot.verbose=2;
    lim_holt.saturation();

    test_mgr tm;
    tm.set_output_level(2);
    tm.test_abs(lim_holt.n0,0.16,0.02,"n0");
    tm.test_abs(lim_holt.eoa*hc_mev_fm,-16.0,2.0,"eoa");
    tm.test_abs(lim_holt.comp*hc_mev_fm,240.0,40.0,"comp");
    tm.test_abs(lim_holt.esym*hc_mev_fm,30.0,5.0,"esym");
    tm.test_abs(lim_holt.msom,0.8,0.2,"msom");
    if (true) {
      cout << "n0: " << lim_holt.n0 << endl;
      cout << "EoA: " << lim_holt.eoa*hc_mev_fm << endl;
      cout << "K: " << lim_holt.comp*hc_mev_fm << endl;
      cout << "Esym: " << lim_holt.esym*hc_mev_fm << endl;
      cout << "msom: " << lim_holt.msom << endl;
      exit(-1);
    }
    if (tm.report()==false) {
      O2SCL_ERR("sk_Tcorr failed.",o2scl::exc_esanity);
    }
    
  }
  return 0;
}

int eos_nuclei::fit_frdm(std::vector<std::string> &sv,
			 bool itive_com) {
  

  ubvector fit_params(10);
  
  cout << "Fitting nuclear mass formula." << endl;
  nucdist_set(nm_fit.dist,m95);
  nm_fit.def_mmin.ntrial*=10;
  double res;
  nm_fit.eval(frdm,res);
  cout << res << endl;
  nm_fit.fit(frdm,res);
  cout << res << endl;
  nm_fit.fit(frdm,res);
  cout << res << endl;
  cout << "Done fitting nuclear mass formula." << endl;
  
  frdm.guess_fun(10,fit_params);
  cout.precision(12);
  for(size_t i=0;i<10;i++) {
    cout << "  fit_params[" << i << "]=" << fit_params[i]
	 << ";" << endl;
  }
  
  return 0;
}

int eos_nuclei::write_nuclei(std::vector<std::string> &sv,
			     bool itive_com) {
  if (sv.size()<2) {
    cerr << "No filename specified in write_nuclei()." << endl;
    return 1;
  }
  write_nuclei(sv[1]);
  return 0;
}

void eos_nuclei::write_nuclei(std::string fname) {

  cout << "Function write_nuclei() file " << fname << endl;
  
  hdf_file hf;
  
  wordexp_single_file(fname);
  
  hf.open_or_create(fname);

  nucleus nuc_temp;
  
  table<> t;
  t.line_of_names("Z N g m be Sn Sp mass_type spin_type");
  for(size_t i=0;i<5;i++) {

    ame.get_nucleus(nuclei[i].Z,nuclei[i].N-1,nuc_temp);
    double Sn=-(nuclei[i].be-nuc_temp.be)*hc_mev_fm;
    ame.get_nucleus(nuclei[i].Z-1,nuclei[i].N,nuc_temp);
    double Sp=-(nuclei[i].be-nuc_temp.be)*hc_mev_fm;
    
    double line[9]={((double)nuclei[i].Z),
      ((double)nuclei[i].N),
      nuclei[i].g,nuclei[i].m*hc_mev_fm,
      nuclei[i].be*hc_mev_fm,Sn,Sp,1,1};
    t.line_of_data(9,line);
  }
  
  for(int Z=5;Z<=200;Z++) {
    for(int N=5;N<=200;N++) {
      
      if (N<=max_ratio*Z && Z<=max_ratio*N) {
	
	double line[9];
	nucleus nuc;
        bool skip_nucleus=false;
	
	if (ame.is_included(Z,N) &&
	    ame.is_included(Z-1,N) &&
	    ame.is_included(Z,N-1)) {
	  
	  ame.get_nucleus(Z,N,nuc);
	  ame.get_nucleus(Z,N-1,nuc_temp);
	  double Sn=-(nuc.be-nuc_temp.be)*hc_mev_fm;
	  ame.get_nucleus(Z-1,N,nuc_temp);
	  double Sp=-(nuc.be-nuc_temp.be)*hc_mev_fm;
	  
	  line[0]=nuc.Z;
	  line[1]=nuc.N;
	  line[2]=nuc.g;
	  line[3]=nuc.m*hc_mev_fm;
	  line[4]=nuc.be*hc_mev_fm;
	  line[5]=Sn;
	  line[6]=Sp;
	  line[7]=2;
	  
	} else if (m95.is_included(Z,N) &&
		   m95.is_included(Z-1,N) &&
		   m95.is_included(Z,N-1)) {
	  
	  m95.get_nucleus(Z,N,nuc);
	  m95.get_nucleus(Z,N-1,nuc_temp);
	  double Sn=-(nuc.be-nuc_temp.be)*hc_mev_fm;
	  m95.get_nucleus(Z-1,N,nuc_temp);
	  double Sp=-(nuc.be-nuc_temp.be)*hc_mev_fm;
	  
	  line[0]=nuc.Z;
	  line[1]=nuc.N;
	  line[2]=nuc.g;
	  line[3]=nuc.m*hc_mev_fm;
	  line[4]=nuc.be*hc_mev_fm;
	  line[5]=Sn;
	  line[6]=Sp;
	  line[7]=3;
	  
	} else if (extend_frdm) {
	  
	  frdm.get_nucleus(Z,N,nuc);
	  frdm.get_nucleus(Z,N-1,nuc_temp);
	  double Sn=-(nuc.be-nuc_temp.be)*hc_mev_fm;
	  frdm.get_nucleus(Z-1,N,nuc_temp);
	  double Sp=-(nuc.be-nuc_temp.be)*hc_mev_fm;
	  
	  line[0]=nuc.Z;
	  line[1]=nuc.N;
	  line[2]=nuc.g;
	  line[3]=nuc.m*hc_mev_fm;
	  line[4]=nuc.be*hc_mev_fm;
	  line[5]=Sn;
	  line[6]=Sp;
	  line[7]=4;
	  
	} else {
          skip_nucleus=true;
        }

        if (skip_nucleus==false) {
          if (hfb.is_included(Z,N)) {
            if (hfb.get_ZN(Z,N).Jexp<99) {
              nuc.g=2.0*hfb.get_ZN(Z,N).Jexp+1.0;
              line[8]=2;
            } else {
              nuc.g=2.0*hfb.get_ZN(Z,N).Jth+1.0;
              line[8]=3;
            }
          } else {
            line[8]=4;
            if (Z%2==0 && N%2==0) {
              nuc.g=1.0;
              //} else if (Z%2==1 && N%2==1) {
              //nuc.g=3.0;
            } else {
              nuc.g=2.0;
            }
          }
          
          t.line_of_data(9,line);
        }
      }
    }
  }

  hdf_output(hf,t,"mass_table");

  hf.close();
  
  return;
}
  
int eos_nuclei::write_results(std::string fname) {

  int mpi_rank=0, mpi_size=1;
#ifndef NO_MPI
  // Get MPI rank, etc.
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);

  /*
  // Ensure that multiple MPI ranks aren't reading from the
  // filesystem at the same time
  int tag=0, buffer=0;
  if (mpi_size>1 && mpi_rank>=1) {
  MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,
  tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
  */
  
  if (mpi_size>1 && mpi_rank>0) {
    cerr << "Shouldn't output multiple ranks to the same file "
         << "in write_results." << endl;
    return 2;
  }
#endif
  
  cout << "Function write_results(): rank " << mpi_rank
       << " writing file " << fname << endl;
  
  hdf_file hf;
  
  wordexp_single_file(fname);
  
  hf.open_or_create(fname);

  if (model_selected) {
    bool matched=false;
    std::string eos_str;
    if (use_alt_eos==false) {
      eos_str=((string)"0 ");
      eos_str+=o2scl::itos(i_ns)+" ";
      eos_str+=o2scl::itos(i_skyrme)+" ";
      eos_str+=o2scl::dtos(qmc_a,-1,true)+" ";
      eos_str+=o2scl::dtos(qmc_alpha,-1,true)+" ";
      eos_str+=o2scl::dtos(eos_S,-1,true)+" ";
      eos_str+=o2scl::dtos(eos_L,-1,true)+" ";
      eos_str+=o2scl::dtos(phi,-1,true);
      matched=true;
    } else {
      eos_str=((string)"1 ");
      eos_had_skyrme *skp=dynamic_cast<eos_had_skyrme *>(eosp_alt);
      if (skp!=0) {
        matched=true;
        if (alt_name.length()==0) {
          eos_str+="skyrme "+alt_name;
        } else {
          eos_str+="skyrme ";
          eos_str+=o2scl::dtos(sk.t0*hc_mev_fm,-1,true)+" ";
          eos_str+=o2scl::dtos(sk.t1*hc_mev_fm,-1,true)+" ";
          eos_str+=o2scl::dtos(sk.t2*hc_mev_fm,-1,true)+" ";
          eos_str+=o2scl::dtos(sk.t3*hc_mev_fm,-1,true)+" ";
          eos_str+=o2scl::dtos(sk.x0,-1,true)+" ";
          eos_str+=o2scl::dtos(sk.x1,-1,true)+" ";
          eos_str+=o2scl::dtos(sk.x2,-1,true)+" ";
          eos_str+=o2scl::dtos(sk.x3,-1,true)+" ";
          eos_str+=o2scl::dtos(sk.alpha,-1,true);
        }
      }
      eos_had_rmf *rmfp=dynamic_cast<eos_had_rmf *>(eosp_alt);
      if (rmfp!=0) {
        matched=true;
        if (alt_name.length()==0) {
          eos_str+="rmf "+alt_name;
        } else {
          eos_str+="rmf ";
          eos_str+=o2scl::dtos(rmf.cs,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.cw,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.cr,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.b,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.c,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.ms*hc_mev_fm,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.mw*hc_mev_fm,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.mr*hc_mev_fm,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.zeta,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.xi,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.a1,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.a2,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.a3,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.a4,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.a5,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.a6,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.b1,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.b2,-1,true)+" ";
          eos_str+=o2scl::dtos(rmf.b3,-1,true);
        }
      }
    }
    if (matched==false) {
      eos_str="<unknown>";
    }
    cout << "write_results(), model string: " << eos_str << endl;
    hf.sets("model",eos_str);
  }
  
  hf.set_szt("n_nB",n_nB2);
  hf.set_szt("n_Ye",n_Ye2);
  hf.set_szt("n_T",n_T2);
  hf.setd_vec("nB_grid",nB_grid2);
  hf.setd_vec("Ye_grid",Ye_grid2);
  hf.setd_vec("T_grid",T_grid2);

  hdf_output(hf,tg_log_xn,"log_xn");
  hdf_output(hf,tg_log_xp,"log_xp");
  hdf_output(hf,tg_Z,"Z");
  hdf_output(hf,tg_A,"A");
  hdf_output(hf,tg_flag,"flag");
  hdf_output(hf,tg_Fint,"Fint");

  // Note that we write Sint and Eint even if derivs_computed is
  // false, because the entropy derivative is analytical
  hdf_output(hf,tg_Sint,"Sint");
  hdf_output(hf,tg_Eint,"Eint");
  
  hdf_output(hf,tg_Xn,"Xn");
  hdf_output(hf,tg_Xp,"Xp");
  hdf_output(hf,tg_Xalpha,"Xalpha");
  hdf_output(hf,tg_Xnuclei,"Xnuclei");
  hdf_output(hf,tg_Xd,"Xd");
  hdf_output(hf,tg_Xt,"Xt");
  hdf_output(hf,tg_XHe3,"XHe3");
  hdf_output(hf,tg_XLi4,"XLi4");

  if (rmf_fields) {
    hdf_output(hf,tg_sigma,"sigma");
    hdf_output(hf,tg_omega,"omega");
    hdf_output(hf,tg_rho,"rho");
  }
  
  if (alg_mode==2 || alg_mode==3) {
    hdf_output(hf,tg_A_min,"A_min");
    hdf_output(hf,tg_A_max,"A_max");
    hdf_output(hf,tg_NmZ_min,"NmZ_min");
    hdf_output(hf,tg_NmZ_max,"NmZ_max");
  }
  
  if (true && alg_mode==4) {
    hdf_output(hf,tg_A_min,"A_min");
    hdf_output(hf,tg_A_max,"A_max");
    hdf_output(hf,tg_NmZ_min,"NmZ_min");
    hdf_output(hf,tg_NmZ_max,"NmZ_max");
  }
  
  hf.seti("baryons_only",1);
  hf.seti("rmf_fields",rmf_fields);
  hf.seti("alg_mode",alg_mode);
  if (derivs_computed) {
    hf.seti("derivs_computed",1);
  } else {
    hf.seti("derivs_computed",0);
  }
  if (with_leptons) {
    hf.seti("with_leptons",1);
  } else {
    hf.seti("with_leptons",0);
  }
  // The parent class has an include_muons bool but this
  // child class doesn't support muons yet
  hf.seti("include_muons",this->include_muons);

  hf.setd("m_neut",neutron.m*hc_mev_fm);
  hf.setd("m_prot",proton.m*hc_mev_fm);
  hf.setd("hc",hc_mev_fm);
  hf.setd("alpha_em",o2scl_const::fine_structure_f<double>());

  if (with_leptons || include_muons) {
    hdf_output(hf,tg_mue,"mue");
    hdf_output(hf,tg_Ymu,"Ymu");
  }
  
  if (derivs_computed) {
    hdf_output(hf,tg_Pint,"Pint");
    hdf_output(hf,tg_mun,"mun");
    hdf_output(hf,tg_mup,"mup");
    if (with_leptons) {
      hdf_output(hf,tg_F,"F");
      hdf_output(hf,tg_E,"E");
      hdf_output(hf,tg_P,"P");
      hdf_output(hf,tg_S,"S");
    }      
  }

  vector<string> oth_names={"Xd","Xt","XHe3","XLi4","flag",
    "log_xn","log_xp"};
  if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
    oth_names.push_back("A_min");
    oth_names.push_back("A_max");
    oth_names.push_back("NmZ_min");
    oth_names.push_back("NmZ_max");
  }
  vector<string> oth_units={"","","","","","",""};
  if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
    oth_units.push_back("");
    oth_units.push_back("");
    oth_units.push_back("");
    oth_units.push_back("");
  }

  hf.seti("detail",include_detail);
  if (include_detail) {
    hdf_output(hf,tg_zn,"zn");
    hdf_output(hf,tg_zp,"zp");
    hdf_output(hf,tg_msn,"msn");
    hdf_output(hf,tg_msp,"msp");
    hdf_output(hf,tg_F1,"F1");
    hdf_output(hf,tg_F2,"F2");
    hdf_output(hf,tg_F3,"F3");
    hdf_output(hf,tg_F4,"F4");
    hdf_output(hf,tg_Un,"Un");
    hdf_output(hf,tg_Up,"Up");
    hdf_output(hf,tg_g,"g");
    hdf_output(hf,tg_dgdT,"dgdT");
    oth_names.push_back("zn");
    oth_units.push_back("");
    oth_names.push_back("zp");
    oth_units.push_back("");
    oth_names.push_back("msn");
    oth_units.push_back("MeV");
    oth_names.push_back("msp");
    oth_units.push_back("MeV");
    oth_names.push_back("F1");
    oth_units.push_back("MeV");
    oth_names.push_back("F2");
    oth_units.push_back("MeV");
    oth_names.push_back("F3");
    oth_units.push_back("MeV");
    oth_names.push_back("F4");
    oth_units.push_back("MeV");
    oth_names.push_back("Un");
    oth_units.push_back("MeV");
    oth_names.push_back("Up");
    oth_units.push_back("MeV");
    oth_names.push_back("g");
    oth_units.push_back("");
    oth_names.push_back("dgdT");
    oth_units.push_back("1/MeV");
  }

  size_t n_oth=oth_names.size();
  hf.set_szt("n_oth",n_oth);
  hf.sets_vec_copy("oth_names",oth_names);
  hf.sets_vec_copy("oth_units",oth_units);
  
  hf.close();
  
  cout << "Function write_results(): rank " << mpi_rank
       << " done writing file." << endl;

#ifndef NO_MPI
  // Send a message to the next MPI rank
  /*
    if (mpi_size>1 && mpi_rank<mpi_size-1) {
    MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
    tag,MPI_COMM_WORLD);
    }
  */
#endif

  return 0;
}

int eos_nuclei::read_results(std::string fname) {

  int mpi_rank=0, mpi_size=1;
#ifndef NO_MPI
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
  
  cout << "Function read_results(): rank " << mpi_rank
       << " reading file " << fname << endl;
  
  hdf_file hf;
  string type;

  wordexp_single_file(fname);

  hf.open(fname);

  vector<string> name_list;
  hf.list_objects_by_type("tensor_grid",name_list);

  size_t n_oth;
  vector<string> oth_names, oth_units;
  hf.get_szt("n_oth",n_oth);
  hf.gets_vec_copy("oth_names",oth_names);
  hf.gets_vec_copy("oth_units",oth_units);

  //
  if (oth_names.size()!=n_oth) {
    //oth_units.size()!=n_oth) {
    cerr << "Value of n_oth is not equal to size of oth_names or "
         << "oth_units." << endl;
    return 1;
  }
  cout << "n_oth,oth_names.size(),oth_units.size(): "
       << n_oth << " " << oth_names.size() << " "
       << oth_units.size() << endl;

  vector<string> tg_list;
  vector<string> standard_list={"A","E","Eint","F","Fint","P",
    "Pint","S","Sint","Xalpha","Xn","Xp","Ymu","Z","mue",
    "mun","mup"};
  hf.list_objects_by_type("tensor_grid",tg_list,false,1);
  for(size_t i=0;i<n_oth;i++) {
    size_t j;
    if (vector_search(tg_list,oth_names[i],j)==false) {
      cout << "Entry " << oth_names[i] << " in oth_names does not "
           << "correspond to a tensor_grid object." << endl;
    }
  }
  for(size_t i=0;i<tg_list.size();i++) {
    size_t j;
    if (vector_search(oth_names,tg_list[i],j)==false &&
        vector_search(standard_list,tg_list[i],j)==false) {
      cout << "Object " << tg_list[i] << " of type tensor_grid "
           << "not found in oth_names." << endl;
    }
  }
  
  if (true) {
    
    std::string mod_str;
    hf.gets_def("model","",mod_str);
    
    if (mod_str.length()>0) {
      
      vector<string> vs2;
      split_string(mod_str,vs2);
      cout << "read_results(): model string: " << mod_str << endl;
      
      if (false && vs2[0]=="0") {
        
        i_ns=o2scl::stoi(vs2[1]);
        i_skyrme=o2scl::stoi(vs2[2]);
        qmc_a=o2scl::stod(vs2[3]);
        qmc_alpha=o2scl::stod(vs2[4]);
        eos_S=o2scl::stod(vs2[5]);
        eos_L=o2scl::stod(vs2[6]);
        phi=o2scl::stod(vs2[7]);

        cout << "read_results(): selecting model: " << mod_str << endl;
        select_internal(i_ns,i_skyrme,qmc_alpha,qmc_a,
                        eos_L,eos_S,phi);
        
      }
    }
  }
                

  // ----------------------------------------------------------------
  // nB, Ye, T grid (the strangeness grid is taken care of later)

  if (verbose>2) cout << "Reading n_nB." << endl;
  hf.get_szt("n_nB",n_nB2);
  if (verbose>2) cout << "Reading n_Ye." << endl;
  hf.get_szt("n_Ye",n_Ye2);
  if (verbose>2) cout << "Reading n_T." << endl;
  hf.get_szt("n_T",n_T2);
  if (n_nB2==0 || n_Ye2==0 || n_T2==0) {
    O2SCL_ERR("One of the grid counts is zero.",o2scl::exc_efailed);
  }
  
  
  if (verbose>2) cout << "Reading nB_grid." << endl;
  hf.getd_vec("nB_grid",nB_grid2);
  if (verbose>2) cout << "Reading Ye_grid." << endl;
  hf.getd_vec("Ye_grid",Ye_grid2);
  if (verbose>2) cout << "Reading T_grid." << endl;
  hf.getd_vec("T_grid",T_grid2);
  
  // ----------------------------------------------------------------
  // Flags

  int itmp;
  if (verbose>2) cout << "Reading strange_axis." << endl;
  hf.geti_def("strange_axis",0,itmp);
  if (itmp==1) strange_axis=true;
  else strange_axis=false;

  if (strange_axis==true) {
    hf.get_szt_def("n_S",0,n_S2);
    if (n_S2>0) {
      hf.getd_vec("S_grid",S_grid2);
    } else {
      cerr << "Variable strange_axis is 1 but n_S2=0." << endl;
      exit(-1);
    }
  } else {
    n_S2=0;
  }
  
  if (verbose>2) cout << "Reading baryons_only." << endl;
  hf.geti_def("baryons_only",1,itmp);
  if (itmp==1) baryons_only=true;
  else baryons_only=false;
  
  if (verbose>2) cout << "Reading with_leptons." << endl;
  hf.geti_def("with_leptons",0,itmp);
  if (itmp==1) with_leptons=true;
  else with_leptons=false;
  
  if (verbose>2) cout << "Reading derivs_computed." << endl;
  hf.geti_def("derivs_computed",0,itmp);
  if (itmp==1) derivs_computed=true;
  else derivs_computed=false;

  if (with_leptons && !derivs_computed) {
    O2SCL_ERR("File indicates leptons but no derivatives.",
              o2scl::exc_eunimpl);
  }
  
  if (verbose>2) cout << "Reading include_detail." << endl;
  hf.geti_def("detail",0,itmp);
  if (itmp==1) include_detail=true;
  else include_detail=false;
  
  if (verbose>2) cout << "Reading include_muons." << endl;
  hf.geti_def("include_muons",0,itmp);
  if (itmp==1) include_muons=true;
  else include_muons=false;

  if (verbose>2) cout << "Reading alg_mode." << endl;
  hf.geti_def("alg_mode",4,alg_mode);

  if (verbose>2) cout << "Reading rmf_fields." << endl;
  hf.geti_def("rmf_fields",0,itmp);
  if (itmp==1) rmf_fields=true;
  else rmf_fields=false;
  
  // ----------------------------------------------------------------
  // Main data

  if (hf.find_object_by_name("log_xn",type)!=0 || type!="tensor_grid") {
    O2SCL_ERR("Couldn't find tensor log_xn in file.",
	      o2scl::exc_enotfound);
  }
  if (verbose>2) cout << "Reading log_xn." << endl;
  hdf_input(hf,tg_log_xn,"log_xn");
  if (verbose>2) cout << "Reading log_xp." << endl;
  hdf_input(hf,tg_log_xp,"log_xp");
  if (verbose>2) cout << "Reading Z." << endl;
  hdf_input(hf,tg_Z,"Z");
  if (verbose>2) cout << "Reading A." << endl;
  hdf_input(hf,tg_A,"A");
  if (verbose>2) cout << "Reading flag." << endl;
  hdf_input(hf,tg_flag,"flag");
  if (verbose>2) cout << "Reading Fint." << endl;
  hdf_input(hf,tg_Fint,"Fint");
  
  // Note that we read Sint and Eint even if derivs_computed is
  // false, because the entropy derivative is analytical
  
  if (verbose>2) cout << "Reading Eint." << endl;
  hdf_input(hf,tg_Eint,"Eint");
  if (verbose>2) cout << "Reading Sint." << endl;
  hdf_input(hf,tg_Sint,"Sint");
  
  if (verbose>2) cout << "Reading Xn." << endl;
  hdf_input(hf,tg_Xn,"Xn");
  if (verbose>2) cout << "Reading Xp." << endl;
  hdf_input(hf,tg_Xp,"Xp");
  if (verbose>2) cout << "Reading Xalpha." << endl;
  hdf_input(hf,tg_Xalpha,"Xalpha");
  if (verbose>2) cout << "Reading Xnuclei." << endl;
  hdf_input(hf,tg_Xnuclei,"Xnuclei");
  if (verbose>2) cout << "Reading Xd." << endl;
  hdf_input(hf,tg_Xd,"Xd");
  if (verbose>2) cout << "Reading Xt." << endl;
  hdf_input(hf,tg_Xt,"Xt");
  if (verbose>2) cout << "Reading XHe3." << endl;
  hdf_input(hf,tg_XHe3,"XHe3");
  if (verbose>2) cout << "Reading XLi4." << endl;
  hdf_input(hf,tg_XLi4,"XLi4");

  if (rmf_fields) {
    if (verbose>2) cout << "Reading σ." << endl;
    hdf_input(hf,tg_sigma,"sigma");
    if (verbose>2) cout << "Reading ω." << endl;
    hdf_input(hf,tg_omega,"omega");
    if (verbose>2) cout << "Reading ρ." << endl;
    hdf_input(hf,tg_rho,"rho");
  }
  
  // ----------------------------------------------------------------
  // Nuclear distribution
  
  if (alg_mode==2 || alg_mode==3 || (alg_mode==4 && true)) {
    if (hf.find_object_by_name("A_min",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor A_min in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading A_min." << endl;
    hdf_input(hf,tg_A_min,"A_min");
    if (hf.find_object_by_name("A_max",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor A_max in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading A_max." << endl;
    hdf_input(hf,tg_A_max,"A_max");
    if (hf.find_object_by_name("NmZ_min",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor NmZ_min in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading NmZ_min." << endl;
    hdf_input(hf,tg_NmZ_min,"NmZ_min");
    if (hf.find_object_by_name("NmZ_max",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor NmZ_max in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading NmZ_max." << endl;
    hdf_input(hf,tg_NmZ_max,"NmZ_max");
  }

  // ----------------------------------------------------------------
  // Derivatives of free energy and leptons

  if (with_leptons && !derivs_computed) {
    O2SCL_ERR2("File says leptons are present but derivatives are not ",
	       "computed.",o2scl::exc_esanity);
  }

  if (with_leptons || include_muons) {
    if (verbose>2) cout << "Reading mue." << endl;
    hdf_input(hf,tg_mue,"mue");
  } else {
    tg_mue.clear();
  }
  
  if (include_muons) {
    if (verbose>2) cout << "Reading Ymu." << endl;
    hdf_input(hf,tg_Ymu,"Ymu");
  } else {
    tg_Ymu.clear();
  }
  
  if (derivs_computed) {
    if (hf.find_object_by_name("Pint",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor Pint in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading Pint." << endl;
    hdf_input(hf,tg_Pint,"Pint");
    if (hf.find_object_by_name("mun",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor mun in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading mun." << endl;
    hdf_input(hf,tg_mun,"mun");
    if (hf.find_object_by_name("mup",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor mup in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading mup." << endl;
    hdf_input(hf,tg_mup,"mup");
    
    if (with_leptons) {
      if (verbose>2) cout << "Reading F." << endl;
      hdf_input(hf,tg_F,"F");
      if (verbose>2) cout << "Reading E." << endl;
      hdf_input(hf,tg_E,"E");
      if (verbose>2) cout << "Reading P." << endl;
      hdf_input(hf,tg_P,"P");
      if (verbose>2) cout << "Reading S." << endl;
      hdf_input(hf,tg_S,"S");
    } else {
      tg_F.clear();
      tg_E.clear();
      tg_P.clear();
      tg_S.clear();
    }
  } else {
    tg_Pint.clear();
    tg_mun.clear();
    tg_mup.clear();
    tg_F.clear();
    tg_E.clear();
    tg_P.clear();
    tg_S.clear();
  }

  // ----------------------------------------------------------------
  // Detail
  
  if (include_detail) {
    if (verbose>2) cout << "Reading zn." << endl;
    hdf_input(hf,tg_zn,"zn");
    if (verbose>2) cout << "Reading zp." << endl;
    hdf_input(hf,tg_zp,"zp");
    if (verbose>2) cout << "Reading msn." << endl;
    hdf_input(hf,tg_msn,"msn");
    if (verbose>2) cout << "Reading msp." << endl;
    hdf_input(hf,tg_msp,"msp");
    if (verbose>2) cout << "Reading F1." << endl;
    hdf_input(hf,tg_F1,"F1");
    if (verbose>2) cout << "Reading F2." << endl;
    hdf_input(hf,tg_F2,"F2");
    if (verbose>2) cout << "Reading F3." << endl;
    hdf_input(hf,tg_F3,"F3");
    if (verbose>2) cout << "Reading F4." << endl;
    hdf_input(hf,tg_F4,"F4");
    if (verbose>2) cout << "Reading Un." << endl;
    hdf_input(hf,tg_Un,"Un");
    if (verbose>2) cout << "Reading Up." << endl;
    hdf_input(hf,tg_Up,"Up");
    if (verbose>2) cout << "Reading g." << endl;
    hdf_input(hf,tg_g,"g");
    if (verbose>2) cout << "Reading dgdT." << endl;
    hdf_input(hf,tg_dgdT,"dgdT");
  }
  
  hf.close();

  loaded=true;
  
  cout << "Function read_results(): rank " << mpi_rank
       << " done reading file." << endl;

#ifndef NO_MPI
  // Send a message to the next MPI rank
  if (mpi_size>1 && mpi_rank<mpi_size-1) {
    MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
	     tag,MPI_COMM_WORLD);
  }
#endif

  return 0;
}

int eos_nuclei::test_random(std::vector<std::string> &sv,
                            bool itive_com) {

  size_t ntests=o2scl::stoszt(sv[1]);

  o2scl::rng<> rg;
  rg.clock_seed();

  int n_changed=0;
  bool lg=false;
  if (sv.size()>=3 && sv[2]=="lg") {
    std::cout << "Liquid-gas only." << std::endl;
    lg=true;
  }

  vector<double> nB_kist, Ye_kist, T_kist;

  size_t ilo=0, ihi=n_nB2;
  if (lg) {
    ilo=vector_lookup(n_nB2,nB_grid2,0.01);
    ihi=vector_lookup(n_nB2,nB_grid2,0.15);
    // The random int below doesn't include ihi, so we increase
    // this by one to make sure we include ihi
    if (ihi<n_nB2) ihi++;
    cout << "inB_lo, inB_hi: " << ilo << " " << ihi << endl;
  }
  
  for(size_t it=0;it<ntests;it++) {

    size_t inB;
    if (lg) {
      inB=rg.random_int(ihi-ilo)+ilo;
    } else {
      inB=rg.random_int(n_nB2);
    }
    size_t iYe=rg.random_int(n_Ye2);
    size_t iT=rg.random_int(n_T2);
    vector<size_t> ix={inB,iYe,iT};

    double nB=nB_grid2[inB];
    double Ye=Ye_grid2[iYe];
    double T=T_grid2[iT]/hc_mev_fm;

    double log_xn=tg_log_xn.get(ix);
    double log_xp=tg_log_xp.get(ix);
    double A_old=tg_A.get(ix);
    
    int A_min=5;
    int A_max=fd_A_max;
    int NmZ_min=-200;
    int NmZ_max=200;

    thermo thx;
    double mun_full, mup_full;
    
    derivs_computed=false;
    with_leptons=false;
    
    double Zbar, Nbar;
    map<string,double> vdet;
    double log_xn_old=log_xn;
    double log_xp_old=log_xp;
    cout << it+1 << "/" << ntests << ": nB,Ye,T[MeV]: "
         << nB << " " << Ye << " " << T*hc_mev_fm << " log_xn,log_xp: "
         << log_xn << " " << log_xp << endl;
    int ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
                          thx,mun_full,mup_full,
                          A_min,A_max,NmZ_min,NmZ_max,vdet,
                          true,false);
    if (fabs(log_xn-log_xn_old)>1.0e-5 ||
        fabs(log_xp-log_xp_old)>1.0e-5) {
      cout << inB << " " << iYe << " " << iT << endl;
      cout << "  ret,log_xn_old,log_xn,log_xp_old,log_xp: "
           << ret << " " << log_xn_old << " " << log_xn << " "
           << log_xp_old << " " << log_xp << endl;
      cout << "  A_old,A: " << A_old << " " << Zbar+Nbar << endl;
      if (ret==0) {
        ubvector X;
        compute_X(nB,X);
        store_point(inB,iYe,iT,nB,Ye,T,thx,log_xn,log_xp,Zbar,Nbar,
                    mun_full,mup_full,X,A_min,A_max,NmZ_min,NmZ_max,
                    10.0,vdet);
        n_changed++;
        if (nB_kist.size()<10) {
          nB_kist.push_back(nB);
          Ye_kist.push_back(Ye);
          T_kist.push_back(T);
        }
      }
        
      //char ch;
      //cin >> ch;
    }
    
  }

  cout << n_changed << " points changed." << endl;
  if (nB_kist.size()>0) {
    cout << "First " << nB_kist.size() << ":" << endl;
    for(size_t i=0;i<nB_kist.size();i++) {
      cout << nB_kist[i] << " " << Ye_kist[i] << " "
           << T_kist[i]*hc_mev_fm << endl;
    }
  }
  
  return 0;
}

int eos_nuclei::point_nuclei_mu(std::vector<std::string> &sv,
                                bool itive_com) {
  
  double mun=o2scl::function_to_double(sv[1])/hc_mev_fm-neutron.m;
  double mup=o2scl::function_to_double(sv[2])/hc_mev_fm-proton.m;
  double T=o2scl::function_to_double(sv[3])/hc_mev_fm;
  
  double nB=o2scl::function_to_double(sv[4]);
  double Ye=o2scl::function_to_double(sv[5]);

  std::cout << "Command 'point-nuclei-mu' computing EOS at (w/o rest mass)\n"
            << "  mun(1/fm),mup(1/fm),T(1/fm): " << mun << " "
            << mup << " " << T << endl;
  std::cout << "Using guess for nB,Ye: " << nB << " " << Ye << std::endl;
  
  double log_xn, log_xp;
  if (sv.size()>=8) {
    log_xn=o2scl::function_to_double(sv[6]);
    log_xp=o2scl::function_to_double(sv[7]);
  } else {
    if (nB<0.16) {
      log_xn=log10(nB*(1.0-Ye)/2.0);
      log_xp=log10(nB*Ye/2.0);
    } else {
      log_xn=log10(nB*(1.0-Ye));
      log_xp=log10(nB*Ye);
    }
  }
  std::cout << "Using guess for log_xn,log_xp: " << log_xn << " "
            << log_xp << std::endl;
  
  thermo thx, th_gas;
  double mun_full, mup_full;
  double Zbar, Nbar;
  map<string,double> vdet;
  int A_min, A_max, NmZ_min, NmZ_max;

  A_min=5;
  A_max=fd_A_max;
  NmZ_min=-200;
  NmZ_max=200;

  cout << "Solving for matter at initial point nB,Ye,T(1/fm): "
       << nB << " " << Ye << " " << T << endl;
  int ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
                        thx,mun_full,mup_full,
                        A_min,A_max,NmZ_min,NmZ_max,vdet,
                        true,false);
  if (ret!=0) {
    cerr << "Initial point failed in 'point-nuclei-mu'." << endl;
    return 1;
  }

  ubvector x(10), y(10);
  x[0]=nB;
  x[1]=Ye;
  x[2]=log_xn;
  x[3]=log_xp;
  x[4]=log_xn;
  x[5]=log_xp;
  x[6]=log_xn;
  x[7]=log_xp;
  x[8]=log_xn;
  x[9]=log_xp;
  
  mm_funct func=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector&,
                     double,double,double,double &,double &,
                     thermo &)>
     (&eos_nuclei::solve_nuclei_mu),this,std::placeholders::_1,
     std::placeholders::_2,std::placeholders::_3,
     mun,mup,T,std::ref(mun_full),std::ref(mup_full),
     std::ref(th_gas));

  if (func(10,x,y)!=0) {
    cerr << "Initial guess to solver failed in 'point-nuclei-mu'."
         << endl;
    return 2;
  }
          
  mroot_hybrids<> mh2;
  mh2.verbose=2;
  mh2.msolve(10,x,func);

  cout << "Solving for matter at final point nB,Ye,T(1/fm): "
       << nB << " " << Ye << " " << T << endl;
  ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
                        thx,mun_full,mup_full,
                        A_min,A_max,NmZ_min,NmZ_max,vdet,
                        true,false);
  if (ret!=0) {
    cerr << "Final point failed in 'point-nuclei-mu'." << endl;
    return 3;
  }
  
  return 0;
}

int eos_nuclei::point_nuclei(std::vector<std::string> &sv,
			     bool itive_com) {

  if (sv.size()<4) {
    cerr << "Not enough arguments in point-nuclei." << endl;
    return 2;
  }

  double nB=o2scl::function_to_double(sv[1]);
  double Ye=o2scl::function_to_double(sv[2]);
  double T=o2scl::function_to_double(sv[3])/hc_mev_fm;

  double log_xn, log_xp;
  size_t nuc_Z1, nuc_N1;
  int A_min, A_max, NmZ_min, NmZ_max;
  bool guess_provided=false;
  int flag=0;
  size_t inB=0, iYe=0, iT=0;
  vector<size_t> ix, ix_left, ix_right;

  // -------------------------------------------------------------------
  // If a table is loaded, then fix the density, electron fraction,
  // and temperature to a grid point
  
  if (loaded) {

    inB=vector_lookup(n_nB2,nB_grid2,nB);
    nB=nB_grid2[inB];
    iYe=vector_lookup(n_Ye2,Ye_grid2,Ye);
    Ye=Ye_grid2[iYe];
    iT=vector_lookup(n_T2,T_grid2,T*hc_mev_fm);
    T=T_grid2[iT]/hc_mev_fm;
    ix={inB,iYe,iT};
    if (inB>0) {
      ix_left={inB-1,iYe,iT};
    }
    if (inB<n_nB2-1) {
      ix_right={inB+1,iYe,iT};
    }
    
    flag=((int)(tg_flag.get(ix)+1.0e-10));
    
    cout << "Found grid point (inB,iYe,iT)=(" << inB << ","
	 << iYe << "," << iT << ")\n\t(nB,Ye,T[MeV])=("
	 << nB << "," << Ye << "," << T*hc_mev_fm << "), flag="
	 << flag << "." << endl;
  } else {
    cout << "No table loaded. Using exact (nB,Ye,T) specified." << endl;
  }

  // -------------------------------------------------------------------
  // Attempt to find a good initial guess, either from the command
  // line or from the table
  
  if ((alg_mode==2 || alg_mode==3 || alg_mode==4) && sv.size()>=10) {
    
    cout << "Function point_nuclei(): "
	 << "Reading guess (log_xn,log_xp,A_min,A_max,NmZ_min,NmZ_max) "
	 << "from command line." << endl;
    
    log_xn=o2scl::function_to_double(sv[4]);
    log_xp=o2scl::function_to_double(sv[5]);
    A_min=((size_t)(o2scl::function_to_double(sv[6])*(1.0+1.0e-12)));
    A_max=((size_t)(o2scl::function_to_double(sv[7])*(1.0+1.0e-12)));
    NmZ_min=((size_t)(o2scl::function_to_double(sv[8])*(1.0+1.0e-12)));
    NmZ_max=((size_t)(o2scl::function_to_double(sv[9])*(1.0+1.0e-12)));
    cout << "  " << log_xn << " " << log_xp << " " << A_min << " "
	 << A_max << " " << NmZ_min << " " << NmZ_max << endl;
    guess_provided=true;
    
  } else if ((alg_mode==0 || alg_mode==1) && sv.size()>=8) {
    
    cout << "Function point_nuclei(): "
	 << "Reading guess (log_xn,log_xp,N,Z) from command line." << endl;
    
    log_xn=o2scl::function_to_double(sv[4]);
    log_xp=o2scl::function_to_double(sv[5]);
    nuc_Z1=((size_t)(o2scl::function_to_double(sv[6])*(1.0+1.0e-12)));
    nuc_N1=((size_t)(o2scl::function_to_double(sv[7])*(1.0+1.0e-12)));
    cout << "  " << log_xn << " "
	 << log_xp << " " << nuc_Z1 << " " << nuc_N1 << endl;
    guess_provided=true;
    
  } else if (loaded) {
    
    // If the flag is 'guess' or 'in progress with guess'
    if (flag>=5 || flag==-5) {
      
      if (alg_mode==0 || alg_mode==1) {
	cout << "Function point_nuclei(): "
	     << "Reading guess (N,Z) from current table." << endl;
	log_xn=tg_log_xn.get(ix);
	log_xp=tg_log_xp.get(ix);
	nuc_Z1=((size_t)(tg_Z.get(ix)));
	nuc_N1=((size_t)(tg_A.get(ix)))-nuc_Z1;
	cout << "  " << log_xn << " "
	     << log_xp << " " << nuc_Z1 << " " << nuc_N1 << endl;
	guess_provided=true;
      } else {
	cout << "Function point_nuclei():\n"
	     << "  Reading guess (log_xn,log_xp,A_min,A_max,NmZ_min,NmZ_max) "
	     << "from current table." << endl;
	log_xn=tg_log_xn.get(ix);
	log_xp=tg_log_xp.get(ix);
	A_min=((int)(tg_A_min.get(ix)));
	A_max=((int)(tg_A_max.get(ix)));
	NmZ_min=((int)(tg_NmZ_min.get(ix)));
	NmZ_max=((int)(tg_NmZ_max.get(ix)));
	cout << "  " << log_xn << " " << log_xp << " " << A_min << " "
	     << A_max << " " << NmZ_min << " " << NmZ_max << endl;
	guess_provided=true;
      }
      
    } else if (inB>0 && tg_flag.get(ix_left)>9.9 && alg_mode==4) {

      // Otherwise if there is a good guess from the adjacent point
      // at lower baryon density
      
      cout << "Reading guess from lower baryon density point." << endl;
      log_xn=tg_log_xn.get(ix_left);
      log_xp=tg_log_xp.get(ix_left);
      A_min=((int)(tg_A_min.get(ix_left)));
      A_max=((int)(tg_A_max.get(ix_left)));
      NmZ_min=((int)(tg_NmZ_min.get(ix_left)));
      NmZ_max=((int)(tg_NmZ_max.get(ix_left)));
      cout << "  log_xn,log_xp,Z,A: "
	   << log_xn << " " << log_xp << " ";
      cout << tg_Z.get(ix_left) << " ";
      cout << tg_A.get(ix_left) << endl;
      cout << "  A_min,A_Max,NmZ_min,NmZ_max: ";
      cout << A_min << " "
	   << A_max << " " << NmZ_min << " " << NmZ_max << endl;
      guess_provided=true;
    }
  }

  // -------------------------------------------------------------------
  // If no guess was found, then just use a default guess
  
  if (guess_provided==false) {
    cout << "Function point_nuclei(): "
	 << "Using hard-coded initial guess." << endl;
    log_xn=-1.0;
    log_xp=-1.0;
    if (alg_mode==2 || alg_mode==3) {
      A_min=20;
      A_max=40;
      NmZ_min=-5;
      NmZ_max=5;
    } else if (alg_mode==4) {
      A_min=5;
      A_max=fd_A_max;
      NmZ_min=-200;
      NmZ_max=200;
    } else {
      nuc_Z1=28;
      nuc_N1=28;
    }
  }

  // -------------------------------------------------------------------
  // If necessary and possible, compute the point
  
  int ret=-1;
  
  if (flag<10 || recompute==true) {

    if (flag<10) {
      cout << "Computing point because flag < 10." << endl;
    } else {
      cout << "Computing point because recompute is true." << endl;
    }

    if (loaded && inB>0 && inB<n_nB2-1 &&
	tg_flag.get(ix_left)>9.9 &&
	tg_flag.get(ix_right)>9.9) {
      cout << "Automatically setting ranges for log_xn and log_xp "
           << "from neighboring points." << endl;
      fd_rand_ranges.resize(4);
      fd_rand_ranges[0]=tg_log_xn.get(ix_left);
      fd_rand_ranges[1]=tg_log_xn.get(ix_right);
      fd_rand_ranges[2]=tg_log_xp.get(ix_left);
      fd_rand_ranges[3]=tg_log_xp.get(ix_right);
      cout << "  ";
      vector_out(cout,fd_rand_ranges,true);
    }
    
    thermo thx;
    double mun_full, mup_full;
    
    double Zbar, Nbar;
    map<string,double> vdet;
    if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
      cout << "alg mode: " << alg_mode << endl;
      ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
			thx,mun_full,mup_full,
			A_min,A_max,NmZ_min,NmZ_max,vdet,
			true,true); // setting this true for no_nuclei to use with pion code
    } else {
      ret=eos_vary_ZN(nB,Ye,T,log_xn,log_xp,nuc_Z1,nuc_N1,
		      thx,mun_full,mup_full,false);
      Nbar=nuc_N1;
      Zbar=nuc_Z1;
    }

    ubvector X;
    if (ret==0) {
      compute_X(nB,X);
    }

    if (ret!=0) {
      cout << "Point failed." << endl;
    } else if (loaded) {
      cout << "Point succeeded. Storing." << endl;

      store_point(inB,iYe,iT,nB,Ye,T,thx,log_xn,log_xp,Zbar,Nbar,
		  mun_full,mup_full,X,A_min,A_max,NmZ_min,NmZ_max,
		  10.0,vdet);
    } else {
      cout << "Point succeeded." << endl;
    }

    if (ret==0) {
      cout << "log_xn: " << log_xn << endl;
      cout << "log_xp: " << log_xp << endl;
      cout << "fint: " << thx.ed-T*thx.en << " 1/fm^4" << endl;
      cout << "Fint: " << (thx.ed-T*thx.en)/nB*hc_mev_fm << " MeV" << endl;
      cout << "A: " << Zbar+Nbar << endl;
      cout << "Z: " << Zbar << endl;
      
      if (true) {
        
        double n0=0.16;
        double κ=1.0-nB/n0;
        double n_fraction, p_fraction, xn, xp, ξ;
        
        if (nB<0.16) {
          
          xn=0.0;
          if (log_xn>-300.0) {
            xn=pow(10.0,log_xn);
          }
          
          xp=0.0;
          if (log_xp>-300.0) {
            xp=pow(10.0,log_xp);
          }
          
          ξ=κ/(1.0-nB*xn/n0-nB*xp/n0);
          n_fraction=xn*ξ;
          p_fraction=xp*ξ;
          
        } else {
          
          xn=1.0-Ye;
          xp=Ye;
          ξ=1.0;
          n_fraction=1.0-Ye;
          p_fraction=Ye;
          
        }

        cout << "ξ: " << ξ << endl;
        cout << "n_{n,gas} [in 1/fm^3; called n_n^{\\prime} in "
             << "Du et al. (2022)]: " << xn*nB << endl;
        cout << "n_{p,gas} [in 1/fm^3; called n_n^{\\prime} in "
             << "Du et al. (2022)]: " << xp*nB << endl;
        cout << "n_{n,avg} [in 1/fm^3: equal to Xn*nB]: "
             << xn*nB*ξ << endl;
        cout << "n_{n,avg} [in 1/fm^3: equal to Xn*nB]: "
             << xp*nB*ξ << endl;
        
        if (inc_hrg) {
          
          table_units<> t;
          store_hrg(vdet["mun_gas"]+neutron.m,
                    vdet["mup_gas"]+proton.m,xn*nB,xp*nB,T,t);
          hdf_file hf;
          hf.open_or_create("hrg.o2");
          hdf_output(hf,t,"hrg");
          hf.close();
          
        }

      }
      
      
      if (include_detail) {
	cout << "zn: " << vdet["zn"] << endl;
	cout << "zp: " << vdet["zp"] << endl;
	cout << "F1: " << vdet["F1"] << " "
             << vdet_units.find("F1")->second << endl;
	cout << "F2: " << vdet["F2"] << " "
             << vdet_units.find("F2")->second << endl;
	cout << "F3: " << vdet["F3"] << " "
             << vdet_units.find("F3")->second << endl;
	cout << "F4: " << vdet["F4"] << " "
             << vdet_units.find("F4")->second << endl;
	cout << "msn: " << vdet["msn"] << " "
             << vdet_units.find("msn")->second << endl;
	cout << "msp: " << vdet["msp"] << " "
             << vdet_units.find("msp")->second << endl;
	cout << "Un: " << vdet["Un"] << " "
             << vdet_units.find("Un")->second << endl;
	cout << "Up: " << vdet["Up"] << " "
             << vdet_units.find("Up")->second << endl;
	cout << "g: " << vdet["g"] << endl;
	cout << "dgdT: " << vdet["dgdT"] << " "
             << vdet_units.find("dgdT")->second << endl;
      }
    }      

    if (fd_rand_ranges.size()>0) fd_rand_ranges.clear();
    
  }

  // (At this point in the code, the integer 'flag' is the old
  // flag, not the new flag, which is stored in tg_flag)

  if (flag>=10 && recompute==false) {
    cout << "Point already computed and recompute is false." << endl;
  }

  // -------------------------------------------------------------------
  // Print out the results when tables have been loaded and if either
  // the computation worked (ret==0) or if the old flag was marked
  // done (flag==10).
  
  if ((ret==0 || flag==10) && loaded) {
    cout << "Results stored in table:" << endl;
    cout << "nB,Ye,T[MeV]: " << nB << " " << Ye << " " << T*hc_mev_fm << endl;
    cout << "flag: " << tg_flag.get(ix) << endl;
    cout << "log_xn: " << tg_log_xn.get(ix) << endl;
    cout << "log_xp: " << tg_log_xp.get(ix) << endl;
    cout << "Z: " << tg_Z.get(ix) << endl;
    cout << "A: " << tg_A.get(ix) << endl;
    cout << "Fint: " << tg_Fint.get(ix) << " MeV" << endl;
    cout << "Sint: " << tg_Sint.get(ix) << endl;
    cout << "Eint: " << tg_Eint.get(ix) << " MeV" << endl;
    cout << "A_min: " << tg_A_min.get(ix) << endl;
    cout << "A_max: " << tg_A_max.get(ix) << endl;
    cout << "NmZ_min: " << tg_NmZ_min.get(ix) << endl;
    cout << "NmZ_max: " << tg_NmZ_max.get(ix) << endl;
    cout << "Xn: " << tg_Xn.get(ix) << endl;
    cout << "Xp: " << tg_Xp.get(ix) << endl;
    cout << "Xd: " << tg_Xd.get(ix) << endl;
    cout << "Xt: " << tg_Xt.get(ix) << endl;
    cout << "XHe3: " << tg_XHe3.get(ix) << endl;
    cout << "XLi4: " << tg_XLi4.get(ix) << endl;
    cout << "Xalpha: " << tg_Xalpha.get(ix) << endl;
    cout << "Xnuclei: " << tg_Xnuclei.get(ix) << endl;
    if (rmf_fields) {
      cout << "sigma: " << tg_sigma.get(ix) << endl;
      cout << "omega: " << tg_omega.get(ix) << endl;
      cout << "rho: " << tg_rho.get(ix) << endl;
    }
    if (include_muons) {
      cout << "Ymu: " << tg_Ymu.get(ix) << endl;
    }
    if (include_muons || with_leptons) {
      cout << "mue: " << tg_mue.get(ix) << " MeV" << endl;
    }
    if (derivs_computed) {
      cout << "Pint: " << tg_Pint.get(ix) << endl;
      cout << "mun: " << tg_mun.get(ix) << endl;
      cout << "mup: " << tg_mup.get(ix) << endl;
      if (with_leptons) {
        cout << "F: " << tg_F.get(ix) << endl;
        cout << "E: " << tg_E.get(ix) << endl;
        cout << "S: " << tg_S.get(ix) << endl;
        cout << "P: " << tg_P.get(ix) << endl;
      }
    }
    if (include_detail) {
      cout << "zn: " << tg_zn.get(ix) << endl;
      cout << "zp: " << tg_zp.get(ix) << endl;
      cout << "F1: " << tg_F1.get(ix) << " MeV" << endl;
      cout << "F2: " << tg_F2.get(ix) << " MeV" << endl;
      cout << "F3: " << tg_F3.get(ix) << " MeV" << endl;
      cout << "F3: " << tg_F3.get(ix) << " MeV" << endl;
      cout << "msn: " << tg_msn.get(ix) << " MeV" << endl;
      cout << "msp: " << tg_msp.get(ix) << " MeV" << endl;
      cout << "Un: " << tg_Un.get(ix) << " MeV" << endl;
      cout << "Up: " << tg_Up.get(ix) << " MeV" << endl;
      cout << "g: " << tg_g.get(ix) << endl;
      cout << "dgdT: " << tg_dgdT.get(ix) << " 1/MeV" << endl;
    }
  }

  // -------------------------------------------------------------------
  // Output the nuclear distribution
  
  if (alg_mode>=2 && show_all_nuclei) {
    
    cout << "Writing distribution to dist.o2." << endl;
    
    table3d t3d;
    A_max=((int)(tg_A_max.get(ix)+1.0e-10));
    t3d.set_xy("N",uniform_grid_end_width<double>(0.0,A_max,1.0),
	       "Z",uniform_grid_end_width<double>(0.0,A_max,1.0));
    
    // New slices
    t3d.new_slice("n_nuc");
    t3d.set_slice_all("n_nuc",0.0);
    t3d.new_slice("X_nuc");
    t3d.set_slice_all("X_nuc",0.0);
    
    // Initialize these to zero to avoid warnings, but they
    // are always set in the i=0 case below
    double n_nuc_min=0.0, X_nuc_min=0.0;
    for(size_t i=0;i<nuclei.size();i++) {
      // Compute minimums
      if (i==0) {
	n_nuc_min=nuclei[i].n;
	X_nuc_min=nuclei[i].n*nuclei[i].A/nB;
      } else {
	if (nuclei[i].n<n_nuc_min) {
	  n_nuc_min=nuclei[i].n;
	}
	if (nuclei[i].n*nuclei[i].A/nB<X_nuc_min) {
	  X_nuc_min=nuclei[i].n*nuclei[i].A/nB;
	}
      }
      // Set n_nuc and X_nuc
      t3d.set(nuclei[i].N,nuclei[i].Z,"n_nuc",nuclei[i].n);
      t3d.set(nuclei[i].N,nuclei[i].Z,"X_nuc",nuclei[i].n*nuclei[i].A/nB);
    }
    
    // Now take logs, replacing zeros with the minimum value
    t3d.new_slice("log10_n_nuc");
    t3d.set_slice_all("log10_n_nuc",log10(n_nuc_min));
    t3d.new_slice("log10_X_nuc");
    t3d.set_slice_all("log10_X_nuc",log10(X_nuc_min));
    for(size_t i=0;i<nuclei.size();i++) {
      t3d.set(nuclei[i].N,nuclei[i].Z,"log10_n_nuc",log10(nuclei[i].n));
      t3d.set(nuclei[i].N,nuclei[i].Z,"log10_X_nuc",log10(nuclei[i].n/nB));
    }
    
    // Now create a table adding up isotopes
    table<> tab;
    tab.line_of_names("Z X_nuc");
    for(int i=1;i<=A_max;i++) {
      double line[2]={((double)(i+1)),0.0};
      tab.line_of_data(2,line);
    }
    for(size_t i=0;i<nuclei.size();i++) {
      int iZ=nuclei[i].Z;
      if (iZ>=1 && iZ<=A_max) {
	tab.set("X_nuc",iZ-1,tab.get("X_nuc",iZ-1)+nuclei[i].n*nuclei[i].A);
      }
    }
    for(int i=1;i<=A_max;i++) {
      tab.set("X_nuc",i-1,tab.get("X_nuc",i-1)/nB);
    }
    
    // Output to file
    hdf_file hfx;
    hfx.open_or_create("dist.o2");
    hdf_output(hfx,(const table3d &)t3d,"dist");
    hdf_output(hfx,tab,"dist_sum");
    hfx.close();
  }
  
  return 0;
}

int eos_nuclei::muses(std::vector<std::string> &sv,
			     bool itive_com){
  if (sv.size()<4) {
    cerr << "Not enough arguments in point-nuclei." << endl;
    return 2;
  }

  YAML::Node node, _baseNode = YAML::LoadFile("api/output/output.yaml"); 

  double nB=o2scl::function_to_double(sv[1]);
  _baseNode["Output_UTK"]["n_B"] = nB; 
  double Ye=o2scl::function_to_double(sv[2]);
  _baseNode["Output_UTK"]["Y_e"] = Ye;
  double T=o2scl::function_to_double(sv[3])/hc_mev_fm;
  _baseNode["Output_UTK"]["T"] = T*hc_mev_fm;

  double log_xn, log_xp;
  size_t nuc_Z1, nuc_N1;
  int A_min, A_max, NmZ_min, NmZ_max;
  bool guess_provided=false;
  int flag=0;
  size_t inB=0, iYe=0, iT=0;
  vector<size_t> ix, ix_left, ix_right;

  // -------------------------------------------------------------------
  // If a table is loaded, then fix the density, electron fraction,
  // and temperature to a grid point
  
  if (loaded) {

    inB=vector_lookup(n_nB2,nB_grid2,nB);
    nB=nB_grid2[inB];
    iYe=vector_lookup(n_Ye2,Ye_grid2,Ye);
    Ye=Ye_grid2[iYe];
    iT=vector_lookup(n_T2,T_grid2,T*hc_mev_fm);
    T=T_grid2[iT]/hc_mev_fm;
    ix={inB,iYe,iT};
    if (inB>0) {
      ix_left={inB-1,iYe,iT};
    }
    if (inB<n_nB2-1) {
      ix_right={inB+1,iYe,iT};
    }
    
    flag=((int)(tg_flag.get(ix)+1.0e-10));
    
    cout << "Found grid point (inB,iYe,iT)=(" << inB << ","
	 << iYe << "," << iT << ")\n\t(nB,Ye,T[MeV])=("
	 << nB << "," << Ye << "," << T*hc_mev_fm << "), flag="
	 << flag << "." << endl;
  } else {
    cout << "No table loaded. Using exact (nB,Ye,T) specified." << endl;
  }

  // -------------------------------------------------------------------
  // Attempt to find a good initial guess, either from the command
  // line or from the table
  
  if ((alg_mode==2 || alg_mode==3 || alg_mode==4) && sv.size()>=10) {
    
    cout << "Function point_nuclei(): "
	 << "Reading guess (log_xn,log_xp,A_min,A_max,NmZ_min,NmZ_max) "
	 << "from command line." << endl;
    
    log_xn=o2scl::function_to_double(sv[4]);
    log_xp=o2scl::function_to_double(sv[5]);
    A_min=((size_t)(o2scl::function_to_double(sv[6])*(1.0+1.0e-12)));
    A_max=((size_t)(o2scl::function_to_double(sv[7])*(1.0+1.0e-12)));
    NmZ_min=((size_t)(o2scl::function_to_double(sv[8])*(1.0+1.0e-12)));
    NmZ_max=((size_t)(o2scl::function_to_double(sv[9])*(1.0+1.0e-12)));
    cout << "  " << log_xn << " " << log_xp << " " << A_min << " "
	 << A_max << " " << NmZ_min << " " << NmZ_max << endl;
    guess_provided=true;
    
  } else if ((alg_mode==0 || alg_mode==1) && sv.size()>=8) {
    
    cout << "Function point_nuclei(): "
	 << "Reading guess (log_xn,log_xp,N,Z) from command line." << endl;
    
    log_xn=o2scl::function_to_double(sv[4]);
    log_xp=o2scl::function_to_double(sv[5]);
    nuc_Z1=((size_t)(o2scl::function_to_double(sv[6])*(1.0+1.0e-12)));
    nuc_N1=((size_t)(o2scl::function_to_double(sv[7])*(1.0+1.0e-12)));
    cout << "  " << log_xn << " "
	 << log_xp << " " << nuc_Z1 << " " << nuc_N1 << endl;
    guess_provided=true;
    
  } else if (loaded) {
    
    // If the flag is 'guess' or 'in progress with guess'
    if (flag>=5 || flag==-5) {
      
      if (alg_mode==0 || alg_mode==1) {
	cout << "Function point_nuclei(): "
	     << "Reading guess (N,Z) from current table." << endl;
	log_xn=tg_log_xn.get(ix);
	log_xp=tg_log_xp.get(ix);
	nuc_Z1=((size_t)(tg_Z.get(ix)));
	nuc_N1=((size_t)(tg_A.get(ix)))-nuc_Z1;
	cout << "  " << log_xn << " "
	     << log_xp << " " << nuc_Z1 << " " << nuc_N1 << endl;
	guess_provided=true;
      } else {
	cout << "Function point_nuclei():\n"
	     << "  Reading guess (log_xn,log_xp,A_min,A_max,NmZ_min,NmZ_max) "
	     << "from current table." << endl;
	log_xn=tg_log_xn.get(ix);
	log_xp=tg_log_xp.get(ix);
	A_min=((int)(tg_A_min.get(ix)));
	A_max=((int)(tg_A_max.get(ix)));
	NmZ_min=((int)(tg_NmZ_min.get(ix)));
	NmZ_max=((int)(tg_NmZ_max.get(ix)));
	cout << "  " << log_xn << " " << log_xp << " " << A_min << " "
	     << A_max << " " << NmZ_min << " " << NmZ_max << endl;
	guess_provided=true;
      }
      
    } else if (inB>0 && tg_flag.get(ix_left)>9.9 && alg_mode==4) {

      // Otherwise if there is a good guess from the adjacent point
      // at lower baryon density
      
      cout << "Reading guess from lower baryon density point." << endl;
      log_xn=tg_log_xn.get(ix_left);
      log_xp=tg_log_xp.get(ix_left);
      A_min=((int)(tg_A_min.get(ix_left)));
      A_max=((int)(tg_A_max.get(ix_left)));
      NmZ_min=((int)(tg_NmZ_min.get(ix_left)));
      NmZ_max=((int)(tg_NmZ_max.get(ix_left)));
      cout << "  log_xn,log_xp,Z,A: "
	   << log_xn << " " << log_xp << " ";
      cout << tg_Z.get(ix_left) << " ";
      cout << tg_A.get(ix_left) << endl;
      cout << "  A_min,A_Max,NmZ_min,NmZ_max: ";
      cout << A_min << " "
	   << A_max << " " << NmZ_min << " " << NmZ_max << endl;
      guess_provided=true;
    }
  }

  // -------------------------------------------------------------------
  // If no guess was found, then just use a default guess
  
  if (guess_provided==false) {
    cout << "Function point_nuclei(): "
	 << "Using hard-coded initial guess." << endl;
    log_xn=-1.0;
    log_xp=-1.0;
    if (alg_mode==2 || alg_mode==3) {
      A_min=20;
      A_max=40;
      NmZ_min=-5;
      NmZ_max=5;
    } else if (alg_mode==4) {
      A_min=5;
      A_max=fd_A_max;
      NmZ_min=-200;
      NmZ_max=200;
    } else {
      nuc_Z1=28;
      nuc_N1=28;
    }
  }

  // -------------------------------------------------------------------
  // If necessary and possible, compute the point
  
  int ret=-1;
  
  if (flag<10 || recompute==true) {

    if (flag<10) {
      cout << "Computing point because flag < 10." << endl;
    } else {
      cout << "Computing point because recompute is true." << endl;
    }

    if (loaded && inB>0 && inB<n_nB2-1 &&
	tg_flag.get(ix_left)>9.9 &&
	tg_flag.get(ix_right)>9.9) {
      cout << "Automatically setting ranges for log_xn and log_xp "
           << "from neighboring points." << endl;
      fd_rand_ranges.resize(4);
      fd_rand_ranges[0]=tg_log_xn.get(ix_left);
      fd_rand_ranges[1]=tg_log_xn.get(ix_right);
      fd_rand_ranges[2]=tg_log_xp.get(ix_left);
      fd_rand_ranges[3]=tg_log_xp.get(ix_right);
      cout << "  ";
      vector_out(cout,fd_rand_ranges,true);
    }
    
    thermo thx;
    double mun_full, mup_full;
    
    double Zbar, Nbar;
    map<string,double> vdet;
    if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
      cout << "alg mode: " << alg_mode << endl;
      ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
			thx,mun_full,mup_full,
			A_min,A_max,NmZ_min,NmZ_max,vdet,
			true,false);
    } else {
      ret=eos_vary_ZN(nB,Ye,T,log_xn,log_xp,nuc_Z1,nuc_N1,
		      thx,mun_full,mup_full,false);
      Nbar=nuc_N1;
      Zbar=nuc_Z1;
    }

    ubvector X;
    if (ret==0) {
      compute_X(nB,X);
    }

    if (ret!=0) {
      cout << "Point failed." << endl;
    } else if (loaded) {
      cout << "Point succeeded. Storing." << endl;

      store_point(inB,iYe,iT,nB,Ye,T,thx,log_xn,log_xp,Zbar,Nbar,
		  mun_full,mup_full,X,A_min,A_max,NmZ_min,NmZ_max,
		  10.0,vdet);
    } else {
      cout << "Point succeeded." << endl;
    }

    if (ret==0) {
      cout << "log_xn: " << log_xn << endl;
      cout << "log_xp: " << log_xp << endl;
      _baseNode["Output_UTK"]["fint"] = thx.ed-T*thx.en;
      cout << "Fint: " << (thx.ed-T*thx.en)/nB*hc_mev_fm << " MeV" << endl;
      _baseNode["Output_UTK"]["A"] = Zbar+Nbar;
      _baseNode["Output_UTK"]["Z"] = Zbar;
      
      if (true) {
        
        double n0=0.16;
        double κ=1.0-nB/n0;
        double n_fraction, p_fraction, xn, xp, ξ;
        
        if (nB<0.16) {
          
          xn=0.0;
          if (log_xn>-300.0) {
            xn=pow(10.0,log_xn);
          }
          
          xp=0.0;
          if (log_xp>-300.0) {
            xp=pow(10.0,log_xp);
          }
          
          ξ=κ/(1.0-nB*xn/n0-nB*xp/n0);
          n_fraction=xn*ξ;
          p_fraction=xp*ξ;
          
        } else {
          
          xn=1.0-Ye;
          xp=Ye;
          ξ=1.0;
          n_fraction=1.0-Ye;
          p_fraction=Ye;
          
        }

        _baseNode["Output_UTK"]["ξ"] =  ξ;
        cout << "n_{n,gas} [in 1/fm^3; called n_n^{\\prime} in "
             << "Du et al. (2022)]: " << xn*nB << endl;
        cout << "n_{p,gas} [in 1/fm^3; called n_n^{\\prime} in "
             << "Du et al. (2022)]: " << xp*nB << endl;
        cout << "n_{n,avg} [in 1/fm^3: equal to Xn*nB]: "
             << xn*nB*ξ << endl;
        cout << "n_{n,avg} [in 1/fm^3: equal to Xn*nB]: "
             << xp*nB*ξ << endl;
        
        if (inc_hrg) {
          
          table_units<> t;
          store_hrg(vdet["mun_gas"]+neutron.m,
                    vdet["mup_gas"]+proton.m,xn*nB,xp*nB,T,t);
          hdf_file hf;
          hf.open_or_create("hrg.o2");
          hdf_output(hf,t,"hrg");
          hf.close();
          
        }

      }
      
      
      if (include_detail) {
	_baseNode["Output_UTK"]["zn"] =  vdet["zn"];
	_baseNode["Output_UTK"]["zp"] =  vdet["zp"];
	_baseNode["Output_UTK"]["F1"] =  vdet["F1"];
	_baseNode["Output_UTK"]["F2"] =  vdet["F2"];
	_baseNode["Output_UTK"]["F3"] =  vdet["F3"];
	_baseNode["Output_UTK"]["F4"] =  vdet["F4"];
	_baseNode["Output_UTK"]["msn"] =  vdet["msn"];
	_baseNode["Output_UTK"]["msp"] =  vdet["msp"];
	_baseNode["Output_UTK"]["Un"] =  vdet["Un"];
	_baseNode["Output_UTK"]["Up"] =  vdet["Up"];
	_baseNode["Output_UTK"]["g"] =  vdet["g"];
	_baseNode["Output_UTK"]["dgdT"] =  vdet["dgdT"];
      }
    }      

    if (fd_rand_ranges.size()>0) fd_rand_ranges.clear();
    
  }

  // (At this point in the code, the integer 'flag' is the old
  // flag, not the new flag, which is stored in tg_flag)

  if (flag>=10 && recompute==false) {
    cout << "Point already computed and recompute is false." << endl;
  }

  // -------------------------------------------------------------------
  // Print out the results when tables have been loaded and if either
  // the computation worked (ret==0) or if the old flag was marked
  // done (flag==10).
  
  if ((ret==0 || flag==10) && loaded) {
    
    _baseNode["Output_UTK"]["flag"] = tg_flag.get(ix);
    _baseNode["Output_UTK"]["log_xn"] = tg_log_xn.get(ix);
    _baseNode["Output_UTK"]["log_xp"] = tg_log_xp.get(ix);
    _baseNode["Output_UTK"]["Z"] = tg_Z.get(ix);
    _baseNode["Output_UTK"]["A"] = tg_A.get(ix);
    _baseNode["Output_UTK"]["Fint"] = tg_Fint.get(ix);
    _baseNode["Output_UTK"]["Sint"] = tg_Sint.get(ix);
    _baseNode["Output_UTK"]["Eint"] = tg_Eint.get(ix);
    _baseNode["Output_UTK"]["A_min"] = tg_A_min.get(ix);
    _baseNode["Output_UTK"]["A_max"] = tg_A_max.get(ix);
    _baseNode["Output_UTK"]["Z_min"] = tg_NmZ_min.get(ix);
    _baseNode["Output_UTK"]["Z_max"] = tg_NmZ_max.get(ix);
    _baseNode["Output_UTK"]["Xn"] = tg_Xn.get(ix);
    _baseNode["Output_UTK"]["Xp"] = tg_Xp.get(ix);
    _baseNode["Output_UTK"]["Xd"] = tg_Xd.get(ix);
    _baseNode["Output_UTK"]["Xt"] = tg_Xt.get(ix);
    _baseNode["Output_UTK"]["XHe3"] = tg_XHe3.get(ix);
    _baseNode["Output_UTK"]["XLi4"] = tg_XLi4.get(ix);
    _baseNode["Output_UTK"]["Xalpha"] = tg_Xalpha.get(ix);
    _baseNode["Output_UTK"]["Xnuclei"] = tg_Xnuclei.get(ix);
    if (rmf_fields) {
      _baseNode["Output_UTK"]["sigma"] = tg_sigma.get(ix);
      _baseNode["Output_UTK"]["omega"] = tg_omega.get(ix);
      _baseNode["Output_UTK"]["rho"] = tg_rho.get(ix);
    }
    if (include_muons) {
      _baseNode["Output_UTK"]["Ymu"] = tg_Ymu.get(ix);
    }
    if (include_muons || with_leptons) {
      _baseNode["Output_UTK"]["mue"] = tg_mue.get(ix);
    }
    if (derivs_computed) {
      _baseNode["Output_UTK"]["Pint"] = tg_Pint.get(ix);
      _baseNode["Output_UTK"]["mun"] = tg_mun.get(ix);
      _baseNode["Output_UTK"]["mup"] = tg_mup.get(ix);
      if (with_leptons) {
        _baseNode["Output_UTK"]["F"] = tg_F.get(ix);
        _baseNode["Output_UTK"]["E"] = tg_E.get(ix);
        _baseNode["Output_UTK"]["S"] = tg_S.get(ix);
        _baseNode["Output_UTK"]["P"] = tg_P.get(ix);
      }
    }
    if (include_detail) {
      _baseNode["Output_UTK"]["zn"] = tg_zn.get(ix);
      _baseNode["Output_UTK"]["zp"] = tg_zp.get(ix);
      _baseNode["Output_UTK"]["F1"] = tg_F1.get(ix);
      _baseNode["Output_UTK"]["F2"] = tg_F2.get(ix);
      _baseNode["Output_UTK"]["F3"] = tg_F3.get(ix);
      _baseNode["Output_UTK"]["F4"] = tg_F3.get(ix);
      _baseNode["Output_UTK"]["msn"] = tg_msn.get(ix);
      _baseNode["Output_UTK"]["msp"] = tg_msp.get(ix);
      _baseNode["Output_UTK"]["Un"] = tg_Un.get(ix);
      _baseNode["Output_UTK"]["Up"] = tg_Up.get(ix);
      _baseNode["Output_UTK"]["g"] = tg_g.get(ix);
      _baseNode["Output_UTK"]["dgdT"] = tg_dgdT.get(ix);
    }
  }

  // -------------------------------------------------------------------
  // Output the nuclear distribution
  
  if (alg_mode>=2 && show_all_nuclei) {
    
    cout << "Writing distribution to dist.o2." << endl;
    
    table3d t3d;
    A_max=((int)(tg_A_max.get(ix)+1.0e-10));
    t3d.set_xy("N",uniform_grid_end_width<double>(0.0,A_max,1.0),
	       "Z",uniform_grid_end_width<double>(0.0,A_max,1.0));
    
    // New slices
    t3d.new_slice("n_nuc");
    t3d.set_slice_all("n_nuc",0.0);
    t3d.new_slice("X_nuc");
    t3d.set_slice_all("X_nuc",0.0);
    
    // Initialize these to zero to avoid warnings, but they
    // are always set in the i=0 case below
    double n_nuc_min=0.0, X_nuc_min=0.0;
    for(size_t i=0;i<nuclei.size();i++) {
      // Compute minimums
      if (i==0) {
	n_nuc_min=nuclei[i].n;
	X_nuc_min=nuclei[i].n*nuclei[i].A/nB;
      } else {
	if (nuclei[i].n<n_nuc_min) {
	  n_nuc_min=nuclei[i].n;
	}
	if (nuclei[i].n*nuclei[i].A/nB<X_nuc_min) {
	  X_nuc_min=nuclei[i].n*nuclei[i].A/nB;
	}
      }
      // Set n_nuc and X_nuc
      t3d.set(nuclei[i].N,nuclei[i].Z,"n_nuc",nuclei[i].n);
      t3d.set(nuclei[i].N,nuclei[i].Z,"X_nuc",nuclei[i].n*nuclei[i].A/nB);
    }
    
    // Now take logs, replacing zeros with the minimum value
    t3d.new_slice("log10_n_nuc");
    t3d.set_slice_all("log10_n_nuc",log10(n_nuc_min));
    t3d.new_slice("log10_X_nuc");
    t3d.set_slice_all("log10_X_nuc",log10(X_nuc_min));
    for(size_t i=0;i<nuclei.size();i++) {
      t3d.set(nuclei[i].N,nuclei[i].Z,"log10_n_nuc",log10(nuclei[i].n));
      t3d.set(nuclei[i].N,nuclei[i].Z,"log10_X_nuc",log10(nuclei[i].n/nB));
    }
    
    // Now create a table adding up isotopes
    table<> tab;
    tab.line_of_names("Z X_nuc");
    for(int i=1;i<=A_max;i++) {
      double line[2]={((double)(i+1)),0.0};
      tab.line_of_data(2,line);
    }
    for(size_t i=0;i<nuclei.size();i++) {
      int iZ=nuclei[i].Z;
      if (iZ>=1 && iZ<=A_max) {
	tab.set("X_nuc",iZ-1,tab.get("X_nuc",iZ-1)+nuclei[i].n*nuclei[i].A);
      }
    }
    for(int i=1;i<=A_max;i++) {
      tab.set("X_nuc",i-1,tab.get("X_nuc",i-1)/nB);
    }
    
    // Output to file
    hdf_file hfx;
    hfx.open_or_create("dist.o2");
    hdf_output(hfx,(const table3d &)t3d,"dist");
    hdf_output(hfx,tab,"dist_sum");
    hfx.close();
  }
  std::ofstream fout("api/output/output.yaml"); 
  fout << _baseNode; 

  std::ofstream myfile;
  myfile.open("example.csv", std::ofstream::out | std::ofstream::app);
  myfile << T*hc_mev_fm << "," << tg_mun.get(ix)+neutron.m*hc_mev_fm << "," << 0 
      << "," << tg_mup.get(ix)+proton.m*hc_mev_fm-tg_mun.get(ix)-neutron.m*hc_mev_fm << "," 
      << nB <<"," << 0 << "," << Ye << "," << tg_E.get(ix) << "," << tg_P.get(ix) << "," << tg_S.get(ix) << std::endl;
  myfile.close();
  return 0;
}

int eos_nuclei::muses_table(std::vector<std::string> &sv,
			     bool itive_com) {
  std::vector<double> nB_grid; size_t n_nB;
  hdf_file hf;
  hf.open("data/fid_3_5_22.o2");
  hf.getd_vec("nB_grid",nB_grid);
  hf.get_szt("n_nB",n_nB);
  hf.close();
  double Ye=0.5;
  double T=0.1;
  std::ofstream myfile;
  myfile.open("example.csv");
  myfile.clear();
  myfile.close();
  for (size_t i=0;i<n_nB;i++){
    double nB=nB_grid[i];
    vector<string> sv2={"",o2scl::dtos(nB),
          o2scl::dtos(Ye),o2scl::dtos(T)};
    muses(sv2,false);
  }

  return 0;
}

int eos_nuclei::create_new_table(std::vector<std::string> &sv,
				 bool itive_com){
  std::vector<double> nB_grid; size_t n_nB;
  hdf_file hf;
  hf.open("data/fid_3_5_22.o2");
  hf.getd_vec("nB_grid",nB_grid);
  hf.get_szt("n_nB",n_nB);
  hf.close();
  double Ye=0.5;
  double T=0.1/hc_mev_fm;
  double log_xn=-4.661752e+01;
  double log_xp=-1.509323e+01;

  thermo thx, th_gas;
  // Chemical potentials for homogeneous matter
  double mun_gas, mup_gas;
  double mun_full, mup_full;
  double Zbar, Nbar;
  map<string,double> vdet;
  int A_min, A_max, NmZ_min, NmZ_max;
  ubvector xt(2), yt(2);
  double fnb1, fnb2, fye1, fye2;

  A_min=5;
  A_max=fd_A_max;
  NmZ_min=-200;
  NmZ_max=200;

  std::ofstream myfile;
  myfile.open("example.csv");

  for (size_t i=0;i<n_nB;i++){
    double nB=nB_grid[i];
    std::cout << "Command 'create-new-table' computing EOS at (w/o rest mass)\n"
            << "  nB(1/fm),Ye(1/fm),T(1/fm): " << nB << " "
            << Ye << " " << T << endl;

    cout << "Solving for matter at initial point nB,Ye,T(1/fm): "
       << nB << " " << Ye << " " << T << endl;
    int ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
                        thx,mun_full,mup_full,
                        A_min,A_max,NmZ_min,NmZ_max,vdet,
                        true,false);
    if (ret!=0) {
      cerr << "Initial point failed in 'create-new-table'." << endl;
      return 1;
    }

    ubvector x(2), y(2);
    x[0]=log_xn;
    x[1]=log_xp;
  
    mm_funct func=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector&,double,double,
		     double,int,double &,double &,thermo &,
		     map<string,double> &)>
      (&eos_nuclei::solve_nuclei),this,std::placeholders::_1,
      std::placeholders::_2,std::placeholders::_3,nB,Ye,T,0,
      std::ref(mun_gas),std::ref(mup_gas),std::ref(th_gas),
      std::ref(vdet));

    /*if (func(10,x,y)!=0) {
      cerr << "Initial guess to solver failed in 'point-nuclei-mu'."
         << endl;
      return 2;
    }*/
          
    mroot_hybrids<> mh2;
    mh2.verbose=1;
    mh2.msolve(2,x,func);

    cout << "Solving for matter at final point nB,Ye,T(1/fm): "
       << nB << " " << Ye << " " << T << endl;
    ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
                        thx,mun_full,mup_full,
                        A_min,A_max,NmZ_min,NmZ_max,vdet,
                        true,false);
    if (ret!=0) {
      cerr << "Final point failed in 'point-nuclei-mu'." << endl;
      return 3;
    }
    
    xt[0]=log_xn;
    xt[1]=log_xp;
    ret=solve_nuclei(2,xt,yt,nB*(1.0-1.0e-4),Ye,T,0,mun_gas,mup_gas,
                   th_gas,vdet);
    if (ret!=0) {
      cerr << "First solve_nuclei function failed." << endl;
      return ret;
    }
    fnb1=compute_fr_nuclei(nB,Ye,T,xt[0],xt[1],thx,th_gas);

    ret=solve_nuclei(2,xt,yt,nB*(1.0+1.0e-4),Ye,T,0,mun_gas,mup_gas,
                   th_gas,vdet);
    if (ret!=0) {
      cerr << "Second solve_nuclei function failed." << endl;
      return ret;
    }
    fnb2=compute_fr_nuclei(nB,Ye,T,xt[0],xt[1],thx,th_gas);

    ret=solve_nuclei(2,xt,yt,nB,Ye*(1.0-1.0e-4),T,0,mun_gas,mup_gas,
                   th_gas,vdet);
    if (ret!=0) {
      cerr << "Third solve_nuclei function failed." << endl;
      return ret;
    }
    fye1=compute_fr_nuclei(nB,Ye,T,xt[0],xt[1],thx,th_gas);

    ret=solve_nuclei(2,xt,yt,nB,Ye*(1.0+1.0e-4),T,0,mun_gas,mup_gas,
                   th_gas,vdet);
    if (ret!=0) {
      cerr << "Fourth solve_nuclei function failed." << endl;
      return ret;
    }
    fye2=compute_fr_nuclei(nB,Ye,T,xt[0],xt[1],thx,th_gas);
  
    double dfdnB=(fnb2-fnb1)/(nB*(1.0+2.0e-4));
    double dfdYe=(fye2-fye1)/(Ye*(1.0+2.0e-4));
  
    double mun2=dfdnB-dfdYe*Ye/nB;
    double mup2=dfdnB-dfdYe*(Ye-1.0)/nB;

    double Fint=(thx.ed-T*thx.en)/nB*hc_mev_fm;
    double Eint=thx.ed/nB*hc_mev_fm;
    double Sint=thx.en/nB;
    double Pint=-Fint*nB+nB*(1-Ye)*mun2+nB*Ye*mup2;

    eos_sn_base esb;
    esb.include_muons=include_muons;     
	
	  thermo lep;
	  double mue;
	  esb.compute_eg_point(nB,Ye,T,lep,mue);

	  double mu_e=mue*hc_mev_fm;
	  double E=Eint+lep.ed/nB*hc_mev_fm;
	  double P=Pint+lep.pr*hc_mev_fm;
	  double S=Sint+lep.en/nB;
      myfile << T << "," << (mun_gas+neutron.m)*hc_mev_fm << "," << 0 
      << "," << (mup_gas+proton.m-mun_gas-neutron.m)*hc_mev_fm << "," 
      << nB <<"," << 0 << "," << Ye << "," << E << "," << P << "," << S << std::endl;

  }
  myfile.close();
  return 0;
}

int eos_nuclei::increase_density(std::vector<std::string> &sv,
				 bool itive_com) {

  if (sv.size()<8) {
    cerr << "Need nB1 nB2 Ye1 Ye2 T1 T2 output_file." << endl;
    return 2;
  }
  
  double nB_start=o2scl::function_to_double(sv[1]);
  double nB_end=o2scl::function_to_double(sv[2]);
  double Ye_start=o2scl::function_to_double(sv[3]);
  double Ye_end=o2scl::function_to_double(sv[4]);
  double T_start=o2scl::function_to_double(sv[5])/hc_mev_fm;
  double T_end=o2scl::function_to_double(sv[6])/hc_mev_fm;
  std::string out_file=sv[7];

  size_t inB_start=vector_lookup(n_nB2,nB_grid2,nB_start);
  size_t iYe_start=vector_lookup(n_Ye2,Ye_grid2,Ye_start);
  size_t iT_start=vector_lookup(n_T2,T_grid2,T_start*hc_mev_fm);

  size_t inB_end=vector_lookup(n_nB2,nB_grid2,nB_end);
  size_t iYe_end=vector_lookup(n_Ye2,Ye_grid2,Ye_end);
  size_t iT_end=vector_lookup(n_T2,T_grid2,T_end*hc_mev_fm);

  double log_xn, log_xp;
  size_t nuc_Z1, nuc_N1;
  int A_min, A_max, NmZ_min, NmZ_max;

  derivs_computed=false;
  with_leptons=false;

  for(size_t iT=iT_start;iT<=iT_end;iT++) {
    
    for(size_t iYe=iYe_start;iYe<=iYe_end;iYe++) {

      vector<size_t> ix_start={inB_start,iYe,iT};
      
      log_xn=tg_log_xn.get(ix_start);
      log_xp=tg_log_xp.get(ix_start);
      A_min=tg_A_min.get(ix_start);
      A_max=tg_A_max.get(ix_start);
      NmZ_min=tg_NmZ_min.get(ix_start);
      NmZ_max=tg_NmZ_max.get(ix_start);

      bool no_nuclei=false;

      vector<double> Zarr, Narr;
      
      for(size_t inB=inB_start;inB<=inB_end;inB++) {
	
	double nB=nB_grid2[inB];
	double Ye=Ye_grid2[iYe];
	double T=T_grid2[iT]/hc_mev_fm;
        vector<size_t> ix={inB,iYe,iT};
	
	thermo thx;
	double mun_full, mup_full;
	
	if (alg_mode!=4) {
	  cout << "Only works for alg_mode=4." << endl;
	  return 1;
	}

	// If the next point will take us outside the limiting N/Z,
	// then just set no_nuclei to true
	if (no_nuclei==false && Zarr.size()>=2 && Narr.size()>=2) {
	  double last_Z=Zarr[Zarr.size()-1];
	  double last_N=Narr[Narr.size()-1];
	  double last_ratio=last_N/last_Z;
	  double sec_last_Z=Zarr[Zarr.size()-2];
	  double sec_last_N=Narr[Narr.size()-2];
	  double sec_last_ratio=sec_last_N/sec_last_Z;
	  double dratio=last_ratio-sec_last_ratio;
	  if (last_ratio>sec_last_ratio && last_ratio+dratio>max_ratio-1.0) {
	    cout << "Predicted ratio: " << last_ratio+dratio << " "
		 << "max_ratio: " << max_ratio << endl;
	    if (last_ratio+dratio>max_ratio-0.1) {
	      cout << "Predicted ratio large, setting no_nuclei to true."
		   << endl;
	      no_nuclei=true;
	    }
	  }
	}
	
	double Zbar, Nbar;
	map<string,double> vdet;
	int ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
			      thx,mun_full,mup_full,
			      A_min,A_max,NmZ_min,NmZ_max,
			      vdet,true,no_nuclei);
	if (Zbar>0.0 && Nbar>0.0) {
	  cout << "ret,nB,Ye,T_MeV,Z,N,N/Z: "
	       << ret << " " << nB << " " << Ye << " " << T*hc_mev_fm << " "
	       << Zbar << " " << Nbar << " " << Nbar/Zbar << endl;
	} else {
	  cout << "ret,nB,Ye,T_MeV,Z,N: "
	       << ret << " " << nB << " " << Ye << " " << T*hc_mev_fm << " "
	       << Zbar << " " << Nbar << endl;
	}

	Zarr.push_back(Zbar);
	Narr.push_back(Nbar);
	
	if (no_nuclei==false && Zbar<1.0e-6 && Nbar<1.0e-6) {
	  cout << "Function eos_vary_dist() found no nuclei preferred.\n"
	       << "  Setting no-nuclei to true for this Ye and T."
	       << endl;
	  no_nuclei=true;
	}
	
	if (ret==0) {
	  
	  ubvector X;
	  compute_X(nB,X);

	  if (tg_A.get(ix)>0.0) {
	    cout << "before: " << tg_Z.get(ix) << " "
		 << tg_A.get(ix) << " "
		 << tg_A.get(ix)/
	      tg_Z.get(ix) << endl;
	  } else {
	    cout << "before: " << tg_Z.get(ix) << " "
		 << tg_A.get(ix) << " " << 0.0 << endl;
	  }
	  
	  store_point(inB,iYe,iT,nB,Ye,T,thx,log_xn,log_xp,Zbar,Nbar,
		      mun_full,mup_full,X,A_min,A_max,NmZ_min,
		      NmZ_max,10.0,vdet);
	  
	  if (tg_A.get(ix)>0.0) {
	    cout << "after: " << tg_log_xn.get(ix) << " "
		 << tg_log_xp.get(ix) << " "
		 << tg_Xnuclei.get(ix) << " "
		 << tg_Z.get(ix) << " "
		 << tg_A.get(ix) << " "
		 << tg_A.get(ix)/
	      tg_Z.get(ix) << endl;
	  } else {
	    cout << "after: " << tg_log_xn.get(ix) << " "
		 << tg_log_xp.get(ix) << " "
		 << tg_Xnuclei.get(ix) << " "
		 << tg_Z.get(ix) << " "
		 << tg_A.get(ix) << " " << 0.0 << endl;
	  }
	  cout << endl;
	  
	} else {
	  // If we failed, presume that we're past the transition
	  if (nB>0.04) {
	    no_nuclei=true;
	  } else {
	    inB=inB_end;
	  }
	}
	
      }
    }

    // We require an output file as it allows the user to specify
    // different files when multiple runs are occurring at the
    // same time
    write_results(out_file);
    
  }
  
  return 0;
}

int eos_nuclei::fix_cc(std::vector<std::string> &sv,
		       bool itive_com) {

  if (sv.size()<2) {
    cerr << "Need output_file." << endl;
    return 2;
  }
  
  std::string out_file=sv[1];

  size_t inB_start=vector_lookup(n_nB2,nB_grid2,0.01);
  size_t inB_end=vector_lookup(n_nB2,nB_grid2,0.16);

  double log_xn, log_xp;
  size_t nuc_Z1, nuc_N1;
  int A_min, A_max, NmZ_min, NmZ_max;

  derivs_computed=false;
  with_leptons=false;
  bool no_nuclei;
  
  for(size_t iT=0;iT<n_T2;iT++) {

    bool found_point=false;
    
    for(size_t iYe=0;iYe<n_Ye2;iYe++) {
  
      for(size_t inB=inB_start;inB<=inB_end;inB++) {
	
	double nB=nB_grid2[inB];
	double Ye=Ye_grid2[iYe];
	double T=T_grid2[iT]/hc_mev_fm;

        vector<size_t> ix={inB,iYe,iT};
        
	if (tg_flag.get(ix)<9.9) {
	  found_point=true;
	  
	  cout << "Found uncomputed point at (nB,Ye,T):\n"
	       << "\t" << nB << " " << Ye << " "
	       << T*hc_mev_fm << endl;
	  
	  no_nuclei=false;
	  
	  if (inB>20 && no_nuclei==false) {
	    for(size_t j=1;j<20 && no_nuclei==false;j++) {
	      cout << "j: " << j << endl;
              vector<size_t> ixmj={inB-j,iYe,iT};
	      if (tg_flag.get(ixmj)>9.9 && tg_Z.get(ixmj)<1.0e-4 &&
		  tg_A.get(ixmj)<1.0e-4) {
		cout << "Found nuclear matter at a lower density."
		     << endl;
		no_nuclei=true;
	      }
	    }
	  }
	  
          vector<size_t> ixm1={inB-1,iYe,iT};
          
	  if (no_nuclei==false && inB>2) {
            vector<size_t> ixm2={inB-2,iYe,iT};
	    double Z_prev=tg_Z.get(ixm1);
	    double Z_prev2=tg_Z.get(ixm2);
	    if (Z_prev>0 && Z_prev2>0) {
	      double A_prev=tg_A.get(ixm1);
	      double A_prev2=tg_A.get(ixm2);
	      double NoZ_prev=(A_prev-Z_prev)/Z_prev;
	      double NoZ_prev2=(A_prev2-Z_prev2)/Z_prev2;
	      cout << "prev: " << NoZ_prev << " " << NoZ_prev2 << endl;
	      if (NoZ_prev>NoZ_prev2) {
		double dNoZ=NoZ_prev-NoZ_prev2;
		if (NoZ_prev+dNoZ>max_ratio) {
		  no_nuclei=true;
		  cout << "Extrapolated Z/N is too large." << endl;
		}
	      }
	    }
	  }
	  
	  // Get initial guess 
	  log_xn=tg_log_xn.get(ixm1);
	  log_xp=tg_log_xp.get(ixm1);
	  A_min=tg_A_min.get(ixm1);
	  A_max=tg_A_max.get(ixm1);
	  NmZ_min=tg_NmZ_min.get(ixm1);
	  NmZ_max=tg_NmZ_max.get(ixm1);
	  
	  thermo thx;
	  double mun_full, mup_full;
	  
	  if (alg_mode!=4) {
	    cout << "Only works for alg_mode=4." << endl;
	    return 1;
	  }
	  
	  double Zbar, Nbar;
	  map<string,double> vdet;
	  int ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
				thx,mun_full,mup_full,
				A_min,A_max,NmZ_min,NmZ_max,
				vdet,true,no_nuclei);
	  if (ret!=0) {
	    cout << "Function eos_vary_dist() failed, going "
		 << "without nuclei." << endl;
	    no_nuclei=true;
	    ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
			      thx,mun_full,mup_full,
			      A_min,A_max,NmZ_min,NmZ_max,
			      vdet,true,no_nuclei);
	  }
	  
	  if (ret==0) {
	    
	    ubvector X;
	    compute_X(nB,X);
	    
	    store_point(inB,iYe,iT,nB,Ye,T,thx,log_xn,log_xp,Zbar,Nbar,
			mun_full,mup_full,X,A_min,A_max,NmZ_min,
			NmZ_max,10.0,vdet);
	    
	    if (tg_A.get(ix)>0.0) {
	      cout << "Success (log_xn,log_xp,Xnuc,Z,N,N/Z): "
		   << tg_log_xn.get(ix) << " "
		   << tg_log_xp.get(ix) << " "
		   << tg_Xnuclei.get(ix) << " "
		   << tg_Z.get(ix) << " "
		   << tg_A.get(ix)-tg_Z.get(ix) << " "
		   << (tg_A.get(ix)-tg_Z.get(ix))/
		tg_Z.get(ix) << endl;
	    } else {
	      cout << "Success (log_xn,log_xp,Xnuc,Z,N): "
		   << tg_log_xn.get(ix) << " "
		   << tg_log_xp.get(ix) << " "
		   << tg_Xnuclei.get(ix) << " "
		   << tg_Z.get(ix) << " "
		   << tg_A.get(ix) << endl;
	    }
	    cout << endl;
	    
	  }

	}
	
      }
      
    }

    if (found_point==true) {
      // We require an output file as it allows the user to specify
      // different files when multiple runs are occurring at the
      // same time
      write_results(out_file);
    }

  }
  return 0;
}

int eos_nuclei::verify(std::vector<std::string> &sv,
		       bool itive_com) {

  if (sv.size()<2) {
    cerr << "Need mode." << endl;
    return 2;
  }
  
  std::string mode=sv[1], out_file=sv[2];
  size_t n_tot=n_nB2*n_Ye2*n_T2;
  if (mode=="random" || mode=="random_lg") {
    n_tot=stoszt(sv[2]);
    out_file=sv[3];
  }
  if (mode=="point") n_tot=1;
  
  double log_xn, log_xp;
  size_t nuc_Z1, nuc_N1;
  int A_min, A_max, NmZ_min, NmZ_max;

  derivs_computed=false;
  with_leptons=false;
  bool no_nuclei;

  size_t inB_lo=0, inB_hi=0;
  if (mode=="random_lg" || mode=="all_lg") {
    inB_lo=vector_lookup(nB_grid2,0.01);
    inB_hi=vector_lookup(nB_grid2,0.15);
    // The random int below doesn't include ihi, so we increase
    // this by one to make sure we include ihi
    if (inB_hi<n_nB2) inB_hi++;
    cout << "inB_lo, inB_hi: " << inB_lo << " " << inB_hi << endl;
  }

  if (mode=="all_lg") {
    n_tot=(inB_hi-inB_lo+1)*n_Ye2*n_T2;
  }
  cout << "n_tot: " << n_tot << endl;

  bool verify_success=true;
  
  for(size_t j=0;j<n_tot;j++) {
    
    size_t inB=0, iYe=0, iT=0;

    if (mode=="all") {
      vector<size_t> vix(tg_A.get_rank());
      tg_A.unpack_index(j,vix);
      inB=vix[0];
      iYe=vix[1];
      iT=vix[2];
    } else if (mode=="all_lg") {
      size_t n_slice=n_Ye2*n_T2;
      inB=j/n_slice;
      size_t j2=j-n_slice*inB;
      iYe=j2/n_T2;
      iT=j2-n_T2*iYe;
      inB+=inB_lo;
    } else if (mode=="point") {
      inB=vector_lookup(nB_grid2,o2scl::function_to_double(sv[2]));
      iYe=vector_lookup(Ye_grid2,o2scl::function_to_double(sv[3]));
      iT=vector_lookup(T_grid2,o2scl::function_to_double(sv[4]));
    } else {
      if (mode=="random_lg") {
	inB=rng.random_int(inB_hi-inB_lo+1)+inB_lo;
      } else {
	inB=rng.random_int(n_nB2);
      }
      iYe=rng.random_int(n_Ye2);
      iT=rng.random_int(n_T2);
    }

    double nB=nB_grid2[inB];
    double Ye=Ye_grid2[iYe];
    double T=T_grid2[iT]/hc_mev_fm;
    vector<size_t> ix={inB,iYe,iT};

    // Get initial guess 
    log_xn=tg_log_xn.get(ix);
    log_xp=tg_log_xp.get(ix);
    A_min=tg_A_min.get(ix);
    A_max=tg_A_max.get(ix);
    NmZ_min=tg_NmZ_min.get(ix);
    NmZ_max=tg_NmZ_max.get(ix);
    no_nuclei=false;
    if (tg_A.get(ix)<1.0e-6 && tg_Z.get(ix)<1.0e-6) {
      no_nuclei=true;
    }
    
    thermo thx;
    double mun_full, mup_full;
    
    if (alg_mode!=4) {
      cout << "Only works for alg_mode=4." << endl;
      return 1;
    }
    
    double Zbar, Nbar;
    map<string,double> vdet;
    
    verify_only=true;
    int ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
			  thx,mun_full,mup_full,
			  A_min,A_max,NmZ_min,NmZ_max,
			  vdet,true,no_nuclei);
    verify_only=false;

    bool computed=false;
    if (ret!=0) {
      tg_flag.get(ix)=5.0;
      cout << "Verification failed. Computing." << endl;
      ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
                        thx,mun_full,mup_full,
                        A_min,A_max,NmZ_min,NmZ_max,
                        vdet,true,no_nuclei);

      if (ret==0) {
	// Note that Zbar and Nbar aren't computed if verify_only is true
	if (fabs(Zbar-tg_Z.get(ix))>1.0e-6) {
	  cout << "Z changed from "
	       << tg_Z.get(ix) << " to " << Zbar << endl;
	}
	if (fabs(Zbar+Nbar-tg_A.get(ix))>1.0e-6) {
	  cout << "A changed from "
	       << tg_A.get(ix) << " to " << Zbar+Nbar << endl;
	}
      }

      computed=true;
    }

    cout << "verify() j,nB,Ye,T,ret: ";
    cout.width(((size_t)log10(n_tot*(10.0-1.0e-8))));
    cout << j << " " << nB << " " << Ye << " " << T*hc_mev_fm
	 << " " << ret << endl;
    
    if (ret!=0) {
      cout << "Verification failed. Stopping." << endl;
      j=n_tot;
      verify_success=false;
    }

    if (ret==0 && computed) {
      ubvector X;
      compute_X(nB,X);
      store_point(inB,iYe,iT,nB,Ye,T,thx,log_xn,log_xp,Zbar,Nbar,
		  mun_full,mup_full,X,A_min,A_max,NmZ_min,NmZ_max,
		  10.0,vdet);
    }
    
    if (((int)j)%file_update_iters==file_update_iters-1) {
      cout << "Updating file." << endl;
      write_results(out_file);
    }

  }

  if (verify_success) {
    cout << "Verification succeeded." << endl;
  }
  
  return 0;
}

int eos_nuclei::stats(std::vector<std::string> &sv,
		      bool itive_com) {

  if (loaded==false) {
    cerr << "No table loaded in stats()." << endl;
    return 3;
  }
  
  cout << endl;
  cout << "derivs_computed: " << derivs_computed << endl;
  cout << "with_leptons: " << with_leptons << endl;
  cout << "alg_mode: " << alg_mode << endl;
  cout << endl;
  
  const vector<double> &data=tg_flag.get_data();
  vector<size_t> flags(22);
  
  size_t nb_frac_count=0, S_neg_count=0, Fint_count=0, F_count=0;

  size_t ti_int_count=0, ti_count=0;
  
  vector<o2scl::tensor_grid<> *> ptrs;
  ptrs.push_back(&tg_log_xn);
  ptrs.push_back(&tg_log_xp);
  ptrs.push_back(&tg_flag);
  if (with_leptons) {
    ptrs.push_back(&tg_F);
    ptrs.push_back(&tg_E);
    ptrs.push_back(&tg_P);
    ptrs.push_back(&tg_S);
    ptrs.push_back(&tg_mue);
  }
  ptrs.push_back(&tg_Fint);
  ptrs.push_back(&tg_Eint);
  ptrs.push_back(&tg_Sint);
  if (derivs_computed) {
    ptrs.push_back(&tg_Pint);
    ptrs.push_back(&tg_mun);
    ptrs.push_back(&tg_mup);
  }
  ptrs.push_back(&tg_Z);
  ptrs.push_back(&tg_A);
  ptrs.push_back(&tg_Xn);
  ptrs.push_back(&tg_Xp);
  ptrs.push_back(&tg_Xalpha);
  ptrs.push_back(&tg_Xnuclei);
  if (include_muons) {
    ptrs.push_back(&tg_Ymu);
  }
  ptrs.push_back(&tg_Xd);
  ptrs.push_back(&tg_Xt);
  ptrs.push_back(&tg_XHe3);
  ptrs.push_back(&tg_XLi4);
  if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
    ptrs.push_back(&tg_A_min);
    ptrs.push_back(&tg_A_max);
    ptrs.push_back(&tg_NmZ_min);
    ptrs.push_back(&tg_NmZ_max);
  }
  
  for(size_t i=0;i<data.size();i++) {
    
    int iflag=((int)(data[i]*(1.0+1.0e-12)));
    if (iflag<-10 || iflag>10) {
      iflag=11;
    }
    flags[iflag+10]++;
    
    if (iflag==10) {
      
      size_t ix[3];
      tg_Xn.unpack_index(i,ix);
      
      double nB=nB_grid2[ix[0]];
      double Ye=Ye_grid2[ix[1]];
      double T_MeV=T_grid2[ix[2]];

      // Check grid
      if (true) {
        for(size_t kk=0;kk<ptrs.size();kk++) {
          if (fabs(ptrs[kk]->get_grid(0,ix[0])-nB)/nB>1.0e-6) {
            std::cout << "Grid broken (nB): " << ptrs[kk]->get_grid(0,ix[0])
                      << " " << kk << " " << nB << " " << ix[0]
                      << std::endl;
            O2SCL_ERR("Baryon density grid broken in "
                      "eos_nuclei::stats().",o2scl::exc_efailed);
          }
          if (fabs(ptrs[kk]->get_grid(1,ix[1])-Ye)/Ye>1.0e-6) {
            std::cout << "Grid broken (Ye): " << ptrs[kk]->get_grid(1,ix[1])
                      << " " << kk << " " << Ye << " " << ix[1]
                      << std::endl;
            O2SCL_ERR("Electron fraction grid broken in "
                      "eos_nuclei::stats().",o2scl::exc_efailed);
          }
          if (fabs(ptrs[kk]->get_grid(2,ix[2])-T_MeV)/T_MeV>1.0e-6) {
            std::cout << "Grid broken (T): " << ptrs[kk]->get_grid(2,ix[2])
                      << " " << kk << " " << T_MeV << " " << ix[2]
                      << std::endl;
            O2SCL_ERR("Temperature grid broken in "
                      "eos_nuclei::stats().",o2scl::exc_efailed);
          }
        }
      }
      
      // Check that X's add up to 1
      if (true) {
	double check_X=tg_Xn.get_data()[i]+tg_Xp.get_data()[i]+
	  tg_Xalpha.get_data()[i]+tg_Xnuclei.get_data()[i]+
	  tg_Xd.get_data()[i]+tg_Xt.get_data()[i]+
	  tg_XHe3.get_data()[i]+tg_XLi4.get_data()[i];
	
	if (fabs(check_X-1.0)>1.0e-6) {
	  nb_frac_count++;
          if (nb_frac_count<1000) {
            cout << "Nuclear fractions do not add up to 1 "
                 << "(i,nB,Ye,T,X_total):\n  "
                 << i << " " << nB << " " << Ye << " "
                 << T_MeV << " " << check_X << endl;
          } else if (nb_frac_count==1000) {
	    cout << "Further nuclear fractions warnings suppressed."
		 << endl;
	  }
	}
      }

      // Check that Fint=Eint-T*Sint
      if (tg_Eint.total_size()>0) {

	double scale=fabs(tg_Fint.get_data()[i]);
	if (scale<10.0) scale=10.0;
	double Fint_check=(tg_Fint.get_data()[i]-
			   tg_Eint.get_data()[i]+
			   T_MeV*tg_Sint.get_data()[i])/scale;
	if (fabs(Fint_check)>1.0e-9) {
	  Fint_count++;
          if (Fint_count<1000) {
            cout << "Fint doesn't match Eint-T*Sint:"
                 << "(i,nB,Ye,T,Fint_check):\n  "
                 << i << " " << nB << " " << Ye << " " << T_MeV
                 << " " << Fint_check << endl;
            cout << "  (Fint,Eint,T,Sint,Eint-T*Sint): "
                 << tg_Fint.get_data()[i] << " ";
            cout << tg_Eint.get_data()[i] << " ";
            cout << T_MeV << " " << tg_Sint.get_data()[i] << " ";
            cout << tg_Eint.get_data()[i]-T_MeV*tg_Sint.get_data()[i] << endl;
          } else if (Fint_count==1000) {
	    cout << "Further Fint warnings suppressed." << endl;
	  }
	}
      }

      // Check that Eint*nB=-Pint+T*Sint*nB+nn*mun+np*mup
      if (derivs_computed) {
	double nn=nB*(1.0-Ye);
	double np=nB*Ye;
	double scale=fabs(tg_Eint.get_data()[i]);
	if (scale<100.0) scale=100.0;
	double ti_int_check=(tg_Eint.get_data()[i]*nB+
			     tg_Pint.get_data()[i]-
			     T_MeV*tg_Sint.get_data()[i]*nB-
			     nn*tg_mun.get_data()[i]-
			     np*tg_mup.get_data()[i])/scale/nB;
	if (fabs(ti_int_check)>1.0e-9) {
	  ti_int_count++;
          if (ti_int_count<1000) {
            cout << "Thermodynamic identity violated for EOS w/o leptons:\n"
                 << "  i,nB,Ye,T [MeV]: " 
                 << i << " " << nB << " " << Ye << " " << T_MeV << endl;
            cout << "  ti_int_check: " << ti_int_check << endl;
	  } else if (ti_int_count==1000) {
	    cout << "Further thermodynamic identity w/o leptons "
                 << "warnings suppressed." << endl;
	  }
	}
      }

      // Check that E*nB=-P+T*S*nB+nn*mun+np*mup+ne*mue
      if (derivs_computed) {
	double nn=nB*(1.0-Ye);
	double np=nB*Ye;
	double scale=fabs(tg_E.get_data()[i]);
	if (scale<10.0) scale=10.0;
	double ti_check=(tg_E.get_data()[i]*nB+
			 tg_P.get_data()[i]-
			 T_MeV*tg_S.get_data()[i]*nB-
			 nn*tg_mun.get_data()[i]-
			 np*tg_mup.get_data()[i]-
			 np*tg_mue.get_data()[i])/scale/nB;
	if (fabs(ti_check)>1.0e-9) {
	  ti_count++;
          if (ti_count<1000) {
            cout << "Thermodynamic identity violated for EOS w/leptons:\n"
                 << "  ti_count,i,nB,Ye,T [MeV]: "
                 << ti_count << " " << i << " "
                 << nB << " " << Ye << " " << T_MeV << endl;
            cout << "  ti_check,scale: "
                 << ti_check << " " << scale << endl;
            if (false) {
              cout.precision(5);
              cout << tg_E.get_data()[i]*nB << " "
                   << tg_P.get_data()[i] << " "
                   << T_MeV*tg_S.get_data()[i]*nB << " "
                   << nn*tg_mun.get_data()[i] << " "
                   << np*tg_mup.get_data()[i] << " "
                   << np*tg_mue.get_data()[i] << endl;
              cout.precision(6);
              char ch;
              cin >> ch;
            }
          } else if (ti_count==1000) {
	    cout << "Further thermodynamic idenity w/leptons "
                 << "warnings suppressed." << endl;
	  }
	}
      }

      // Check that F=E-T*S
      if (with_leptons) {

	double scale=fabs(tg_F.get_data()[i]);
	if (scale<10.0) scale=10.0;
	double F_check=(tg_F.get_data()[i]-
			tg_E.get_data()[i]+
			T_MeV*tg_S.get_data()[i])/scale;
	if (fabs(F_check)>1.0e-9) {
	  F_count++;
          if (F_count<1000) {
            cout << "F doesn't match E-T*S (F_count,i,nB,Ye,T,F_check,scale): "
                 << F_count << " " << i << "\n  " << nB << " " << Ye
                 << " " << T_MeV
                 << " " << F_check << " " << scale << endl;
          } else if (F_count==1000) {
	    cout << "Further Fint warnings suppressed." << endl;
	  }
          exit(-1);
	}
      
      }

      // Check that Sint>0
      if (false && S_neg_count<1000 && tg_Sint.get_data()[i]<0.0) {
	cout << "Entropy per baryon negative (i,nB,Ye,T,Sint): "
	     << i << " " << nB << " " << Ye << " "
	     << T_MeV << " " << tg_Sint.get_data()[i] << endl;
	S_neg_count++;
	if (S_neg_count==1000) {
	  cout << "Further negative entropy warnings suppressed."
	       << endl;
	}
      }
      
    }
  }
  cout << "nb_frac_count: " << nb_frac_count << endl;
  //cout << "Sint_neg_count: " << S_neg_count << endl;
  if (tg_Eint.total_size()>0) {
    cout << "Fint_count: " << Fint_count << endl;
  }
  if (derivs_computed) {
    cout << "ti_int_count: " << ti_int_count << endl;
    if (with_leptons) {
      cout << "F_count: " << F_count << endl;
      cout << "ti_count: " << ti_count << endl;
    }
  }
  cout << endl;
  cout << "flag counts:" << endl;
  size_t count=0;
  for(int j=0;j<22;j++) {
    if (j!=20) count+=flags[j];
    cout.width(3);
    cout << j-10 << ": " << flags[j] << endl;
  }

  // Less than 200 points missing, so output the missing points
  if (count>0 && count<200) {
    cout << "\nOnly " << count << " non-converged points: " << endl;
    size_t temp_ctr=0;
    cout << "i nB Ye T_MeV flag: " << endl;
    for(size_t i=0;i<data.size();i++) {
      int iflag=((int)(data[i]*(1.0+1.0e-12)));
      if (iflag!=10) {
	size_t ix[3];
	tg_Xn.unpack_index(i,ix);
	double nB=nB_grid2[ix[0]];
	double Ye=Ye_grid2[ix[1]];
	double T_MeV=T_grid2[ix[2]];
	cout.width(2);
	cout << temp_ctr << " ";
	cout.width(8);
	cout << i << " " << nB << " " << Ye << " "
	     << T_MeV << " " << iflag << endl;
	temp_ctr++;
      }
    }
    
  }
  
  return 0;
}

int eos_nuclei::merge_tables(std::vector<std::string> &sv,
			     bool itive_com) {

  if (sv.size()<4) {
    cerr << "Command 'merge-tables' needs at least 3 arguments." << endl;
    cerr << "<filename> <N> <nB func> <Ye func> <T func> "
	 << "[log_xn_0 log_xp_0 Z_0 N_0]" << endl;
    return 1;
  }
  
  string in1=sv[1];
  string in2=sv[2];
  string out=sv[3];

  eos_nuclei en2;
  
  vector<size_t> counts(22);

  // Read the first file
  
  read_results(in1);

  const vector<double> &d=tg_flag.get_data();
  for(size_t i=0;i<tg_flag.total_size();i++) {
    size_t szt_tmp=((size_t)((d[i]+10.0)*(1.0+1.0e-12)));
    if (szt_tmp>20) szt_tmp=21;
    counts[szt_tmp]++;
  }
  cout << "Counts for file 1: ";
  vector_out(cout,counts,true);

  // Read the second file

  en2.read_results(in2);

  for(size_t i=0;i<22;i++) counts[i]=0;
  const vector<double> &d2=en2.tg_flag.get_data();
  for(size_t i=0;i<en2.tg_flag.total_size();i++) {
    size_t szt_tmp=((size_t)((d2[i]+10.0)*(1.0+1.0e-12)));
    if (szt_tmp>20) szt_tmp=21;
    counts[szt_tmp]++;
  }
  cout << "Counts for file 2: ";
  vector_out(cout,counts,true);

  // Check if two tables match

  if (en2.derivs_computed!=derivs_computed) {
    cerr << "derivs_computed doesn't match." << endl;
    return 1;
  }
  if (en2.with_leptons!=with_leptons) {
    cerr << "with_leptons doesn't match." << endl;
    return 2;
  }
  if (n_nB2!=en2.n_nB2 || n_Ye2!=en2.n_Ye2 || n_T2!=en2.n_T2) {
    cerr << "Grid sizes don't match." << endl;
    return 3;
  }
  if (nB_grid2!=en2.nB_grid2 || Ye_grid2!=en2.Ye_grid2 ||
      T_grid2!=en2.T_grid2) {
    cerr << "Grid don't match." << endl;
    return 4;
  }
  
  // Now perform the merge
  
  size_t nx=tg_flag.get_size(0);
  size_t ny=tg_flag.get_size(1);
  size_t nz=tg_flag.get_size(2);

  size_t c1=0, c2=0;
  
  for(size_t i=0;i<nx;i++) {
    double nB=nB_grid2[i];
    if (i%10==9) {
      cout << (i+1) << "/" << nx << endl;
    }
    for(size_t j=0;j<ny;j++) {
      double Ye=Ye_grid2[j];
      for(size_t k=0;k<nz;k++) {
        vector<size_t> ix={i,j,k};
	double T_MeV=T_grid2[k];
	double flag1=tg_flag.get(ix);
	double flag2=en2.tg_flag.get(ix);
	double Fint1=tg_Fint.get(ix);
	double Fint2=en2.tg_Fint.get(ix);

	if (fabs(flag1)<0.5 && (Fint1>1.0e90 ||
				!std::isfinite(Fint1) || Fint1<(-1.0e10))) {
	  flag1=0.0;
	  tg_flag.set(ix,0.0);
	  tg_Fint.set(ix,0.0);
	  tg_Z.set(ix,0.0);
	  tg_A.set(ix,0.0);
	  tg_log_xn.set(ix,0.0);
	  tg_log_xp.set(ix,0.0);
	  c1++;
	  cout << "Invalid free energy in table 1 (" << nB << ","
	       << Ye << "," << T_MeV << ")." << endl;
          O2SCL_ERR("Invalid Fint in table 1.",o2scl::exc_einval);
          
	}
        
	if (fabs(flag2)<0.5 && (Fint2>1.0e90 ||
				!std::isfinite(Fint2)|| Fint2<(-1.0e10))) {
          
	  flag2=0.0;
	  en2.tg_flag.set(ix,0.0);
	  en2.tg_Fint.set(ix,0.0);
	  en2.tg_Z.set(ix,0.0);
	  en2.tg_A.set(ix,0.0);
	  en2.tg_log_xn.set(ix,0.0);
	  en2.tg_log_xp.set(ix,0.0);
	  cout << "Invalid free energy in table 2 (" << nB << ","
	       << Ye << "," << T_MeV << ")." << endl;
          O2SCL_ERR("Invalid Fint in table 2.",o2scl::exc_einval);
          
	} else if ((flag1>=10.0 && flag2>=10.0 && Fint2<Fint1) ||
		   (flag1<10.0 && flag2>=10.0) ||
		   (flag1>0.0 && flag2>0.0 && flag1<10.0 && flag2<0.0 &&
		    Fint2<Fint1) ||
		   (flag2>0.0 && flag1==0.0)) {
	  /*
	    Four cases when we need to copy "2" results to "1":
	    1: They both have flags=10 but "2" has smaller free energy
	    2: "2" has a flag=10 but "1" does not
	    3: They both have flags smaller than 10 but larger than
	    zero and "2" has a smaller free energy
	    4: "2" has a non-zero flag but "1" does not
	  */

	  tg_flag.set(ix,en2.tg_flag.get(ix));

	  tg_log_xn.set(ix,en2.tg_log_xn.get(ix));
	  tg_log_xp.set(ix,en2.tg_log_xp.get(ix));
	  tg_Z.set(ix,en2.tg_Z.get(ix));
	  tg_A.set(ix,en2.tg_A.get(ix));
	  
	  tg_A_min.set(ix,en2.tg_A_min.get(ix));
	  tg_A_max.set(ix,en2.tg_A_max.get(ix));
	  tg_NmZ_min.set(ix,en2.tg_NmZ_min.get(ix));
	  tg_NmZ_max.set(ix,en2.tg_NmZ_max.get(ix));
	  
	  tg_Fint.set(ix,en2.tg_Fint.get(ix));
	  tg_Sint.set(ix,en2.tg_Sint.get(ix));
	  tg_Eint.set(ix,en2.tg_Eint.get(ix));

	  tg_Xn.set(ix,en2.tg_Xn.get(ix));
	  tg_Xp.set(ix,en2.tg_Xp.get(ix));
	  tg_Xalpha.set(ix,en2.tg_Xalpha.get(ix));
	  tg_Xnuclei.set(ix,en2.tg_Xnuclei.get(ix));
	  tg_Xd.set(ix,en2.tg_Xd.get(ix));
	  tg_Xt.set(ix,en2.tg_Xt.get(ix));
	  tg_XHe3.set(ix,en2.tg_XHe3.get(ix));
	  tg_XLi4.set(ix,en2.tg_XLi4.get(ix));

	  if (derivs_computed) {
	    tg_Pint.set(ix,en2.tg_Pint.get(ix));
	    tg_mun.set(ix,en2.tg_mun.get(ix));
	    tg_mup.set(ix,en2.tg_mup.get(ix));
	    if (with_leptons) {
	      tg_F.set(ix,en2.tg_F.get(ix));
	      tg_E.set(ix,en2.tg_E.get(ix));
	      tg_P.set(ix,en2.tg_P.get(ix));
	      tg_S.set(ix,en2.tg_S.get(ix));
	      tg_mue.set(ix,en2.tg_mue.get(ix));
	    }
	  }

	  c2++;
	}
      }
    }
  }

  cout << "Points removed from table 1: " << c1 << endl;
  cout << "Points copied from table 2 to table 1: " << c2 << endl;
  
  for(size_t i=0;i<22;i++) counts[i]=0;
  const vector<double> &d3=tg_flag.get_data();
  for(size_t i=0;i<tg_flag.total_size();i++) {
    size_t szt_tmp=((size_t)((d3[i]+10.0)*(1.0+1.0e-12)));
    if (szt_tmp>20) szt_tmp=21;
    counts[szt_tmp]++;
  }
  cout << "Counts for merged file: ";
  vector_out(cout,counts,true);
  
  cout << "Writing merged results." << endl;
  write_results(out);
  cout << "Done." << endl;
  
  return 0;
}

int eos_nuclei::compare_tables(std::vector<std::string> &sv,
			       bool itive_com) {
  
  if (sv.size()<2) {
    cerr << "Command 'compare-tables' needs 2 arguments." << endl;
    cerr << "<file 1> <file 2> [quantity]" << endl;
    return 1;
  }

  string in1=sv[1];
  string in2=sv[2];
  string quantity;
  if (sv.size()>=4) {
    quantity=sv[3];
  }

  // Read the first file
  
  read_results(in1);

  // Read the second file
  
  eos_nuclei en2;
  en2.read_results(in2);

  bool grids_same=true;
  if (en2.n_nB2!=n_nB2 || en2.n_Ye2!=n_Ye2 || en2.n_T2!=n_T2) {
    cout << "Grids have different size." << endl;
    grids_same=false;
  } else {
    for(size_t i=0;i<n_nB2;i++) {
      if (fabs(en2.nB_grid2[i]-nB_grid2[i])/en2.nB_grid2[i]>1.0e-8) {
	grids_same=false;
      }
    }
    for(size_t i=0;i<n_Ye2;i++) {
      if (fabs(en2.Ye_grid2[i]-Ye_grid2[i])/en2.Ye_grid2[i]>1.0e-8) {
	grids_same=false;
      }
    }
    for(size_t i=0;i<n_T2;i++) {
      if (fabs(en2.T_grid2[i]-T_grid2[i])/en2.T_grid2[i]>1.0e-8) {
	grids_same=false;
      }
    }
    if (grids_same==false) {
      cout << "Grids have same size but different values." << endl;
    }
  }
  if (grids_same) {
    cout << "The grids match to within 1.0e-8." << endl;
  }

  vector<string> names;
  names.push_back("log_xn");
  names.push_back("log_xp");
  names.push_back("Z");
  names.push_back("A");
  names.push_back("Fint");
  names.push_back("Sint");
  names.push_back("Eint");
  names.push_back("Xn");
  names.push_back("Xp");
  names.push_back("Xalpha");
  names.push_back("Xnuclei");
  names.push_back("Xd");
  names.push_back("Xt");
  names.push_back("XHe3");
  names.push_back("XLi4");
  if (derivs_computed) {
    names.push_back("Pint");
    names.push_back("mun");
    names.push_back("mup");
    if (with_leptons) {
      names.push_back("F");
      names.push_back("E");
      names.push_back("P");
      names.push_back("S");
      names.push_back("mue");
    }
  }
  
  vector<o2scl::tensor_grid<> *> ptrs;
  ptrs.push_back(&tg_log_xn);
  ptrs.push_back(&tg_log_xp);
  ptrs.push_back(&tg_Z);
  ptrs.push_back(&tg_A);
  ptrs.push_back(&tg_Fint);
  ptrs.push_back(&tg_Xn);
  ptrs.push_back(&tg_Xp);
  ptrs.push_back(&tg_Xalpha);
  ptrs.push_back(&tg_Xnuclei);
  ptrs.push_back(&tg_Xd);
  ptrs.push_back(&tg_Xt);
  ptrs.push_back(&tg_XHe3);
  ptrs.push_back(&tg_XLi4);
  ptrs.push_back(&tg_Sint);
  ptrs.push_back(&tg_Eint);
  if (derivs_computed) {
    ptrs.push_back(&tg_Pint);
    ptrs.push_back(&tg_mun);
    ptrs.push_back(&tg_mup);
    if (with_leptons) {
      ptrs.push_back(&tg_F);
      ptrs.push_back(&tg_E);
      ptrs.push_back(&tg_P);
      ptrs.push_back(&tg_S);
      ptrs.push_back(&tg_mue);
    }
  }
  
  vector<o2scl::tensor_grid<> *> ptrs2;
  ptrs2.push_back(&en2.tg_log_xn);
  ptrs2.push_back(&en2.tg_log_xp);
  ptrs2.push_back(&en2.tg_Z);
  ptrs2.push_back(&en2.tg_A);
  ptrs2.push_back(&en2.tg_Fint);
  ptrs2.push_back(&en2.tg_Xn);
  ptrs2.push_back(&en2.tg_Xp);
  ptrs2.push_back(&en2.tg_Xalpha);
  ptrs2.push_back(&en2.tg_Xnuclei);
  ptrs2.push_back(&en2.tg_Xd);
  ptrs2.push_back(&en2.tg_Xt);
  ptrs2.push_back(&en2.tg_XHe3);
  ptrs2.push_back(&en2.tg_XLi4);
  ptrs2.push_back(&en2.tg_Sint);
  ptrs2.push_back(&en2.tg_Eint);
  if (derivs_computed) {
    ptrs2.push_back(&en2.tg_Pint);
    ptrs2.push_back(&en2.tg_mun);
    ptrs2.push_back(&en2.tg_mup);
    if (with_leptons) {
      ptrs2.push_back(&en2.tg_F);
      ptrs2.push_back(&en2.tg_E);
      ptrs2.push_back(&en2.tg_P);
      ptrs2.push_back(&en2.tg_S);
      ptrs2.push_back(&en2.tg_mue);
    }
  }
  
  if (grids_same) {
    
    // Visually separate from the file reading output
    cout << endl;
    
    bool found_points=true;
    for(size_t ell=0;ell<ptrs.size() && found_points;ell++) {
      if (quantity.length()==0 || names[ell]==quantity) {
	double max_dev=0.0;
	size_t imax=0, jmax=0, kmax=0;
	bool found=false;
	for(size_t i=0;i<n_nB2;i++) {
	  for(size_t j=0;j<n_Ye2;j++) {
	    for(size_t k=0;k<n_T2;k++) {
              vector<size_t> ix={i,j,k};
	      if (tg_flag.get(ix)>9.9 && en2.tg_flag.get(ix)>9.9) {
		found=true;
		double v1=ptrs[ell]->get(ix);
		double v2=ptrs2[ell]->get(ix);
		if (fabs(v1-v2)/fabs(v1)>max_dev) {
		  max_dev=fabs(v1-v2)/fabs(v1);
		  imax=i;
		  jmax=j;
		  kmax=k;
		}
	      }
	    }
	  }
	}
	if (found==true) {
          vector<size_t> ix_max={imax,jmax,kmax};
	  cout << "The maximum deviation for " << names[ell]
	       << " is at point (" << imax << "," << jmax << ","
	       << kmax << ")" << endl;
	  cout.precision(5);
	  cout << "  (" << en2.nB_grid2[imax] << ","
	       << en2.Ye_grid2[jmax] << "," << en2.T_grid2[kmax]
	       << ") "
	       << "with relative deviation " << max_dev << endl;
	  cout.precision(6);
	  cout << "Value in file " << in1 << " is "
	       << ptrs[ell]->get(ix_max) << endl;
	  cout << "Value in file " << in2 << " is "
	       << ptrs2[ell]->get(ix_max) << endl;
	  if (quantity.length()==0) cout << endl;
	} else {
	  cout << "Could not find any points to compare." << endl;
	  found_points=false;
	}
      }
    }
  } else {
    cout << "Compare tables doesn't yet support different grids." << endl;
  }
  
  return 0;
}

void eos_nuclei::new_table() {

  n_nB2=301;
  n_Ye2=70;
  n_T2=160;

  if (strange_axis) {
    n_S2=5;
  } else {
    n_S2=0;
  }
  cout << "Beginning new table." << endl;

  calc_utf8<> calc;
  std::map<std::string,double> vars;
    
  vector<double> packed;
  vector<std::string> split_res;
  
  split_string_delim(nB_grid_spec,split_res,',');
  n_nB2=stoszt(split_res[0]);
  
  calc.compile(split_res[1].c_str());
  for(size_t i=0;i<n_nB2;i++) {
    vars["i"]=((double)i);
    nB_grid2.push_back(calc.eval(&vars));
    packed.push_back(nB_grid2[i]);
  }
    
  split_string_delim(Ye_grid_spec,split_res,',');
  n_Ye2=stoszt(split_res[0]);

  calc.compile(split_res[1].c_str());
  for(size_t i=0;i<n_Ye2;i++) {
    vars["i"]=((double)i);
    Ye_grid2.push_back(calc.eval(&vars));
    packed.push_back(Ye_grid2[i]);
  }
    
  split_string_delim(T_grid_spec,split_res,',');
  n_T2=stoszt(split_res[0]);

  calc.compile(split_res[1].c_str());
  for(size_t i=0;i<n_T2;i++) {
    vars["i"]=((double)i);
    T_grid2.push_back(calc.eval(&vars));
    packed.push_back(T_grid2[i]);
  }

  if (strange_axis) {
    split_string_delim(S_grid_spec,split_res,',');
    n_S2=stoszt(split_res[0]);
    
    calc.compile(split_res[1].c_str());
    for(size_t i=0;i<n_S2;i++) {
      vars["i"]=((double)i);
      S_grid2.push_back(calc.eval(&vars));
      packed.push_back(S_grid2[i]);
    }
  }
  
  if (strange_axis) {
    
    size_t st[4]={n_nB2,n_Ye2,n_T2,n_S2};
    tg_log_xn.resize(4,st);
    tg_log_xp.resize(4,st);
    tg_Z.resize(4,st);
    tg_A.resize(4,st);
    tg_flag.resize(4,st);
    tg_Fint.resize(4,st);
    tg_Sint.resize(4,st);
    tg_Eint.resize(4,st);
  
    tg_Xn.resize(4,st);
    tg_Xp.resize(4,st);
    tg_Xalpha.resize(4,st);
    tg_Xnuclei.resize(4,st);
    tg_Xd.resize(4,st);
    tg_Xt.resize(4,st);
    tg_XHe3.resize(4,st);
    tg_XLi4.resize(4,st);

  } else {
    
    size_t st[3]={n_nB2,n_Ye2,n_T2};
    tg_log_xn.resize(3,st);
    tg_log_xp.resize(3,st);
    tg_Z.resize(3,st);
    tg_A.resize(3,st);
    tg_flag.resize(3,st);
    tg_Fint.resize(3,st);
    tg_Sint.resize(3,st);
    tg_Eint.resize(3,st);
  
    tg_Xn.resize(3,st);
    tg_Xp.resize(3,st);
    tg_Xalpha.resize(3,st);
    tg_Xnuclei.resize(3,st);
    tg_Xd.resize(3,st);
    tg_Xt.resize(3,st);
    tg_XHe3.resize(3,st);
    tg_XLi4.resize(3,st);

  }
  
  tg_log_xn.set_grid_packed(packed);
  tg_log_xp.set_grid_packed(packed);
  tg_Z.set_grid_packed(packed);
  tg_A.set_grid_packed(packed);
  tg_flag.set_grid_packed(packed);
  tg_Fint.set_grid_packed(packed);
  tg_Sint.set_grid_packed(packed);
  tg_Eint.set_grid_packed(packed);
  
  tg_Xn.set_grid_packed(packed);
  tg_Xp.set_grid_packed(packed);
  tg_Xalpha.set_grid_packed(packed);
  tg_Xnuclei.set_grid_packed(packed);
  tg_Xd.set_grid_packed(packed);
  tg_Xt.set_grid_packed(packed);
  tg_XHe3.set_grid_packed(packed);
  tg_XLi4.set_grid_packed(packed);

  tg_log_xn.set_all(0.0);
  tg_log_xp.set_all(0.0);
  tg_Z.set_all(0.0);
  tg_A.set_all(0.0);
  tg_flag.set_all(0.0);
  tg_Fint.set_all(0.0);
  tg_Sint.set_all(0.0);
  tg_Eint.set_all(0.0);

  tg_Xn.set_all(0.0);
  tg_Xp.set_all(0.0);
  tg_Xalpha.set_all(0.0);
  tg_Xnuclei.set_all(0.0);
  tg_Xd.set_all(0.0);
  tg_Xt.set_all(0.0);
  tg_XHe3.set_all(0.0);
  tg_XLi4.set_all(0.0);
  
  if (alg_mode==2 || alg_mode==3 || alg_mode==4) {

    if (strange_axis) {
    
      size_t st[4]={n_nB2,n_Ye2,n_T2,n_S2};
      tg_A_min.resize(4,st);
      tg_A_max.resize(4,st);
      tg_NmZ_min.resize(4,st);
      tg_NmZ_max.resize(4,st);

    } else {

      size_t st[3]={n_nB2,n_Ye2,n_T2};
      tg_A_min.resize(3,st);
      tg_A_max.resize(3,st);
      tg_NmZ_min.resize(3,st);
      tg_NmZ_max.resize(3,st);
    
    }
    
    tg_A_min.set_grid_packed(packed);
    tg_A_max.set_grid_packed(packed);
    tg_NmZ_min.set_grid_packed(packed);
    tg_NmZ_max.set_grid_packed(packed);
    
    tg_A_min.set_all(0.0);
    tg_A_max.set_all(0.0);
    tg_NmZ_min.set_all(0.0);
    tg_NmZ_max.set_all(0.0);

  }

  if ((derivs_computed || include_muons) && with_leptons) {
    
    if (strange_axis) {
    
      size_t st[4]={n_nB2,n_Ye2,n_T2,n_S2};
      tg_mue.resize(4,st);
      tg_mue.set_grid_packed(packed);
      tg_mue.set_all(0.0);
      
    } else {
      
      size_t st[3]={n_nB2,n_Ye2,n_T2};
      tg_mue.resize(3,st);
      tg_mue.set_grid_packed(packed);
      tg_mue.set_all(0.0);

    }
    
  }
  
  if (derivs_computed) {
    
    if (strange_axis) {
      
      size_t st[4]={n_nB2,n_Ye2,n_T2,n_S2};
      tg_Pint.resize(4,st);
      tg_Pint.set_grid_packed(packed);
      tg_Pint.set_all(0.0);
      tg_mun.resize(4,st);
      tg_mun.set_grid_packed(packed);
      tg_mun.set_all(0.0);
      tg_mup.resize(4,st);
      tg_mup.set_grid_packed(packed);
      tg_mup.set_all(0.0);

    } else {

      size_t st[3]={n_nB2,n_Ye2,n_T2};
      tg_Pint.resize(3,st);
      tg_Pint.set_grid_packed(packed);
      tg_Pint.set_all(0.0);
      tg_mun.resize(3,st);
      tg_mun.set_grid_packed(packed);
      tg_mun.set_all(0.0);
      tg_mup.resize(3,st);
      tg_mup.set_grid_packed(packed);
      tg_mup.set_all(0.0);

    }
    
    if (with_leptons) {
      
      if (strange_axis) {
      
        size_t st[4]={n_nB2,n_Ye2,n_T2,n_S2};
        tg_E.resize(4,st);
        tg_E.set_grid_packed(packed);
        tg_E.set_all(0.0);
        tg_P.resize(4,st);
        tg_P.set_grid_packed(packed);
        tg_P.set_all(0.0);
        tg_S.resize(4,st);
        tg_S.set_grid_packed(packed);
        tg_S.set_all(0.0);
        tg_F.resize(4,st);
        tg_F.set_grid_packed(packed);
        tg_F.set_all(0.0);
      
      } else {
      
        size_t st[3]={n_nB2,n_Ye2,n_T2};
        tg_E.resize(3,st);
        tg_E.set_grid_packed(packed);
        tg_E.set_all(0.0);
        tg_P.resize(3,st);
        tg_P.set_grid_packed(packed);
        tg_P.set_all(0.0);
        tg_S.resize(3,st);
        tg_S.set_grid_packed(packed);
        tg_S.set_all(0.0);
        tg_F.resize(3,st);
        tg_F.set_grid_packed(packed);
        tg_F.set_all(0.0);
      
      }
      
    }
  }

  loaded=true;
  
  return;
}  

int eos_nuclei::edit_data(std::vector<std::string> &sv,
			  bool itive_com) {

  string select_func=sv[1];
  string tensor_to_change, value_func, out_file;
  
  if (sv.size()>2) {
    if (sv.size()<4) {
      cerr << "Not enough arguments to edit-data." << endl;
      return 2;
    }
    tensor_to_change=sv[2];
    value_func=sv[3];
  } else if (sv.size()<2) {
    cerr << "Not enough arguments to edit-data." << endl;
    return 1;
  }
  
  size_t count=0;
  
  calc_utf8<> calc, calc2;
  std::map<std::string,double> vars;
  calc.compile(select_func.c_str());
  if (sv.size()>2) {
    calc2.compile(value_func.c_str());
  }

  for(int inB=0;inB<((int)n_nB2);inB++) {
    for(int iYe=0;iYe<((int)n_Ye2);iYe++) {
      for(int iT=0;iT<((int)n_T2);iT++) {

        vector<size_t> ix={((size_t)inB),((size_t)iYe),((size_t)iT)};
	
	vars["inB"]=inB;
	vars["iYe"]=iYe;
	vars["iT"]=iT;
	vars["nnB"]=nB_grid2.size();
	vars["nYe"]=Ye_grid2.size();
	vars["nT"]=T_grid2.size();
	vars["nB"]=nB_grid2[inB];
	vars["Ye"]=Ye_grid2[iYe];
	vars["T"]=T_grid2[iT];
	
	vars["flag"]=tg_flag.get(ix);
	vars["Fint"]=tg_Fint.get(ix);
	vars["Sint"]=tg_Sint.get(ix);
	vars["Eint"]=tg_Eint.get(ix);
	vars["Z"]=tg_Z.get(ix);
	vars["A"]=tg_A.get(ix);
	vars["log_xn"]=tg_log_xn.get(ix);
	vars["log_xp"]=tg_log_xp.get(ix);

        if (include_muons && with_leptons) {
          vars["Ymu"]=tg_Ymu.get(ix);
        }
	
	vars["Xn"]=tg_Xn.get(ix);
	vars["Xp"]=tg_Xp.get(ix);
	vars["Xalpha"]=tg_Xalpha.get(ix);
	vars["Xnuclei"]=tg_Xnuclei.get(ix);
	vars["Xd"]=tg_Xd.get(ix);
	vars["Xt"]=tg_Xt.get(ix);
	vars["XHe3"]=tg_XHe3.get(ix);
	vars["XLi4"]=tg_XLi4.get(ix);
	vars["dAdnB"]=0.0;
        
        vector<size_t> ixm1={((size_t)(inB-1)),((size_t)iYe),
          ((size_t)iT)};
        vector<size_t> ixp1={((size_t)(inB+1)),((size_t)iYe),
          ((size_t)iT)};
        
	if (inB>0 && inB<((int)n_nB2)-1) {
          vars["dPdnB"]=(tg_P.get(ixp1)-tg_P.get(ixm1))/
            (nB_grid2[inB+1]-nB_grid2[inB-1]);
        } else if (inB==0) {
          vars["dPdnB"]=(tg_P.get(ixp1)-tg_P.get(ix))/
            (nB_grid2[inB+1]-nB_grid2[inB]);
        } else {
          vars["dPdnB"]=(tg_P.get(ix)-tg_P.get(ixm1))/
            (nB_grid2[inB]-nB_grid2[inB-1]);
        }
 
	if (inB>0 && inB<((int)n_nB2)-1) {
          vars["dAdnB"]=(tg_A.get(ixp1)-tg_A.get(ixm1))/
            (nB_grid2[inB+1]-nB_grid2[inB-1]);
        } else if (inB==0) {
          vars["dAdnB"]=(tg_A.get(ixp1)-tg_A.get(ix))/
            (nB_grid2[inB+1]-nB_grid2[inB]);
        } else {
          vars["dAdnB"]=(tg_A.get(ix)-tg_A.get(ixm1))/
            (nB_grid2[inB]-nB_grid2[inB-1]);
        }
 
	if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
	  vars["A_min"]=tg_A_min.get(ix);
	  vars["A_max"]=tg_A_max.get(ix);
	  vars["NmZ_min"]=tg_NmZ_min.get(ix);
	  vars["NmZ_max"]=tg_NmZ_max.get(ix);
	}

	if (derivs_computed) {
	  vars["Pint"]=tg_Pint.get(ix);
	  vars["mun"]=tg_mun.get(ix);
	  vars["mup"]=tg_mup.get(ix);
          if (with_leptons) {
            vars["F"]=tg_F.get(ix);
            vars["E"]=tg_E.get(ix);
            vars["P"]=tg_P.get(ix);
            vars["S"]=tg_S.get(ix);
          }
	}

	double val=calc.eval(&vars);
	if (val>0.5) {

	  count++;

          if (count<10) {
            cout << count << "/? (nB,Ye,T[MeV]) "
                 << nB_grid2[inB] << " "
                 << Ye_grid2[iYe] << " "
                 << T_grid2[iT] << endl;
          }
          
	  if (sv.size()>3) {
	    double val2=calc2.eval(&vars);
	    if (tensor_to_change=="flag") {
	      tg_flag.get(ix)=val2;
	    } else if (tensor_to_change=="Fint") {
	      tg_Fint.get(ix)=val2;
	    } else if (tensor_to_change=="Z") {
	      tg_Z.get(ix)=val2;
	    } else if (tensor_to_change=="A") {
	      tg_A.get(ix)=val2;
	    } else if (tensor_to_change=="log_xn") {
	      tg_log_xn.get(ix)=val2;
	    } else if (tensor_to_change=="log_xp") {
	      tg_log_xp.get(ix)=val2;
	    } else if (tensor_to_change=="Xn") {
	      tg_Xn.get(ix)=val2;
	    } else if (tensor_to_change=="Xp") {
	      tg_Xp.get(ix)=val2;
	    } else if (tensor_to_change=="Xalpha") {
	      tg_Xalpha.get(ix)=val2;
	    } else if (tensor_to_change=="Xnuclei") {
	      tg_Xnuclei.get(ix)=val2;
	    } else if (tensor_to_change=="Xd") {
	      tg_Xd.get(ix)=val2;
	    } else if (tensor_to_change=="Xt") {
	      tg_Xt.get(ix)=val2;
	    } else if (tensor_to_change=="XHe3") {
	      tg_XHe3.get(ix)=val2;
	    } else if (tensor_to_change=="XLi4") {
	      tg_XLi4.get(ix)=val2;
	    } else if (tensor_to_change=="A_min") {
	      tg_A_min.get(ix)=val2;
	    } else if (tensor_to_change=="A_max") {
	      tg_A_max.get(ix)=val2;
	    } else if (tensor_to_change=="NmZ_min") {
	      tg_NmZ_min.get(ix)=val2;
	    } else if (tensor_to_change=="NmZ_max") {
	      tg_NmZ_max.get(ix)=val2;
	    } else if (tensor_to_change=="Ymu") {
	      tg_Ymu.get(ix)=val2;
	    } else if (tensor_to_change=="Eint") {
	      tg_Eint.get(ix)=val2;
	    } else if (tensor_to_change=="Pint") {
	      tg_Pint.get(ix)=val2;
	    } else if (tensor_to_change=="Sint") {
	      tg_Sint.get(ix)=val2;
	    } else if (tensor_to_change=="F") {
	      tg_F.get(ix)=val2;
	    } else if (tensor_to_change=="E") {
	      tg_E.get(ix)=val2;
	    } else if (tensor_to_change=="P") {
	      tg_P.get(ix)=val2;
	    } else if (tensor_to_change=="S") {
	      tg_S.get(ix)=val2;
	    } else if (tensor_to_change=="mue") {
	      tg_mue.get(ix)=val2;
	    } else {
	      cerr << "Invalid value for tensor to change." << endl;
	      return 3;
	    }
	  }
	}
      }
    }
  }
  
  cout << "Matched " << count << "/" << n_nB2*n_Ye2*n_T2
       << " entries." << endl;
  
  return 0;
}

int eos_nuclei::generate_table(std::vector<std::string> &sv,
			       bool itive_com) {

  if (derivs_computed) {
    cout << "Function generate-table setting derivs_computed to false."
	 << endl;
    derivs_computed=false;
  }
  if (with_leptons) {
    cout << "Function generate-table setting with_leptons to false."
	 << endl;
    with_leptons=false;
  }
  
  int mpi_rank=0, mpi_size=1;
  double mpi_start_time;
  
#ifndef NO_MPI

  // Compute rank and number of processors
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  mpi_start_time=MPI_Wtime();

  // MPI status object
  MPI_Status status;

#else

  mpi_start_time=time(0);
  
#endif

  static const size_t buf_size=100;
  int gt_verbose=2;

  // -----------------------------------------------------
  // External guess
  
#ifndef NO_MPI
  // Get MPI rank, etc.
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  
  // Ensure that multiple MPI ranks aren't reading from the
  // filesystem at the same time
  int tag=0, buffer2=0;
  if (mpi_size>1 && mpi_rank>=1) {
    MPI_Recv(&buffer2,1,MPI_INT,mpi_rank-1,
	     tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
#endif

  eos_nuclei external;
  if (ext_guess.length()>0) {
    external.read_results(ext_guess);
    if (nB_grid2!=external.nB_grid2 ||
	Ye_grid2!=external.Ye_grid2 ||
	T_grid2!=external.T_grid2) {
      O2SCL_ERR("Grids don't match for external guess.",
		o2scl::exc_einval);
    }
  }
  
#ifndef NO_MPI
  // Send a message to the next MPI rank
  if (mpi_size>1 && mpi_rank<mpi_size-1) {
    MPI_Send(&buffer2,1,MPI_INT,mpi_rank+1,
	     tag,MPI_COMM_WORLD);
  }
#endif
  
  // -----------------------------------------------------
  // All processors read input file in turn

  string out_file;
  if (sv.size()>=2) {
    out_file=sv[1];
  } else {
    cout << "Automatically setting output file to \"eos_nuclei.o2\"."
	 << endl;
    out_file="eos_nuclei.o2";
  }
  
  if (mpi_rank==0) {

    // -----------------------------------------------------
    // Parent MPI process or sole process for non-MPI version

    double last_file_time;
#ifdef NO_MPI
    last_file_time=time(0);
#else
    last_file_time=MPI_Wtime();
#endif

    // -----------------------------------------------------
    // Read the input file
    
    if (loaded==false) {

      cout << "No data. Creating new table." << endl;
      
      new_table();
      
      double nB=nB_grid2[250];
      double Ye=Ye_grid2[50];
      double T=T_grid2[80]/hc_mev_fm;
      double log_xn=-2.46;
      double log_xp=-1.64;
      if (alg_mode==2 || alg_mode==4) {
	log_xn=-2.5;
	log_xp=-1.64;
      }
      double Zbar, Nbar;
      size_t nuc_Z1=50, nuc_N1=50;
      int NmZ_min=-10;
      int NmZ_max=9;
      int A_min=8;
      int A_max=126;
      double f_min;
      thermo thx;
      double mun_full=0.0, mup_full=0.0;
      int first_ret=0;
      map<string,double> vdet;
      if (include_muons) {
        vdet["mue"]=electron.m;
      }
      if (alg_mode==1) {
	cerr << "Mode alg_mode=1 no longer supported in generate-table."
	     << endl;
	return 1;
      } else if (alg_mode==0) {
	first_ret=eos_vary_ZN(nB,Ye,T,log_xn,log_xp,nuc_Z1,nuc_N1,
			      thx,mun_full,mup_full,false);
	Zbar=nuc_Z1;
	Nbar=nuc_N1;
      } else if (alg_mode==2 || alg_mode==4) {
	first_ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,thx,
				mun_full,mup_full,A_min,A_max,
				NmZ_min,NmZ_max,vdet,true,false);
      } else if (alg_mode==3) {
	first_ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,thx,
				mun_full,mup_full,A_min,A_max,
				NmZ_min,NmZ_max,vdet,true,false);
      }
      if (first_ret!=0) {
	cerr << "Initial point for blank table failed." << endl;
	return 1;
      }
      
      ubvector X;
      compute_X(nB,X);
      
      store_point(250,50,80,nB,Ye,T,thx,log_xn,log_xp,
		  Zbar,Nbar,mun_full,mup_full,X,A_min,A_max,
		  NmZ_min,NmZ_max,10.0,vdet);

      // End of 'else' for 'if (in_file!="none")'
    }

    // Adjust flags in input file if necessary
    for(int inB=0;inB<((int)n_nB2);inB++) {
      for(int iYe=0;iYe<((int)n_Ye2);iYe++) {
	for(int iT=0;iT<((int)n_T2);iT++) {
          vector<size_t> ix={((size_t)inB),((size_t)iYe),((size_t)iT)};
	  int iflag=((int)(tg_flag.get(ix)*(1.0+1.0e-12)));
	  // If recompute is true, set all finished points to "guess"
	  if (recompute==true) {
	    if (iflag==iflag_done) {
	      tg_flag.get(ix)=iflag_guess;
	    }
	  }
	  // Get rid of previous in progress markings
	  if (iflag==iflag_in_progress_empty) {
	    tg_flag.get(ix)=iflag_empty;
	  }
	  if (iflag==iflag_in_progress_with_guess) {
	    tg_flag.get(ix)=iflag_guess;
	  }
	}
      }
    }

    if (edge_list.length()>0) {
      vector<string> edge_vec;
      split_string(edge_list,edge_vec);
      for(size_t k=0;k<edge_vec.size();k++) {
	o2scl::tensor_grid<> *ptr;
	if (edge_vec[k]=="A") {
	  ptr=&tg_A;
	} else if (edge_vec[k]=="Xnuclei") {
	  ptr=&tg_Xnuclei;
	} else if (edge_vec[k]=="Xalpha") {
	  ptr=&tg_Xalpha;
	} else {
	  cerr << "Invalid name in edge_list." << endl;
	  return 5;
	}
	for(int inB=1;inB<((int)n_nB2-1);inB++) {
	  for(int iYe=1;iYe<((int)n_Ye2-1);iYe++) {
	    for(int iT=1;iT<((int)n_T2-1);iT++) {
              vector<size_t> ix={((size_t)inB),((size_t)iYe),((size_t)iT)};
              vector<size_t> ixm1={((size_t)(inB-1)),((size_t)iYe),
                ((size_t)iT)};
              vector<size_t> ixp1={((size_t)(inB+1)),((size_t)iYe),
                ((size_t)iT)};
              vector<size_t> ixm2={((size_t)inB),((size_t)(iYe-1)),
                ((size_t)iT)};
              vector<size_t> ixp2={((size_t)inB),((size_t)(iYe+1)),
                ((size_t)iT)};
              vector<size_t> ixm3={((size_t)inB),((size_t)iYe),
                ((size_t)(iT-1))};
              vector<size_t> ixp3={((size_t)inB),((size_t)iYe),
                ((size_t)(iT+1))};
              
	      int iflag=((int)(tg_flag.get(ix)*(1.0+1.0e-12)));
	      if (iflag==iflag_done) {
		double X0=ptr->get(ix);

		if (X0>1.0e-4) {

		  if (tg_flag.get(ixm1)>9.9 &&
		      tg_flag.get(ixp1)>9.9) {
		    double X1=ptr->get(ixm1);
		    double X2=ptr->get(ixp1);
		
		    if (X0<std::min(X1,X2) || X0>std::max(X1,X2)) {
		      tg_flag.get(ixm1)=5.0;
		      tg_flag.get(ixp1)=5.0;
		      tg_flag.get(ix)=5.0;
		    }
		  }
		  if (false) {
		    if (tg_flag.get(ixm2)>9.9 &&
			tg_flag.get(ixp2)>9.9) {
		      double X1=ptr->get(ixm2);
		      double X2=ptr->get(ixp2);
		    
		      if (X0<std::min(X1,X2) || X0>std::max(X1,X2)) {
			tg_flag.get(ix)=5.0;
			tg_flag.get(ixm2)=5.0;
			tg_flag.get(ixp2)=5.0;
		      }
		    }
		  }
		  if (tg_flag.get(ixm3)>9.9 &&
		      tg_flag.get(ixp3)>9.9) {
		    double X1=ptr->get(ixm3);
		    double X2=ptr->get(ixp3);
		  
		    if (X0<std::min(X1,X2) || X0>std::max(X1,X2)) {
		      tg_flag.get(ix)=5.0;
		      tg_flag.get(ixm3)=5.0;
		      tg_flag.get(ixp3)=5.0;
		    }
		  }
		
		}
	      }
	    }
	  }
	}
	
      }
      
      write_results("temp.o2");

      // End of 'if (edge_list.length()>0)'
    }

    // Task counters
    size_t total_tasks=n_nB2*n_Ye2*n_T2;
    size_t current_tasks=0;

    // -----------------------------------------------------
    // Main calculation loop

    bool done=false;
    while (done==false) {

      bool one_success=false;
      
      // -----------------------------------------------------
      // Compute tasks (either w/o MPI or w/MPI on rank 0)
    
      // Setup task list as a two sets of triplets, source first,
      // destination second
      vector<size_t> tasks;

      // Store initial guesses in a table, one row for each task
      table<> gtab;
      gtab.line_of_names(((string)"log_xn log_xp Z A A_min ")+
                         "A_max NmZ_min NmZ_max mue");

      size_t nB_step=0;
      size_t Ye_step=0;
      size_t T_step=0;

      // Interpret variable Ye_list as an array of type size_t
      vector<size_t> Ye_list_sizet;
      if (Ye_list.length()==0) {
	for(size_t i=0;i<n_Ye2;i++) {
	  Ye_list_sizet.push_back(i);
	}
      } else {
	int lret=string_to_uint_list(Ye_list,Ye_list_sizet);
	if (lret!=0) {
	  cerr << "Could not interpret Ye list: " << Ye_list << endl;
	  return 4;
	}
      }

      // Output Ye_list
      string sout=((string)"Ye_list: ")+Ye_list+" ";
      for(size_t jk=0;jk<Ye_list_sizet.size();jk++) {
	sout+=o2scl::szttos(Ye_list_sizet[jk])+" ";
      }
      vector<string> sv2;
      rewrap(sout,sv2);
      for(size_t jk=0;jk<sv2.size();jk++) {
	cout << sv2[jk] << endl;
      }

      // Loop over all points to compute task list
      for(int inB=0;inB<((int)n_nB2);inB++) {
	for(int iYe_list=0;iYe_list<((int)Ye_list_sizet.size());
	    iYe_list++) {
	  int iYe=Ye_list_sizet[iYe_list];
	  for(int iT=0;iT<((int)n_T2);iT++) {
            vector<size_t> ix={((size_t)inB),((size_t)iYe),((size_t)iT)};
            vector<size_t> ixm1={((size_t)(inB-1)),((size_t)iYe),
              ((size_t)iT)};
            vector<size_t> ixp1={((size_t)(inB+1)),((size_t)iYe),
              ((size_t)iT)};
            vector<size_t> ixm2={((size_t)inB),((size_t)(iYe-1)),
              ((size_t)iT)};
            vector<size_t> ixp2={((size_t)inB),((size_t)(iYe+1)),
              ((size_t)iT)};
            vector<size_t> ixm3={((size_t)inB),((size_t)iYe),
              ((size_t)(iT-1))};
            vector<size_t> ixp3={((size_t)inB),((size_t)iYe),
              ((size_t)(iT+1))};
          
	    int iflag=((int)(tg_flag.get(ix)*(1.0+1.0e-12)));
	    bool guess_found=false;

	    if (iflag==iflag_guess) {
	      
	      // A point which is not finished, is not yet being
	      // computed, but has an initial guess
	      tasks.push_back(inB);
	      tasks.push_back(iYe);
	      tasks.push_back(iT);
	      tasks.push_back(inB);
	      tasks.push_back(iYe);
	      tasks.push_back(iT);

	      if (alg_mode>=2) {
                double mue=0.0;
                if (include_muons) {
                  mue=tg_mue.get(ix)/hc_mev_fm;
                }
		vector<double> line={tg_log_xn.get(ix),
                  tg_log_xp.get(ix),0.0,0.0,
                  tg_A_min.get(ix),tg_A_max.get(ix),tg_NmZ_min.get(ix),
                  tg_NmZ_max.get(ix),mue};
		gtab.line_of_data(line.size(),line);
	      } else {
		vector<double> line={tg_log_xn.get(ix),
                  tg_log_xp.get(ix),tg_Z.get(ix),tg_A.get(ix),
                  0.0,0.0,0.0,0.0,0.0};
		gtab.line_of_data(line.size(),line);
	      }	
	      
	      tg_flag.get(ix)=(double)iflag_in_progress_with_guess;
	      guess_found=true;
	    }

	    if (ext_guess.length()>0 &&
		external.tg_flag.get(ix)>9.9 &&
		(iflag==iflag_guess || iflag==iflag_empty)) {
	      tasks.push_back(inB);
	      tasks.push_back(iYe);
	      tasks.push_back(iT);
	      tasks.push_back(inB);
	      tasks.push_back(iYe);
	      tasks.push_back(iT);
	      if (alg_mode>=2) {
                double mue=0.0;
                if (include_muons) {
                  mue=tg_mue.get(ix)/hc_mev_fm;
                }
		vector<double> line={external.tg_log_xn.get(ix),
                  external.tg_log_xp.get(ix),
                  0.0,0.0,
                  external.tg_A_min.get(ix),
                  external.tg_A_max.get(ix),
                  external.tg_NmZ_min.get(ix),
                  external.tg_NmZ_max.get(ix),mue};
		gtab.line_of_data(line.size(),line);
	      } else {
		vector<double> line={external.tg_log_xn.get(ix),
                  external.tg_log_xp.get(ix),
                  external.tg_Z.get(ix),
                  external.tg_A.get(ix),
                  0.0,0.0,0.0,0.0,0.0};
		gtab.line_of_data(line.size(),line);
	      }
	      guess_found=true;
	    }
	    
	    // If six_neighbors is true, set up six additional tasks
	    // for neighbors with useful initial guesses
	    if (six_neighbors>0 &&
		(iflag==iflag_guess ||
		 (propagate_points && iflag!=iflag_done))) {

	      if (inB>0) {
		int iflag2=((int)(tg_flag.get(ixm1)*
				  (1.0+1.0e-12)));
		if (iflag2==iflag_in_progress_with_guess ||
		    iflag2==iflag_guess || iflag2==iflag_done) {
		  tasks.push_back(inB-1);
		  tasks.push_back(iYe);
		  tasks.push_back(iT);
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT);
		  
		  if (alg_mode>=2) {
                    double mue=0.0;
		    // If we're at small densities and we've converged
		    // at both a higher density and a lower density,
		    // then we can linearly interpolate to get a good
		    // guess. At higher densities this is not
		    // necessarily a good guess because of the
		    // liquid-gas phase transition
		    if (nB_grid2[inB]<1.0e-3 && inB<((int)n_nB2)-1 &&
			tg_flag.get(ixp1)>9.9) {
                      if (include_muons) {
                        mue=(tg_mue.get(ixm1)+
                             tg_mue.get(ixp1))/2.0/hc_mev_fm;
                      }
		      vector<double> line={(tg_log_xn.get(ixm1)+
                                            tg_log_xn.get(ixp1))/2.0,
                        (tg_log_xp.get(ixm1)+
                         tg_log_xp.get(ixp1))/2.0,
                        0.0,0.0,
                        (tg_A_min.get(ixm1)+
                         tg_A_min.get(ixp1))/2.0,
                        (tg_A_max.get(ixm1)+
                         tg_A_max.get(ixp1))/2.0,
                        (tg_NmZ_min.get(ixm1)+
                         tg_NmZ_min.get(ixp1))/2.0,
                        (tg_NmZ_max.get(ixm1)+
                         tg_NmZ_max.get(ixp1))/2.0,mue};
		      gtab.line_of_data(line.size(),line);
		    } else {
                      if (include_muons) {
                        mue=tg_mue.get(ixm1)/hc_mev_fm;
                      }
		      vector<double> line={tg_log_xn.get(ixm1),
                        tg_log_xp.get(ixm1),
                        0.0,0.0,
                        tg_A_min.get(ixm1),
                        tg_A_max.get(ixm1),
                        tg_NmZ_min.get(ixm1),
                        tg_NmZ_max.get(ixm1),mue};
		      gtab.line_of_data(line.size(),line);
		    }
		  } else {
                    double mue=0.0;
                    if (include_muons) {
                      mue=tg_mue.get(ixm1)/hc_mev_fm;
                    }
		    vector<double> line={tg_log_xn.get(ixm1),
                      tg_log_xp.get(ixm1),
                      tg_Z.get(ixm1),
                      tg_A.get(ixm1),
                      0.0,0.0,0.0,0.0,mue};
		    gtab.line_of_data(line.size(),line);
		  }	
		  
		  guess_found=true;
		}
	      }
	      
	      if (six_neighbors>4 && iYe>0) {
		int iflag2=((int)(tg_flag.get(ixm2)*
				  (1.0+1.0e-12)));
		if (iflag2==iflag_in_progress_with_guess ||
		    iflag2==iflag_guess || iflag2==iflag_done) {
		  tasks.push_back(inB);
		  tasks.push_back(iYe-1);
		  tasks.push_back(iT);
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT);
		  
		  if (alg_mode>=2) {
                    double mue=0.0;
                    if (include_muons) {
                      mue=tg_mue.get(ixm2)/hc_mev_fm;
                    }
		    vector<double> line={tg_log_xn.get(ixm2),
                      tg_log_xp.get(ixm2),
                      0.0,0.0,
                      tg_A_min.get(ixm2),
                      tg_A_max.get(ixm2),
                      tg_NmZ_min.get(ixm2),
                      tg_NmZ_max.get(ixm2),mue};
		    gtab.line_of_data(line.size(),line);
		  } else {
                    double mue=0.0;
                    if (include_muons) {
                      mue=tg_mue.get(ixm2)/hc_mev_fm;
                    }
		    vector<double> line={tg_log_xn.get(ixm2),
                      tg_log_xp.get(ixm2),
                      tg_Z.get(ixm2),
                      tg_A.get(ixm2),
                      0.0,0.0,0.0,0.0,mue};
		    gtab.line_of_data(line.size(),line);
		  }	
		  
		  guess_found=true;
		}
	      }
	      if (six_neighbors>2 && iT>0) {
		int iflag2=((int)(tg_flag.get(ixm3)*
				  (1.0+1.0e-12)));
		if (iflag2==iflag_in_progress_with_guess ||
		    iflag2==iflag_guess || iflag2==iflag_done) {
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT-1);
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT);

                  double mue=0.0;
                  if (include_muons) {
                    mue=tg_mue.get(ixm3)/hc_mev_fm;
                  }
                  
		  if (alg_mode>=2) {
		    vector<double> line={tg_log_xn.get(ixm3),
                      tg_log_xp.get(ixm3),
                      0.0,0.0,
                      tg_A_min.get(ixm3),
                      tg_A_max.get(ixm3),
                      tg_NmZ_min.get(ixm3),
                      tg_NmZ_max.get(ixm3),mue};
		    gtab.line_of_data(line.size(),line);
		  } else {
		    vector<double> line={tg_log_xn.get(ixm3),
                      tg_log_xp.get(ixm3),
                      tg_Z.get(ixm3),
                      tg_A.get(ixm3),
                      0.0,0.0,0.0,0.0,mue};
		    gtab.line_of_data(line.size(),line);
		  }	
		  
		  guess_found=true;
		}
	      }
	      if (six_neighbors>1 && inB<((int)n_nB2)-1) {
		int iflag2=((int)(tg_flag.get(ixp1)*
				  (1.0+1.0e-12)));
		if (iflag2==iflag_in_progress_with_guess ||
		    iflag2==iflag_guess || iflag2==iflag_done) {
		  tasks.push_back(inB+1);
		  tasks.push_back(iYe);
		  tasks.push_back(iT);
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT);

		  if (alg_mode>=2) {
		    if (nB_grid2[inB]<1.0e-3 && inB>0 &&
			tg_flag.get(ixm1)>9.9) {
                      double mue=0.0;
                      if (include_muons) {
                        mue=(tg_mue.get(ixm1)+
                             tg_mue.get(ixp1))/2.0/hc_mev_fm;
                      }
		      vector<double> line={(tg_log_xn.get(ixm1)+
                                            tg_log_xn.get(ixp1))/2.0,
                        (tg_log_xp.get(ixm1)+
                         tg_log_xp.get(ixp1))/2.0,
                        0.0,0.0,
                        (tg_A_min.get(ixm1)+
                         tg_A_min.get(ixp1))/2.0,
                        (tg_A_max.get(ixm1)+
                         tg_A_max.get(ixp1))/2.0,
                        (tg_NmZ_min.get(ixm1)+
                         tg_NmZ_min.get(ixp1))/2.0,
                        (tg_NmZ_max.get(ixm1)+
                         tg_NmZ_max.get(ixp1))/2.0,mue};
		      gtab.line_of_data(line.size(),line);
		    } else {
                      double mue=0.0;
                      if (include_muons) {
                        mue=tg_mue.get(ixp1)/hc_mev_fm;
                      }
		      vector<double> line={tg_log_xn.get(ixp1),
                        tg_log_xp.get(ixp1),
                        0.0,0.0,
                        tg_A_min.get(ixp1),
                        tg_A_max.get(ixp1),
                        tg_NmZ_min.get(ixp1),
                        tg_NmZ_max.get(ixp1),mue};
		      gtab.line_of_data(line.size(),line);
		    }
		  } else {
                    double mue=0.0;
                    if (include_muons) {
                      mue=tg_mue.get(ixp1)/hc_mev_fm;
                    }
		    vector<double> line={tg_log_xn.get(ixp1),
                      tg_log_xp.get(ixp1),
                      tg_Z.get(ixp1),
                      tg_A.get(ixp1),
                      0.0,0.0,0.0,0.0,mue};
		    gtab.line_of_data(line.size(),line);
		  }	
		  
		  guess_found=true;
		}
	      }
	      if (six_neighbors>4 && iYe<((int)n_Ye2)-1) {
		int iflag2=((int)(tg_flag.get(ixp2)*
				  (1.0+1.0e-12)));
		if (iflag2==iflag_in_progress_with_guess ||
		    iflag2==iflag_guess || iflag2==iflag_done) {
		  tasks.push_back(inB);
		  tasks.push_back(iYe+1);
		  tasks.push_back(iT);
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT);

                  double mue=0.0;
                  if (include_muons) {
                    mue=tg_mue.get(ixp2)/hc_mev_fm;
                  }
                  
		  if (alg_mode>=2) {
		    vector<double> line={tg_log_xn.get(ixp2),
                      tg_log_xp.get(ixp2),
                      0.0,0.0,
                      tg_A_min.get(ixp2),
                      tg_A_max.get(ixp2),
                      tg_NmZ_min.get(ixp2),
                      tg_NmZ_max.get(ixp2),mue};
		    gtab.line_of_data(line.size(),line);
		  } else {
		    vector<double> line={tg_log_xn.get(ixp2),
                      tg_log_xp.get(ixp2),
                      tg_Z.get(ixp2),
                      tg_A.get(ixp2),
                      0.0,0.0,0.0,0.0,mue};
		    gtab.line_of_data(line.size(),line);
		  }	
		  
		  guess_found=true;
		}
	      }
	      if (six_neighbors>2 && iT<((int)n_T2)-1) {
		int iflag2=((int)(tg_flag.get(ixp3)*
				  (1.0+1.0e-12)));
		if (iflag2==iflag_in_progress_with_guess ||
		    iflag2==iflag_guess || iflag2==iflag_done) {
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT+1);
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT);

                  double mue=0.0;
                  if (include_muons) {
                    mue=tg_mue.get(ixp3)/hc_mev_fm;
                  }
                  
		  if (alg_mode>=2) {
		    vector<double> line={tg_log_xn.get(ixp3),
                      tg_log_xp.get(ixp3),
                      0.0,0.0,
                      tg_A_min.get(ixp3),
                      tg_A_max.get(ixp3),
                      tg_NmZ_min.get(ixp3),
                      tg_NmZ_max.get(ixp3),mue};
		    gtab.line_of_data(line.size(),line);
		  } else {
		    vector<double> line={tg_log_xn.get(ixp3),
                      tg_log_xp.get(ixp3),
                      tg_Z.get(ixp3),
                      tg_A.get(ixp3),
                      0.0,0.0,0.0,0.0,mue};
		    gtab.line_of_data(line.size(),line);
		  }	
		  
		  guess_found=true;
		}
	      }
	    }
	    
	    if (propagate_points==true && guess_found==false && 
		iflag!=iflag_done &&
		iflag!=iflag_in_progress_empty &&
		iflag!=iflag_in_progress_with_guess) {
	      
	      // A point which is not finished, is not yet being
	      // computed, and has no initial guess, so we look
	      // for a guess nearby
	      
	      bool guess_found2=false;
	      for(int jnB=inB-nB_step;guess_found2==false && 
		    jnB<=inB+((int)nB_step);jnB++) {
		for(int jYe=iYe-Ye_step;guess_found2==false && 
		      jYe<=iYe+((int)Ye_step);jYe++) {
		  for(int jT=iT-T_step;guess_found2==false && 
			jT<=iT+((int)T_step);jT++) {
                    vector<size_t> jx={((size_t)jnB),
                      ((size_t)jYe),((size_t)jT)};
                    
		    // Ensure that the j indices are in bounds
		    if ((jnB!=inB || jYe!=iYe || jT!=iT) && jnB>=0 &&
			jYe>=0 && jT>=0 && jnB<((int)n_nB2) &&
			jYe<((int)n_Ye2) && jT<((int)n_T2)) {

		      int jflag=((int)(tg_flag.get(jx)*
				       (1.0+1.0e-12)));
		      
		      if (jflag==iflag_guess || jflag==iflag_done) {
			guess_found2=true;
			// Tasks are stored with j first and i
			// second, and the calculation goes j -> i
			tasks.push_back(jnB);
			tasks.push_back(jYe);
			tasks.push_back(jT);
			tasks.push_back(inB);
			tasks.push_back(iYe);
			tasks.push_back(iT);

                        double mue=0.0;
                        if (include_muons) {
                          mue=tg_mue.get(jx)/hc_mev_fm;
                        }
                        
			if (alg_mode>=2) {
			  vector<double> line={tg_log_xn.get(jx),
                            tg_log_xp.get(jx),
                            0.0,0.0,
                            tg_A_min.get(jx),
                            tg_A_max.get(jx),
                            tg_NmZ_min.get(jx),
                            tg_NmZ_max.get(jx),mue};
			  gtab.line_of_data(line.size(),line);
			} else {
			  vector<double> line={tg_log_xn.get(jx),
                            tg_log_xp.get(jx),
                            tg_Z.get(jx),
                            tg_A.get(jx),
                            0.0,0.0,0.0,0.0,mue};
			  gtab.line_of_data(line.size(),line);
			}	
		  
			if (iflag==iflag_guess) {
			  tg_flag.get(ix)=
			    (double)iflag_in_progress_with_guess;
			} else {
			  tg_flag.get(ix)=
			    (double)iflag_in_progress_empty;
			}
		      }
		    }
		    
		    // End of three loops over j indices
		  }
		}
	      }

	      // End over conditional for propagate points
	    }
	      
	    // End of three loops over i indices
	  }
	}
      }

      size_t ntasks=tasks.size()/6;
      if (gt_verbose>1) {
	if (ntasks>0) {
	  cout << "Rank " << mpi_rank 
	       << " tasks " << ntasks << endl;
	} else {
	  cout << "Rank " << mpi_rank << " tasks " << ntasks << endl;
	}
      }
      if (ntasks==0) {
	cout << "Found no tasks to complete." << endl;
      }

#ifndef NO_MPI

      // -------------------------------------------------------------
      // If we're running MPI with more than one rank, then send tasks
      // to the child MPI ranks until we run out of tasks.
      
      if (mpi_size>1) {
	
	// Setup buffers and MPI objects for each processor
	vector<MPI_Request> request(mpi_size-1);
	MPI_Request blank_request;
	vector<vector<double> > output_buffers(mpi_size-1);
	vector<vector<double> > input_buffers(mpi_size-1);
	
	for(size_t i=0;i<((size_t)(mpi_size-1));i++) {
	  output_buffers[i].resize(buf_size);
	  input_buffers[i].resize(buf_size);
	}
	
	// -----------------------------------------------------
	// Task loop
	
	for(size_t i=0;i<ntasks+((size_t)(mpi_size-1));i++) {
	  int proc_index; 
	  
	  if (i<((size_t)(mpi_size-1))) {
	    
	    // Automatically assign processors for the first set of
	    // tasks
	    proc_index=((int)i);
	    
	  } else {
	    
	    // Read the message from a processor which is finished
	    MPI_Waitany(mpi_size-1,&request[0],&proc_index,&status);
	    
	    int iflag=((int)(input_buffers[proc_index][vi["flag"]]*
			     (1.0+1.0e-12)));
	    size_t task_index=
	      ((size_t)(input_buffers[proc_index][vi["index"]]+1.0e-6));
	    
	    // Index of the point which was computed
	    size_t i2=task_index*6+3;
	    size_t j2=task_index*6+4;
	    size_t k2=task_index*6+5;
	    
	    // If it succeeded, store the results, otherwise, just
	    // leave the point alone
	    if (iflag>=10) {
              one_success=true;
	      thermo thy;
	      thy.ed=input_buffers[proc_index][vi["ed"]];
	      thy.pr=input_buffers[proc_index][vi["pr"]];
	      thy.en=input_buffers[proc_index][vi["en"]];
	      if (false && gt_verbose>1) {
		double nBt=input_buffers[proc_index][vi["nB"]];
		double Tt=input_buffers[proc_index][vi["T"]]/hc_mev_fm;
		cout << mpi_rank << " storing point "
		     << tasks[i2] << " " << tasks[j2] << " "
		     << tasks[k2] << " with fr=" << thy.ed-Tt*thy.en
		     << " F=" << (thy.ed-Tt*thy.en)/nBt*hc_mev_fm << endl;
	      }
	      ubvector X(6);
	      map<string,double> vdet;
	      X[0]=input_buffers[proc_index][vi["Xalpha"]];
	      X[1]=input_buffers[proc_index][vi["Xd"]];
	      X[2]=input_buffers[proc_index][vi["Xt"]];
	      X[3]=input_buffers[proc_index][vi["XHe3"]];
	      X[4]=input_buffers[proc_index][vi["XLi4"]];
	      X[5]=input_buffers[proc_index][vi["Xnuclei"]];
	      if (include_detail) {
		vdet["zn"]=input_buffers[proc_index][vi["zn"]];
		vdet["zp"]=input_buffers[proc_index][vi["zp"]];
	      }
              if (include_muons) {
		vdet["mue"]=input_buffers[proc_index][vi["mue"]];
		vdet["Ymu"]=input_buffers[proc_index][vi["Ymu"]];
              }
	      store_point(tasks[i2],tasks[j2],tasks[k2],
			  input_buffers[proc_index][vi["nB"]],
			  input_buffers[proc_index][vi["Ye"]],
			  input_buffers[proc_index][vi["T"]]/hc_mev_fm,
			  thy,
			  input_buffers[proc_index][vi["log_xn"]],
			  input_buffers[proc_index][vi["log_xp"]],
			  input_buffers[proc_index][vi["Z"]],
			  input_buffers[proc_index][vi["N"]],
			  input_buffers[proc_index][vi["mun"]],
			  input_buffers[proc_index][vi["mup"]],X,
			  input_buffers[proc_index][vi["A_min"]],
			  input_buffers[proc_index][vi["A_max"]],
			  input_buffers[proc_index][vi["NmZ_min"]],
			  input_buffers[proc_index][vi["NmZ_max"]],
			  input_buffers[proc_index][vi["flag"]],vdet);
	    } else {
	      if (gt_verbose>1) {
		cout << "Rank " << mpi_rank
		     << " point at (" << tasks[i2] << "," << tasks[j2]
		     << "," << tasks[k2] << ") (nB,Ye,T)=(";
		cout.precision(4);
		cout << input_buffers[proc_index][vi["nB"]] << ","
		     << input_buffers[proc_index][vi["Ye"]] << ","
		     << input_buffers[proc_index][vi["T"]] 
		     << ") failed." << endl;
		cout.precision(6);
	      }
	    }
	    
	  }
	  
	  if (i<ntasks) {

	    // Send the ith task to a processor to be computed
	    
	    // Index of the point which contains the initial guess
	    size_t i2=i*6, j2=i*6+1, k2=i*6+2;
	    
	    if (gt_verbose>1) {
	      cout << "Rank " << mpi_rank << " sending " 
		   << i << "/" << ntasks << " (" << tasks[i*6] << ","
		   << tasks[i*6+1] << "," 
		   << tasks[i*6+2] << ") -> ("
		   << tasks[i*6+3] << "," 
		   << tasks[i*6+4] << ","
		   << tasks[i*6+5] << ") to " << proc_index+1 << endl;
	    }

	    output_buffers[proc_index][vi["msg"]]=((double)message_continue);
	    
	    output_buffers[proc_index][vi["index"]]=i;
	    output_buffers[proc_index][vi["inB"]]=tasks[i*6+3];
	    output_buffers[proc_index][vi["iYe"]]=tasks[i*6+4];
	    output_buffers[proc_index][vi["iT"]]=tasks[i*6+5];
	    output_buffers[proc_index][vi["nB"]]=nB_grid2[tasks[i*6+3]];
	    output_buffers[proc_index][vi["Ye"]]=Ye_grid2[tasks[i*6+4]];
	    output_buffers[proc_index][vi["T"]]=T_grid2[tasks[i*6+5]];
	    output_buffers[proc_index][vi["log_xn"]]=gtab.get("log_xn",i);
	    output_buffers[proc_index][vi["log_xp"]]=gtab.get("log_xp",i);
	    output_buffers[proc_index][vi["Z"]]=gtab.get("Z",i);
	    output_buffers[proc_index][vi["N"]]=gtab.get("A",i)-
	      gtab.get("Z",i);
	    if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
	      output_buffers[proc_index][vi["A_min"]]=
		gtab.get("A_min",i);
	      output_buffers[proc_index][vi["A_max"]]=
		gtab.get("A_max",i);
	      output_buffers[proc_index][vi["NmZ_min"]]=
		gtab.get("NmZ_min",i);
	      output_buffers[proc_index][vi["NmZ_max"]]=
		gtab.get("NmZ_max",i);
	    }
            if (include_muons) {
	      output_buffers[proc_index][vi["mue"]]=
		gtab.get("mue",i);
	      output_buffers[proc_index][vi["Ymu"]]=
		gtab.get("Ymu",i);
            }

	    // Determine if the "no_nuclei" flag should be set
	    output_buffers[proc_index][vi["no_nuclei"]]=0.0;
            // AWS, 1/16/21: I'm taking this out as it is confusing
            // and can result in incorrect results if not used
            // carefully
	    if (false) {
	      size_t inB_dest=tasks[i*6+3];
	      size_t iYe_dest=tasks[i*6+4];
	      size_t iT_dest=tasks[i*6+5];
              vector<size_t> ix_dest_m3={inB_dest,iYe_dest,iT_dest-1};
	      if (iT_dest>0) {
		if (tg_flag.get(ix_dest_m3)>9.9) {
		  double X_all_nuclei=0.0;
		  X_all_nuclei+=tg_Xalpha.get(ix_dest_m3);
		  X_all_nuclei+=tg_Xnuclei.get(ix_dest_m3);
		  X_all_nuclei+=tg_Xd.get(ix_dest_m3);
		  X_all_nuclei+=tg_Xt.get(ix_dest_m3);
		  X_all_nuclei+=tg_XHe3.get(ix_dest_m3);
		  X_all_nuclei+=tg_XLi4.get(ix_dest_m3);
		  if (X_all_nuclei<1.0e-20 && nB_grid2[tasks[i*6+3]]<0.02) {
		    cout << "X_all_nuclei small (1): " << X_all_nuclei << " "
			 << T_grid2[iT_dest] << " "
                         << nB_grid2[tasks[i*6+3]] << endl;
		    output_buffers[proc_index][vi["no_nuclei"]]=1.0;
		  }
		}
	      }
	    }
	    
	    MPI_Isend(&(output_buffers[proc_index][0]),buf_size,
		      MPI_DOUBLE,proc_index+1,0,
		      MPI_COMM_WORLD,&blank_request);
	    
	    MPI_Irecv(&(input_buffers[proc_index][0]),buf_size,
		      MPI_DOUBLE,proc_index+1,1,
		      MPI_COMM_WORLD,&(request[proc_index]));
	  }

	  // Update file if necessary
	  if (((int)i)%file_update_iters==file_update_iters-1 ||
	      MPI_Wtime()-last_file_time>file_update_time) {
	    
	    cout << "Updating file." << endl;
	    write_results(out_file);
	    last_file_time=MPI_Wtime();
	    
	    size_t tc=0,conv2_count=0;
	    for(size_t ii=0;ii<n_nB2;ii++) {
	      for(size_t j=0;j<n_Ye2;j++) {
		for(size_t k=0;k<n_T2;k++) {
                  vector<size_t> iix={ii,j,k};
		  if (tg_flag.get(iix)>=10.0) conv2_count++;
		  tc++;
		}
	      }
	    }
	    
	    string msg="Table "+out_file+" is "+
	      o2scl::dtos(((double)conv2_count)/((double)tc)*100.0)+
	      " percent completed.";
	    slack.send(msg);
	    
	  }
	  
	  // End of i loop 
	}
	
	// End of if (mpi_size>1) 
      }

#endif

      // If we're running without MPI, or if we're running MPI with
      // only one rank, then complete the tasks on this MPI rank
      
      if (mpi_size==1) {

	for(size_t i=0;i<ntasks;i++) {

	  size_t jnB=tasks[i*6];
	  size_t jYe=tasks[i*6+1];
	  size_t jT=tasks[i*6+2];
	  size_t inB=tasks[i*6+3];
	  size_t iYe=tasks[i*6+4];
	  size_t iT=tasks[i*6+5];
	  cout << "Computing " 
	       << i << "/" << ntasks << " (" << tasks[i*6] << ","
	       << tasks[i*6+1] << "," 
	       << tasks[i*6+2] << ") -> ("
	       << tasks[i*6+3] << "," 
	       << tasks[i*6+4] << ","
	       << tasks[i*6+5] << ")" << endl;

	  double nB=nB_grid2[inB];
	  double Ye=Ye_grid2[iYe];
	  double T=T_grid2[iT]/hc_mev_fm;
	  cout << "  nB,Ye,T[MeV]: " << nB << " " << Ye << " " << T*hc_mev_fm
	       << endl;
	  double log_xn=gtab.get("log_xn",i);
	  double log_xp=gtab.get("log_xp",i);
	  size_t nuc_Z1=((size_t)(gtab.get("Z",i)+1.0e-12));
	  size_t nuc_N1=((size_t)(gtab.get("A",i)+1.0e-12))-nuc_Z1;
	  bool no_nuclei=false;
          
          // AWS, 1/16/21: I'm taking this out as it is confusing and
          // can result in incorrect results if not used carefully
	  if (false) {
	    if (iT>0) {
              vector<size_t> ixm3={inB,iYe,iT-1};
	      if (tg_flag.get(ixm3)>9.9) {
		double X_all_nuclei=0.0;
		X_all_nuclei+=tg_Xalpha.get(ixm3);
		X_all_nuclei+=tg_Xnuclei.get(ixm3);
		X_all_nuclei+=tg_Xd.get(ixm3);
		X_all_nuclei+=tg_Xt.get(ixm3);
		X_all_nuclei+=tg_XHe3.get(ixm3);
		X_all_nuclei+=tg_XLi4.get(ixm3);
		if (X_all_nuclei<1.0e-20) {
		  cout << "X_all_nuclei small (2): " << X_all_nuclei << " "
		       << T_grid2[iT] << endl;
		  no_nuclei=true;
		}
	      }
	    }
	  }
	    
	  thermo thx;
	  double mun_full, mup_full, Zbar, Nbar;
	  int A_min=0, A_max=0, NmZ_min=0, NmZ_max=0;
	  map<string,double> vdet;
          vdet["mue"]=gtab.get("mue",i);
	  
	  int ret=0;
	  if (alg_mode==1) {
	    cerr << "Mode alg_mode=1 no longer supported in generate-table."
		 << endl;
	    return 1;
	  } else if (alg_mode==0) {
	    ret=eos_vary_ZN(nB,Ye,T,log_xn,log_xp,nuc_Z1,nuc_N1,
			    thx,mun_full,mup_full,no_nuclei);
	    Zbar=nuc_Z1;
	    Nbar=nuc_N1;
	  } else if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
	    A_min=gtab.get("A_min",i);
	    A_max=gtab.get("A_max",i);
	    NmZ_min=gtab.get("NmZ_min",i);
	    NmZ_max=gtab.get("NmZ_max",i);
	    ret=eos_vary_dist
	      (nB,Ye,T,log_xn,log_xp,Zbar,Nbar,thx,mun_full,mup_full,
	       A_min,A_max,NmZ_min,NmZ_max,vdet,true,no_nuclei);    
	  }
	  if (gt_verbose>1) {
            cout.precision(5);
	    cout << "Point at (nB,Ye,T[MeV])=("
		 << nB << "," << Ye << "," << T*hc_mev_fm << "), ret="
		 << ret << ", " << i << "/" << ntasks << endl;
            cout.precision(6);
	  }
	  if (ret==0) {
	    
            one_success=true;
	    ubvector X;
	    compute_X(nB,X);
	    
	    store_point(inB,iYe,iT,nB,Ye,T,thx,log_xn,log_xp,
			Zbar,Nbar,mun_full,mup_full,X,A_min,A_max,
			NmZ_min,NmZ_max,10.0,vdet);
	  }

	  
#ifdef NO_MPI
	  double curr_time=time(0);
#else
	  double curr_time=MPI_Wtime();
#endif
	  if (((int)i)%(file_update_iters)==file_update_iters-1 ||
	      curr_time-last_file_time>file_update_time) {
	    
	    cout << "Updating file." << endl;
            write_results(out_file);
#ifdef NO_MPI
	    last_file_time=time(0);
#else
	    last_file_time=MPI_Wtime();
#endif
            
            size_t tc=0,conv2_count=0;
            for(size_t ii=0;ii<n_nB2;ii++) {
              for(size_t j=0;j<n_Ye2;j++) {
                for(size_t k=0;k<n_T2;k++) {
                  vector<size_t> iix={ii,j,k};
                  if (tg_flag.get(iix)>=10.0) conv2_count++;
                  tc++;
                }
              }
            }
            
            string msg="Table "+out_file+" is "+
              o2scl::dtos(((double)conv2_count)/((double)tc)*100.0)+
              " percent completed.";
            slack.send(msg);
            

	  }
	    
	  // End of loop over tasks
	}

	// End of 'if (mpi_size==1)'
      }

      // -------------------------------------------------------
      // The end of the main loop for mpi_rank = 0
      
      current_tasks+=ntasks;
      if (gt_verbose>0) {
	cout << "Rank " << mpi_rank << " computed " << current_tasks
	     << " out of " << total_tasks << ". "
	     << ((int)(((double)current_tasks)/
		       ((double)total_tasks)*100.0))
	     << " percent done." << endl;
      }

      if (one_success==false) {
        if (gt_verbose>0) {
          cout << "Rank " << mpi_rank << " found no successes. Stopping."
               << endl;
        }
        done=true;
      }
      
      double elapsed, write_elapsed;
#ifdef NO_MPI
      elapsed=time(0)-mpi_start_time;
      write_elapsed=time(0)-last_file_time;
#else
      elapsed=MPI_Wtime()-mpi_start_time;
      write_elapsed=MPI_Wtime()-last_file_time;
#endif
      
      if (ntasks==0 || (max_time>0.0 && elapsed>max_time)) {
	
	cout << "Finished. " << ntasks << " " << max_time << " "
	     << elapsed << endl;
	done=true;
	
      } else if (write_elapsed>file_update_time) {
	
	cout << "Updating file." << endl;
	write_results(out_file);
#ifdef NO_MPI	
	last_file_time=time(0);
#else
	last_file_time=MPI_Wtime();
#endif

	// Go through the entire table and count how many points have
	// finished
	
	size_t tc=0, conv2_count=0;
	for(size_t i=0;i<n_nB2;i++) {
	  for(size_t j=0;j<n_Ye2;j++) {
	    for(size_t k=0;k<n_T2;k++) {
              vector<size_t> ix={i,j,k};
	      if (tg_flag.get(ix)>=10.0) conv2_count++;
	      tc++;
	    }
	  }
	}
	
	string msg="Table "+out_file+" is "+
	  o2scl::dtos(((double)conv2_count)/((double)tc)*100.0)+
	  " percent completed.";
	slack.send(msg);
      }
      
      // For debugging
      // done=true;

      // Go back to see if we have any remaining tasks
      
      // End of 'while (done==false)'
    }

#ifndef NO_MPI
    
    // -----------------------------------------------------
    // Send exit to other processors

    for(int iproc=0;iproc<mpi_size-1;iproc++) {
      double buffer[buf_size];
      buffer[vi["msg"]]=((double)message_done);
      MPI_Send(&buffer[0],buf_size,MPI_DOUBLE,iproc+1,0,MPI_COMM_WORLD);
    }
    
#endif

    // -----------------------------------------------------
    // Count number of points finished

    size_t conv2_count=0;
    for(size_t i=0;i<n_nB2;i++) {
      for(size_t j=0;j<n_Ye2;j++) {
	for(size_t k=0;k<n_T2;k++) {
          vector<size_t> ix={i,j,k};
	  if (tg_flag.get(ix)>=10.0) conv2_count++;
	}
      }
    }
    cout << "There are " << conv2_count << " total points finished "
	 << "out of " << total_tasks << endl;

    // -----------------------------------------------------
    // Output file
    
    write_results(out_file);
    
    if (gt_verbose>0) {
      cout << "Rank " << mpi_rank << " sending exit to children." << endl;
    }
      
    string msg="Function generate_table() done. Wrote file "+
      out_file+" . There are "+
      o2scl::szttos(conv2_count)+" points finished out of "+
      o2scl::szttos(total_tasks)+".";
    slack.send(msg);
    
    // End of else for 'if (mpi_rank==0)'
    
  } else {

#ifndef NO_MPI
    
    // -----------------------------------------------------
    // Child processors

    double input_buffer[buf_size], output_buffer[buf_size];
    MPI_Recv(&input_buffer[0],buf_size,MPI_DOUBLE,0,0,
	     MPI_COMM_WORLD,&status);
      
    int message=((int)(input_buffer[vi["msg"]]+1.0e-6));

    if (false && gt_verbose>1) {
      cout << "Rank " << mpi_rank << " got message "
	   << message << endl;
    }
    
    thermo thx;
    double mun_full;
    double mup_full;
    
    while (message!=message_done) {

      // Obtain the task number and initial guess from the input buffer
      double task_index_d=input_buffer[vi["index"]];
      double nB=input_buffer[vi["nB"]];
      double Ye=input_buffer[vi["Ye"]];
      double T=input_buffer[vi["T"]]/hc_mev_fm;
      double log_xn=input_buffer[vi["log_xn"]];
      double log_xp=input_buffer[vi["log_xp"]];
      size_t nuc_Z1=((size_t)(input_buffer[vi["Z"]]+1.0e-12));
      size_t nuc_N1=((size_t)(input_buffer[vi["N"]]+1.0e-12));
      bool no_nuclei=input_buffer[vi["no_nuclei"]];
      int A_min, A_max, NmZ_min, NmZ_max;
      double Zbar, Nbar;
      if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
	A_min=((int)(input_buffer[vi["A_min"]]*(1.0+1.0e-12)));
	A_max=((int)(input_buffer[vi["A_max"]]*(1.0+1.0e-12)));
	NmZ_min=((int)(input_buffer[vi["NmZ_min"]]*(1.0+1.0e-12)));
	NmZ_max=((int)(input_buffer[vi["NmZ_max"]]*(1.0+1.0e-12)));
	Zbar=input_buffer[vi["Z"]];
	Nbar=input_buffer[vi["N"]];
      }
      
      int ret=-10;
      double elapsed=MPI_Wtime()-mpi_start_time;
      map<string,double> vdet;
      
      if (max_time==0.0 || elapsed<max_time) {

	if (gt_verbose>1) {
	  cout << "Rank " << mpi_rank
	       << " computing point at nB,Ye,T[MeV]: " << nB << " " 
	       << Ye << " " << T*hc_mev_fm << endl;
	}
	
        if (include_muons) {
          vdet["mue"]=input_buffer[vi["mue"]];
          vdet["Ymu"]=input_buffer[vi["Ymu"]];
        }
	
	if (alg_mode==1) {
	  cerr << "Mode alg_mode=1 no longer supported in generate-table."
	       << endl;
	  return 1;
	} else if (alg_mode==0) {
	  ret=eos_vary_ZN(nB,Ye,T,log_xn,log_xp,nuc_Z1,nuc_N1,
			  thx,mun_full,mup_full,no_nuclei);
	} else if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
	  ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
			    thx,mun_full,mup_full,
			    A_min,A_max,NmZ_min,NmZ_max,vdet,true,no_nuclei);
	}

	if (gt_verbose>1) {
	  cout.precision(4);
	  if (ret==0) {
	    cout << "Rank " << mpi_rank
		 << " done. nB,Ye,T[MeV],F[MeV]: " << nB << " " 
		 << Ye << " " << T*hc_mev_fm << " "
		 << (thx.ed-thx.en*T)/nB*hc_mev_fm << endl;
	    cout << "Rank " << mpi_rank
		 << " done. (inB,iYe,iT), Z,A: ("
		 << ((size_t)(input_buffer[vi["inB"]]*1.0+1.0e-12)) << ","
		 << ((size_t)(input_buffer[vi["iYe"]]*1.0+1.0e-12)) << ","
		 << ((size_t)(input_buffer[vi["iT"]]*1.0+1.0e-12)) << "), ";
	    if (alg_mode<2) {
	      cout << nuc_Z1 << " "
		   << nuc_Z1+nuc_N1 << endl;
	    } else {
	      cout << Zbar << " "
		   << Zbar+Nbar << endl;
	    }
	  } else {
	    cout << "Rank " << mpi_rank
		 << " failed. nB,Ye,T[MeV],ret: " << nB << " " 
		 << Ye << " " << T*hc_mev_fm << " "
		 << ret << endl;
	  }
	  cout.precision(6);
	}
	
	// End of 'if (max_time==0.0 || elapsed<max_time)'
      } else {
	// If we're out of time, don't compute anything, just return
	ret=o2scl::exc_efailed;
      }

      // Send results back to parent processor
      output_buffer[vi["index"]]=task_index_d;
      output_buffer[vi["nB"]]=nB;
      output_buffer[vi["Ye"]]=Ye;
      output_buffer[vi["T"]]=T*hc_mev_fm;
      if (ret==0) {
	output_buffer[vi["flag"]]=10.0;
	output_buffer[vi["log_xn"]]=log_xn;
	output_buffer[vi["log_xp"]]=log_xp;
	if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
	  output_buffer[vi["Z"]]=Zbar;
	  output_buffer[vi["N"]]=Nbar;
	  output_buffer[vi["A_min"]]=A_min;
	  output_buffer[vi["A_max"]]=A_max;
	  output_buffer[vi["NmZ_min"]]=NmZ_min;
	  output_buffer[vi["NmZ_max"]]=NmZ_max;
	} else {
	  output_buffer[vi["Z"]]=nuc_Z1;
	  output_buffer[vi["N"]]=nuc_N1;
	}
        if (include_muons) {
	  output_buffer[vi["mue"]]=vdet["mue"];
	  output_buffer[vi["Ymu"]]=vdet["Ymu"];
        }
	if (include_detail) {
	  //output_buffer[vi["zn"]]=vdet["zn"];
	}
	//cout << "Sending back: " << f_min << " "
	//<< f_min/nB*hc_mev_fm << endl;
	output_buffer[vi["fr"]]=thx.ed-T*thx.en;
	output_buffer[vi["ed"]]=thx.ed;
	output_buffer[vi["pr"]]=thx.pr;
	output_buffer[vi["en"]]=thx.en;
	if (nuclei.size()<5) {
	  O2SCL_ERR("Nuclei array improperly sized",o2scl::exc_esanity);
	}
	output_buffer[vi["Xalpha"]]=nuclei[0].n*4.0/nB;
	output_buffer[vi["Xd"]]=nuclei[1].n*2.0/nB;
	output_buffer[vi["Xt"]]=nuclei[2].n*3.0/nB;
	output_buffer[vi["XHe3"]]=nuclei[3].n*3.0/nB;
	output_buffer[vi["XLi4"]]=nuclei[4].n*4.0/nB;
	if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
	  double tmp=0.0;
	  for(size_t i=5;i<nuclei.size();i++) {
	    tmp+=nuclei[i].n*(nuclei[i].Z+nuclei[i].N)/nB;
	  }
	  output_buffer[vi["Xnuclei"]]=tmp;
	} else {
	  output_buffer[vi["Xnuclei"]]=nuclei[5].n*
	    (nuclei[5].Z+nuclei[5].N)/nB;
	}
      } else {
	output_buffer[vi["flag"]]=0.0;
      }
      
      MPI_Send(&output_buffer[0],buf_size,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
	
      // Obtain next instruction
      MPI_Recv(&input_buffer[0],buf_size,MPI_DOUBLE,0,0,
	       MPI_COMM_WORLD,&status);
      message=((int)(input_buffer[vi["msg"]]+1.0e-6));
    }

    if (gt_verbose>0) {
      cout << "Rank " << mpi_rank << " received done message." << endl;
    }
    
#endif
    
  }    
  
  return 0;

}

int eos_nuclei::check_virial(std::vector<std::string> &sv,
			     bool itive_com) {

  size_t n_nB=0, n_Ye=0, n_T=0;
  vector<double> nB_grid, Ye_grid, T_grid;
  
  n_nB=301;
  n_Ye=70;
  n_T=160;

  calc_utf8<> calc;
  std::map<std::string,double> vars;

  vector<double> packed, packed2;
  vector<std::string> split_res;
  
  o2scl::tensor_grid<> tg_zn2;
  o2scl::tensor_grid<> tg_zp2;

  split_string_delim(nB_grid_spec,split_res,',');
  n_nB2=stoszt(split_res[0]);
  
  calc.compile(split_res[1].c_str());
  for(size_t i=0;i<n_nB;i++) {
    vars["i"]=((double)i);
    nB_grid.push_back(calc.eval(&vars));
    packed.push_back(nB_grid[i]);
  }
  for(size_t i=0;i<21;i++) {
    packed2.push_back(pow(10.0,-50.0+i*2));
  }
    
  split_string_delim(Ye_grid_spec,split_res,',');
  n_Ye2=stoszt(split_res[0]);
  
  calc.compile(split_res[1].c_str());
  for(size_t i=0;i<n_Ye;i++) {
    vars["i"]=((double)i);
    Ye_grid.push_back(calc.eval(&vars));
    packed.push_back(Ye_grid[i]);
    packed2.push_back(Ye_grid[i]);
  }
    
  split_string_delim(T_grid_spec,split_res,',');
  n_T2=stoszt(split_res[0]);
  
  calc.compile(split_res[1].c_str());
  for(size_t i=0;i<n_T;i++) {
    vars["i"]=((double)i);
    T_grid.push_back(calc.eval(&vars));
    packed.push_back(T_grid[i]);
    packed2.push_back(T_grid[i]);
  }

  size_t st[3]={n_nB,n_Ye,n_T};
  tg_zn.resize(3,st);
  tg_zp.resize(3,st);

  tg_zn.set_grid_packed(packed);
  tg_zp.set_grid_packed(packed);

  tg_zn.set_all(0.0);
  tg_zp.set_all(0.0);
  
  size_t st2[3]={21,n_Ye,n_T};
  tg_zn2.resize(3,st2);
  tg_zp2.resize(3,st2);

  tg_zn2.set_grid_packed(packed2);
  tg_zp2.set_grid_packed(packed2);

  tg_zn2.set_all(0.0);
  tg_zp2.set_all(0.0);
  
  for(int inB=0;inB<((int)n_nB);inB++) {
    double nB=nB_grid[inB];
    for(int iYe=0;iYe<((int)n_Ye);iYe++) {
      double Ye=Ye_grid[iYe];
      for(int iT=0;iT<((int)n_T);iT++) {
        vector<size_t> ix={((size_t)inB),((size_t)iYe),((size_t)iT)};
	double T=T_grid[iT];
	double lambda=sqrt(4.0*o2scl_const::pi/(neutron.m+proton.m)/T);
	
	double b_n=ecv.bn_f(T);
	double b_pn=ecv.bpn_f(T);
	
	double zn, zp;
	//cout << "A: " << nB << " " << Ye << " " << T << " " << endl;
	vsd.solve_fugacity(nB*(1.0-Ye),nB*Ye,lambda,
			   lambda,b_n,b_pn,zn,zp);
	
	tg_zn.get(ix)=zn;
	tg_zp.get(ix)=zp;
      }
    }
  }

  for(int inB=0;inB<21;inB++) {
    double nB=pow(10.0,-50.0+inB*2);
    for(int iYe=0;iYe<((int)n_Ye);iYe++) {
      double Ye=Ye_grid[iYe];
      for(int iT=0;iT<((int)n_T);iT++) {
        vector<size_t> ix={((size_t)inB),((size_t)iYe),((size_t)iT)};
        
	double T=T_grid[iT];
	double lambda=sqrt(4.0*o2scl_const::pi/(neutron.m+proton.m)/T);
	
	double b_n=ecv.bn_f(T);
	double b_pn=ecv.bpn_f(T);
	
	double zn, zp;
	//cout << "B: " << nB << " " << Ye << " " << T << " " << endl;
	vsd.solve_fugacity(nB*(1.0-Ye),nB*Ye,lambda,
			   lambda,b_n,b_pn,zn,zp);

	tg_zn2.get(ix)=log10(zn);
	tg_zp2.get(ix)=log10(zp);
      }
    }
  }

  hdf_file hf;
  hf.open_or_create("check_virial.o2");
  hdf_output(hf,tg_zn,"zn");
  hdf_output(hf,tg_zp,"zp");
  hdf_output(hf,tg_zn2,"zn2");
  hdf_output(hf,tg_zp2,"zp2");
  hf.close();
  cout << "Created check_virial.o2" << endl;
  
  return 0;
}

void eos_nuclei::setup_cli(o2scl::cli &cl) {
  
  eos::setup_cli(cl,false);
  
  static const int nopt=30;

  o2scl::comm_option_s options[nopt]=
    {{0,"eos-deriv","",0,0,"","",
       new o2scl::comm_option_mfptr<eos_nuclei>
       (this,&eos_nuclei::eos_deriv),o2scl::cli::comm_option_both,
       1,"","eos_nuclei","eos_deriv","doc/xml/classeos__nuclei.xml"},
     {0,"add-eg","",0,0,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::add_eg),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","add_eg","doc/xml/classeos__nuclei.xml"},
     {0,"eg-table","",1,1,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::eg_table),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","eg_table","doc/xml/classeos__nuclei.xml"},
     {0,"maxwell","",0,0,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::maxwell),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","maxwell","doc/xml/classeos__nuclei.xml"},
     {0,"fit-frdm","",0,0,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::fit_frdm),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","fit_frdm","doc/xml/classeos__nuclei.xml"},
     {0,"check-virial","",0,0,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::check_virial),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","check_virial","doc/xml/classeos__nuclei.xml"},
     {0,"generate-table","",0,1,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::generate_table),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","generate_table","doc/xml/classeos__nuclei.xml"},
     {0,"test-random","",1,2,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::test_random),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","test_random","doc/xml/classeos__nuclei.xml"},
     {0,"load","",0,1,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::load),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","load","doc/xml/classeos__nuclei.xml"},
     {0,"output","",0,1,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::output),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","output","doc/xml/classeos__nuclei.xml"},
     {0,"edit-data","",1,4,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::edit_data),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","edit_data","doc/xml/classeos__nuclei.xml"},
     {0,"merge-tables","",3,3,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::merge_tables),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","merge_tables","doc/xml/classeos__nuclei.xml"},
     {0,"write-nuclei","",1,1,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::write_nuclei),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","write_nuclei","doc/xml/classeos__nuclei.xml"},
     {0,"compare-tables","",2,3,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::compare_tables),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","compare_tables","doc/xml/classeos__nuclei.xml"},
     {0,"interp-point","",4,4,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::interp_point),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","interp_point","doc/xml/classeos__nuclei.xml"},
     {0,"stats","",0,0,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::stats),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","stats","doc/xml/classeos__nuclei.xml"},
     {0,"mcarlo-nuclei","",0,0,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::mcarlo_nuclei),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","mcarlo_nuclei","doc/xml/classeos__nuclei.xml"},
     {0,"mcarlo-nuclei2","",5,5,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::mcarlo_nuclei2),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","mcarlo_nuclei2","doc/xml/classeos__nuclei.xml"},
     {0,"mcarlo-beta","",1,2,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::mcarlo_beta),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","mcarlo_beta","doc/xml/classeos__nuclei.xml"},
     {0,"mcarlo-neutron","",1,2,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::mcarlo_neutron),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","mcarlo_neutron","doc/xml/classeos__nuclei.xml"},
     {0,"point-nuclei-mu","",-1,-1,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::point_nuclei_mu),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","point_nuclei_mu","doc/xml/classeos__nuclei.xml"},
      {0,"muses-table","",1,1,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::muses_table),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","muses_table","doc/xml/classeos__nuclei.xml"},
      {0,"create-new-table","",-1,-1,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::create_new_table),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","create_new_table","doc/xml/classeos__nuclei.xml"},
     {0,"point-nuclei","",-1,-1,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::point_nuclei),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","point_nuclei","doc/xml/classeos__nuclei.xml"},
      {0,"muses","",-1,-1,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::muses),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","muses","doc/xml/classeos__nuclei.xml"},
     {0,"increase-density","",7,7,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::increase_density),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","increase_density","doc/xml/classeos__nuclei.xml"},
     {0,"fix-cc","",1,1,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::fix_cc),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","fix_cc","doc/xml/classeos__nuclei.xml"},
     {0,"stability","",1,3,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::stability),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","stability","doc/xml/classeos__nuclei.xml"},
     {0,"verify","",2,4,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::verify),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","verify","doc/xml/classeos__nuclei.xml"},
     {0,"select-high-T","",1,1,"","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::select_high_T),o2scl::cli::comm_option_both,
      1,"","eos_nuclei","select_high_T","doc/xml/classeos__nuclei.xml"}
    };

  cl.doc_o2_file="data/eos_nuclei_docs.o2";

  p_nB_grid_spec.str=&nB_grid_spec;
  p_nB_grid_spec.help="";
  p_nB_grid_spec.doc_class="eos_nuclei";
  p_nB_grid_spec.doc_name="nB_grid_spec";
  p_nB_grid_spec.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("nB_grid_spec",&p_nB_grid_spec));
  
  p_Ye_grid_spec.str=&Ye_grid_spec;
  p_Ye_grid_spec.help="";
  p_Ye_grid_spec.doc_class="eos_nuclei";
  p_Ye_grid_spec.doc_name="Ye_grid_spec";
  p_Ye_grid_spec.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("Ye_grid_spec",&p_Ye_grid_spec));
  
  p_T_grid_spec.str=&T_grid_spec;
  p_T_grid_spec.help="";
  p_T_grid_spec.doc_class="eos_nuclei";
  p_T_grid_spec.doc_name="T_grid_spec";
  p_T_grid_spec.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("T_grid_spec",&p_T_grid_spec));
  
  p_show_all_nuclei.b=&show_all_nuclei;
  p_show_all_nuclei.help="";
  p_show_all_nuclei.doc_class="eos_nuclei";
  p_show_all_nuclei.doc_name="show_all_nuclei";
  p_show_all_nuclei.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("show_all_nuclei",&p_show_all_nuclei));
  
  p_extend_frdm.b=&extend_frdm;
  p_extend_frdm.help="";
  p_extend_frdm.doc_class="eos_nuclei";
  p_extend_frdm.doc_name="extend_frdm";
  p_extend_frdm.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("extend_frdm",&p_extend_frdm));
  
  p_survey_eqs.b=&survey_eqs;
  p_survey_eqs.help="";
  p_survey_eqs.doc_class="eos_nuclei";
  p_survey_eqs.doc_name="survey_eqs";
  p_survey_eqs.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("survey_eqs",&p_survey_eqs));
  
  p_fd_A_max.i=&fd_A_max;
  p_fd_A_max.help="";
  p_fd_A_max.doc_class="eos_nuclei";
  p_fd_A_max.doc_name="fd_A_max";
  p_fd_A_max.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("fd_A_max",&p_fd_A_max));
  
  p_recompute.b=&recompute;
  p_recompute.help="";
  p_recompute.doc_class="eos_nuclei";
  p_recompute.doc_name="recompute";
  p_recompute.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("recompute",&p_recompute));
  
  p_verify_only.b=&verify_only;
  p_verify_only.help="";
  p_verify_only.doc_class="eos_nuclei";
  p_verify_only.doc_name="verify_only";
  p_verify_only.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("verify_only",&p_verify_only));
  
  p_edge_list.str=&edge_list;
  p_edge_list.help="";
  p_edge_list.doc_class="eos_nuclei";
  p_edge_list.doc_name="edge_list";
  p_edge_list.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("edge_list",&p_edge_list));
  
  p_six_neighbors.i=&six_neighbors;
  p_six_neighbors.help="";
  p_six_neighbors.doc_class="eos_nuclei";
  p_six_neighbors.doc_name="six_neighbors";
  p_six_neighbors.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("six_neighbors",&p_six_neighbors));
  
  p_rnuc_less_rws.b=&rnuc_less_rws;
  p_rnuc_less_rws.help="";
  p_rnuc_less_rws.doc_class="eos_nuclei";
  p_rnuc_less_rws.doc_name="rnuc_less_rws";
  p_rnuc_less_rws.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("rnuc_less_rws",&p_rnuc_less_rws));

  p_propagate_points.b=&propagate_points;
  p_propagate_points.help="";
  p_propagate_points.doc_class="eos_nuclei";
  p_propagate_points.doc_name="propagate_points";
  p_propagate_points.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("propagate_points",&p_propagate_points));
  
  p_mh_tol_rel.d=&mh_tol_rel;
  p_mh_tol_rel.help="";
  p_mh_tol_rel.doc_class="eos_nuclei";
  p_mh_tol_rel.doc_name="mh_tol_rel";
  p_mh_tol_rel.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("mh_tol_rel",&p_mh_tol_rel));
  
  p_max_time.d=&max_time;
  p_max_time.help="";
  p_max_time.doc_class="eos_nuclei";
  p_max_time.doc_name="max_time";
  p_max_time.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("max_time",&p_max_time));
  
  p_nucleon_func.str=&nucleon_func;
  p_nucleon_func.help="";
  p_nucleon_func.doc_class="eos_nuclei";
  p_nucleon_func.doc_name="nucleon_func";
  p_nucleon_func.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("nucleon_func",&p_nucleon_func));

  p_ext_guess.str=&ext_guess;
  p_ext_guess.help="";
  p_ext_guess.doc_class="eos_nuclei";
  p_ext_guess.doc_name="ext_guess";
  p_ext_guess.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("ext_guess",&p_ext_guess));

  p_Ye_list.str=&Ye_list;
  p_Ye_list.help="";
  p_Ye_list.doc_class="eos_nuclei";
  p_Ye_list.doc_name="Ye_list";
  p_Ye_list.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("Ye_list",&p_Ye_list));

  p_alg_mode.i=&alg_mode;
  p_alg_mode.help="";
  p_alg_mode.doc_class="eos_nuclei";
  p_alg_mode.doc_name="alg_mode";
  p_alg_mode.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("alg_mode",&p_alg_mode));
  
  p_cs2_verbose.i=&cs2_verbose;
  p_cs2_verbose.help="";
  p_cs2_verbose.doc_class="eos_nuclei";
  p_cs2_verbose.doc_name="cs2_verbose";
  p_cs2_verbose.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("cs2_verbose",&p_cs2_verbose));
  
  p_fixed_dist_alg.i=&fixed_dist_alg;
  p_fixed_dist_alg.help="";
  p_fixed_dist_alg.doc_class="eos_nuclei";
  p_fixed_dist_alg.doc_name="fixed_dist_alg";
  p_fixed_dist_alg.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("fixed_dist_alg",&p_fixed_dist_alg));
  
  p_function_verbose.i=&function_verbose;
  p_function_verbose.help="";
  p_function_verbose.doc_class="eos_nuclei";
  p_function_verbose.doc_name="function_verbose";
  p_function_verbose.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("function_verbose",&p_function_verbose));
  
  p_max_ratio.d=&max_ratio;
  p_max_ratio.help="";
  p_max_ratio.doc_class="eos_nuclei";
  p_max_ratio.doc_name="max_ratio";
  p_max_ratio.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("max_ratio",&p_max_ratio));

  p_file_update_time.d=&file_update_time;
  p_file_update_time.help="";
  p_file_update_time.doc_class="eos_nuclei";
  p_file_update_time.doc_name="file_update_time";
  p_file_update_time.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("file_update_time",&p_file_update_time));

  p_file_update_iters.i=&file_update_iters;
  p_file_update_iters.help="";
  p_file_update_iters.doc_class="eos_nuclei";
  p_file_update_iters.doc_name="file_update_iters";
  p_file_update_iters.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("file_update_iters",&p_file_update_iters));

  p_include_detail.b=&include_detail;
  p_include_detail.help="";
  p_include_detail.doc_class="eos_nuclei";
  p_include_detail.doc_name="include_detail";
  p_include_detail.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("include_detail",&p_include_detail));

  p_strange_axis.b=&strange_axis;
  p_strange_axis.help="";
  p_strange_axis.doc_class="eos_nuclei";
  p_strange_axis.doc_name="strange_axis";
  p_strange_axis.doc_xml_file="doc/xml/classeos.xml";
  cl.par_list.insert(make_pair("strange_axis",&p_strange_axis));
  
  p_inc_hrg.b=&inc_hrg;
  p_inc_hrg.help="If true, include the hadron resonance gas (default false)";
  p_inc_hrg.doc_class="eos_nuclei";
  p_inc_hrg.doc_name="inc_hrg";
  p_inc_hrg.doc_xml_file="doc/xml/classeos__nuclei.xml";
  cl.par_list.insert(make_pair("inc_hrg",&p_inc_hrg));
  
  cl.set_comm_option_vec(nopt,options);
  
  cl.xml_subs.push_back("<formula> $ 10^{-6} $ </formula>");
  cl.xml_subs.push_back("10⁻⁶");
  cl.xml_subs.push_back("<formula> $ N/Z $ </formula>");
  cl.xml_subs.push_back("N/Z");
  cl.xml_subs.push_back("<formula> $ Z/N $ </formula>");
  cl.xml_subs.push_back("Z/N");
  cl.xml_subs.push_back("<computeroutput> dist.o2 </computeroutput>");
  cl.xml_subs.push_back("dist.o2");
  // Just ignore the function references
  cl.xml_subs.push_back("<ref>");
  cl.xml_subs.push_back("");
  cl.xml_subs.push_back("</ref>");
  cl.xml_subs.push_back("");
  
  if (file_exists(cl.doc_o2_file)) {
    //cl.set_verbose(3);
    cl.read_docs();
  } else {
    cerr << "Couldn't find file " << cl.doc_o2_file
         << " for run-time documentation." << endl;
  }

  return;
}
