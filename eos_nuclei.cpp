/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2020, Xingfu Du and Andrew W. Steiner
  
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

#include <o2scl/root_brent_gsl.h>
#include <o2scl/classical.h>
#include <o2scl/interp.h>
#include <o2scl/svd.h>
#include <o2scl/vector.h>
#include <o2scl/mmin_bfgs2.h>
#include <o2scl/diff_evo_adapt.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

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

  show_all_nuclei=false;
  recompute=false;
  verify_only=false;
  edge_list="";
  six_neighbors=0;
  propagate_points=true;
  
  mh.ntrial=10000;
  mh.err_nonconv=false;
  mh.def_jac.err_nonconv=false;
  mh.tol_rel=1.0e-6;
  mh_tol_rel=1.0e-6;
  rbg.err_nonconv=false;

  baryons_only_loaded=true;
  derivs_computed=false;
  with_leptons_loaded=false;

  //nucleon_func="(nb<0.04)*((i<100)*10+(i>=100)*sqrt(i))+(nb>=0.04)*100";
  nucleon_func="(i<100)*10+(i>=100)*sqrt(i)";
  
  vi.append({"msg","index","flag","nB","Ye","T","log_xn","log_xp",
	     "Z","N","fr","ed","pr","en","mun","mup","no_nuclei",
	     "inB","iYe","iT","A_min","A_max","NmZ_min","NmZ_max",
	     "Xalpha","Xd","Xt","XHe3","XLi4","Xnuclei",
	     "zn","zp","F1","F2","F3","F4","Un","Up","msn",
	     "msp","g","dgdT"});

  select_high_T(6);
  
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

  baryons_only_loaded=true;

  n_nB2=0;
  n_Ye2=0;
  n_T2=0;

  ext_guess="";
  include_detail=false;
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
  
#ifdef O2SCL_CORI
  // Load the nuclear masses
  cout << "Rank " << mpi_rank << " loading nuclear masses." << endl;  
  o2scl_hdf::ame_load_ext(ame,"data/ame16.o2","ame16.o2");
  o2scl_hdf::mnmsk_load(m95,"data/mnmsk.o2");
  o2scl_hdf::hfb_sp_load(hfb,27,"data");
  cout << "Rank " << mpi_rank << " finished loading nuclear masses."
       << endl;  
#else
  // Load the nuclear masses
  o2scl_hdf::ame_load(ame);
  o2scl_hdf::mnmsk_load(m95);
  o2scl_hdf::hfb_sp_load(hfb,27);
#endif
  
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
    ame.get_nucleus(nuclei[i].Z,nuclei[i].N-1,nuc_temp);
    Sneut[i]=-(nuclei[i].be-nuc_temp.be)*hc_mev_fm;
    ame.get_nucleus(nuclei[i].Z-1,nuclei[i].N,nuc_temp);
    Sprot[i]=-(nuclei[i].be-nuc_temp.be)*hc_mev_fm;  
  }
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
    Ec1[i]=-0.6*nuclei[i].Z*nuclei[i].Z*fine_structure/rA*
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
  double Ecalpha=-0.6*nuclei[0].Z*nuclei[0].Z*fine_structure
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
  double Echeavy=-0.6*nuclei[5].Z*nuclei[5].Z*fine_structure
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

int eos_nuclei::maxwell_test(std::vector<std::string> &sv,
			     bool itive_com) {

  if (sv.size()<5) {
    cerr << "Need input, Ye, T, and output files for \"maxwell-test\"" << endl;
    return 1;
  }
  
  read_results(sv[1]);
  
  if (with_leptons_loaded==false || derivs_computed==false) {
    cout << "maxwell_test only works with leptons and derivatives." << endl;
    return 0;
  }

  size_t iYe=vector_lookup(Ye_grid2,o2scl::function_to_double(sv[2]));
  size_t iT=vector_lookup(T_grid2,o2scl::function_to_double(sv[3]));
  cout << "iYe, iT: " << iYe << " " << iT << endl;
  cout << Ye_grid2[iYe] << " " << T_grid2[iT] << endl;

  table_units<> tu;
  tu.line_of_names("nB P mun");
  for(size_t i=0;i<n_nB2;i++) {
    double line[3]={nB_grid2[i],tg3_P.get(i,iYe,iT),
		    tg3_mun.get(i,iYe,iT)};
    tu.line_of_data(3,line);
  }

  cout << "Writing file: " << sv[4] << endl;
  hdf_file hf;
  hf.open_or_create(sv[4]);
  hdf_output(hf,tu,"mt");
  hf.close();
  
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

  if (tg3_E.total_size()==0) {
    tg3_E.resize(3,st);
    tg3_E.set_grid_packed(packed);
    tg3_E.set_all(0.0);
  }
  if (tg3_P.total_size()==0) {
    tg3_P.resize(3,st);
    tg3_P.set_grid_packed(packed);
    tg3_P.set_all(0.0);
  }
  if (tg3_S.total_size()==0) {
    tg3_S.resize(3,st);
    tg3_S.set_grid_packed(packed);
    tg3_S.set_all(0.0);
  }
  if (tg3_F.total_size()==0) {
    tg3_F.resize(3,st);
    tg3_F.set_grid_packed(packed);
    tg3_F.set_all(0.0);
  }
  if (tg3_mue.total_size()==0) {
    tg3_mue.resize(3,st);
    tg3_mue.set_grid_packed(packed);
    tg3_mue.set_all(0.0);
  }
  
  eos_sn_base eso;
  eso.include_muons=false;
  
  for (size_t i=0;i<n_nB2;i++) {
    for (size_t j=0;j<n_Ye2;j++) {
      for (size_t k=0;k<n_T2;k++) {
	
	thermo lep;
	double mue;
	eso.compute_eg_point(nB_grid2[i],Ye_grid2[j],T_grid2[k],lep,mue);
	
	tg3_F.set(i,j,k,tg3_Fint.get(i,j,k)+
		  (hc_mev_fm*lep.ed-T_grid2[k]*lep.en)/nB_grid2[i]);
	tg3_E.set(i,j,k,tg3_Eint.get(i,j,k)+
		  hc_mev_fm*lep.ed/nB_grid2[i]);
	tg3_P.set(i,j,k,tg3_Pint.get(i,j,k)+
		  hc_mev_fm*lep.pr);
	tg3_S.set(i,j,k,tg3_Sint.get(i,j,k)+
		  lep.en/nB_grid2[i]);
	tg3_mue.set(i,j,k,hc_mev_fm*mue);

      }
    }
    cout << "add_eg(): " << i+1 << "/" << n_nB2 << endl;
  }

  with_leptons_loaded=true;

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
  
  tensor_grid3<> mu_e;
  tensor_grid3<> E;
  tensor_grid3<> n_mu;
  tensor_grid3<> P;
  tensor_grid3<> S;
  tensor_grid3<> F;
  
  size_t st[3]={n_nB2,n_Ye2,n_T2};
  vector<vector<double> > grid={nB_grid2,Ye_grid2,T_grid2};

  mu_e.resize(3,st);
  n_mu.resize(3,st);
  E.resize(3,st);
  P.resize(3,st);
  S.resize(3,st);
  F.resize(3,st);

  mu_e.set_grid(grid);
  n_mu.set_grid(grid);
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
	double T_MeV=T_grid2[k];
	
	thermo lep;
	double mue;
	esb.compute_eg_point(nB_grid2[i],Ye_grid2[j],T_grid2[k],lep,mue);

	mu_e.set(i,j,k,mue*o2scl_const::hc_mev_fm);
	E.set(i,j,k,lep.ed/nB*o2scl_const::hc_mev_fm);
	P.set(i,j,k,lep.pr*o2scl_const::hc_mev_fm);
	S.set(i,j,k,lep.en/nB);
	F.set(i,j,k,(lep.ed-T_MeV*lep.en)/nB*o2scl_const::hc_mev_fm);
	if (include_muons) {
	  n_mu.set(i,j,k,muon.n*o2scl_const::hc_mev_fm);
	}
      }
    }
    cout << "eg_table(): " << i+1 << "/" << n_nB2 << endl;
  }

  hdf_file hf;
  hf.open_or_create(sv[1]);
  hdf_output(hf,mu_e,"mue");
  if (include_muons) {
    hdf_output(hf,n_mu,"nmu");
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
    
    calculator calc;
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
    
    tg3_Eint.resize(3,st);
    tg3_Pint.resize(3,st);
    tg3_Sint.resize(3,st);
    tg3_mun.resize(3,st);
    tg3_mup.resize(3,st);
    
    tg3_Eint.set_grid_packed(packed);
    tg3_Pint.set_grid_packed(packed);
    tg3_Sint.set_grid_packed(packed);
    tg3_mun.set_grid_packed(packed);
    tg3_mup.set_grid_packed(packed);
    
    tg3_Eint.set_all(0.0);
    tg3_Pint.set_all(0.0);
    tg3_Sint.set_all(0.0);
    tg3_mun.set_all(0.0);
    tg3_mup.set_all(0.0);
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
	fint_vec_T[k]=tg3_Fint.get(i,j,k)*nB_grid2[i];
      }
      itp_stf.set(n_T2,T_grid2,fint_vec_T);
      for(size_t k=0;k<n_T2;k++) {
	// Note that fint_vec_T above is stored in MeV/fm^3, so
	// when we take a derivative wrt to temperature (stored
	// in MeV), we get the correct units of 1/fm^3, and then
	// we divide by the baryon density in units of 1/fm^3
	tg3_Sint.set(i,j,k,-itp_stf.deriv(T_grid2[k])/nB_grid2[i]);
      }
    }
  }

  // Second, compute the electron fraction derivative, which we
  // temporarily store in tg3_mup
  for (size_t i=0;i<n_nB2;i++) {
    for(size_t k=0;k<n_T2;k++) {
      for (size_t j=0;j<n_Ye2;j++) {
	fint_vec_Ye[j]=tg3_Fint.get(i,j,k)*nB_grid2[i];
      }
      itp_stf.set(n_Ye2,Ye_grid2,fint_vec_Ye);
      for (size_t j=0;j<n_Ye2;j++) {
	tg3_mup.set(i,j,k,itp_stf.deriv(Ye_grid2[j]));
      }
    }
  }

  // Third, compute the baryon density derivative, which we
  // temporarily store in tg3_mun
  for (size_t j=0;j<n_Ye2;j++) {
    for(size_t k=0;k<n_T2;k++) {
      for (size_t i=0;i<n_nB2;i++) {
	fint_vec_nB[i]=tg3_Fint.get(i,j,k)*nB_grid2[i];
      }
      itp_stf.set(n_nB2,nB_grid2,fint_vec_nB);
      for (size_t i=0;i<n_nB2;i++) {
	tg3_mun.set(i,j,k,itp_stf.deriv(nB_grid2[i]));
      }
    }
  }

  // Now go through every point and compute the remaining
  // quantities
  for (size_t i=0;i<n_nB2;i++) {
    for (size_t j=0;j<n_Ye2;j++) {
      for (size_t k=0;k<n_T2;k++) {
	
	double en=tg3_Sint.get(i,j,k)*nB_grid2[i];
	double dfdnB=tg3_mun.get(i,j,k);
	double dfdYe=tg3_mup.get(i,j,k);

	if (false && nB_grid2[i]>0.3) {
	  cout << tg3_mun.get(i,j,k) << endl;
	  cout << tg3_mup.get(i,j,k) << endl;
	}
	
	tg3_mun.get(i,j,k)=dfdnB-dfdYe*Ye_grid2[j]/nB_grid2[i];
	tg3_mup.get(i,j,k)=dfdnB-dfdYe*(Ye_grid2[j]-1.0)/nB_grid2[i];

	if (false && nB_grid2[i]>0.3) {
	  cout << tg3_mun.get(i,j,k) << endl;
	  cout << tg3_mup.get(i,j,k) << endl;
	  exit(-1);
	}
	
	// E = F + T S
	tg3_Eint.get(i,j,k)=tg3_Fint.get(i,j,k)+
	  T_grid2[k]*tg3_Sint.get(i,j,k);
	
	// P = - F + mun * nn + mup * np
	tg3_Pint.get(i,j,k)=-tg3_Fint.get(i,j,k)*nB_grid2[i]+
	  nB_grid2[i]*(1.0-Ye_grid2[j])*tg3_mun.get(i,j,k)+
	  nB_grid2[i]*Ye_grid2[j]*tg3_mup.get(i,j,k);
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

  tensor_grid3<> tg3_dmundnB;
  tensor_grid3<> tg3_dmundYe;
  tensor_grid3<> tg3_dmundT;
  tensor_grid3<> tg3_dmupdYe;
  tensor_grid3<> tg3_dmupdT;
  tensor_grid3<> tg3_dsdT;
  
  derivs_computed=true;
  
  size_t st[3]={n_nB2,n_Ye2,n_T2};
  vector<vector<double> > grid={nB_grid2,Ye_grid2,T_grid2};
  
  tg3_dmundnB.resize(3,st);
  tg3_dmundYe.resize(3,st);
  tg3_dmundT.resize(3,st);
  tg3_dmupdYe.resize(3,st);
  tg3_dmupdT.resize(3,st);
  tg3_dsdT.resize(3,st);
  
  tg3_dmundnB.set_grid(grid);
  tg3_dmundYe.set_grid(grid);
  tg3_dmundT.set_grid(grid);
  tg3_dmupdYe.set_grid(grid);
  tg3_dmupdT.set_grid(grid);
  tg3_dsdT.set_grid(grid);
  
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
	mun_vec_T[k]=tg3_mun.get(i,j,k);
	mup_vec_T[k]=tg3_mup.get(i,j,k);
	s_vec_T[k]=tg3_S.get(i,j,k)/nB;
      }
      itp_stf.set(n_T2,T_grid2,mun_vec_T);
      for(size_t k=0;k<n_T2;k++) {
	tg3_dmundT.set(i,j,k,itp_stf.deriv(T_grid2[k]));
      }
      itp_stf.set(n_T2,T_grid2,mup_vec_T);
      for(size_t k=0;k<n_T2;k++) {
	tg3_dmupdT.set(i,j,k,itp_stf.deriv(T_grid2[k]));
      }
      itp_stf.set(n_T2,T_grid2,s_vec_T);
      for(size_t k=0;k<n_T2;k++) {
	tg3_dsdT.set(i,j,k,itp_stf.deriv(T_grid2[k]));
      }
    }
  }

  // Second, compute the Ye derivatives
  for (size_t i=0;i<n_nB2;i++) {
    for(size_t k=0;k<n_T2;k++) {
      for (size_t j=0;j<n_Ye2;j++) {
	mun_vec_Ye[j]=tg3_mun.get(i,j,k);
	mup_vec_Ye[j]=tg3_mup.get(i,j,k);
      }
      itp_stf.set(n_Ye2,Ye_grid2,mun_vec_Ye);
      for (size_t j=0;j<n_Ye2;j++) {
	tg3_dmundYe.set(i,j,k,itp_stf.deriv(Ye_grid2[j]));
      }
      itp_stf.set(n_Ye2,Ye_grid2,mup_vec_Ye);
      for (size_t j=0;j<n_Ye2;j++) {
	tg3_dmupdYe.set(i,j,k,itp_stf.deriv(Ye_grid2[j]));
      }
    }
  }

  // Third, compute the baryon density derivative
  for (size_t j=0;j<n_Ye2;j++) {
    for(size_t k=0;k<n_T2;k++) {
      for (size_t i=0;i<n_nB2;i++) {
	mun_vec_nB[i]=tg3_mun.get(i,j,k);
      }
      itp_stf.set(n_nB2,nB_grid2,mun_vec_nB);
      for (size_t i=0;i<n_nB2;i++) {
	tg3_dmundnB.set(i,j,k,itp_stf.deriv(nB_grid2[i]));
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

	double dmupdnB=tg3_dmundnB.get(i,j,k)+
	  tg3_dmundYe.get(i,j,k)*(1.0-Ye)/nB+
	  tg3_dmupdYe.get(i,j,k)*Ye/nB;
	  
	double f_nnnn=tg3_dmundnB.get(i,j,k)-
	  Ye*tg3_dmundYe.get(i,j,k)/nB;
	double f_nnnp=tg3_dmundnB.get(i,j,k)+
	  (1.0-Ye)*tg3_dmundYe.get(i,j,k)/nB;
	double f_npnp=dmupdnB+
	  (1.0-Ye)*tg3_dmupdYe.get(i,j,k)/nB;
	double f_nnT=tg3_dmundT.get(i,j,k);
	double f_npT=tg3_dmupdT.get(i,j,k);
	double f_TT=-tg3_dsdT.get(i,j,k);
	double nn=nB*(1.0-Ye);
	double np=nB*Ye;
	double en=tg3_Sint.get(i,j,k)*nB;

	double mue=0.0;
	double den=en*T_MeV+(tg3_mun.get(i,j,k)+neutron.m)*nn+
	  (tg3_mup.get(i,j,k)+proton.m+mue)*np;
	
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

  string outfile=sv[1];

  // The six second derivatives we need to compute
  o2scl::tensor_grid3<> dmundYe, dmundnB, dmupdYe, dsdT, dsdnB, dsdYe;
  o2scl::tensor_grid3<> egv[4], cs2;

  interp_vec<vector<double> > itp_sta_a, itp_sta_b, itp_sta_c;

  std::string in_file;
  read_results(in_file);

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
  dmundYe.resize(3,st);
  dmupdYe.resize(3,st);
  dsdT.resize(3,st);
  dsdnB.resize(3,st);
  dsdYe.resize(3,st);
  for(size_t i=0;i<4;i++) {
    egv[i].resize(3,st);
  }
  cs2.resize(3,st);

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

  // The baryon density derivatives
  for (size_t j=0;j<n_Ye2;j++) {
    for (size_t k=0;k<n_T2;k++) {
      vector<double> mun_of_nB, s_of_nB;
      for (size_t i=0;i<n_nB2;i++) {
	mun_of_nB.push_back(tg3_mun.get(i,j,k)/o2scl_const::hc_mev_fm);
	s_of_nB.push_back(tg3_Sint.get(i,j,k)*nB_grid2[i]);
      }
      itp_sta_a.set(n_nB2,nB_grid2,mun_of_nB,itp_steffen);
      itp_sta_b.set(n_nB2,nB_grid2,s_of_nB,itp_steffen);
      for (size_t i=0;i<n_nB2;i++) {
	dmundnB.get(i,j,k)=itp_sta_a.eval(nB_grid2[i]);
	dsdnB.get(i,j,k)=itp_sta_b.eval(nB_grid2[i]);
      }
    }
  }
  
  // The electron fraction derivatives
  for (size_t i=0;i<n_nB2;i++) {
    double nB=nB_grid2[i];
    for (size_t k=0;k<n_T2;k++) {
      vector<double> mun_of_Ye, mup_of_Ye, s_of_Ye;
      for (size_t j=0;j<n_Ye2;j++) {
	mun_of_Ye.push_back(tg3_mun.get(i,j,k)/o2scl_const::hc_mev_fm);
	mup_of_Ye.push_back(tg3_mup.get(i,j,k)/o2scl_const::hc_mev_fm);
	s_of_Ye.push_back(tg3_Sint.get(i,j,k)*nB);
      }
      itp_sta_a.set(n_Ye2,Ye_grid2,mun_of_Ye,itp_steffen);
      itp_sta_b.set(n_Ye2,Ye_grid2,mup_of_Ye,itp_steffen);
      itp_sta_c.set(n_Ye2,Ye_grid2,s_of_Ye,itp_steffen);
      for (size_t j=0;j<n_Ye2;j++) {
	dmundYe.get(i,j,k)=itp_sta_a.eval(Ye_grid2[j]);
	dmupdYe.get(i,j,k)=itp_sta_b.eval(Ye_grid2[j]);
	dsdYe.get(i,j,k)=itp_sta_b.eval(Ye_grid2[j]);
      }
    }
  }
  
  // The temperature derivative
  for (size_t i=0;i<n_nB2;i++) {
    double nB=nB_grid2[i];
    for (size_t j=0;j<n_Ye2;j++) {
      vector<double> s_of_T;
      for (size_t k=0;k<n_T2;k++) {
	s_of_T.push_back(tg3_Sint.get(i,j,k)*nB);
      }
      itp_sta_a.set(n_T2,T_grid2,s_of_T,itp_steffen);
      for (size_t k=0;k<n_T2;k++) {
	dsdT.get(i,j,k)=itp_sta_a.eval(T_grid2[k]);
      }
    }
  }

  // Storage for the matrix and SVD
  ubmatrix mat(4,4), V(4,4);
  ubvector sing(4), work(4);
  
  /// Compute the stability matrix and its eigenvalues at each point
  for (size_t i=0;i<n_nB2;i++) {
    double nB=nB_grid2[i];
    for (size_t j=0;j<n_Ye2;j++) {
      double Ye=Ye_grid2[j];
      for (size_t k=0;k<n_T2;k++) {
	double T_MeV=T_grid2[k];

	// Entropy and densities
	double en=tg3_S.get(i,j,k)*nB;
	double nn2=(1.0-Ye)*nB;
	double np2=Ye*nB;
	
	// Temporarily store the six derivatives
	double dmundnBv=dmundnB.get(i,j,k);
	double dmundYev=dmundYe.get(i,j,k);
	double dmupdYev=dmupdYe.get(i,j,k);
	double dsdnBv=dsdnB.get(i,j,k);
	double dsdYev=dsdYe.get(i,j,k);
	double dsdTv=dsdT.get(i,j,k);

	// Compute dmupdnB
	double dmupdnBv=Ye/nB*dmupdYev+dmundnBv+(1.0-Ye)/nB*dmundYev;

	// Transform to nn and np
	double dmundnn=dmundnBv-Ye/nB*dmundYev;
	double dmundnp=dmundnBv+(1.0-Ye)/nB*dmundYev;
	double dmupdnn=dmundnp;
	double dmupdnp=dmupdnBv+(1.0-Ye)/nB*dmupdYev;
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

	// Compute eigenvalues using SVD
	o2scl_linalg::SV_decomp(4,4,mat,V,sing,work);

	for(size_t ij=0;ij<4;ij++) {
	  egv[ij].get(i,j,k)=sing[ij];
	}

	// Compute squared speed of sound
	double f_nnnn=dmundnn;
	double f_nnnp=dmundnp;
	double f_npnp=dmupdnp;
	double f_nnT=dmundT;
	double f_npT=dmupdT;
	double f_TT=-dsdTv;
	double den=en*T_MeV+(tg3_mun.get(i,j,k)/hc_mev_fm+neutron.m)*nn2+
	  (tg3_mup.get(i,j,k)/hc_mev_fm+proton.m)*np2+electron.mu*electron.n;
	double cs_sq=(nn2*nn2*(f_nnnn-f_nnT*f_nnT/f_TT)+
		      2.0*nn2*np2*(f_nnnp-f_nnT*f_npT/f_TT)+
		      np2*np2*(f_npnp-f_npT*f_npT/f_TT)-
		      2.0*en*(nn2*f_nnT/f_TT+np2*f_npT/f_TT)-en*en/f_TT)/den;
	cs2.get(i,j,k)=cs_sq;
      }
    }
  }
  
  hdf_file hf;
  hf.open(outfile);
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
  hf.close();
  
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
    Ec[i]=-0.6*nuclei[i].Z*nuclei[i].Z*fine_structure/rA*
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

    if (use_skalt) {
      eosp_alt->calc_temp_e(neutron,proton,T,th2);
    } else {
      free_energy_density(neutron,proton,T,th2);
    }
    double mup_shift=proton.mu-
      T*log(1.0/proton.g*pow(2.0*pi/proton.ms/T,1.5))+100.0*T*log(10.0);

    // Use shift to compute correct proton chemical potential
    proton.mu=mup_shift+T*log(1.0/proton.g*pow(2.0*pi/proton.ms/T,1.5))+
      log_xp*T*log(10.0);
    
  } else if (nn_small) {

    // Compute chemical potential shift at a fiducial neutron density
    neutron.n=1.0e-100;

    if (use_skalt) {
      eosp_alt->calc_temp_e(neutron,proton,T,th2);
    } else {
      free_energy_density(neutron,proton,T,th2);
    }
    double mun_shift=neutron.mu-
      T*log(1.0/neutron.g*pow(2.0*pi/neutron.ms/T,1.5))+100.0*T*log(10.0);

    // Use shift to compute correct neutron chemical potential
    neutron.mu=mun_shift+T*log(1.0/neutron.g*pow(2.0*pi/neutron.ms/T,1.5))+
      log_xp*T*log(10.0);
    
  } else {
    
    if (use_skalt) {
      eosp_alt->calc_temp_e(neutron,proton,T,th_gas);
    } else {
      free_energy_density_detail(neutron,proton,T,th_gas,vdet["zn"],
				 vdet["zp"],vdet["f1"],vdet["f2"],
				 vdet["f3"],vdet["f4"],
				 vdet["g"],vdet["dgdT"]);
    }
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
      if (arg>900.0) {
	return 7;
      } else if (arg<-900.0) {
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

  y[0]=nn_tilde/nB/(1.0-Ye)-1.0;
  y[1]=np_tilde/nB/Ye-1.0;

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
    
  } else {
    
    frdm.get_nucleus(nucZ1,nucN1,*nuc_heavy);
    frdm.get_nucleus(nucZ1,nucN1-1,nuc_temp);
    Sneut[5]=-(nuclei[5].be-nuc_temp.be)*hc_mev_fm;
    frdm.get_nucleus(nucZ1-1,nucN1,nuc_temp);
    Sprot[5]=-(nuclei[5].be-nuc_temp.be)*hc_mev_fm;

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
      double zEc=(nuclei[i].Z-1.0)*fine_structure/zR*hc_mev_fm;
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
	  cout << "Nuclear free energy not finite." << endl;
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
    p_c2+=-0.6*nuclei[i].n*nuclei[i].Z*nuclei[i].Z*fine_structure/rA*
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

int eos_nuclei::nuc_matter(double nB, double Ye, double T,
			   double &log_xn, double &log_xp,
			   double &Zbar, double &Nbar, thermo &thx,
			   double &mun_full, double &mup_full,
			   int &A_min, int &A_max,
			   int &NmZ_min, int &NmZ_max,
			   map<string,double> &vdet) {
  
  neutron.n=nB*(1.0-Ye);
  proton.n=nB*Ye;
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

  // Now compute the full EOS
  if (use_skalt) {
    eosp_alt->calc_temp_e(neutron,proton,T,thx);
  } else {
    double f1, f2, f3, f4;
    free_energy_density_detail(neutron,proton,T,thx,vdet["zn"],
			       vdet["zp"],f1,f2,f3,f4,vdet["g"],
			       vdet["dgdT"]);
    if (include_detail) {
      vdet["F1"]=f1/nB*hc_mev_fm;
      vdet["F2"]=f2/nB*hc_mev_fm;
      vdet["F3"]=f3/nB*hc_mev_fm;
      vdet["F4"]=f4/nB*hc_mev_fm;
      vdet["msn"]=neutron.ms*hc_mev_fm;
      vdet["msp"]=proton.ms*hc_mev_fm;
      vdet["Un"]=neutron.nu*hc_mev_fm;
      vdet["Up"]=proton.nu*hc_mev_fm;
      vdet["dgdT"]/=hc_mev_fm;
    }
  }
  
  mun_full=neutron.mu;
  mup_full=proton.mu;
  
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
    
    if (thx2.ed-T*thx2.en<thx.ed-T*thx.en) {
      cout << "Nuclear matter preferred." << endl;
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
	  count++;
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
      ame.get_nucleus(nuclei[i].Z,nuclei[i].N-1,nuc_temp);
      Sneut[i]=-(nuclei[i].be-nuc_temp.be)*hc_mev_fm;
      ame.get_nucleus(nuclei[i].Z-1,nuclei[i].N,nuc_temp);
      Sprot[i]=-(nuclei[i].be-nuc_temp.be)*hc_mev_fm;
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

	  if (index>=count) {
	    O2SCL_ERR2("Indexing problem in ",
		       "eos_nuclei::eos_fixed_dist().",o2scl::exc_esanity);
	  }
	  
	  // Obtain heavy nucleus and neutron and proton separation
	  // energies. Ensuring that all three nuclei are in the same
	  // mass formula avoids problems with the neutron and
	  // separation energies at the boundary between the mass
	  // formulae.

	  if (index>=nuclei.size()) {
	    O2SCL_ERR("Nuclei size problem point 1.",
		      o2scl::exc_esanity);
	  }	    
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

	  } else {
	    
	    frdm.get_nucleus(Z,N,nuclei[index]);
	    frdm.get_nucleus(Z,N-1,nuc_temp);
	    Sneut[index]=-(nuclei[index].be-nuc_temp.be)*hc_mev_fm;
	    frdm.get_nucleus(Z-1,N,nuc_temp);
	    Sprot[index]=-(nuclei[index].be-nuc_temp.be)*hc_mev_fm;
	    
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
	    } else {
	      nuclei[index].g=2.0;
	    }
	  }
	  
	  index++;
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
      double zEc=(nuclei[i].Z-1.0)*fine_structure/zR*hc_mev_fm;
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
    cout << "Initial guess: " << x1[0] << " " << x1[1] << endl;
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
      cout << "  nB,Ye,T: " << nB << " " << Ye << " "
	   << T*hc_mev_fm << endl;
      cout << "  y[0],y[1],tol: " << y1[0] << " " << y1[1] << " "
	   << mh.tol_rel << endl;
      return 1;
    }
  }
  
  if (loc_verbose>2) mh.verbose=1;

  double qual_best=1.0e100;

  if (mh_ret!=0) {
    mh_ret=mh.msolve(2,x1,sn_func);
    if (mh_ret==0 && mpi_size==1) {
      cout << "Rank " << mpi_rank << " finished after initial solve." << endl;
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
	cout << "Rank " << mpi_rank << " finished solve " << k+1
	     << "/" << n_solves << " " << x1[0] << " " << x1[1] << " "
	     << qual_best << endl;
      }
      if (mh_ret==0 && mpi_size==1) {
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
	cout << "Rank " << mpi_rank << " finished bracket " << k+1
	     << "/" << n_brackets << " " << x1[0] << " " << x1[1] << endl;
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
	  if (mpi_size==1) {
	    cout << "Rank " << mpi_rank << " finished after "
		 << k << " brackets." << endl;
	  }
	} else if (k>0 && qual_last==qual) {
	  // If our quality hasn't improved, stop
	  k=10;	  
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
	  cout << "Rank " << mpi_rank << " finished minimization " << jk+1
	       << "/" << n_minimizes << " " << mret << " " << qual_best
	       << " " << x1[0] << " " << x1[1] << endl;
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
	if (mpi_size==1) {
	  cout << "Rank " << mpi_rank << " finished after "
	       << jk << " minimizer calls." << endl;
	}
      }
      
    }
    
    // ---------------------------------------------------------------
    // Survey the solution space

    if (false) {
      table3d surv;
      uniform_grid_log_end<double> ug(1.0e-8,0.01,99);
      surv.set_xy("dx",ug,"dy",ug);
      surv.line_of_names("pp1 pp2 pps ppret pm1 pm2 pms pmret");
      surv.line_of_names("mp1 mp2 mps mpret mm1 mm2 mms mmret");
      ubvector x2(2);
      for(size_t sj=0;sj<surv.get_nx();sj++) {
	for(size_t suk=0;suk<surv.get_ny();suk++) {
	  x2[0]=x1[0]+surv.get_grid_x(sj);
	  x2[1]=x1[1]+surv.get_grid_y(suk);
	  int siret=sn_func(2,x2,y1);
	  cout << "Survey: " << sj << " " << suk << " ";
	  cout << siret << " " << y1[0] << " " << y1[1] << " ";
	  if (siret==0) {
	    surv.set(sj,suk,"pp1",y1[0]);
	    surv.set(sj,suk,"pp2",y1[1]);
	    surv.set(sj,suk,"pps",fabs(y1[0])+fabs(y1[1]));
	  } else {
	    surv.set(sj,suk,"pp1",0.0);
	    surv.set(sj,suk,"pp2",0.0);
	    surv.set(sj,suk,"pps",0.0);
	  }
	  surv.set(sj,suk,"ppret",siret);
	  x2[0]=x1[0]+surv.get_grid_x(sj);
	  x2[1]=x1[1]-surv.get_grid_y(suk);
	  siret=sn_func(2,x2,y1);
	  cout << siret << " " << y1[0] << " " << y1[1] << " ";
	  if (siret==0) {
	    surv.set(sj,suk,"pm1",y1[0]);
	    surv.set(sj,suk,"pm2",y1[1]);
	    surv.set(sj,suk,"pms",fabs(y1[0])+fabs(y1[1]));
	  } else {
	    surv.set(sj,suk,"pm1",0.0);
	    surv.set(sj,suk,"pm2",0.0);
	    surv.set(sj,suk,"pms",0.0);
	  }
	  surv.set(sj,suk,"pmret",siret);
	  x2[0]=x1[0]-surv.get_grid_x(sj);
	  x2[1]=x1[1]+surv.get_grid_y(suk);
	  siret=sn_func(2,x2,y1);
	  cout << siret << " " << y1[0] << " " << y1[1] << " ";
	  if (siret==0) {
	    surv.set(sj,suk,"mp1",y1[0]);
	    surv.set(sj,suk,"mp2",y1[1]);
	    surv.set(sj,suk,"mps",fabs(y1[0])+fabs(y1[1]));
	  } else {
	    surv.set(sj,suk,"mp1",0.0);
	    surv.set(sj,suk,"mp2",0.0);
	    surv.set(sj,suk,"mps",0.0);
	  }
	  surv.set(sj,suk,"mpret",siret);
	  x2[0]=x1[0]-surv.get_grid_x(sj);
	  x2[1]=x1[1]-surv.get_grid_y(suk);
	  siret=sn_func(2,x2,y1);
	  cout << siret << " " << y1[0] << " " << y1[1] << endl;
	  if (siret==0) {
	    surv.set(sj,suk,"mm1",y1[0]);
	    surv.set(sj,suk,"mm2",y1[1]);
	    surv.set(sj,suk,"mms",fabs(y1[0])+fabs(y1[1]));
	  } else {
	    surv.set(sj,suk,"mm1",0.0);
	    surv.set(sj,suk,"mm2",0.0);
	    surv.set(sj,suk,"mms",0.0);
	  }
	  surv.set(sj,suk,"mmret",siret);
	}
      }
      double pps_max=matrix_max_value<ubmatrix,double>
	(surv.get_slice("pps"));
      for(size_t sj=0;sj<surv.get_nx();sj++) {
	for(size_t suk=0;suk<surv.get_ny();suk++) {
	  if (fabs(surv.get(sj,suk,"ppret"))>1.0e-4) {
	    surv.set(sj,suk,"pps",pps_max);
	  }
	  if (surv.get(sj,suk,"pps")>1.0) {
	    surv.set(sj,suk,"pps",1.0);
	  }
	}
      }
      double pms_max=matrix_max_value<ubmatrix,double>
	(surv.get_slice("pms"));
      for(size_t sj=0;sj<surv.get_nx();sj++) {
	for(size_t suk=0;suk<surv.get_ny();suk++) {
	  if (fabs(surv.get(sj,suk,"pmret"))>1.0e-4) {
	    surv.set(sj,suk,"pms",pms_max);
	  }
	  if (surv.get(sj,suk,"pms")>1.0) {
	    surv.set(sj,suk,"pms",1.0);
	  }
	}
      }
      double mps_max=matrix_max_value<ubmatrix,double>
	(surv.get_slice("mps"));
      for(size_t sj=0;sj<surv.get_nx();sj++) {
	for(size_t suk=0;suk<surv.get_ny();suk++) {
	  if (fabs(surv.get(sj,suk,"mpret"))>1.0e-4) {
	    surv.set(sj,suk,"mps",mps_max);
	  }
	  if (surv.get(sj,suk,"mps")>1.0) {
	    surv.set(sj,suk,"mps",1.0);
	  }
	}
      }
      double mms_max=matrix_max_value<ubmatrix,double>
	(surv.get_slice("mms"));
      for(size_t sj=0;sj<surv.get_nx();sj++) {
	for(size_t suk=0;suk<surv.get_ny();suk++) {
	  if (fabs(surv.get(sj,suk,"mmret"))>1.0e-4) {
	    surv.set(sj,suk,"mms",mms_max);
	  }
	  if (surv.get(sj,suk,"mms")>1.0) {
	    surv.set(sj,suk,"mms",1.0);
	  }
	}
      }
      cout << pps_max << " " << pms_max << " " << mps_max << " "
	   << mms_max << endl;
      double val;
      size_t suoi, suoj;
      matrix_min_index<ubmatrix,double>(surv.get_slice("pps"),
					suoi,suoj,val);
      cout << val << " " << surv.get_grid_x(suoi) << " "
	   << surv.get_grid_y(suoj) << " " 
	   << x1[0]+surv.get_grid_x(suoi) << " "
	   << x1[1]+surv.get_grid_y(suoj) << endl;
      matrix_min_index<ubmatrix,double>(surv.get_slice("pms"),
					suoi,suoj,val);
      cout << val << " " << surv.get_grid_x(suoi) << " "
	   << surv.get_grid_y(suoj) << " " 
	   << x1[0]+surv.get_grid_x(suoi) << " "
	   << x1[1]-surv.get_grid_y(suoj) << endl;
      matrix_min_index<ubmatrix,double>(surv.get_slice("mps"),
					suoi,suoj,val);
      cout << val << " " << surv.get_grid_x(suoi) << " "
	   << surv.get_grid_y(suoj) << " " 
	   << x1[0]-surv.get_grid_x(suoi) << " "
	   << x1[1]+surv.get_grid_y(suoj) << endl;
      matrix_min_index<ubmatrix,double>(surv.get_slice("mms"),
					suoi,suoj,val);
      cout << val << " " << surv.get_grid_x(suoi) << " "
	   << surv.get_grid_y(suoj) << " " 
	   << x1[0]-surv.get_grid_x(suoi) << " "
	   << x1[1]-surv.get_grid_y(suoj) << endl;
      hdf_file hf;
      hf.open_or_create("survey.o2");
      hdf_output(hf,(const table3d &)surv,"survey");
      hf.close();
      exit(-1);
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
	cout << "x1,ranges: " << x1[0] << " " << x1[1] << " "
	     << ranges[0] << " " << ranges[1] << " "
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
	  cout << "Rank " << mpi_rank << " finished random solve " << kk+1
	       << "/" << n_randoms << " " << x1[0] << " " << x1[1] << " "
	       << qual_best << endl;
	}
	if (mh_ret==0 && mpi_size==1) {
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

    if (mh_ret!=0) {
      if (loc_verbose>1) {
	cout << "Failed in eos_fixed_dist(), "
	     << "x1[0], x1[1], y1[0], y1[1]:\n  " 
	     << x1[0] << " " << x1[1] << " " << y1[0] << " "
	     << y1[1] << endl;
	cout << "qual_best: " << qual_best << endl;
	cout << "nB,Ye,T: " << nB << " " << Ye << " " << T << endl;
      }
      if (loc_verbose==8) {
	char ch;
	cin >> ch;
      } else if (loc_verbose>=9) {
	exit(-1);
      }
      return 5;
    }
    
  }

  // Perform a final call to solve_nuclei() to ensure
  // the chemical potentials are updated
  sn_func(2,x1,y1);
  
  if (mpi_size==1 && loc_verbose>1) {
    sn_func(2,x1,y1);
    cout << "Success in eos_fixed_dist(), "
	 << "x1[0], x1[1], y1[0], y1[1]:\n  " 
	 << x1[0] << " " << x1[1] << " " << y1[0] << " "
	 << y1[1] << endl;
    cout << "nB,Ye,T: " << nB << " " << Ye << " " << T << endl;
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
      en_nuc=nuclei[i].n*(log(vomega[i]/nuclei[i].n/pow(lambda,3.0))+
			  5.0/2.0+vomega_prime[i]*T/vomega[i]);
      if (!std::isfinite(fr_nuc)) {
	if (nuclei[i].n<1.0e-200) {
	  nuclei[i].n=0.0;
	  fr_nuc=0.0;
	  en_nuc=0.0;
	} else {
	  cout << "Nuclear free energy not finite." << endl;
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

  calculator calc;
  std::map<std::string,double> vars;
  calc.compile(nucleon_func.c_str());

  vars["nB"]=nB;
  vars["Ye"]=Ye;
  vars["T"]=T*o2scl_const::hc_mev_fm;
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

  int iflag=((int)(tg3_flag.get(i_nB,i_Ye,i_T)*(1.0+1.0e-12)));

  if (iflag==iflag_done) {
    double fr_old=tg3_Fint.get(i_nB,i_Ye,i_T)/hc_mev_fm*nB;
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
      cout << "  nB,Ye,T: " << nB << " " << Ye << " " << T << endl;
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
    
  tg3_log_xn.set(i_nB,i_Ye,i_T,log_xn);
  tg3_log_xp.set(i_nB,i_Ye,i_T,log_xp);
  tg3_Z.set(i_nB,i_Ye,i_T,Zbar);
  tg3_A.set(i_nB,i_Ye,i_T,Zbar+Nbar);
  tg3_flag.set(i_nB,i_Ye,i_T,loc_flag);
  tg3_Fint.set(i_nB,i_Ye,i_T,(th.ed-T*th.en)/nB*hc_mev_fm);
  tg3_Eint.set(i_nB,i_Ye,i_T,th.ed/nB*hc_mev_fm);
  tg3_Sint.set(i_nB,i_Ye,i_T,th.en/nB);

  if (loc_verbose>1) {
    cout << "store_point() output:\n  nB, Ye, T (MeV): "
	 << nB << " " << Ye << " " << T*hc_mev_fm << endl;
    cout << "  log_xn, log_xp, Z, A: "
	 << log_xn << " " << log_xp << " " << Zbar << " "
	 << Zbar+Nbar << endl;
    cout << "  flag Fint Xn Xp: " << loc_flag << " "
	 << (th.ed-T*th.en)/nB*hc_mev_fm << " "
	 << n_fraction << " " << p_fraction << endl;
  }

  //cout << "Storing: " << (th.ed-T*th.en)/nB*hc_mev_fm << endl;
  tg3_Xn.set(i_nB,i_Ye,i_T,n_fraction);
  tg3_Xp.set(i_nB,i_Ye,i_T,p_fraction);

  tg3_Xalpha.set(i_nB,i_Ye,i_T,X[0]);
  tg3_Xd.set(i_nB,i_Ye,i_T,X[1]);
  tg3_Xt.set(i_nB,i_Ye,i_T,X[2]);
  tg3_XHe3.set(i_nB,i_Ye,i_T,X[3]);
  tg3_XLi4.set(i_nB,i_Ye,i_T,X[4]);
  tg3_Xnuclei.set(i_nB,i_Ye,i_T,X[5]);

  if (loc_verbose>1) {
    cout << "  Xalpha, Xd, Xt: " << X[0] << " " << X[1] << " "
	 << X[2] << endl;
    cout << "  XHe3, XLi4, Xnuclei: " << X[3] << " " << X[4] << " "
	 << X[5] << endl;
  }
  
  if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
    tg3_A_min.set(i_nB,i_Ye,i_T,A_min);
    tg3_A_max.set(i_nB,i_Ye,i_T,A_max);
    tg3_NmZ_min.set(i_nB,i_Ye,i_T,NmZ_min);
    tg3_NmZ_max.set(i_nB,i_Ye,i_T,NmZ_max);

    if (loc_verbose>1) {
      cout << "  A_min, A_max, NmZ_min, NmZ_max: ";
      cout.precision(4);
      cout << A_min << " " << A_max << " " << NmZ_min << " "
	   << NmZ_max << endl;
      cout.precision(6);
    }
  }

  tg3_Sint.set(i_nB,i_Ye,i_T,th.en/nB);

  if (include_detail) {
    tg3_zn.set(i_nB,i_Ye,i_T,vdet["zn"]);
    tg3_zp.set(i_nB,i_Ye,i_T,vdet["zp"]);
    tg3_F1.set(i_nB,i_Ye,i_T,vdet["f1"]);
    tg3_F2.set(i_nB,i_Ye,i_T,vdet["f2"]);
    tg3_F3.set(i_nB,i_Ye,i_T,vdet["f3"]);
    tg3_F4.set(i_nB,i_Ye,i_T,vdet["f4"]);
    tg3_Un.set(i_nB,i_Ye,i_T,vdet["Un"]);
    tg3_Up.set(i_nB,i_Ye,i_T,vdet["Up"]);
    tg3_msn.set(i_nB,i_Ye,i_T,vdet["msn"]);
    tg3_msp.set(i_nB,i_Ye,i_T,vdet["msp"]);
    tg3_g.set(i_nB,i_Ye,i_T,vdet["g"]);
    tg3_dgdT.set(i_nB,i_Ye,i_T,vdet["dgdT"]);
  }
  
  // AWS 8/4/2020: This section is commented out because the code does
  // not correctly analytically compute mun, mup, Eint and Pint
  if (false && derivs_computed) {
    
    tg3_Pint.set(i_nB,i_Ye,i_T,th.pr*hc_mev_fm);

    if (loc_verbose>1) {
      cout << "  Eint, Pint, Sint: " << th.ed/nB*hc_mev_fm << " "
	   << th.pr*hc_mev_fm << " " << th.en/nB << endl;
    }
    tg3_mun.set(i_nB,i_Ye,i_T,mun_full*hc_mev_fm);
    tg3_mup.set(i_nB,i_Ye,i_T,mup_full*hc_mev_fm);

    if (loc_verbose>1) {
      cout << "  mun, mup: " << mun_full*hc_mev_fm << " "
	   << mup_full*hc_mev_fm << endl;
    }
    
    if (with_leptons_loaded) {
      
      electron.n=nB*Ye;
      electron.mu=electron.m;
      relf.pair_density(electron,T);
      photon.massless_calc(T);
      

      tg3_F.set(i_nB,i_Ye,i_T,(fr+electron.ed+photon.ed-
			       T*(electron.en+photon.en))/nB*hc_mev_fm);
      tg3_E.set(i_nB,i_Ye,i_T,(th.ed+electron.ed+photon.ed)/nB*hc_mev_fm);
      tg3_P.set(i_nB,i_Ye,i_T,(th.pr+electron.pr+photon.pr)*hc_mev_fm);
      tg3_S.set(i_nB,i_Ye,i_T,(th.en+electron.en+photon.en)/nB);
      tg3_mue.set(i_nB,i_Ye,i_T,electron.mu*hc_mev_fm);
      
    }
    
  }
  
  if (loc_verbose>8) {
    exit(-1);
  }

  return 0;
}

int eos_nuclei::select_high_T_cl(std::vector<std::string> &sv,
				 bool itive_com) {
  if (sv.size()<2) {
    cerr << "Command \"select-high-T\" needs an integer argument."
	 << endl;
    return 1;
  }
  return select_high_T(o2scl::stoi(sv[1]));
}

int eos_nuclei::select_high_T(int option) {
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
    tm.test_abs(skyrme_ext.msom,0.8,0.2,"msom");
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
	  
	} else {
	  
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
	  
	}
	
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
	  } else {
	    nuc.g=2.0;
	  }
	}

	t.line_of_data(9,line);
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
  
  // Ensure that multiple MPI ranks aren't reading from the
  // filesystem at the same time
  int tag=0, buffer=0;
  if (mpi_size>1 && mpi_rank>=1) {
    MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,
	     tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
#endif
  
  cout << "Function write_results(): rank " << mpi_rank
       << " writing file " << fname << endl;
  
  hdf_file hf;
  
  wordexp_single_file(fname);
  
  hf.open_or_create(fname);
  
  hf.set_szt("n_nB",n_nB2);
  hf.set_szt("n_Ye",n_Ye2);
  hf.set_szt("n_T",n_T2);
  hf.setd_vec("nB_grid",nB_grid2);
  hf.setd_vec("Ye_grid",Ye_grid2);
  hf.setd_vec("T_grid",T_grid2);

  hdf_output(hf,tg3_log_xn,"log_xn");
  hdf_output(hf,tg3_log_xp,"log_xp");
  hdf_output(hf,tg3_Z,"Z");
  hdf_output(hf,tg3_A,"A");
  hdf_output(hf,tg3_flag,"flag");
  hdf_output(hf,tg3_Fint,"Fint");
  hdf_output(hf,tg3_Sint,"Sint");
  hdf_output(hf,tg3_Eint,"Eint");
  
  hdf_output(hf,tg3_Xn,"Xn");
  hdf_output(hf,tg3_Xp,"Xp");
  hdf_output(hf,tg3_Xalpha,"Xalpha");
  hdf_output(hf,tg3_Xnuclei,"Xnuclei");
  hdf_output(hf,tg3_Xd,"Xd");
  hdf_output(hf,tg3_Xt,"Xt");
  hdf_output(hf,tg3_XHe3,"XHe3");
  hdf_output(hf,tg3_XLi4,"XLi4");

  if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
    hdf_output(hf,tg3_A_min,"A_min");
    hdf_output(hf,tg3_A_max,"A_max");
    hdf_output(hf,tg3_NmZ_min,"NmZ_min");
    hdf_output(hf,tg3_NmZ_max,"NmZ_max");
  }
  
  hf.seti("baryons_only",1);
  hf.seti("alg_mode",alg_mode);
  if (derivs_computed) {
    hf.seti("derivs_computed",1);
  } else {
    hf.seti("derivs_computed",0);
  }
  if (with_leptons_loaded) {
    hf.seti("with_leptons",1);
  } else {
    hf.seti("with_leptons",0);
  }
  // The parent class has an include_muons bool but this
  // child class doesn't support muons yet
  hf.seti("include_muons",0);

  hf.setd("m_neut",neutron.m*o2scl_const::hc_mev_fm);
  hf.setd("m_prot",proton.m*o2scl_const::hc_mev_fm);
  hf.setd("hc",o2scl_const::hc_mev_fm);
  hf.setd("alpha_em",o2scl_const::fine_structure);
  
  if (derivs_computed) {
    hdf_output(hf,tg3_Pint,"Pint");
    hdf_output(hf,tg3_mun,"mun");
    hdf_output(hf,tg3_mup,"mup");
    if (with_leptons_loaded) {
      hdf_output(hf,tg3_F,"F");
      hdf_output(hf,tg3_E,"E");
      hdf_output(hf,tg3_P,"P");
      hdf_output(hf,tg3_S,"S");
      hdf_output(hf,tg3_mup,"mue");
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
    hdf_output(hf,tg3_zn,"zn");
    hdf_output(hf,tg3_zp,"zp");
    hdf_output(hf,tg3_msn,"msn");
    hdf_output(hf,tg3_msp,"msp");
    hdf_output(hf,tg3_F1,"F1");
    hdf_output(hf,tg3_F2,"F2");
    hdf_output(hf,tg3_F3,"F3");
    hdf_output(hf,tg3_F4,"F4");
    hdf_output(hf,tg3_Un,"Un");
    hdf_output(hf,tg3_Up,"Up");
    hdf_output(hf,tg3_g,"g");
    hdf_output(hf,tg3_dgdT,"dgdT");
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
  hf.sets_vec("oth_names",oth_names);
  hf.sets_vec("oth_units",oth_units);
  
  hf.close();
  
  cout << "Function write_results(): rank " << mpi_rank
       << " done writing file." << endl;

#ifndef NO_MPI
  // Send a message to the next MPI rank
  if (mpi_size>1 && mpi_rank<mpi_size-1) {
    MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,
	     tag,MPI_COMM_WORLD);
  }
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

  // nB, Ye, T grid

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
  
  if (hf.find_object_by_name("log_xn",type)!=0 || type!="tensor_grid") {
    O2SCL_ERR("Couldn't find tensor log_xn in file.",
	      o2scl::exc_enotfound);
  }

  // Flags

  int itmp;
  if (verbose>2) cout << "Reading baryons_only." << endl;
  hf.geti_def("baryons_only",1,itmp);
  if (itmp==1) baryons_only_loaded=true;
  else baryons_only_loaded=false;
  
  if (verbose>2) cout << "Reading with_leptons." << endl;
  hf.geti_def("with_leptons",0,itmp);
  if (itmp==1) with_leptons_loaded=true;
  else with_leptons_loaded=false;
  
  if (verbose>2) cout << "Reading derivs_computed." << endl;
  hf.geti_def("derivs_computed",0,itmp);
  if (itmp==1) derivs_computed=true;
  else derivs_computed=false;
  
  if (verbose>2) cout << "Reading include_detail." << endl;
  hf.geti_def("detail",0,itmp);
  if (itmp==1) include_detail=true;
  else include_detail=false;
  
  include_muons=false;

  if (verbose>2) cout << "Reading alg_mode." << endl;
  hf.geti_def("alg_mode",4,alg_mode);

  // Main data
  
  if (verbose>2) cout << "Reading log_xn." << endl;
  hdf_input(hf,tg3_log_xn,"log_xn");
  if (verbose>2) cout << "Reading log_xp." << endl;
  hdf_input(hf,tg3_log_xp,"log_xp");
  if (verbose>2) cout << "Reading Z." << endl;
  hdf_input(hf,tg3_Z,"Z");
  if (verbose>2) cout << "Reading A." << endl;
  hdf_input(hf,tg3_A,"A");
  if (verbose>2) cout << "Reading flag." << endl;
  hdf_input(hf,tg3_flag,"flag");
  if (verbose>2) cout << "Reading Fint." << endl;
  hdf_input(hf,tg3_Fint,"Fint");
  
  if (hf.find_object_by_name("Sint",type)==0 && type=="tensor_grid") {
    if (verbose>2) cout << "Reading Sint." << endl;
    hdf_input(hf,tg3_Sint,"Sint");
  }

  // Fix old tables which don't have Eint
  if (hf.find_object_by_name("Eint",type)==0 && type=="tensor_grid") {    
    if (verbose>2) cout << "Reading Eint." << endl;
    hdf_input(hf,tg3_Eint,"Eint");
  } else if (tg3_Sint.total_size()>0) {
    cout << "Fixing file with missing Eint." << endl;
    vector<vector<double> > grid={nB_grid2,Ye_grid2,
				  T_grid2};
    vector<size_t> sz={n_nB2,n_Ye2,n_T2};
    tg3_Eint.resize(3,sz);
    tg3_Eint.set_grid(grid);
    for(size_t i=0;i<tg3_Eint.total_size();i++) {
      //std::cout << "J: " << i << " " << tg3_Eint.total_size() << endl;
      size_t ix[3];
      tg3_Eint.unpack_index(i,ix);
      //std::cout << "H: " << i << " " << ix[0] << std::endl;
      
      double T_MeV=T_grid2[ix[2]];
      
      double Eint_val=tg3_Fint.get_data()[i]+
	T_MeV*tg3_Sint.get_data()[i];
      tg3_Eint.get(ix[0],ix[1],ix[2])=Eint_val;
    }
  }
  
  if (verbose>2) cout << "Reading Xn." << endl;
  hdf_input(hf,tg3_Xn,"Xn");
  if (verbose>2) cout << "Reading Xp." << endl;
  hdf_input(hf,tg3_Xp,"Xp");
  if (verbose>2) cout << "Reading Xalpha." << endl;
  hdf_input(hf,tg3_Xalpha,"Xalpha");
  if (verbose>2) cout << "Reading Xnuclei." << endl;
  hdf_input(hf,tg3_Xnuclei,"Xnuclei");
  if (verbose>2) cout << "Reading Xd." << endl;
  hdf_input(hf,tg3_Xd,"Xd");
  if (verbose>2) cout << "Reading Xt." << endl;
  hdf_input(hf,tg3_Xt,"Xt");
  if (verbose>2) cout << "Reading XHe3." << endl;
  hdf_input(hf,tg3_XHe3,"XHe3");
  if (verbose>2) cout << "Reading XLi4." << endl;
  hdf_input(hf,tg3_XLi4,"XLi4");

  // Nuclear distribution
  
  if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
    if (hf.find_object_by_name("A_min",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor A_min in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading A_min." << endl;
    hdf_input(hf,tg3_A_min,"A_min");
    if (hf.find_object_by_name("A_max",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor A_max in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading A_max." << endl;
    hdf_input(hf,tg3_A_max,"A_max");
    if (hf.find_object_by_name("NmZ_min",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor NmZ_min in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading NmZ_min." << endl;
    hdf_input(hf,tg3_NmZ_min,"NmZ_min");
    if (hf.find_object_by_name("NmZ_max",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor NmZ_max in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading NmZ_max." << endl;
    hdf_input(hf,tg3_NmZ_max,"NmZ_max");
  }

  // Derivatives of free energy

  if (with_leptons_loaded && !derivs_computed) {
    O2SCL_ERR2("File says leptons are present but derivatives are not ",
	       "computed.",o2scl::exc_esanity);
  }
  
  if (derivs_computed) {
    if (hf.find_object_by_name("Pint",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor Pint in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading Pint." << endl;
    hdf_input(hf,tg3_Pint,"Pint");
    if (hf.find_object_by_name("mun",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor mun in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading mun." << endl;
    hdf_input(hf,tg3_mun,"mun");
    if (hf.find_object_by_name("mup",type)!=0 || type!="tensor_grid") {
      O2SCL_ERR("Couldn't find tensor mup in file.",
		o2scl::exc_enotfound);
    }
    if (verbose>2) cout << "Reading mup." << endl;
    hdf_input(hf,tg3_mup,"mup");
    
    if (with_leptons_loaded) {
      if (verbose>2) cout << "Reading F." << endl;
      hdf_input(hf,tg3_F,"F");
      if (verbose>2) cout << "Reading E." << endl;
      hdf_input(hf,tg3_E,"E");
      if (verbose>2) cout << "Reading P." << endl;
      hdf_input(hf,tg3_P,"P");
      if (verbose>2) cout << "Reading S." << endl;
      hdf_input(hf,tg3_S,"S");
      if (verbose>2) cout << "Reading mue." << endl;
      hdf_input(hf,tg3_mue,"mue");
    } else {
      tg3_F.clear();
      tg3_E.clear();
      tg3_P.clear();
      tg3_S.clear();
      tg3_mue.clear();
    }
  } else {
    tg3_Pint.clear();
    tg3_mun.clear();
    tg3_mup.clear();
    tg3_F.clear();
    tg3_E.clear();
    tg3_P.clear();
    tg3_S.clear();
    tg3_mue.clear();
  }

  if (include_detail) {
    if (verbose>2) cout << "Reading zn." << endl;
    hdf_input(hf,tg3_zn,"zn");
    if (verbose>2) cout << "Reading zp." << endl;
    hdf_input(hf,tg3_zp,"zp");
    if (verbose>2) cout << "Reading msn." << endl;
    hdf_input(hf,tg3_msn,"msn");
    if (verbose>2) cout << "Reading msp." << endl;
    hdf_input(hf,tg3_msp,"msp");
    if (verbose>2) cout << "Reading F1." << endl;
    hdf_input(hf,tg3_F1,"F1");
    if (verbose>2) cout << "Reading F2." << endl;
    hdf_input(hf,tg3_F2,"F2");
    if (verbose>2) cout << "Reading F3." << endl;
    hdf_input(hf,tg3_F3,"F3");
    if (verbose>2) cout << "Reading F4." << endl;
    hdf_input(hf,tg3_F4,"F4");
    if (verbose>2) cout << "Reading Un." << endl;
    hdf_input(hf,tg3_Un,"Un");
    if (verbose>2) cout << "Reading Up." << endl;
    hdf_input(hf,tg3_Up,"Up");
    if (verbose>2) cout << "Reading g." << endl;
    hdf_input(hf,tg3_g,"g");
    if (verbose>2) cout << "Reading dgdT." << endl;
    hdf_input(hf,tg3_dgdT,"dgdT");
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
    
    flag=((int)(tg3_flag.get(inB,iYe,iT)+1.0e-10));
    
    cout << "Found grid point (inB,iYe,iT)=(" << inB << ","
	 << iYe << "," << iT << ")\n\t(nB,Ye,T)=("
	 << nB << "," << Ye << "," << T*hc_mev_fm << "), flag="
	 << flag << "." << endl;
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
    cout << log_xn << " " << log_xp << " " << A_min << " "
	 << A_max << " " << NmZ_min << " " << NmZ_max << endl;
    guess_provided=true;
    
  } else if ((alg_mode==0 || alg_mode==1) && sv.size()>=8) {
    
    cout << "Function point_nuclei(): "
	 << "Reading guess (log_xn,log_xp,N,Z) from command line." << endl;
    log_xn=o2scl::function_to_double(sv[4]);
    log_xp=o2scl::function_to_double(sv[5]);
    nuc_Z1=((size_t)(o2scl::function_to_double(sv[6])*(1.0+1.0e-12)));
    nuc_N1=((size_t)(o2scl::function_to_double(sv[7])*(1.0+1.0e-12)));
    cout << log_xn << " "
	 << log_xp << " " << nuc_Z1 << " " << nuc_N1 << endl;
    guess_provided=true;
    
  } else if (loaded) {
    
    // If the flag is 'guess' or 'in progress with guess'
    if (flag>=5 || flag==-5) {
      
      if (alg_mode==0 || alg_mode==1) {
	cout << "Function point_nuclei(): "
	     << "Reading guess (N,Z) from current table." << endl;
	log_xn=tg3_log_xn.get(inB,iYe,iT);
	log_xp=tg3_log_xp.get(inB,iYe,iT);
	nuc_Z1=((size_t)(tg3_Z.get(inB,iYe,iT)));
	nuc_N1=((size_t)(tg3_A.get(inB,iYe,iT)))-nuc_Z1;
	cout << log_xn << " "
	     << log_xp << " " << nuc_Z1 << " " << nuc_N1 << endl;
	guess_provided=true;
      } else {
	cout << "Function point_nuclei(): "
	     << "Reading guess (log_xn,log_xp,A_min,A_max,NmZ_min,NmZ_max) "
	     << "from current table." << endl;
	log_xn=tg3_log_xn.get(inB,iYe,iT);
	log_xp=tg3_log_xp.get(inB,iYe,iT);
	A_min=((int)(tg3_A_min.get(inB,iYe,iT)));
	A_max=((int)(tg3_A_max.get(inB,iYe,iT)));
	NmZ_min=((int)(tg3_NmZ_min.get(inB,iYe,iT)));
	NmZ_max=((int)(tg3_NmZ_max.get(inB,iYe,iT)));
	cout << log_xn << " " << log_xp << " " << A_min << " "
	     << A_max << " " << NmZ_min << " " << NmZ_max << endl;
	guess_provided=true;
      }
      
    } else if (inB>0 && tg3_flag.get(inB-1,iYe,iT)>9.9 && alg_mode==4) {

      // Otherwise if there is a good guess from the adjacent point
      // at lower baryon density
      
      cout << "Reading guess from lower baryon density point." << endl;
      log_xn=tg3_log_xn.get(inB-1,iYe,iT);
      log_xp=tg3_log_xp.get(inB-1,iYe,iT);
      A_min=((int)(tg3_A_min.get(inB-1,iYe,iT)));
      A_max=((int)(tg3_A_max.get(inB-1,iYe,iT)));
      NmZ_min=((int)(tg3_NmZ_min.get(inB-1,iYe,iT)));
      NmZ_max=((int)(tg3_NmZ_max.get(inB-1,iYe,iT)));
      cout << "  log_xn,log_xp,Z,A: "
	   << log_xn << " " << log_xp << " ";
      cout << tg3_Z.get(inB-1,iYe,iT) << " ";
      cout << tg3_A.get(inB-1,iYe,iT) << endl;
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
	tg3_flag.get(inB-1,iYe,iT)>9.9 &&
	tg3_flag.get(inB+1,iYe,iT)>9.9) {
      cout << "Automatically setting ranges." << endl;
      fd_rand_ranges.resize(4);
      fd_rand_ranges[0]=tg3_log_xn.get(inB-1,iYe,iT);
      fd_rand_ranges[1]=tg3_log_xn.get(inB+1,iYe,iT);
      fd_rand_ranges[2]=tg3_log_xp.get(inB-1,iYe,iT);
      fd_rand_ranges[3]=tg3_log_xp.get(inB+1,iYe,iT);
      cout << "  ";
      vector_out(cout,fd_rand_ranges,true);
    }
    
    thermo thx;
    double mun_full, mup_full;
    
    double Zbar, Nbar;
    map<string,double> vdet;
    if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
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

    if (ret!=0) {
      cout << "Point failed." << endl;
    } else if (loaded) {
      cout << "Point succeeded. Storing." << endl;

      ubvector X;
      compute_X(nB,X);

      store_point(inB,iYe,iT,nB,Ye,T,thx,log_xn,log_xp,Zbar,Nbar,
		  mun_full,mup_full,X,A_min,A_max,NmZ_min,NmZ_max,
		  10.0,vdet);
    } else {
      cout << "Point succeeded." << endl;
    }

    if (ret==0) {
      cout << "log_xn: " << log_xn << endl;
      cout << "log_xp: " << log_xp << endl;
      cout << "fint: " << thx.ed-T*thx.en << endl;
      cout << "Fint: " << (thx.ed-T*thx.en)/nB*hc_mev_fm << " MeV" << endl;
      cout << "A: " << Zbar+Nbar << endl;
      cout << "Z: " << Zbar << endl;
      if (include_detail) {
	cout << "zn: " << vdet["zn"] << endl;
	cout << "zp: " << vdet["zp"] << endl;
	cout << "F1: " << vdet["F1"] << " MeV" << endl;
	cout << "F2: " << vdet["F2"] << " MeV" << endl;
	cout << "F3: " << vdet["F3"] << " MeV" << endl;
	cout << "F3: " << vdet["F3"] << " MeV" << endl;
	cout << "msn: " << vdet["msn"] << " MeV" << endl;
	cout << "msp: " << vdet["msp"] << " MeV" << endl;
	cout << "Un: " << vdet["Un"] << " MeV" << endl;
	cout << "Up: " << vdet["Up"] << " MeV" << endl;
	cout << "g: " << vdet["g"] << endl;
	cout << "dgdT: " << vdet["dgdT"] << " 1/MeV" << endl;
      }
    }      

    if (fd_rand_ranges.size()>0) fd_rand_ranges.clear();
    
  }

  // (At this point in the code, the integer 'flag' is the old
  // flag, not the new flag, which is stored in tg3_flag)

  if (flag>=10 && recompute==false) {
    cout << "Point already computed and recompute is false." << endl;
  }

  // -------------------------------------------------------------------
  // Print out the results when tables have been loaded and if either
  // the computation worked (ret==0) or if the old flag was marked
  // done (flag==10).
  
  if ((ret==0 || flag==10) && loaded) {
    cout << "Results stored in table:" << endl;
    cout << "nB,Ye,T: " << nB << " " << Ye << " " << T*hc_mev_fm << endl;
    cout << "flag: " << tg3_flag.get(inB,iYe,iT) << endl;
    cout << "log_xn: " << tg3_log_xn.get(inB,iYe,iT) << endl;
    cout << "log_xp: " << tg3_log_xp.get(inB,iYe,iT) << endl;
    cout << "Z: " << tg3_Z.get(inB,iYe,iT) << endl;
    cout << "A: " << tg3_A.get(inB,iYe,iT) << endl;
    cout << "Fint: " << tg3_Fint.get(inB,iYe,iT) << " MeV" << endl;
    cout << "Sint: " << tg3_Sint.get(inB,iYe,iT) << endl;
    cout << "Eint: " << tg3_Eint.get(inB,iYe,iT) << " MeV" << endl;
    cout << "A_min: " << tg3_A_min.get(inB,iYe,iT) << endl;
    cout << "A_max: " << tg3_A_max.get(inB,iYe,iT) << endl;
    cout << "NmZ_min: " << tg3_NmZ_min.get(inB,iYe,iT) << endl;
    cout << "NmZ_max: " << tg3_NmZ_max.get(inB,iYe,iT) << endl;
    cout << "Xn: " << tg3_Xn.get(inB,iYe,iT) << endl;
    cout << "Xp: " << tg3_Xp.get(inB,iYe,iT) << endl;
    cout << "Xd: " << tg3_Xd.get(inB,iYe,iT) << endl;
    cout << "Xt: " << tg3_Xt.get(inB,iYe,iT) << endl;
    cout << "XHe3: " << tg3_XHe3.get(inB,iYe,iT) << endl;
    cout << "XLi4: " << tg3_XLi4.get(inB,iYe,iT) << endl;
    cout << "Xalpha: " << tg3_Xalpha.get(inB,iYe,iT) << endl;
    cout << "Xnuclei: " << tg3_Xnuclei.get(inB,iYe,iT) << endl;
    if (include_detail) {
      cout << "zn: " << tg3_zn.get(inB,iYe,iT) << endl;
      cout << "zp: " << tg3_zp.get(inB,iYe,iT) << endl;
      cout << "F1: " << tg3_F1.get(inB,iYe,iT) << " MeV" << endl;
      cout << "F2: " << tg3_F2.get(inB,iYe,iT) << " MeV" << endl;
      cout << "F3: " << tg3_F3.get(inB,iYe,iT) << " MeV" << endl;
      cout << "F3: " << tg3_F3.get(inB,iYe,iT) << " MeV" << endl;
      cout << "msn: " << tg3_msn.get(inB,iYe,iT) << " MeV" << endl;
      cout << "msp: " << tg3_msp.get(inB,iYe,iT) << " MeV" << endl;
      cout << "Un: " << tg3_Un.get(inB,iYe,iT) << " MeV" << endl;
      cout << "Up: " << tg3_Up.get(inB,iYe,iT) << " MeV" << endl;
      cout << "g: " << tg3_g.get(inB,iYe,iT) << endl;
      cout << "dgdT: " << tg3_dgdT.get(inB,iYe,iT) << " 1/MeV" << endl;
    }
  }

  // -------------------------------------------------------------------
  // Output the nuclear distribution
  
  if (alg_mode>=2 && show_all_nuclei) {
    
    cout << "Writing distribution to dist.o2." << endl;
    
    table3d t3d;
    A_max=((int)(tg3_A_max.get(inB,iYe,iT)+1.0e-10));
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
  with_leptons_loaded=false;
  
  for(size_t iT=iT_start;iT<=iT_end;iT++) {
    
    for(size_t iYe=iYe_start;iYe<=iYe_end;iYe++) {
  
      log_xn=tg3_log_xn.get(inB_start,iYe,iT);
      log_xp=tg3_log_xp.get(inB_start,iYe,iT);
      A_min=tg3_A_min.get(inB_start,iYe,iT);
      A_max=tg3_A_max.get(inB_start,iYe,iT);
      NmZ_min=tg3_NmZ_min.get(inB_start,iYe,iT);
      NmZ_max=tg3_NmZ_max.get(inB_start,iYe,iT);

      bool no_nuclei=false;

      vector<double> Zarr, Narr;
      
      for(size_t inB=inB_start;inB<=inB_end;inB++) {
	
	double nB=nB_grid2[inB];
	double Ye=Ye_grid2[iYe];
	double T=T_grid2[iT]/hc_mev_fm;
	
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

	  if (tg3_A.get(inB,iYe,iT)>0.0) {
	    cout << "before: " << tg3_Z.get(inB,iYe,iT) << " "
		 << tg3_A.get(inB,iYe,iT) << " "
		 << tg3_A.get(inB,iYe,iT)/
	      tg3_Z.get(inB,iYe,iT) << endl;
	  } else {
	    cout << "before: " << tg3_Z.get(inB,iYe,iT) << " "
		 << tg3_A.get(inB,iYe,iT) << " " << 0.0 << endl;
	  }
	  
	  store_point(inB,iYe,iT,nB,Ye,T,thx,log_xn,log_xp,Zbar,Nbar,
		      mun_full,mup_full,X,A_min,A_max,NmZ_min,
		      NmZ_max,10.0,vdet);
	  
	  if (tg3_A.get(inB,iYe,iT)>0.0) {
	    cout << "after: " << tg3_log_xn.get(inB,iYe,iT) << " "
		 << tg3_log_xp.get(inB,iYe,iT) << " "
		 << tg3_Xnuclei.get(inB,iYe,iT) << " "
		 << tg3_Z.get(inB,iYe,iT) << " "
		 << tg3_A.get(inB,iYe,iT) << " "
		 << tg3_A.get(inB,iYe,iT)/
	      tg3_Z.get(inB,iYe,iT) << endl;
	  } else {
	    cout << "after: " << tg3_log_xn.get(inB,iYe,iT) << " "
		 << tg3_log_xp.get(inB,iYe,iT) << " "
		 << tg3_Xnuclei.get(inB,iYe,iT) << " "
		 << tg3_Z.get(inB,iYe,iT) << " "
		 << tg3_A.get(inB,iYe,iT) << " " << 0.0 << endl;
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
  with_leptons_loaded=false;
  bool no_nuclei;
  
  for(size_t iT=0;iT<n_T2;iT++) {

    bool found_point=false;
    
    for(size_t iYe=0;iYe<n_Ye2;iYe++) {
  
      for(size_t inB=inB_start;inB<=inB_end;inB++) {
	
	double nB=nB_grid2[inB];
	double Ye=Ye_grid2[iYe];
	double T=T_grid2[iT]/hc_mev_fm;

	if (tg3_flag.get(inB,iYe,iT)<9.9) {
	  found_point=true;
	  
	  cout << "Found uncomputed point at (nB,Ye,T):\n"
	       << "\t" << nB << " " << Ye << " "
	       << T*hc_mev_fm << endl;
	  
	  no_nuclei=false;
	  
	  if (inB>20 && no_nuclei==false) {
	    for(size_t j=1;j<20 && no_nuclei==false;j++) {
	      cout << "j: " << j << endl;
	      if (tg3_flag.get(inB-j,iYe,iT)>9.9 &&
		  tg3_Z.get(inB-j,iYe,iT)<1.0e-4 &&
		  tg3_A.get(inB-j,iYe,iT)<1.0e-4) {
		cout << "Found nuclear matter at a lower density."
		     << endl;
		no_nuclei=true;
	      }
	    }
	  }
	  
	  if (no_nuclei==false && inB>2) {
	    double Z_prev=tg3_Z.get(inB-1,iYe,iT);
	    double Z_prev2=tg3_Z.get(inB-2,iYe,iT);
	    if (Z_prev>0 && Z_prev2>0) {
	      double A_prev=tg3_A.get(inB-1,iYe,iT);
	      double A_prev2=tg3_A.get(inB-2,iYe,iT);
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
	  log_xn=tg3_log_xn.get(inB-1,iYe,iT);
	  log_xp=tg3_log_xp.get(inB-1,iYe,iT);
	  A_min=tg3_A_min.get(inB-1,iYe,iT);
	  A_max=tg3_A_max.get(inB-1,iYe,iT);
	  NmZ_min=tg3_NmZ_min.get(inB-1,iYe,iT);
	  NmZ_max=tg3_NmZ_max.get(inB-1,iYe,iT);
	  
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
	    
	    if (tg3_A.get(inB,iYe,iT)>0.0) {
	      cout << "Success (log_xn,log_xp,Xnuc,Z,N,N/Z): "
		   << tg3_log_xn.get(inB,iYe,iT) << " "
		   << tg3_log_xp.get(inB,iYe,iT) << " "
		   << tg3_Xnuclei.get(inB,iYe,iT) << " "
		   << tg3_Z.get(inB,iYe,iT) << " "
		   << tg3_A.get(inB,iYe,iT)-tg3_Z.get(inB,iYe,iT) << " "
		   << (tg3_A.get(inB,iYe,iT)-tg3_Z.get(inB,iYe,iT))/
		tg3_Z.get(inB,iYe,iT) << endl;
	    } else {
	      cout << "Success (log_xn,log_xp,Xnuc,Z,N): "
		   << tg3_log_xn.get(inB,iYe,iT) << " "
		   << tg3_log_xp.get(inB,iYe,iT) << " "
		   << tg3_Xnuclei.get(inB,iYe,iT) << " "
		   << tg3_Z.get(inB,iYe,iT) << " "
		   << tg3_A.get(inB,iYe,iT) << endl;
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
  with_leptons_loaded=false;
  bool no_nuclei;

  size_t inB_lo=0, inB_hi=0;
  if (mode=="random_lg" || mode=="all_lg") {
    inB_lo=vector_lookup(nB_grid2,0.01);
    inB_hi=vector_lookup(nB_grid2,0.15);
  }

  if (mode=="all_lg") {
    n_tot=(inB_hi-inB_lo+1)*n_Ye2*n_T2;
  }

  bool verify_success=true;
  
  for(size_t j=0;j<n_tot;j++) {
    
    size_t inB=0, iYe=0, iT=0;

    if (mode=="all") {
      vector<size_t> vix;
      tg3_A.unpack_index(j,vix);
      inB=vix[0];
      iYe=vix[1];
      iT=vix[2];
    } else if (mode=="all_lg") {
      size_t n_slice=n_Ye2*n_T2;
      inB=j/n_slice;
      size_t j2=j-n_slice*inB;
      iYe=j2/n_Ye2;
      iT=j2-n_T2*iYe;
      cout << j << " " << inB << " " << iYe << " " << iT << endl;
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

    // Get initial guess 
    log_xn=tg3_log_xn.get(inB,iYe,iT);
    log_xp=tg3_log_xp.get(inB,iYe,iT);
    A_min=tg3_A_min.get(inB,iYe,iT);
    A_max=tg3_A_max.get(inB,iYe,iT);
    NmZ_min=tg3_NmZ_min.get(inB,iYe,iT);
    NmZ_max=tg3_NmZ_max.get(inB,iYe,iT);
    no_nuclei=false;
    if (tg3_A.get(inB,iYe,iT)<1.0e-6 &&
	tg3_Z.get(inB,iYe,iT)<1.0e-6) {
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
    if (ret!=0) {
      cout << "Verification failed. Computing." << endl;
      ret=eos_vary_dist(nB,Ye,T,log_xn,log_xp,Zbar,Nbar,
			  thx,mun_full,mup_full,
			  A_min,A_max,NmZ_min,NmZ_max,
			  vdet,true,no_nuclei);
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
  cout << "with_leptons: " << with_leptons_loaded << endl;
  cout << "alg_mode: " << alg_mode << endl;
  cout << endl;
  
  const vector<double> &data=tg3_flag.get_data();
  vector<size_t> flags(22);
  
  size_t nb_frac_count=0, S_neg_count=0, Fint_count=0, F_count=0;

  size_t ti_int_count=0, ti_count=0;
  
  for(size_t i=0;i<data.size();i++) {
    
    int iflag=((int)(data[i]*(1.0+1.0e-12)));
    if (iflag<-10 || iflag>10) {
      iflag=11;
    }
    flags[iflag+10]++;
    
    if (iflag==10) {
      
      size_t ix[3];
      tg3_Xn.unpack_index(i,ix);
      
      double nB=nB_grid2[ix[0]];
      double Ye=Ye_grid2[ix[1]];
      double T_MeV=T_grid2[ix[2]];

      // Check that X's add up to 1
      if (nb_frac_count<1000) {
	double check_X=tg3_Xn.get_data()[i]+tg3_Xp.get_data()[i]+
	  tg3_Xalpha.get_data()[i]+tg3_Xnuclei.get_data()[i]+
	  tg3_Xd.get_data()[i]+tg3_Xt.get_data()[i]+
	  tg3_XHe3.get_data()[i]+tg3_XLi4.get_data()[i];
	
	if (fabs(check_X-1.0)>1.0e-6) {
	  cout << "Nuclear fractions do not add up to 1 "
	       << "(i,nB,Ye,T,X_total):\n  "
	       << i << " " << nB << " " << Ye << " "
	       << T_MeV << " " << check_X << endl;
	  nb_frac_count++;
	  if (nb_frac_count==1000) {
	    cout << "Further nuclear fractions warnings suppressed."
		 << endl;
	  }
	}
      }

      // Check that Fint=Eint-T*Sint
      if (tg3_Eint.total_size()>0 && Fint_count<1000) {

	double scale=tg3_Fint.get_data()[i];
	if (scale<10.0) scale=10.0;
	double Fint_check=(tg3_Fint.get_data()[i]-
			   tg3_Eint.get_data()[i]+
			   T_MeV*tg3_Sint.get_data()[i])/scale;
	if (fabs(Fint_check)>1.0e-9) {
	  cout << "Fint doesn't match Eint-T*Sint:"
	       << "(i,nB,Ye,T,Fint_check):\n  "
	       << i << " " << nB << " " << Ye << " " << T_MeV
	       << " " << Fint_check << endl;
	  cout << "  (Fint,Eint,T,Sint,Eint-T*Sint): "
	       << tg3_Fint.get_data()[i] << " ";
	  cout << tg3_Eint.get_data()[i] << " ";
	  cout << T_MeV << " " << tg3_Sint.get_data()[i] << " ";
	  cout << tg3_Eint.get_data()[i]-T_MeV*tg3_Sint.get_data()[i] << endl;
	  Fint_count++;
	  if (Fint_count==1000) {
	    cout << "Further Fint warnings suppressed." << endl;
	  }
	}
      }

      // Check that Eint*nB=-Pint+T*Sint*nB+nn*mun+np*mup
      if (derivs_computed && ti_int_count<1000) {
	double nn=nB*(1.0-Ye);
	double np=nB*Ye;
	double scale=tg3_Eint.get_data()[i];
	if (scale<100.0) scale=100.0;
	double ti_int_check=(tg3_Eint.get_data()[i]*nB+
			     tg3_Pint.get_data()[i]-
			     T_MeV*tg3_Sint.get_data()[i]*nB-
			     nn*tg3_mun.get_data()[i]-
			     np*tg3_mup.get_data()[i])/scale/nB;
	if (fabs(ti_int_check)>1.0e-9) {
	  cout << "Therm. ident. doens't hold (i,nB,Ye,T,ti_int_check): "
	       << i << " " << nB << " " << Ye << " " << T_MeV
	       << " " << ti_int_check << endl;
	  ti_int_count++;
	  if (ti_int_count==1000) {
	    cout << "Further Fint warnings suppressed." << endl;
	  }
	}
      }

      // Check that E*nB=-P+T*S*nB+nn*mun+np*mup+ne*mue
      if (derivs_computed && ti_count<1000) {
	double nn=nB*(1.0-Ye);
	double np=nB*Ye;
	double scale=tg3_E.get_data()[i];
	if (scale<10.0) scale=10.0;
	double ti_check=(tg3_E.get_data()[i]*nB+
			 tg3_P.get_data()[i]-
			 T_MeV*tg3_S.get_data()[i]*nB-
			 nn*tg3_mun.get_data()[i]-
			 np*tg3_mup.get_data()[i]-
			 np*tg3_mue.get_data()[i])/scale/nB;
	if (fabs(ti_check)>1.0e-9) {
	  cout << "Therm. ident. doens't hold (i,nB,Ye,T,ti_check): "
	       << i << " " << nB << " " << Ye << " " << T_MeV
	       << " " << ti_check << endl;
	  ti_count++;
	  if (ti_count==1000) {
	    cout << "Further Fint warnings suppressed." << endl;
	  }
	}
      }

      // Check that F=E-T*S
      if (with_leptons_loaded && F_count<1000) {

	double scale=tg3_F.get_data()[i];
	if (scale<10.0) scale=10.0;
	double F_check=(tg3_F.get_data()[i]-
			tg3_E.get_data()[i]+
			T_MeV*tg3_S.get_data()[i])/scale;
	if (fabs(F_check)>1.0e-9) {
	  cout << "F doesn't match E-T*S (i,nB,Ye,T,F_check): "
	       << i << " " << nB << " " << Ye << " " << T_MeV
	       << " " << F_check << endl;
	  F_count++;
	  if (F_count==1000) {
	    cout << "Further Fint warnings suppressed." << endl;
	  }
	}
      
      }

      // Check that Sint>0
      if (false && S_neg_count<1000 && tg3_Sint.get_data()[i]<0.0) {
	cout << "Entropy per baryon negative (i,nB,Ye,T,Sint): "
	     << i << " " << nB << " " << Ye << " "
	     << T_MeV << " " << tg3_Sint.get_data()[i] << endl;
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
  if (tg3_Eint.total_size()>0) {
    cout << "Fint_count: " << Fint_count << endl;
  }
  if (derivs_computed) {
    cout << "ti_int_count: " << ti_int_count << endl;
    if (with_leptons_loaded) {
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
	tg3_Xn.unpack_index(i,ix);
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

  const vector<double> &d=tg3_flag.get_data();
  for(size_t i=0;i<tg3_flag.total_size();i++) {
    size_t szt_tmp=((size_t)((d[i]+10.0)*(1.0+1.0e-12)));
    if (szt_tmp>20) szt_tmp=21;
    counts[szt_tmp]++;
  }
  cout << "Counts for file 1: ";
  vector_out(cout,counts,true);

  // Read the second file

  en2.read_results(in2);

  for(size_t i=0;i<22;i++) counts[i]=0;
  const vector<double> &d2=en2.tg3_flag.get_data();
  for(size_t i=0;i<en2.tg3_flag.total_size();i++) {
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
  if (en2.with_leptons_loaded!=with_leptons_loaded) {
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
  
  size_t nx=tg3_flag.get_size(0);
  size_t ny=tg3_flag.get_size(1);
  size_t nz=tg3_flag.get_size(2);

  size_t c1=0, c2=0;
  
  for(size_t i=0;i<nx;i++) {
    double nB=nB_grid2[i];
    if (i%10==9) {
      cout << (i+1) << "/" << nx << endl;
    }
    for(size_t j=0;j<ny;j++) {
      double Ye=Ye_grid2[j];
      for(size_t k=0;k<nz;k++) {
	double T_MeV=T_grid2[k];
	double flag1=tg3_flag.get(i,j,k);
	double flag2=en2.tg3_flag.get(i,j,k);
	double Fint1=tg3_Fint.get(i,j,k);
	double Fint2=en2.tg3_Fint.get(i,j,k);

	if (fabs(flag1)<0.5 && (Fint1>1.0e90 ||
				!std::isfinite(Fint1) || Fint1<(-1.0e10))) {
	  flag1=0.0;
	  tg3_flag.set(i,j,k,0.0);
	  tg3_Fint.set(i,j,k,0.0);
	  tg3_Z.set(i,j,k,0.0);
	  tg3_A.set(i,j,k,0.0);
	  tg3_log_xn.set(i,j,k,0.0);
	  tg3_log_xp.set(i,j,k,0.0);
	  c1++;
	  cout << "Invalid free energy in table 1 (" << nB << ","
	       << Ye << "," << T_MeV << ")." << endl;
	  exit(-1);
	}
	if (fabs(flag2)<0.5 && (Fint2>1.0e90 ||
				!std::isfinite(Fint2)|| Fint2<(-1.0e10))) {
	  flag2=0.0;
	  en2.tg3_flag.set(i,j,k,0.0);
	  en2.tg3_Fint.set(i,j,k,0.0);
	  en2.tg3_Z.set(i,j,k,0.0);
	  en2.tg3_A.set(i,j,k,0.0);
	  en2.tg3_log_xn.set(i,j,k,0.0);
	  en2.tg3_log_xp.set(i,j,k,0.0);
	  cout << "Invalid free energy in table 2 (" << nB << ","
	       << Ye << "," << T_MeV << ")." << endl;
	  exit(-1);
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

	  tg3_flag.set(i,j,k,en2.tg3_flag.get(i,j,k));

	  tg3_log_xn.set(i,j,k,en2.tg3_log_xn.get(i,j,k));
	  tg3_log_xp.set(i,j,k,en2.tg3_log_xp.get(i,j,k));
	  tg3_Z.set(i,j,k,en2.tg3_Z.get(i,j,k));
	  tg3_A.set(i,j,k,en2.tg3_A.get(i,j,k));
	  
	  tg3_A_min.set(i,j,k,en2.tg3_A_min.get(i,j,k));
	  tg3_A_max.set(i,j,k,en2.tg3_A_max.get(i,j,k));
	  tg3_NmZ_min.set(i,j,k,en2.tg3_NmZ_min.get(i,j,k));
	  tg3_NmZ_max.set(i,j,k,en2.tg3_NmZ_max.get(i,j,k));
	  
	  tg3_Fint.set(i,j,k,en2.tg3_Fint.get(i,j,k));
	  tg3_Sint.set(i,j,k,en2.tg3_Sint.get(i,j,k));
	  tg3_Eint.set(i,j,k,en2.tg3_Eint.get(i,j,k));

	  tg3_Xn.set(i,j,k,en2.tg3_Xn.get(i,j,k));
	  tg3_Xp.set(i,j,k,en2.tg3_Xp.get(i,j,k));
	  tg3_Xalpha.set(i,j,k,en2.tg3_Xalpha.get(i,j,k));
	  tg3_Xnuclei.set(i,j,k,en2.tg3_Xnuclei.get(i,j,k));
	  tg3_Xd.set(i,j,k,en2.tg3_Xd.get(i,j,k));
	  tg3_Xt.set(i,j,k,en2.tg3_Xt.get(i,j,k));
	  tg3_XHe3.set(i,j,k,en2.tg3_XHe3.get(i,j,k));
	  tg3_XLi4.set(i,j,k,en2.tg3_XLi4.get(i,j,k));

	  if (derivs_computed) {
	    tg3_Pint.set(i,j,k,en2.tg3_Pint.get(i,j,k));
	    tg3_mun.set(i,j,k,en2.tg3_mun.get(i,j,k));
	    tg3_mup.set(i,j,k,en2.tg3_mup.get(i,j,k));
	    if (with_leptons_loaded) {
	      tg3_F.set(i,j,k,en2.tg3_F.get(i,j,k));
	      tg3_E.set(i,j,k,en2.tg3_E.get(i,j,k));
	      tg3_P.set(i,j,k,en2.tg3_P.get(i,j,k));
	      tg3_S.set(i,j,k,en2.tg3_S.get(i,j,k));
	      tg3_mue.set(i,j,k,en2.tg3_mue.get(i,j,k));
	    }
	  }

	  c2++;
	}
	/*
	  flag1=tg3_flag1.get(i,j,k);
	  flag2=en2.tg3_flag.get(i,j,k);
	  Fint1=tg3_Fint1.get(i,j,k);
	  Fint2=en2.tg3_Fint.get(i,j,k);
	  if (i==99 && j==39) {
	  cout << i << " " << j << " " << k << " "
	  << flag1 << " " << flag2 << " " << Fint1 << " "
	  << Fint2 << endl;
	  exit(-1);
	  }
	*/
      }
    }
  }

  cout << "Points removed from table 1: " << c1 << endl;
  cout << "Points copied from table 2 to table 1: " << c2 << endl;
  
  for(size_t i=0;i<22;i++) counts[i]=0;
  const vector<double> &d3=tg3_flag.get_data();
  for(size_t i=0;i<tg3_flag.total_size();i++) {
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
    if (with_leptons_loaded) {
      names.push_back("F");
      names.push_back("E");
      names.push_back("P");
      names.push_back("S");
      names.push_back("mue");
    }
  }
  
  vector<tensor_grid3<> *> ptrs;
  ptrs.push_back(&tg3_log_xn);
  ptrs.push_back(&tg3_log_xp);
  ptrs.push_back(&tg3_Z);
  ptrs.push_back(&tg3_A);
  ptrs.push_back(&tg3_Fint);
  ptrs.push_back(&tg3_Xn);
  ptrs.push_back(&tg3_Xp);
  ptrs.push_back(&tg3_Xalpha);
  ptrs.push_back(&tg3_Xnuclei);
  ptrs.push_back(&tg3_Xd);
  ptrs.push_back(&tg3_Xt);
  ptrs.push_back(&tg3_XHe3);
  ptrs.push_back(&tg3_XLi4);
  ptrs.push_back(&tg3_Sint);
  ptrs.push_back(&tg3_Eint);
  if (derivs_computed) {
    ptrs.push_back(&tg3_Pint);
    ptrs.push_back(&tg3_mun);
    ptrs.push_back(&tg3_mup);
    if (with_leptons_loaded) {
      ptrs.push_back(&tg3_F);
      ptrs.push_back(&tg3_E);
      ptrs.push_back(&tg3_P);
      ptrs.push_back(&tg3_S);
      ptrs.push_back(&tg3_mue);
    }
  }
  
  vector<tensor_grid3<> *> ptrs2;
  ptrs2.push_back(&en2.tg3_log_xn);
  ptrs2.push_back(&en2.tg3_log_xp);
  ptrs2.push_back(&en2.tg3_Z);
  ptrs2.push_back(&en2.tg3_A);
  ptrs2.push_back(&en2.tg3_Fint);
  ptrs2.push_back(&en2.tg3_Xn);
  ptrs2.push_back(&en2.tg3_Xp);
  ptrs2.push_back(&en2.tg3_Xalpha);
  ptrs2.push_back(&en2.tg3_Xnuclei);
  ptrs2.push_back(&en2.tg3_Xd);
  ptrs2.push_back(&en2.tg3_Xt);
  ptrs2.push_back(&en2.tg3_XHe3);
  ptrs2.push_back(&en2.tg3_XLi4);
  ptrs2.push_back(&en2.tg3_Sint);
  ptrs2.push_back(&en2.tg3_Eint);
  if (derivs_computed) {
    ptrs2.push_back(&en2.tg3_Pint);
    ptrs2.push_back(&en2.tg3_mun);
    ptrs2.push_back(&en2.tg3_mup);
    if (with_leptons_loaded) {
      ptrs2.push_back(&en2.tg3_F);
      ptrs2.push_back(&en2.tg3_E);
      ptrs2.push_back(&en2.tg3_P);
      ptrs2.push_back(&en2.tg3_S);
      ptrs2.push_back(&en2.tg3_mue);
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
	      if (tg3_flag.get(i,j,k)>9.9 && en2.tg3_flag.get(i,j,k)>9.9) {
		found=true;
		double v1=ptrs[ell]->get(i,j,k);
		double v2=ptrs2[ell]->get(i,j,k);
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
	       << ptrs[ell]->get(imax,jmax,kmax) << endl;
	  cout << "Value in file " << in2 << " is "
	       << ptrs2[ell]->get(imax,jmax,kmax) << endl;
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
  cout << "Beginning new table." << endl;

  calculator calc;
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

  size_t st[3]={n_nB2,n_Ye2,n_T2};
  tg3_log_xn.resize(3,st);
  tg3_log_xp.resize(3,st);
  tg3_Z.resize(3,st);
  tg3_A.resize(3,st);
  tg3_flag.resize(3,st);
  tg3_Fint.resize(3,st);
  tg3_Sint.resize(3,st);
  tg3_Eint.resize(3,st);
  
  tg3_Xn.resize(3,st);
  tg3_Xp.resize(3,st);
  tg3_Xalpha.resize(3,st);
  tg3_Xnuclei.resize(3,st);
  tg3_Xd.resize(3,st);
  tg3_Xt.resize(3,st);
  tg3_XHe3.resize(3,st);
  tg3_XLi4.resize(3,st);

  tg3_log_xn.set_grid_packed(packed);
  tg3_log_xp.set_grid_packed(packed);
  tg3_Z.set_grid_packed(packed);
  tg3_A.set_grid_packed(packed);
  tg3_flag.set_grid_packed(packed);
  tg3_Fint.set_grid_packed(packed);
  tg3_Sint.set_grid_packed(packed);
  tg3_Eint.set_grid_packed(packed);
  
  tg3_Xn.set_grid_packed(packed);
  tg3_Xp.set_grid_packed(packed);
  tg3_Xalpha.set_grid_packed(packed);
  tg3_Xnuclei.set_grid_packed(packed);
  tg3_Xd.set_grid_packed(packed);
  tg3_Xt.set_grid_packed(packed);
  tg3_XHe3.set_grid_packed(packed);
  tg3_XLi4.set_grid_packed(packed);

  tg3_log_xn.set_all(0.0);
  tg3_log_xp.set_all(0.0);
  tg3_Z.set_all(0.0);
  tg3_A.set_all(0.0);
  tg3_flag.set_all(0.0);
  tg3_Fint.set_all(0.0);
  tg3_Sint.set_all(0.0);
  tg3_Eint.set_all(0.0);

  tg3_Xn.set_all(0.0);
  tg3_Xp.set_all(0.0);
  tg3_Xalpha.set_all(0.0);
  tg3_Xnuclei.set_all(0.0);
  tg3_Xd.set_all(0.0);
  tg3_Xt.set_all(0.0);
  tg3_XHe3.set_all(0.0);
  tg3_XLi4.set_all(0.0);

  if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
    
    tg3_A_min.resize(3,st);
    tg3_A_max.resize(3,st);
    tg3_NmZ_min.resize(3,st);
    tg3_NmZ_max.resize(3,st);
    
    tg3_A_min.set_grid_packed(packed);
    tg3_A_max.set_grid_packed(packed);
    tg3_NmZ_min.set_grid_packed(packed);
    tg3_NmZ_max.set_grid_packed(packed);
    
    tg3_A_min.set_all(0.0);
    tg3_A_max.set_all(0.0);
    tg3_NmZ_min.set_all(0.0);
    tg3_NmZ_max.set_all(0.0);
  
  }
  
  if (derivs_computed) {
    tg3_Pint.resize(3,st);
    tg3_Pint.set_grid_packed(packed);
    tg3_Pint.set_all(0.0);
    tg3_mun.resize(3,st);
    tg3_mun.set_grid_packed(packed);
    tg3_mun.set_all(0.0);
    tg3_mup.resize(3,st);
    tg3_mup.set_grid_packed(packed);
    tg3_mup.set_all(0.0);
    if (with_leptons_loaded) {
      tg3_E.resize(3,st);
      tg3_E.set_grid_packed(packed);
      tg3_E.set_all(0.0);
      tg3_P.resize(3,st);
      tg3_P.set_grid_packed(packed);
      tg3_P.set_all(0.0);
      tg3_S.resize(3,st);
      tg3_S.set_grid_packed(packed);
      tg3_S.set_all(0.0);
      tg3_F.resize(3,st);
      tg3_F.set_grid_packed(packed);
      tg3_F.set_all(0.0);
      tg3_mue.resize(3,st);
      tg3_mue.set_grid_packed(packed);
      tg3_mue.set_all(0.0);
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
  
  calculator calc, calc2;
  std::map<std::string,double> vars;
  calc.compile(select_func.c_str());
  if (sv.size()>2) {
    calc2.compile(value_func.c_str());
  }

  for(int inB=0;inB<((int)n_nB2);inB++) {
    for(int iYe=0;iYe<((int)n_Ye2);iYe++) {
      for(int iT=0;iT<((int)n_T2);iT++) {
	
	vars["inB"]=inB;
	vars["iYe"]=iYe;
	vars["iT"]=iT;
	vars["nnB"]=nB_grid2.size();
	vars["nYe"]=Ye_grid2.size();
	vars["nT"]=T_grid2.size();
	vars["nB"]=nB_grid2[inB];
	vars["Ye"]=Ye_grid2[iYe];
	vars["T"]=T_grid2[iT];
	
	vars["flag"]=tg3_flag.get(inB,iYe,iT);
	vars["Fint"]=tg3_Fint.get(inB,iYe,iT);
	vars["Sint"]=tg3_Sint.get(inB,iYe,iT);
	vars["Eint"]=tg3_Eint.get(inB,iYe,iT);
	vars["Z"]=tg3_Z.get(inB,iYe,iT);
	vars["A"]=tg3_A.get(inB,iYe,iT);
	vars["log_xn"]=tg3_log_xn.get(inB,iYe,iT);
	vars["log_xp"]=tg3_log_xp.get(inB,iYe,iT);
	
	vars["Xn"]=tg3_Xn.get(inB,iYe,iT);
	vars["Xp"]=tg3_Xp.get(inB,iYe,iT);
	vars["Xalpha"]=tg3_Xalpha.get(inB,iYe,iT);
	vars["Xnuclei"]=tg3_Xnuclei.get(inB,iYe,iT);
	vars["Xd"]=tg3_Xd.get(inB,iYe,iT);
	vars["Xt"]=tg3_Xt.get(inB,iYe,iT);
	vars["XHe3"]=tg3_XHe3.get(inB,iYe,iT);
	vars["XLi4"]=tg3_XLi4.get(inB,iYe,iT);
	vars["dAdnB"]=0.0;
	if (inB>0 && inB<((int)n_nB2)-1) {
	  double a_low=tg3_A.get(inB-1,iYe,iT);
	  double a_mid=tg3_A.get(inB,iYe,iT);
	  double a_high=tg3_A.get(inB+1,iYe,iT);
	  if (fabs(a_mid-a_low)>fabs(a_mid-a_high)) {
	    vars["dAdnB"]=fabs(a_mid-a_low);
	  } else {
	    vars["dAdnB"]=fabs(a_mid-a_high);
	  }
	}

	if (alg_mode==2 || alg_mode==3 || alg_mode==4) {
	  vars["A_min"]=tg3_A_min.get(inB,iYe,iT);
	  vars["A_max"]=tg3_A_max.get(inB,iYe,iT);
	  vars["NmZ_min"]=tg3_NmZ_min.get(inB,iYe,iT);
	  vars["NmZ_max"]=tg3_NmZ_max.get(inB,iYe,iT);
	}

	if (derivs_computed) {
	  vars["Pint"]=tg3_Pint.get(inB,iYe,iT);
	  vars["mun"]=tg3_mun.get(inB,iYe,iT);
	  vars["mup"]=tg3_mup.get(inB,iYe,iT);
	}
	
	double val=calc.eval(&vars);
	if (val>0.5) {
	  count++;
	  if (sv.size()>3) {
	    double val2=calc2.eval(&vars);
	    if (tensor_to_change=="flag") {
	      tg3_flag.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="Fint") {
	      tg3_Fint.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="Z") {
	      tg3_Z.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="A") {
	      tg3_A.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="log_xn") {
	      tg3_log_xn.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="log_xp") {
	      tg3_log_xp.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="Xn") {
	      tg3_Xn.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="Xp") {
	      tg3_Xp.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="Xalpha") {
	      tg3_Xalpha.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="Xnuclei") {
	      tg3_Xnuclei.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="Xd") {
	      tg3_Xd.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="Xt") {
	      tg3_Xt.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="XHe3") {
	      tg3_XHe3.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="XLi4") {
	      tg3_XLi4.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="A_min") {
	      tg3_A_min.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="A_max") {
	      tg3_A_max.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="NmZ_min") {
	      tg3_NmZ_min.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="NmZ_max") {
	      tg3_NmZ_max.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="Eint") {
	      tg3_Eint.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="Pint") {
	      tg3_Pint.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="Sint") {
	      tg3_Sint.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="F") {
	      tg3_F.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="E") {
	      tg3_E.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="P") {
	      tg3_P.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="S") {
	      tg3_S.get(inB,iYe,iT)=val2;
	    } else if (tensor_to_change=="mue") {
	      tg3_mue.get(inB,iYe,iT)=val2;
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
  if (with_leptons_loaded) {
    cout << "Function generate-table setting with_leptons to false."
	 << endl;
    with_leptons_loaded=false;
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
	  int iflag=((int)(tg3_flag.get(inB,iYe,iT)*(1.0+1.0e-12)));
	  // If recompute is true, set all finished points to "guess"
	  if (recompute==true) {
	    if (iflag==iflag_done) {
	      tg3_flag.get(inB,iYe,iT)=iflag_guess;
	    }
	  }
	  // Get rid of previous in progress markings
	  if (iflag==iflag_in_progress_empty) {
	    tg3_flag.get(inB,iYe,iT)=iflag_empty;
	  }
	  if (iflag==iflag_in_progress_with_guess) {
	    tg3_flag.get(inB,iYe,iT)=iflag_guess;
	  }
	}
      }
    }

    if (edge_list.length()>0) {
      vector<string> edge_vec;
      split_string(edge_list,edge_vec);
      for(size_t k=0;k<edge_vec.size();k++) {
	tensor_grid3<> *ptr;
	if (edge_vec[k]=="A") {
	  ptr=&tg3_A;
	} else if (edge_vec[k]=="Xnuclei") {
	  ptr=&tg3_Xnuclei;
	} else if (edge_vec[k]=="Xalpha") {
	  ptr=&tg3_Xalpha;
	} else {
	  cerr << "Invalid name in edge_list." << endl;
	  return 5;
	}
	for(int inB=1;inB<((int)n_nB2-1);inB++) {
	  for(int iYe=1;iYe<((int)n_Ye2-1);iYe++) {
	    for(int iT=1;iT<((int)n_T2-1);iT++) {
	      int iflag=((int)(tg3_flag.get(inB,iYe,iT)*(1.0+1.0e-12)));
	      if (iflag==iflag_done) {
		double X0=ptr->get(inB,iYe,iT);

		if (X0>1.0e-4) {

		  if (tg3_flag.get(inB-1,iYe,iT)>9.9 &&
		      tg3_flag.get(inB+1,iYe,iT)>9.9) {
		    double X1=ptr->get(inB-1,iYe,iT);
		    double X2=ptr->get(inB+1,iYe,iT);
		
		    if (X0<std::min(X1,X2) || X0>std::max(X1,X2)) {
		      tg3_flag.get(inB-1,iYe,iT)=5.0;
		      tg3_flag.get(inB+1,iYe,iT)=5.0;
		      tg3_flag.get(inB,iYe,iT)=5.0;
		    }
		  }
		  if (false) {
		    if (tg3_flag.get(inB,iYe-1,iT)>9.9 &&
			tg3_flag.get(inB,iYe+1,iT)>9.9) {
		      double X1=ptr->get(inB,iYe-1,iT);
		      double X2=ptr->get(inB,iYe+1,iT);
		    
		      if (X0<std::min(X1,X2) || X0>std::max(X1,X2)) {
			tg3_flag.get(inB,iYe,iT)=5.0;
			tg3_flag.get(inB,iYe-1,iT)=5.0;
			tg3_flag.get(inB,iYe+1,iT)=5.0;
		      }
		    }
		  }
		  if (tg3_flag.get(inB,iYe,iT-1)>9.9 &&
		      tg3_flag.get(inB,iYe,iT+1)>9.9) {
		    double X1=ptr->get(inB,iYe,iT-1);
		    double X2=ptr->get(inB,iYe,iT+1);
		  
		    if (X0<std::min(X1,X2) || X0>std::max(X1,X2)) {
		      tg3_flag.get(inB,iYe,iT)=5.0;
		      tg3_flag.get(inB,iYe,iT-1)=5.0;
		      tg3_flag.get(inB,iYe,iT+1)=5.0;
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
      
      // -----------------------------------------------------
      // Compute tasks (either w/o MPI or w/MPI on rank 0)
    
      // Setup task list as a two sets of triplets, source first,
      // destination second
      vector<size_t> tasks;

      // Store initial guesses in a table, one row for each task
      table<> gtab;
      gtab.line_of_names("log_xn log_xp Z A A_min A_max NmZ_min NmZ_max");

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
	    
	    int iflag=((int)(tg3_flag.get(inB,iYe,iT)*(1.0+1.0e-12)));
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
		double line[8]={tg3_log_xn.get(inB,iYe,iT),
				tg3_log_xp.get(inB,iYe,iT),
				0.0,0.0,
				tg3_A_min.get(inB,iYe,iT),
				tg3_A_max.get(inB,iYe,iT),
				tg3_NmZ_min.get(inB,iYe,iT),
				tg3_NmZ_max.get(inB,iYe,iT)};
		gtab.line_of_data(8,line);
	      } else {
		double line[8]={tg3_log_xn.get(inB,iYe,iT),
				tg3_log_xp.get(inB,iYe,iT),
				tg3_Z.get(inB,iYe,iT),
				tg3_A.get(inB,iYe,iT),
				0.0,0.0,0.0,0.0};
		gtab.line_of_data(8,line);
	      }	
	      
	      tg3_flag.get(inB,iYe,iT)=(double)iflag_in_progress_with_guess;
	      guess_found=true;
	    }

	    if (ext_guess.length()>0 &&
		external.tg3_flag.get(inB,iYe,iT)>9.9 &&
		(iflag==iflag_guess || iflag==iflag_empty)) {
	      tasks.push_back(inB);
	      tasks.push_back(iYe);
	      tasks.push_back(iT);
	      tasks.push_back(inB);
	      tasks.push_back(iYe);
	      tasks.push_back(iT);
	      if (alg_mode>=2) {
		double line[8]={external.tg3_log_xn.get(inB,iYe,iT),
				external.tg3_log_xp.get(inB,iYe,iT),
				0.0,0.0,
				external.tg3_A_min.get(inB,iYe,iT),
				external.tg3_A_max.get(inB,iYe,iT),
				external.tg3_NmZ_min.get(inB,iYe,iT),
				external.tg3_NmZ_max.get(inB,iYe,iT)};
		gtab.line_of_data(8,line);
	      } else {
		double line[8]={external.tg3_log_xn.get(inB,iYe,iT),
				external.tg3_log_xp.get(inB,iYe,iT),
				external.tg3_Z.get(inB,iYe,iT),
				external.tg3_A.get(inB,iYe,iT),
				0.0,0.0,0.0,0.0};
		gtab.line_of_data(8,line);
	      }
	      guess_found=true;
	    }
	    
	    // If six_neighbors is true, set up six additional tasks
	    // for neighbors with useful initial guesses
	    if (six_neighbors>0 &&
		(iflag==iflag_guess ||
		 (propagate_points && iflag!=iflag_done))) {

	      if (inB>0) {
		int iflag2=((int)(tg3_flag.get(inB-1,iYe,iT)*
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
		    // If we're at small densities and we've converged
		    // at both a higher density and a lower density,
		    // then we can linearly interpolate to get a good
		    // guess. At higher densities this is not
		    // necessarily a good guess because of the
		    // liquid-gas phase transition
		    if (nB_grid2[inB]<1.0e-3 && inB<((int)n_nB2)-1 &&
			tg3_flag.get(inB+1,iYe,iT)>9.9) {
		      double line[8]={(tg3_log_xn.get(inB-1,iYe,iT)+
				       tg3_log_xn.get(inB+1,iYe,iT))/2.0,
				      (tg3_log_xp.get(inB-1,iYe,iT)+
				       tg3_log_xp.get(inB-1,iYe,iT))/2.0,
				      0.0,0.0,
				      (tg3_A_min.get(inB-1,iYe,iT)+
				       tg3_A_min.get(inB-1,iYe,iT))/2.0,
				      (tg3_A_max.get(inB-1,iYe,iT)+
				       tg3_A_max.get(inB-1,iYe,iT))/2.0,
				      (tg3_NmZ_min.get(inB-1,iYe,iT)+
				       tg3_NmZ_min.get(inB-1,iYe,iT))/2.0,
				      (tg3_NmZ_max.get(inB-1,iYe,iT)+
				       tg3_NmZ_max.get(inB-1,iYe,iT))/2.0};
		      gtab.line_of_data(8,line);
		    } else {
		      double line[8]={tg3_log_xn.get(inB-1,iYe,iT),
				      tg3_log_xp.get(inB-1,iYe,iT),
				      0.0,0.0,
				      tg3_A_min.get(inB-1,iYe,iT),
				      tg3_A_max.get(inB-1,iYe,iT),
				      tg3_NmZ_min.get(inB-1,iYe,iT),
				      tg3_NmZ_max.get(inB-1,iYe,iT)};
		      gtab.line_of_data(8,line);
		    }
		  } else {
		    double line[8]={tg3_log_xn.get(inB-1,iYe,iT),
				    tg3_log_xp.get(inB-1,iYe,iT),
				    tg3_Z.get(inB-1,iYe,iT),
				    tg3_A.get(inB-1,iYe,iT),
				    0.0,0.0,0.0,0.0};
		    gtab.line_of_data(8,line);
		  }	
		  
		  guess_found=true;
		}
	      }
	      
	      if (six_neighbors>4 && iYe>0) {
		int iflag2=((int)(tg3_flag.get(inB,iYe-1,iT)*
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
		    double line[8]={tg3_log_xn.get(inB,iYe-1,iT),
				    tg3_log_xp.get(inB,iYe-1,iT),
				    0.0,0.0,
				    tg3_A_min.get(inB,iYe-1,iT),
				    tg3_A_max.get(inB,iYe-1,iT),
				    tg3_NmZ_min.get(inB,iYe-1,iT),
				    tg3_NmZ_max.get(inB,iYe-1,iT)};
		    gtab.line_of_data(8,line);
		  } else {
		    double line[8]={tg3_log_xn.get(inB,iYe-1,iT),
				    tg3_log_xp.get(inB,iYe-1,iT),
				    tg3_Z.get(inB,iYe-1,iT),
				    tg3_A.get(inB,iYe-1,iT),
				    0.0,0.0,0.0,0.0};
		    gtab.line_of_data(8,line);
		  }	
		  
		  guess_found=true;
		}
	      }
	      if (six_neighbors>2 && iT>0) {
		int iflag2=((int)(tg3_flag.get(inB,iYe,iT-1)*
				  (1.0+1.0e-12)));
		if (iflag2==iflag_in_progress_with_guess ||
		    iflag2==iflag_guess || iflag2==iflag_done) {
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT-1);
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT);

		  if (alg_mode>=2) {
		    double line[8]={tg3_log_xn.get(inB,iYe,iT-1),
				    tg3_log_xp.get(inB,iYe,iT-1),
				    0.0,0.0,
				    tg3_A_min.get(inB,iYe,iT-1),
				    tg3_A_max.get(inB,iYe,iT-1),
				    tg3_NmZ_min.get(inB,iYe,iT-1),
				    tg3_NmZ_max.get(inB,iYe,iT-1)};
		    gtab.line_of_data(8,line);
		  } else {
		    double line[8]={tg3_log_xn.get(inB,iYe,iT-1),
				    tg3_log_xp.get(inB,iYe,iT-1),
				    tg3_Z.get(inB,iYe,iT-1),
				    tg3_A.get(inB,iYe,iT-1),
				    0.0,0.0,0.0,0.0};
		    gtab.line_of_data(8,line);
		  }	
		  
		  guess_found=true;
		}
	      }
	      if (six_neighbors>1 && inB<((int)n_nB2)-1) {
		int iflag2=((int)(tg3_flag.get(inB+1,iYe,iT)*
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
			tg3_flag.get(inB-1,iYe,iT)>9.9) {
		      double line[8]={(tg3_log_xn.get(inB-1,iYe,iT)+
				       tg3_log_xn.get(inB+1,iYe,iT))/2.0,
				      (tg3_log_xp.get(inB-1,iYe,iT)+
				       tg3_log_xp.get(inB-1,iYe,iT))/2.0,
				      0.0,0.0,
				      (tg3_A_min.get(inB-1,iYe,iT)+
				       tg3_A_min.get(inB-1,iYe,iT))/2.0,
				      (tg3_A_max.get(inB-1,iYe,iT)+
				       tg3_A_max.get(inB-1,iYe,iT))/2.0,
				      (tg3_NmZ_min.get(inB-1,iYe,iT)+
				       tg3_NmZ_min.get(inB-1,iYe,iT))/2.0,
				      (tg3_NmZ_max.get(inB-1,iYe,iT)+
				       tg3_NmZ_max.get(inB-1,iYe,iT))/2.0};
		      gtab.line_of_data(8,line);
		    } else {
		      double line[8]={tg3_log_xn.get(inB+1,iYe,iT),
				      tg3_log_xp.get(inB+1,iYe,iT),
				      0.0,0.0,
				      tg3_A_min.get(inB+1,iYe,iT),
				      tg3_A_max.get(inB+1,iYe,iT),
				      tg3_NmZ_min.get(inB+1,iYe,iT),
				      tg3_NmZ_max.get(inB+1,iYe,iT)};
		      gtab.line_of_data(8,line);
		    }
		  } else {
		    double line[8]={tg3_log_xn.get(inB+1,iYe,iT),
				    tg3_log_xp.get(inB+1,iYe,iT),
				    tg3_Z.get(inB+1,iYe,iT),
				    tg3_A.get(inB+1,iYe,iT),
				    0.0,0.0,0.0,0.0};
		    gtab.line_of_data(8,line);
		  }	
		  
		  guess_found=true;
		}
	      }
	      if (six_neighbors>4 && iYe<((int)n_Ye2)-1) {
		int iflag2=((int)(tg3_flag.get(inB,iYe+1,iT)*
				  (1.0+1.0e-12)));
		if (iflag2==iflag_in_progress_with_guess ||
		    iflag2==iflag_guess || iflag2==iflag_done) {
		  tasks.push_back(inB);
		  tasks.push_back(iYe+1);
		  tasks.push_back(iT);
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT);

		  if (alg_mode>=2) {
		    double line[8]={tg3_log_xn.get(inB,iYe+1,iT),
				    tg3_log_xp.get(inB,iYe+1,iT),
				    0.0,0.0,
				    tg3_A_min.get(inB,iYe+1,iT),
				    tg3_A_max.get(inB,iYe+1,iT),
				    tg3_NmZ_min.get(inB,iYe+1,iT),
				    tg3_NmZ_max.get(inB,iYe+1,iT)};
		    gtab.line_of_data(8,line);
		  } else {
		    double line[8]={tg3_log_xn.get(inB,iYe+1,iT),
				    tg3_log_xp.get(inB,iYe+1,iT),
				    tg3_Z.get(inB,iYe+1,iT),
				    tg3_A.get(inB,iYe+1,iT),
				    0.0,0.0,0.0,0.0};
		    gtab.line_of_data(8,line);
		  }	
		  
		  guess_found=true;
		}
	      }
	      if (six_neighbors>2 && iT<((int)n_T2)-1) {
		int iflag2=((int)(tg3_flag.get(inB,iYe,iT+1)*
				  (1.0+1.0e-12)));
		if (iflag2==iflag_in_progress_with_guess ||
		    iflag2==iflag_guess || iflag2==iflag_done) {
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT+1);
		  tasks.push_back(inB);
		  tasks.push_back(iYe);
		  tasks.push_back(iT);

		  if (alg_mode>=2) {
		    double line[8]={tg3_log_xn.get(inB,iYe,iT+1),
				    tg3_log_xp.get(inB,iYe,iT+1),
				    0.0,0.0,
				    tg3_A_min.get(inB,iYe,iT+1),
				    tg3_A_max.get(inB,iYe,iT+1),
				    tg3_NmZ_min.get(inB,iYe,iT+1),
				    tg3_NmZ_max.get(inB,iYe,iT+1)};
		    gtab.line_of_data(8,line);
		  } else {
		    double line[8]={tg3_log_xn.get(inB,iYe,iT+1),
				    tg3_log_xp.get(inB,iYe,iT+1),
				    tg3_Z.get(inB,iYe,iT+1),
				    tg3_A.get(inB,iYe,iT+1),
				    0.0,0.0,0.0,0.0};
		    gtab.line_of_data(8,line);
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

		    // Ensure that the j indices are in bounds
		    if ((jnB!=inB || jYe!=iYe || jT!=iT) && jnB>=0 &&
			jYe>=0 && jT>=0 && jnB<((int)n_nB2) &&
			jYe<((int)n_Ye2) && jT<((int)n_T2)) {

		      int jflag=((int)(tg3_flag.get(jnB,jYe,jT)*
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
			
			if (alg_mode>=2) {
			  double line[8]={tg3_log_xn.get(jnB,jYe,jT),
					  tg3_log_xp.get(jnB,jYe,jT),
					  0.0,0.0,
					  tg3_A_min.get(jnB,jYe,jT),
					  tg3_A_max.get(jnB,jYe,jT),
					  tg3_NmZ_min.get(jnB,jYe,jT),
					  tg3_NmZ_max.get(jnB,jYe,jT)};
			  gtab.line_of_data(8,line);
			} else {
			  double line[8]={tg3_log_xn.get(jnB,jYe,jT),
					  tg3_log_xp.get(jnB,jYe,jT),
					  tg3_Z.get(jnB,jYe,jT),
					  tg3_A.get(jnB,jYe,jT),
					  0.0,0.0,0.0,0.0};
			  gtab.line_of_data(8,line);
			}	
		  
			if (iflag==iflag_guess) {
			  tg3_flag.get(inB,iYe,iT)=
			    (double)iflag_in_progress_with_guess;
			} else {
			  tg3_flag.get(inB,iYe,iT)=
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

	    // Determine if the "no_nuclei" flag should be set
	    output_buffers[proc_index][vi["no_nuclei"]]=0.0;
	    {
	      size_t inB_dest=tasks[i*6+3];
	      size_t iYe_dest=tasks[i*6+4];
	      size_t iT_dest=tasks[i*6+5];
	      if (iT_dest>0) {
		if (tg3_flag.get(inB_dest,iYe_dest,iT_dest-1)>9.9) {
		  double X_all_nuclei=0.0;
		  X_all_nuclei+=tg3_Xalpha.get(inB_dest,iYe_dest,iT_dest-1);
		  X_all_nuclei+=tg3_Xnuclei.get(inB_dest,iYe_dest,iT_dest-1);
		  X_all_nuclei+=tg3_Xd.get(inB_dest,iYe_dest,iT_dest-1);
		  X_all_nuclei+=tg3_Xt.get(inB_dest,iYe_dest,iT_dest-1);
		  X_all_nuclei+=tg3_XHe3.get(inB_dest,iYe_dest,iT_dest-1);
		  X_all_nuclei+=tg3_XLi4.get(inB_dest,iYe_dest,iT_dest-1);
		  if (X_all_nuclei<1.0e-20 && nB_grid2[tasks[i*6+3]]<0.02) {
		    cout << "X_all_nuclei small: " << X_all_nuclei << " "
			 << T_grid2[iT_dest] << endl;
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
		  if (tg3_flag.get(ii,j,k)>=10.0) conv2_count++;
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
	  cout << "  nB,Ye,T: " << nB << " " << Ye << " " << T*hc_mev_fm
	       << endl;
	  double log_xn=gtab.get("log_xn",i);
	  double log_xp=gtab.get("log_xp",i);
	  size_t nuc_Z1=((size_t)(gtab.get("Z",i)+1.0e-12));
	  size_t nuc_N1=((size_t)(gtab.get("A",i)+1.0e-12))-nuc_Z1;
	  bool no_nuclei=false;
	  
	  {
	    if (iT>0) {
	      if (tg3_flag.get(inB,iYe,iT-1)>9.9) {
		double X_all_nuclei=0.0;
		X_all_nuclei+=tg3_Xalpha.get(inB,iYe,iT-1);
		X_all_nuclei+=tg3_Xnuclei.get(inB,iYe,iT-1);
		X_all_nuclei+=tg3_Xd.get(inB,iYe,iT-1);
		X_all_nuclei+=tg3_Xt.get(inB,iYe,iT-1);
		X_all_nuclei+=tg3_XHe3.get(inB,iYe,iT-1);
		X_all_nuclei+=tg3_XLi4.get(inB,iYe,iT-1);
		if (X_all_nuclei<1.0e-20) {
		  cout << "X_all_nuclei small: " << X_all_nuclei << " "
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
	    cout << "Point at (nB,Ye,T)=("
		 << nB << "," << Ye << "," << T*hc_mev_fm << "), ret="
		 << ret << ", " << i << "/" << ntasks << endl;
	  }
	  if (ret==0) {
	    
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
                  if (tg3_flag.get(ii,j,k)>=10.0) conv2_count++;
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
	
      } else if (true || write_elapsed>file_update_time) {
	
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
	      if (tg3_flag.get(i,j,k)>=10.0) conv2_count++;
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
	  if (tg3_flag.get(i,j,k)>=10.0) conv2_count++;
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
      if (max_time==0.0 || elapsed<max_time) {

	if (gt_verbose>1) {
	  cout << "Rank " << mpi_rank
	       << " computing point at nB,Ye,T(MeV): " << nB << " " 
	       << Ye << " " << T*hc_mev_fm << endl;
	}
	
	map<string,double> vdet;
	
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
		 << " done. nB,Ye,T(MeV),F(MeV): " << nB << " " 
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
		 << " failed. nB,Ye,T(MeV),ret: " << nB << " " 
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

  calculator calc;
  std::map<std::string,double> vars;

  vector<double> packed, packed2;
  vector<std::string> split_res;

  //o2scl::tensor_grid3<> tg3_zn;
  //o2scl::tensor_grid3<> tg3_zp;
  o2scl::tensor_grid3<> tg3_zn2;
  o2scl::tensor_grid3<> tg3_zp2;

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
  tg3_zn.resize(3,st);
  tg3_zp.resize(3,st);

  tg3_zn.set_grid_packed(packed);
  tg3_zp.set_grid_packed(packed);

  tg3_zn.set_all(0.0);
  tg3_zp.set_all(0.0);
  
  size_t st2[3]={21,n_Ye,n_T};
  tg3_zn2.resize(3,st2);
  tg3_zp2.resize(3,st2);

  tg3_zn2.set_grid_packed(packed2);
  tg3_zp2.set_grid_packed(packed2);

  tg3_zn2.set_all(0.0);
  tg3_zp2.set_all(0.0);
  
  for(int inB=0;inB<((int)n_nB);inB++) {
    double nB=nB_grid[inB];
    for(int iYe=0;iYe<((int)n_Ye);iYe++) {
      double Ye=Ye_grid[iYe];
      for(int iT=0;iT<((int)n_T);iT++) {
	double T=T_grid[iT];
	double lambda=sqrt(4.0*o2scl_const::pi/(neutron.m+proton.m)/T);
	
	double b_n=ecv.bn_f(T);
	double b_pn=ecv.bpn_f(T);
	
	double zn, zp;
	//cout << "A: " << nB << " " << Ye << " " << T << " " << endl;
	vsd.solve_fugacity(nB*(1.0-Ye),nB*Ye,lambda,
			   lambda,b_n,b_pn,zn,zp);
	
	tg3_zn.get(inB,iYe,iT)=zn;
	tg3_zp.get(inB,iYe,iT)=zp;
      }
    }
  }

  for(int inB=0;inB<21;inB++) {
    double nB=pow(10.0,-50.0+inB*2);
    for(int iYe=0;iYe<((int)n_Ye);iYe++) {
      double Ye=Ye_grid[iYe];
      for(int iT=0;iT<((int)n_T);iT++) {
	double T=T_grid[iT];
	double lambda=sqrt(4.0*o2scl_const::pi/(neutron.m+proton.m)/T);
	
	double b_n=ecv.bn_f(T);
	double b_pn=ecv.bpn_f(T);
	
	double zn, zp;
	//cout << "B: " << nB << " " << Ye << " " << T << " " << endl;
	vsd.solve_fugacity(nB*(1.0-Ye),nB*Ye,lambda,
			   lambda,b_n,b_pn,zn,zp);

	tg3_zn2.get(inB,iYe,iT)=log10(zn);
	tg3_zp2.get(inB,iYe,iT)=log10(zp);
      }
    }
  }

  hdf_file hf;
  hf.open_or_create("check_virial.o2");
  hdf_output(hf,tg3_zn,"zn");
  hdf_output(hf,tg3_zp,"zp");
  hdf_output(hf,tg3_zn2,"zn2");
  hdf_output(hf,tg3_zp2,"zp2");
  hf.close();
  cout << "Created check_virial.o2" << endl;
  
  return 0;
}

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
      if (use_skalt) {
	line.push_back(eosp_alt->calc_temp_e(neutron,proton,T,th2)/
		       nB_arr[k]*hc_mev_fm);
      } else {
	line.push_back(free_energy_density(neutron,proton,T,th2)/
		       nB_arr[k]*hc_mev_fm);
      }
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

  if (tg3_flag.get(inB,iYe,iT)<9.9) {
    cerr << "Point not converged." << endl;
    return 3;
  }

  table<> tab;
  tab.line_of_names("log_xn log_xp Z N A ZoA Xnuclei");

  static const int N=1000;
  for(int j=0;j<N;j++) {
    std::cout << "j: " << j << endl;
    
    std::vector<std::string> obj;
    random(obj,false);
      
    double log_xn, log_xp;
    int A_min=((int)(tg3_A_min.get(inB,iYe,iT)));
    int A_max=((int)(tg3_A_max.get(inB,iYe,iT)));
    int NmZ_min=((int)(tg3_NmZ_min.get(inB,iYe,iT)));
    int NmZ_max=((int)(tg3_NmZ_max.get(inB,iYe,iT)));
    double mun_full, mup_full;
    thermo thx;
    bool dist_changed=true;
    bool no_nuclei=false;

    double log_xn_best=0.0, log_xp_best=0.0;
    double fr_best=1.0e10;
    double Zbar, Nbar;
    
    for(size_t k=0;k<10;k++) {
      
      log_xn=tg3_log_xn.get(inB,iYe,iT)+(rng.random()*6.0-3.0);
      log_xp=tg3_log_xp.get(inB,iYe,iT)+(rng.random()*6.0-3.0);
      
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
      double line[7]={log_xn,log_xp,Zbar,Nbar,Zbar+Nbar,
		      Zbar/(Zbar+Nbar),X[5]};
      tab.line_of_data(7,line);
    }

    if (j%100==0) {
      hdf_file hf;
      hf.open_or_create(sv[4]);
      hdf_output(hf,tab,"mn2");
      hf.close();
    }
    
  }

  return 0;
}

void eos_nuclei::setup_cli(o2scl::cli &cl) {
  
  eos::setup_cli(cl);
  
  static const int nopt=21;
  
  o2scl::comm_option_s options[nopt]=
    {{0,"eos-deriv","compute derivatives",
      0,0,"<file in> <file out>","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::eos_deriv),
      o2scl::cli::comm_option_both},
     {0,"add-eg","Add electrons and photons.",
      0,0,"","",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::add_eg),
      o2scl::cli::comm_option_both},
     {0,"eg-table","Make electron and photon table.",
      1,1,"<output file>","",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::eg_table),
      o2scl::cli::comm_option_both},
     {0,"maxwell-test","Experimental test of Maxwell construction.",
      4,4,"<file in> <file out>","",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::maxwell_test),
      o2scl::cli::comm_option_both},
     {0,"fit-frdm","Fit the FRDM mass model.",
      0,0,"","",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::fit_frdm),
      o2scl::cli::comm_option_both},
     {0,"check-virial","Check the virial EOS.",
      0,0,"","",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::check_virial),
      o2scl::cli::comm_option_both},
     {0,"generate-table","Generate an EOS table.",
      0,1,"[out file]",
      "",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::generate_table),o2scl::cli::comm_option_both},
     {0,"load","Load an EOS table.",0,1,"<filename>",
      "",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::load),o2scl::cli::comm_option_both},
     {0,"output","Output an EOS table.",
      0,1,"<filename>",
      "",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::output),o2scl::cli::comm_option_both},
     {0,"edit-data","Edit data in the EOS tables.",1,4,
      "<select func.> [tensor to modify] [value func.]",
      ((string)"The \"edit-data\" command counts the number of ")+
      "(nB,Ye,T) points matching the "+
      "criteria specified in <select func.>. If the remaining two "+
      "arguments are given, then the values of [tensor to modify] "+
      "for the selected points "+
      "are changed to the result of the function [value func.].",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::edit_data),o2scl::cli::comm_option_both},
     {0,"merge-tables","Merge two output tables to create a third.",
      3,3,"<input file 1> <input file 2> <output file>",
      "",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::merge_tables),o2scl::cli::comm_option_both},
     {0,"write-nuclei","Output nuclear masses.",
      1,1,"<output file>",
      "",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::write_nuclei),o2scl::cli::comm_option_both},
     {0,"compare-tables","Compare two output tables.",
      2,3,"<input file 1> <input file 2> [quantity]",
      ((std::string)"Compare two EOS tables. If the optional argument ")+
      "is unspecified, then all quantities are compared. If [quantity] "+
      "is specified, then only that particular quantitiy is compared. "+
      "Only points for which flag=10 in both tables are compared. "+
      "If derivs_computed is true, then Pint, mun, and "+
      "mup are available for comparisons. If with_leptons_loaded is "+
      "true, then "+
      "F, E, P, and S, are also available for comparisons. Any current "+
      "EOS data stored is cleared before the comparison. If the "+
      "nB, Ye, or T grids do not match, then no comparison is performed.",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::compare_tables),o2scl::cli::comm_option_both},
     {0,"stats","Output convergence statistics and simple checks.",0,0,"",
      "",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::stats),o2scl::cli::comm_option_both},
     {0,"mcarlo-nuclei","",
      0,0,"<file>",
      "",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::mcarlo_nuclei),o2scl::cli::comm_option_both},
     {0,"mcarlo-nuclei2","",
      4,4,"<nB> <Ye> <T> <filename>",
      "",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::mcarlo_nuclei2),o2scl::cli::comm_option_both},
     {0,"point-nuclei",
      "Compute and/or show EOS results at one (n_B,Y_e,T) point.",
      -1,-1,((string)"<n_B> <Y_e> <T (MeV)> [log(xn) log(xp) Z N] ")+
      "[alg_mode 2-3: log(xn) log(xp) A_min A_max NmZ_min NmZ_max] [fname]",
      ((std::string)"If an EOS is loaded, then the n_B, Y_e, and T ")+
      "values are modified to ensure that they lie on a grid point. "+
      "If an initial guess is specified on the command line, it is "+
      "used even if there is a good guess already in the table. "+
      "If the flag is not 10 or if \"recompute\" is true, then the EOS is "+
      "recomputed. If an EOS is loaded or the recompute was successful, "+
      "then the results are output to the screen. If the point was "+
      "successful it is stored in the current tables. If \"show_all"+
      "_nuclei\" is true, then a file named \"dist.o2\" is created "+
      "which holds the full nuclear distribution.",
      new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::point_nuclei),o2scl::cli::comm_option_both},
     {0,"increase-density",
      "Increase nB to optimize the phase transition.",7,7,
      "<nB low> <nB high> <Ye low> <Ye high> <T low> <T high> <output file>",
      "",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::increase_density),o2scl::cli::comm_option_both},
     {0,"fix-cc",
      "Increase nB to optimize the phase transition.",1,1,
      "<output file>",
      "",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::fix_cc),o2scl::cli::comm_option_both},
     {0,"verify","",1,4,"",
      "",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::verify),o2scl::cli::comm_option_both},
     {0,"select-high-T",
      "Choose the Skyrme model for the finite T corrections.",
      1,1,"<index>","",new o2scl::comm_option_mfptr<eos_nuclei>
      (this,&eos_nuclei::select_high_T_cl),o2scl::cli::comm_option_both}
    };
  cl.set_comm_option_vec(nopt,options);

  p_nB_grid_spec.str=&nB_grid_spec;
  p_nB_grid_spec.help=((string)"The function for default baryon ")+
    "density grid. Used by the new_table() function, and the "+
    "check-virial and eos-deriv commands.";
  cl.par_list.insert(make_pair("nB_grid_spec",&p_nB_grid_spec));

  p_Ye_grid_spec.str=&Ye_grid_spec;
  p_Ye_grid_spec.help=((string)"The function for default electron ")+
    "fraction grid. Used by the new_table() function, and the "+
    "check-virial and eos-deriv commands.";
  cl.par_list.insert(make_pair("Ye_grid_spec",&p_Ye_grid_spec));
  
  p_T_grid_spec.str=&T_grid_spec;
  p_T_grid_spec.help=((string)"The function for default temperature ")+
    "grid. Used by the new_table() function, and the "+
    "check-virial and eos-deriv commands.";
  cl.par_list.insert(make_pair("T_grid_spec",&p_T_grid_spec));
  
  p_show_all_nuclei.b=&show_all_nuclei;
  p_show_all_nuclei.help=((string)"If true, show all nuclei considered at ")+
    "every point (default false). This applies to the point-nuclei "+
    "command and the eos_vary_ZN() function.";
  cl.par_list.insert(make_pair("show_all_nuclei",&p_show_all_nuclei));
  
  p_fd_A_max.i=&fd_A_max;
  p_fd_A_max.help=((string)"The maximum value of A for the alg_mode=4 ")+
    "fixed distribution (default 200).";
  cl.par_list.insert(make_pair("fd_A_max",&p_fd_A_max));
  
  p_recompute.b=&recompute;
  p_recompute.help=((string)"If true, recompute points in the table ")+
    "and ignore the flag value (default false).";
  cl.par_list.insert(make_pair("recompute",&p_recompute));
  
  p_verify_only.b=&verify_only;
  p_verify_only.help=((string)"If true, don't call the solver, just ")+
    "verify the initial guess (default false).";
  cl.par_list.insert(make_pair("verify_only",&p_verify_only));
  
  p_edge_list.str=&edge_list;
  p_edge_list.help=((string)"List of names");
  cl.par_list.insert(make_pair("edge_list",&p_edge_list));
  
  p_six_neighbors.i=&six_neighbors;
  p_six_neighbors.help=((string)"If greater than 0, use nearest ")+
    "neighbors as guesses (default 0). Values greater than 0 use "+
    "the point at the next smallest density, values greater than 1 "+
    "use the point at the next largest density, values greater than 2 "+
    "use points at the next largest and next smallest temperature, "+
    "and values greater than 4 use the next largest and smallest "+
    "electron fraction.";
  cl.par_list.insert(make_pair("six_neighbors",&p_six_neighbors));
  
  p_rnuc_less_rws.b=&rnuc_less_rws;
  p_rnuc_less_rws.help=((string)"If true, restrict rnuc to ")+
    "be smaller than rws (default true).";
  cl.par_list.insert(make_pair("rnuc_less_rws",&p_rnuc_less_rws));

  p_propagate_points.b=&propagate_points;
  p_propagate_points.help=
    "If true, use previously computed points as guesses for neighbors";
  cl.par_list.insert(make_pair("propagate_points",&p_propagate_points));
  
  p_mh_tol_rel.d=&mh_tol_rel;
  p_mh_tol_rel.help=((std::string)"Relative tolerance for the solver ")+
    "in the eos_fixed_dist() function (default 10^(-6)).";
  cl.par_list.insert(make_pair("mh_tol_rel",&p_mh_tol_rel));
  
  p_max_time.d=&max_time;
  p_max_time.help=((string)"Maximum time, in seconds, for the ")+
    "generate-table command. The default is 0.0 which is "+
    "no maximum time.";
  cl.par_list.insert(make_pair("max_time",&p_max_time));
  
  p_nucleon_func.str=&nucleon_func;
  p_nucleon_func.help="Function for delta Z and delta N in the SNA.";
  cl.par_list.insert(make_pair("nucleon_func",&p_nucleon_func));

  p_ext_guess.str=&ext_guess;
  p_ext_guess.help=((string)"Filename containing separate table ")+
    "to use as a guess for the generate-table command.";
  cl.par_list.insert(make_pair("ext_guess",&p_ext_guess));

  p_Ye_list.str=&Ye_list;
  p_Ye_list.help=((string)"The list of electron fractions to consider ")+
    "for the generate-table command. Can be several comma-separated "+
    "ranges e.g. \"1-3,5-7,59-60\". Zero-indexed.";
  cl.par_list.insert(make_pair("Ye_list",&p_Ye_list));

  p_alg_mode.i=&alg_mode;
  p_alg_mode.help=((string)"Algorithm mode (default 4; ")+
    "0 for SNA, 1 for old SNA method, 2 for vary dist., "+
    "3 for old vary dist., and 4 for fixed dist.)";
  cl.par_list.insert(make_pair("alg_mode",&p_alg_mode));
  
  p_fixed_dist_alg.i=&fixed_dist_alg;
  p_fixed_dist_alg.help=((string)"Modify the algorithm for the ")+
    "eos_fixed_dist() function. The 1s digit is the number "+
    "of solves minus 1, the 10s digit is the number of brackets "+
    "divided by 10, the "+
    "100s digit is the number of minimizes, and the 1000s digit "+
    "is the number of random guesses to try divided by 1000. "+
    "The default is 1111. Other good options are 1319 and 9999.";
  cl.par_list.insert(make_pair("fixed_dist_alg",&p_fixed_dist_alg));
  
  p_function_verbose.i=&function_verbose;
  p_function_verbose.help=((string)"Verbose for individual functions ")+
    "(default value 11111).\n\t1s digit: fixed_ZN()\n\t10s digit: "+
    "vary_ZN()\n\t100s digit: "+
    "fixed_dist()\n\t1000s digit: vary_dist()\n\t10000s digit: "+
    "store_point().";
  cl.par_list.insert(make_pair("function_verbose",&p_function_verbose));
  
  p_max_ratio.d=&max_ratio;
  p_max_ratio.help="The maximum of N/Z or Z/N (default 2.25).";
  cl.par_list.insert(make_pair("max_ratio",&p_max_ratio));

  p_file_update_time.d=&file_update_time;
  p_file_update_time.help=((string)"The time (in seconds) between ")+
    "output file updates for the generate-table command. Default is "+
    "1800 or 30 minutes.";
  cl.par_list.insert(make_pair("file_update_time",&p_file_update_time));

  p_file_update_iters.i=&file_update_iters;
  p_file_update_iters.help=((string)"The number of points between ")+
    "output file updates for the generate-table command. Default is "+
    "100000.";
  cl.par_list.insert(make_pair("file_update_iters",&p_file_update_iters));

  p_include_detail.b=&include_detail;
  p_include_detail.help="";
  cl.par_list.insert(make_pair("include_detail",&p_include_detail));

}
