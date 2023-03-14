/*
  -------------------------------------------------------------------
  
  Copyright (C) 2022-2023, Andrew W. Steiner
  
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

#ifdef O2SCL_EIGEN
#include <eigen3/Eigen/Dense>
#endif

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

#ifdef O2SCL_NEVER_DEFINED

/** \brief
    
    Compute the free energy per baryon at (jnB,jYe,jT) assuming a
    correction centered at (inB,iYe,iT) given parameters in \c pars
*/
double eos_nuclei::F_interp
(size_t inB, size_t iYe, size_t iT, double r_inner,
 double r_outer, vector<size_t> nB_list,
 vector<size_t> Ye_list, vector<size_t> T_list,
 vector<double> F_list, size_t np, ubvector &pars,
 size_t i) {
    
  size_t jnB=nB_list[i];
  size_t jYe=Ye_list[i];
  size_t jT=T_list[i];
  
  double dnB=((double)jnB)-((double)inB);
  double dYe=((double)jYe)-((double)iYe);
  double dT=((double)jT)-((double)iT);
  
  double F_func=pars[0]+
    pars[1]*dnB+pars[2]*dnB*dnB+
    pars[3]*dYe+pars[4]*dYe*dYe+
    pars[5]*dT+pars[6]*dT*dT+
    pars[7]*dnB*dYe+pars[8]*dnB*dT+pars[9]*dYe*dT;
  vector<size_t> ix={jnB,jYe,jT};
  double F_orig=tg_F.get(ix);
  double rad=sqrt(dnB*dnB+dYe*dYe+dT*dT);
  double ret=(F_func-F_orig)/(1.0+exp(rad-r_inner)/(r_outer-r_inner)*8.0)+
    F_orig;
  
  return ret;
}

/** \brief Using parameters in \c p, compute the 
    relative deviations in \c f
*/
double eos_nuclei::interp_min
(size_t inB, size_t iYe, size_t iT, double r_inner,
 double r_outer, vector<size_t> nB_list,
 vector<size_t> Ye_list, vector<size_t> T_list,
 vector<double> F_list, size_t np, const vec_t &p) {

  for(size_t i=0;i<nd;i++) {
    double yi=F_new(np,p,i);
    f[i]=(yi-Flist[i])/0.1;
  }
  return;
}

#endif

int eos_nuclei::interp_point(std::vector<std::string> &sv,
                             bool itive_com) {

  if (sv.size()<5) {
    cerr << "Not enough arguments interp-point." << endl;
    return 1;
  }
  
  double nB_cent=o2scl::function_to_double(sv[1]);
  double Ye_cent=o2scl::function_to_double(sv[2]);
  double T_cent=o2scl::function_to_double(sv[3])/hc_mev_fm;

  if (!loaded) {
    O2SCL_ERR("No EOS loaded in interp_point.",o2scl::exc_einval);
  }
  if (with_leptons==false) {
    O2SCL_ERR("No EOS leptons in interp_point.",o2scl::exc_einval);
  }
  
  size_t inB=0, iYe=0, iT=0;
  inB=vector_lookup(n_nB2,nB_grid2,nB_cent);
  nB_cent=nB_grid2[inB];
  iYe=vector_lookup(n_Ye2,Ye_grid2,Ye_cent);
  Ye_cent=Ye_grid2[iYe];
  iT=vector_lookup(n_T2,T_grid2,T_cent*hc_mev_fm);
  T_cent=T_grid2[iT]/hc_mev_fm;
  cout << "Adjusted to grid, nB, Ye, T[MeV]: " << nB_cent << " "
       << Ye_cent << " " << T_cent*hc_mev_fm << endl;

  int window=o2scl::stoi(sv[4]);
  cout << "Using window size: " << window << endl;

  int count=0;
  for(int dnB=-window;dnB<=window;dnB++) {
    for(int dYe=-window;dYe<=window;dYe++) {
      for(int dT=-window;dT<=window;dT++) {
        if (abs(dnB)+abs(dYe)+abs(dT)<=window) {
          count++;
        }
      }
    }
  }
  
  cout << "Using " << count << " points to interpolate." << endl;
  
  ubmatrix ix(count,3);
  ubmatrix iy(1,count);

  count=0;
  for(int dnB=-window;dnB<=window;dnB++) {
    for(int dYe=-window;dYe<=window;dYe++) {
      for(int dT=-window;dT<=window;dT++) {
        
        if (abs(dnB)+abs(dYe)+abs(dT)<=window &&
            inB+dnB>=0 && inB+dnB<n_nB2 &&
            iYe+dYe>=0 && iYe+dYe<n_Ye2 &&
            iT+dT>=0 && iT+dT<n_T2) {
          
          ix(count,0)=inB+dnB;
          ix(count,1)=iYe+dYe;
          ix(count,2)=iT+dT;
          vector<size_t> index={inB+dnB,iYe+dYe,iT+dT};
          iy(0,count)=tg_F.get(index);

          if (true) {
            cout << ix(count,0) << " " << ix(count,1) << " "
                 << ix(count,2) << " " << iy(0,count) << endl;
          }
          
          count++;
          
        }
        
      }
    }
  }

  ubmatrix ix2=ix;

  interpm_krige_eos ike;
  ike.mode=ike.mode_loo_cv_bf;
  ike.full_min=true;
  ike.def_mmin.verbose=1;
  vector<mcovar_funct_rbf_noise> mfr(1);
  mfr[0].len.resize(3);
  
  //vector<double> len_list={0.1,0.3,0.5,0.75,1.0,2.0,3.0,4.0,5.0,6.0,
  //8.0,10.0};
  vector<double> len_list={2.0,3.0};
  vector<double> l10_list={-15,-13,-11,-9};  
  vector<vector<double>> ptemp;
  ptemp.push_back(len_list);
  ptemp.push_back(len_list);
  ptemp.push_back(len_list);
  ptemp.push_back(l10_list);
  vector<vector<vector<double>>> param_lists;
  param_lists.push_back(ptemp);
  cout << "Going to set_covar()." << endl;
  ike.set_covar(mfr,param_lists);

  cout << "Going to set_data()." << endl;
  ike.verbose=2;
  ike.set_data(3,1,count,ix,iy);

  cout << "Herex." << endl;
  
  for(int k=0;k<count;k++) {

    cout << "k: " << k << endl;
    
    ix2(k,0)=inB;
    ix2(k,1)=iYe;
    ix2(k,2)=iT;
    
    vector<size_t> index={((size_t)(ix2(k,0)*1.0001)),
      ((size_t)(ix2(k,1)*1.0001)),((size_t)(ix2(k,2)*1.0001))};
    vector<double> point={((double)inB),((double)iYe),((double)iT)};
    vector<double> out(1);

    double nB=nB_grid2[ix2(k,0)];
    double Ye=Ye_grid2[ix2(k,1)];
    double T_MeV=T_grid2[ix2(k,2)];

    cout << nB_cent << " " << Ye_cent << " " << T_cent*hc_mev_fm << endl;
    cout << inB << " " << iYe << " " << iT << endl;
    cout << nB << " " << Ye << " " << T_MeV << endl;
    
    // Derivatives of the physical coordinates with respect to the indices
    
    double dnBdi=2.0*0.04*log(10.0)*pow(10.0,((double)inB)*0.04-12.0);
    double dYedj=0.01;
    double dTdk=0.1*log(1.046)*pow(1.046,iT);

    // Evaluate the free energy and its derivatives analytically
    // using the interpolator
    
    ike.eval(index,out);
    double Fintp=out[0];
    
    ike.deriv(index,out,0);
    double dF_dnB=out[0]/hc_mev_fm/dnBdi;
    ike.deriv(index,out,1);
    double dF_dYe=out[0]/hc_mev_fm/dYedj;
    ike.deriv(index,out,2);
    // No hbarc here b/c dTdk has units of MeV as does out[0]
    double dF_dT=out[0]/dTdk;

    ike.deriv2(index,out,0,0);
    double F_nBnB=out[0]/hc_mev_fm;
    ike.deriv2(index,out,0,1);
    double F_nBYe=out[0]/hc_mev_fm;
    ike.deriv2(index,out,1,1);
    double F_YeYe=out[0]/hc_mev_fm;
    ike.deriv2(index,out,0,2);
    double F_nBT=out[0]/hc_mev_fm;
    ike.deriv2(index,out,1,2);
    double F_YeT=out[0]/hc_mev_fm;
    ike.deriv2(index,out,2,2);
    double F_TT=out[0]/hc_mev_fm;
    
    // Use those derivatives to compute the chemical potentials and
    // the entropy density
    
    double mun=Fintp/hc_mev_fm-Ye*dF_dYe+nB*dF_dnB;
    double mue=tg_mue.get(index)/hc_mev_fm;
    double mup=Fintp/hc_mev_fm+(1.0-Ye)*dF_dYe+nB*dF_dnB-mue;
    double en=-nB*dF_dT;

    // Compare theose derivatives with the stored values

    cout << "Stored  : mun[MeV],mup[MeV],mue[MeV],S: ";
    cout << tg_mun.get(index) << " ";
    cout << tg_mup.get(index) << " ";
    cout << tg_mue.get(index) << " ";
    cout << tg_S.get(index) << endl;
    cout << "Computed: mun[MeV],mup[MeV],mue[MeV],S: " << mun*hc_mev_fm << " "
         << mup*hc_mev_fm << " " << mue*hc_mev_fm << " ";
    cout << en/nB << endl;

    // Now compute the second derivatives and the speed of sound
    
    double f_nnnn=(Ye*Ye*F_YeYe+nB*(2.0*dF_dnB-2.0*Ye*F_nBYe+nB*F_nBnB))/nB;
    double f_nnnp=((Ye-1.0)*Ye*F_YeYe+nB*(2.0*dF_dnB+(1.0-2.0*Ye)*
                                          F_nBYe+nB*F_nBnB))/nB;
    double f_npnp=((Ye-1.0)*(Ye-1.0)*F_YeYe+nB*(2.0*dF_dnB-2.0*(Ye-1.0)*
                                                F_nBYe+nB*F_nBnB))/nB;
    double f_nnT=dF_dT-Ye*F_YeT+nB*F_nBT;
    double f_npT=dF_dT-(Ye-1.0)*F_YeT+nB*F_nBT;
    double f_TT=nB*F_TT;
    
    double den=en*T_MeV/hc_mev_fm+(mun+neutron.m)*nB*(1.0-Ye)+
      (mup+proton.m)*nB*Ye+mue*nB*Ye;
    double nn2=nB*(1.0-Ye);
    double np2=nB*Ye;
    double cs_sq=(nn2*nn2*(f_nnnn-f_nnT*f_nnT/f_TT)+
                  2.0*nn2*np2*(f_nnnp-f_nnT*f_npT/f_TT)+
                  np2*np2*(f_npnp-f_npT*f_npT/f_TT)-
                  2.0*en*(nn2*f_nnT/f_TT+np2*f_npT/f_TT)-en*en/f_TT)/den;

    cout << "Here: " << cs_sq << endl;

    // Also compute dPdnB

    vector<size_t> im1={index[0]-1,index[1],index[2]};
    
    double dmun_dnB=2.0*dF_dnB-Ye*F_nBYe+nB*F_nBnB;
    double dmup_dnB=2.0*dF_dnB-(Ye-1.0)*F_nBYe+nB*F_nBnB;
    cout << dmun_dnB*hc_mev_fm << " " << dmup_dnB*hc_mev_fm << endl;
    cout << (tg_mun.get(index)-tg_mun.get(im1))/hc_mev_fm/
      (nB_grid2[index[0]]-nB_grid2[index[0]-1]) << " ";
    cout << (tg_mup.get(index)-tg_mup.get(im1))/hc_mev_fm/
      (nB_grid2[index[0]]-nB_grid2[index[0]-1]) << endl;
    
    double dmuden_dnB=((1.0-Ye)*mun+nn2*dmun_dnB+
                       Ye*(mup+mue)+np2*dmup_dnB);
    double dPdnB=dmuden_dnB-dF_dnB*nB-Fintp/hc_mev_fm;
    cout << dPdnB << endl;
    cout << (tg_P.get(index)-tg_P.get(im1))/hc_mev_fm/
      (nB_grid2[index[0]]-nB_grid2[index[0]-1]) << endl;
    
    exit(-1);
  }

  return 0;
}

void interpm_krige_eos::set(std::vector<double> &nB_grid2,
                            std::vector<double> &Ye_grid2,
                            std::vector<double> &T_grid2,
                            o2scl::tensor_grid<> &tg_F,
                            o2scl::tensor_grid<> &tg_P,
                            o2scl::tensor_grid<> &tg_S,
                            o2scl::tensor_grid<> &tg_mun,
                            o2scl::tensor_grid<> &tg_mup,
                            o2scl::tensor_grid<> &tg_mue,
                            double mn, double mpx,
                            int window) {

  mneut=mn;
  mprot=mpx;
  
  nB_grid=nB_grid2;
  Ye_grid=Ye_grid2;
  T_grid=T_grid2;

  tgp_F=&tg_F;
  tgp_P=&tg_P;
  tgp_S=&tg_S;
  tgp_mun=&tg_mun;
  tgp_mup=&tg_mup;
  tgp_mue=&tg_mue;
    
  size_t n_nB=nB_grid2.size();
  size_t n_Ye=Ye_grid2.size();
  size_t n_T=T_grid2.size();

  ubmatrix ix, iy;

  int count;
    
  for(size_t irun=0;irun<2;irun++) {
      
    count=0;
    
    for(int ilist=0;ilist<((int)index_list.size());ilist++) {
        
      int inB=index_list[ilist*3];
      int iYe=index_list[ilist*3+1];
      int iT=index_list[ilist*3+2];
        
      for(int dnB=-window;dnB<=window;dnB++) {
        for(int dYe=-window;dYe<=window;dYe++) {
          for(int dT=-window;dT<=window;dT++) {
            if (inB+dnB>=0 && iYe+dYe>=0 && iT+dT>=0 &&
                inB+dnB<((int)n_nB) && iYe+dYe<((int)n_Ye) &&
                iT+dT<((int)n_T) && abs(dnB)+abs(dYe)+abs(dT)<=window) {

              if (irun==1) {
                ix(count,0)=inB+dnB;
                ix(count,1)=iYe+dYe;
                ix(count,2)=iT+dT;
                std::vector<size_t> index={inB+dnB,iYe+dYe,iT+dT};
                iy(0,count)=tg_F.get(index);
                  
                if (true) {
                  std::cout << ix(count,0) << " " << ix(count,1) << " "
                            << ix(count,2) << " " << iy(0,count) << std::endl;
                }
              }

              count++;
            }
          }
        }
      }
    }

    if (irun==0) {
        
      std::cout << "Using " << count << " points to interpolate."
                << std::endl;
        
      ix.resize(count,3);
      iy.resize(1,count);
    }
      
  }
    
  mode=mode_loo_cv_bf;
  full_min=true;
  def_mmin.verbose=1;
  std::vector<o2scl::mcovar_funct_rbf_noise> mfr(1);
  mfr[0].len.resize(3);
    
  //vector<double> len_list={0.1,0.3,0.5,0.75,1.0,2.0,3.0,4.0,5.0,6.0,
  //8.0,10.0};
  std::vector<double> len_list={2.0,3.0};
  std::vector<double> l10_list={-15,-13,-11,-9};  
  std::vector<std::vector<double>> ptemp;
  ptemp.push_back(len_list);
  ptemp.push_back(len_list);
  ptemp.push_back(len_list);
  ptemp.push_back(l10_list);
  std::vector<std::vector<std::vector<double>>> param_lists;
  param_lists.push_back(ptemp);
  std::cout << "Going to set_covar()." << std::endl;
  set_covar(mfr,param_lists);
    
  std::cout << "Going to set_data()." << std::endl;
  verbose=2;
  set_data(3,1,count,ix,iy);
    
  return;
}
  
int interpm_krige_eos::addl_const(double &ret) {

  ret=0.0;
  bool compare=true;
    
  for(int ilist=0;ilist<index_list.size();ilist++) {
      
    std::cout << "ilist: " << ilist << std::endl;
      
    size_t inB=index_list[ilist*3];
    size_t iYe=index_list[ilist*3+1];
    size_t iT=index_list[ilist*3+2];
      
    std::vector<size_t> index={((size_t)inB),((size_t)iYe),((size_t)iT)};
        
    double nB=nB_grid[inB];
    double Ye=Ye_grid[iYe];
    double T_MeV=T_grid[iT];

    // Derivatives of the physical coordinates with respect to the indices
    
    double dnBdi=2.0*0.04*log(10.0)*pow(10.0,((double)inB)*0.04-12.0);
    double dYedj=0.01;
    double dTdk=0.1*log(1.046)*pow(1.046,iT);

    // Evaluate the free energy and its derivatives analytically
    // using the interpolator
        
    std::vector<double> out(1);
    eval(index,out);
    double Fintp=out[0];
        
    deriv(index,out,0);
    double dF_dnB=out[0]/hc_mev_fm/dnBdi;
    deriv(index,out,1);
    double dF_dYe=out[0]/hc_mev_fm/dYedj;
    deriv(index,out,2);
    // No hbarc here b/c dTdk has units of MeV as does out[0]
    double dF_dT=out[0]/dTdk;
        
    deriv2(index,out,0,0);
    double F_nBnB=out[0]/hc_mev_fm;
    deriv2(index,out,0,1);
    double F_nBYe=out[0]/hc_mev_fm;
    deriv2(index,out,1,1);
    double F_YeYe=out[0]/hc_mev_fm;
    deriv2(index,out,0,2);
    double F_nBT=out[0]/hc_mev_fm;
    deriv2(index,out,1,2);
    double F_YeT=out[0]/hc_mev_fm;
    deriv2(index,out,2,2);
    double F_TT=out[0]/hc_mev_fm;
        
    // Use those derivatives to compute the chemical potentials and
    // the entropy density
        
    double mun=Fintp/hc_mev_fm-Ye*dF_dYe+nB*dF_dnB;
    double mue=tgp_mue->get(index)/hc_mev_fm;
    double mup=Fintp/hc_mev_fm+(1.0-Ye)*dF_dYe+nB*dF_dnB-mue;
    double en=-nB*dF_dT;
        
    // Compare theose derivatives with the stored values

    if (compare) {
      std::cout << "Stored  : mun[MeV],mup[MeV],mue[MeV],S: ";
      std::cout << tgp_mun->get(index) << " ";
      std::cout << tgp_mup->get(index) << " ";
      std::cout << tgp_mue->get(index) << " ";
      std::cout << tgp_S->get(index) << std::endl;
      std::cout << "Computed: mun[MeV],mup[MeV],mue[MeV],S: "
                << mun*hc_mev_fm << " " << mup*hc_mev_fm << " "
                << mue*hc_mev_fm << " ";
      std::cout << en/nB << std::endl;
    }
        
    // Now compute the second derivatives and the speed of sound
        
    double f_nnnn=(Ye*Ye*F_YeYe+nB*(2.0*dF_dnB-2.0*Ye*F_nBYe+nB*F_nBnB))/nB;
    double f_nnnp=((Ye-1.0)*Ye*F_YeYe+nB*(2.0*dF_dnB+(1.0-2.0*Ye)*
                                          F_nBYe+nB*F_nBnB))/nB;
    double f_npnp=((Ye-1.0)*(Ye-1.0)*F_YeYe+nB*(2.0*dF_dnB-2.0*(Ye-1.0)*
                                                F_nBYe+nB*F_nBnB))/nB;
    double f_nnT=dF_dT-Ye*F_YeT+nB*F_nBT;
    double f_npT=dF_dT-(Ye-1.0)*F_YeT+nB*F_nBT;
    double f_TT=nB*F_TT;
    
    double den=en*T_MeV/hc_mev_fm+(mun+mneut)*nB*(1.0-Ye)+
      (mup+mprot)*nB*Ye+mue*nB*Ye;
    double nn2=nB*(1.0-Ye);
    double np2=nB*Ye;
    double cs_sq=(nn2*nn2*(f_nnnn-f_nnT*f_nnT/f_TT)+
                  2.0*nn2*np2*(f_nnnp-f_nnT*f_npT/f_TT)+
                  np2*np2*(f_npnp-f_npT*f_npT/f_TT)-
                  2.0*en*(nn2*f_nnT/f_TT+np2*f_npT/f_TT)-en*en/f_TT)/den;
        
    std::cout << "Here: " << cs_sq << std::endl;
    if (cs_sq>1.0 || cs_sq<0.0 || !std::isfinite(cs_sq)) {
      return 1;
    }
        
    // Also compute dPdnB
        
    std::vector<size_t> im1={index[0]-1,index[1],index[2]};
        
    double dmun_dnB=2.0*dF_dnB-Ye*F_nBYe+nB*F_nBnB;
    double dmup_dnB=2.0*dF_dnB-(Ye-1.0)*F_nBYe+nB*F_nBnB;
    if (compare) {
      std::cout << dmun_dnB*hc_mev_fm << " " << dmup_dnB*hc_mev_fm << std::endl;
      std::cout << (tgp_mun->get(index)-tgp_mun->get(im1))/hc_mev_fm/
        (nB_grid[index[0]]-nB_grid[index[0]-1]) << " ";
      std::cout << (tgp_mup->get(index)-tgp_mup->get(im1))/hc_mev_fm/
        (nB_grid[index[0]]-nB_grid[index[0]-1]) << std::endl;
    }
        
    double dmuden_dnB=((1.0-Ye)*mun+nn2*dmun_dnB+
                       Ye*(mup+mue)+np2*dmup_dnB);
    double dPdnB=dmuden_dnB-dF_dnB*nB-Fintp/hc_mev_fm;
    if (compare) {
      std::cout << dPdnB << std::endl;
      std::cout << (tgp_P->get(index)-tgp_P->get(im1))/hc_mev_fm/
        (nB_grid[index[0]]-nB_grid[index[0]-1]) << std::endl;
    }
    if (dPdnB<=0.0) return 2;
      
  }
      
        
  return 0;
}
