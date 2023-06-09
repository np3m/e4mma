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
  if (sv.size()<6) {
    cerr << "Not enough arguments interp-point." << endl;
    return 1;
  }

  double nB_cent=o2scl::function_to_double(sv[1]);
  double Ye_cent=o2scl::function_to_double(sv[2]);
  double T_cent=o2scl::function_to_double(sv[3])/hc_mev_fm;
  
  int window=o2scl::stoi(sv[4]);

  std::string st_o2=sv[5];
  hdf_file hff;
  o2scl::tensor_grid<> tgp_cs2;
  hff.open(st_o2);
  hdf_input(hff, tgp_cs2);

  eos_nuclei::interpolate(nB_cent, Ye_cent, T_cent, window, st_o2, tgp_cs2, itive_com);
  
  hff.close();
  return 0;
}


void eos_nuclei::interpolate(double nB_p,
                                            double Ye_p,
                                            double T_p,
                                            int window,
                                            std::string st_o2,
                                            o2scl::tensor_grid<> &tg_cs2,
                                            bool itive_com) {

  if (!loaded) {
    O2SCL_ERR("No EOS loaded in interp_point.",o2scl::exc_einval);
  }
  if (with_leptons==false) {
    O2SCL_ERR("No EOS leptons in interp_point.",o2scl::exc_einval);
  } 

//  std::vector<double> results;
//  results.assign(0.0, 5);
  double nB_cent=nB_p;
  double Ye_cent=Ye_p;
  double T_cent=T_p/hc_mev_fm;
  size_t inB=0, iYe=0, iT=0;
  inB=vector_lookup(n_nB2,nB_grid2,nB_cent);
  nB_cent=nB_grid2[inB];
  iYe=vector_lookup(n_Ye2,Ye_grid2,Ye_cent);
  Ye_cent=Ye_grid2[iYe];
  iT=vector_lookup(n_T2,T_grid2,T_cent*hc_mev_fm);
  T_cent=T_grid2[iT]/hc_mev_fm;
  cout << "Adjusted to grid, nB, Ye, T[MeV]: " << nB_cent << " "
       << Ye_cent << " " << T_cent*hc_mev_fm << endl;
  cout << "At grid point: " << inB << " " << iYe << " " << iT << endl;

  cout << "Using window size: " << window << endl;

  // Create interpolation object
  interpm_krige_eos ike;
  ike.mode=ike.mode_loo_cv_bf;
  ike.full_min=true;
  ike.def_mmin.verbose=1;
  
  /// Load cs2 from a file
  hdf_file hff;
  hff.open(st_o2);
  hdf_input(hff,ike.tgp_cs2);
  hff.close();

  for(int dnB=-window;dnB<=window;dnB++) {
    for(int dYe=-window;dYe<=window;dYe++) {
      for(int dT=-window;dT<=window;dT++) {
        vector<size_t> index={inB+dnB,iYe+dYe,iT+dT};
        if (abs(dnB)+abs(dYe)+abs(dT)<=window &&
            inB+dnB>=0 && inB+dnB<n_nB2 &&
            iYe+dYe>=0 && iYe+dYe<n_Ye2 &&
            iT+dT>=0 && iT+dT<n_T2) {
          if (ike.tgp_cs2.get_rank()>=3 &&
              (ike.tgp_cs2.get(index)>1.0 ||
               !std::isfinite(ike.tgp_cs2.get(index)) ||
               ike.tgp_cs2.get(index)<0.0)) {
            ike.fix_list.push_back(index[0]);
            ike.fix_list.push_back(index[1]);
            ike.fix_list.push_back(index[2]);
            //cout << "fix: " << index[0] << " " << index[1] << " "
            //<< index[2] << endl;
          } else {
            ike.calib_list.push_back(index[0]);
            ike.calib_list.push_back(index[1]);
            ike.calib_list.push_back(index[2]);
            //cout << "cal: " << index[0] << " " << index[1] << " "
            //<< index[2] << endl;
          }
        }
      }
    }
  }
  /*
  if (ike.fix_list.size()==0) {
    ike.fix_list.push_back(nB_p);
    ike.fix_list.push_back(Ye_p);
    ike.fix_list.push_back(T_p);
  }*/

  cout << "Using " << ike.calib_list.size()/3 << " points to calibrate"
       << endl;
  cout << "  and attempting to fix " << ike.fix_list.size()/3 << " points."
       << endl;
  size_t count=ike.calib_list.size()/3;

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
  vector<mcovar_funct_rbf_noise> mfr(1);
  mfr[0].len.resize(3);
  ike.set_covar(mfr,param_lists);
 
  ike.skip_optim=true;
  ike.set(nB_grid2,Ye_grid2,T_grid2,tg_F,tg_P,tg_S,
          tg_mun,tg_mup,tg_mue,neutron.m,proton.m);


  double min_qual=1.0e99;
  vector<double> pnt(4), min_p;  
  for(pnt[0]=2.0;pnt[0]<20.0;pnt[0]*=1.5) {
    for(pnt[1]=2.0;pnt[1]<20.0;pnt[1]*=1.5) {
      for(pnt[2]=2.0;pnt[2]<20.0;pnt[2]*=1.5) {
        for(pnt[3]=-15.0;pnt[3]<-8.99;pnt[3]+=2.0) {
          vector_out(cout,pnt,true);
          (*ike.cf)[0].set_params(pnt);
          int success;
          double q=ike.qual_fun(0,success);
          cout << q << " " << success << endl;
          if (q<min_qual) {
            min_p=pnt;
          }
        }
      }
    }
  }
  vector_out(cout,min_p,true);

  // Use the interpolation results to fix points 
  std::vector<double> out(1);
  for(size_t j=0;j<ike.fix_list.size();j+=3) {

    vector<size_t> index={ike.fix_list[j],ike.fix_list[j+1],
      ike.fix_list[j+2]};
    
    double nB=nB_grid2[inB];
    double Ye=Ye_grid2[iYe];
    double T_MeV=T_grid2[iT];
    
    // Derivatives of the physical coordinates with respect to the indices
    
    double dnBdi=2.0*0.04*log(10.0)*pow(10.0,((double)inB)*0.04-12.0);
    double dYedj=0.01;
    double dTdk=0.1*log(1.046)*pow(1.046,iT);
    
    ike.eval(index,out);
    double Fintp=out[0];
    tg_F.get(index)=Fintp;
    
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

    double mun=Fintp/hc_mev_fm-Ye*dF_dYe+nB*dF_dnB;
    double mue=tg_mue.get(index)/hc_mev_fm;
    double mup=Fintp/hc_mev_fm+(1.0-Ye)*dF_dYe+nB*dF_dnB-mue;
    double en=-nB*dF_dT;
    tg_mun.get(index)=mun;
    tg_mup.get(index)=mup;
    tg_mue.get(index)=mue;
    
    tg_S.get(index)=en/nB;

    // unverified
    tg_E.get(index)=tg_F.get(index)+T_MeV*tg_S.get(index);
    tg_P.get(index)=tg_F.get(index)+mun*neutron.m+mup*proton.m+
      mue*electron.m;
//  code for speed of sound squared taken from main/eos_interp.cpp
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
    tg_cs2.get(index)=cs_sq;

    // Attempt at calculating dP/dnB
    // code taken from main/eos_interp.cpp
    cout << "starting to calc dP/dnB\n";
    vector<size_t> im1={index[0]-1,index[1],index[2]};
    double dmundnB=(f_nnnn*(1-Ye))+(f_nnnp*Ye);
    double dmupmuednB=(f_nnnp*(1-Ye))+(f_npnp*Ye);
    //double dmundnB=F_nBnB-(Ye*(((1/nB)*F_nBYe)-((1/(nB*nB))*dF_dYe)));
    //double dmupmuednB=F_nBnB-((1-Ye)*(((1/nB)*F_nBYe)-((1/(nB*nB))*dF_dYe)));
    double dPdnB=(dmundnB*nB*(1-Ye))+(mun*(1-Ye))+(dmupmuednB*Ye*nB)+(Ye*(mup+mue))-((mun*(1-Ye))+((mup+mue)*Ye));

    cout << "Calculated: mun[MeV],mup[MeV],mue[MeV]: " << mun << " " << mup << " " << mue << " \n";
    cout << "dmundnB: " << dmundnB << endl;
    cout << "dmupmuednB: " << dmupmuednB << endl;
    cout << "dmundnB*hc: " << dmundnB*hc_mev_fm << endl;
    cout << "dmupmuednB*hc: " << dmupmuednB*hc_mev_fm << endl;
    cout << "dPdnB: " << dPdnB << endl;
    if (!(index[0]==0)){
        cout << "dmundnB from table: " << (tg_mun.get(index)-tg_mun.get(im1))/hc_mev_fm/
        (nB_grid2[index[0]]-nB_grid2[index[0]-1]) << " ";
        cout << "dmupmuednB from table: " << ((tg_mup.get(index)+tg_mue.get(index))-(tg_mup.get(im1)+tg_mue.get(im1)))/hc_mev_fm/
        (nB_grid2[index[0]]-nB_grid2[index[0]-1]) << endl;
        cout << "dPdnB from table: " << (tg_P.get(index)-tg_P.get(im1))/hc_mev_fm/
      (nB_grid2[index[0]]-nB_grid2[index[0]-1]) << endl;
    }
    if (dPdnB>0.0) {
        cout<<"success\n";
    }
    else {
        cout<<"failure\n";
    }
  }
}

int eos_nuclei::interp_file(std::vector<std::string> &sv,
                            bool itive_com) {
  if (sv.size()<5) {
    cerr << "Not enough arguments interp-file." << endl;
    return 1;
  }
  
  double results=0.0;
//  std::vector<double> results;
  std::vector<double> point;
  point.assign(0.0, 3);

  std::string csv_path = sv[1];
  int window=o2scl::stoi(sv[2]);
  std::string st_o2=sv[3];
  std::string stfix_o2=sv[4];
  hdf_file hff;
  hdf_file hff2;
  o2scl::tensor_grid<> tgp_cs2;
  o2scl::tensor_grid<> tgp_cs2_old;
  hff.open(st_o2);
  hdf_input(hff, tgp_cs2);
  hdf_input(hff, tgp_cs2_old);
  std::filesystem::copy_file(st_o2, stfix_o2, std::filesystem::copy_options::overwrite_existing);
  hff2.open(stfix_o2);
  std::ifstream csvfile;

  std::string line;
  int col = 0;
  int x;
  std::string val;
  while (std::getline(csvfile, line)) {
      std::stringstream s(line);
      if (true) {
          x=0;
          while (s>>val) {
              point[x] = std::stod(val);
              x++;
              if (s.peek()==',') {
                  s.ignore();
              }
          }
        eos_nuclei::interpolate(point[0], point[1], point[2], window, st_o2, tgp_cs2, itive_com);
      }
      else {
          cout<<"csv file of superliminal points must have 3 terms in each line";
      }
  }

  //solution for reseting file to beginning taken from stack exchange.
  csvfile.clear();
  csvfile.seekg(0, ios::beg);
  double fixed = 0.0;
  col = 0;
  while (std::getline(csvfile, line)) {
      std::stringstream s(line);
      if (true) {
          x=0;
          while (s>>val) {
              point[x] = std::stod(val);
              x++;
              if (s.peek()==',') {
                  s.ignore();
              }
          }
        //fixed= some function of tgp_cs2.get(point) and tgp_cs2_old.get(point);
        tgp_cs2.get(point)=fixed;
      }
      else {
          cout<<"csv file of superliminal points must have 3 terms in each line";
      }
  }
  hdf_output(hff2,tgp_cs2,"cs_sq");
  hff.close();
  hff2.close();
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
                            double mn, double mpx) {

  ubmatrix ix(calib_list.size()/3,3);
  ubmatrix iy(1,calib_list.size()/3);

  for(size_t j=0;j<calib_list.size();j+=3) {
    ix(j/3,0)=calib_list[j];
    ix(j/3,1)=calib_list[j+1];
    ix(j/3,2)=calib_list[j+2];
    vector<size_t> index={calib_list[j],calib_list[j+1],
      calib_list[j+2]};
    iy(0,j/3)=tg_F.get(index);
    if (true) {
      cout << ix(j/3,0) << " " << ix(j/3,1) << " "
           << ix(j/3,2) << " " << iy(0,j/3) << endl;
    }
  }

  // Make a copy because interpm_krige_eos will keep it
  ubmatrix ix2=ix;

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
    
  std::cout << "Going to set_data()." << std::endl;
  verbose=2;
  set_data(3,1,calib_list.size()/3,ix,iy);

  return;
}
  
int interpm_krige_eos::addl_const(size_t iout, double &ret) {

  // First, we need to ensure the interpolator has been
  // setup to be able to use the eval() and deriv()
  // functions
  
  // Select the row of the data matrix
  mat_y_row_t yiout2(this->y,iout);
  
  // Construct the KXX matrix
  size_t size=this->x.size1();
  
  mat_inv_kxx_t KXX(size,size);
  for(size_t irow=0;irow<size;irow++) {
    mat_x_row_t xrow(this->x,irow);
    for(size_t icol=0;icol<size;icol++) {
      mat_x_row_t xcol(this->x,icol);
      if (irow>icol) {
        KXX(irow,icol)=KXX(icol,irow);
      } else {
        KXX(irow,icol)=(*cf)[iout](xrow,xcol);
      }
    }
  }
  
  if (verbose>2) {
    std::cout << "Done creating covariance matrix with size "
              << size << std::endl;
  }
  
  // Perform the matrix inversion and compute the determinant
  
  double lndet;
  
  // Construct the inverse of KXX
  if (verbose>2) {
    std::cout << "Performing matrix inversion with size "
              << size << std::endl;
  }
  this->inv_KXX[iout].resize(size,size);
  int cret=this->mi.invert_det(size,KXX,this->inv_KXX[iout],lndet);
  if (cret!=0) {
    cout << "Return failed inversion." << endl;
    return 3;
  }
  
  lndet=log(lndet);
  
  if (verbose>2) {
    std::cout << "Done performing matrix inversion with size "
              << size << std::endl;
  }
        
  // Inverse covariance matrix times function vector
  this->Kinvf[iout].resize(size);
  o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                     o2scl_cblas::o2cblas_NoTrans,
                     size,size,1.0,this->inv_KXX[iout],
                     yiout2,0.0,this->Kinvf[iout]);

  // Done setting interpolator
  // ----------
  
  ret=0.0;
  bool compare=true;

  cout << "i nB Ye T cs2_itp cs2_tab d2FdnB2_itp d2FdnB2_tab "
       << "d2FdYe2_itp d2FdYe2_tab "
       << "dPdnB_itp dPdnB_tab1 dPdnB_tab2" << endl;
  for(size_t ilist=0;ilist<(calib_list.size()+fix_list.size())/3;
      ilist++) {
      
    std::cout << ilist << " ";

    size_t inB, iYe, iT;
    if (ilist<calib_list.size()/3) {
      inB=calib_list[ilist*3];
      iYe=calib_list[ilist*3+1];
      iT=calib_list[ilist*3+2];
    } else {
      inB=fix_list[(ilist-calib_list.size()/3)*3];
      iYe=fix_list[(ilist-calib_list.size()/3)*3+1];
      iT=fix_list[(ilist-calib_list.size()/3)*3+2];
    }
    //inB=8;
    //iYe=48;
    //iT=1;
      
    std::vector<size_t> index={((size_t)inB),((size_t)iYe),((size_t)iT)};
        
    double nB=nB_grid[inB];
    double Ye=Ye_grid[iYe];
    double T_MeV=T_grid[iT];

    cout << nB << " " << Ye << " " << T_MeV << " ";
    
    // Derivatives of the physical coordinates with respect to the indices
    
    double dnBdi=2.0*0.04*log(10.0)*pow(10.0,((double)inB)*0.04-12.0);
    double dYedj=0.01;
    double dTdk=0.1*log(1.046)*pow(1.046,iT);

    double didnB=25.0/nB/log(10.0);
    double d2idnB2=-25.0/nB/nB/log(10.0);
    double djdYe=100.0;
    double d2jdYe2=0.0;
    double dkdT=1.0/T_MeV/log(1.046);
    double d2kdT2=-1.0/T_MeV/T_MeV/log(1.046);
  
    // Evaluate the free energy and its derivatives analytically
    // using the interpolator

    std::vector<double> out(1);
    eval(index,out);
    double Fintp=out[0];
        
    deriv(index,out,0);
    double dFdi=out[0]/hc_mev_fm;
    double dF_dnB=dFdi*didnB;
    deriv(index,out,1);
    double dFdj=out[0]/hc_mev_fm;
    double dF_dYe=dFdj*djdYe;
    deriv(index,out,2);
    double dFdk=out[0]/hc_mev_fm;
    double dF_dT=dFdk*dkdT*hc_mev_fm;
        
    deriv2(index,out,0,0);
    double d2Fdi2=out[0]/hc_mev_fm;
    double F_nBnB=d2Fdi2*didnB*didnB+dFdi*d2idnB2;
    
    deriv2(index,out,0,1);
    double d2Fdidj=out[0]/hc_mev_fm;
    double F_nBYe=d2Fdidj*didnB*djdYe;
    
    deriv2(index,out,1,1);
    double d2Fdj2=out[0]/hc_mev_fm;
    double F_YeYe=d2Fdj2*djdYe*djdYe+dFdj*d2jdYe2;
    
    deriv2(index,out,0,2);
    double d2Fdidk=out[0]/hc_mev_fm;
    double F_nBT=d2Fdidk*didnB*dkdT*hc_mev_fm;
    
    deriv2(index,out,1,2);
    double d2Fdjdk=out[0]/hc_mev_fm;
    double F_YeT=d2Fdjdk*djdYe*dkdT*hc_mev_fm;
    
    deriv2(index,out,2,2);
    double d2Fdk2=out[0]/hc_mev_fm;
    double F_TT=(d2Fdk2*dkdT*dkdT+dFdk*d2kdT2)*hc_mev_fm*hc_mev_fm;

    // Use those derivatives to compute the chemical potentials and
    // the entropy density

    double mun=Fintp/hc_mev_fm-Ye*dF_dYe+nB*dF_dnB;
    double mue=tgp_mue->get(index)/hc_mev_fm;
    double mup=Fintp/hc_mev_fm+(1.0-Ye)*dF_dYe+nB*dF_dnB-mue;
    double en=-nB*dF_dT;
        
    // Compare theose derivatives with the stored values

    if (false && compare) {
      std::cout << "Indices: " << index[0] << " " << index[1] << " "
                << index[2] << std::endl;
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

    cout.setf(ios::showpos);
    std::cout << cs_sq << " ";
    cout << tgp_cs2.get(index) << " ";
        
    // Also compute dPdnB

    if (index[0]>0) {
      std::vector<size_t> im1={index[0]-1,index[1],index[2]};
      std::vector<size_t> ip1={index[0]+1,index[1],index[2]};
      std::vector<size_t> jm1={index[0],index[1]-1,index[2]};
      std::vector<size_t> jp1={index[0],index[1]+1,index[2]};
      std::vector<size_t> km1={index[0],index[1],index[2]-1};
      std::vector<size_t> kp1={index[0],index[1],index[2]+1};
      //cout << "cs22: " << tgp_cs2.get(im1) << endl;
      
      /*
        cout << "6: " << dF_dnB << " "
        << (tgp_F->get(index)-tgp_F->get(im1))/hc_mev_fm/
        (nB_grid[index[0]]-nB_grid[index[0]-1]) << " "
        << (tgp_F->get(ip1)-tgp_F->get(index))/hc_mev_fm/
        (nB_grid[index[0]+1]-nB_grid[index[0]]) << " ";
        cout << "7: " << dF_dYe << " "
        << (tgp_F->get(index)-tgp_F->get(jm1))/hc_mev_fm/
        (Ye_grid[index[1]]-Ye_grid[index[1]-1]) << " "
        << (tgp_F->get(jp1)-tgp_F->get(index))/hc_mev_fm/
        (Ye_grid[index[1]+1]-Ye_grid[index[1]]) << endl;
        cout << "8: " << dF_dT << " "
        << (tgp_F->get(index)-tgp_F->get(km1))/
        (T_grid[index[2]]-T_grid[index[2]-1]) << " "
        << (tgp_F->get(kp1)-tgp_F->get(index))/
        (T_grid[index[2]+1]-T_grid[index[2]]) << endl;
      */
      
      double t1=(tgp_F->get(index)-tgp_F->get(im1))/hc_mev_fm/
        (nB_grid[index[0]]-nB_grid[index[0]-1]);
      double t2=(tgp_F->get(ip1)-tgp_F->get(index))/hc_mev_fm/
        (nB_grid[index[0]+1]-nB_grid[index[0]]);
      double t3=(t2-t1)*2.0/(nB_grid[index[0]+1]-nB_grid[index[0]-1]);
      cout << F_nBnB << " " << t3 << " ";
     if (!(index[1]==0)) { 
        t1=(tgp_F->get(index)-tgp_F->get(jm1))/hc_mev_fm/
          (Ye_grid[index[1]]-Ye_grid[index[1]-1]);
        t2=(tgp_F->get(jp1)-tgp_F->get(index))/hc_mev_fm/
          (Ye_grid[index[1]+1]-Ye_grid[index[1]]);
        t3=(t2-t1)*2.0/(Ye_grid[index[1]+1]-Ye_grid[index[1]-1]);
        cout << F_YeYe << " " << t3 << " ";
     }
      
      if (false) {
        t1=(tgp_F->get(index)-tgp_F->get(km1))/
          (T_grid[index[2]]-T_grid[index[2]-1]);
        t2=(tgp_F->get(kp1)-tgp_F->get(index))/
          (T_grid[index[2]+1]-T_grid[index[2]]);
        t3=hc_mev_fm*(t1-t2)*2.0/(T_grid[index[2]+1]-T_grid[index[2]-1]);
        cout << "c: " << F_TT << " " << t3 << endl;
      }
      
      //double mun=Fintp/hc_mev_fm-Ye*dF_dYe+nB*dF_dnB;
      //double mup=Fintp/hc_mev_fm+(1.0-Ye)*dF_dYe+nB*dF_dnB-mue;
     // double dmun_dnB=2*dF_dnB-Ye*F_nBYe+nB*F_nBnB;
      //cout << "9: " << 2*dF_dnB << " " << -Ye*F_nBYe << " "
      //<< nB*F_nBnB << endl;
        
      //double dmup_dnB=2.0*dF_dnB+(Ye-1.0)*F_nBYe+nB*F_nBnB;

      double dmundnB=(f_nnnn*(1-Ye))+(f_nnnp*Ye);
      double dmupmuednB=(f_nnnp*(1-Ye))+(f_npnp*Ye);
      //double dmundnB=F_nBnB-(Ye*(((1/nB)*F_nBYe)-((1/(nB*nB))*dF_dYe)));
      //double dmupmuednB=F_nBnB-((1-Ye)*(((1/nB)*F_nBYe)-((1/(nB*nB))*dF_dYe)));
      double dPdnB=(dmundnB*nB*(1-Ye))+(mun*(1-Ye))+(dmupmuednB*Ye*nB)+(Ye*(mup+mue))-((mun*(1-Ye))+((mup+mue)*Ye));
      if (compare) {
        //std::cout << "dmun_dnB dmup_dnB" << endl;
        std::cout << "dmun_dnB d(mup+mue)_dnB" << endl;
        cout << "Computed:   ";
        //std::cout << dmun_dnB << " " << dmup_dnB << std::endl;
        std::cout << dmundnB << " " << dmupmuednB << std::endl;
        cout << "From table: ";
        std::cout << (tgp_mun->get(index)-tgp_mun->get(im1))/hc_mev_fm/
          (nB_grid[index[0]]-nB_grid[index[0]-1]) << " ";
        //std::cout << (tgp_mup->get(index)-tgp_mup->get(im1))/hc_mev_fm/
         // (nB_grid[index[0]]-nB_grid[index[0]-1]) << std::endl;

        std::cout << ((tgp_mup->get(index)+tgp_mue->get(index))-(tgp_mup->get(im1)+tgp_mue->get(im1)))/hc_mev_fm/
          (nB_grid[index[0]]-nB_grid[index[0]-1]) << std::endl;
        cout << "From table: ";
        std::cout << (tgp_mun->get(ip1)-tgp_mun->get(index))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]) << " ";
        std::cout << (tgp_mup->get(ip1)-tgp_mup->get(index))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]) << std::endl;
      }
      
      //double dmuden_dnB=((1.0-Ye)*mun+nn2*dmun_dnB+
      //                   Ye*(mup+mue)+np2*dmup_dnB);
      //double dPdnB=dmuden_dnB-dF_dnB*nB-Fintp/hc_mev_fm;
      if (compare) {
        cout << dPdnB << " ";
        /*
          cout << "dP_dnB: " << endl;
          std::cout << "Computed: " << dPdnB << std::endl;
          std::cout << "From table: "
          << (tgp_P->get(index)-tgp_P->get(im1))/hc_mev_fm/
          (nB_grid[index[0]]-nB_grid[index[0]-1]) << " "
          << (tgp_P->get(ip1)-tgp_P->get(index))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]) << std::endl;
        */
        cout << (tgp_P->get(index)-tgp_P->get(im1))/hc_mev_fm/
          (nB_grid[index[0]]-nB_grid[index[0]-1]) << " "
             << (tgp_P->get(ip1)-tgp_P->get(index))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]) << std::endl;
        
      }
      cout.unsetf(ios::showpos);
      if (dPdnB<=0.0) {
        cout << "Return failure for dPdnB<=0.0." << endl;
        return 2;
      }
    }
    if (cs_sq>1.0 || cs_sq<0.0 || !std::isfinite(cs_sq)) {
      cout << "Return failure for unphysical cs_sq." << endl;
      return 1;
    }
      
  }
        
return 0;
}
