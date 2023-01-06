/*
  -------------------------------------------------------------------
  
  Copyright (C) 2022, Andrew W. Steiner
  
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

          if (false) {
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
  mcovar_funct_rbf mfr;
  mfr.len.resize(3);
  
  //vector<double> len_list={0.1,0.3,0.5,0.75,1.0,2.0,3.0,4.0,5.0,6.0,
  //8.0,10.0};
  vector<double> len_list={2.0,3.0};
  vector<vector<double>> plist;
  plist.push_back(len_list);
  plist.push_back(len_list);
  plist.push_back(len_list);
  cout << "Going to set_covar()." << endl;
  ike.set_covar(mfr,plist);

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
    
    // This is the derivative with respect to the index
    // df/di*didnB

    // Derivatives of the physical coordinates with respect to the indices
    double dnBdi=2.0*0.04*log(10.0)*pow(10.0,((double)inB)*0.04-12.0);
    double dYedj=0.01;
    double dTdk=0.1*log(1.046)*pow(1.046,iT);

    ike.eval(point,out);
    double Fintp=out[0];
    
    ike.deriv(point,out,0);
    double F_nB=out[0]/hc_mev_fm/dnBdi;
    ike.deriv(point,out,1);
    double F_Ye=out[0]/hc_mev_fm/dYedj;
    ike.deriv(point,out,2);
    double F_T=out[0]/hc_mev_fm/dTdk;

    {
      vector<double> pointx={((double)inB+1),((double)iYe),((double)iT)};
      double F_nB2=(tg_F.get(pointx)-tg_F.get(point))/hc_mev_fm/
        (nB_grid2[inB+1]-nB_grid2[inB]);
      cout << "F_nB: " << F_nB << " " << F_nB2 << " "
           << F_nB2/F_nB << endl;
    }
    {
      vector<double> pointx={((double)inB),((double)iYe+1),((double)iT)};
      double F_Ye2=(tg_F.get(pointx)-tg_F.get(point))/hc_mev_fm/
        (Ye_grid2[iYe+1]-Ye_grid2[iYe]);
      cout << "F_Ye: " << F_Ye << " " << F_Ye2 << " "
           << F_Ye2/F_Ye << endl;
    }
    {
      vector<double> pointx={((double)iT),((double)iYe),((double)iT+1)};
      double F_T2=(tg_F.get(pointx)-tg_F.get(point))/
        (T_grid2[iT+1]-T_grid2[iT])/pow(hc_mev_fm,2.0);
      cout << "F_T: " << F_T << " " << F_T2 << endl;
    }
    
    double mun=Fintp/hc_mev_fm-Ye*F_Ye+nB*F_nB;
    double mup=Fintp/hc_mev_fm+(1.0-Ye)*F_Ye+nB*F_nB;
    double mue=tg_mue.get(index)/hc_mev_fm;

    // The entropy density
    double en=-nB*F_T;

    vector<double> point2={nB_cent,Ye_cent,T_cent*hc_mev_fm};

    cout << "F[MeV]: " << tg_F.get(point) << " " << Fintp << endl;
    cout << "mun,mup,mue: ";
    cout << tg_mun.interp_linear(point2) << " ";
    cout << tg_mup.interp_linear(point2) << " ";
    cout << tg_mue.interp_linear(point2) << endl;
    cout << "A: " << Fintp << " " << (-Ye*F_Ye)*hc_mev_fm << endl;
    cout << "A: " << Fintp << " " << (nB*F_nB)*hc_mev_fm << endl;
    cout << "B: " << Fintp << " " << ((1.0-Ye)*F_Ye)*hc_mev_fm << endl;
    cout << "B: " << Fintp << " " << (nB*F_nB)*hc_mev_fm << endl;
    cout << "mun,mup,mue: " << mun*hc_mev_fm << " "
         << mup*hc_mev_fm << " " << mue*hc_mev_fm << endl;
    cout << "en: " << tg_S.interp_linear(point2)*nB << " ";
    cout << en*hc_mev_fm << endl;
    exit(-1);
    
    ike.deriv2(point,out,0,0);
    double F_nBnB=out[0]/hc_mev_fm;
    ike.deriv2(point,out,0,1);
    double F_nBYe=out[0]/hc_mev_fm;
    ike.deriv2(point,out,1,1);
    double F_YeYe=out[0]/hc_mev_fm;
    ike.deriv2(point,out,0,2);
    double F_nBT=out[0]/hc_mev_fm;
    ike.deriv2(point,out,1,2);
    double F_YeT=out[0]/hc_mev_fm;
    ike.deriv2(point,out,2,2);
    double F_TT=out[0]/hc_mev_fm;
    
    double f_nnnn=(Ye*Ye*F_YeYe+nB*(2.0*F_nB-2.0*Ye*F_nBYe+nB*F_nBnB))/nB;
    double f_nnnp=((Ye-1.0)*Ye*F_YeYe+nB*(2.0*F_nB+(1.0-2.0*Ye)*
                                          F_nBYe+nB*F_nBnB))/nB;
    double f_npnp=((Ye-1.0)*(Ye-1.0)*F_YeYe+nB*(2.0*F_nB-2.0*(Ye-1.0)*
                                                F_nBYe+nB*F_nBnB))/nB;
    double f_nnT=F_T-Ye*F_YeT+nB*F_nBT;
    double f_npT=F_T-(Ye-1.0)*F_YeT+nB*F_nBT;
    double f_TT=nB*F_TT;
    
    double den=en*T_MeV/hc_mev_fm+mun*nB*(1.0-Ye)+mup*nB*Ye+mue*nB*Ye;
    double nn2=nB*(1.0-Ye);
    double np2=nB*Ye;
    double cs_sq=(nn2*nn2*(f_nnnn-f_nnT*f_nnT/f_TT)+
                  2.0*nn2*np2*(f_nnnp-f_nnT*f_npT/f_TT)+
                  np2*np2*(f_npnp-f_npT*f_npT/f_TT)-
                  2.0*en*(nn2*f_nnT/f_TT+np2*f_npT/f_TT)-en*en/f_TT)/den;

    exit(-1);
  }
  
  return 0;
}

