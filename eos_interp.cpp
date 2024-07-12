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

#ifdef NEVER_DEFINED
int eos_nuclei::interp_fix_table(std::vector<std::string> &sv,
                                 bool itive_com) {

  std::string st_in, st_out;
  std::string table_out;
  st_in=sv[1];
  size_t window=o2scl::stoszt(sv[2]);
  table_out=sv[3];
  st_out=sv[4];
  
  // Create interpolation object
  interpm_krige_eos ike;
  ike.mode=ike.mode_loo_cv_bf;
  ike.full_min=true;
  ike.def_mmin.verbose=1;
  ike.enp=this;
  
  /// Load cs2 from a file
  cout << "eos_nuclei::interp_fix_table() reading stability file: "
       << st_in << endl;
  hdf_file hff;
  hff.open(st_in);
  hdf_input(hff,tg_cs2,"cs2");
  hff.close();

  int ipx_count=0;

  cout << "Starting loop over entire grid:." << endl;
  for(size_t i=0;i<nB_grid2.size();i++) {
    for(size_t j=0;j<Ye_grid2.size();j++) {
      for(size_t k=0;k<T_grid2.size();k++) {
        
        vector<size_t> ix={i,j,k};

        double dPdnB;
        if (i>0 && i<nB_grid2.size()-1) {
          vector<size_t> ixp1={i+1,j,k};
          vector<size_t> ixm1={i-1,j,k};
          dPdnB=(tg_P.get(ixp1)-tg_P.get(ixm1))/
            (nB_grid2[i+1]-nB_grid2[i-1])/2;
        } else if (i>0) {
          vector<size_t> ixm1={i-1,j,k};
          dPdnB=(tg_P.get(ix)-tg_P.get(ixm1))/
            (nB_grid2[i]-nB_grid2[i-1]);
        } else {
          vector<size_t> ixp1={i+1,j,k};
          dPdnB=(tg_P.get(ixp1)-tg_P.get(ix))/
            (nB_grid2[i+1]-nB_grid2[i]);
        }
        
        if (tg_cs2.get(ix)>1.0 ||
            tg_cs2.get(ix)<0.0 || 
            !std::isfinite(tg_cs2.get(ix)) || 
            dPdnB<=0.0 || !std::isfinite(dPdnB)) {
          
          size_t i_fix=i, j_fix=j, k_fix=k;
          cout << "Found point to fix at (" << i << "," << j << ","
               << k << ") = (" << nB_grid2[i] << ","
               << Ye_grid2[j] << "," << T_grid2[k] << ")\n  cs2: "
               << tg_cs2.get(ix) << " dPdnB: " << dPdnB << endl;

          if (nB_grid2[i]<1.0e-3) {
          
            bool done=false;
            double delta_F=tg_Fint.get(ix)/1.0e5;
            double Fint_0=tg_Fint.get(ix);
            
            for(size_t ell=1;ell<20 && done==false;ell++) {
              
              for(size_t em=0;em<2 && done==false;em++) {
                
                if (em==0) {
                  tg_Fint.get(ix)=Fint_0-delta_F*ell;
                } else {
                  tg_Fint.get(ix)=Fint_0+delta_F*ell;
                }
                
                std::vector<std::string> sv2;
                
                // Computing the derivatives is fast, so we just do
                // the full table
                eos_deriv(sv2,false);
                
                size_t i_min, i_max, j_min, j_max, k_min, k_max;
                if (i_fix>0) i_min=i_fix-1;
                else i_min=i_fix;
                if (i_fix>=n_nB2-1) i_max=n_nB2-1;
                else i_max=i_fix+1;
                if (j_fix>0) j_min=j_fix-1;
                else j_min=j_fix;
                if (j_fix>=n_Ye2-1) j_max=n_Ye2-1;
                else j_max=j_fix+1;
                if (k_fix>0) k_min=k_fix-1;
                else k_min=k_fix;
                if (k_fix>=n_T2-1) k_max=n_T2-1;
                else k_max=k_fix+1;
                
                sv2={"stability",o2scl::szttos(i_min),o2scl::szttos(i_max),
                     o2scl::szttos(j_min),o2scl::szttos(j_max),
                     o2scl::szttos(k_min),o2scl::szttos(k_max)};
                
                // Update the electron-photon EOS for the specified points
                add_eg(sv2,false);
                
                // Update the stability and second derivatives
                stability(sv2,false);
                
                cout << "XX " << ell << " " << i_fix << " " << j_fix
                     << " " << k_fix << " " << tg_Fint.get(ix) << " "
                     << n_stability_fail << endl;
                
                if (n_stability_fail==0) {
                  done=true;
                  ipx_count++;
                }
              }
              
            }
            
          } else if (false) {
            
            // Create a copy of the free energy for temporary storage
            tg_Fint_old=tg_Fint;
            tg_F_old=tg_F;
            
            int ii_ret=interp_internal(i_fix,j_fix,k_fix,window,ike);
            cout << "Herexx." << endl;
            exit(-1);
            if (ii_ret!=0) {
              cerr << "Interpolation failed." << endl;
              O2SCL_ERR("Interpolation failed.",o2scl::exc_efailed);
            }
            
            std::vector<std::string> sv2;
            eos_deriv(sv2,itive_com);
            add_eg(sv2,itive_com);
            
            size_t i_min, i_max, j_min, j_max, k_min, k_max;
            i_min=ike.fix_list[0];
            i_max=ike.fix_list[0];
            j_min=ike.fix_list[1];
            j_max=ike.fix_list[1];
            k_min=ike.fix_list[2];
            k_max=ike.fix_list[2];
            
            for(size_t ifx=3;ifx<ike.fix_list.size();ifx+=3) {
              if (ike.fix_list[ifx]<i_min) {
                i_min=ike.fix_list[ifx];
              }
              if (ike.fix_list[ifx]>i_max) {
                i_max=ike.fix_list[ifx];
              }
              if (ike.fix_list[ifx+1]<j_min) {
                j_min=ike.fix_list[ifx+1];
              }
              if (ike.fix_list[ifx+1]>j_max) {
                j_max=ike.fix_list[ifx+1];
              }
              if (ike.fix_list[ifx+2]<k_min) {
                k_min=ike.fix_list[ifx+2];
              }
              if (ike.fix_list[ifx+2]>k_max) {
                k_max=ike.fix_list[ifx+2];
              }
            }
            for(size_t icx=3;icx<ike.calib_list.size();icx+=3) {
              if (ike.calib_list[icx]<i_min) {
                i_min=ike.calib_list[icx];
              }
              if (ike.calib_list[icx]>i_max) {
                i_max=ike.calib_list[icx];
              }
              if (ike.calib_list[icx+1]<j_min) {
                j_min=ike.calib_list[icx+1];
              }
              if (ike.calib_list[icx+1]>j_max) {
                j_max=ike.calib_list[icx+1];
              }
              if (ike.calib_list[icx+2]<k_min) {
                k_min=ike.calib_list[icx+2];
              }
              if (ike.calib_list[icx+2]>k_max) {
                k_max=ike.calib_list[icx+2];
              }
            }
            for(size_t ii=0;ii<window;ii++) {
              if (i_min>0) i_min--;
              if (i_max<nB_grid2.size()-1) i_max++;
              if (j_min>0) j_min--;
              if (j_max<Ye_grid2.size()-1) j_max++;
              if (k_min>0) k_min--;
              if (k_max<T_grid2.size()-1) k_max++;
            }
            
            cout << "Computed i_min, i_max: " << i_min << " " << i_max << endl;
            cout << "Computed j_min, j_max: " << j_min << " " << j_max << endl;
            cout << "Computed k_min, k_max: " << k_min << " " << k_max << endl;
            
            std::vector<std::string> sv3;
            sv3={"stability",o2scl::szttos(i_min),o2scl::szttos(i_max),
                 o2scl::szttos(j_min),o2scl::szttos(j_max),
                 o2scl::szttos(k_min),o2scl::szttos(k_max)};
            stability(sv3,itive_com);
            
            k+=window*2;
            
            ipx_count++;
          
          }
          
          if (false) {
            hdf_file hf;
            hf.open_or_create(st_out);
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
            hdf_output(hf,tg_cs2,"cs2");
            hdf_output(hf,tg_cs2_hom,"cs2_hom");
            hf.close();
            
            write_results(table_out);
          }            
          
          cout << "ipx_count: " << ipx_count << endl;
          
          //if (ipx_count==10) {
          //exit(-1);
          //}

          if (ipx_count==100) {
            i=nB_grid2.size();
            j=Ye_grid2.size();
            k=T_grid2.size();
          }
          
        }
      }
      
    }
  }
  
  if (true) {
    hdf_file hf;
    hf.open_or_create(st_out);
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
    hdf_output(hf,tg_cs2,"cs2");
    hdf_output(hf,tg_cs2_hom,"cs2_hom");
    hf.close();
    
    write_results(table_out);
  }            
          
  return 0;
}

int eos_nuclei::interp_internal(size_t i_fix, size_t j_fix, size_t k_fix,
                                size_t window, interpm_krige_eos &ike) {
  
  // Using the specified window, compute the list of points to fix,
  // and the list of calibration points

  if (tg_cs2.get_rank()<3) {
    O2SCL_ERR("No cs2 in interp_internal().",o2scl::exc_einval);
  }
  
  int iwindow=((int)window);
  for(int dnB=-iwindow;dnB<=iwindow;dnB++) {
    for(int dYe=-iwindow;dYe<=iwindow;dYe++) {
      for(int dT=-iwindow;dT<=iwindow;dT++) {
        vector<size_t> index={i_fix+dnB,j_fix+dYe,k_fix+dT};
        if (abs(dnB)+abs(dYe)+abs(dT)<=iwindow &&
            i_fix+dnB>=0 && i_fix+dnB<n_nB2 &&
            j_fix+dYe>=0 && j_fix+dYe<n_Ye2 &&
            k_fix+dT>=0 && k_fix+dT<n_T2) {
          double dPdnB;
          size_t i=((size_t)(((int)i_fix)+dnB));
          size_t j=((size_t)(((int)j_fix)+dYe));
          size_t k=((size_t)(((int)k_fix)+dT));
          if (i>0 && i<nB_grid2.size()-1) {
            vector<size_t> ixp1={i+1,j,k};
            vector<size_t> ixm1={i-1,j,k};
            dPdnB=(tg_P.get(ixp1)-tg_P.get(ixm1))/
              (nB_grid2[i+1]-nB_grid2[i-1])/2;
          } else if (i>0) {
            vector<size_t> ixm1={i-1,j,k};
            dPdnB=(tg_P.get(index)-tg_P.get(ixm1))/
              (nB_grid2[i]-nB_grid2[i-1]);
          } else {
            vector<size_t> ixp1={i+1,j,k};
            dPdnB=(tg_P.get(ixp1)-tg_P.get(index))/
              (nB_grid2[i+1]-nB_grid2[i]);
          }
          
          if (dPdnB<=0.0 ||
              !std::isfinite(dPdnB) ||
              tg_cs2.get(index)>1.0 ||
              !std::isfinite(tg_cs2.get(index)) ||
              tg_cs2.get(index)<0.0) {
            ike.fix_list.push_back(index[0]);
            ike.fix_list.push_back(index[1]);
            ike.fix_list.push_back(index[2]);
          } else {
            ike.calib_list.push_back(index[0]);
            ike.calib_list.push_back(index[1]);
            ike.calib_list.push_back(index[2]);
          }
        }
      }
    }
  }
  
  cout << "In eos_nuclei:interp_internal(), using "
       << ike.calib_list.size()/3 << " points to calibrate"
       << endl;
  cout << "  and attempting to fix " << ike.fix_list.size()/3 << " points."
       << endl;
  size_t count=ike.calib_list.size()/3;

  if (ike.fix_list.size()==0) {
    cerr << "No points to fix." << endl;
    return 1;
  }

  // Compute the distanes between the calibration points and the
  // points to fix
  
  ike.compute_dists();
  if (ike.calib_list.size()/3!=ike.calib_dists.size()) {
    O2SCL_ERR("Error in calibration distance math.",o2scl::exc_esanity);
  }

  // Initialize the covariance object 
  
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

  // Use the covariance object to set the data, skipping the
  // optimization for now and performing it manually
  
  ike.skip_optim=true;
  ike.set();

  ike.addl_verbose=1;

  // Manually optimize the Gaussian process interpolation by
  // exhaustively searching
  
  double min_qual=1.0e99;
  vector<double> p(4), min_p;
  
  //for(p[0]=80.0;p[0]>8.0;p[0]/=1.4) {
  //for(p[1]=80.0;p[1]>8.0;p[1]/=1.4) {
  //for(p[2]=80.0;p[2]>8.0;p[2]/=1.4) {
  //for(p[3]=-15.0;p[3]<-2.99;p[3]+=1.0) {
  std::cout << "In eos_nuclei::interp_internal(), loop over "
            << "hyperparameters." << std::endl;
  for(p[0]=10.0;p[0]>1.99;p[0]/=1.4) {
    for(p[1]=10.0;p[1]>1.99;p[1]/=1.4) {
      for(p[2]=10.0;p[2]>1.99;p[2]/=1.4) {
        for(p[3]=-15.0;p[3]<-14.99;p[3]+=1.0) {

          if (ike.addl_verbose>=1) {
            cout << "Covariance parameters: ";
            vector_out(cout,p,true);
          }
          
          (*ike.cf)[0].set_params(p);
          int success;
          double q=ike.qual_fun(0,success);
          if (ike.addl_verbose>=1) {
            cout << "q,min_qual,success: "
                 << q << " " << min_qual << " " << success << endl;
            cout << endl;
          }
	  exit(-1);
          
          if (success==0 && q<min_qual) {
            min_p=p;
            min_qual=q;
          }
        }
      }
    }
  }

  // Choose the best point
  
  if (min_qual>0.9e99) {
    cerr << "All points failed in Gaussian process interpolation." << endl;
    return 2;
  } else {
    (*ike.cf)[0].set_params(min_p);
    int st;
    cout << "Last run:\n" << endl;
    vector_out(cout,min_p,true);
    double qt=ike.qual_fun(0,st);
    cout << "min_qual,qt,st,min_p: "
         << min_qual << " " << qt << " " << st << " ";
    vector_out(cout,min_p,true);
    cout << "Success." << endl;
  }

  // Use the interpolation results to modify the points to be fixed
  
  std::vector<double> out(1);
  
  for(size_t j=0;j<ike.fix_list.size();j+=3) {

    vector<size_t> index={ike.fix_list[j],ike.fix_list[j+1],
      ike.fix_list[j+2]};
    
    double nB=nB_grid2[index[0]];
    double Ye=Ye_grid2[index[1]];
    double T_MeV=T_grid2[index[2]];

    // Electron photon contribution in MeV
    double F_eg=tg_F.get(index)-tg_Fint.get(index);
    
    ike.eval(index,out);
    if (ike.interp_Fint==false) {
      cout << "Change (5) from " << tg_F.get(index) << " to "
           << out[0] << endl;
      cout << "Change (5b) from " << tg_Fint.get(index) << " to "
           << out[0]-F_eg << endl;
      tg_F.get(index)=out[0];
      tg_Fint.get(index)=out[0]-F_eg;
    } else {
      cout << "Change (6) from " << tg_Fint.get(index) << " to "
           << out[0] << endl;
      tg_Fint.get(index)=out[0];
      tg_F.get(index)=out[0]+F_eg;
    }
    
  }

  // Use the interpolation results to modify the calibration points nearby

  for(size_t j=0;j<ike.calib_list.size();j+=3) {

    vector<size_t> index={ike.calib_list[j],ike.calib_list[j+1],
      ike.calib_list[j+2]};
    
    double nB=nB_grid2[index[0]];
    double Ye=Ye_grid2[index[1]];
    double T_MeV=T_grid2[index[2]];
    
    double fact=1.0/(1.0+exp(2.0*(ike.calib_dists[j/3]-window/2.0)));

    // Electron photon contribution in MeV
    double F_eg=tg_F.get(index)-tg_Fint.get(index);
    
    ike.eval(index,out);
    if (ike.interp_Fint==false) {
      double corr=out[0]-tg_F.get(index);
      cout << "Change (7) from " << tg_F.get(index) << " to "
           << tg_F.get(index)+fact*corr << " at dist: "
           << ike.calib_dists[j/3] << endl;
      cout << "Change (7b) from " << tg_Fint.get(index) << " to "
           << tg_F.get(index)-F_eg << " at dist: "
           << ike.calib_dists[j/3] << endl;
      tg_F.get(index)+=fact*corr;
      tg_Fint.get(index)=tg_F.get(index)-F_eg;
    } else {
      double corr=out[0]-tg_Fint.get(index);
      cout << "Change (8) from " << tg_Fint.get(index) << " to "
           << tg_Fint.get(index)+fact*corr << " at dist: "
           << ike.calib_dists[j/3] << endl;
      tg_Fint.get(index)+=fact*corr;
      tg_F.get(index)=out[0]+F_eg;
    }
    
  }
  
  // Make sure to set the derivative and lepton flags so that they
  // can be recomputed later
  
  derivs_computed=false;
  with_leptons=false;

  return 0;
}

int eos_nuclei::interp_point(std::vector<std::string> &sv,
                             bool itive_com) {

  if (sv.size()<6) {
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
  cout << "At grid point: " << inB << " " << iYe << " " << iT << endl;

  size_t window=o2scl::stoszt(sv[4]);
  cout << "Using window size: " << window << endl;

  // Create interpolation object
  interpm_krige_eos ike;
  ike.mode=ike.mode_loo_cv_bf;
  ike.full_min=true;
  ike.def_mmin.verbose=1;
  ike.enp=this;
  
  cout << "Making a copy." << endl;
  tg_Fint_old=tg_Fint;
  tg_F_old=tg_F;
  cout << "Done making a copy." << endl;
  
  /// Load cs2 from a file
  std::string st_o2="";
  st_o2=sv[5];
  hdf_file hff;
  hff.open(st_o2);
  hdf_input(hff,tg_cs2);
  hff.close();

  interp_internal(inB,iYe,iT,window,ike);

  return 0;
}

double interpm_krige_eos::dist_cf(size_t i_calib, size_t i_fix) {

  // These are size_t's so we have to convert to double before
  // we subtract
  double dist1=((double)calib_list[i_calib*3])-
    ((double)fix_list[i_fix*3]);
  double dist2=((double)calib_list[i_calib*3+1])-
    ((double)fix_list[i_fix*3+1]);
  double dist3=((double)calib_list[i_calib*3+2])-
    ((double)fix_list[i_fix*3+2]);
  
  return sqrt(dist1*dist1+dist2*dist2+dist3*dist3);
}

void interpm_krige_eos::compute_dists() {

  // Make sure we start with an empty array
  calib_dists.clear();
  
  // Collect counts
  size_t calib_count=calib_list.size()/3;
  size_t fix_count=fix_list.size()/3;

  // Compute a distance for each calibration points
  for(size_t i_calib=0;i_calib<calib_count;i_calib++) {

    // Compute the minimum distance over all points to fix
    double dist_min=dist_cf(i_calib,0);
    for(size_t i_fix=1;i_fix<fix_count;i_fix++) {
      double dist=dist_cf(i_calib,i_fix);
      if (dist<dist_min) dist_min=dist;
    }
    
    calib_dists.push_back(dist_min);
  }
  
  return;
}

void interpm_krige_eos::set() {

  // Set the grids and the pointers to the tensor_grid objects

  eos_nuclei *enp2=(eos_nuclei *)enp;
  nB_grid=enp2->nB_grid2;
  Ye_grid=enp2->Ye_grid2;
  T_grid=enp2->T_grid2;

  tgp_F=&enp2->tg_F;
  tgp_P=&enp2->tg_P;
  tgp_S=&enp2->tg_S;
  tgp_mun=&enp2->tg_mun;
  tgp_mup=&enp2->tg_mup;
  tgp_mue=&enp2->tg_mue;
  tgp_Fint=&enp2->tg_Fint;
  tgp_Sint=&enp2->tg_Sint;
  tgp_cs2=&enp2->tg_cs2;
  tgp_Fint_old=&enp2->tg_Fint_old;
  tgp_F_old=&enp2->tg_F_old;

  mneut=enp2->neutron.m;
  mprot=enp2->proton.m;
  
  // Reformat the data into the format which interpm_krige_eos
  // uses
  
  ubmatrix ix(calib_list.size()/3,3);
  ubmatrix iy(1,calib_list.size()/3);

  cout << "In interpm_krige_eos::set(): Fix list: " << endl;
  if (interp_Fint) {
    cout << "ix[0] ix[1] ix[2] Fint" << endl;
  } else {
    cout << "ix[0] ix[1] ix[2] F" << endl;
  }
  for(size_t j=0;j<fix_list.size();j+=3) {
    ix(j/3,0)=fix_list[j];
    ix(j/3,1)=fix_list[j+1];
    ix(j/3,2)=fix_list[j+2];
    vector<size_t> index={fix_list[j],fix_list[j+1],
      fix_list[j+2]};
    if (interp_Fint) {
      iy(0,j/3)=tgp_Fint->get(index);
    } else {
      iy(0,j/3)=tgp_F->get(index);
    }
    if (true) {
      cout << ((int)ix(j/3,0)) << " " << ((int)ix(j/3,1)) << " "
           << ((int)ix(j/3,2)) << " "
           << nB_grid[fix_list[j]] << " "
           << Ye_grid[fix_list[j+1]] << " "
           << T_grid[fix_list[j+2]] << " "
           << iy(0,j/3) << " ";
      cout << endl;
    }
  }
  cout << "In interpm_krige_eos::set(): Calibration list: " << endl;
  if (interp_Fint) {
    cout << "ix[0] ix[1] ix[2] Fint dist" << endl;
  } else {
    cout << "ix[0] ix[1] ix[2] F dist" << endl;
  }
  for(size_t j=0;j<calib_list.size();j+=3) {
    ix(j/3,0)=calib_list[j];
    ix(j/3,1)=calib_list[j+1];
    ix(j/3,2)=calib_list[j+2];
    vector<size_t> index={calib_list[j],calib_list[j+1],
      calib_list[j+2]};
    if (interp_Fint) {
      iy(0,j/3)=tgp_Fint->get(index);
    } else {
      iy(0,j/3)=tgp_F->get(index);
    }
    if (true) {
      cout << ((int)ix(j/3,0)) << " " << ((int)ix(j/3,1)) << " "
           << ((int)ix(j/3,2)) << " "
           << nB_grid[calib_list[j]] << " "
           << Ye_grid[calib_list[j+1]] << " "
           << T_grid[calib_list[j+2]] << " "
           << iy(0,j/3) << " ";
      if (fix_list.size()>3) {
        cout << calib_dists[j/3];
      }
      cout << endl;
    }
  }
  
  // Make a copy because interpm_krige_eos will keep it
  ubmatrix ix2=ix;

  size_t n_nB=nB_grid.size();
  size_t n_Ye=Ye_grid.size();
  size_t n_T=T_grid.size();

  // Call the Gaussian process set_data() function
  
  std::cout << "Going to interp_krige_eos::set_data()." << std::endl;
  verbose=2;
  set_data(3,1,calib_list.size()/3,ix,iy);

  return;
}
  
int interpm_krige_eos::addl_const(size_t iout, double &ret) {

  std::cout << "In interpm_krige_eos::addl_const()." << std::endl;
  
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
  double window=5.0;

  if (true) {
    
    // Use the interpolation results to modify the points to be fixed
    
    std::vector<double> out(1);
    
    for(size_t j=0;j<fix_list.size();j+=3) {
      
      vector<size_t> index={fix_list[j],fix_list[j+1],
        fix_list[j+2]};
      
      double nB=nB_grid[index[0]];
      double Ye=Ye_grid[index[1]];
      double T_MeV=T_grid[index[2]];
      
      // Electron photon contribution in MeV
      double F_eg=tgp_F->get(index)-tgp_Fint->get(index);
      
      eval(index,out);

      if (true) {
        if (interp_Fint==false) {
          if (false) {
            cout << "Change F (fix) from " << tgp_F_old->get(index) << " "
                 << tgp_F->get(index) << " to "
                 << out[0] << endl;
            cout << "Change Fint (fix) from "
                 << tgp_Fint_old->get(index) << " "
                 << tgp_Fint->get(index) << " to "
                 << out[0]-F_eg << endl;
          } else {
            cout << index[0] << " " << index[1] << " "
                 << index[2] << " " << nB <<  " " << Ye << " "
                 << T_MeV << " " << tgp_F_old->get(index)
                 << " " << out[0] << endl;
          }
          tgp_F->get(index)=out[0];
          tgp_Fint->get(index)=out[0]-F_eg;
        } else {
          cout << "Change Fint (fix) from " << tgp_Fint->get(index)
               << " to " << out[0] << endl;
          tgp_Fint->get(index)=out[0];
          tgp_F->get(index)=out[0]+F_eg;
        }
      }
      
      
    }

    // Use the interpolation results to modify the calibration points nearby
    
    for(size_t j=0;j<calib_list.size();j+=3) {
      
      vector<size_t> index={calib_list[j],calib_list[j+1],
        calib_list[j+2]};
      
      double nB=nB_grid[index[0]];
      double Ye=Ye_grid[index[1]];
      double T_MeV=T_grid[index[2]];

      double fact=1.0/(1.0+exp(2.0*(calib_dists[j/3]-window/2.0)));
      
      // Lepton photon contribution in MeV
      double F_eg=tgp_F->get(index)-tgp_Fint->get(index);
      
      eval(index,out);
      //cout << "dist,fact,F,out: " << calib_dists[j/3] << " " << fact << " "
      //<< tgp_F->get(index) << " " << out[0] << endl;

      if (true) {
        if (interp_Fint==false) {
          double corr=out[0]-tgp_F->get(index);
          if (false) {
            cout << "Change F (calib) from " << tgp_F->get(index) << " to "
                 << tgp_F->get(index)+fact*corr << " at dist: "
                 << calib_dists[j/3] << endl;
            cout << "Change Fint (calib) from "
                 << tgp_Fint->get(index) << " to "
                 << tgp_F->get(index)-F_eg << " at dist: "
                 << calib_dists[j/3] << endl;
          } else {
            cout << index[0] << " " << index[1] << " "
                 << index[2] << " " << nB <<  " " << Ye << " "
                 << T_MeV << " " << tgp_F_old->get(index) << " "
                 << tgp_F_old->get(index)+fact*corr << " "
                 << fact*corr << endl;
          }
          tgp_F->get(index)+=fact*corr;
          tgp_Fint->get(index)=tgp_F->get(index)-F_eg;
        } else {
          double corr=out[0]-tgp_Fint->get(index);
          cout << "Change Fint (calib) from " << tgp_Fint->get(index)
               << " to "
               << tgp_Fint->get(index)+fact*corr << " at dist: "
               << calib_dists[j/3] << endl;
          tgp_Fint->get(index)+=fact*corr;
          tgp_F->get(index)=out[0]+F_eg;
        }
      }
      
    }

    size_t i_min, i_max, j_min, j_max, k_min, k_max;
    i_min=fix_list[0];
    i_max=fix_list[0];
    j_min=fix_list[1];
    j_max=fix_list[1];
    k_min=fix_list[2];
    k_max=fix_list[2];
    
    for(size_t ifx=3;ifx<fix_list.size();ifx+=3) {
      if (fix_list[ifx]<i_min) {
        i_min=fix_list[ifx];
      }
      if (fix_list[ifx]>i_max) {
        i_max=fix_list[ifx];
      }
      if (fix_list[ifx+1]<j_min) {
        j_min=fix_list[ifx+1];
      }
      if (fix_list[ifx+1]>j_max) {
        j_max=fix_list[ifx+1];
      }
      if (fix_list[ifx+2]<k_min) {
        k_min=fix_list[ifx+2];
      }
      if (fix_list[ifx+2]>k_max) {
        k_max=fix_list[ifx+2];
      }
    }
    for(size_t icx=0;icx<calib_list.size();icx+=3) {
      if (calib_list[icx]<i_min) {
        i_min=calib_list[icx];
      }
      if (calib_list[icx]>i_max) {
        i_max=calib_list[icx];
      }
      if (calib_list[icx+1]<j_min) {
        j_min=calib_list[icx+1];
      }
      if (calib_list[icx+1]>j_max) {
        j_max=calib_list[icx+1];
      }
      if (calib_list[icx+2]<k_min) {
        k_min=calib_list[icx+2];
      }
      if (calib_list[icx+2]>k_max) {
        k_max=calib_list[icx+2];
      }
    }

    // Expand by a factor of three if we're not at the boundary
    //for(size_t ii=0;ii<window;ii++) {

    // Expand by one unit to make sure to catch any interface problems
    if (true) {
      if (i_min>0) i_min--;
      if (i_max<nB_grid.size()-1) i_max++;
      if (j_min>0) j_min--;
      if (j_max<Ye_grid.size()-1) j_max++;
      if (k_min>0) k_min--;
      if (k_max<T_grid.size()-1) k_max++;
    }
    
    cout << "Computed i_min, i_max: " << i_min << " " << i_max << endl;
    cout << "Computed j_min, j_max: " << j_min << " " << j_max << endl;
    cout << "Computed k_min, k_max: " << k_min << " " << k_max << endl;

    if (false) {
      tensor_grid3 tg3x;

      // Create a mapping between index and baryon density
      vector<double> gi, gnb;
      for(size_t i=i_min;i<=i_max;i++) {
        gi.push_back(i);
        gnb.push_back(nB_grid[i]);
      }
      interp_vec<vector<double>> iv(gi.size(),gi,gnb,itp_steffen);
      
      vector<vector<double>> g(3);
      for(double id=i_min;id<i_max+1.0e-6;id+=0.1) {
        g[0].push_back(iv.eval(id));
      }
      for(size_t j=j_min;j<=j_max;j++) {
        g[1].push_back(Ye_grid[j]);
      }
      for(size_t k=k_min;k<=k_max;k++) {
        g[2].push_back(T_grid[k]);
      }
      vector<size_t> sz={g[0].size(),g[1].size(),g[2].size()};
      tg3x.resize(3,sz);
      
      tg3x.set_grid(g);
      for(size_t j=j_min;j<=j_max;j++) {
        for(size_t k=k_min;k<=k_max;k++) {
          size_t i=0;
          for(double id=i_min;id<i_max+1.0e-6;id+=0.1) {
            vector<size_t> ix={i,j,k};
            vector<double> iix={id,((double)j),((double)k)};
            //cout << i << " " << id << " " << j << " " << k << endl;
            eval(iix,out);
            tg3x.set(i,j-j_min,k-k_min,out[0]);
            i++;
          }
        }
      }

      table_units<> tux;
      tux.line_of_names("i j k nB Ye T F out");
      for(size_t i=i_min;i<=i_max;i++) {
        for(size_t j=j_min;j<=j_max;j++) {
          for(size_t k=k_min;k<=k_max;k++) {
            vector<size_t> ix={i,j,k};
            vector<double> iix={((double)i),((double)j),((double)k)};
            eval(iix,out);
            vector<double> line={((double)i),((double)j),((double)k),
                                 nB_grid[i],Ye_grid[j],T_grid[k],
                                 tgp_F_old->get(ix),out[0]};
            tux.line_of_data(line.size(),line);
          }
        }
      }
      
      hdf_file hf;
      hf.open_or_create("itest.o2");
      hdf_output(hf,tg3x,"tg");
      hdf_output(hf,tux,"tu");
      hf.close();

      cout << "Done." << endl;
      exit(-1);
    }
    
    std::vector<std::string> sv2;

    eos_nuclei *enp2=(eos_nuclei *)enp;
  
    // Computing the derivatives is fast, so we just do the full table
    enp2->eos_deriv(sv2,false);

    sv2={"stability",o2scl::szttos(i_min),o2scl::szttos(i_max),
      o2scl::szttos(j_min),o2scl::szttos(j_max),
      o2scl::szttos(k_min),o2scl::szttos(k_max)};
    
    // Update the electron-photon EOS for the specified points
    enp2->add_eg(sv2,false);

    // Update the stability and second derivatives
    enp2->stability(sv2,false);
    
    // Change free energies back to original

    for(size_t j=0;j<fix_list.size();j+=3) {
      
      vector<size_t> index={fix_list[j],fix_list[j+1],
        fix_list[j+2]};

      tgp_Fint->get(index)=tgp_Fint_old->get(index);
      tgp_F->get(index)=tgp_F_old->get(index);
      
    }    
    
    for(size_t j=0;j<calib_list.size();j+=3) {
      
      vector<size_t> index={calib_list[j],calib_list[j+1],
        calib_list[j+2]};

      tgp_Fint->get(index)=tgp_Fint_old->get(index);
      tgp_F->get(index)=tgp_F_old->get(index);
      
    }    
    
    for(size_t j=0;j<fix_list.size();j+=3) {
      
      vector<size_t> ix={fix_list[j],fix_list[j+1],
        fix_list[j+2]};
      
      double dPdnB;
      if (ix[0]>0 && ix[0]<nB_grid.size()-1) {
        vector<size_t> ixp1={ix[0]+1,ix[1],ix[2]};
        vector<size_t> ixm1={ix[0]-1,ix[1],ix[2]};
        dPdnB=(tgp_P->get(ixp1)-tgp_P->get(ixm1))/
          (nB_grid[ix[0]+1]-nB_grid[ix[0]-1])/2;
      } else if (ix[0]>0) {
        vector<size_t> ixm1={ix[0]-1,ix[1],ix[2]};
        dPdnB=(tgp_P->get(ix)-tgp_P->get(ixm1))/
          (nB_grid[ix[0]]-nB_grid[ix[0]-1]);
      } else {
        vector<size_t> ixp1={ix[0]+1,ix[1],ix[2]};
        dPdnB=(tgp_P->get(ixp1)-tgp_P->get(ix))/
          (nB_grid[ix[0]+1]-nB_grid[ix[0]]);
      }
      
      if (tgp_cs2->get(ix)>1.0 ||
          tgp_cs2->get(ix)<0.0 || 
          !std::isfinite(tgp_cs2->get(ix)) || 
          dPdnB<0.0 || !std::isfinite(dPdnB)) {
        if (addl_verbose>=1) {
          cout << "Return failure at fix point." << endl;
        }
        if (addl_verbose>=2) {
          cout << endl;
        }
        return 1;
      }
      
    }    
    
    size_t n_calib_fail=0;
    
    for(size_t j=0;j<calib_list.size();j+=3) {
      
      vector<size_t> ix={calib_list[j],calib_list[j+1],
        calib_list[j+2]};

      double dPdnB;
      if (ix[0]>0 && ix[0]<nB_grid.size()-1) {
        vector<size_t> ixp1={ix[0]+1,ix[1],ix[2]};
        vector<size_t> ixm1={ix[0]-1,ix[1],ix[2]};
        dPdnB=(tgp_P->get(ixp1)-tgp_P->get(ixm1))/
          (nB_grid[ix[0]+1]-nB_grid[ix[0]-1])/2;
      } else if (ix[0]>0) {
        vector<size_t> ixm1={ix[0]-1,ix[1],ix[2]};
        dPdnB=(tgp_P->get(ix)-tgp_P->get(ixm1))/
          (nB_grid[ix[0]]-nB_grid[ix[0]-1]);
      } else {
        vector<size_t> ixp1={ix[0]+1,ix[1],ix[2]};
        dPdnB=(tgp_P->get(ixp1)-tgp_P->get(ix))/
          (nB_grid[ix[0]+1]-nB_grid[ix[0]]);
      }
      
      if (tgp_cs2->get(ix)>1.0 ||
          tgp_cs2->get(ix)<0.0 || 
          !std::isfinite(tgp_cs2->get(ix)) || 
          dPdnB<0.0 || !std::isfinite(dPdnB)) {
        n_calib_fail++;
      }
      
    }
    
    if (n_calib_fail>fix_list.size()/9) {
      if (addl_verbose>=1) {
        cout << n_calib_fail << " calibration failures and "
             << fix_list.size()/3 << " points to fix." << endl;
        cout << "Return failure." << endl;
      }
      if (addl_verbose>=2) {
        cout << endl;
      }
      return 2;
    }

    /*
    if (enp2->n_stability_fail>0) {
      if (addl_verbose>=1) {
        cout << "Return failure." << endl;
      }
      if (addl_verbose>=2) {
        cout << endl;
      }
      return 1;
    }
    */

    cout << "Success." << endl;
    
    exit(-1);
  }
  
  for(size_t ilist=0;ilist<(calib_list.size()+fix_list.size())/3;
      ilist++) {
    
    if (addl_verbose>=1) {
      std::cout << "i,nB,Ye,T: " << ilist << " ";
    }
    
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
      
    std::vector<size_t> index={((size_t)inB),((size_t)iYe),((size_t)iT)};
        
    double nB=nB_grid[inB];
    double Ye=Ye_grid[iYe];
    double T_MeV=T_grid[iT];
    
    if (addl_verbose>=1) {
      cout << inB << " " << iYe <<  " " << iT << " "
           << nB << " " << Ye << " " << T_MeV << " "
           << interp_Fint << endl;
    }
    
    // Derivatives of the physical coordinates with respect to the indices

    double didnB, djdYe, dkdT, dnBdi, dYedj, dTdk;
    double d2idnB2, d2jdYe2, d2kdT2;
    if (inB>0 && inB<nB_grid.size()-1) {
      dnBdi=(nB_grid[inB+1]-nB_grid[inB-1])/2;
      double t2=(nB_grid[inB+1]-nB_grid[inB]);
      double t1=(nB_grid[inB]-nB_grid[inB-1]);
      d2idnB2=(1.0/t2-1.0/t1)/(nB_grid[inB+1]-nB_grid[inB-1])*2.0;
    } else if (inB==0) {
      dnBdi=(nB_grid[1]-nB_grid[0]);
      double t2=(nB_grid[2]-nB_grid[1]);
      double t1=(nB_grid[1]-nB_grid[0]);
      d2idnB2=(1.0/t2-1.0/t1)/(nB_grid[2]-nB_grid[0])*2.0;
    } else {
      dnBdi=(nB_grid[nB_grid.size()-1]-nB_grid[nB_grid.size()-2]);
      double t2=(nB_grid[nB_grid.size()-1]-nB_grid[nB_grid.size()-2]);
      double t1=(nB_grid[nB_grid.size()-2]-nB_grid[nB_grid.size()-3]);
      d2idnB2=(1.0/t2-1.0/t1)/(nB_grid[nB_grid.size()-1]-
                               nB_grid[nB_grid.size()-3])*2.0;
    }
    if (iYe>0 && iYe<Ye_grid.size()-1) {
      dYedj=(Ye_grid[iYe+1]-Ye_grid[iYe-1])/2;
      double t2=(Ye_grid[iYe+1]-Ye_grid[iYe]);
      double t1=(Ye_grid[iYe]-Ye_grid[iYe-1]);
      d2jdYe2=(1.0/t2-1.0/t1)/(Ye_grid[iYe+1]-Ye_grid[iYe-1])*2.0;
    } else if (iYe==0) {
      dYedj=(Ye_grid[1]-Ye_grid[0]);
      double t2=(Ye_grid[2]-Ye_grid[1]);
      double t1=(Ye_grid[1]-Ye_grid[0]);
      d2jdYe2=(1.0/t2-1.0/t1)/(Ye_grid[2]-Ye_grid[0])*2.0;
    } else {
      dYedj=(Ye_grid[Ye_grid.size()-1]-Ye_grid[Ye_grid.size()-1]);
      double t2=(Ye_grid[Ye_grid.size()-1]-Ye_grid[Ye_grid.size()-2]);
      double t1=(Ye_grid[Ye_grid.size()-2]-Ye_grid[Ye_grid.size()-3]);
      d2jdYe2=(1.0/t2-1.0/t1)/(Ye_grid[Ye_grid.size()-1]-
                               Ye_grid[Ye_grid.size()-3])*2.0;
    }
    if (iT>0 && iT<T_grid.size()-1) {
      dTdk=(T_grid[iT+1]-T_grid[iT-1])/2;
      double t2=(T_grid[iT+1]-T_grid[iT]);
      double t1=(T_grid[iT]-T_grid[iT-1]);
      d2kdT2=(1.0/t2-1.0/t1)/(T_grid[iT+1]-T_grid[iT-1])*2.0;
    } else if (iT==0) {
      dTdk=(T_grid[1]-T_grid[0]);
      double t2=(T_grid[2]-T_grid[1]);
      double t1=(T_grid[1]-T_grid[0]);
      d2kdT2=(1.0/t2-1.0/t1)/(T_grid[2]-T_grid[0])*2.0;
    } else {
      dTdk=(T_grid[T_grid.size()-1]-T_grid[T_grid.size()-1]);
      double t2=(T_grid[T_grid.size()-1]-T_grid[T_grid.size()-2]);
      double t1=(T_grid[T_grid.size()-2]-T_grid[T_grid.size()-3]);
      d2kdT2=(1.0/t2-1.0/t1)/(T_grid[T_grid.size()-1]-
                              T_grid[T_grid.size()-3])*2.0;
    }
    didnB=1.0/dnBdi;
    djdYe=1.0/dYedj;
    dkdT=1.0/dTdk;
    
    // Evaluate the free energy and its derivatives analytically
    // using the interpolator, and then remove a factor of hbar c

    std::vector<double> out(1);
    eval(index,out);

    double Seg=tgp_S->get(index)-tgp_Sint->get(index);
    double Feg=(tgp_F->get(index)-tgp_Fint->get(index))/hc_mev_fm;
    
    double Fint, F;

    if (interp_Fint) {
      // Interacting free energy per baryon in units of 1/fm
      Fint=out[0]/hc_mev_fm;
      // Total free energy
      F=Fint+Feg;
    } else {
      // Total free energy
      F=out[0]/hc_mev_fm;
      // Interacting free energy per baryon in units of 1/fm
      Fint=F-Feg;
    }

    // Electron chemical potential in 1/fm
    double mue=tgp_mue->get(index)/hc_mev_fm;

    elep.include_deriv=true;
    elep.include_muons=false;
    elep.include_photons=true;
    elep.e.mu=elep.e.m;
    elep.e.n=Ye*nB;
    elep.pair_density_eq(Ye*nB,T_MeV/hc_mev_fm);
    if (addl_verbose>=2) {
      cout << "Feg[elep],Feg[table]: "
           << (elep.th.ed-T_MeV/hc_mev_fm*elep.th.en)/nB << " "
           << Feg << endl;
      cout << "Seg,Se,Sg,Seg[table]: " << elep.th.en/nB << " "
           << elep.e.en/nB << " " 
           << elep.ph.en/nB << " " << Seg << endl;
      cout << "mue,mue[table]: " << elep.e.mu << " "
           << tgp_mue->get(index)/hc_mev_fm << endl;
    }

    double dFint_dnB, dFint_dYe, dFint_dT, dF_dnB, dF_dYe, dF_dT;
    double Fint_nBnB, Fint_nBYe, Fint_YeYe, Fint_nBT, Fint_YeT, Fint_TT;
    double F_nBnB, F_nBYe, F_YeYe, F_nBT, F_YeT, F_TT;

    if (interp_Fint) {
      
      // First derivatives of the free energy
      deriv(index,out,0);
      double dFdi=out[0]/hc_mev_fm;
      dFint_dnB=dFdi*didnB;
      
      deriv(index,out,1);
      double dFdj=out[0]/hc_mev_fm;
      dFint_dYe=dFdj*djdYe;
      
      deriv(index,out,2);
      double dFdk=out[0]/hc_mev_fm;
      dFint_dT=dFdk*dkdT*hc_mev_fm;
      
      // Convert derivatives of Fint to derivatives of F
      dF_dnB=dFint_dnB+Ye*mue/nB-Feg/nB;
      dF_dYe=dFint_dYe+mue;
      dF_dT=dFint_dT-Seg;
      
      deriv2(index,out,0,0);
      double d2Fdi2=out[0]/hc_mev_fm;
      // d2FdnB2 in fm^5
      Fint_nBnB=d2Fdi2*didnB*didnB+dFdi*d2idnB2;
      
      deriv2(index,out,0,1);
      double d2Fdidj=out[0]/hc_mev_fm;
      // d2FdnBdYe in fm^2
      Fint_nBYe=d2Fdidj*didnB*djdYe;
      
      deriv2(index,out,1,1);
      double d2Fdj2=out[0]/hc_mev_fm;
      // d2FdYe2 in fm^{-1}
      Fint_YeYe=d2Fdj2*djdYe*djdYe+dFdj*d2jdYe2;
      
      deriv2(index,out,0,2);
      double d2Fdidk=out[0]/hc_mev_fm;
      Fint_nBT=d2Fdidk*didnB*dkdT*hc_mev_fm;
      
      deriv2(index,out,1,2);
      double d2Fdjdk=out[0]/hc_mev_fm;
      Fint_YeT=d2Fdjdk*djdYe*dkdT*hc_mev_fm;
      
      deriv2(index,out,2,2);
      double d2Fdk2=out[0]/hc_mev_fm;
      Fint_TT=(d2Fdk2*dkdT*dkdT+dFdk*d2kdT2)*hc_mev_fm*hc_mev_fm;
      
      // Convert second derivatives of Fint to second derivatives of F
      F_nBnB=Fint_nBnB+Ye*Ye/nB/elep.ed.dndmu-2.0*Ye*mue/nB/nB+
        2.0*Feg/nB/nB;
      F_nBYe=Fint_nBYe+Ye/elep.ed.dndmu;
      F_YeYe=Fint_YeYe+nB/elep.ed.dndmu;
      F_nBT=Fint_nBT+Ye/nB/elep.ed.dndmu*elep.ed.dndT+Seg/nB;
      F_YeT=Fint_YeT+1.0/elep.ed.dndmu*elep.ed.dndT;
      F_TT=Fint_TT-elep.ed.dsdT/nB;

    } else {

      // First derivatives of the free energy
      deriv(index,out,0);
      double dFdi=out[0]/hc_mev_fm;
      dF_dnB=dFdi*didnB;
      
      deriv(index,out,1);
      double dFdj=out[0]/hc_mev_fm;
      dF_dYe=dFdj*djdYe;
      
      deriv(index,out,2);
      double dFdk=out[0]/hc_mev_fm;
      dF_dT=dFdk*dkdT*hc_mev_fm;
      
      // Convert derivatives of Fint to derivatives of F
      dFint_dnB=dF_dnB-Ye*mue/nB+Feg/nB;
      dFint_dYe=dF_dYe-mue;
      dFint_dT=dF_dT+Seg;
      
      deriv2(index,out,0,0);
      double d2Fdi2=out[0]/hc_mev_fm;
      // d2FdnB2 in fm^5
      F_nBnB=d2Fdi2*didnB*didnB+dFdi*d2idnB2;
      
      deriv2(index,out,0,1);
      double d2Fdidj=out[0]/hc_mev_fm;
      // d2FdnBdYe in fm^2
      F_nBYe=d2Fdidj*didnB*djdYe;
      
      deriv2(index,out,1,1);
      double d2Fdj2=out[0]/hc_mev_fm;
      // d2FdYe2 in fm^{-1}
      F_YeYe=d2Fdj2*djdYe*djdYe+dFdj*d2jdYe2;
      
      deriv2(index,out,0,2);
      double d2Fdidk=out[0]/hc_mev_fm;
      F_nBT=d2Fdidk*didnB*dkdT*hc_mev_fm;
      
      deriv2(index,out,1,2);
      double d2Fdjdk=out[0]/hc_mev_fm;
      F_YeT=d2Fdjdk*djdYe*dkdT*hc_mev_fm;
      
      deriv2(index,out,2,2);
      double d2Fdk2=out[0]/hc_mev_fm;
      F_TT=(d2Fdk2*dkdT*dkdT+dFdk*d2kdT2)*hc_mev_fm*hc_mev_fm;
      
      // Convert second derivatives of Fint to second derivatives of F
      Fint_nBnB=F_nBnB-Ye*Ye/nB/elep.ed.dndmu+2.0*Ye*mue/nB/nB-
        2.0*Feg/nB/nB;
      Fint_nBYe=F_nBYe-Ye/elep.ed.dndmu;
      Fint_YeYe=F_YeYe-nB/elep.ed.dndmu;
      Fint_nBT=F_nBT-Ye/nB/elep.ed.dndmu*elep.ed.dndT-Seg/nB;
      Fint_YeT=F_YeT-1.0/elep.ed.dndmu*elep.ed.dndT;
      Fint_TT=F_TT+elep.ed.dsdT/nB;

    }

    // Use those derivatives to compute the chemical potentials and
    // the entropy density
    double mun=F-Ye*dF_dYe+nB*dF_dnB;
    double mupmue=F+(1.0-Ye)*dF_dYe+nB*dF_dnB;
    double mup=mupmue-mue;
    double en=-nB*dF_dT;
        
    // Compare theose derivatives with the stored values

    if (addl_verbose>=2) {
      std::cout << "  mun table [1/fm], mun interp [1/fm]: "
                << tgp_mun->get(index)/hc_mev_fm << " "
                << mun << endl;
      cout << "    " << F << " " << -Ye*dF_dYe << " "
           << nB*dF_dnB << endl;
      std::cout << "  mup table [1/fm], mup interp [1/fm]: "
                << tgp_mup->get(index)/hc_mev_fm << " "
                << mup << endl;
    }

    std::vector<size_t> im1, ip1, jm1, jp1, km1, kp1;
    std::vector<size_t> ip1jp1, ip1kp1, jp1kp1;
    
    if (index[0]>0) im1={index[0]-1,index[1],index[2]};
    else im1={0,0,0};
    if (index[1]>0) jm1={index[0],index[1]-1,index[2]};
    else jm1={0,0,0};
    if (index[2]>0) km1={index[0],index[1],index[2]-1};
    else km1={0,0,0};
    if (index[0]<nB_grid.size()-1) ip1={index[0]+1,index[1],index[2]};
    else ip1={0,0,0};
    if (index[1]<Ye_grid.size()-1) jp1={index[0],index[1]+1,index[2]};
    else jp1={0,0,0};
    if (index[2]<T_grid.size()-1) kp1={index[0],index[1],index[2]+1};
    else kp1={0,0,0};
    if (index[0]<nB_grid.size()-1 && index[1]<Ye_grid.size()-1) {
      ip1jp1={index[0]+1,index[1]+1,index[2]};
    }
    if (index[0]<nB_grid.size()-1 && index[2]<T_grid.size()-1) {
      ip1kp1={index[0]+1,index[1],index[2]+1};
    }
    if (index[1]<Ye_grid.size()-1 && index[2]<T_grid.size()-1) {
      jp1kp1={index[0],index[1]+1,index[2]+1};
    }
        
    else jp1={0,0,0};

    // Now compute the second derivatives and the speed of sound
        
    double f_nnnn=(Ye*Ye*F_YeYe+nB*(2.0*dF_dnB-2.0*Ye*F_nBYe+nB*F_nBnB))/nB;
    double f_nnnp=((Ye-1.0)*Ye*F_YeYe+nB*(2.0*dF_dnB+(1.0-2.0*Ye)*
                                          F_nBYe+nB*F_nBnB))/nB;
    double f_npnp=((Ye-1.0)*(Ye-1.0)*F_YeYe+nB*(2.0*dF_dnB-2.0*(Ye-1.0)*
                                                F_nBYe+nB*F_nBnB))/nB;
    double f_nnT=dF_dT-Ye*F_YeT+nB*F_nBT;
    double f_npT=dF_dT+(Ye-1.0)*F_YeT+nB*F_nBT;
    double f_TT=nB*F_TT;
    
    double den=en*T_MeV/hc_mev_fm+(mun+mneut)*nB*(1.0-Ye)+
      (mup+mprot)*nB*Ye+mue*nB*Ye;
    if (addl_verbose>=2) {
      std::cout << "mun,mup,mn,mp: " << mun << " " << mup << " " << mneut << " "
                << mprot << std::endl;
      std::cout << "den,sT,mun*nn,mup*np,mue*np: "
                << den << " " << en*T_MeV/hc_mev_fm << " "
                << (mun+mneut)*nB*(1.0-Ye) << " "
                << (mup+mprot)*nB*Ye << " " << mue*nB*Ye << std::endl;
    }
    double nn2=nB*(1.0-Ye);
    double np2=nB*Ye;
    double cs_sq=(nn2*nn2*(f_nnnn-f_nnT*f_nnT/f_TT)+
                  2.0*nn2*np2*(f_nnnp-f_nnT*f_npT/f_TT)+
                  np2*np2*(f_npnp-f_npT*f_npT/f_TT)-
                  2.0*en*(nn2*f_nnT/f_TT+np2*f_npT/f_TT)-en*en/f_TT)/den;
    if (addl_verbose>=2) {
      cout << "f_nnnn, f_nnnp, f_npnp, f_nnT, f_npT, f_TT, den:\n  "
           << f_nnnn << " " << f_nnnp << " " << f_npnp << " "
           << f_nnT << "\n  " << f_npT << " " << f_TT << " "
           << den << endl;
      
      cout.setf(ios::showpos);
      cout << "en,f_TT,den: " << en << " " << f_TT << " " << den << endl;
      std::cout << "  t1,t2,t3,t4,t5,t6,cs2 (interp,table): " 
                << nn2*nn2*(f_nnnn-f_nnT*f_nnT/f_TT)/den << " "
                << 2.0*nn2*np2*(f_nnnp-f_nnT*f_npT/f_TT)/den << " "
                << np2*np2*(f_npnp-f_npT*f_npT/f_TT)/den << " "
                << 2.0*en*nn2*f_nnT/f_TT/den << " " 
                << 2.0*en*np2*f_npT/f_TT/den << " "
                << -en*en/f_TT/den << " " << cs_sq << " ";
      cout << tgp_cs2->get(index) << endl;
    }
        
    // Also compute dPdnB

    if (addl_verbose>=2) {
      if (index[0]>0 && index[0]<nB_grid.size()-1) {
        cout << "  dF_dnB (interp, table lo, table hi): " << dF_dnB << " "
             << (tgp_F->get(index)-tgp_F->get(im1))/hc_mev_fm/
          (nB_grid[index[0]]-nB_grid[index[0]-1]) << " "
             << (tgp_F->get(ip1)-tgp_F->get(index))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]) << endl;
      }
      if (index[1]>0 && index[1]<Ye_grid.size()-1) {
        cout << "  dF_dYe (interp, table lo, table hi): " << dF_dYe << " "
             << (tgp_F->get(index)-tgp_F->get(jm1))/hc_mev_fm/
          (Ye_grid[index[1]]-Ye_grid[index[1]-1]) << " "
             << (tgp_F->get(jp1)-tgp_F->get(index))/hc_mev_fm/
          (Ye_grid[index[1]+1]-Ye_grid[index[1]]) << endl;
      }
      if (index[2]>0 && index[2]<T_grid.size()-1) {
        cout << "  dF_dT (interp, table lo, table hi): " << dF_dT << " "
             << (tgp_F->get(index)-tgp_F->get(km1))/
          (T_grid[index[2]]-T_grid[index[2]-1]) << " "
             << (tgp_F->get(kp1)-tgp_F->get(index))/
          (T_grid[index[2]+1]-T_grid[index[2]]) << endl;
      }

      if (index[0]>0 && index[0]<nB_grid.size()-1) {
        cout << "  dFint_dnB (interp, table lo, table hi): " << dFint_dnB << " "
             << (tgp_Fint->get(index)-tgp_Fint->get(im1))/hc_mev_fm/
          (nB_grid[index[0]]-nB_grid[index[0]-1]) << " "
             << (tgp_Fint->get(ip1)-tgp_Fint->get(index))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]) << endl;
      }
      if (index[1]>0 && index[1]<Ye_grid.size()-1) {
        cout << "  dFint_dYe (interp, table lo, table hi): " << dFint_dYe << " "
             << (tgp_Fint->get(index)-tgp_Fint->get(jm1))/hc_mev_fm/
          (Ye_grid[index[1]]-Ye_grid[index[1]-1]) << " "
             << (tgp_Fint->get(jp1)-tgp_Fint->get(index))/hc_mev_fm/
          (Ye_grid[index[1]+1]-Ye_grid[index[1]]) << endl;
      }
      if (index[2]>0 && index[2]<T_grid.size()-1) {
        cout << "  dFint_dT (interp, table lo, table hi): " << dFint_dT << " "
             << (tgp_Fint->get(index)-tgp_Fint->get(km1))/
          (T_grid[index[2]]-T_grid[index[2]-1]) << " "
             << (tgp_Fint->get(kp1)-tgp_Fint->get(index))/
          (T_grid[index[2]+1]-T_grid[index[2]]) << endl;
      }
 
      cout << "\nSecond derivatives of Fint and F:" << endl;
      if (index[0]>0 && index[0]<nB_grid.size()-1) {
        double t1=(tgp_Fint->get(index)-tgp_Fint->get(im1))/hc_mev_fm/
          (nB_grid[index[0]]-nB_grid[index[0]-1]);
        double t2=(tgp_Fint->get(ip1)-tgp_Fint->get(index))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]);
        double t3=(t2-t1)*2.0/(nB_grid[index[0]+1]-nB_grid[index[0]-1]);
        cout << "  Fint_nBnB,Fint_nBnB_intp: "
             << Fint_nBnB << " " << t3 << endl;
        t1=(tgp_F->get(index)-tgp_F->get(im1))/hc_mev_fm/
          (nB_grid[index[0]]-nB_grid[index[0]-1]);
        t2=(tgp_F->get(ip1)-tgp_F->get(index))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]);
        t3=(t2-t1)*2.0/(nB_grid[index[0]+1]-nB_grid[index[0]-1]);
        cout << "  F_nBnB,F_nBnB_intp: "
             << F_nBnB << " " << t3 << endl;
      }

      if (index[0]>0 && index[0]<nB_grid.size()-1 &&
          index[1]>0 && index[1]<Ye_grid.size()-1) {
        double t1=(tgp_Fint->get(ip1)-tgp_Fint->get(index))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]);
        double t2=(tgp_Fint->get(ip1jp1)-tgp_Fint->get(jp1))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]);
        double t3=(t2-t1)/(Ye_grid[index[1]+1]-Ye_grid[index[1]]);
        cout << "  Fint_nBYe,Fint_nBYe_intp: "
             << Fint_nBYe << " " << t3 << endl;
        t1=(tgp_F->get(ip1)-tgp_F->get(index))/hc_mev_fm/
          (nB_grid[index[0]]-nB_grid[index[0]-1]);
        t2=(tgp_F->get(ip1jp1)-tgp_F->get(jp1))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]);
        t3=(t2-t1)/(Ye_grid[index[1]+1]-Ye_grid[index[1]]);
        cout << "  F_nBYe,F_nBYe_intp: "
             << F_nBYe << " " << t3 << endl;
      }

      if (index[1]>0 && index[1]<Ye_grid.size()-1) {
        double t1=(tgp_Fint->get(index)-tgp_Fint->get(jm1))/hc_mev_fm/
          (Ye_grid[index[1]]-Ye_grid[index[1]-1]);
        double t2=(tgp_Fint->get(jp1)-tgp_Fint->get(index))/hc_mev_fm/
          (Ye_grid[index[1]+1]-Ye_grid[index[1]]);
        double t3=(t2-t1)*2.0/(Ye_grid[index[1]+1]-Ye_grid[index[1]-1]);
        cout << "  Fint_YeYe,Fint_YeYe_intp: " << Fint_YeYe << " "
             << t3 << endl;
        t1=(tgp_F->get(index)-tgp_F->get(jm1))/hc_mev_fm/
          (Ye_grid[index[1]]-Ye_grid[index[1]-1]);
        t2=(tgp_F->get(jp1)-tgp_F->get(index))/hc_mev_fm/
          (Ye_grid[index[1]+1]-Ye_grid[index[1]]);
        t3=(t2-t1)*2.0/(Ye_grid[index[1]+1]-Ye_grid[index[1]-1]);
        cout << "  F_YeYe,F_YeYe_intp: " << F_YeYe << " " << t3 << endl;
      }
        
      if (index[0]>0 && index[0]<nB_grid.size()-1 &&
          index[2]>0 && index[2]<T_grid.size()-1) {
        double t1=(tgp_Fint->get(ip1)-tgp_Fint->get(index))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]);
        double t2=(tgp_Fint->get(ip1kp1)-tgp_Fint->get(kp1))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]);
        double t3=(t2-t1)*hc_mev_fm/
          (T_grid[index[2]+1]-T_grid[index[2]]);
        cout << "  Fint_nBT,Fint_nBT_intp: "
             << Fint_nBT << " " << t3 << endl;
        t1=(tgp_F->get(ip1)-tgp_F->get(index))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]);
        t2=(tgp_F->get(ip1kp1)-tgp_F->get(kp1))/hc_mev_fm/
          (nB_grid[index[0]+1]-nB_grid[index[0]]);
        t3=(t2-t1)*hc_mev_fm/
          (T_grid[index[2]+1]-T_grid[index[2]]);
        cout << "  F_nBT,F_nBT_intp: "
             << F_nBT << " " << t3 << endl;
      }

      if (index[1]>0 && index[1]<Ye_grid.size()-1 &&
          index[2]>0 && index[2]<T_grid.size()-1) {
        double t1=(tgp_Fint->get(jp1)-tgp_Fint->get(index))/hc_mev_fm/
          (Ye_grid[index[1]+1]-Ye_grid[index[1]]);
        double t2=(tgp_Fint->get(jp1kp1)-tgp_Fint->get(kp1))/hc_mev_fm/
          (Ye_grid[index[1]+1]-Ye_grid[index[1]]);
        double t3=(t2-t1)*hc_mev_fm/
          (T_grid[index[2]+1]-T_grid[index[2]]);
        cout << "  Fint_YeT,Fint_YeT_intp: "
             << Fint_YeT << " " << t3 << endl;
        t1=(tgp_F->get(jp1)-tgp_F->get(index))/hc_mev_fm/
          (Ye_grid[index[1]+1]-Ye_grid[index[1]]);
        t2=(tgp_F->get(jp1kp1)-tgp_F->get(kp1))/hc_mev_fm/
          (Ye_grid[index[1]+1]-Ye_grid[index[1]]);
        t3=(t2-t1)*hc_mev_fm/
          (T_grid[index[2]+1]-T_grid[index[2]]);
        cout << "  F_YeT,F_YeT_intp: "
             << F_YeT << " " << t3 << endl;
      }

      if (index[2]>0 && index[2]<T_grid.size()-1) {
        double t1=(tgp_Fint->get(index)-tgp_Fint->get(km1))/
          (T_grid[index[2]]-T_grid[index[2]-1]);
        double t2=(tgp_Fint->get(kp1)-tgp_Fint->get(index))/
          (T_grid[index[2]+1]-T_grid[index[2]]);
        double t3=hc_mev_fm*(t2-t1)*2.0/(T_grid[index[2]+1]-T_grid[index[2]-1]);
        cout << "  Fint_TT,Fint_TT_intp: " << Fint_TT << " " << t3 << endl;
        t1=(tgp_F->get(index)-tgp_F->get(km1))/
          (T_grid[index[2]]-T_grid[index[2]-1]);
        t2=(tgp_F->get(kp1)-tgp_F->get(index))/
          (T_grid[index[2]+1]-T_grid[index[2]]);
        t3=hc_mev_fm*(t2-t1)*2.0/(T_grid[index[2]+1]-T_grid[index[2]-1]);
        cout << "  F_TT,F_TT_intp: " << F_TT << " " << t3 << endl;
      }
      cout << endl;
      
      cout << "  TI1,TI2: " << nB*(1.0-Ye)*mun+(mup+mue)*nB*Ye-
        tgp_P->get(index)/hc_mev_fm << " "
           << F*nB/hc_mev_fm << endl;
    }
      
    double dmun_dnB=f_nnnn*(1.0-Ye)+f_nnnp*Ye;
    double dmupmue_dnB=f_nnnp*(1.0-Ye)+f_npnp*Ye;
    
    double dmun_dYe=nB*(f_nnnp-f_nnnn);
    double dmupmue_dYe=nB*(f_npnp-f_nnnp);
    double ds_dnB=-f_nnT*(1.0-Ye)-f_npT*Ye;
    double ds_dYe=nB*(f_nnT-f_npT);
    
    if (addl_verbose>=2 && index[0]>0 && index[0]<nB_grid.size()-1) {
      std::cout << "  dmun_dnB,dmupmue_dnB,dmun_dnB_intp,"
                << "dmupmue_dnB_intp,dmun_dnB_intp2,dmupmue_dnB_intp2:\n  "
                << dmun_dnB << " " << dmupmue_dnB << " ";
      std::cout << (tgp_mun->get(index)-tgp_mun->get(im1))/hc_mev_fm/
        (nB_grid[index[0]]-nB_grid[index[0]-1]) << " ";
      std::cout << (tgp_mup->get(index)+tgp_mue->get(index)-
                    tgp_mup->get(im1)-tgp_mue->get(im1))/hc_mev_fm/
        (nB_grid[index[0]]-nB_grid[index[0]-1]) << " ";
      std::cout << (tgp_mun->get(ip1)-tgp_mun->get(index))/hc_mev_fm/
        (nB_grid[index[0]+1]-nB_grid[index[0]]) << " ";
      std::cout << (tgp_mup->get(ip1)+tgp_mue->get(ip1)-
                    tgp_mup->get(index)-tgp_mue->get(index))/hc_mev_fm/
        (nB_grid[index[0]+1]-nB_grid[index[0]]) << endl;
    }
    
    double dfdnB=mun*(1.0-Ye)+(mup+mue)*Ye;
    double dPdnB=dmun_dnB*nB*(1.0-Ye)+mun*(1.0-Ye)+
      dmupmue_dnB*nB*Ye+(mup+mue)*Ye-dfdnB;
    if (addl_verbose>=2 && index[0]>0 && index[0]<nB_grid.size()-1) {
      cout << "  dPdnB (interp, table lo, table hi): " << dPdnB << " ";
      cout << (tgp_P->get(index)-tgp_P->get(im1))/hc_mev_fm/
        (nB_grid[index[0]]-nB_grid[index[0]-1]) << " "
           << (tgp_P->get(ip1)-tgp_P->get(index))/hc_mev_fm/
        (nB_grid[index[0]+1]-nB_grid[index[0]]) << std::endl;
      cout.unsetf(ios::showpos);
    }

    if (true) {
      double expr1=nB*nB*dmun_dnB-nB*nB*ds_dnB*ds_dnB/f_TT+
        nB*Ye*(1.0-Ye)*dmun_dYe+nB*Ye*Ye*dmupmue_dYe;
      double expr2=-nB*ds_dnB/f_TT;
      cs_sq=(expr1-2.0*en*expr2-en*en/f_TT)/den;
    }
    
    if (dPdnB<=0.0) {
      if (addl_verbose>=1) {
        cout << "Return failure for dPdnB<=0.0." << endl;
      }
      if (addl_verbose>=2) {
        cout << endl;
      }
      return 2;
    }
    
    if (cs_sq>1.0 || cs_sq<0.0 || !std::isfinite(cs_sq)) {
      if (addl_verbose>=1) {
        cout << "Return failure for unphysical cs_sq." << endl;
      }
      if (addl_verbose>=2) {
        cout << endl;
      }
      return 1;
    }
    
    if (addl_verbose>=2) {
      cout << endl;
    }

    // End of loop over ilist
  }
  
  if (addl_verbose>=1) {
    cout << "Return success." << endl;
  }
  return 0;
}
#endif
