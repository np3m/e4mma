/*
  -------------------------------------------------------------------
  
  Copyright (C) 2022-2023, Andrew W. Steiner, Josue Bautista
  
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

int eos_nuclei::test_hdf5io () {
    std::string cs2in = "/home/jbaut001/st.o2";
    std::string Fintin = "/home/awsteiner/wcs/eos/fid_3_14_23.o2";
    std::string cs2o2 = "/home/jbaut001/test/sttest.o2";
    std::string Finto2 = "/home/jbaut001/test/fidtest.o2";
    std::filesystem::copy_file(cs2in, cs2o2, std::filesystem::copy_options::overwrite_existing);
    std::filesystem::copy_file(Fintin, Finto2, std::filesystem::copy_options::overwrite_existing);

    hdf_file hff, hff2;
    o2scl::tensor_grid<> tgp_cs2, tgp_Fint;
    hff.open_or_create(cs2o2);
    hff2.open_or_create(Finto2);
    hdf_input(hff, tgp_cs2, "cs2");
    hdf_input(hff2, tgp_Fint, "Fint");
    eos_nuclei::change_tgp(tgp_cs2, 3.0);
    eos_nuclei::change_tgp(tgp_Fint, 50.0);
    hdf_output(hff, tgp_cs2, "cs2");
    hdf_output(hff2, tgp_Fint, "Fint");
    hff.seti("derivs_computed", 0);
    hff.seti("with_leptons", 0);
    hff2.seti("derivs_computed", 0);
    hff.seti("with_leptons", 0);
    hff.close();
    hff2.close();
    return 0;
}

void eos_nuclei::change_tgp(o2scl::tensor_grid<>& tg_file, double value) {
    for (size_t inB=210;inB<=213;inB++) {
        for (size_t iYe=2;iYe<=5;iYe++) {
            for (size_t iT=64;iT<=65;iT++) {
                vector<size_t> index={inB, iYe, iT};
                cout << tg_file.get(index) << endl;
                tg_file.get(index)=value;
                cout << tg_file.get(index) << endl;
            }
        }
    }
}

int eos_nuclei::interp_point(std::vector<std::string> &sv,
                             bool itive_com) {
  if (sv.size()<6) {
    cerr << "Not enough arguments interp-point." << endl;
    return 1;
  }

  double nB_cent=o2scl::function_to_double(sv[1]);
  double Ye_cent=o2scl::function_to_double(sv[2]);
  double T_cent=o2scl::function_to_double(sv[3]);

  int window=o2scl::stoi(sv[4]);

  std::string st_o2=sv[5];
  std::string stfix_o2="";
  hdf_file hff, hff2;
  o2scl::tensor_grid<> tgp_cs2, tgp_file;
  hff.open(st_o2);
  if (sv.size()>=7) {
      stfix_o2=sv[6];
      std::filesystem::copy_file(table_path, stfix_o2, std::filesystem::copy_options::overwrite_existing);
      hff2.open_or_create(stfix_o2);
  }
  hdf_input(hff, tgp_cs2, "cs2");
  tgp_file = tg_Fint;

  eos_nuclei::interpolate(nB_cent, Ye_cent, T_cent, window, (window*2), st_o2, tgp_cs2, tgp_file, itive_com);

  if (sv.size()>=7) {
      hdf_output(hff2, tgp_file, "Fint");
      hff2.seti("derivs_computed", 0);
      hff2.seti("with_leptons", 0);
      hff2.close();
  }

  hff.close();
  return 0;
}


void eos_nuclei::interpolate(double nB_p,
                                            double Ye_p,
                                            double T_p,
                                            int window,
                                            int neighborhood,
                                            std::string st_o2,
                                            o2scl::tensor_grid<> &tg_cs2,
                                            o2scl::tensor_grid<> &tg_file,
                                            bool itive_com) {

  if (!loaded) {
    O2SCL_ERR("No EOS loaded in interp_point.",o2scl::exc_einval);
  }
  if (with_leptons==false) {
    O2SCL_ERR("No EOS leptons in interp_point.",o2scl::exc_einval);
  } 

  std::map<std::vector<size_t>, std::vector<double>> results_no_mue;
  std::map<std::vector<size_t>, std::vector<double>> results_with_mue;
  std::map<std::vector<size_t>, std::vector<double>> results_table;
  std::map<std::vector<size_t>, double> results;
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
          if ((ike.tgp_cs2.get_rank()>=3 &&
              (ike.tgp_cs2.get(index)>1.0 ||
               !std::isfinite(ike.tgp_cs2.get(index)) ||
               ike.tgp_cs2.get(index)<0.0)) &&
               (tg_cs2.get_rank()>=3 &&
               (tg_cs2.get(index)>1.0 ||
               !std::isfinite(tg_cs2.get(index)) ||
               tg_cs2.get(index)<0.0))) {
            ike.fix_list.push_back(index[0]);
            ike.fix_list.push_back(index[1]);
            ike.fix_list.push_back(index[2]);
            //cout << "fix: " << index[0] << " " << index[1] << " "
            //<< index[2] << endl;
          } else if (ike.tgp_cs2.get_rank()>=3 &&
              (ike.tgp_cs2.get(index)<=1.0 &&
               std::isfinite(ike.tgp_cs2.get(index)) &&
               ike.tgp_cs2.get(index)>=0.0)) {
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
  ptemp.push_back(len_list);
  ptemp.push_back(len_list);
  ptemp.push_back(len_list);
  ptemp.push_back(len_list);
  ptemp.push_back(len_list);
  ptemp.push_back(len_list);
  ptemp.push_back(l10_list);
  std::vector<std::vector<std::vector<double>>> param_lists;
  param_lists.push_back(ptemp);
  std::cout << "Going to set_covar()." << std::endl;
  vector<mcovar_funct_quad_correl> mfr(1);
  mfr[0].len.resize(9);
  ike.set_covar(mfr,param_lists);
 
  ike.skip_optim=true;
  ike.set(nB_grid2,Ye_grid2,T_grid2,tg_Fint,tg_P,tg_Sint,
          tg_mun,tg_mup,tg_mue,tg_F,tg_S,neutron.m,proton.m);

  latin_hypercube_sampling(10,2000,ike);
  //minimize_parameters(ike);

  // Use the interpolation results to fix points
  cout << "\nstarting to interpolate here: \n";
  if (ike.fix_list.size() != 0) {
    results_no_mue=ike.interpolate_points(ike.fix_list);
    results_table=calculate_table_values(ike.fix_list, ike);
  }
  else {
    results_no_mue=ike.interpolate_points(ike.calib_list);
    results_table=calculate_table_values(ike.calib_list, ike);
  }
  std::map<std::vector<size_t>, double> results_final;

  std::map<std::vector<size_t>, std::vector<double>>::iterator it, it2;
  for (it=results_no_mue.begin();it!=results_no_mue.end();++it) {
    results[it->first]=it->second.at(5);
    results_final[it->first]=it->second.at(0);
  }
  it2=results_table.begin();
  std::ofstream latex{"table.tex", std::ios_base::app};
  for (it=results_no_mue.begin(); it!=results_no_mue.end();++it) {
  cout << "Point: " << it->first.at(0) << " " << it->first.at(1) << " " << it->first.at(2) << endl;
  latex << "Point: " << it->first.at(0) << " " << it->first.at(1) << " " << it->first.at(2) << "\\\\\n";
  latex << nB_grid2[it->first.at(0)] << " " << Ye_grid2[it->first.at(1)] << " " << T_grid2[it->first.at(2)] << endl;
  latex << "\\begin{center}" << endl;
  latex << " \\begin{tabular}{|c|c|c|c|}" << endl;
  latex << "  \\hline" << endl;
  latex << "  Variable & Fint & error & Table \\\\" << endl;
  latex << "  \\hline" << endl;
  double diff = std::abs(std::abs(it2->second.at(5)-it->second.at(5))/it2->second.at(5));
  latex << "  $\\mathrm{cs}^2$ & " << it->second.at(5) << " & " << diff << " & ";
    latex << it2->second.at(5) << " \\\\ " << endl;
    cout << "Here: " << it->second.at(5) << endl;
    cout << "$\\mathrm{cs}^2\\_\\mathrm{tab}$ " << it2->second.at(5) << endl;

    cout << "F\\_intp " << it->second.at(0) << endl;
    cout << "F\\_tab " << it2->second.at(0) << endl;
    diff = std::abs(std::abs(it2->second.at(0)-it->second.at(0))/it2->second.at(0));
    cout << diff << endl;
    if (diff < 0.001) {  
        cout << "success\n";
    }
    else {
        cout << "failure\n";
    }
    cout << "Calculated: mun[MeV],mup[MeV],mue[MeV],S: " << it->second.at(1)*hc_mev_fm << " " << it->second.at(3)*hc_mev_fm << " " << it->second.at(2)*hc_mev_fm << " " << it->second.at(4) << " \n";
    cout << "Stored: mun[MeV],mup[MeV],mue[MeV],S: " << it2->second.at(1) << " ";
    cout << it2->second.at(3) << " ";
    cout << it2->second.at(2) << " " << it2->second.at(4) << " \n";
    diff = (std::abs(std::abs((it2->second.at(1)-(it->second.at(1)*hc_mev_fm)))/it2->second.at(1)));
    cout << diff << endl;
    if (diff<0.01) {
        cout<<"mun: success\n";
    }
    else {
        cout<<"mun: failure\n";
    }
    cout << "dmundnB: " << it->second.at(6) << endl;
    cout << "dmupmuednB: " << it->second.at(7) << endl;
    cout << "dPdnB: " << it->second.at(8) << endl;
    cout << "dmundnB from table: " << it2->second.at(6) << endl;
    cout << "dmupmuednB from table: " << it2->second.at(7) << endl;
    cout << "dPdnB from table: " << it2->second.at(8) << endl;
    latex << "  \\hline" << endl;
    diff = std::abs(std::abs(it2->second.at(8)-it2->second.at(8))/it2->second.at(8));
    latex << "  dPdnB & " << it2->second.at(8) << " & " << diff << " & ";
    latex << it->second.at(8) << it2->second.at(8) << " \\\\ " << endl;
    latex << "  \\hline" << endl;
    diff = std::abs(std::abs(it2->second.at(0)-it->second.at(0))/it2->second.at(0));
    latex << " Fint & " << it->second.at(0) << " & " << diff << " & ";
    latex << it2->second.at(0) << " \\\\ " << endl;
    latex << "  \\hline" << endl;
    latex << " \\end{tabular}" << endl;
    latex << "\\end{center}" << endl;
    cout << diff << endl;
    if (it->second.at(8)>0.0 && diff<0.10) {
        cout<<"success\n";
    }
    else {
        cout<<"failure\n";
    }
    ++it2;
  }
  //begin fixing points near acausal cs_sq.
  if (ike.fix_list.size() != 0) {
    std::map<std::vector<size_t>, double> external_acausal_points;
    std::map<std::vector<size_t>, std::pair<double, double>> fix_list;
    std::pair<double, double> closest;
    double nearest_internal=0.0;
    double nearest_external=0.0;
    bool isEmpty=false;

    for (int dnB=-(neighborhood*2); dnB<(neighborhood*2); dnB++) {
      for (int dYe=-(neighborhood*2); dYe<(neighborhood*2); dYe++) {
        for (int dT=-(neighborhood*2); dT<(neighborhood*2); dT++) {
          std::vector<size_t> index = {inB+dnB,iYe+dYe,iT+dT};
          if (std::abs(dnB)+std::abs(dYe)+std::abs(dT)<=(neighborhood*2) &&
              inB+dnB>=0 && inB+dnB<n_nB2 &&
              iYe+dYe>=0 && iYe+dYe<n_Ye2 &&
              iT+dT>=0 && iT+dT<n_T2) {
                  if (std::abs(dnB)+std::abs(dYe)+std::abs(dT)>window &&
                    inB+dnB>=0 && inB+dnB<n_nB2 &&
                    iYe+dYe>=0 && iYe+dYe<n_Ye2 &&
                    iT+dT>=0 && iT+dT<n_T2) {
                    if ((ike.tgp_cs2.get_rank()>=3 &&
                    (ike.tgp_cs2.get(index)>1.0 ||
                    !std::isfinite(ike.tgp_cs2.get(index)) ||
                    ike.tgp_cs2.get(index)<0.0))) {
                      external_acausal_points[index]=0.0;
                    }
                  }
                  else if (std::abs(dnB)+std::abs(dYe)+std::abs(dT)<=window &&
                    inB+dnB>=0 && inB+dnB<n_nB2 &&
                    iYe+dYe>=0 && iYe+dYe<n_Ye2 &&
                    iT+dT>=0 && iT+dT<n_T2) {
                    if ((ike.tgp_cs2.get_rank()>=3 &&
                    (ike.tgp_cs2.get(index)>1.0 ||
                    !std::isfinite(ike.tgp_cs2.get(index)) ||
                    ike.tgp_cs2.get(index)<0.0)) &&
                    (tg_cs2.get_rank()>=3 &&
                    (tg_cs2.get(index)<=1.0 &&
                    std::isfinite(tg_cs2.get(index)) &&
                    tg_cs2.get(index)>=0.0))) {
                      external_acausal_points[index]=0.0;
                    }
                  }
                }
            }
        }
    }

    for (int dnB=-neighborhood; dnB<neighborhood; dnB++) {
      for (int dYe=-neighborhood; dYe<neighborhood; dYe++) {
        for (int dT=-neighborhood; dT<neighborhood; dT++) {
          std::vector<size_t> index = {inB+dnB,iYe+dYe,iT+dT};
          if (std::abs(dnB)+std::abs(dYe)+std::abs(dT)<=neighborhood &&
              inB+dnB>=0 && inB+dnB<n_nB2 &&
              iYe+dYe>=0 && iYe+dYe<n_Ye2 &&
              iT+dT>=0 && iT+dT<n_T2) {
                isEmpty=false;
                if ((ike.tgp_cs2.get_rank()>=3 &&
                    (ike.tgp_cs2.get(index)<=1.0 ||
                    std::isfinite(ike.tgp_cs2.get(index)) ||
                    ike.tgp_cs2.get(index)>=0.0))) {
                    if (!external_acausal_points.empty()) {
                      closest=vector_distance<double>(index, external_acausal_points);
                      nearest_external=closest.first;
                    }
                    else {
                      nearest_external=0.0;
                      isEmpty=true;
                    }
                    closest=vector_distance<double>(index, results_final);
                    nearest_internal=closest.first;
                    if ((nearest_internal<=nearest_external) || (isEmpty==true)) {
                        fix_list[index]=std::make_pair(closest.second,nearest_internal);
                    }
                }
          }
        }
      }
    }
    
    std::map<std::vector<size_t>,double>::iterator it3;
    for (it3=results_final.begin(); it3 != results_final.end(); ++it3) {
      fix_list[it3->first]=std::make_pair(it3->second,0.0);
    }

    double fixed = 0.0;
    double eta = 0.2;
    std::map<std::vector<size_t>, std::pair<double, double>>::iterator it4;
    for (it4=fix_list.begin(); it4 != fix_list.end(); ++it4) {
      fixed=ike.tgp_F->get(it4->first)+((it4->second.first-ike.tgp_F->get(it4->first))*std::exp(-std::pow(it4->second.second, 2.0)/std::pow(eta, 2.0)));
      cout << tg_file.get(it4->first) << " " << fixed << " Original: " << ike.tgp_F->get(it4->first) << " Closest: " << it4->second.first << " Distance: " << it4->second.second << endl;
      tg_file.get(it4->first)=fixed;
      cout << tg_file.get(it4->first) << endl;
    }
    std::map<std::vector<size_t>, double>::iterator it5;
    for (it5=results.begin(); it5 != results.end(); ++it5) {
      tg_cs2.get(it5->first)=it5->second;
    }
  }
}

int eos_nuclei::interp_file(std::vector<std::string> &sv,
                            bool itive_com) {
  if (sv.size()<4) {
    cerr << "Not enough arguments interp-file." << endl;
    return 1;
  }
  
  double nB = 0.0;
  double Ye = 0.0;
  double T_MeV = 0.0;

  std::string st_o2=sv[1];
  std::string stfix_o2=sv[2];
  int window=o2scl::stoi(sv[3]);
  std::string csv_path = "";
  if (sv.size()>=5) {
    csv_path = sv[4];
  }
  hdf_file hff;
  hdf_file hff2;
  o2scl::tensor_grid<> tgp_cs2, tgp_cs2_old, tgp_file;
  hff.open(st_o2);
  hdf_input(hff, tgp_cs2, "cs2");
  hdf_input(hff, tgp_cs2_old, "cs2");
  tgp_file=tg_Fint;
  std::filesystem::copy_file(table_path, stfix_o2, std::filesystem::copy_options::overwrite_existing);
  hff2.open_or_create(stfix_o2);
  std::ifstream csvfile;
  if (sv.size()==4) {
    for (size_t inB=0; inB<nB_grid2.size(); inB++) {
      for (size_t iYe=0; iYe<Ye_grid2.size(); iYe++) {
        for (size_t iT=0; iT<T_grid2.size(); iT++) {
          std::vector<size_t> index = {inB, iYe, iT};
          if (tgp_cs2_old.get_rank()>=3 &&
              (tgp_cs2_old.get(index)>1.0 ||
              !std::isfinite(tgp_cs2_old.get(index)) ||
              tgp_cs2_old.get(index)<0.0) &&
              nB_grid2[inB]<max_nB_inter) {
              std::vector<double> point = {nB_grid2[inB], Ye_grid2[iYe], T_grid2[iT]};
            eos_nuclei::interpolate(point[0], point[1], point[2], window, (window*2), st_o2, tgp_cs2, tgp_file, itive_com);
          }
        }
      }
    }
  }
  else {
    csvfile.open(csv_path);
    std::string line;
    int x;
    std::string val;
    while (std::getline(csvfile, line)) {
        std::stringstream s(line);
        if (true) {
            std::vector<double> file_point;
            x=0;
            while (std::getline(s, val, ',')) {
                file_point.push_back(std::stod(val));
                x++;
                if (s.peek()==',') {
                    s.ignore();
                }
            }
            eos_nuclei::interpolate(file_point[0], file_point[1], file_point[2], window, (window*2), st_o2, tgp_cs2, tgp_file, itive_com);
        }
        else {
            cout<<"csv file of superliminal points must have 3 terms in each line";
        }
    }
  }

  hdf_output(hff2,tgp_file,"Fint");
  hff2.seti("derivs_computed", 0);
  hff2.seti("with_leptons", 0);
  hff.close();
  hff2.close();
  return 0;
}

  template<typename T>
  std::pair<double, T> eos_nuclei::vector_distance(std::vector<size_t> start, 
                                                  std::map<std::vector<size_t>, T> points) {
    std::pair<double, T> results;
    double distance=0.0;
    typename std::map<std::vector<size_t>, T>::iterator i;
    for (i=points.begin(); i != points.end(); ++i) {
        if ((start.size() != i->first.size()) || (start.size()==0)) {
            results = std::make_pair(0.0,T());
            return results;
        }
        for (size_t x=0; x<start.size(); x++) {
            if (i->first.at(x)>start.at(x)) {
                distance+=std::pow((i->first.at(x)-start.at(x)), 2.0);
            }
            else {
                distance+=std::pow((start.at(x)-i->first.at(x)), 2.0);
            }
        }
        distance=std::sqrt(distance);
        if ((results.first==0.0) || (distance<results.first)) {
            results.first=distance;
            results.second=i->second;
        }
    }
    return results;
}

//to make more generic consider returning list of strata and handling generating hyperparameters in another function.
void eos_nuclei::latin_hypercube_sampling(int dim, size_t samples, interpm_krige_eos& ike) {
    cout << "Attempting Latin Hypercube Sampling" << std::endl;
    double min_qual=1.0e99;
    bool unique=true;
    vector<double> pnt(dim), min_p;
    vector<vector<int>> strata(dim);
    int iter=0;
    double rangehp1=20.0-2.0;
    double rangehp2=std::abs(-15.0-(-8.99));
    double mar_prob=1.0/samples;
    std::random_device rd;
    std::mt19937 gen(rd());
    //std::srand(time(NULL));
    //consider changing vectors to strings with values separated by a delimiting character here
    //Works as long as there are no repeating numbers for each dimension
    for (int x=0;x<strata.size();x++) {
        vector<int> numrange(samples);
        for (int y=0;y<numrange.size();y++) {
            numrange[y]=y;
        }
        shuffle (numrange.begin(),numrange.end(),gen);
        strata[x]=numrange;
        cout << "Done with dim " << x << std::endl;
    }
    /*
    while (strata.size() != samples) {
        vector<int> entry(dim);
        for (int x=0;x<entry.size();x++) {
            entry[x]=gen() % samples;
        }
        for (int y=0;y<strata.size();y++) {
            for (int z=0;z<strata[y].size();z++) {
                if (entry[z] == strata[y][z]) {
                    unique=false;
                }
            }
        }
        if (unique==true) {
            strata.push_back(entry);
            cout << "Added entry: " << strata.size() << std::endl;
        }
        unique=true;
    }*/
    //maybe shuffle strata
    while (iter<samples) {
        //assign strata to point
        //randomly sample each hparam for the strata assigned to each dim in the point
        //vector<int> lhstrata = strata.back();
        //strata.pop_back();
        for (int x=0;x<pnt.size();x++) {
            if (x<(pnt.size()-1)) {
                std::uniform_real_distribution<double> dist((2.0+(rangehp1*mar_prob*strata[x].back())),(2.0+(rangehp1*mar_prob*(strata[x].back()+1))));
                pnt[x]=dist(gen);
            }
            else {
                std::uniform_real_distribution<double> dist((-15.0+(rangehp2*mar_prob*strata[x].back())),(-15.0+(rangehp2*mar_prob*(strata[x].back()+1))));
                pnt[x]=dist(gen);
            }
            cout << pnt[x] << " " << strata[x].back() << std::endl;
            strata[x].pop_back();
        }
        //test with min_qual func
        vector_out(cout,pnt,true);
        (*ike.cf)[0].set_params(pnt);
        int success;
        double q=ike.qual_fun(0,success);
        cout << q << " " << success << endl;
        if (q<min_qual) {
            min_p=pnt;
        }
        iter++;
    }
    vector_out(cout,min_p,true);
}

void eos_nuclei::minimize_parameters(interpm_krige_eos& ike) {
  double min_qual=1.0e99;
  vector<double> pnt(10), min_p;  
  for(pnt[0]=2.0;pnt[0]<20.0;pnt[0]*=1.5) {
    for(pnt[1]=2.0;pnt[1]<20.0;pnt[1]*=1.5) {
      for(pnt[2]=2.0;pnt[2]<20.0;pnt[2]*=1.5) {
        for(pnt[3]=2.0;pnt[3]<20.0;pnt[3]*=1.5) {
          for(pnt[4]=2.0;pnt[4]<20.0;pnt[4]*=1.5) {
            for(pnt[5]=2.0;pnt[5]<20.0;pnt[5]*=1.5) {
              for(pnt[6]=2.0;pnt[6]<20.0;pnt[6]*=1.5) {
                for(pnt[7]=2.0;pnt[7]<20.0;pnt[7]*=1.5) {
                  for(pnt[8]=2.0;pnt[8]<20.0;pnt[8]*=1.5) {
                    for(pnt[9]=-15.0;pnt[9]<-8.99;pnt[9]+=2.0) {
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
              }
            }
          }
        }
      }
    }
  vector_out(cout,min_p,true);
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
                            o2scl::tensor_grid<> &tg_Fall,
                            o2scl::tensor_grid<> &tg_Sall,
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
  tgp_Fall=&tg_Fall;
  tgp_Sall=&tg_Sall;

  size_t n_nB=nB_grid2.size();
  size_t n_Ye=Ye_grid2.size();
  size_t n_T=T_grid2.size();

  includeMue=false;
    
  std::cout << "Going to set_data()." << std::endl;
  verbose=2;
  set_data(3,1,calib_list.size()/3,ix,iy);

  return;
}

std::map<std::vector<size_t>, std::vector<double>> eos_nuclei::calculate_table_values(std::vector<size_t> points_list, interpm_krige_eos& ike) {
  std::map<std::vector<size_t>, std::vector<double>> results_map;
  for(size_t j=0;j<points_list.size();j+=3) {
    std::vector<double> results;
    size_t pnB=points_list[j];
    size_t pYe=points_list[j+1];
    size_t pT=points_list[j+2];

    vector<size_t> index={(size_t) pnB,(size_t) pYe,(size_t) pT};

    double nB=nB_grid2[pnB];
    double Ye=Ye_grid2[pYe];
    double T_MeV=T_grid2[pT];
    
    results.push_back(tg_Fint.get(index));
    results.push_back(tg_mun.get(index));
    results.push_back(tg_mue.get(index));
    results.push_back(tg_mup.get(index));
    results.push_back(tg_S.get(index));
    results.push_back(ike.tgp_cs2.get(index));

    vector<size_t> im1={index[0]-1,index[1],index[2]};
    if (!(index[0]==0)){
        double dmundnB = (tg_mun.get(index)-tg_mun.get(im1))/hc_mev_fm/(nB_grid2[index[0]]-nB_grid2[index[0]-1]);
        double dmupmuednB = ((tg_mup.get(index)+tg_mue.get(index))-(tg_mup.get(im1)+tg_mue.get(im1)))/hc_mev_fm/(nB_grid2[index[0]]-nB_grid2[index[0]-1]);
        double dPdnB=(tg_P.get(index)-tg_P.get(im1))/hc_mev_fm/(nB_grid2[index[0]]-nB_grid2[index[0]-1]);
        results.push_back(dmundnB);
        results.push_back(dmupmuednB);
        results.push_back(dPdnB);
    }
    else {
        results.push_back(0.0);
        results.push_back(0.0);
        results.push_back(0.0);
    }
    results_map[index]=results;
  }
  return results_map;
}

std::map<std::vector<size_t>, std::vector<double>> interpm_krige_eos::interpolate_points(std::vector<size_t> points_list) {
  std::vector<double> out(1);
  std::map<std::vector<size_t>, std::vector<double>> results_map;
  for(size_t j=0;j<points_list.size();j+=3) {
    std::vector<double> results;
    size_t pnB=points_list[j];
    size_t pYe=points_list[j+1];
    size_t pT=points_list[j+2];

    vector<size_t> index={(size_t) pnB,(size_t) pYe,(size_t) pT};

    double mue=tgp_mue->get(index)/hc_mev_fm;
    double mue_keep=mue;
    if (includeMue==false) {
        mue=0.0;
    }
    double Feg=tgp_Fall->get(index)-tgp_F->get(index);
    double Seg=tgp_Sall->get(index)-tgp_S->get(index);

    double nB=nB_grid[pnB];
    double Ye=Ye_grid[pYe];
    double T_MeV=T_grid[pT];
    
    // Derivatives of the physical coordinates with respect to the indices
    // new derivation of these variables taken from Dr. Steiner's branch.
    double didnB, djdYe, dkdT, dnBdi, dYedj, dTdk;
    double d2idnB2, d2jdYe2, d2kdT2;
    if (pnB>0 && pnB<nB_grid.size()-1) {
      dnBdi=(nB_grid[pnB+1]-nB_grid[pnB-1])/2;
      double t2=(nB_grid[pnB+1]-nB_grid[pnB]);
      double t1=(nB_grid[pnB]-nB_grid[pnB-1]);
      d2idnB2=(1.0/t2-1.0/t1)/(nB_grid[pnB+1]-nB_grid[pnB-1])*2.0;
    } else if (pnB==0) {
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
    if (pYe>0 && pYe<Ye_grid.size()-1) {
      dYedj=(Ye_grid[pYe+1]-Ye_grid[pYe-1])/2;
      double t2=(Ye_grid[pYe+1]-Ye_grid[pYe]);
      double t1=(Ye_grid[pYe]-Ye_grid[pYe-1]);
      d2jdYe2=(1.0/t2-1.0/t1)/(Ye_grid[pYe+1]-Ye_grid[pYe-1])*2.0;
    } else if (pYe==0) {
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
    if (pT>0 && pT<T_grid.size()-1) {
      dTdk=(T_grid[pT+1]-T_grid[pT-1])/2;
      double t2=(T_grid[pT+1]-T_grid[pT]);
      double t1=(T_grid[pT]-T_grid[pT-1]);
      d2kdT2=(1.0/t2-1.0/t1)/(T_grid[pT+1]-T_grid[pT-1])*2.0;
    } else if (pT==0) {
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
    // using the interpolator
    eval(index,out);
    double Fintp=out[0];
    results.push_back(Fintp);
        
    deriv(index,out,0);
    double dFdi=out[0]/hc_mev_fm;
    double dF_dnB=dFdi*didnB;
    if (interp_Fint) {
        dF_dnB=dF_dnB+((Ye*mue)/nB)-(Feg/nB);
    }
    deriv(index,out,1);
    double dFdj=out[0]/hc_mev_fm;
    double dF_dYe=dFdj*djdYe;
    if (interp_Fint) {
        dF_dYe=dF_dYe+mue;
    }
    deriv(index,out,2);
    double dFdk=out[0]/hc_mev_fm;
    double dF_dT=dFdk*dkdT*hc_mev_fm;
    if (interp_Fint) {
        dF_dT=dF_dT-Seg;
    }

    //This section of code calculates second derivatives used to add contribtion from photons and electrons back to Fint. Taken from Dr. Steiner's code.
    if (true) {
      //cout << "Here." << endl;
      elep.include_deriv=true;
      elep.include_muons=false;
      elep.include_photons=true;
      elep.e.mu=elep.e.m;
      elep.e.n=Ye*nB;
      //cout << elep.e.n << " " << elep.e.mu << endl;
      elep.pair_density_eq(Ye*nB,T_MeV/hc_mev_fm);
      //cout << elep.e.n << " " << elep.e.mu << " "
      //<< elep.e.ed << " " << elep.ed.dndmu << endl;
      cout << "Feg[elep],Feg[table]: "
           << (elep.th.ed-T_MeV/hc_mev_fm*elep.th.en)/nB << " "
           << Feg << endl;
      cout << "Seg,Se,Sg,Seg[table]: " << elep.th.en/nB << " "
           << elep.e.en/nB << " " 
           << elep.ph.en/nB << " " << Seg << endl;
      cout << "mue,mue[table]: " << elep.e.mu << " "
           << tgp_mue->get(index)/hc_mev_fm << endl;
    }

    deriv2(index,out,0,0);
    double d2Fdi2=out[0]/hc_mev_fm;
    double F_nBnB=d2Fdi2*didnB*didnB+dFdi*d2idnB2;
    //added !interp_Fint originally for some reason?
    if (interp_Fint) {
        F_nBnB=F_nBnB+(((Ye*Ye)/nB)*(1/elep.ed.dndmu))-((2*Ye*mue)/(nB*nB))+((2*Feg)/(nB*nB));
    } 
    deriv2(index,out,0,1);
    double d2Fdidj=out[0]/hc_mev_fm;
    double F_nBYe=d2Fdidj*didnB*djdYe;
    if (interp_Fint) {
        F_nBYe=F_nBYe+(Ye*(1/elep.ed.dndmu));
    } 
    deriv2(index,out,1,1);
    double d2Fdj2=out[0]/hc_mev_fm;
    double F_YeYe=d2Fdj2*djdYe*djdYe+dFdj*d2jdYe2;
    if (interp_Fint) {
        F_YeYe=F_YeYe+(nB*(1/elep.ed.dndmu));
    }
    deriv2(index,out,0,2);
    double d2Fdidk=out[0]/hc_mev_fm;
    double F_nBT=d2Fdidk*didnB*dkdT*hc_mev_fm;
    if (interp_Fint) {
        F_nBT=F_YeYe+((Ye/nB)*(elep.ed.dndT/elep.ed.dndmu))+(Seg/nB);
    }
    deriv2(index,out,1,2);
    double d2Fdjdk=out[0]/hc_mev_fm;
    double F_YeT=d2Fdjdk*djdYe*dkdT*hc_mev_fm;
    if (interp_Fint) {
        F_YeT=F_YeT+(elep.ed.dndT/elep.ed.dndmu);
    }
    deriv2(index,out,2,2);
    double d2Fdk2=out[0]/hc_mev_fm;
    double F_TT=(d2Fdk2*dkdT*dkdT+dFdk*d2kdT2)*hc_mev_fm*hc_mev_fm;
    if (interp_Fint) {
        F_TT=F_TT-((1/nB)*(elep.ed.dsdT));
    }
    //changed order of results added to vector in original for some reason? 
    double mun=Fintp/hc_mev_fm+(nB*(dF_dnB+((Ye*mue)/nB)-(Feg/nB)))-(Ye*(dF_dYe+mue));
    results.push_back(mun);
    results.push_back(mue_keep);
    double mup=Fintp/hc_mev_fm+(1.0-Ye)*dF_dYe+nB*dF_dnB-mue;
    results.push_back(mup);
    double en=-nB*dF_dT;
    results.push_back(en);
    //tg_mun.get(index)=mun;
    //tg_mup.get(index)=mup;
    //tg_mue.get(index)=mue;
    
    //tg_S.get(index)=en/nB;

    // unverified
    //tg_E.get(index)=tg_F.get(index)+T_MeV*tg_S.get(index);
    //tg_P.get(index)=tg_F.get(index)+mun*neutron.m+mup*proton.m+
      //mue*electron.m;
//  code for speed of sound squared taken from main/eos_interp.cpp
    double f_nnnn=(Ye*Ye*F_YeYe+nB*(2.0*dF_dnB-2.0*Ye*F_nBYe+nB*F_nBnB))/nB;
    double f_nnnp=((Ye-1.0)*Ye*F_YeYe+nB*(2.0*dF_dnB+(1.0-2.0*Ye)*
                                          F_nBYe+nB*F_nBnB))/nB;
    double f_npnp=((Ye-1.0)*(Ye-1.0)*F_YeYe+nB*(2.0*dF_dnB-2.0*(Ye-1.0)*
                                                F_nBYe+nB*F_nBnB))/nB;
    //alternate f_nnnn, f_nnnp, f_npnp
    //double f_npnp=F_nBnB-(((1-Ye)/(nB*nB))*dF_dYe)-((((1-Ye)*(1-Ye))/(nB*nB))*F_YeYe);
    //double f_nnnn=F_nBnB-((Ye/(nB*nB))*dF_dYe)-(((Ye*Ye)/(nB*nB))*F_YeYe);
    //double f_nnnp=0.0;
    double f_nnT=dF_dT-Ye*F_YeT+nB*F_nBT;
    double f_npT=dF_dT+(Ye-1.0)*F_YeT+nB*F_nBT;
    double f_TT=nB*F_TT;

    double den=en*T_MeV/hc_mev_fm+(mun+mneut)*nB*(1.0-Ye)+
        (mup+mprot)*nB*Ye+mue*nB*Ye;

    double nn2=nB*(1.0-Ye);
    double np2=nB*Ye;

    //maybe change back?
    double dmundnB=(f_nnnn*(1-Ye))+(f_nnnp*Ye);
    double dmundYe=nB*(f_nnnp-f_nnnn);
    double dmupmuedYe=nB*(f_npnp-f_nnnp);
    double dsdnB=-f_nnT*(1.0-Ye)-f_npT*Ye;
    double dsdYe=nB*(f_nnT-f_npT);
    double e1=(nB*nB*dmundnB)+((nB*nB*dsdnB*dsdnB)/f_TT)+(nB*Ye*(1-Ye)*dmundYe)+(nB*Ye*Ye*dmupmuedYe);
    double e2=(nB*dsdnB)/f_TT;
    double cs_sq=e1-(2*en*e2)-((en*en)/f_TT);
/*    double cs_sq=(nn2*nn2*(f_nnnn-f_nnT*f_nnT/f_TT)+
                  2.0*nn2*np2*(f_nnnp-f_nnT*f_npT/f_TT)+
                  np2*np2*(f_npnp-f_npT*f_npT/f_TT)-
                  2.0*en*(nn2*f_nnT/f_TT+np2*f_npT/f_TT)-en*en/f_TT)/den;*/
    results.push_back(cs_sq);

    // Attempt at calculating dP/dnB
    // code taken from main/eos_interp.cpp
    //double dmundnB=(f_nnnn*(1-Ye))+(f_nnnp*Ye);
    results.push_back(dmundnB);
    double dmupmuednB=(f_nnnp*(1-Ye))+(f_npnp*Ye);
    results.push_back(dmupmuednB);
    //double dmundnB=F_nBnB-(Ye*(((1/nB)*F_nBYe)-((1/(nB*nB))*dF_dYe)));
    //double dmupmuednB=F_nBnB-((1-Ye)*(((1/nB)*F_nBYe)-((1/(nB*nB))*dF_dYe)));
    double dPdnB=(dmundnB*nB*(1-Ye))+(mun*(1-Ye))+(dmupmuednB*Ye*nB)+(Ye*(mup+mue))-((mun*(1-Ye))+((mup+mue)*Ye));
    results.push_back(dPdnB);
    results_map[index]=results;
  }
  return results_map;
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
      ilist+=3) {
      
    std::cout << ilist << " ";

    size_t inB, iYe, iT;
    if (ilist<calib_list.size()/3) {
      inB=calib_list[ilist*3];
      iYe=calib_list[ilist*3+1];
      iT=calib_list[ilist*3+2];
      cout << "calib" << endl;
    } else {
      inB=fix_list[(ilist-calib_list.size()/3)*3];
      iYe=fix_list[(ilist-calib_list.size()/3)*3+1];
      iT=fix_list[(ilist-calib_list.size()/3)*3+2];
      cout << "fix" << endl;
    }
    //inB=8;
    //iYe=48;
    //iT=1;
      
    std::vector<size_t> index={((size_t)inB),((size_t)iYe),((size_t)iT)};
    cout << index[0] << " " << index[1] << " "
                << index[2] << std::endl;

        
    double nB=nB_grid[inB];
    double Ye=Ye_grid[iYe];
    double T_MeV=T_grid[iT];

    cout << nB << " " << Ye << " " << T_MeV << " ";

    double mue=tgp_mue->get(index)/hc_mev_fm;
    //not in originally
    if (includeMue==false) {
        mue=0.0;
    }
    double Feg=tgp_Fall->get(index)-tgp_F->get(index);
    double Seg=tgp_Sall->get(index)-tgp_S->get(index);
    
    // Derivatives of the physical coordinates with respect to the indices
    // new derivation of these variables taken from Dr. Steiner's branch.
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
    // using the interpolator

    std::vector<double> out(1);
    eval(index,out);
    double Fintp=out[0];
        
    deriv(index,out,0);
    double dFdi=out[0]/hc_mev_fm;
    double dF_dnB=dFdi*didnB;
    if (interp_Fint) {
        dF_dnB=dF_dnB+((Ye*mue)/nB)-(Feg/nB);
    }
    deriv(index,out,1);
    double dFdj=out[0]/hc_mev_fm;
    double dF_dYe=dFdj*djdYe;
    if (interp_Fint) {
        dF_dYe=dF_dYe+mue;
    }
    deriv(index,out,2);
    double dFdk=out[0]/hc_mev_fm;
    double dF_dT=dFdk*dkdT*hc_mev_fm;
    if (interp_Fint) {
        dF_dT=dF_dT-Seg;
    }

    //This section of code calculates second derivatives used to add contribtion from photons and electrons back to Fint. Taken from Dr. Steiner's code.
    if (true) {
      //cout << "Here." << endl;
      elep.include_deriv=true;
      elep.include_muons=false;
      elep.include_photons=true;
      elep.e.mu=elep.e.m;
      elep.e.n=Ye*nB;
      //cout << elep.e.n << " " << elep.e.mu << endl;
      elep.pair_density_eq(Ye*nB,T_MeV/hc_mev_fm);
      //cout << elep.e.n << " " << elep.e.mu << " "
      //<< elep.e.ed << " " << elep.ed.dndmu << endl;
      cout << "Feg[elep],Feg[table]: "
           << (elep.th.ed-T_MeV/hc_mev_fm*elep.th.en)/nB << " "
           << Feg << endl;
      cout << "Seg,Se,Sg,Seg[table]: " << elep.th.en/nB << " "
           << elep.e.en/nB << " " 
           << elep.ph.en/nB << " " << Seg << endl;
      cout << "mue,mue[table]: " << elep.e.mu << " "
           << tgp_mue->get(index)/hc_mev_fm << endl;
    }

    deriv2(index,out,0,0);
    double d2Fdi2=out[0]/hc_mev_fm;
    double F_nBnB=d2Fdi2*didnB*didnB+dFdi*d2idnB2;
    if (interp_Fint) {
        F_nBnB=F_nBnB+(((Ye*Ye)/nB)*(1/elep.ed.dndmu))-((2*Ye*mue)/(nB*nB))+((2*Feg)/(nB*nB));
    }
    deriv2(index,out,0,1);
    double d2Fdidj=out[0]/hc_mev_fm;
    double F_nBYe=d2Fdidj*didnB*djdYe;
    if (interp_Fint) {
        F_nBYe=F_nBYe+(Ye*(1/elep.ed.dndmu));
    }
    deriv2(index,out,1,1);
    double d2Fdj2=out[0]/hc_mev_fm;
    double F_YeYe=d2Fdj2*djdYe*djdYe+dFdj*d2jdYe2;
    if (interp_Fint) {
        F_YeYe=F_YeYe+(nB*(1/elep.ed.dndmu));
    }
    deriv2(index,out,0,2);
    double d2Fdidk=out[0]/hc_mev_fm;
    double F_nBT=d2Fdidk*didnB*dkdT*hc_mev_fm;
    if (interp_Fint) {
        F_nBT=F_nBT+((Ye/nB)*(elep.ed.dndT/elep.ed.dndmu))+(Seg/nB);
        //F_nBT=F_YeYe+((Ye/nB)*(elep.ed.dndT/elep.ed.dndmu))+(Seg/nB);
    }
    deriv2(index,out,1,2);
    double d2Fdjdk=out[0]/hc_mev_fm;
    double F_YeT=d2Fdjdk*djdYe*dkdT*hc_mev_fm;
    if (interp_Fint) {
        F_YeT=F_YeT+(elep.ed.dndT/elep.ed.dndmu);
    } 
    deriv2(index,out,2,2);
    double d2Fdk2=out[0]/hc_mev_fm;
    double F_TT=(d2Fdk2*dkdT*dkdT+dFdk*d2kdT2)*hc_mev_fm*hc_mev_fm;
    if (interp_Fint) {
        F_TT=F_TT-((1/nB)*(elep.ed.dsdT));
    }

    // Use those derivatives to compute the chemical potentials and
    // the entropy density

    double mun=Fintp/hc_mev_fm+(nB*(dF_dnB+((Ye*mue)/nB)-(Feg/nB)))-(Ye*(dF_dYe+mue));
    double mup=Fintp/hc_mev_fm+(1.0-Ye)*dF_dYe+nB*dF_dnB-mue;
    double en=-nB*dF_dT;
        
    // Compare theose derivatives with the stored values

    if (true && compare) {
      //std::cout << "Indices: " << index[0] << " " << index[1] << " "
      //          << index[2] << std::endl;
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
    /*
    double f_nnnn=(Ye*Ye*F_YeYe+nB*(2.0*dF_dnB-2.0*Ye*F_nBYe+nB*F_nBnB))/nB;
    double f_nnnp=((1.0-Ye)*Ye*F_YeYe+nB*(2.0*dF_dnB+(1.0-2.0*Ye)*
                                          F_nBYe+nB*F_nBnB))/nB;
    double f_npnp=((1.0-Ye)*(1.0-Ye)*F_YeYe+nB*(2.0*dF_dnB-2.0*(1.0-Ye)*
                                                F_nBYe+nB*F_nBnB))/nB;*/
    //alternate f_nnnn, f_nnnp, f_npnp
    //double f_npnp=F_nBnB-(((1-Ye)/(nB*nB))*dF_dYe)-((((1-Ye)*(1-Ye))/(nB*nB))*F_YeYe);
    //double f_nnnn=F_nBnB-((Ye/(nB*nB))*dF_dYe)-(((Ye*Ye)/(nB*nB))*F_YeYe);
    //double f_nnnp=0.0;
    double f_nnnn=(Ye*Ye*F_YeYe+nB*(2.0*dF_dnB-2.0*Ye*F_nBYe+nB*F_nBnB))/nB;
    double f_nnnp=((Ye-1.0)*Ye*F_YeYe+nB*(2.0*dF_dnB+(1.0-2.0*Ye)*
                                          F_nBYe+nB*F_nBnB))/nB;
    double f_npnp=((Ye-1.0)*(Ye-1.0)*F_YeYe+nB*(2.0*dF_dnB-2.0*(Ye-1.0)*
                                                F_nBYe+nB*F_nBnB))/nB;
    double f_nnT=dF_dT-Ye*F_YeT+nB*F_nBT;
    double f_npT=dF_dT+(Ye-1.0)*F_YeT+nB*F_nBT;
    double f_TT=nB*F_TT; 
    /*
    double f_nnT=dF_dT-Ye*F_YeT+nB*F_nBT;
    double f_npT=dF_dT-(1.0-Ye)*F_YeT+nB*F_nBT;
    double f_TT=nB*F_TT;*/
    
    double den=en*T_MeV/hc_mev_fm+(mun+mneut)*nB*(1.0-Ye)+
      (mup+mprot)*nB*Ye+mue*nB*Ye;
    double nn2=nB*(1.0-Ye);
    double np2=nB*Ye;
    //maybe change back?
    double dmundnB=(f_nnnn*(1-Ye))+(f_nnnp*Ye);
    double dmundYe=nB*(f_nnnp-f_nnnn);
    double dmupmuedYe=nB*(f_npnp-f_nnnp);
    double dsdnB=-f_nnT*(1.0-Ye)-f_npT*Ye;
    double dsdYe=nB*(f_nnT-f_npT);
    double e1=(nB*nB*dmundnB)+((nB*nB*dsdnB*dsdnB)/f_TT)+(nB*Ye*(1-Ye)*dmundYe)+(nB*Ye*Ye*dmupmuedYe);
    double e2=(nB*dsdnB)/f_TT;
    double cs_sq=e1-(2*en*e2)-((en*en)/f_TT);
/*    double cs_sq=(nn2*nn2*(f_nnnn-f_nnT*f_nnT/f_TT)+
                  2.0*nn2*np2*(f_nnnp-f_nnT*f_npT/f_TT)+
                  np2*np2*(f_npnp-f_npT*f_npT/f_TT)-
                  2.0*en*(nn2*f_nnT/f_TT+np2*f_npT/f_TT)-en*en/f_TT)/den;*/

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
      cout << "F_nB: " << dF_dnB << " " << t1 << " " << t2 << endl;
      double tab_dF_dnB = t1;
      cout << F_nBnB << " " << t3 << " ";
     if (!(index[1]==0) && !(index[1]==69)) { 
        t1=(tgp_F->get(index)-tgp_F->get(jm1))/hc_mev_fm/
          (Ye_grid[index[1]]-Ye_grid[index[1]-1]);
        t2=(tgp_F->get(jp1)-tgp_F->get(index))/hc_mev_fm/
          (Ye_grid[index[1]+1]-Ye_grid[index[1]]);
        t3=(t2-t1)*2.0/(Ye_grid[index[1]+1]-Ye_grid[index[1]-1]);
        cout << "F_Ye: " << dF_dYe << " " << t1 << " " << t2 << endl;
        cout << F_YeYe << " " << t3 << " ";

        double tab_mun=tgp_F->get(index)/hc_mev_fm-Ye*t1+nB*tab_dF_dnB;
        double tab_mup=tgp_F->get(index)/hc_mev_fm+(1.0-Ye)*t1+nB*tab_dF_dnB;
        cout << "tab mun, mup: " << tab_mun << " " << tab_mup << endl;
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

      double dmundnB=(f_nnnn*(1.0-Ye))+(f_nnnp*Ye);
      double dmupmuednB=(f_nnnp*(1.0-Ye))+(f_npnp*Ye);
      //double dmundnB=F_nBnB-(Ye*(((1/nB)*F_nBYe)-((1/(nB*nB))*dF_dYe)));
      //double dmupmuednB=F_nBnB-((1-Ye)*(((1/nB)*F_nBYe)-((1/(nB*nB))*dF_dYe)));
      double dPdnB=(dmundnB*nB*(1.0-Ye))+(mun*(1.0-Ye))+(dmupmuednB*Ye*nB)+(Ye*(mup+mue))-((mun*(1.0-Ye))+((mup+mue)*Ye));
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

    cout << "\nF_intp " << Fintp << endl;
    cout << "F_tab " << tgp_F->get(index) << endl;
    double diff = std::abs(tgp_F->get(index)-Fintp)/std::abs(tgp_F->get(index));
    cout << diff << endl;
    if (diff < 0.001) {
        cout << "success\n";
    }
    else {
        cout << "failure\n";
    }
    diff = std::abs(tgp_mun->get(index)-(mun*hc_mev_fm))/std::abs(tgp_mun->get(index));
    cout << diff << endl;
    if (diff < 0.01) {
        cout<<"mun: success\n";
    }
    else {
        cout<<"mun: failure\n";
    }
    double tab_dPdnB=(tgp_P->get(index)-tgp_P->get(im1))/hc_mev_fm/
      (nB_grid[index[0]]-nB_grid[index[0]-1]);
    
    diff = std::abs(tab_dPdnB-dPdnB)/std::abs(tab_dPdnB);
    cout << diff << endl;
    if (dPdnB>0.0 && diff<0.10) {
        cout<<"dPdnB: success\n";
    }
    else {
        cout<<"dPdnB: failure\n";
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
