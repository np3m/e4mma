/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018, Xingfu Du and Andrew W. Steiner
  
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

using namespace std;
using namespace o2scl;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  eos eph;
  cli cl;
  
  static const int nopt=11;
  o2scl::comm_option_s options[nopt]={
    {0,"test_deriv","Desc.",0,0,"","",
     new o2scl::comm_option_mfptr<eos>
     (&eph,&eos::test_deriv),o2scl::cli::comm_option_both},
    {0,"table_Ye","Desc.",2,2,"<fname> <Ye>","",
     new o2scl::comm_option_mfptr<eos>
     (&eph,&eos::table_Ye),o2scl::cli::comm_option_both},
    {0,"table_full","Desc.",1,1,"<fname>","",
     new o2scl::comm_option_mfptr<eos>
     (&eph,&eos::table_full),o2scl::cli::comm_option_both},
    {0,"vir_fit","Desc.",0,0,"","",
     new o2scl::comm_option_mfptr<eos>
     (&eph,&eos::vir_fit),o2scl::cli::comm_option_both},
    {0,"eos_sn","Desc.",0,0,"","",
     new o2scl::comm_option_mfptr<eos>
     (&eph,&eos::eos_sn),o2scl::cli::comm_option_both},
    {0,"mcarlo_data","Desc.",0,1,"Monte Carlo function","",
     new o2scl::comm_option_mfptr<eos>
     (&eph,&eos::mcarlo_data),o2scl::cli::comm_option_both},
    {0,"point","Desc.",0,3,"","",
     new o2scl::comm_option_mfptr<eos>
     (&eph,&eos::point),o2scl::cli::comm_option_both},
    {0,"random","Desc.",0,0,"","",
     new o2scl::comm_option_mfptr<eos>
     (&eph,&eos::random),o2scl::cli::comm_option_both},
    {0,"select_model","Desc.",7,7,"","",
     new o2scl::comm_option_mfptr<eos>
     (&eph,&eos::select_model),o2scl::cli::comm_option_both},
    {0,"teg","Desc.",0,0,"","",
     new o2scl::comm_option_mfptr<eos>
     (&eph,&eos::test_eg),o2scl::cli::comm_option_both},
    {0,"vir_comp","Desc.",0,0,"","",
     new o2scl::comm_option_mfptr<eos>
     (&eph,&eos::vir_comp),o2scl::cli::comm_option_both}
  };
  cl.set_comm_option_vec(nopt,options);
  cl.gnu_intro=false;

  o2scl::cli::parameter_int p_verbose;
  p_verbose.i=&eph.verbose;
  p_verbose.help="Verbose parameter (default 1)";
  cl.par_list.insert(make_pair("verbose",&p_verbose));

  o2scl::cli::parameter_bool p_old_ns_fit;
  p_old_ns_fit.b=&eph.old_ns_fit;
  p_old_ns_fit.help="Old NS fit (default 0)";
  cl.par_list.insert(make_pair("old_ns_fit",&p_old_ns_fit));

  o2scl::cli::parameter_bool p_ns_record;
  p_ns_record.b=&eph.ns_record;
  p_ns_record.help="Record NS fit (default 0)";
  cl.par_list.insert(make_pair("ns_record",&p_ns_record));

  o2scl::cli::parameter_bool p_include_muons;
  p_include_muons.b=&eph.include_muons;
  p_include_muons.help="If true, include muons (default true)";
  cl.par_list.insert(make_pair("include_muons",&p_include_muons));

  o2scl::cli::parameter_bool p_select_cs2_test;
  p_select_cs2_test.b=&eph.select_cs2_test;
  p_select_cs2_test.help="Test cs2 in select_internal() (default 1)";
  cl.par_list.insert(make_pair("select_cs2_test",&p_select_cs2_test));

  o2scl::cli::parameter_bool p_test_ns_cs2;
  p_test_ns_cs2.b=&eph.test_ns_cs2;
  p_test_ns_cs2.help="Test neutron star cs2 (default 0)";
  cl.par_list.insert(make_pair("test_ns_cs2",&p_test_ns_cs2));

  o2scl::cli::parameter_double p_a_virial;
  p_a_virial.d=&eph.a_virial;
  p_a_virial.help="Virial coefficient a (default 3.0)";
  cl.par_list.insert(make_pair("a_virial",&p_a_virial));

  o2scl::cli::parameter_double p_b_virial;
  p_b_virial.d=&eph.b_virial;
  p_b_virial.help="Virial coefficient b (default 0.0)";
  cl.par_list.insert(make_pair("b_virial",&p_b_virial));

  cl.run_auto(argc,argv);
  
  return 0;
}
