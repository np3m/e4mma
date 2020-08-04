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
#include "eos.h"
#include "eos_had_skyrme_ext.h"
#include "virial_solver_deriv.h"
#include <o2scl/nucmass_fit.h>
#include <o2scl/slack_messenger.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

/** \brief Compute partition functions using Fowler prescription
 */
class partition_func {
  
public:

  /// Temperature in MeV
  double T_MeV;
    
  /// Nuclear level spacing in 1/MeV
  double a;
    
  /// Backshift parameter in MeV
  double delta;
    
  /// Coefficient for large \f$ \delta \f$ in 1/MeV
  double C;
    
  /// Critical temperature connecting low and high energies in MeV
  double Tc;
    
  /** \brief Integrand for partition function when \f$ \delta \f$
      is small

      From eq. (25) & (27) in Shen 10.
  */
  double delta_small_iand(double x);
    
  /** \brief Integrand for derivative of partition function
      with respect to temperature when \f$ \delta \f$
      is small

      From eq. (26) & (27) in Shen 10.
  */
  double delta_small_iand_prime(double x);
    
  /** \brief Integrand when \f$ \delta \f$ is large

      From eq. (25) and (30) in Shen 10.
   */
  double delta_large_iand(double x);
    
  /** \brief Integrand for temperature derivative when \f$ \delta
      \f$ is large

      From eq. (26) and (30) in Shen 10.
  */
  double delta_large_iand_prime(double x);
    
};

/** \brief Solve for the EOS including nuclei

    \todo Rename n_nB2 to n_nB, etc.
    \todo Make child of eos_sn_base
*/
class eos_nuclei : public eos {

public:

  eos_nuclei();

  virtual ~eos_nuclei();
  
  size_t n_nB2;
  size_t n_Ye2;
  size_t n_T2;
  std::vector<double> nB_grid2;
  std::vector<double> Ye_grid2;
  std::vector<double> T_grid2;
  bool loaded;

  /** \brief Object for sending slack messages
   */
  o2scl::slack_messenger slack;

  /// Virial solver used in \ref check_virial()
  virial_solver_deriv vsd;
  
  /// Integrator
  o2scl::inte_qag_gsl<> iqg;

  /** \brief List of electron fractions to examine for computing
   */
  std::string Ye_list;
  
  /// Extended Skyrme model for finite-temperature corrections
  o2scl::eos_had_skyrme_ext skyrme_ext;
  
  /** \brief Select high-temperature EOS
      
      \note This function is called by the constructor and thus
      cannot be virtual
  */
  int select_high_T(int option);

  /** \brief Compute the second derivatives and the
      eigenvalues of the stability matrix
  */
  int stability(std::vector<std::string> &sv,
		bool itive_com);
  
  /** \brief Command-line interface for selection of the high-temperature EOS
   */
  int select_high_T_cl(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Fit nuclear masses
   */
  int fit_frdm(std::vector<std::string> &sv,
	       bool itive_com);
  
  /** \brief Theoretical nuclear masses
   */
  o2scl::nucmass_mnmsk m95;

  /** \brief HFB masses for spin predictions
   */
  o2scl::nucmass_hfb_sp hfb;

  /** \brief Experimental nuclear masses
   */
  o2scl::nucmass_ame ame;

  /** \brief Theoretical nuclear masses far from stability
   */
  o2scl::nucmass_frdm frdm;

  /** \brief Fit theory masses
   */
  o2scl::nucmass_fit nm_fit;

  /// Desc (default 2.25)
  double max_ratio;
  
  /** \brief The number of neutrons in the nucleus
   */
  size_t nucN;

  /** \brief The number of protons in the nucleus
   */
  size_t nucZ;

  /// \name Nucleus objects
  //@{
  std::vector<o2scl::nucleus> nuclei;
  o2scl::nucleus *nuc_alpha;
  o2scl::nucleus *nuc_deut;
  o2scl::nucleus *nuc_trit;
  o2scl::nucleus *nuc_he3;
  o2scl::nucleus *nuc_li4;
  o2scl::nucleus *nuc_heavy;
  //@}

  /** \brief Solver 
   */
  o2scl::mroot_hybrids<> mh;

  /** \brief Bracketing solver
   */
  o2scl::root_brent_gsl<> rbg;
  
  /** \brief Neutron separation energies (in MeV)
   */
  ubvector Sneut;
  
  /** \brief Proton separation energies (in MeV)
   */
  ubvector Sprot;

  /** \brief Partition function
   */
  ubvector vomega;
  
  /** \brief Derivative of partition function with respect to 
      temperature
  */
  ubvector vomega_prime;
  
  /** \brief Coulomb energy (in \f$ \mathrm{fm}^{-1} \f$ )
   */
  ubvector Ec;
  
  /** \ Store the results of solve_nuclei
   */
  double mun_old, mup_old, mun, mup, nn, np;
  double p_c, p_c_test;
  double f0_nuc, p0_nuc, ed0_nuc, en0_nuc;
  double f_total, Pr_total, en_total, ed_total;
  double flag;
  bool heavy;

  double file_update_time;
  int file_update_iters;
  o2scl::cli::parameter_double p_file_update_time;
  o2scl::cli::parameter_int p_file_update_iters;
  
  /// \name MPI message values
  //@{
  static const int message_continue=0;
  static const int message_done=1;
  //@}

  /// \name Grid specification
  //@{
  std::string nB_grid_spec;
  std::string Ye_grid_spec;
  std::string T_grid_spec;
  o2scl::cli::parameter_string p_nB_grid_spec;
  o2scl::cli::parameter_string p_Ye_grid_spec;
  o2scl::cli::parameter_string p_T_grid_spec;
  //@}

  /// Dictionary for mapping buffers to physical quantities
  o2scl::vec_index vi;

  /** \brief Function to describe the Z and N range 

      This is the function used to specify the range of Z and N which
      is varied to find the smallest free energy
  */
  std::string nucleon_func;

  /// Maximum MPI time (default 0)
  double max_time;
  
  /** \brief Show all nuclei considered at every point (default false)
   */
  bool show_all_nuclei;

  /** \brief If true, ensure that the nuclear radius is less than the
      Wigner-Seitz radius (default true)
  */
  bool rnuc_less_rws;
  
  /** \brief If true, recompute all points, irrespective of the
      value of the convergence flag (default false)
  */
  bool recompute;

  /** \brief If true, recompute points where Z, A, log_xn, 
      or log_xp reach a maximum or minimum (default false)
  */
  std::string edge_list;

  /** \brief Algorithm mode (default 1)

      0 for AWS SNA, 1 for XD SNA, 2 for AWS dist, 3 for XD dist
  */
  int alg_mode;

  /** \brief Algorithm for \ref eos_fixed_dist()
   */
  int fixed_dist_alg;
  
  /** \brief If true, when computing a point, perform the calculation
      also using the six neighboring points as an initial guess
      (default false)
  */
  bool six_neighbors;

  /** \brief A new function verbose parameter
   */
  int function_verbose;
  
  /** \brief If true, use previously computed points (or guesses) as
      an initial guess to compute adjacent points (default true)
  */
  bool propagate_points;

  /** \brief If true, output all of the data necessary for a full EOS

      If true, include Eint, Pint, Sint, mun, and mup. If include_eg
      is additionally true, then add F, E, P, and S.
  */
  bool full_results;

  /** \brief If true, include electrons and photons
   */
  bool include_eg;

  /// \name Other parameter objects
  //@{
  o2scl::cli::parameter_bool p_show_all_nuclei;
  o2scl::cli::parameter_bool p_recompute;
  o2scl::cli::parameter_string p_edge_list;
  o2scl::cli::parameter_bool p_six_neighbors;
  o2scl::cli::parameter_bool p_full_results;
  o2scl::cli::parameter_bool p_rnuc_less_rws;
  o2scl::cli::parameter_bool p_include_eg;
  o2scl::cli::parameter_bool p_propagate_points;
  o2scl::cli::parameter_double p_mh_tol_rel;
  o2scl::cli::parameter_double p_max_time;
  o2scl::cli::parameter_string p_nucleon_func;
  o2scl::cli::parameter_int p_alg_mode;
  o2scl::cli::parameter_int p_fixed_dist_alg;
  o2scl::cli::parameter_int p_function_verbose;
  o2scl::cli::parameter_string p_Ye_list;
  o2scl::cli::parameter_double p_max_ratio;
  //@}

  /** \brief Compute derivatives analytically
   */
  int eos_deriv(std::vector<std::string> &sv, bool itive_com);

  /** \brief Desc
   */
  int add_eg(std::vector<std::string> &sv, bool itive_com);

  /** \brief Desc
   */
  int maxwell_test(std::vector<std::string> &sv, bool itive_com);

  /** \brief Check the virial solver by using it to compute
      the EOS over a wide range of densities and temperatures
      
      This function is particularly good for checking to make
      sure that the virial part of the EOS does not lead to 
      any discontinuities. 

      For example, 

      o2graph -read check_virial.o2 zn -Ye-slice 0.3 -set logz 1 \
      -den-plot slice -show
  */
  int check_virial(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Construct an equation to solve for matter at low 
      densities
  */
  double solve_nuclei_ld(double x2, size_t nv, const ubvector &x, 
			 double nb, double ye, double T,
			 int ix, double &mun_gas, double &mup_gas,
			 o2scl::thermo &th_gas);

  /** \brief Desc
   */
  double solve_nuclei_min(size_t nv, const ubvector &x, 
			  double nb, double ye, double T,
			  double &mun_gas, double &mup_gas,
			  o2scl::thermo &th_gas);

  /** \brief Write results to an HDF5 file

      \todo Eventually replace this with eos_sn_base::output()
  */
  int write_results(std::string fname);
  
  /** \brief Read results from an HDF5 file

      \todo Eventually replace this with eos_sn_base::load()
  */
  int read_results(std::string fname);
  
  /** \brief Construct equations to solve for a fixed baryon
      density and electron fraction (AWS version)
  */
  int solve_nuclei(size_t nv, const ubvector &x, ubvector &y, double nb,
		   double ye, double T, 
		   int loc_verbose, double &mun_gas, double &mup_gas,
		   o2scl::thermo &th_gas);
  
  /** \brief Determine the EOS presuming a fixed single heavy nucleus
      and solving for the log (base 10) of the
      free neutron and proton abundances (AWS version)
  */
  int eos_fixed_ZN(double nb, double ye, double T,
		   double &log_xn, double &log_xp,
		   size_t nuc_Z1, size_t nuc_N1,
		   o2scl::thermo &thx,
		   double &mun_full, double &mup_full);

  /** \brief Store data in the tensor objects
   */
  int store_point(size_t i_nB, size_t i_Ye, size_t i_T,
		  double nB, double Ye, double T,
		  o2scl::thermo &th, double log_xn, double log_xp,
		  double Zbar, double Nbar, 
		  double mun_full, double mup_full, ubvector &X,
		  double A_min, double A_max, double NmZ_min, double NmZ_max,
		  double flag=10.0);


  /** \brief Determine the EOS allowing the Z and N of the nucleus
      to vary
  */
  int eos_vary_ZN(double nb, double ye, double T,
		  double &log_xn, double &log_xp,
		  size_t &nuc_Z1, size_t &nuc_N1,
		  o2scl::thermo &thx,
		  double &mun_full, double &mup_full,
		  bool nu_nuclei=false);

  /** \brief Determine the EOS presuming a distribution of nuclei
      with fixed limits in A and \f$ N-Z \f$
  */
  int eos_fixed_dist
  (double nB, double Ye, double T, double &log_xn, double &log_xp,
   o2scl::thermo &thx, double &mun_full, double &mup_full, int &A_min,
   int &A_max, int &NmZ_min, int &NmZ_max, bool dist_changed,
   bool no_nuclei);


  /** \brief Determine the EOS presuming a distribution of nuclei
      with fixed limits in A and \f$ N-Z \f$ but used to fix table only
  */
  int eos_fixed_dist_fix_table
  (double nB, double Ye, double T, double &log_xn, double &log_xp,
   o2scl::thermo &thx, double &mun_full, double &mup_full, int &A_min,
   int &A_max, int &NmZ_min, int &NmZ_max, bool dist_changed,
   bool no_nuclei);
  
  /** \brief Determine the EOS presuming a distribution of nuclei
      and optimizing the limits in A and \f$ N-Z \f$
  */
  int eos_vary_dist
  (double nB, double Ye, double T, double &log_xn, double &log_xp,
   double &Zbar, double &Nbar, 
   o2scl::thermo &thx, double &mun_full, double &mup_full, int &A_min,
   int &A_max, int &NmZ_min, int &NmZ_max, bool dist_changed,
   bool no_nuclei);
   
  /** \brief Generate a table (MPI version)
   */
  int generate_table(std::vector<std::string> &sv, bool itive_com);

  /** \brief Desc
   */
  int load(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Desc
   */
  int output(std::vector<std::string> &sv, bool itive_com);

  /** \brief Edit an EOS table
   */
  int edit_data(std::vector<std::string> &sv, bool itive_com);

  /// \name Flag values
  //@{
  /// Point is empty
  static const int iflag_empty=0;
  /// Point is in progress, and empty
  static const int iflag_in_progress_empty=-1;
  /// Point is in progress, and has initial guess 
  static const int iflag_in_progress_with_guess=-5;
  /// Point has initial guess
  static const int iflag_guess=5;
  /// Point is finished
  static const int iflag_done=10;
  //@}
  
  /** \brief Merge two tables
   */
  int merge_tables(std::vector<std::string> &sv, bool itive_com);

  /** \brief Compare two tables
   */
  int compare_tables(std::vector<std::string> &sv, bool itive_com);

  /** \brief Output the statistics on flag values for a table
   */
  int stats(std::vector<std::string> &sv, bool itive_com);

  /** \brief Compute the EOS at one point
   */
  int point_nuclei(std::vector<std::string> &sv, bool itive_com);

  /** \brief Desc
   */
  int mcarlo_nuclei(std::vector<std::string> &sv, bool itive_com);

  /** \brief Setup the command-line interface
   */
  virtual void setup_cli(o2scl::cli &cli); 

  /** \brief Initialize tensors for a new EOS table
   */
  void new_table();
  
  /// \name Tensors for full output
  //@{
  o2scl::tensor_grid3<> tg3_log_xn;
  o2scl::tensor_grid3<> tg3_log_xp;
  o2scl::tensor_grid3<> tg3_flag;
  o2scl::tensor_grid3<> tg3_F;
  o2scl::tensor_grid3<> tg3_E;
  o2scl::tensor_grid3<> tg3_P;
  o2scl::tensor_grid3<> tg3_S;
  o2scl::tensor_grid3<> tg3_Fint;
  o2scl::tensor_grid3<> tg3_Eint;
  o2scl::tensor_grid3<> tg3_Pint;
  o2scl::tensor_grid3<> tg3_Sint;
  o2scl::tensor_grid3<> tg3_mun;
  o2scl::tensor_grid3<> tg3_mup;
  o2scl::tensor_grid3<> tg3_Z;
  o2scl::tensor_grid3<> tg3_A;
  o2scl::tensor_grid3<> tg3_Xn;
  o2scl::tensor_grid3<> tg3_Xp;
  o2scl::tensor_grid3<> tg3_Xalpha;
  o2scl::tensor_grid3<> tg3_Xnuclei;
  o2scl::tensor_grid3<> tg3_Xd;
  o2scl::tensor_grid3<> tg3_Xt;
  o2scl::tensor_grid3<> tg3_XHe3;
  o2scl::tensor_grid3<> tg3_XLi4;
  o2scl::tensor_grid3<> tg3_A_min;
  o2scl::tensor_grid3<> tg3_A_max;
  o2scl::tensor_grid3<> tg3_NmZ_min;
  o2scl::tensor_grid3<> tg3_NmZ_max;
  //@}

  /// Object which integrates partition functions
  partition_func part_func;
  
  /** \brief Compute eos with nuclei by searching minimum
   */
  double f_min_search(size_t nvar,const ubvector &x,
		      double nb, double ye, double T);
  
  /** \brief initialization for differential evolution approach
   */		      
  int init_function(size_t dim, const ubvector &x, ubvector &y);
  
};

