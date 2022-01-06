/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2022, Xingfu Du, Zidu Lin, and Andrew W. Steiner
  
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
#include <map>

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
      is smaller than \f$ E_d \f$

      From eq. (25) & (27) in Shen 10.
  */
  double delta_small_iand(double E);
    
  /** \brief Integrand for derivative of partition function
      with respect to temperature when \f$ \delta \f$
      is smaller than \f$ E_d \f$

      From eq. (26) & (27) in Shen 10.
  */
  double delta_small_iand_prime(double E);
    
  /** \brief Integrand when \f$ \delta \f$ is greater than \f$ E_d \f$

      From eq. (25) and (30) in Shen 10.
   */
  double delta_large_iand(double E);
    
  /** \brief Integrand for temperature derivative when \f$ \delta
      \f$ is greater than \f$ E_d \f$

      From eq. (26) and (30) in Shen 10.
  */
  double delta_large_iand_prime(double E);
    
};

/** \brief Solve for the EOS including nuclei

    \verbatim embed:rst

    .. todo::

       Class eos_nuclei:

       - Rename n_nB2 to n_nB, etc.
       - Make child of eos_sn_base
       - Use \c loaded instead of testing n_nB2==0
       - Move select_high_T to the parent
       - Future: Allow different form for the NS fit
       - Future: Allow user to specify where data files are
         located or to manually specify the nuclear model?
       - Future: Allow RMF rather than just Skyrme.

    \endverbatim
*/
class eos_nuclei : public eos {

public:

  /// \name Constructor and destructor
  //@{
  eos_nuclei();
  
  virtual ~eos_nuclei();
  //@}

  /// \name Grid specification
  //@{
  size_t n_nB2;
  size_t n_Ye2;
  size_t n_T2;
  std::vector<double> nB_grid2;
  std::vector<double> Ye_grid2;
  std::vector<double> T_grid2;
  std::string nB_grid_spec;
  std::string Ye_grid_spec;
  std::string T_grid_spec;
  //@}

  /// \name Nuclear masses and their fits
  //@{
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
  //@}

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

  /// \name Mathematical algorithm objects
  //@{
  /// Integrator
  o2scl::inte_qag_gsl<> iqg;

  /** \brief Solver 
   */
  o2scl::mroot_hybrids<> mh;

  /** \brief Bracketing solver
   */
  o2scl::root_brent_gsl<> rbg;
  //@}

  /// \name The nuclear partition function
  //@{
  /// Object which integrates partition functions
  partition_func part_func;
  
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
  //@}

  /// \name MPI message values
  //@{
  static const int message_continue=0;
  static const int message_done=1;
  //@}

  /// \name Parameters modifiable by the CLI user
  //@{
  /// The time between file updates in seconds (default 1800)
  double file_update_time;
  /// The number of iterations between file updates (default 1000)
  int file_update_iters;
    
  /** \brief The maximum value of A for a fixed distribution
      (when alg_mode is 4)
  */
  int fd_A_max;

  /// If true, extend FRDM beyond neutron/proton drip (default false)
  bool extend_frdm;
  
  /// The maximum value of \f$ N/Z \f$ or \f$ Z/N \f$ (default 2.25)
  double max_ratio;
  
  /// Fiducial value for solver tolerance (default \f$ 10^{-6} \f$)
  double mh_tol_rel;

  /// File containing external guess
  std::string ext_guess;
  
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

  /** \brief If true, verify points only and do not attempt to solve
      (default false)
  */
  bool verify_only;

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
  int six_neighbors;

  /** \brief A new function verbose parameter
   */
  int function_verbose;
  
  /** \brief If true, use previously computed points (or guesses) as
      an initial guess to compute adjacent points (default true)
  */
  bool propagate_points;

  /** \brief If true survey the nB and Ye equations (default false)
   */
  bool survey_eqs;

  /** \brief If true, output all of the data necessary for a full EOS

      If true, include Eint, Pint, Sint, mun, and mup. If include_eg
      is additionally true, then add F, E, P, and S.
  */
  bool derivs_computed;

  /** \brief Always true, for consistency with o2scl::eos_sn_base
   */
  bool baryons_only_loaded;
  
  /** \brief If true, include electrons and photons

      Requires that derivs_computed is also true
   */
  bool with_leptons_loaded;
  //@}

  /// \name Other internal objects
  //@{
  /** \brief Coulomb energy (in \f$ \mathrm{fm}^{-1} \f$ )
      
      This quantity is computed in \ref solve_nuclei() and then
      used later in \ref eos_fixed_dist().
  */
  ubvector Ec;
  
  /// Dictionary for mapping buffers to physical quantities
  o2scl::vec_index vi;

  /// True if an EOS is currently loaded
  bool loaded;

  /** \brief Ranges for randomly selected ranges in 
      \ref eos_fixed_dist()
   */
  std::vector<double> fd_rand_ranges;

  /** \brief Object for sending slack messages
   */
  o2scl::slack_messenger slack;

  /** \brief List of electron fractions to examine for computing
   */
  std::string Ye_list;
  //@}

  /// \name Other internal physics objects
  //@{
  /// Virial solver used in \ref check_virial()
  virial_solver_deriv vsd;

  /// Extended Skyrme model for finite-temperature corrections
  eos_had_skyrme_ext skyrme_ext;
  //@}
  
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
  
  /// \name Main EOS table storage
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
  o2scl::tensor_grid3<> tg3_mue;
  o2scl::tensor_grid3<> tg3_Z;
  o2scl::tensor_grid3<> tg3_A;
  o2scl::tensor_grid3<> tg3_Xn;
  o2scl::tensor_grid3<> tg3_Xp;
  o2scl::tensor_grid3<> tg3_Xalpha;
  o2scl::tensor_grid3<> tg3_Xnuclei;
  o2scl::tensor_grid3<> tg3_Ymu;
  o2scl::tensor_grid3<> tg3_Xd;
  o2scl::tensor_grid3<> tg3_Xt;
  o2scl::tensor_grid3<> tg3_XHe3;
  o2scl::tensor_grid3<> tg3_XLi4;
  o2scl::tensor_grid3<> tg3_A_min;
  o2scl::tensor_grid3<> tg3_A_max;
  o2scl::tensor_grid3<> tg3_NmZ_min;
  o2scl::tensor_grid3<> tg3_NmZ_max;
  //@}

  /// \name Detail storage
  //@{
  bool include_detail;
  o2scl::tensor_grid3<> tg3_zn;
  o2scl::tensor_grid3<> tg3_zp;
  o2scl::tensor_grid3<> tg3_F1;
  o2scl::tensor_grid3<> tg3_F2;
  o2scl::tensor_grid3<> tg3_F3;
  o2scl::tensor_grid3<> tg3_F4;
  o2scl::tensor_grid3<> tg3_Un;
  o2scl::tensor_grid3<> tg3_Up;
  o2scl::tensor_grid3<> tg3_msn;
  o2scl::tensor_grid3<> tg3_msp;
  o2scl::tensor_grid3<> tg3_g;
  o2scl::tensor_grid3<> tg3_dgdT;
  //@}
  
  /// \name Other parameter objects
  //@{
  o2scl::cli::parameter_bool p_survey_eqs;
  o2scl::cli::parameter_bool p_extend_frdm;
  o2scl::cli::parameter_bool p_show_all_nuclei;
  o2scl::cli::parameter_int p_fd_A_max;
  o2scl::cli::parameter_bool p_recompute;
  o2scl::cli::parameter_bool p_verify_only;
  o2scl::cli::parameter_string p_edge_list;
  o2scl::cli::parameter_string p_ext_guess;
  o2scl::cli::parameter_int p_six_neighbors;
  o2scl::cli::parameter_bool p_full_results;
  o2scl::cli::parameter_bool p_rnuc_less_rws;
  o2scl::cli::parameter_bool p_include_eg;
  o2scl::cli::parameter_bool p_propagate_points;
  o2scl::cli::parameter_bool p_include_detail;
  o2scl::cli::parameter_double p_mh_tol_rel;
  o2scl::cli::parameter_double p_max_time;
  o2scl::cli::parameter_string p_nucleon_func;
  o2scl::cli::parameter_int p_alg_mode;
  o2scl::cli::parameter_int p_fixed_dist_alg;
  o2scl::cli::parameter_int p_function_verbose;
  o2scl::cli::parameter_string p_Ye_list;
  o2scl::cli::parameter_double p_max_ratio;
  o2scl::cli::parameter_double p_file_update_time;
  o2scl::cli::parameter_int p_file_update_iters;
  o2scl::cli::parameter_string p_nB_grid_spec;
  o2scl::cli::parameter_string p_Ye_grid_spec;
  o2scl::cli::parameter_string p_T_grid_spec;
  //@}

  /// \name Functions for the main algorithm
  //@{
  /** \brief Use only one of the two equations for a 
      function for a root bracketing algorithm
  */
  double solve_nuclei_ld(double x2, size_t nv, const ubvector &x, 
			 double nb, double ye, double T,
			 int ix, double &mun_gas, double &mup_gas,
			 o2scl::thermo &th_gas);

  /** \brief Solve for the nuclei via minimization
   */
  double solve_nuclei_min(size_t nv, const ubvector &x, 
			  double nb, double ye, double T,
			  double &mun_gas, double &mup_gas,
			  o2scl::thermo &th_gas);

  /** \brief Construct equations to solve for a fixed baryon
      density and electron fraction
  */
  int solve_nuclei(size_t nv, const ubvector &x, ubvector &y, double nb,
		   double ye, double T, 
		   int loc_verbose, double &mun_gas, double &mup_gas,
		   o2scl::thermo &th_gas, std::map<std::string,double> &vdet);
  
  /** \brief Determine the EOS presuming a fixed single heavy nucleus
      and solving for the log (base 10) of the
      free neutron and proton abundances
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
		  double flag, std::map<std::string,double> &vdet);

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
   int &A_max, int &NmZ_min, int &NmZ_max,
   std::map<std::string,double> &vdet,
   bool dist_changed, bool no_nuclei);

  /** \brief Determine the EOS presuming a distribution of nuclei
      with fixed limits in A and \f$ N-Z \f$ but used to fix table only
  */
  int eos_fixed_dist_fix_table
  (double nB, double Ye, double T, double &log_xn, double &log_xp,
   o2scl::thermo &thx, double &mun_full, double &mup_full, int &A_min,
   int &A_max, int &NmZ_min, int &NmZ_max, bool dist_changed,
   bool no_nuclei);

  /** \brief Compute the EOS presuming homogeneous nuclear matter
   */
  int nuc_matter(double nB, double Ye, double T,
		 double &log_xn, double &log_xp,
		 double &Zbar, double &Nbar, o2scl::thermo &thx,
		 double &mun_full, double &mup_full,
		 int &A_min, int &A_max,
		 int &NmZ_min, int &NmZ_max,
		 std::map<std::string,double> &vdet);

  /** \brief Determine the EOS presuming a distribution of nuclei
      and optimizing the limits in A and \f$ N-Z \f$
  */
  int eos_vary_dist(double nB, double Ye, double T, double &log_xn,
		    double &log_xp, double &Zbar, double &Nbar, 
		    o2scl::thermo &thx, double &mun_full, double &mup_full,
		    int &A_min, int &A_max, int &NmZ_min, int &NmZ_max,
		    std::map<std::string,double> &vdet, bool dist_changed,
		    bool no_nuclei);
   
  /** \brief Generate a table (MPI version)
   */
  int generate_table(std::vector<std::string> &sv, bool itive_com);
  //@}

  /// \name EOS post-processing functions
  //@{
  /** \brief Compute derivatives numerically
   */
  int eos_deriv(std::vector<std::string> &sv, bool itive_com);
  int eos_deriv_v2(std::vector<std::string> &sv, bool itive_com);

  /** \brief Compute second derivatives numerically

      The derivatives 

   */
  int eos_second_deriv(std::vector<std::string> &sv, bool itive_com);

  /** \brief Add electrons and photons
   */
  int add_eg(std::vector<std::string> &sv, bool itive_com);

  /** \brief Construct an electrons and photon table
   */
  int eg_table(std::vector<std::string> &sv, bool itive_com);

  /** \brief Edit an EOS table
   */
  int edit_data(std::vector<std::string> &sv, bool itive_com);

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
  
  int test_random(std::vector<std::string> &sv, bool itive_com);
  //@}

  /// \name File I/O functions
  //@{
  /** \brief Load an EOS table (CLI wrapper)
   */
  int load(std::vector<std::string> &sv, bool itive_com);

  /** \brief Output an EOS table to a file (CLI wrapper)
   */
  int output(std::vector<std::string> &sv, bool itive_com);

  /** \brief Write results to an HDF5 file

      \future Eventually replace this with eos_sn_base::output()?
  */
  int write_results(std::string fname);
  
  /** \brief Read results from an HDF5 file

      \future Eventually replace this with eos_sn_base::load()?
  */
  int read_results(std::string fname);
  
  /** \brief Write the nuclear masses to an HDF5 file
   */
  int write_nuclei(std::vector<std::string> &sv,
			       bool itive_com);

  /** \brief Load nuclear masses

      This function is called in <tt>main()</tt>.
   */
  void load_nuclei();
  
  /** \brief Write nuclear masses to a file
   */
  void write_nuclei(std::string fname);
  //@}
  
  /// \name Miscellaneous functions
  //@{
  /** \brief Setup the command-line interface
   */
  virtual void setup_cli(o2scl::cli &cli); 

  /** \brief Initialize tensors for a new EOS table
   */
  void new_table();

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
  
  /** \brief Use results lower densities to provide initial 
      guesses to higher densities
   */
  int increase_density(std::vector<std::string> &sv, bool itive_com);

  /** \brief Fix the core-crust transition
   */
  int fix_cc(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Verify the EOS
   */
  int verify(std::vector<std::string> &sv, bool itive_com);

  /** \brief Monte Carlo results with nuclei
   */
  int mcarlo_nuclei(std::vector<std::string> &sv, bool itive_com);

  /** \brief Monte Carlo results with nuclei (version 2)
   */
  int mcarlo_nuclei2(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Monte Carlo results in beta equilibrium
   */
  int mcarlo_beta(std::vector<std::string> &sv, bool itive_com);

  /** \brief Test the neutrino opacities
   */
  int test_neutrino(std::vector<std::string> &sv, bool itive_com);
  
  /// Compute the baryon number fractions and put them in \c X
  void compute_X(double nB, ubvector &X);
  
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

  /** \brief Compute eos with nuclei by searching minimum
   */
  double f_min_search(size_t nvar,const ubvector &x,
		      double nb, double ye, double T);
  
  /** \brief Initialization for differential evolution approach
   */		      
  int init_function(size_t dim, const ubvector &x, ubvector &y);
  
  /** \brief Old Maxwell construction test
   */
  int maxwell_test(std::vector<std::string> &sv, bool itive_com);
  //@}
  
};

