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
#include <fstream>

#ifndef NO_MPI
#include <mpi.h>
#endif

#include <o2scl/test_mgr.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/cli.h>
#include <o2scl/lib_settings.h>
#include <o2scl/eos_crust_virial.h>
#include <o2scl/eos_sn.h>
#include <o2scl/rng.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/eos_had_rmf_hyp.h>
#include <o2scl/eos_had_virial.h>
#include <o2scl/part_pdg.h>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

/** \brief An updated version of \ref o2scl::eos_crust_virial
    with a better fit for the virial coefficients
*/
class eos_crust_virial_v2 : public o2scl::eos_crust_virial {
  
 public:

  /** \brief If true, include the deuteron contribution in the
      virial coefficients
   */
  bool include_deuteron;
  
  /** \brief Temperature grid for alpha-nucleon virial coefficients
   */
  std::vector<double> ba_T;
  
  /** \brief Alpha-nucleon virial coefficient
   */
  std::vector<double> vba;
  
  /** \brief Isospin alpha-nucleon virial coefficient
   */
  std::vector<double> vbpna;
  
  /** \brief Interpolator for alpha-nucleon virial coefficient
   */
  o2scl::interp_vec<std::vector<double> > iv_ba;
  
  /** \brief Interpolator for isospin alpha-nucleon virial coefficient
   */
  o2scl::interp_vec<std::vector<double> > iv_bpna;
  
  /** \brief Free nucleon virial coefficient
   */
  double bn0_free();

  /** \brief Free isospin nucleon virial coefficient
   */
  double bpn0_free();

  /** \brief Free spin nucleon virial coefficient
   */
  double bn1_free();

  /** \brief Free spin-isospin nucleon virial coefficient
   */
  double bpn1_free();

  /** \brief Nucleon virial coefficient
   */
  double bn0(double T) {
    return (bn_f(T)-bna(T))/2.0;
  }
  
  /** \brief Spin nucleon virial coefficient
   */
  double bn1(double T) {
    return (bn_f(T)+bna(T))/2.0;
  }
  
  /** \brief Isopin nucleon virial coefficient
   */
  double bpn0(double T) {
    return (bpn_f(T)-bpna(T))/2.0;
  }
  
  /** \brief Spin-isospin nucleon virial coefficient
   */
  double bpn1(double T) {
    return (bpn_f(T)+bpna(T))/2.0;
  }
  
  /** \brief Nucleon-alpha virial coefficient

      The temperature must be specified in MeV
   */
  double bna(double T);

  /** \brief Nucleon-alpha isospin virial coefficient

      The temperature must be specified in MeV
   */
  double bpna(double T);

  /** \brief Fermi-liquid parameter \f$ F_0 \f$

      The value of lambda should be in 1/MeV and the 
      temperature should be specified in MeV. The result
      is returned in units of 1/MeV^2.
   */
  double f0(double lambda, double T);

  /** \brief Fermi-liquid parameter \f$ F_0^{\prime} \f$

      The value of lambda should be in 1/MeV and the 
      temperature should be specified in MeV. The result
      is returned in units of 1/MeV^2.
   */
  double f0p(double lambda, double T);

  /** \brief Fermi-liquid parameter \f$ G_0 \f$

      The value of lambda should be in 1/MeV and the 
      temperature should be specified in MeV. The result
      is returned in units of 1/MeV^2.
   */
  double g0(double lambda, double T);
  
  /** \brief Fermi-liquid parameter \f$ G_0^{\prime} \f$
      
      The value of lambda should be in 1/MeV and the 
      temperature should be specified in MeV. The result
      is returned in units of 1/MeV^2.
  */
  double g0p(double lambda, double T);

  /** \brief The neutron-neutron virial coefficient given the
      function parameters specified in \c par
   */
  double bn_func(size_t np, const std::vector<double> &par, double T);
  
  /** \brief The neutron-proton virial coefficient given the
      function parameters specified in \c par
   */
  double bpn_func(size_t np, const std::vector<double> &par, double T);

  /** \brief The neutron-neutron virial coefficient

      The temperature should be specified in MeV.
   */
  double bn_f(double T);
  
  /** \brief The neutron-proton virial coefficient

      The temperature should be specified in MeV.
   */
  double bpn_f(double T);
  
  /** \brief The temperature derivative of the
      neutron-neutron virial coefficient

      The temperature should be specified in MeV.
   */
  double dbndT_f(double T);
  
  /** \brief The temperature derivative of the
      neutron-proton virial coefficient
      
      The temperature should be specified in MeV.
  */
  double dbpndT_f(double T);
  
  /** \brief The current neutron-neutron virial coefficient parameters
   */
  std::vector<double> bn_params;
  
  /** \brief The current neutron-proton virial coefficient parameters
   */
  std::vector<double> bpn_params;

  /** \brief The number of neutron-neutron virial coefficient parameters
   */
  static const size_t bn_np=10;

  /** \brief The number of neutron-proton virial coefficient parameters
   */
  static const size_t bpn_np=6;
  
  /** \brief Perform the fit to the scattering data
   */
  virtual void fit(bool show_fit=false);

  /// \name Constructor and destructor
  //@{
  eos_crust_virial_v2();

  virtual ~eos_crust_virial_v2() {}
  //@}
  
}; 

/** \brief Phenomenological EOS for homogeneous nucleonic matter
 */
class eos {
  
public:

  /** \brief If true, increase the verbosity for cs2 (default false)
   */
  int cs2_verbose;
  
  /** \brief Particle database
   */
  std::vector<o2scl::part_pdg_db::pdg_entry> part_db;

  /** \brief Fermionic resonances
   */
  std::vector<o2scl::fermion> res_f;
  
  /** \brief Bosonic resonances
   */
  std::vector<o2scl::boson> res_b;
  
 protected:

  /// Desc
  void read_data_files();  

  /** \brief Process the grid specification strings and store
      the values in the arrays
   */
  int process_grid_spec();

public:
  
  /// \name Grid specification
  //@{
  size_t n_nB2;
  size_t n_Ye2;
  size_t n_T2;
  size_t n_S2;
  std::vector<double> nB_grid2;
  std::vector<double> Ye_grid2;
  std::vector<double> T_grid2;
  std::vector<double> S_grid2;
  
  /** \brief The function for default baryon density grid. 
      
      This parameter is used by the new_table() function, and the
      \c check-virial and \c eos-deriv commands.
  */
  std::string nB_grid_spec;
  
  /** \brief The function for default electron fraction grid. 
      
      This parameter is used by the new_table() function, and the
      \c check-virial and eos-deriv \c commands.
  */
  std::string Ye_grid_spec;
  
  /** \brief The function for default temperature grid. 
      
      This parameter is used by the new_table() function, and the
      \c check-virial and \c eos-deriv commands.
  */
  std::string T_grid_spec;
  
  /** \brief The function for default strangeness grid

      This parameter is used by the new_table() function, and the
      \c check-virial and \c eos-deriv commands.
   */
  std::string S_grid_spec;
  //@}

  /// \name Main EOS parameters [protected]
  //@{
  /// The first exponent for density in the QMC EOS (unitless)
  double qmc_alpha;

  /// The first coefficient for the QMC EOS (in MeV)
  double qmc_a;
  
  /** \brief The speed of sound in neutron star matter at 
      \f$ n_B=2~\mathrm{fm}^{-3} \f$
   */
  double phi;

  /// The symmetry energy
  double eos_S;

  /// The slope of the symmetry energy
  double eos_L;

  /// The index of the neutron star model
  int i_ns;

  /// The index of the Skyrme model
  int i_skyrme;
  //@}

  /// \name Basic EOS functions [protected]
  //@{
  /** \brief Return the total free energy density of matter
      (without the rest mass contribution for the nucleons)
   */
  double free_energy_density
    (o2scl::fermion &n, o2scl::fermion &p, double T,
     o2scl::thermo &th);

  /** \brief Compute the free energy returning several 
      details as output parameters
      
      f1 is g*f_virial
      f2 is (1-g)*f_skyrme
      f3 is (1-g)*delta^2*esym
      f4 is (1-g)*delta f_hot

      so that the total homogeneous free energy is f1+f2+f3+f4
   */
  double free_energy_density_detail
  (o2scl::fermion &n, o2scl::fermion &p, double T, o2scl::thermo &th,
   std::map<std::string,double> &vdet);

  /** \brief Compute the free energy density using the virial 
      expansion including derivative information
  */
  virtual double free_energy_density_virial
    (o2scl::fermion &n, o2scl::fermion &p, double T,
     o2scl::thermo &th, double &dmundnn, double &dmundnp,
     double &dmupdnn, double &dmupdnp, double &dmundT,
     double &dmupdT);    

  /** \brief Compute the free energy density using the virial 
      expansion
  */
  virtual double free_energy_density_virial
    (o2scl::fermion &n, o2scl::fermion &p, double T,
     o2scl::thermo &th) {
    double x1, x2, x3, x4, x5, x6;
    return free_energy_density_virial(n,p,T,th,x1,x2,x3,x4,x5,x6);
  }
  
  /** \brief Alternate form of \ref free_energy_density() for
      computing derivatives

      This function does not include electrons or photons.
  */
  double free_energy_density_alt(o2scl::fermion &n, o2scl::fermion &p,
				 double nn, double np, double T,
				 o2scl::thermo &th);

  /** \brief Alternate form of \ref free_energy_density() 
      which includes electrons, positrons, and photons.
  */
  double free_energy_density_ep(double nn, double np, double T);
  
  /** \brief Compute the entropy density including photons,
      electrons, and positrons

      This function is used in \ref cs2_func() .
  */
  double entropy(o2scl::fermion &n, o2scl::fermion &p,
		 double nn, double np, double T, o2scl::thermo &th);

  /** \brief Compute energy density including photons and electons
      (without the rest mass energy density for the nucleons)
  */
  double ed(o2scl::fermion &n, o2scl::fermion &p,
	    double nn, double np, double T, o2scl::thermo &th);

  /** \brief Compute the squared speed of sound 

      The temperature should be in \f$ 1/\mathrm{fm} \f$.
   */
  double cs2_func(o2scl::fermion &n, o2scl::fermion &p, double T,
		  o2scl::thermo &th);
  //@}

  /// \name Internal variables [protected]
  //@{
  /// The table which stores the neutron star EOS results
  o2scl::table_units<> nstar_tab;

  /// The table which stores the Skyrme fits
  o2scl::table_units<> UNEDF_tab;
  
  /** \brief If true, a model has been selected (default false)
   */
  bool model_selected;
  
  /// Random number generator
  o2scl::rng<> rng;
  //@}
  
  /// \name EOS outputs
  //@{
  /// The free energy of degenerate matter
  double f_deg;
  
  /// The virial free energy
  double f_virial;
  
  /// The virial entropy
  double s_virial;

  /** \brief The value of \f$ \bar{\Lambda} \f$ for a 1.4 solar mass
      neutron star
  */
  double Lambda_bar_14;
  //@}
  
  /// \name The fit to the neutron star EOS [protected]
  //@{
  /** \brief Compute the energy density (in \f$ \mathrm{fm}^{-4} \f$)
      of neutron matter at high density from the neutron star data
      using the most recent fit (without the rest mass contribution)

      \note Currently this just returns the value of
      \ref ed_fit() .
  */
  double energy_density_ns(double nn);

  /// Parameters for the function which fits the neutron star EOS
  std::vector<double> ns_fit_parms;

  /** \brief The fit function for the energy per particle
      in units of MeV as a function of the baryon density
      (in \f$ \mathrm{fm}^{-3} \f$ )

      Note this function does not include the rest mass 
      energy density for the nucleons. 
  */
  double fit_fun(size_t np, const std::vector<double> &parms,
		 double nb);

  /** \brief The energy density (in \f$ \mathrm{fm}^{-4} \f$ )
      as a function of baryon density (in \f$ \mathrm{fm}^{-3} \f$ )

      Note this function does not include the rest mass 
      energy density for the nucleons. 
  */
  double ed_fit(double nb);
  
  /** \brief The inverse susceptibility (in \f$ \mathrm{fm}^{2} \f$ )
      as a function of baryon density (in \f$ \mathrm{fm}^{-3} \f$ )
   */
  double dmudn_fit(double nb);
  
  /** \brief The speed of sound
      as a function of baryon density (in \f$ \mathrm{fm}^{-3} \f$ )
   */
  double cs2_fit(double nb);
  
  /** \brief Compute the minimum and maximum speed of sound
      between 0.08 and \ref ns_nb_max
  */
  void min_max_cs2(double &cs2_min, double &cs2_max);

  /** \brief Fit neutron star data from Bamr to an analytical 
      expression 
  */
  void ns_fit(int row);

  /// The chi-squared for the neutron star fit
  double chi2_ns;

  /** \brief The maximum baryon density at which the neutron star
      EOS is causal

      This quantity is determined by \ref ns_fit()
  */
  double ns_nb_max;
  
  /** \brief The baryon number chemical potential (in \f$
      \mathrm{fm}^{-1} \f$ ) as a function of number density (in \f$
      \mathrm{fm}^{-3} \f$ )

      Note this function does not include the rest mass 
      for the nucleons. 
  */
  double mu_fit(double nb);
  //@}

  /// \name Parameter objects
  //@{
  o2scl::cli::parameter_int p_verbose;
  o2scl::cli::parameter_bool p_old_ns_fit;
  o2scl::cli::parameter_bool p_ns_record;
  o2scl::cli::parameter_bool p_include_muons;
  o2scl::cli::parameter_bool p_select_cs2_test;
  o2scl::cli::parameter_bool p_test_ns_cs2;
  o2scl::cli::parameter_bool p_use_alt_eos;
  o2scl::cli::parameter_double p_a_virial;
  o2scl::cli::parameter_double p_b_virial;
  o2scl::cli::parameter_int p_cs2_verbose;
  o2scl::cli::parameter_string p_nB_grid_spec;
  o2scl::cli::parameter_string p_Ye_grid_spec;
  o2scl::cli::parameter_string p_T_grid_spec;
  o2scl::cli::parameter_string p_S_grid_spec;
  o2scl::cli::parameter_string p_data_dir;
  //@}

  /// If true, then RMF fields are included
  bool rmf_fields;
  
  /// \name Other EOS functions [protected]
  //@{
  /** \brief Compute the energy density (in \f$ \mathrm{fm}^{-4} \f$)
      of neutron matter from quantum Monte Carlo (without the rest
      mass contribution)
  */
  double energy_density_qmc(double nn, double pn);

  /** \brief Construct a new neutron star EOS which ensures
      causality at high densities
  */
  int new_ns_eos(double nb, o2scl::fermion &n, double &e_ns,
		 double &densdnn);

  /** \brief Compute dfdnn including photons and electons

      This function is used in \ref cs2_func() .
   */
  double dfdnn_total(o2scl::fermion &n, o2scl::fermion &p,
		     double nn, double pn, double T, o2scl::thermo &th);
  
  /** \brief Compute dfdnp including photons and electons

      This function is used in \ref cs2_func() .
   */
  double dfdnp_total(o2scl::fermion &n, o2scl::fermion &p,
		     double nn, double pn, double T, o2scl::thermo &th);
  
  /** \brief Solve for Ye to ensure a specified value of muL at fixed T
   */
  int solve_Ye(size_t nv, const ubvector &x, ubvector &y,
		double nb, double T, double muL);
  
  /** \brief solve for a1 and a2 when cs_ns(2.0)>cs_ns(1.28)
  */
  int solve_coeff_big(size_t nv, const ubvector &x, ubvector &y, 
        double nb_last, double cs_ns_2, double cs_ns_last);

  /** \brief solve for a1 and a2 when cs_ns(2.0)<cs_ns(1.28)
   */
  int solve_coeff_small(size_t nv, const ubvector &x, ubvector &y, 
			double nb_last, double cs_ns_2, double cs_ns_last);
  
  /** \brief Internal select function
   */
  int select_internal(int i_ns_loc, int i_skyrme_loc,
		      double qmc_alpha_loc, double qmc_a_loc,
		      double eos_L_loc, double eos_S_loc, double phi_loc);
  //@}

  /// \name Particle objects [protected]
  //@{
  /// New lepton object
#ifdef O2SCL_NO_BOOST_MULTIPRECISION
  o2scl::eos_leptons elep;
#else
  o2scl::eos_leptons_multip elep;
#endif
  
  /** \brief Electron/positron
   */
  o2scl::fermion electron;

  /** \brief Muon/anti-muon
   */
  o2scl::fermion muon;

  /** \brief Photon
   */
  o2scl::boson photon;

public:
  
  /// Neutron
  o2scl::fermion neutron;

  /// Proton
  o2scl::fermion proton;

protected:

  /// Neutron for chiral part
  o2scl::fermion n_chiral;

  /// Proton for chiral part
  o2scl::fermion p_chiral;

  /// Neutrino
  o2scl::fermion neutrino;  
  //@}

  /// \name Base physics objects [protected]
  //@{
  /// The virial equation solver (now part of O2scl)
  o2scl::eos_had_virial vsd;

  /** \brief Object for computing thermodynamic integrals for leptons
   */
  o2scl::fermion_rel relf;

  /// Thermodynamic quantities
  o2scl::thermo th2;
  
  /// Thermodynamic quantities for chiral part
  o2scl::thermo th_chiral;
  
  /// Base EOS model
  o2scl::eos_had_skyrme sk;

  /// Base RMF model
  o2scl::eos_had_rmf rmf;

  /// Base RMF model
  o2scl::eos_had_rmf_hyp rmf_hyp;

  /// Skyrme model for finite-temperature correction
  o2scl::eos_had_skyrme sk_Tcorr;

  /// Pointer to EOS for finite-temperature corrections
  o2scl::eos_had_temp_eden_base *eos_Tcorr;

  /// The virial EOS
  eos_crust_virial_v2 ecv;

  /// Alternative skryme model
  o2scl::eos_had_skyrme sk_alt;

  /// Pointer to alternative model
  o2scl::eos_had_temp_base *eosp_alt;

  /// Name of the alternate EOS
  std::string alt_name;
  //@}

  /// \name The parameters for the QMC energy density [protected]
  //@{
  /// The second exponent for density in the QMC EOS (unitless)
  double qmc_beta;
  /// The second coefficient for the QMC EOS (in MeV)
  double qmc_b;
  /** \brief The saturation density of the QMC EOS, equal to 
      \f$ 0.16~\mathrm{fm}^{-3} \f$
  */
  double qmc_n0;
  //@}
  
  /// \name Output saturation properties [protected]
  //@{
  /// The binding energy per particle
  double eos_EoA;
  
  /// The incompressibility
  double eos_K;
  
  /// The saturation density
  double eos_n0;
  //@}

 public:

  /// \name Constructor and destructor
  //@{
  eos();

  virtual ~eos() {
  }
  //@}

  /// \name Settings [public]
  //@{
  /** \brief If true, use the EOS from the Du et al. (2019) paper
      instead of the Du et al. (2022) update (default false)
   */
  bool old_version;
  
  /** \brief Use an alternate EOS rather than the Du et al. 
      combined EOS (default false)
  */
  bool use_alt_eos;

  /** \brief If true, strangeness is included as an additional
      tensor rank (default false)
  */
  bool strange_axis;
  
  /** \brief If true, test the neutron star speed of sound 
      (default true)
   */
  bool test_ns_cs2;
  
  /** \brief If true, save the results of the neutron star fit to
      a file, and immediately exit (default false)
   */
  bool ns_record;

  /** \brief If true, use the old neutron star fit (default true)

      This defaults to true because the old fit performs a bit
      better than the new one. The new fit was never used
      in a publication. 
   */
  bool old_ns_fit;

  /** \brief Generic verbosity parameter (default 0)

      See also \c cs2_verbose and \c function_verbose .
   */
  int verbose;

  /// If true, create output files for individual EOSs
  bool output_files;

  /** \brief Coefficient for modulation of virial EOS (default 10.0)
   */
  double a_virial;

  /** \brief Coefficient for modulation of virial EOS (default 10.0)
   */
  double b_virial;
  
  /** \brief If true, include muons (default false)
   */
  bool include_muons;

  /** \brief If true, test cs2 in the \ref select_internal() function
      (default true)
  */
  bool select_cs2_test;
  
  /// Directory containing data files, default "data"
  std::string data_dir;
  //@}

  /// \name Command-line interface functions [public]
  //@{
  /** \brief Construct a table at fixed electron fraction

      <filename> <Ye>

      Constructs the EOS as a table3d object and outputs to 
      <filename>. Currently stores Fint, Pint, Sint, g, 
      msn, and msp.
   */
  int table_Ye(std::vector<std::string> &sv,
	       bool itive_com);

  /** \brief Load the HRG table

      <filename>

      Loads a list of resonances from a text file.
   */
  int hrg_load(std::vector<std::string> &sv,
	       bool itive_com);
  
  /** \brief Use alternate, rather than the Du et al. EOS

      <"Skyrme"> <name> or <"RMF"> <name>, etc.

      For a Skyrme model, the first argument should be the word 
      Skyrme, and the second should be the name of the desired Skyrme
      model. Similarly for an RMF model. (Hyperons are not yet
      fully supported.)
   */
  int alt_model(std::vector<std::string> &sv,
	       bool itive_com);

  /** \brief Construct a table at fixed baryon density

      <filename> <nB>

      Constructs the EOS as a table3d object and outputs to 
      <filename>. Currently stores Fint, Pint, Sint, g, 
      msn, and msp.
   */
  int table_nB(std::vector<std::string> &sv,
	       bool itive_com);

  /** \brief Construct the PNS EOS and the M-R curve

      <entropy per baryon> <lepton fraction> <output filename>

      Use YL=0 for beta equilibrium. Currently always uses a cold
      crust.
   */
  int pns_eos(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Construct a full 3D EOS table without nuclei

      <filename>

      This constructs a full 3D EOS table without nuclei using the
      specified model. The resulting file has several tensor_grid
      objects including Fint, Eint, Pint, Sint, mun, mup, cs2, mue, F,
      E, P, and S. This function does not yet support muons or
      strangeness.
   */
  int table_full(std::vector<std::string> &sv, bool itive_com);

  /** \brief Test the first derivatives of the free energy (no nuclei)

      (no parameters)

      This function tests the first derivatives of the homogeneous
      matter EOS without nuclei. The model must be selected before
      running this function. It tests derivatives over a range of
      densities for three temperatures and two electron fractions.
   */
  int test_deriv(std::vector<std::string> &sv, bool itive_com);

  /** \brief Select an EOS model

      <i_ns> <i_skyrme> <alpha> <a> <L> <S> <phi>

      Select an EOS model given the 7 specified parameters.
   */
  int select_model(std::vector<std::string> &sv, bool itive_com);

  /** \brief Compare the virial and full EOS

      Params.

      Help.
  */
  int vir_comp(std::vector<std::string> &sv, bool itive_com);

  /** \brief Evaluate the EOS at one (nB,Ye,T) point

      <nB> <Ye> <TMeV>

      Compute the EOS (without nuclei, leptons, or photons) at one
      point and output the results to the screen.
  */
  int point(std::vector<std::string> &sv, bool itive_com);

  /** \brief Select a random EOS model

      (No arguments.)

      Select a random EOS, checking several physical constraints
      and re-selecting a new random EOS until all the constraints
      are met.
   */
  int random(std::vector<std::string> &sv, bool itive_com);

  /** \brief Compute the data for the comparison figures

      Params.

      Help.
   */
  int comp_figs(std::vector<std::string> &sv, bool itive_com);

  /** \brief Compute the data for the Monte Carlo figures

      Params.

      Help.
   */
  int mcarlo_data(std::vector<std::string> &sv, bool itive_com);

  /** \brief Fit the virial EOS

      <no parameters>

      Fit the virial EOS table to the functional form specified
      in Du et al. (2019) using \ref eos_crust_virial_v2::fit() .
   */
  int vir_fit(std::vector<std::string> &sv, bool itive_com);

  /** \brief Test the electron and photon EOS

      [filename]

      This function tests the electron and photon EOS to ensure that
      it does not call the error handler (it does not test accuracy).
      It uses a larger grid than the default EOS grid and stores the
      results in tensor_grid objects. If a file is specified, these
      tensor_grid objects are output to the specified file.

      This function does not require an EOS table or any hadronic
      model specification.
   */
  int test_eg(std::vector<std::string> &sv, bool itive_com);

  /** \brief Compare to other EOSs?

      Params.

      Help.
   */
  int eos_sn(std::vector<std::string> &sv, bool itive_com);
  //@}

  /** \brief Free energy density as a function of the 
      strangeness fraction
   */
  double free_energy_density_s
  (o2scl::fermion &n, o2scl::fermion &p, double Y_s, double T,
   o2scl::thermo &th);
  
  /** \brief \brief Free energy density as a function of the 
      strangeness fraction
   */
  double free_energy_density_detail_s
  (o2scl::fermion &n, o2scl::fermion &p, double Y_s, double T,
   o2scl::thermo &th,
   std::map<std::string,double> &vdet);
  
  /// \name Miscellaneous functions [public]
  //@{
  /** \brief Solve for fixed entropy per baryon and fixed
      lepton fraction
  */
  int solve_fixed_sonb_YL(size_t nv, const ubvector &x, ubvector &y,
			  double nB, double sonb, double YL);

  /** \brief Solve for T to ensure a specified value of sonb at fixed Ye
   */
  int solve_T(size_t nv, const ubvector &x, ubvector &y,
	      double nb, double Ye, double sonb);
  
  
  /** \brief Setup the command-line interface with commands and
      parameters
   */
  virtual void setup_cli(o2scl::cli &cl, bool read_docs=true);

  /// Desc
  int comm_set(std::vector<std::string> &sv, bool itive_com);
  
  //@}
  
};

