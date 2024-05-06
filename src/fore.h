#include <iostream>
#include <string>
#include <vector>

#include <o2scl/constants.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/inte_gauss_cern.h>
#include <o2scl/inte_qagiu_gsl.h>
#include <o2scl/min_cern.h>
#include <o2scl/root_cern.h>
#include <o2scl/root_brent_gsl.h>

#include <functional>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/test_mgr.h>
#include <o2scl/boson.h>
#include <o2scl/fermion.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef std::function<std::vector<double>> funcv;

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

class fore {

protected:
  
public:

  int verbose;

  double hbar, hbar_c, n_0, n_0_MeV3, m_n, m_p, m_e, m_mu, m_pi;
  double Y_p, n_B, mu_n, mu_p, meff_n, meff_p, U_n, U_p,mu_pi, n_e, n_pi;
  int flag = 1; // 1 when values are computed, 0 for condensation
  double Y_pi; // Pion abundance n_pi/n_B
  double e_pi; // Energy density of pions
  double s_pi; // entropy of pions
  double press_pi; // pressure of pions

  ubmatrix nucleon_mod;

  ubvector com_momentum, com_energy, se_list, p_list_New, d3_sum, d3d1_sum;

  // Phase shifts with isospin=n/2 and angular momentum=l
  // m => -1/2 spin instead of +1/2
  // (\alpha_0^{1/2}, \alpha_1^{1/2}, \alpha_0^{3/2}, \alpha_1^{3/2})
  vector<double> d30, d31, d31m, d10, d11, d11m, p_list, pseudo_pot_params;

  // functions with same names and utility as the python variant
  // returns total phase shifts for neutrons or protons for a given com_mom
  funct neutron_phase_shift_sum, proton_phase_shift_sum;

  o2scl::interp_vec<> neutron_sum, proton_sum, interp_se;

  bool include_p_wave, const_after_data;

  size_t n_nB, n_Ye, n_T;
  bool con_exist;
  o2scl::tensor_grid<> pi_con;

  // The integration methods
  // Fixed-order Gaussian quadrature of order 8-16
  o2scl::inte_gauss_cern<> fixed_quad;
  // Basic Gauss-Kronrod integration class (GSL)
  o2scl::inte_kronrod_boost<> kron;
  o2scl::inte_qagiu_gsl<> gu;
  
  // One-dimensional minimization (CERNLIB)
  o2scl::min_cern<> mcn;

  // heaviside function
  int heaviside(double x);

  // Nucleon distribution function for a non-rel fermion.
  // Depends on the mass, temperature and chemical potential
  // of the fermion.
  funct NR_fermion_distro(double m, double T, double mu);
  
  // Read pion-nucleon phase shifts from data files
  // and compute center of mass energy and momentum
  void get_phase_shifts();

  void interp_phase_shift_sum(std::vector<double> params);
  
  double inner_integral(double p, double k, funct phase_shift_sum,
                        double m_N, double last_val);
  
  // Computes self energy for a pion momentum p
  double self_energy(double p, std::vector<double> params, ubmatrix nuc_mod, double T);

  // Sets up the interpolation vector and returns self energy for a given pion 
  // momentum through interpolation
  void self_energy_interp(vector<double> params, ubmatrix nuc_mod, double T);
  double se_func(double p);
  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------
  // Nucleon distribution function for a non-rel boson.
  funct rel_boson_distro(double m, double T, double mu, funct sigma);

  //Return True if the given self energy has a minimum less than $\mu$.

  //The current self energy tends to have a single minimum in the
  //hundreds of MeV. This translates to a minimum in the pion energy
  //either at p=0 or near the minimum of $\Sigma_\pi$. So we will search
  //for the location of the potential minimum at finite momentum
  //starting the minimum function near the minimum of $\Sigma_\pi$ which
  //should only have one minima and another to search near zero.
  double condensation_exists(funct sigma_pi, double mu);

  // Computes the number density of interacting pions
  double rel_pion_number_density(double T, double mu, vector<double> params, ubmatrix nuc_mod);

  // Computes the number density of non-interacting pions
  double non_int_rel_pion_number_density(double T, double mu);

  // Calculate entropy of bosons given their distribution function.
  // Input is unitless distribution function which takes in momentum in
  // MeV Output is MeV^3
  double pion_entropy(double T, double mu, vector<double> params, ubmatrix nuc_mod);

  // Energy of pions in medium given by nucleon_mod using 
  // pseudo-potential.
  //  Function includes interaction energy.
  double pion_energy(double T, double mu, vector<double> params, ubmatrix nuc_mod);

  // Computes non-interacting boson entropy without self-energy contributions
  // given their distribution function.
  double boson_entropy(funct distro);

  void single_point_data(double Y_p, double T, double n_B, double mu_n, 
                  double mu_p, double meff_n, double meff_p);

  int calc_mu(boson &b, fermion &n, fermion &p, double T, double n_B);

  void load_pion();

  void con_chk();
};          