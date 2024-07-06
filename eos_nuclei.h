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
#include "eos.h"
#include "eos_had_skyrme_ext.h"
#include <o2scl/nucmass_fit.h>
#include <o2scl/slack_messenger.h>
#include <o2scl/part_funcs.h>
#include <o2scl/interpm_krige.h>
#include <map>
#include <filesystem>
#include <random>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;
typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_col;

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

/** New kernel for GPI object
 */
class mcovar_funct_rbf_new {
    public:
        //length parameters
        std::vector<double> len;
        std::vector<double> slope;
        std::vector<double> pos;

        //noise
        double noise;

        mcovar_funct_rbf_new() {
            noise=0.0;
        }

        //get number of parameters
        size_t get_n_params() {
            return len.size();
        }

        //Set the parameters
        template<class vec_t> void set_params(vec_t &p) {
            for(size_t j=0;j<(len.size()/3);j++) {
                len[j]=p[j];
                slope[j]=p[j+3];
                pos[j]=p[j+6];
            }
            return;
        } 

        //Covarience function
        template<class vec_t, class vec2_t>
        double operator()(const vec_t &x1, const vec2_t &x2) {
            double sum=0.0;
            bool equal=true;
            for (size_t j=0;j<len.size();j++) {
                for (size_t k=0;k<slope.size();k++) {
                    for (size_t l=0;l<pos.size();l++) {
                        if (x1[j]!=x2[j]) {
                            equal=false;
                        }
                        sum+=-(((x1[j]-x2[j])*(x1[j]-x2[j]))/((len[j]*len[j])+((slope[k]*slope[k])*(((x1[j]+x2[j])-pos[l])*((x1[j]+x2[j])-pos[l])))*2.0));
                    }
                }
            }
            if (equal) {
                return exp(sum)+noise;
            }
            return exp(sum);
        }

        /** first derivative of the covariance function wrt first argument
         */
        template<class vec_t, class vec2_t>
        double deriv(const vec_t &x1, const vec2_t &x2, size_t ix) {
            double sum=0.0;
            for (size_t j=0;j<len.size();j++) {
                sum+=-(((x1[j]-x2[j])*(x1[j]-x2[j]))/((len[j]*len[j])+((slope[j]*slope[j])*(((x1[j]+x2[j])-pos[j])*((x1[j]+x2[j])-pos[j])))*2.0));
            }
            return exp(sum)*((-(len[ix]*len[ix]+slope[ix]*slope[ix]*((x1[ix]+x2[ix])-pos[ix])*((x1[ix]-x2[ix])-pos[ix]))*(x1[ix]-x2[ix])
                    +pow((x1[ix]-x2[ix]),2.0)*slope[ix]*slope[ix]*((x1[ix]-x2[ix])-pos[ix]))
                    /pow((len[ix]*len[ix]+slope[ix]*slope[ix]*pow(((x1[ix]+x2[ix])-pos[ix]),2.0)),2.0));
        }

        /** The second derivative of the covariance function.
         */
        template<class vec_t, class vec2_t>
        double deriv2(vec_t &x1, vec2_t &x2, size_t ix, size_t iy) {
            double sum=0.0;
            for (size_t j=0;j<len.size();j++) {
                sum+=-(((x1[j]-x2[j])*(x1[j]-x2[j]))/((len[j]*len[j])+((slope[j]*slope[j])*(((x1[j]+x2[j])-pos[j])*((x1[j]+x2[j])-pos[j])))*2.0));
            }
            double diffx2=(x1[ix]-x2[ix])*(x1[ix]-x2[ix]);
            double len2=len[ix]*len[ix];
            double m2=slope[ix]*slope[ix];
            double diffxn2=((x1[ix]+x2[ix])-pos[ix])*((x1[ix]+x2[ix])-pos[ix]);
            double denom2=pow((len2+m2*diffxn2),2.0);
            if (ix==iy) {
            return exp(sum)*(((-(len2+m2*diffxn2)*(x1[ix]-x2[ix])
                    +diffx2*m2*((x1[ix]-x2[ix])-pos[ix]))
                    /denom2)
                    *((-(len2+m2*diffxn2)*(x1[ix]-x2[ix])
                    +diffx2*m2*((x1[ix]-x2[ix])-pos[ix]))
                    /denom2)
                    +((-2.0*slope[ix]*((x1[ix]-x2[ix])-pos[ix])*(x1[ix]-x2[ix])-(len2+m2*diffxn2))
                    +2.0*(x1[ix]-x2[ix])*m2*((x1[ix]-x2[ix])-pos[ix])
                    +diffx2*m2)*(1/denom2)
                    +((-(len2+m2*diffxn2)*(x1[ix]-x2[ix])
                    +diffx2*m2*((x1[ix]-x2[ix])-pos[ix]))
//                    /denom2)
                    *-2.0*pow((len2+m2*diffxn2),-3.0)*m2*2.0*((x1[ix]-x2[ix]-pos[ix]))));
            }
            return exp(sum)*(((-(len2+m2*diffxn2)*(x1[ix]-x2[ix])
                    +diffx2*m2*((x1[ix]-x2[ix])-pos[ix]))
                    /denom2)
                    *((-(len[iy]*len[iy]+slope[iy]*slope[iy]*((x1[iy]+x2[iy])-pos[iy])*((x1[iy]-x2[iy])-pos[iy]))*(x1[iy]-x2[iy])
                    +pow((x1[iy]-x2[iy]),2.0)*slope[iy]*slope[iy]*((x1[iy]-x2[iy])-pos[iy]))
                    /pow((len[iy]*len[iy]+slope[iy]*slope[iy]*pow(((x1[iy]+x2[iy])-pos[iy]),2.0)),2.0)));
        }
};


/** New kernel for GPI object with noise
 */
class mcovar_funct_quad_correl {
    public:
        //length parameters
        std::vector<double> len;
//        std::vector<double> length;
//        std::vector<double> slope;
//        std::vector<double> pos;

        //noise
        double log10_noise;

        //get number of parameters
        size_t get_n_params() {
            return len.size()+1;
        }
/*        //constructor
        mcovar_funct_quad_correl() {
            length.push_back(3);
            slope.push_back(3);
            pos.push_back(3);
        }*/

        //Set the parameters
        template<class vec_t> void set_params(vec_t &p) {
            for(size_t j=0;j<(len.size()/3);j++) {
//                length[j]=p[j];
//                slope[j]=p[j+3];
//                pos[j]=p[j+6];
                len[j]=p[j];
                len[j+3]=p[j+3];
                len[j+6]=p[j+6];
            }
            log10_noise=p[len.size()];
            return;
        } 

        //Covarience function
        template<class vec_t, class vec2_t>
        double operator()(const vec_t &x1, const vec2_t &x2) {
            double sum=0.0;
            bool equal=true;
            for (size_t j=0;j<3;j++) {
                for (size_t k=0;k<3;k++) {
                    for (size_t l=0;l<3;l++) {
                        if (x1[j]!=x2[j]) {
                            equal=false;
                        }
//                        sum+=-(((x1[j]-x2[j])*(x1[j]-x2[j]))/((length[j]*length[j])+((slope[k]*slope[k])*(((x1[j]+x2[j])-pos[l])*((x1[j]+x2[j])-pos[l])))*2.0));
//                        sum+=-(((x1[j]-x2[j])*(x1[j]-x2[j]))/((length[j]*length[j])+((slope[j]*slope[j])*(((x1[j]+x2[j])-pos[j])*((x1[j]+x2[j])-pos[j])))*2.0));
                        sum+=-(((x1[j]-x2[j])*(x1[j]-x2[j]))/((len[j]*len[j])+((len[k+3]*len[k+3])*(((x1[j]+x2[j])-len[l+6])*((x1[j]+x2[j])-len[l+6])))*2.0));
                    }
                }
            }
            if (equal) {
                return exp(sum)+pow(10.0,log10_noise);
            }
            return exp(sum);
        }

        /** first derivative of the covariance function wrt first argument
         */
        template<class vec_t, class vec2_t>
        double deriv(const vec_t &x1, const vec2_t &x2, size_t ix) {
            double sum=0.0;
            for (size_t j=0;j<3;j++) {
//                sum+=-(((x1[j]-x2[j])*(x1[j]-x2[j]))/((length[j]*length[j])+((slope[j]*slope[j])*(((x1[j]+x2[j])-pos[j])*((x1[j]+x2[j])-pos[j])))*2.0));
                sum+=-(((x1[j]-x2[j])*(x1[j]-x2[j]))/((len[j]*len[j])+((len[j+3]*len[j+3])*(((x1[j]+x2[j])-len[j+6])*((x1[j]+x2[j])-len[j+6])))*2.0));
            }
/*            return exp(sum)*((-(length[ix]*length[ix]+slope[ix]*slope[ix]*((x1[ix]+x2[ix])-pos[ix])*((x1[ix]-x2[ix])-pos[ix]))*(x1[ix]-x2[ix])
                    +pow((x1[ix]-x2[ix]),2.0)*slope[ix]*slope[ix]*((x1[ix]-x2[ix])-pos[ix]))
                    /pow((length[ix]*length[ix]+slope[ix]*slope[ix]*pow(((x1[ix]+x2[ix])-pos[ix]),2.0)),2.0));*/
            return exp(sum)*((-(len[ix]*len[ix]+len[ix+3]*len[ix+3]*((x1[ix]+x2[ix])-len[ix+6])*((x1[ix]-x2[ix])-len[ix+6]))*(x1[ix]-x2[ix])
                    +pow((x1[ix]-x2[ix]),2.0)*len[ix+3]*len[ix+3]*((x1[ix]-x2[ix])-len[ix+6]))
                    /pow((len[ix]*len[ix]+len[ix+3]*len[ix+3]*pow(((x1[ix]+x2[ix])-len[ix+6]),2.0)),2.0));
        }

        /** The second derivative of the covariance function.
         */
        template<class vec_t, class vec2_t>
        double deriv2(vec_t &x1, vec2_t &x2, size_t ix, size_t iy) {
            double sum=0.0;
            for (size_t j=0;j<3;j++) {
                for (size_t k=0;k<3;k++) {
                    for (size_t l=0;l<3;l++) {
                        if ((k == l) || (j == k) || (j == l)) {
//                            sum+=-(((x1[j]-x2[j])*(x1[j]-x2[j]))/((length[j]*length[j])+((slope[k]*slope[k])*(((x1[j]+x2[j])-pos[l])*((x1[j]+x2[j])-pos[l])))*2.0));
                            sum+=-(((x1[j]-x2[j])*(x1[j]-x2[j]))/((len[j]*len[j])+((len[k+3]*len[k+3])*(((x1[j]+x2[j])-len[l+6])*((x1[j]+x2[j])-len[l+6])))*2.0));
                        }
//                        sum+=-(((x1[j]-x2[j])*(x1[j]-x2[j]))/((length[j]*length[j])+((slope[j]*slope[j])*(((x1[j]+x2[j])-pos[j])*((x1[j]+x2[j])-pos[j])))*2.0));
                    }
                }
            }
/*            double diffx2=(x1[ix]-x2[ix])*(x1[ix]-x2[ix]);
            double len2=length[ix]*length[ix];
            double m2=slope[ix]*slope[ix];
            double diffxn2=((x1[ix]+x2[ix])-pos[ix])*((x1[ix]+x2[ix])-pos[ix]);
            double denom2=pow((len2+m2*diffxn2),2.0);
            if (ix==iy) {
            return exp(sum)*(((-(len2+m2*diffxn2)*(x1[ix]-x2[ix])
                    +diffx2*m2*((x1[ix]-x2[ix])-pos[ix]))
                    /denom2)
                    *((-(len2+m2*diffxn2)*(x1[ix]-x2[ix])
                    +diffx2*m2*((x1[ix]-x2[ix])-pos[ix]))
                    /denom2)
                    +((-2.0*slope[ix]*((x1[ix]-x2[ix])-pos[ix])*(x1[ix]-x2[ix])-(len2+m2*diffxn2))
                    +2.0*(x1[ix]-x2[ix])*m2*((x1[ix]-x2[ix])-pos[ix])
                    +diffx2*m2)*(1/denom2)
                    +((-(len2+m2*diffxn2)*(x1[ix]-x2[ix])
                    +diffx2*m2*((x1[ix]-x2[ix])-pos[ix]))
//                    /denom2)
                    *-2.0*pow((len2+m2*diffxn2),-3.0)*m2*2.0*((x1[ix]-x2[ix]-pos[ix]))));
            }
            return exp(sum)*(((-(len2+m2*diffxn2)*(x1[ix]-x2[ix])
                    +diffx2*m2*((x1[ix]-x2[ix])-pos[ix]))
                    /denom2)
                    *((-(length[iy]*length[iy]+slope[iy]*slope[iy]*((x1[iy]+x2[iy])-pos[iy])*((x1[iy]-x2[iy])-pos[iy]))*(x1[iy]-x2[iy])
                    +pow((x1[iy]-x2[iy]),2.0)*slope[iy]*slope[iy]*((x1[iy]-x2[iy])-pos[iy]))
                    /pow((length[iy]*length[iy]+slope[iy]*slope[iy]*pow(((x1[iy]+x2[iy])-pos[iy]),2.0)),2.0)));*/
            double diffx2=(x1[ix]-x2[ix])*(x1[ix]-x2[ix]);
            double len2=len[ix]*len[ix];
            double m2=len[ix+3]*len[ix+3];
            double diffxn2=((x1[ix]+x2[ix])-len[ix+6])*((x1[ix]+x2[ix])-len[ix+6]);
            double denom2=pow((len2+m2*diffxn2),2.0);
            if (ix==iy) {
            return exp(sum)*(((-(len2+m2*diffxn2)*(x1[ix]-x2[ix])
                    +diffx2*m2*((x1[ix]-x2[ix])-len[ix+6]))
                    /denom2)
                    *((-(len2+m2*diffxn2)*(x1[ix]-x2[ix])
                    +diffx2*m2*((x1[ix]-x2[ix])-len[ix+6]))
                    /denom2)
                    +((-2.0*len[ix+3]*((x1[ix]-x2[ix])-len[ix+6])*(x1[ix]-x2[ix])-(len2+m2*diffxn2))
                    +2.0*(x1[ix]-x2[ix])*m2*((x1[ix]-x2[ix])-len[ix+6])
                    +diffx2*m2)*(1/denom2)
                    +((-(len2+m2*diffxn2)*(x1[ix]-x2[ix])
                    +diffx2*m2*((x1[ix]-x2[ix])-len[ix+6]))
//                    /denom2)
                    *-2.0*pow((len2+m2*diffxn2),-3.0)*m2*2.0*((x1[ix]-x2[ix]-len[ix+6]))));
            }
            return exp(sum)*(((-(len2+m2*diffxn2)*(x1[ix]-x2[ix])
                    +diffx2*m2*((x1[ix]-x2[ix])-len[ix+6]))
                    /denom2)
                    *((-(len[iy]*len[iy]+len[iy+3]*len[iy+3]*((x1[iy]+x2[iy])-len[iy+6])*((x1[iy]-x2[iy])-len[iy+6]))*(x1[iy]-x2[iy])
                    +pow((x1[iy]-x2[iy]),2.0)*len[iy+3]*len[iy+3]*((x1[iy]-x2[iy])-len[iy+6]))
                    /pow((len[iy]*len[iy]+len[iy+3]*len[iy+3]*pow(((x1[iy]+x2[iy])-len[iy+6]),2.0)),2.0)));
        }
};

/** \brief Specialized Gaussian process interpolation object
 */
class interpm_krige_eos :
  public o2scl::interpm_krige_optim
<std::vector<mcovar_funct_quad_correl>,ubvector,ubmatrix,ubmatrix_row,
 ubmatrix,ubmatrix_row,Eigen::MatrixXd,
 o2scl_linalg::matrix_invert_det_eigen<Eigen::MatrixXd>> {
  
public:
  
  // Typedefs from the specialized template parameters
  typedef ubmatrix_row mat_y_row_t;
  typedef ubmatrix_row mat_x_row_t;
  typedef Eigen::MatrixXd mat_inv_kxx_t;
  typedef o2scl_linalg::matrix_invert_det_eigen<Eigen::MatrixXd>
  mat_inv_t;

  /** \brief List of calibration points
   */
  std::vector<size_t> calib_list;

  /** \brief List of points to fix
   */
  std::vector<size_t> fix_list;

  /// \name Grids
  //@{
  std::vector<double> nB_grid;
  std::vector<double> Ye_grid;
  std::vector<double> T_grid;
  //@}

  /// \name Tensor grid data
  //@{
  o2scl::tensor_grid<> *tgp_F;
  o2scl::tensor_grid<> *tgp_P;
  o2scl::tensor_grid<> *tgp_S;
  o2scl::tensor_grid<> *tgp_mun;
  o2scl::tensor_grid<> *tgp_mup;
  o2scl::tensor_grid<> *tgp_mue;
  o2scl::tensor_grid<> *tgp_Fint;
  o2scl::tensor_grid<> *tgp_Sint;
  o2scl::tensor_grid<> *tgp_Fall;
  o2scl::tensor_grid<> *tgp_Sall;
  o2scl::tensor_grid<> tgp_cs2;
  //@}

  /// Neutron mass
  double mneut;
  /// Proton mass
  double mprot;
  /// Include electrons and positrons in calculation of the free energy
  bool includeMue;

  /// Lepton EOS object
  o2scl::eos_leptons elep;

  /// If true, interpolate Fint, otherwise, interpolate F
  bool interp_Fint;
  
  interpm_krige_eos() {
    elep.include_photons=true;
    interp_Fint=false;
    //full_min=true;
  } 

  /** \brief Set the interpolator given the specified EOS
      objects
   */
  virtual void set(std::vector<double> &nB_grid2,
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
                   double mn, double mpx);
  
  /** \brief Additional constraints for the interpolation
   */
  virtual int addl_const(size_t iout, double &ret);

  /** \brief calculate interpolated values for points in list
   */
  std::map<std::vector<size_t>, std::vector<double>> interpolate_points(std::vector<size_t> points_list);

};

/** \brief Solve for the EOS including nuclei

    \verbatim embed:rst

    .. todo::

    Class eos_nuclei:

    - Rename n_nB2 to n_nB, etc.
    - Make child of eos_sn_base
    - Use \c loaded instead of testing n_nB2==0
    - Move select_high_T to the parent
    - 1/15/22: I'm not sure if fix_cc() is really useful or 
    not anymore?
    - 1/15/22: I'm not sure the non-derivative virial solver
    is needed anymore?
    - Future: Allow different form for the NS fit
    - Future: Allow user to specify where data files are
    located or to manually specify the nuclear model?
    - Future: Allow RMF rather than just Skyrme.

    \endverbatim

    \comment
    MUSES list:
    - convert all 3-rank tensors to multi-rank tensors (done)
    - convert get3 to get and set3 to set everywhere (done)
    - add sigma, omega, rho, delta? fields
    - Fix tensor_list HDF5 output
    - add beta-equilibrium functions
    - add RMF models
    - finish implementing 'strangeness' parameter
    - add RMFH models
    - Add muons out of equilibrium
    - improve fermion pair_density in non-degenerate limit
    \endcomment
*/
class eos_nuclei : public eos {

public:

  /// \name Constructor and destructor
  //@{
  eos_nuclei();
  
  virtual ~eos_nuclei();
  //@}

  /// Partition functions for the nuclei
  o2scl::part_funcs pfuncs;

  /// Particle object for bosons
  o2scl::boson_rel relb;
  
  /// Particle object for bosons
  o2scl::boson_eff effb;

  /** \brief If true, include a hadron resonance gas
   */
  bool inc_hrg;

  /** \brief If true, will save points with a baryon density under a value
   * set in max_nB_inter to csv file for later correction with the interpolation
   * function.
   */
  bool save_to_csv;

  /** \brief Maximum baryon density of points with acausal speeds of sound squared
   * to be saved in a csv file for later corrections with the interpolation function.
   */
  double max_nB_inter;

  /** \brief Location where csv file will be saved
   */
  std::string csv_location;

  /** \brief Location of eos table file
   */
  std::string table_path;

  /** \brief Solve for charge neutrality and fixed baryon fraction
      with a hadron resonance gas
   */
  int solve_hrg(size_t nv, const ubvector &x,
                ubvector &y, double nB, double Ye, double T);

  /** \brief Used to change values in hdf5 i/o test function
   */
  void change_tgp(o2scl::tensor_grid<>& tg_file, double value);

  /** \brief Test input and output to o2 files.
   */
  int test_hdf5io ();

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
   */
  std::string S_grid_spec;
  //@}

  /// \name Nuclear masses and their fits
  //@{
  /** \brief Fit the FRDM mass model

      <no parameters>

      Fit the FRDM mass model.
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
      temperature, in units of \f$ \mathrm{fm} \f$
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
  /** \brief The time (in seconds) between output file updates for
      \ref generate_table() (default 1800)
  */
  double file_update_time;
  
  /// The number of iterations between file updates (default 100000)
  int file_update_iters;
    
  /** \brief The maximum value of A for a fixed distribution when
      alg_mode is 4 (default 600)
  */
  int fd_A_max;

  /** \brief If true, attempt to extend FRDM beyond
      the drip lines (default false).
  */
  bool extend_frdm;
  
  /// The maximum value of \f$ N/Z \f$ or \f$ Z/N \f$ (default 7)
  double max_ratio;
  
  /** \brief Relative tolerance for the solver in the \ref
      eos_fixed_dist() function (default \f$ 10^{-6} \f$)
  */ 
  double mh_tol_rel;

  /** \brief Filename containing separate table to use as a guess for
      the \c generate-table command (default is the empty string)
  */
  std::string ext_guess;
  
  /** \brief Function for delta Z and delta N in the single 
      nucleus approximation
      
      This is the function used to specify the range of Z and N which
      is varied to find the smallest free energy
  */
  std::string nucleon_func;

  /** \brief Maximum time, in seconds, for the ``generate-table``
      command. The default is 0.0 which is interpreted as no maximum
      time
  */
  double max_time;
  
  /** \brief If true, show all nuclei considered at every point (default 
      false)
      
      This applies to the \c point-nuclei command and the eos_vary_ZN()
      function.
  */
  bool show_all_nuclei;

  /** \brief If true, ensure that the nuclear radius is less than the
      Wigner-Seitz radius (default true)
  */
  bool rnuc_less_rws;

  /// Units of objects in the vdet arrays
  std::map<std::string,std::string> vdet_units;
  
  /** \brief If true, recompute all points, irrespective of the
      value of the convergence flag (default false)

      This setting is used in \c point-nuclei and \c generate-table.
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

  /** \brief Algorithm mode (default 4)

      0 for SNA, 1 for old SNA method, 2 for vary dist., 
      3 for old vary dist., and 4 for fixed dist.
  */
  int alg_mode;

  /** \brief Algorithm for \ref eos_fixed_dist()

      Modify the algorithm for the \ref eos_fixed_dist() function. The
      1s digit is the number of solves minus 1, the 10s digit is the
      number of brackets divided by 10, the 100s digit is the number
      of minimizes, and the 1000s digit is the number of random
      guesses to try divided by 1000. The default is 1111. Other good
      options are 1319, 1919, 1999, and 9999.
  */
  int fixed_dist_alg;
  
  /** \brief If true, when computing a point, perform the calculation
      using some of the neighboring points as an initial guess
      (default 0)

      Values greater than 0 use the point at the next smallest
      density, values greater than 1 use the point at the next largest
      density, values greater than 2 use points at the next largest
      and next smallest temperature, and values greater than 4 use the
      next largest and smallest electron fraction.
  */
  int six_neighbors;

  /** \brief A new function verbose parameter

      Verbose for individual functions
      (default value 11111).\n\t1s digit: fixed_ZN()\n\t10s digit:
      vary_ZN()\n\t100s digit:
      fixed_dist()\n\t1000s digit: vary_dist()\n\t10000s digit:
      store_point().
  */
  int function_verbose;
  
  void store_hrg(double mun, double mup,
                 double nn, double np, double T,
                 o2scl::table_units<> &tab);
  
  /** \brief If true, use previously computed points (or guesses) as
      an initial guess to compute adjacent points (default true)
  */
  bool propagate_points;

  /** \brief If true, survey the nB and Ye equations near a failed point
      (default false)
  */
  bool survey_eqs;

  /** \brief If true, output all of the data necessary for a full EOS

      If true, include Eint, Pint, Sint, mun, and mup. If include_eg
      is additionally true, then add F, E, P, and S.
  */
  bool derivs_computed;

  /** \brief Always true, included for consistency with o2scl::eos_sn_base
   */
  bool baryons_only;

  /** \brief If true, include electrons and photons

      Requires that derivs_computed is also true
  */
  bool with_leptons;
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

  /** \brief The list of electron fractions to consider
      for the \c generate-table command. 

      Can be several comma-separated ranges e.g. "1-3,5-7,59-60".
      Zero-indexed.
  */
  std::string Ye_list;
  //@}

  /// \name Other internal physics objects
  //@{
  /// Extended Skyrme model for finite-temperature corrections
  eos_had_skyrme_ext skyrme_ext;

  /// Extended Skyrme model for finite-temperature corrections
  eos_had_lim_holt lim_holt;
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
  o2scl::tensor_grid<> tg_log_xn;
  o2scl::tensor_grid<> tg_log_xp;
  o2scl::tensor_grid<> tg_flag;
  o2scl::tensor_grid<> tg_F;
  o2scl::tensor_grid<> tg_E;
  o2scl::tensor_grid<> tg_P;
  o2scl::tensor_grid<> tg_S;
  o2scl::tensor_grid<> tg_Fint;
  o2scl::tensor_grid<> tg_Eint;
  o2scl::tensor_grid<> tg_Pint;
  o2scl::tensor_grid<> tg_Sint;
  o2scl::tensor_grid<> tg_mun;
  o2scl::tensor_grid<> tg_mup;
  o2scl::tensor_grid<> tg_mue;
  o2scl::tensor_grid<> tg_Z;
  o2scl::tensor_grid<> tg_A;
  o2scl::tensor_grid<> tg_Xn;
  o2scl::tensor_grid<> tg_Xp;
  o2scl::tensor_grid<> tg_Xalpha;
  o2scl::tensor_grid<> tg_Xnuclei;
  o2scl::tensor_grid<> tg_Ymu;
  o2scl::tensor_grid<> tg_Xd;
  o2scl::tensor_grid<> tg_Xt;
  o2scl::tensor_grid<> tg_XHe3;
  o2scl::tensor_grid<> tg_XLi4;
  o2scl::tensor_grid<> tg_A_min;
  o2scl::tensor_grid<> tg_A_max;
  o2scl::tensor_grid<> tg_NmZ_min;
  o2scl::tensor_grid<> tg_NmZ_max;
  //@}

  /// \name Detail storage
  //@{
  bool include_detail;
  o2scl::tensor_grid<> tg_zn;
  o2scl::tensor_grid<> tg_zp;
  o2scl::tensor_grid<> tg_F1;
  o2scl::tensor_grid<> tg_F2;
  o2scl::tensor_grid<> tg_F3;
  o2scl::tensor_grid<> tg_F4;
  o2scl::tensor_grid<> tg_Un;
  o2scl::tensor_grid<> tg_Up;
  o2scl::tensor_grid<> tg_msn;
  o2scl::tensor_grid<> tg_msp;
  o2scl::tensor_grid<> tg_g;
  o2scl::tensor_grid<> tg_dgdT;
  o2scl::tensor_grid<> tg_sigma;
  o2scl::tensor_grid<> tg_omega;
  o2scl::tensor_grid<> tg_rho;
  //@}

  /// \name Other parameter objects
  //@{
  o2scl::cli::parameter_bool p_inc_hrg;
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
  o2scl::cli::parameter_bool p_strange_axis;
  o2scl::cli::parameter_double p_mh_tol_rel;
  o2scl::cli::parameter_double p_max_time;
  o2scl::cli::parameter_string p_nucleon_func;
  o2scl::cli::parameter_int p_alg_mode;
  o2scl::cli::parameter_int p_cs2_verbose;
  o2scl::cli::parameter_int p_fixed_dist_alg;
  o2scl::cli::parameter_int p_function_verbose;
  o2scl::cli::parameter_string p_Ye_list;
  o2scl::cli::parameter_double p_max_ratio;
  o2scl::cli::parameter_double p_file_update_time;
  o2scl::cli::parameter_int p_file_update_iters;
  o2scl::cli::parameter_string p_nB_grid_spec;
  o2scl::cli::parameter_string p_Ye_grid_spec;
  o2scl::cli::parameter_string p_T_grid_spec;
  o2scl::cli::parameter_bool p_save_to_csv;
  o2scl::cli::parameter_double p_max_nB_inter;
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

  /** \brief Compute muons in nuclear matter
   */
  int nuc_matter_muons(size_t nv, const ubvector &x, ubvector &y,
                       double nB, double Ye, double T,
                       std::map<std::string,double> &vdet);
  
  /** \brief Compute muons
   */
  int new_muons(size_t nv, const ubvector &x, ubvector &y,
                double nB, double Ye, double T,
                std::map<std::string,double> &vdet,
                o2scl::eos_sn_base &eso);
  
  /** \brief Determine the EOS presuming a distribution of nuclei
      and optimizing the limits in A and \f$ N-Z \f$
  */
  int eos_vary_dist(double nB, double Ye, double T, double &log_xn,
		    double &log_xp, double &Zbar, double &Nbar, 
		    o2scl::thermo &thx, double &mun_full, double &mup_full,
		    int &A_min, int &A_max, int &NmZ_min, int &NmZ_max,
		    std::map<std::string,double> &vdet, bool dist_changed,
		    bool no_nuclei);
   
  /** \brief Generate an EOS table

      [out file]

      Help.
  */
  int generate_table(std::vector<std::string> &sv, bool itive_com);
  //@}

  /// \name EOS post-processing functions
  //@{
  /** \brief Compute derivatives numerically

      <no parameters>

      Help.
  */
  int eos_deriv(std::vector<std::string> &sv, bool itive_com);

  /** \brief Interpolates the EOS around a point
   */
  void interpolate(double nB_cent, double Ye_cent, double T_cent, int window, int neighborhood, std::string st_o2, o2scl::tensor_grid<> &tg_cs2, o2scl::tensor_grid<> &tg_file, bool itive_com);

  /** \brief Interpolate the EOS around a specified point

      <nB> <Ye> <T MeV> <window> [st.o2] [file.o2]

      This function requires that an EOS with leptons has been
      loaded. 
   */
  int interp_point(std::vector<std::string> &sv, bool itive_com);

  /** \brief Interpolate the EOS around the points specified ina csv file.
   *
   * <st.o2> <fixed file.o2> <window> [csv file]
   *
   * This function requires that an EOS with leptons be loaded
   */
  int interp_file(std::vector<std::string> &sv, bool itive_com);

  /** \brief finds distance between two size_t vectors of arbitrary length as long as both vectors have the same length.
   */
  template<typename T>
  std::pair<double, T> vector_distance(std::vector<size_t> start, std::map<std::vector<size_t>, T> points);
  /** \brief calculate results from table for points in list
   */
  std::map<std::vector<size_t>, std::vector<double>> calculate_table_values(std::vector<size_t> points_list, interpm_krige_eos& ike);

  /** \brief minimize parameters for interpolation object using LHS and random sampling.
   */
  void latin_hypercube_sampling(int dim, size_t samples, interpm_krige_eos& ike);

  /** \brief minimize parameters for interpolation object using grid search.
   */
  void minimize_parameters(interpm_krige_eos& ike);

  /** \brief Desc
   */
  int eos_deriv_v2(std::vector<std::string> &sv, bool itive_com);

  /** \brief Compute second derivatives numerically.
   */
  int eos_second_deriv(std::vector<std::string> &sv, bool itive_com);

  /** \brief Add electrons and photons

      <no parameters>
      
      Help.
  */
  int add_eg(std::vector<std::string> &sv, bool itive_com);

  /** \brief Construct an electrons and photon table

      <output file>

      Help.
  */
  int eg_table(std::vector<std::string> &sv, bool itive_com);

  /** \brief Edit an EOS table

      <select func.> [tensor to modify] [value func.]

      The \c edit-data command counts the number of (nB,Ye,T) points
      matching the criteria specified in <select func.>. If the
      remaining two arguments are given, then the values of [tensor to
      modify] for the selected points are changed to the result of the
      function [value func.].
  */
  int edit_data(std::vector<std::string> &sv, bool itive_com);

  /** \brief Merge two output tables to create a third

      <input file 1> <input file 2> <output file>

      Tables can only be merged if their grids and settings match. If
      the Fint table is anomalously small or large or not-finite, then
      this function calls the error handler. Otherwise, for each point
      in (nB,Ye,T), there are four reasons for which a point is copied
      from the second table to the first: (i) they both have flag=10
      but the second has a smaller Fint, (ii) the second has flag=10
      but the first does not, (iii) they both have flags less than 10
      but the second has a non-zero flag with a smaller Fint, or (iv)
      the second table has a non-zero flag and the first does not.
      After the merge, the number of points modified is reported.
  */
  int merge_tables(std::vector<std::string> &sv, bool itive_com);

  /** \brief Compare two output tables

      <input file 1> <input file 2> [quantity]

      Compare two EOS tables. If the optional argument ")+
      is unspecified, then all quantities are compared. If [quantity] "+
      is specified, then only that particular quantitiy is compared. "+
      Only points for which flag=10 in both tables are compared. "+
      If derivs_computed is true, then Pint, mun, and "+
      mup are available for comparisons. If with_leptons is "+
      true, then "+
      F, E, P, and S, are also available for comparisons. Any current "+
      EOS data stored is cleared before the comparison. If the "+
      nB, Ye, or T grids do not match, then no comparison is performed.
  */
  int compare_tables(std::vector<std::string> &sv, bool itive_com);

  /** \brief Output convergence statistics and simple checks

      <no parameters>

      If an EOS is loaded, this function counts
      the number of points with each flag value, checks that
      the nuclear fractions add up to 1, checks that the free energy
      internal energy, and entropy are consistent, and checks the
      thermodynamic identity.
  */
  int stats(std::vector<std::string> &sv, bool itive_com);

  /** \brief Compute and/or show EOS results at one (n_B,Y_e,T) point

      <n_B> <Y_e> <T (MeV)> [log(xn) log(xp) Z N] [alg_mode 2-4:
      log(xn) log(xp) A_min A_max NmZ_min NmZ_max] [fname]

      If an EOS is loaded, then the n_B, Y_e, and T
      values are modified to ensure that they lie on a grid point.
      If an initial guess is specified on the command line, it is
      used even if there is a good guess already in the table.
      If the flag is not 10 or if \ref recompute is true, then the EOS is
      recomputed. If an EOS is loaded or the recompute was successful,
      then the results are output to the screen. If the point was
      successful it is stored in the current tables. If \ref show_all_nuclei
      is true, then a file named \c dist.o2 is created
      which holds the full nuclear distribution.
  */
  int point_nuclei(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Test an EOS at random points in (nB,Ye,T)

      <n_tests> [\"lg\"]

      This function tests the EOS at randomly chosen points in
      (nB,Ye,T) space. If the new calculation and the stored result
      disagree, then the new result is stored in the table. A table
      must be loaded and the full model must be specified (either with
      select-model or from the table).
      
      If the additional argument "lg" is given, then the random points
      are all selected at densities between 0.01 and 0.15 1/fm^3.
  */
  int test_random(std::vector<std::string> &sv, bool itive_com);
  //@}

  /// \name File I/O functions
  //@{
  /** \brief Load an EOS table

      <filename> 

      Loads an EOS table in to memory. In the case
      where MPI is used, only one MPI rank reads the table at a time.
  */
  int load(std::vector<std::string> &sv, bool itive_com);

  /** \brief Output an EOS table to a file

      <filename>

      Loads an EOS table in to memory. In the case
      where MPI is used, only one MPI rank writes the table at a time.
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

      <output file>

      Help.
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
  virtual void setup_cli_nuclei(o2scl::cli &cli); 

  /** \brief Initialize tensors for a new EOS table
   */
  void new_table();

  /** \brief Check the virial EOS

      <no parameters>

      This function checks the solver by using it to compute the EOS
      over a wide range of densities and temperatures This function is
      particularly good for checking to make sure that the virial part
      of the EOS does not lead to any discontinuities. For example,
      o2graph -read check_virial.o2 zn -Ye-slice 0.3 -set logz 1
      -den-plot slice -show.

      This function creates a file 'check_virial.o2'
      with four tensor_grid objects which store the neutron and
      proton fugacities. This file can be plotted with, e.g.
      o2graph -set logz 1 -read check_virial.o2 zn -set logx 1
      -set logy 1 -set colbar 1 -to-table3d 0 2 slice 0.01
      -den-plot slice -show.
  */
  int check_virial(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Use low densities to improve results at high densities

      <nB low> <nB high> <Ye low> <Ye high> <T low> <T high> <output
      file>

      This function computes the EOS at higher densities
      using initial guess from lower densities. It is particularly
      helpful in getting the phase transition between nuclei and
      nuclear matter correct. The outermost loop is temperature, the
      second loop is electron fraction and the inner loop is density.
      This function requires a table has been loaded and the EOS is
      specified. It has no parallelization support.
  */
  int increase_density(std::vector<std::string> &sv, bool itive_com);

  /** \brief Increase nB to optimize the phase transition

      <output file>

      Help
  */
  int fix_cc(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Verify the EOS

      "random" <n_tests> <output file>, "random_lg" <n_tests>
      <output file>, "all" <output file>, "all_lg" <output
      file>, or "point" <output file> <nB> <Ye> <T>

      Verify the EOS, recompute if a point fails and the write final
      results to the specified output file. An table must be loaded
      and the full model must be specified (either with select-model
      or from the table).

      This function only
      verifies that the baryon density and electron fraction equations
      are solved to within the current tolerance and does not attempt
      to solve them again. The test-random function is different, it
      actually re-solves the equations to show the answer is correct.
      Thus, this function requires a bit less running time at each
      point.

      The first argument is a mode parameter which determines which
      points will be verified. If it is "random", then random points
      will be verified. If it is "random_lg", then random points with
      baryon densities from 0.01 to 0.16 fm^{-3} will be verified. If
      is is "all", then all points will be verified. If it is "all_lg"
      then all points with baryon densities from 0.01 to 0.15 fm^{-3}
      will be verified. If it is "point", then only the specified
      point will be verified. If any point fails, the verification
      procedure stops immediately and no further points are tested.

      The baryon density, electron fraction, and temperature are
      output to the screen, as well as the integer from \ref
      eos_vary_dist(), which is 0 for success.

      This function does not appear to require electrons and 
      photons but only works with alg_mode=4. 
  */
  int verify(std::vector<std::string> &sv, bool itive_com);

  /** \brief Monte Carlo results with nuclei
      
      Params.
      
      Help.
  */
  int mcarlo_nuclei(std::vector<std::string> &sv, bool itive_com);

  /** \brief Monte Carlo results with nuclei (v2)

      <nB> <Ye> <T> <N> <filename>

      Help
  */
  int mcarlo_nuclei2(std::vector<std::string> &sv, bool itive_com);
  
  /** \brief Monte Carlo neutrino opacity in beta equilibrium

      <filename> [n_point]

      Help
  */
  int mcarlo_beta(std::vector<std::string> &sv, bool itive_com);

  /** \brief Monte Carlo neutrino opacity in pure neutron matter

      <filename> [n_point]

      Help
  */
  int mcarlo_neutron(std::vector<std::string> &sv, bool itive_com);

  /// Compute the baryon number fractions and put them in \c X
  void compute_X(double nB, ubvector &X);
  
  /** \brief Select high-temperature EOS
      
      \note This function is called by the constructor and thus
      cannot be virtual
  */
  int select_high_T_internal(int option);

  /** \brief Compute the second derivatives and the stability matrix
      
      <output file> or <nB> <Ye> <T>

      This command computes the second derivatives, speed of sound,
      and stability matrix of the current EOS table. A table with
      electrons must be loaded and the full model must be specified
      (either with select-model or from the table).

      If one output file argument is specified, the \c stability
      command creates an output file with several additional data
      objects for the second derivatives of the free energy, the
      eigenvalues of the curvature matrix, and the squared speed of
      sound.

      Otherwise, if a density, electron fraction, and temperature are
      specified, the \c stability command compares the heterogeneous
      and homogeneous matter sound speeds at the grid point nearest
      to the specified values.
  */
  int stability(std::vector<std::string> &sv,
		bool itive_com);
  
  /** \brief Select the high-temperature Skyrme EOS

      <index>

      Select 0 for the original DSH fit, 1 for NRAPR, 
      2 for Sk chi 414, 3 for Skchi450, 4 for Skchi500, 5 for ?, "+
      and 6 for Sk chi m* (the default).
  */
  int select_high_T(std::vector<std::string> &sv, bool itive_com);

  /** \brief Compute eos with nuclei by searching minimum
   */
  double f_min_search(size_t nvar,const ubvector &x,
		      double nb, double ye, double T);
  
  /** \brief Initialization for differential evolution approach
   */		      
  int init_function(size_t dim, const ubvector &x, ubvector &y);
  
  /** \brief Maxwell construction

      Params.

      Help.
  */
  int maxwell(std::vector<std::string> &sv, bool itive_com);

  /** \brief Desc
   */
  int max_fun(size_t nv, const ubvector &x, ubvector &y,
              o2scl::interp_vec<std::vector<double>,ubvector> &itp_P, 
              o2scl::interp_vec<std::vector<double>,ubvector> &itp_mun);
  
  //@}
  
};

