/*
  -------------------------------------------------------------------
  
  Copyright (C) 2021-2022, Andrew W. Steiner
  
  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#ifndef INTE_CUSTOM
#define INTE_CUSTOM

#include <o2scl/inte.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_qags_gsl.h>
#include <o2scl/inte_qng_gsl.h>
#include <o2scl/mcarlo_vegas.h>
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/inte_kronrod_boost.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

namespace o2scl {

  /** Desc
   */
  template<class func_t=funct> class inte_custom : 
    public inte<func_t> {

  protected:
    
    inte_kronrod_boost<> kb;
    
    inte_qag_gsl<> qag;
    
    inte_qags_gsl<> qags;
    
    inte_qng_gsl<> qng;

    inte_qag_smooth<> iqs;
    
   // inte_adapt_cern<> ac;
    
    mcarlo_vegas<> mv;
    
    func_t *fp;
    
    bool record;
    
    double funct_wrapper(double x) {
      if (fp==0) {
        O2SCL_ERR("Function pointer not set.",o2scl::exc_einval);
      }
      double ret=(*fp)(x);
      if (record) {
        std::vector<double> line={x,ret};
        t.line_of_data(line.size(),line);
      }
      return ret;
    }

    table<> t;
  
  public:

    inte_custom() {
      fp=0;
      record=false;
      qags.tol_rel=1.0e-6;
      qags.tol_abs=0.0;
      qags.err_nonconv=false;
      qags.set_limit(250);
      mv.n_points=3000;
      mv.tol_rel=1.0e-6;
      mv.verbose=0;
    }
      
    virtual ~inte_custom() {
    }

    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    virtual int integ_err(func_t &func, double a, double b, 
                          double &res, double &err) {

      fp=&func;
      funct f=std::bind(std::mem_fn<double(double)>
                        (&inte_custom::funct_wrapper),this,
                        std::placeholders::_1);
      int ret;
      
      std::cout << "1";
      qags.tol_rel=1.0e-6;
      ret=qags.integ_err(f,a,b,res,err);
      if (ret!=0) {
        std::cout << res << " " << err << std::endl;
        iqs.qag.set_limit(100);
        iqs.integ_err(f,a,b,res,err);
        std::cout << res << " " << err << std::endl;
        exit(-1);
        std::cout << "2";
        mv.tol_rel=1.0e-6;
        ret=mv.integ_err(f,a,b,res,err);
      }
      if (ret!=0) {
        std::cout << "3";
        qags.tol_rel=1.0e-4;
        ret=qags.integ_err(f,a,b,res,err);
      }
      if (ret!=0) {
        std::cout << "4";
        mv.tol_rel=1.0e-4;
        ret=mv.integ_err(f,a,b,res,err);
      }
      if (ret!=0) {
        std::cout << "5";
        qng.tol_rel=1.0e-6;
        ret=qng.integ_err(f,a,b,res,err);
      }
      if (ret!=0) {
        std::cout << "6";
        qng.tol_rel=1.0e-4;
        ret=qng.integ_err(f,a,b,res,err);
      }
      if (ret!=0) {
        std::cout << "Failed." << std::endl;
        exit(-1);
        record=true;
        t.clear();
        t.line_of_names("x y");
        ret=qags.integ_err(f,a,b,res,err);
        std::cout << "Integration failed. Writing file." << std::endl;
        o2scl_hdf::hdf_file hf;
        hf.open_or_create("inte_custom.o2");
        hdf_output(hf,t,"table");
        hf.close();
        record=false;
        exit(-1);
      }
      return 0;
    }
           
    /// Return string denoting type ("inte_custom")
    const char *type() { return "inte_custom"; }
  
  };
  
}

#endif
