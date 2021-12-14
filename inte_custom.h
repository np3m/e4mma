/*
  -------------------------------------------------------------------
  
  Copyright (C) 2021, Andrew W. Steiner
  
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
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

namespace o2scl {

  /** Desc
   */
  template<class func_t=funct> class inte_custom : 
    public inte<func_t> {

  protected:

  inte_kronrod_boost kb;
  
  inte_qag_gsl qag;
  
  inte_qags_gsl qags;
  
  inte_qng_gsl qng;
  
  inte_adapt_cern ac;
  
  func_t *fp;

  bool record;
  
  double funct_wrapper(double x) {
    if (fp==0) {
      O2SCL_ERR("Function pointer not set.",o2scl::exc_einval);
    }
    double ret=(*fp)(x);
    if (record) {
      vector<double> line={x,ret};
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
    qags.tol_abs=1.0e-6;
    qags.err_nonconv=false;
    qags.set_limit(100);
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
    int ret1=qags.integ_err(f,a,b,res,err);
    if (ret1!=0) {
      record=true;
      t.clear();
      t.line_of_names("x y");
      ret1=qags.integ_err(f,a,b,res,err);
      cout << "Integration failed. Writing file." << endl;
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
