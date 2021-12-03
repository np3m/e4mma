/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2021, Xingfu Du, Zidu Lin, and Andrew W. Steiner
  
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
#include <cubature.h>

using namespace std;
using namespace o2scl;

int f(unsigned ndim, const double *x, void *fdata,
      unsigned fdim, double *fval) {
  double x1=x[0];
  double y1=x[1];
  fval[0]=exp(-pow(x1-2.0,2.0)-pow(y1-3.0,2.0));
  return 0;
}

int f2(unsigned ndim, const double *x, void *fdata,
       unsigned fdim, double *fval) {
  double t=x[0];
  double x1=1.0+t/(1.0-t);
  double dxdt=1.0/(1.0-t)/(1.0-t);
  double y1=x[1];
  fval[0]=exp(-pow(x1-2.0,2.0)-pow(y1-3.0,2.0))*dxdt;
  return 0;
}

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  double xmin[2]={0,0};
  double xmax[2]={10,10};
  double val, err;
  hcubature(1,f,0,2,xmin,xmax,0,0,1.0e-8,ERROR_INDIVIDUAL,&val,&err);
  cout << val << " " << err << endl;

  xmin[0]=0.0;
  xmax[0]=1.0;
  hcubature(1,f2,0,2,xmin,xmax,0,0,1.0e-8,ERROR_INDIVIDUAL,&val,&err);
  cout << val << " " << err << endl;
  
#ifndef NO_MPI
  // Init MPI
  MPI_Init(&argc,&argv);
#endif
  
  eos_nuclei eph;

  eph.load_nuclei();
  
  cli cl;
  
  eph.setup_cli(cl);

  cl.run_auto(argc,argv);

#ifndef NO_MPI
  // Finalize MPI
  MPI_Finalize();
#endif
  
  return 0;
}
