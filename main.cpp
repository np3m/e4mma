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
#include "eos_nuclei.h"

using namespace std;
using namespace o2scl;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

#ifndef NO_MPI
  // Init MPI
  MPI_Init(&argc,&argv);
#endif
  
  eos_nuclei eph;
  
  // The command line interface needs the data directory to get the
  // help text for all the commands, so we have to set the data
  // directory first, before initializing the cli class.
  for(int i=0;i<argc;i++) {
    
    if (i+2<argc && ((std::string)argv[i])=="-set" &&
        ((std::string)argv[i+1])=="data_dir") {
      eph.data_dir=(std::string)argv[i+2];
      cout << "Setting data_dir to " << eph.data_dir << endl;
    }

  }
    
  cli cl;
  
  eph.setup_cli_nuclei(cl);

  eph.load_nuclei();
  
  cl.run_auto(argc,argv);

#ifndef NO_MPI
  // Finalize MPI
  MPI_Finalize();
#endif
  
  return 0;
}
