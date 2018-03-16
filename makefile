# LIBS is the list of libraries
# LFC is the local fortran compiler
# LCXX is the local C++ compiler
# LCFLAGS are the local C++ compiler flags
# LFFLAGS are the local fortran compiler flags

# ----------------------------------------------------------------
# Various user-specific settings
# ----------------------------------------------------------------

ifeq ($(HOSTNAME),isospin.roam.utk.edu)

ifeq ($(USER),awsteiner)

# On isospin for Andrew
LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
	-lo2scl_eos -lo2scl_part -lo2scl_hdf -lo2scl -lhdf5 -lgsl
LCXX = g++ 
LCFLAGS = -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include \
	-I/usr/include/eigen3 -Wno-deprecated-declarations \
	-O3 -std=c++11 -DNO_MPI -Wshadow -DO2SCL_HDF5_COMP 
LFFLAGS = -O3 
LMPICXX = mpic++

else

# On isospin for Xingfu
LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
	-lo2scl_eos -lo2scl_part -lo2scl_hdf -lo2scl -lhdf5 -lgsl
LCXX = g++ 
LCFLAGS = -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include \
	-I/usr/include/eigen3 -Wno-deprecated-declarations \
	-O3 -std=c++11 -DNO_MPI -Wshadow -DO2SCL_HDF5_COMP
LFFLAGS = -O3 
LMPICXX = mpic++

endif

else
ifeq ($(USER),x5a)

# On mimosa
LCXX = g++ 
LIBS = -L$(O2SCL_LIB) -L$(GSL_LIB) -L$(HDF5_LIB) \
	-lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl -lhdf5 -lgsl
LCFLAGS = -O3 -std=c++11 -I$(O2SCL_INC) -I$(EIGEN_INC) -I$(GSL_INC) \
	-Wno-deprecated-declarations -I$(HDF5_INC) \
	-DO2SCL_HDF5_COMP
LFFLAGS = -O3
LMPICXX = mpic++

else
ifeq ($(USER),awsteiner)

# On Andrew's laptop
LFC = mpif90
LCXX = $(MPI_CXX)
LIBS = -L$(O2SCL_LIB) -L$(GSL_LIB) -L$(HDF5_LIB) \
	-lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl -lhdf5 -lgsl
LCFLAGS = -O3 -std=c++11 -I$(O2SCL_INC) \
	-I$(EIGEN_INC) -I$(GSL_INC) -Wno-ignored-attributes \
	-Wno-deprecated-declarations -I$(HDF5_INC) \
	-Wshadow
LFFLAGS = -O3
LMPICXX = mpic++

else

# Default settings
LFC = $(FC)
LCXX = $(CXX)
LIBS = -L/usr/local/lib -lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl \
	-lhdf5 -lgsl
LCFLAGS = -O3 -std=c++11 
LFFLAGS = -O3
LMPICXX = mpic++

endif
endif
endif

# ----------------------------------------------------------------
# Main targets
# ----------------------------------------------------------------

empty:

doc: empty
	cd doc; cp ~/o2scl/doc/o2scl/o2scl.tag .
	cd doc; cp ~/o2scl/doc/o2scl/part/o2scl_part.tag .
	cd doc; cp ~/o2scl/doc/o2scl/eos/o2scl_eos.tag .
	cd doc; doxygen doxyfile
	cd sphinx; make html
#	cp -r sphinx/build/html/* $(HOME)/data/ecn

eos: eos.o main.o
	$(LCXX) $(LCFLAGS) -o eos eos.o main.o $(LIBS) \
		-lreadline

eos.o: eos.cpp virial_solver.h eos.h
	$(LCXX) $(LCFLAGS) -o eos.o -c eos.cpp

main.o: main.cpp virial_solver.h eos.h
	$(LCXX) $(LCFLAGS) -o main.o -c main.cpp

clean:
	rm -f *.o eos
