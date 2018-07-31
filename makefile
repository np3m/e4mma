# ----------------------------------------------------------------
# Various user-specific settings
# ----------------------------------------------------------------

# LIBS is the list of libraries
# LCXX is the local C++ compiler
# LCFLAGS are the local C++ compiler flags

# Default settings
LCXX = $(CXX)
LIBS = -L/usr/local/lib -lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl \
	-lhdf5 -lgsl
LCFLAGS = -O3 -std=c++11 -DNO_MPI

# ----------------------------------------------------------------
# UTK-specific settings
# ----------------------------------------------------------------

ifdef UTKNA_MAKEFILE

include $(UTKNA_MAKEFILE)

# UTK configuration
LIBS = $(UTKNA_O2SCL_LIBS)
LCXX = $(UTKNA_CXX) 
LCFLAGS = $(UTKNA_O2SCL_INCS) $(UTKNA_CFLAGS) -DNO_MPI

else

ifeq ($(MACHINE),isospin)

# On isospin for Xingfu
LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
	-lo2scl_eos -lo2scl_part -lo2scl_hdf -lo2scl -lhdf5 -lgsl
LCXX = g++ 
LCFLAGS = -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include \
	-I/usr/include/eigen3 -Wno-deprecated-declarations \
	-O3 -std=c++11 -DNO_MPI -Wshadow -DO2SCL_HDF5_COMP

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

sync-doc:
	sudo cp -r sphinx/build/html/* $(STATIC_DOC_DIR)/eos

eos: eos.o main.o
	$(LCXX) $(LCFLAGS) -o eos eos.o main.o $(LIBS) \
		-lreadline

eos.o: eos.cpp virial_solver.h eos.h
	$(LCXX) $(LCFLAGS) -o eos.o -c eos.cpp

main.o: main.cpp virial_solver.h eos.h
	$(LCXX) $(LCFLAGS) -o main.o -c main.cpp

clean:
	rm -f *.o eos
