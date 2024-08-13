help:
	@echo "Makefile targets:"
	@echo "─────────────────────────────────────────────────────────"
	@echo "help:             This makefile documentation."
	@echo "eos:              Main homogeneous matter executable."
	@echo "eos_nompi:        Non-parallel version of eos."
	@echo "eos_nuclei:       Main heterogeneous matter executable."
	@echo "eos_nuclei_nompi: Non-parallel version of eos_nuclei."
	@echo "enn:              Alias for 'eos_nuclei_nompi'."
	@echo "web-doc:          Open online documentation in browser."
	@echo "clean:            Delete executables and .o files."
	@echo "check:            Check the code."
	@echo ""
	@echo "Developer targets:"
	@echo "─────────────────────────────────────────────────────────"
	@echo "doc:              Generate documentation."
	@echo "open-doc:         Open local documentation in browser."
	@echo "sync-doc:         Upload documentation to web."
	@echo "test-sync:         Upload documentation to web."

# ----------------------------------------------------------------
# Various user-specific settings
# ----------------------------------------------------------------

# LIBS is the list of libraries
# LCXX is the local C++ compiler
# LCFLAGS are the local C++ compiler flags

# Default settings
LCXX = $(CXX)
LMPI_CXX = $(MPI_CXX)
LIBS = -L/usr/local/lib -lo2scl -lhdf5 -lgsl -lreadline $(LDFLAGS) 
LMPI_CFLAGS = -O3 -std=c++11 -DO2SCL_MPI -fopenmp $(CFLAGS) $(MPI_CFLAGS)
LCFLAGS = -O3 -std=c++11 -DNO_MPI -fopenmp $(CFLAGS)

# ----------------------------------------------------------------
# UTK-specific settings
# ----------------------------------------------------------------

ifdef UTKNA_MAKEFILE

include $(UTKNA_MAKEFILE)

# UTK configuration

LIBS = $(UTKNA_O2SCL_LIBS) $(UTKNA_PYTHON_LDFLAGS)

LCXX = $(UTKNA_CXX) 
LMPI_CXX = $(UTKNA_MPI_CXX)
LCFLAGS = $(UTKNA_O2SCL_INCS) $(UTKNA_CFLAGS) -DNO_MPI \
        $(UTKNA_OPENMP_FLAGS)
LMPI_CFLAGS = $(UTKNA_O2SCL_INCS) $(UTKNA_CFLAGS) \
	$(UTKNA_OPENMP_FLAGS) $(UTKNA_MPI_CFLAGS)

endif

# ----------------------------------------------------------------
# Main
# ----------------------------------------------------------------

eos.o: eos.cpp eos.h
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos.o -c eos.cpp

eos_nuclei.o: eos_nuclei.cpp eos_nuclei.h
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos_nuclei.o -c eos_nuclei.cpp

eos_interp.o: eos_interp.cpp
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos_interp.o -c eos_interp.cpp

eos_neutrino.o: eos_neutrino.cpp
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos_neutrino.o -c eos_neutrino.cpp

eos_had_skyrme_ext.o: eos_had_skyrme_ext.cpp \
		eos_had_skyrme_ext.h
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos_had_skyrme_ext.o \
		-c eos_had_skyrme_ext.cpp

main.o: main.cpp eos.h 
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o main.o -c main.cpp 

neutrino/Couplings.o: neutrino/Couplings.cpp neutrino/Couplings.hpp
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-DNUOPAC_HAS_GSL -o neutrino/Couplings.o \
		-c neutrino/Couplings.cpp

neutrino/FluidState.o: neutrino/FluidState.cpp neutrino/FluidState.hpp
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-DNUOPAC_HAS_GSL -o neutrino/FluidState.o \
		-c neutrino/FluidState.cpp

neutrino/FunctionIntegrator.o: neutrino/FunctionIntegrator.cpp \
	neutrino/FunctionIntegrator.hpp
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-DNUOPAC_HAS_GSL -o \
		neutrino/FunctionIntegrator.o \
		-c neutrino/FunctionIntegrator.cpp

neutrino/Polarization.o: neutrino/Polarization.cpp \
		neutrino/Polarization.hpp
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-DNUOPAC_HAS_GSL -o neutrino/Polarization.o \
		-c neutrino/Polarization.cpp

neutrino/PolarizationNonRelv2Apr8.o: \
		neutrino/PolarizationNonRelv2Apr8.cpp 
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-DNUOPAC_HAS_GSL \
		-o neutrino/PolarizationNonRelv2Apr8.o \
		-c neutrino/PolarizationNonRelv2Apr8.cpp

neutrino/jacobi_rule.o: neutrino/jacobi_rule.cpp neutrino/jacobi_rule.hpp
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-DNUOPAC_HAS_GSL -o neutrino/jacobi_rule.o \
		-c neutrino/jacobi_rule.cpp

eos_nuclei: eos.o main.o eos_nuclei.o eos_had_skyrme_ext.o eos_interp.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		eos_neutrino.o
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos_nuclei eos.o main.o \
		eos_nuclei.o eos_had_skyrme_ext.o eos_interp.o \
		neutrino/Couplings.o neutrino/FluidState.o eos_neutrino.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		$(LIBS) -lreadline

main_eos.o: main_eos.cpp eos.h 
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o main_eos.o -c main_eos.cpp 

eos: eos.o main_eos.o \
		eos_had_skyrme_ext.o 
	$(LMPI_CXX) $(LMPI_CFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos eos.o \
		main_eos.o eos_had_skyrme_ext.o $(LIBS) \
		-lreadline

# ----------------------------------------------------------------
# Main with multiprecision
# ----------------------------------------------------------------

eos_mp.o: eos.cpp eos.h
	$(LMPI_CXX) $(LMPI_CFLAGS) \
		-o eos_mp.o -c eos.cpp

eos_nuclei_mp.o: eos_nuclei.cpp eos_nuclei.h
	$(LMPI_CXX) $(LMPI_CFLAGS) \
		-o eos_nuclei_mp.o -c eos_nuclei.cpp

eos_interp_mp.o: eos_interp.cpp
	$(LMPI_CXX) $(LMPI_CFLAGS) \
		-o eos_interp_mp.o -c eos_interp.cpp

eos_neutrino_mp.o: eos_neutrino.cpp
	$(LMPI_CXX) $(LMPI_CFLAGS) \
		-o eos_neutrino_mp.o -c eos_neutrino.cpp

eos_had_skyrme_ext_mp.o: eos_had_skyrme_ext.cpp \
		eos_had_skyrme_ext.h
	$(LMPI_CXX) $(LMPI_CFLAGS) -o eos_had_skyrme_ext_mp.o \
		-c eos_had_skyrme_ext.cpp

main_mp.o: main.cpp eos.h 
	$(LMPI_CXX) $(LMPI_CFLAGS) -o main_mp.o -c main.cpp 

neutrino/Couplings_mp.o: neutrino/Couplings.cpp neutrino/Couplings.hpp
	$(LMPI_CXX) $(LMPI_CFLAGS) -DNUOPAC_HAS_GSL \
		-o neutrino/Couplings_mp.o \
		-c neutrino/Couplings.cpp

neutrino/FluidState_mp.o: neutrino/FluidState.cpp neutrino/FluidState.hpp
	$(LMPI_CXX) $(LMPI_CFLAGS) -DNUOPAC_HAS_GSL \
		-o neutrino/FluidState_mp.o \
		-c neutrino/FluidState.cpp

neutrino/FunctionIntegrator_mp.o: neutrino/FunctionIntegrator.cpp \
	neutrino/FunctionIntegrator.hpp
	$(LMPI_CXX) $(LMPI_CFLAGS) -DNUOPAC_HAS_GSL -o \
	neutrino/FunctionIntegrator_mp.o \
	-c neutrino/FunctionIntegrator.cpp

neutrino/Polarization_mp.o: neutrino/Polarization.cpp \
		neutrino/Polarization.hpp
	$(LMPI_CXX) $(LMPI_CFLAGS) -DNUOPAC_HAS_GSL -o \
		neutrino/Polarization_mp.o \
		-c neutrino/Polarization.cpp

neutrino/PolarizationNonRelv2Apr8_mp.o: \
		neutrino/PolarizationNonRelv2Apr8.cpp 
	$(LMPI_CXX) $(LMPI_CFLAGS) -DNUOPAC_HAS_GSL \
		-o neutrino/PolarizationNonRelv2Apr8_mp.o \
		-c neutrino/PolarizationNonRelv2Apr8.cpp

neutrino/jacobi_rule_mp.o: neutrino/jacobi_rule.cpp \
		neutrino/jacobi_rule.hpp
	$(LMPI_CXX) $(LMPI_CFLAGS) -DNUOPAC_HAS_GSL -o \
		neutrino/jacobi_rule_mp.o \
		-c neutrino/jacobi_rule.cpp

eos_nuclei_mp: eos_mp.o main_mp.o eos_nuclei_mp.o \
		eos_had_skyrme_ext_mp.o eos_interp_mp.o \
		neutrino/Couplings_mp.o neutrino/FluidState_mp.o \
		neutrino/FunctionIntegrator_mp.o \
		neutrino/Polarization_mp.o \
		neutrino/PolarizationNonRelv2Apr8_mp.o \
		neutrino/jacobi_rule_mp.o \
		eos_neutrino_mp.o
	$(LMPI_CXX) $(LMPI_CFLAGS) -o eos_nuclei_mp eos_mp.o main_mp.o \
		eos_nuclei_mp.o eos_had_skyrme_ext_mp.o eos_interp_mp.o \
		neutrino/Couplings_mp.o neutrino/FluidState_mp.o \
		eos_neutrino_mp.o neutrino/FunctionIntegrator_mp.o \
		neutrino/Polarization_mp.o \
		neutrino/PolarizationNonRelv2Apr8_mp.o \
		neutrino/jacobi_rule_mp.o \
		$(LIBS) -lreadline

# ----------------------------------------------------------------
# Version without MPI
# ----------------------------------------------------------------

neutrino/Couplings_nompi.o: neutrino/Couplings.cpp neutrino/Couplings.hpp
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-DNUOPAC_HAS_GSL -o neutrino/Couplings_nompi.o \
		-c neutrino/Couplings.cpp

neutrino/FluidState_nompi.o: neutrino/FluidState.cpp neutrino/FluidState.hpp
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-DNUOPAC_HAS_GSL -o neutrino/FluidState_nompi.o \
		-c neutrino/FluidState.cpp

neutrino/FunctionIntegrator_nompi.o: neutrino/FunctionIntegrator.cpp \
		neutrino/FunctionIntegrator.hpp
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-DNUOPAC_HAS_GSL -o neutrino/FunctionIntegrator_nompi.o \
	-c neutrino/FunctionIntegrator.cpp

neutrino/Polarization_nompi.o: neutrino/Polarization.cpp \
		neutrino/Polarization.hpp
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-DNUOPAC_HAS_GSL -o neutrino/Polarization_nompi.o \
		-c neutrino/Polarization.cpp

neutrino/PolarizationNonRelv2Apr8_nompi.o: \
		neutrino/PolarizationNonRelv2Apr8.cpp 
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-DNUOPAC_HAS_GSL \
		-o neutrino/PolarizationNonRelv2Apr8_nompi.o \
		-c neutrino/PolarizationNonRelv2Apr8.cpp

neutrino/jacobi_rule_nompi.o: neutrino/jacobi_rule.cpp \
		neutrino/jacobi_rule.hpp
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-DNUOPAC_HAS_GSL -o neutrino/jacobi_rule_nompi.o \
		-c neutrino/jacobi_rule.cpp

eos_nompi.o: eos.cpp eos.h
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos_nompi.o -c eos.cpp

eos_nuclei_nompi.o: eos_nuclei.cpp eos_nuclei.h
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos_nuclei_nompi.o -c eos_nuclei.cpp

eos_interp_nompi.o: eos_interp.cpp
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos_interp_nompi.o -c eos_interp.cpp

eos_neutrino_nompi.o: eos_neutrino.cpp
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos_neutrino_nompi.o -c eos_neutrino.cpp

eos_had_skyrme_ext_nompi.o: eos_had_skyrme_ext.cpp \
		eos_had_skyrme_ext.h
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos_had_skyrme_ext_nompi.o \
		-c eos_had_skyrme_ext.cpp

main_nompi.o: main.cpp eos.h 
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o main_nompi.o -c main.cpp 

main_eos_nompi.o: main_eos.cpp eos.h 
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o main_eos_nompi.o -c main_eos.cpp 

eos_nuclei_nompi: eos_nompi.o main_nompi.o eos_nuclei_nompi.o \
		eos_had_skyrme_ext_nompi.o eos_interp_nompi.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o \
		neutrino/jacobi_rule.o eos_neutrino_nompi.o
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION \
		-o eos_nuclei_nompi eos_nompi.o main_nompi.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		 eos_nuclei_nompi.o eos_had_skyrme_ext_nompi.o $(LIBS) \
		-lreadline

# A shorthand alias for eos_nuclei_nompi
enn: eos_nompi.o main_nompi.o eos_nuclei_nompi.o eos_interp_nompi.o \
		eos_had_skyrme_ext_nompi.o eos_nompi.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		eos_neutrino_nompi.o
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION -o enn \
		eos_nompi.o main_nompi.o \
		eos_interp_nompi.o neutrino/Couplings.o \
		neutrino/FluidState.o neutrino/FunctionIntegrator.o \
		neutrino/Polarization.o eos_neutrino_nompi.o \
		neutrino/PolarizationNonRelv2Apr8.o \
		neutrino/jacobi_rule.o eos_nuclei_nompi.o \
		eos_had_skyrme_ext_nompi.o $(LIBS)

eos_nompi: eos_nompi.o main_eos_nompi.o \
		eos_had_skyrme_ext_nompi.o 
	$(LCXX) $(LCFLAGS) -DO2SCL_NO_BOOST_MULTIPRECISION -o eos_nompi \
		eos_nompi.o main_eos_nompi.o \
		eos_had_skyrme_ext_nompi.o $(LIBS) -lreadline

# ----------------------------------------------------------------
# Other targets
# ----------------------------------------------------------------

empty:

doc: empty
	cd doc; cp ~/o2scl/doc/o2scl/o2scl.tag .
	cd doc; cp ~/o2scl/doc/o2scl/html/objects.inv \
		o2scl_objects.inv
	cd doc; doxygen doxyfile
	cd doc; make html

doc-auto: enn
	cd doc; cp ~/o2scl/doc/o2scl/o2scl.tag .
	cd doc; cp ~/o2scl/doc/o2scl/html/objects.inv \
		o2scl_objects.inv
	cd doc; doxygen doxyfile
	enn -xml-to-o2
	cd doc; make html

BROWSER = 
UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        BROWSER += xdg-open
    endif
    ifeq ($(UNAME_S),Darwin)
        BROWSER += open
    endif

web-doc:
	$(BROWSER) https://np3m.org/code/e4mma

open-doc:
	$(BROWSER) doc/build/html/index.html

sync-doc:
	rsync -Cavzu doc/build/html/* ~/wcs/np3m/np3m.github.io/code/e4mma

test-sync:
	rsync -Cavzun doc/build/html/* ~/wcs/np3m/np3m.github.io/code/e4mma

clean:
	rm -f *.o eos_nuclei eos_nuclei_nompi eos eos_nompi enn neutrino/*.o

# ----------------------------------------------------------------
# EOS parameter sets from Du et al. (2022)
# ----------------------------------------------------------------

P_FIDUCIAL = 470 738 0.5 13.0 62.4 32.8 0.9
P_LARGE_MMAX = 783 738 0.5 13.0 62.4 32.8 0.9
P_SMALL_R = 214 738 0.5 13.0 62.4 32.8 0.9
P_SMALLER_R = 256 738 0.5 13.0 62.4 32.8 0.9
P_LARGE_R = 0 738 0.5 13.0 62.4 32.8 0.9
P_SMALL_SL = 470 738 0.5 13.0 23.7 29.5 0.9
P_LARGE_SL = 470 738 0.5 13.0 100.0 36.0 0.9

# ----------------------------------------------------------------
# Check target
# ----------------------------------------------------------------

check: eos eos_nuclei
	cd examples && ./A_point_nuc.scr
	cd examples && ./B_point.scr
	cd examples && ./C_help.scr

# This optional file, makefile.user, is an alternate place to store
# the user's makefile targets.
-include makefile.user
