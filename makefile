help:
	@echo Help

# ----------------------------------------------------------------
# Various user-specific settings
# ----------------------------------------------------------------

# LIBS is the list of libraries
# LCXX is the local C++ compiler
# LCFLAGS are the local C++ compiler flags

# Default settings
LCXX = $(CXX)
LMPI_CXX = $(MPI_CXX)
LIBS = -L/usr/local/lib -lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl \
        -lhdf5 -lgsl -lreadline
LMPI_CFLAGS = -O3 -std=c++11 -DTEMP_UPDATES -DO2SCL_MPI \
	-DO2SCL_OPENMP -fopenmp
LCFLAGS = -O3 -std=c++11 -DNO_MPI -DTEMP_UPDATES \
	-DO2SCL_OPENMP -fopenmp

# ----------------------------------------------------------------
# UTK-specific settings
# ----------------------------------------------------------------

ifdef UTKNA_MAKEFILE

include $(UTKNA_MAKEFILE)

# UTK configuration
LIBS = $(UTKNA_O2SCL_LIBS)
LCXX = $(UTKNA_CXX) 
LMPI_CXX = $(UTKNA_MPI_CXX)
EOS_DIR = $(UTKNA_EOS_DIR)
LCFLAGS = -ggdb $(UTKNA_O2SCL_INCS) $(UTKNA_CFLAGS) -DNO_MPI -DTEMP_UPDATES \
        -I$(EOS_DIR) $(UTKNA_OPENMP_FLAGS)
LMPI_CFLAGS = -ggdb $(UTKNA_O2SCL_INCS) $(UTKNA_MPI_CFLAGS) -DTEMP_UPDATES \
        -I$(EOS_DIR) -DO2SCL_MPI $(UTKNA_OPENMP_FLAGS)

endif

# ----------------------------------------------------------------
# Main
# ----------------------------------------------------------------

eos.o: eos.cpp virial_solver.h eos.h
	$(LMPI_CXX) $(LMPI_CFLAGS) -o eos.o -c eos.cpp

eos_nuclei.o: eos_nuclei.cpp virial_solver.h eos_nuclei.h
	$(LMPI_CXX) $(LMPI_CFLAGS) -o eos_nuclei.o -c eos_nuclei.cpp

eos_had_skyrme_ext.o: eos_had_skyrme_ext.cpp virial_solver.h \
		eos_had_skyrme_ext.h
	$(LMPI_CXX) $(LMPI_CFLAGS) -o eos_had_skyrme_ext.o \
		-c eos_had_skyrme_ext.cpp

main.o: main.cpp virial_solver.h eos.h 
	$(LMPI_CXX) $(LMPI_CFLAGS) -o main.o -c main.cpp 

zidu/Couplings.o: zidu/Couplings.cpp zidu/Couplings.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o zidu/Couplings.o \
		-c zidu/Couplings.cpp

zidu/FluidState.o: zidu/FluidState.cpp zidu/FluidState.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o zidu/FluidState.o \
		-c zidu/FluidState.cpp

zidu/FunctionIntegrator.o: zidu/FunctionIntegrator.cpp \
	zidu/FunctionIntegrator.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o zidu/FunctionIntegrator.o \
	-c zidu/FunctionIntegrator.cpp

zidu/Polarization.o: zidu/Polarization.cpp zidu/Polarization.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o zidu/Polarization.o \
		-c zidu/Polarization.cpp

zidu/PolarizationNonRelv2Apr8.o: zidu/PolarizationNonRelv2Apr8.cpp 
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL \
		-o zidu/PolarizationNonRelv2Apr8.o \
		-c zidu/PolarizationNonRelv2Apr8.cpp

zidu/jacobi_rule.o: zidu/jacobi_rule.cpp zidu/jacobi_rule.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o zidu/jacobi_rule.o \
		-c zidu/jacobi_rule.cpp

eos_nuclei: eos.o main.o eos_nuclei.o eos_had_skyrme_ext.o \
		virial_solver_deriv.h \
		zidu/Couplings.o zidu/FluidState.o \
		zidu/FunctionIntegrator.o zidu/Polarization.o \
		zidu/PolarizationNonRelv2Apr8.o zidu/jacobi_rule.o 
	$(LMPI_CXX) $(LMPI_CFLAGS) -o eos_nuclei eos.o main.o \
		eos_nuclei.o eos_had_skyrme_ext.o \
		zidu/Couplings.o zidu/FluidState.o \
		zidu/FunctionIntegrator.o zidu/Polarization.o \
		zidu/PolarizationNonRelv2Apr8.o zidu/jacobi_rule.o \
		$(LIBS) -lreadline

# ----------------------------------------------------------------
# Version without MPI
# ----------------------------------------------------------------

zidu/Couplings_nompi.o: zidu/Couplings.cpp zidu/Couplings.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o zidu/Couplings_nompi.o \
		-c zidu/Couplings.cpp

zidu/FluidState_nompi.o: zidu/FluidState.cpp zidu/FluidState.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o zidu/FluidState_nompi.o \
		-c zidu/FluidState.cpp

zidu/FunctionIntegrator_nompi.o: zidu/FunctionIntegrator.cpp \
	zidu/FunctionIntegrator.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o zidu/FunctionIntegrator_nompi.o \
	-c zidu/FunctionIntegrator.cpp

zidu/Polarization_nompi.o: zidu/Polarization.cpp zidu/Polarization.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o zidu/Polarization_nompi.o \
		-c zidu/Polarization.cpp

zidu/PolarizationNonRelv2Apr8_nompi.o: zidu/PolarizationNonRelv2Apr8.cpp 
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL \
		-o zidu/PolarizationNonRelv2Apr8_nompi.o \
		-c zidu/PolarizationNonRelv2Apr8.cpp

zidu/jacobi_rule_nompi.o: zidu/jacobi_rule.cpp zidu/jacobi_rule.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o zidu/jacobi_rule_nompi.o \
		-c zidu/jacobi_rule.cpp

eos_nompi.o: eos.cpp virial_solver.h eos.h
	$(LCXX) $(LCFLAGS) -o eos_nompi.o -c eos.cpp

eos_nuclei_nompi.o: eos_nuclei.cpp virial_solver.h eos_nuclei.h
	$(LCXX) $(LCFLAGS) -o eos_nuclei_nompi.o -c eos_nuclei.cpp

eos_had_skyrme_ext_nompi.o: eos_had_skyrme_ext.cpp virial_solver.h \
		eos_had_skyrme_ext.h
	$(LCXX) $(LCFLAGS) -o eos_had_skyrme_ext_nompi.o \
		-c eos_had_skyrme_ext.cpp

main_nompi.o: main.cpp virial_solver.h eos.h 
	$(LCXX) $(LCFLAGS) -o main_nompi.o -c main.cpp 

main_eos_nompi.o: main_eos.cpp virial_solver.h eos.h 
	$(LCXX) $(LCFLAGS) -o main_eos_nompi.o -c main_eos.cpp 

eos_nuclei_nompi: eos_nompi.o main_nompi.o eos_nuclei_nompi.o \
		eos_had_skyrme_ext_nompi.o virial_solver_deriv.h
	$(LCXX) $(LCFLAGS) -o eos_nuclei_nompi eos_nompi.o main_nompi.o \
		 eos_nuclei_nompi.o eos_had_skyrme_ext_nompi.o $(LIBS) \
		-lreadline

# A shorthand alias for eos_nuclei_nompi
enn: eos_nompi.o main_nompi.o eos_nuclei_nompi.o \
		eos_had_skyrme_ext_nompi.o eos_nompi.o \
		zidu/Couplings.o zidu/FluidState.o \
		zidu/FunctionIntegrator.o zidu/Polarization.o \
		zidu/PolarizationNonRelv2Apr8.o zidu/jacobi_rule.o 
	$(LCXX) $(LCFLAGS) -o enn eos_nompi.o main_nompi.o \
		zidu/Couplings.o zidu/FluidState.o \
		zidu/FunctionIntegrator.o zidu/Polarization.o \
		zidu/PolarizationNonRelv2Apr8.o zidu/jacobi_rule.o \
		eos_nuclei_nompi.o eos_had_skyrme_ext_nompi.o $(LIBS) \
		-lreadline

eos_nompi: eos_nompi.o main_eos_nompi.o \
		eos_had_skyrme_ext_nompi.o virial_solver_deriv.h
	$(LCXX) $(LCFLAGS) -o eos_nompi eos_nompi.o \
		main_eos_nompi.o eos_had_skyrme_ext_nompi.o $(LIBS) \
		-lreadline

# ----------------------------------------------------------------
# Other targets
# ----------------------------------------------------------------

empty:

doc: empty
	cd doc; cp ~/o2scl/doc/o2scl/o2scl.tag .
	cd doc; cp ~/o2scl/doc/o2scl/part/o2scl_part.tag .
	cd doc; cp ~/o2scl/doc/o2scl/eos/o2scl_eos.tag .
	cd doc; cp ~/o2scl/doc/o2scl/sphinx/build/html/objects.inv \
		o2scl_objects.inv
	cd doc; cp ~/o2scl/doc/o2scl/part/sphinx/build/html/objects.inv \
		o2scl_part_objects.inv
	cd doc; cp ~/o2scl/doc/o2scl/eos/sphinx/build/html/objects.inv \
		o2scl_eos_objects.inv
	cd doc; doxygen doxyfile
	cd doc; make html

sync-doc:
	rsync -Cavzu doc/build/html/* $(STATIC_DOC_DIR)/eos

test-sync:
	rsync -Cavzun doc/build/html/* $(STATIC_DOC_DIR)/eos

clean:
	rm -f *.o eos_nuclei eos_nuclei_nompi eos eos_nompi zidu/*.o

# ----------------------------------------------------------------
# New EOS parameter sets
# ----------------------------------------------------------------

P_FIDUCIAL = 470 738 0.5 13.0 62.4 32.8 0.9
P_LARGE_MMAX = 783 738 0.5 13.0 62.4 32.8 0.9
P_SMALL_R = 214 738 0.5 13.0 62.4 32.8 0.9
P_SMALLER_R = 256 738 0.5 13.0 62.4 32.8 0.9
P_LARGE_R = 0 738 0.5 13.0 62.4 32.8 0.9
P_SMALL_SL = 470 738 0.5 13.0 23.7 29.5 0.9
P_LARGE_SL = 470 738 0.5 13.0 100.0 36.0 0.9

# ----------------------------------------------------------------

mbt:
	enn \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load ~/data/eos/final/fid_6_30_21.o2 \
		-mcarlo-beta 0.14 10.0 mb_temp.o2

nt:
	enn -test-neutrino
