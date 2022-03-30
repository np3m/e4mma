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
LCFLAGS = $(UTKNA_O2SCL_INCS) $(UTKNA_CFLAGS) -DNO_MPI \
        $(UTKNA_OPENMP_FLAGS) -DO2SCL_EIGEN \
	-DO2SCL_NEW_BOOST_INTEGRATION
LMPI_CFLAGS = $(UTKNA_O2SCL_INCS) $(UTKNA_CFLAGS) \
	$(UTKNA_OPENMP_FLAGS) $(UTKNA_MPI_CFLAGS) \
	-DO2SCL_NEW_BOOST_INTEGRATION

endif

# ----------------------------------------------------------------
# Main
# ----------------------------------------------------------------

eos.o: eos.cpp eos.h
	$(LMPI_CXX) $(LMPI_CFLAGS) \
		-o eos.o -c eos.cpp

eos_nuclei.o: eos_nuclei.cpp eos_nuclei.h
	$(LMPI_CXX) $(LMPI_CFLAGS) \
		-o eos_nuclei.o -c eos_nuclei.cpp

eos_had_skyrme_ext.o: eos_had_skyrme_ext.cpp \
		eos_had_skyrme_ext.h
	$(LMPI_CXX) $(LMPI_CFLAGS) -o eos_had_skyrme_ext.o \
		-c eos_had_skyrme_ext.cpp

main.o: main.cpp eos.h 
	$(LMPI_CXX) $(LMPI_CFLAGS) -o main.o -c main.cpp 

neutrino/Couplings.o: neutrino/Couplings.cpp neutrino/Couplings.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o neutrino/Couplings.o \
		-c neutrino/Couplings.cpp

neutrino/FluidState.o: neutrino/FluidState.cpp neutrino/FluidState.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o neutrino/FluidState.o \
		-c neutrino/FluidState.cpp

neutrino/FunctionIntegrator.o: neutrino/FunctionIntegrator.cpp \
	neutrino/FunctionIntegrator.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o neutrino/FunctionIntegrator.o \
	-c neutrino/FunctionIntegrator.cpp

neutrino/Polarization.o: neutrino/Polarization.cpp neutrino/Polarization.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o neutrino/Polarization.o \
		-c neutrino/Polarization.cpp

neutrino/PolarizationNonRelv2Apr8.o: neutrino/PolarizationNonRelv2Apr8.cpp 
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL \
		-o neutrino/PolarizationNonRelv2Apr8.o \
		-c neutrino/PolarizationNonRelv2Apr8.cpp

neutrino/jacobi_rule.o: neutrino/jacobi_rule.cpp neutrino/jacobi_rule.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o neutrino/jacobi_rule.o \
		-c neutrino/jacobi_rule.cpp

eos_nuclei: eos.o main.o eos_nuclei.o eos_had_skyrme_ext.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o 
	$(LMPI_CXX) $(LMPI_CFLAGS) -o eos_nuclei eos.o main.o \
		eos_nuclei.o eos_had_skyrme_ext.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		$(LIBS) -lreadline

# ----------------------------------------------------------------
# Version without MPI
# ----------------------------------------------------------------

neutrino/Couplings_nompi.o: neutrino/Couplings.cpp neutrino/Couplings.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o neutrino/Couplings_nompi.o \
		-c neutrino/Couplings.cpp

neutrino/FluidState_nompi.o: neutrino/FluidState.cpp neutrino/FluidState.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o neutrino/FluidState_nompi.o \
		-c neutrino/FluidState.cpp

neutrino/FunctionIntegrator_nompi.o: neutrino/FunctionIntegrator.cpp \
		neutrino/FunctionIntegrator.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o \
	neutrino/FunctionIntegrator_nompi.o \
	-c neutrino/FunctionIntegrator.cpp

neutrino/Polarization_nompi.o: neutrino/Polarization.cpp \
		neutrino/Polarization.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o neutrino/Polarization_nompi.o \
		-c neutrino/Polarization.cpp

neutrino/PolarizationNonRelv2Apr8_nompi.o: \
		neutrino/PolarizationNonRelv2Apr8.cpp 
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL \
		-o neutrino/PolarizationNonRelv2Apr8_nompi.o \
		-c neutrino/PolarizationNonRelv2Apr8.cpp

neutrino/jacobi_rule_nompi.o: neutrino/jacobi_rule.cpp neutrino/jacobi_rule.hpp
	$(LCXX) $(LCFLAGS) -DNUOPAC_HAS_GSL -o neutrino/jacobi_rule_nompi.o \
		-c neutrino/jacobi_rule.cpp

eos_nompi.o: eos.cpp eos.h
	$(LCXX) $(LCFLAGS) \
		-o eos_nompi.o -c eos.cpp

eos_nuclei_nompi.o: eos_nuclei.cpp eos_nuclei.h
	$(LCXX) $(LCFLAGS) \
		-o eos_nuclei_nompi.o -c eos_nuclei.cpp

eos_had_skyrme_ext_nompi.o: eos_had_skyrme_ext.cpp \
		eos_had_skyrme_ext.h
	$(LCXX) $(LCFLAGS) -o eos_had_skyrme_ext_nompi.o \
		-c eos_had_skyrme_ext.cpp

main_nompi.o: main.cpp eos.h 
	$(LCXX) $(LCFLAGS) -o main_nompi.o -c main.cpp 

main_eos_nompi.o: main_eos.cpp eos.h 
	$(LCXX) $(LCFLAGS) -o main_eos_nompi.o -c main_eos.cpp 

eos_nuclei_nompi: eos_nompi.o main_nompi.o eos_nuclei_nompi.o \
		eos_had_skyrme_ext_nompi.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o 
	$(LCXX) $(LCFLAGS) -o eos_nuclei_nompi eos_nompi.o main_nompi.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		 eos_nuclei_nompi.o eos_had_skyrme_ext_nompi.o $(LIBS) \
		-lreadline

# A shorthand alias for eos_nuclei_nompi
enn: eos_nompi.o main_nompi.o eos_nuclei_nompi.o \
		eos_had_skyrme_ext_nompi.o eos_nompi.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o 
	$(LCXX) $(LCFLAGS) -o enn eos_nompi.o main_nompi.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		eos_nuclei_nompi.o eos_had_skyrme_ext_nompi.o $(LIBS) \
		-lreadline

eos_nompi: eos_nompi.o main_eos_nompi.o \
		eos_had_skyrme_ext_nompi.o 
	$(LCXX) $(LCFLAGS) -o eos_nompi eos_nompi.o \
		main_eos_nompi.o eos_had_skyrme_ext_nompi.o $(LIBS) \
		-lreadline

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

open-doc:
	open doc/build/html/index.html

sync-doc:
	rsync -Cavzu doc/build/html/* $(STATIC_DOC_DIR)/eos

test-sync:
	rsync -Cavzun doc/build/html/* $(STATIC_DOC_DIR)/eos

clean:
	rm -f *.o eos_nuclei eos_nuclei_nompi eos eos_nompi neutrino/*.o

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

mbtest:
	enn \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load ~/data/eos/final/fid_6_30_21.o2 \
		-mcarlo-beta mb_temp1.o2

mbtestd:
	enn \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load ~/data/eos/final/fid_6_30_21.o2 \
		-mcarlo-beta mb_temp1.o2 100

mb1:
	eos_nuclei \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load ~/data/eos/final/fid_6_30_21.o2 \
		-mcarlo-beta mb1.o2 > mb1.out 2>&1 &

mb2:
	eos_nuclei \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load ~/data/eos/final/fid_6_30_21.o2 \
		-mcarlo-beta mb2.o2 > mb2.out 2>&1 &

mb3:
	eos_nuclei \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load ~/data/eos/final/fid_6_30_21.o2 \
		-mcarlo-beta mb3.o2 > mb3.out 2>&1 &

mb4:
	eos_nuclei \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load ~/data/eos/final/fid_6_30_21.o2 \
		-mcarlo-beta mb4.o2 > mb4.out 2>&1 &

mb1d:
	enn \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load ~/data/eos/final/fid_6_30_21.o2 \
		-mcarlo-beta mb1d.o2 100 > mb1d.out 2>&1 &

mb2d:
	enn \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load ~/data/eos/final/fid_6_30_21.o2 \
		-mcarlo-beta mb2d.o2 100 > mb2d.out 2>&1 &

mbt2:
	eos_nuclei \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load ~/data/eos/final/fid_6_30_21.o2 \
		-mcarlo-beta mb_temp2.o2 > mbt2.out 2>&1 &

mbt3:
	eos_nuclei \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load ~/data/eos/final/fid_6_30_21.o2 \
		-mcarlo-beta mb_temp3.o2 > mbt3.out 2>&1 &

nt:
	./enn -test-neutrino

nt2:
	./enn -test-neutrino > tn_new.out
	 head -n 170 tn.out | tail -n 10
	 head -n 40 tn_new.out | tail -n 12

yetest:
	o2graph -set logx 1 \
		-plotv "grid:1.0e-4,0.15,(0.15/1.0e-4)^(1/99),log" \
		"hdf5:~awsteiner/data/21/09/23/mb1d.o2:mb:0:Ye_best_*" \
		-create table x \
		"grid:1.0e-4,0.15,(0.15/1.0e-4)^(1/99),log" \
		-function "log10(x)*31.485+126" i \
		-function "0.05+0.28*exp(-i/24)" ye2 \
		-plot x ye2 \
		-show

imfps:
	o2graph -subplots 2 2 -set logx 1 -set logy 1 \
		-selax 0 \
		-plotv "grid:1.0e-4,0.15,(0.15/1.0e-4)^(1/99),log" \
		"hdf5:mb1d.o2:mb:0-4:cc_vec_imfp" \
		-ttext 0.3 0.8 "CC,vec" \
		-selax 1 \
		-plotv "grid:1.0e-4,0.15,(0.15/1.0e-4)^(1/99),log" \
		"hdf5:mb1d.o2:mb:0-4:cc_axvec_imfp" \
		-ttext 0.2 0.8 "CC,ax" \
		-selax 2 \
		-plotv "grid:1.0e-4,0.15,(0.15/1.0e-4)^(1/99),log" \
		"hdf5:mb1d.o2:mb:0-4:nc_vec_imfp" \
		-ttext 0.65 0.8 "NC,vec" \
		-selax 3 \
		-plotv "grid:1.0e-4,0.15,(0.15/1.0e-4)^(1/99),log" \
		"hdf5:mb1d.o2:mb:0-4:nc_axvec_imfp" \
		-ttext 0.3 0.8 "NC,ax" \
		-subadj "left=0.12,right=0.99,top=0.99,bottom=0.09,wspace=0.27,hspace=0.17" \
		-save imfps.pdf -show
