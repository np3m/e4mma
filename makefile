help:
	@echo "Makefile targets:"
	@echo "─────────────────────────────────────────────────────────"
	@echo "help:             This makefile documentation."
	@echo "eos:              Main homogeneous matter executable."
	@echo "eos_nompi:        Non-parallel version of eos."
	@echo "eos_nuclei:       Main heterogeneous matter executable."
	@echo "eos_nuclei_nompi: Non-parallel version of eos_nuclei."
	@echo "open-doc:         Open local documentation in browser."
	@echo "doc:              Generate documentation."
	@echo "sync-doc:         Upload documentation to web."
	@echo "web-doc:          Open online documentation in browser."

# ----------------------------------------------------------------
# Various user-specific settings
# ----------------------------------------------------------------

# LIBS is the list of libraries
# LCXX is the local C++ compiler
# LCFLAGS are the local C++ compiler flags

# Default settings
LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
	-lo2scl -lhdf5 -lgsl \
	-lreadline
# PLIBS = -L/usr/lib/x86_64-linux-gnu/ 
LCXX = g++
LMPI_CXX = mpic++
LCFLAGS = -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include \
	-DNO_MPI -DNO_OPENMP -DO2SCL_NO_BOOST_MULTIPRECISION 
LFFLAGS = -O3
LMPI_CFLAGS = -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include \
	-DO2SCL_MPI -DO2SCL_OPENMP \
	-fopenmp -DTEMP_UPDATES -DO2SCL_NO_BOOST_MULTIPRECISION
	
COMMENT = "default"
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

eos_interp.o: eos_interp.cpp
	$(LMPI_CXX) $(LMPI_CFLAGS) \
		-o eos_interp.o -c eos_interp.cpp

eos_neutrino.o: eos_neutrino.cpp
	$(LMPI_CXX) $(LMPI_CFLAGS) \
		-o eos_neutrino.o -c eos_neutrino.cpp

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

fore.o: fore.cpp fore.h
	$(LMPI_CXX) $(LMPI_CFLAGS) -o fore.o -c fore.cpp 

fore_nompi.o: fore.cpp fore.h
	$(LCXX) $(LCFLAGS) -o fore_nompi.o -c fore.cpp 

fore_test.o: fore_test.cpp fore.h
	$(LMPI_CXX) $(LMPI_CFLAGS) -o fore_test.o -c fore_test.cpp

eos_nuclei: eos.o main.o eos_nuclei.o fore.o eos_had_skyrme_ext.o eos_interp.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		eos_neutrino.o
	$(LMPI_CXX) $(LMPI_CFLAGS) -I/usr/local/include/yaml-cpp/ -o eos_nuclei eos.o main.o \
		eos_nuclei.o fore.o eos_had_skyrme_ext.o eos_interp.o \
		neutrino/Couplings.o neutrino/FluidState.o eos_neutrino.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		$(LIBS) -lyaml-cpp

fore_test: fore_test.o fore.o 
	$(LMPI_CXX) $(LMPI_CFLAGS)-o fore_test fore_test.o fore.o $(LIBS)

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

eos_interp_nompi.o: eos_interp.cpp
	$(LMPI_CXX) $(LMPI_CFLAGS) \
		-o eos_interp_nompi.o -c eos_interp.cpp

eos_neutrino_nompi.o: eos_neutrino.cpp
	$(LMPI_CXX) $(LMPI_CFLAGS) \
		-o eos_neutrino_nompi.o -c eos_neutrino.cpp

eos_had_skyrme_ext_nompi.o: eos_had_skyrme_ext.cpp \
		eos_had_skyrme_ext.h
	$(LCXX) $(LCFLAGS) -o eos_had_skyrme_ext_nompi.o \
		-c eos_had_skyrme_ext.cpp

main_nompi.o: main.cpp eos.h 
	$(LCXX) $(LCFLAGS) -o main_nompi.o -c main.cpp 

main_eos_nompi.o: main_eos.cpp eos.h 
	$(LCXX) $(LCFLAGS) -o main_eos_nompi.o -c main_eos.cpp 

eos_nuclei_nompi: eos_nompi.o main_nompi.o eos_nuclei_nompi.o \
		eos_had_skyrme_ext_nompi.o eos_interp_nompi.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		eos_neutrino_nompi.o
	$(LCXX) $(LCFLAGS) -o eos_nuclei_nompi eos_nompi.o main_nompi.o \
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
	$(LCXX) $(LCFLAGS) -o enn eos_nompi.o main_nompi.o \
		eos_interp_nompi.o neutrino/Couplings.o \
		neutrino/FluidState.o neutrino/FunctionIntegrator.o \
		neutrino/Polarization.o eos_neutrino_nompi.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		eos_nuclei_nompi.o eos_had_skyrme_ext_nompi.o $(LIBS)

eos_nompi: eos_nompi.o main_eos_nompi.o \
		eos_had_skyrme_ext_nompi.o 
	$(LCXX) $(LCFLAGS) -o eos_nompi eos_nompi.o \
		main_eos_nompi.o eos_had_skyrme_ext_nompi.o $(LIBS) \
		-lreadline

ymltest: 
	$(LCXX) -I/usr/local/include -L/usr/local/lib -lyaml-cpp -o ymltest ymltest.cpp 

test: 
	$(LMPI_CXX) $(LMPI_CFLAGS) fore.o -o test test.cpp $(LIBS)
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

BROWSER = 
UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        BROWSER += xdg-open
    endif
    ifeq ($(UNAME_S),Darwin)
        BROWSER += open
    endif

web-doc:
	$(BROWSER) https://neutronstars.utk.edu/code/eos/

open-doc:
	$(BROWSER) doc/build/html/index.html

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

mn-test:
	enn \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 9999 \
		-set verbose 3 \
		-load ~/data/eos/final/fid_6_30_21.o2 \
		-mcarlo-neutron mn_test.o2

mbnew:
	./eos_nuclei \
		-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-set recompute 1 \
		-point-nuclei 0.16 0.465 0.1 

mbpi:
	./eos_nuclei \
	-set select_cs2_test 0 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-set inc_hrg false \
		-load data/fid_3_5_22.o2 \
		-hrg-load ./pdg_uh_nonp.dat \
		-set recompute 1 \
		-point-nuclei 0.16 0.5 30 

enn_fid_lep:
	./eos_nuclei \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load data/fid_3_5_22.o2 \
		-muses-table create

enn_fid_nolep:
	./eos_nuclei \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 3 \
		-load data/fid_nolep_noderiv_3_4_22.o2 \
		-muses-table create

mbmuses:
	./eos_nuclei \
		-select-model $(P_FIDUCIAL) \
		-set recompute 1 \
		-create-new-table create
		
-include makefile.aws

yml_gen:
	python3 yaml_generator.py \
	--select_model P_FIDUCIAL \
	--a_virial 10.0 --b_virial 10.0 \
	--extend_frdm 0 \
	--fd_A_max 600 --max_ratio 7.0 \
	--fixed_dist_alg 1999 \
	--function_verbose 0 \
	--load data/fid_3_5_22.o2 \
	--output_format HDF5


# Read YAML parameter using yq	
# Define YAML file
#CONFIG_FILE := api/input/config.yaml

# Extract keys from the YAML file excluding 'set' (assuming 'set' is an object)
#ALL_KEYS := $(shell yq eval 'keys | .[]' $(CONFIG_FILE) | grep -v '^set$$')

# Generate variables for each parameter
#$(foreach key,$(ALL_KEYS),$(eval $(key) := $(shell yq eval '.$(key)' $(CONFIG_FILE))))

# Extract keys from the 'set' section
#SET_KEYS := $(shell yq eval '.set | keys | .[]' $(CONFIG_FILE))

# Generate variables for each 'set' parameter
#$(foreach key,$(SET_KEYS),$(eval set_$(key) := $(shell yq eval '.set.$(key)' $(CONFIG_FILE))))

# Default target

# enn_fid_lep target
#enn_fid_lep_yaml:
#	./eos_nuclei \
		-select-model $($(select_model)) \
		$(if $(set_a_virial),-set a_virial $(set_a_virial)) \
		$(if $(set_b_virial),-set b_virial $(set_b_virial)) \
		$(if $(set_extend_frdm),-set extend_frdm $(set_extend_frdm)) \
		$(if $(set_fd_A_max),-set fd_A_max $(set_fd_A_max)) \
		$(if $(set_max_ratio),-set max_ratio $(set_max_ratio)) \
		$(if $(set_fixed_dist_alg),-set fixed_dist_alg $(set_fixed_dist_alg)) \
		$(if $(set_function_verbose),-set function_verbose $(set_function_verbose)) \
		-load $(load) \
		-muses-table create \

#	cp utk_eos.csv api/output/utk_eos.csv

# muses plots

# nB vs T den-plot for A at Ye=0.4
plot1:
	o2graph -read data/fid_3_5_22.o2 mun -set colbar 1 \
	-to-table3d 0 1 slice 0.1 -den-plot slice pcm=True \
	-xtitle "$$ n_B~(\mathrm{fm}^{-3}) $$" -ytitle "$$ Y_e $$" \
	-show

# nB vs Ye den-plot for mup at T=0.1
plot2:
	o2graph -read data/fid_3_5_22.o2 mup -set colbar 1 \
	-to-table3d 0 1 slice 1 -den-plot slice pcm=True \
	-xtitle "$$ n_B~(\mathrm{fm}^{-3}) $$" -ytitle "$$ Y_e $$" \
	-show

# nB vs Ye den-plot for P at T=0.1
plot3:
	o2graph -read data/fid_3_5_22.o2 P -set colbar 1 \
	-to-table3d 0 1 slice 0.1 -den-plot slice pcm=True \
	-xtitle "$$ n_B~(\mathrm{fm}^{-3}) $$" -ytitle "$$ Y_e $$" \
	-show

# nB vs Ye den-plot for E at T=0.1
plot4:
	o2graph -read data/fid_3_5_22.o2 E -set colbar 1 \
	-to-table3d 0 1 slice 0.1 -den-plot slice pcm=True \
	-xtitle "$$ n_B~(\mathrm{fm}^{-3}) $$" -ytitle "$$ Y_e $$" \
	-show

plot5:
	o2graph \
		-read data/fid_3_5_22.o2 mun -to-table 0 nB mun 0.5 0.1 -plot nB mun \
		-read data/fid_3_5_22.o2 mup -to-table 0 nB mup 0.5 0.1 -plot nB mup \
	-xtitle "$$ n_B~(\mathrm{fm}^{-3}) $$" -ytitle "$$ \mu_N~(\mathrm{MeV}) $$" \
	-show

plot6:
	o2graph \
		-read data/fid_3_5_22.o2 P -to-table 0 nB P 0.5 0.1 -plot nB P \
		-read data/fid_3_5_22.o2 S -to-table 0 nB S 0.5 0.1 -plot nB S \
	-xtitle "$$ n_B~(\mathrm{fm}^{-3}) $$" -ytitle "$$ P, S~(\mathrm{MeV/fm^{-3}}) $$" \
	-save ps.png -show

plot7:
	o2graph \
		-read data/fid_3_5_22.o2 P -to-table 2 T P 0.16 0.5 -plot T P \
		-read data/fid_3_5_22.o2 S -to-table 2 T S 0.16 0.5 -plot T S \
	-xtitle "$$ T~(\mathrm{MeV}) $$" -ytitle "$$ P, S~(\mathrm{MeV/fm^{-3}}) $$" \
	-save psvsT.png -show

plot8:
	o2graph \
		-read data/fid_3_5_22.o2 mun -to-table 1 Ye mun 0.16 0.1 -plot Ye mun \
		-read data/fid_3_5_22.o2 mup -to-table 1 Ye mup 0.16 0.1 -plot Ye mup \
		-read data/fid_3_5_22.o2 mue -to-table 1 Ye mue 0.16 0.1 -plot Ye mue \
		-create table Ye "grid:1.0e-2,0.7,(0.7-1.0e-2)/70" -function "mun(Ye)-mup(Ye)-mue(Ye)" mutot -internal data/mutot.o2 \
		-read data/mutot.o2 -plot Ye mutot \
	-xtitle "$$ Ye $$" -ytitle "$$ \mu_n, \mu_p, \mu_e S~(\mathrm{MeV}) $$" \
	 -show
