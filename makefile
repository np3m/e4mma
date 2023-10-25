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

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	LCFLAGS += -g -Wall -Wextra
	LMPI_CFLAGS += -g -Wall -Wextra
endif

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

eos_nuclei: eos.o main.o eos_nuclei.o eos_had_skyrme_ext.o eos_interp.o \
		neutrino/Couplings.o neutrino/FluidState.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		eos_neutrino.o
	$(LMPI_CXX) $(LMPI_CFLAGS) -o eos_nuclei eos.o main.o \
		eos_nuclei.o eos_had_skyrme_ext.o eos_interp.o \
		neutrino/Couplings.o neutrino/FluidState.o eos_neutrino.o \
		neutrino/FunctionIntegrator.o neutrino/Polarization.o \
		neutrino/PolarizationNonRelv2Apr8.o neutrino/jacobi_rule.o \
		$(LIBS) -lreadline

test_program: test.cpp
	$(LMPI_CXX) $(LMPI_CFLAGS) \
		-o test -c test.cpp
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

docker_clean:
	-sudo docker rm \
	`sudo docker ps --all | grep -i -v container | awk '{print $$1}'`
	-sudo docker rmi \
	`sudo docker images --all | grep -i -v container | awk '{print $$3}'`

docker_show:
	- sudo docker ps --all
	- sudo docker images --all

docker_build:
	sudo docker build - < docker \
		> docker.out 2>&1 &

test:
	./eos_nuclei -load /home/awsteiner/wcs/eos/fid_3_14_23.o2 \
	-interp-point 1.0e-10 0.4 10 2 /home/jbaut001/st.o2

test2:
	./eos_nuclei -load /home/awsteiner/wcs/eos/fid_3_14_23.o2 \
	-interp-point 0.06 0.67 10 2 /home/jbaut001/st.o2

test3:
	./eos_nuclei -load /home/awsteiner/wcs/eos/fid_3_14_23.o2 \
	-interp-point 4.178592e-12 5.000000e-01 10 2 /home/jbaut001/st.o2

test4:
	./eos_nuclei -load /home/awsteiner/wcs/eos/fid_3_14_23.o2 \
	-interp-point 0.000550846 0.05 10 2 /home/jbaut001/st.o2

test5:
	./eos_nuclei -load /home/awsteiner/wcs/eos/fid_3_14_23.o2 \
	-interp-point 0.000550846 0.05 0.1 2 /home/jbaut001/st.o2

point-test:
	./eos_nuclei \
		-set select_cs2_test 1 \
		-select-model $(P_FIDUCIAL) \
		-set a_virial 10 -set b_virial 10 \
		-set extend_frdm 0 \
		-set fd_A_max 600 -set max_ratio 7.0 \
		-set fixed_dist_alg 1999 \
		-set function_verbose 0 \
		-load /home/awsteiner/wcs/eos/fid_3_14_23.o2 \
		-point-nuclei 0.000502377 0.04 1.8602
makeimg:
	o2graph -read /home/jbaut001/eos/fid_7_21_23_multi_final.o2 Fint -set logx 1 -set xlo 0.01 -set xhi 0.1 -set ylo "(-1)" -set yhi 1 -to-table 0 nB Fint 0.04 1.77839 -function "Fint/10^3" Fs -plot nB Fs -xtitle "$$ n_B$$" -ytitle "$$ Fint$$" -ttext 0.1 0.45 "$$ Fint$$" -ttext 0.1 0.6 "$$ \frac{\partial Fint}{\partial n_B}$$" -ttext 0.67 0.75 "$$ \frac{\partial^2 Fint}{\partial n_B^2}$$" -ttext 0.5 0.95 " $$ Fint$$, $$ \frac{\partial Fint}{\partial n_B}$$, and $$ \frac{\partial Fint}{\partial n_B}$$ vs $$ n_B$$ " \
	       -deriv nB Fint Fp -function "Fp/10^3" Fs2 -plot nB Fs2\
	       -deriv2 nB Fint Fw -function "Fw/10^3" Fs3 -plot nB Fs3\
       	       -save fint_after.png


makeimga:
	o2graph -read /home/jbaut001/eos/fid_7_21_23_multi_final.o2 Fint -set logx 1 -set xlo 0.01 -set xhi 0.1 -set ylo "(-1)" -set yhi 1 -to-table 0 nB Fint 0.04 1.77839 -function "Fint/10^4" Fs -plot nB Fs -xtitle "$$ n_B$$" -ytitle "$$ Fint$$" -ttext 0.1 0.45 "$$ Fint$$" -ttext 0.1 0.55 "$$ \frac{\partial Fint}{\partial n_B}$$" -ttext 0.67 0.4 "$$ \frac{\partial^2 Fint}{\partial n_B^2}$$" -ttext 0.5 0.95 " $$ Fint$$, $$ \frac{\partial Fint}{\partial n_B}$$, and $$ \frac{\partial Fint}{\partial n_B}$$ vs $$ n_B$$ " \
	       -deriv nB Fint Fp -function "Fp/10^4" Fs2 -plot nB Fs2\
	       -deriv2 nB Fint Fw -function "Fw/10^4" Fs3 -plot nB Fs3\
       	       -save finta_after.png

makeimg2:
	o2graph -read /home/awsteiner/wcs/eos/fid_3_14_23.o2 Fint -set logx 1 -set xlo 0.01 -set xhi 0.1 --set ylo "(-50)" -set yhi 50 -to-table 0 nB Fint 0.05 1.8602 -plot nB Fint -xtitle "$$ n_B$$" -ytitle "$$ Fint$$" -ttext 0.5 1.021 " $$ Fint$$ vs $$ n_B$$ " \
       	       -save fint2_after.png


makeimg3:
	o2graph -read /home/awsteiner/wcs/eos/fid_3_14_23.o2 Fint -set logx 1 -to-table 0 nB Fint 0.05 1.8602 -function "Fint/10^23" Fs -plot nB Fs -xtitle "$$ n_B$$" -ytitle "$$ Fint$$" -ttext 0.5 1.021 " $$ Fint$$ vs $$ n_B$$ " \
       	       -save fint3_after.png

makeimg4:
	o2graph -read /home/awsteiner/wcs/eos/fid_3_14_23.o2 Fint -set logx 1 -set xlo 0.01 -set xhi 0.1 -set ylo "(-1)" -set yhi 1 -to-table 0 nB Fint 0.04 1.77839 -function "Fint/10^3" Fs -plot nB Fs -xtitle "$$ n_B$$" -ytitle "$$ Fint$$" -ttext 0.1 0.45 "$$ Fint$$" -ttext 0.1 0.6 "$$ \frac{\partial Fint}{\partial n_B}$$" -ttext 0.67 0.75 "$$ \frac{\partial^2 Fint}{\partial n_B^2}$$" -ttext 0.5 0.95 " $$ Fint$$, $$ \frac{\partial Fint}{\partial n_B}$$, and $$ \frac{\partial Fint}{\partial n_B}$$ vs $$ n_B$$ " \
	       -deriv nB Fint Fp -function "Fp/10^3" Fs2 -plot nB Fs2\
	       -deriv2 nB Fint Fw -function "Fw/10^3" Fs3 -plot nB Fs3\
       	       -save fint4_after.png


makeimg5:
	o2graph -read /home/awsteiner/wcs/eos/fid_3_14_23.o2 Fint -set logx 1 -set xlo 0.01 -set xhi 0.1 -set ylo "(-1)" -set yhi 1 -to-table 0 nB Fint 0.04 1.77839 -function "Fint/10^4" Fs -plot nB Fs -xtitle "$$ n_B$$" -ytitle "$$ Fint$$" -ttext 0.1 0.45 "$$ Fint$$" -ttext 0.1 0.55 "$$ \frac{\partial Fint}{\partial n_B}$$" -ttext 0.67 0.4 "$$ \frac{\partial^2 Fint}{\partial n_B^2}$$" -ttext 0.5 0.95 " $$ Fint$$, $$ \frac{\partial Fint}{\partial n_B}$$, and $$ \frac{\partial Fint}{\partial n_B}$$ vs $$ n_B$$ " \
	       -deriv nB Fint Fp -function "Fp/10^4" Fs2 -plot nB Fs2\
	       -deriv2 nB Fint Fw -function "Fw/10^4" Fs3 -plot nB Fs3\
       	       -save fint5_after.png

testsegfault:
	./eos_nuclei -load /home/awsteiner/wcs/eos/fid_3_14_23.o2 -interp-file /home/jbaut001/st.o2 ~/eos/fid_10_18_23.o2 4

-include makefile.aws
