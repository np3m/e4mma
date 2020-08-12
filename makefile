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

eos_nuclei: eos.o main.o eos_nuclei.o eos_had_skyrme_ext.o \
		virial_solver_deriv.h
	$(LMPI_CXX) $(LMPI_CFLAGS) -o eos_nuclei eos.o main.o \
		 eos_nuclei.o eos_had_skyrme_ext.o $(LIBS) \
		-lreadline

# ----------------------------------------------------------------
# Version without MPI
# ----------------------------------------------------------------

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
	cd doc; cp ~/o2scl/doc/o2scl/sphinx/build/html/objects.inv o2scl_objects.inv
	cd doc; cp ~/o2scl/doc/o2scl/part/sphinx/build/html/objects.inv o2scl_part_objects.inv
	cd doc; cp ~/o2scl/doc/o2scl/eos/sphinx/build/html/objects.inv o2scl_eos_objects.inv
	cd doc; doxygen doxyfile
	cd doc; make html

sync-doc:
	rsync -Cavzu doc/build/html/* $(STATIC_DOC_DIR)/eos

test-sync:
	rsync -Cavzun doc/build/html/* $(STATIC_DOC_DIR)/eos

clean:
	rm -f *.o eos_nuclei eos_nuclei_nompi eos eos_nompi

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

cv0:
	eos_nuclei_nompi -select-model $(P_FIDUCIAL) -check-virial

md0:
	eos_nuclei -select-model $(P_FIDUCIAL) -mcarlo-data

# ----------------------------------------------------------------
# Test XD's tables
# ----------------------------------------------------------------

testx:
	eos_nuclei -set full_results 1 -select-model $(P_FIDUCIAL) \
		-point-nuclei-aws 0.05 0.5 0.1 \
		~/data/eos/dsh_eos_fiducial.o2

# ----------------------------------------------------------------
# Recompute XD's tables
# ----------------------------------------------------------------

restart:
	cp ~/data/eos/dsh_eos_fiducial.o2 data/check_fiducial.o2
	eos_nuclei_nompi -set full_results 1 -select-model $(P_FIDUCIAL) \
		-edit-data data/check_fiducial.o2 \
		"inB%13==0 && iYe%13==0 && iT%13==0" flag 5 \
		data/check_fiducial.o2

# 35141 is just a random prime number which gives 8144 matches
restart2:
	cp ~/data/eos/dsh_eos_fiducial.o2 data/check_fiducial.o2
	eos_nuclei_nompi -set full_results 1 -select-model $(P_FIDUCIAL) \
		-edit-data data/check_fiducial.o2 \
		"floor(sin(inB*nYe*nT+iYe*nT+iT)*35141)%35141==0" flag 5 \
		data/check_fiducial.o2

restart3:
	cp ~/data/eos/dsh_eos_fiducial.o2 data/check_fiducial.o2
	eos_nuclei_nompi -set full_results 1 -select-model $(P_FIDUCIAL) \
		-edit-data data/check_fiducial.o2 \
		"abs(Ye-0.1)<0.001 && nB>0.04 && nB<0.07" flag 5 \
		data/check_fiducial.o2

recompute:
	eos_nuclei_nompi -set full_results 1 -select-model $(P_FIDUCIAL) \
		-set function_verbose 12211 \
		-generate-table data/check_fiducial.o2 \
		data/check_fiducial.o2

recompute2:
	eos_nuclei_nompi -set full_results 1 -select-model $(P_FIDUCIAL) \
		-set six_neighbors 1 \
		-generate-table data/check_fiducial.o2 \
		data/check_fiducial.o2

testx2:
	eos_nuclei_nompi -set full_results 1 -select-model $(P_FIDUCIAL) \
		-point-nuclei-aws 0.050237 0.19 12.86536 \
		~/data/eos/dsh_eos_fiducial.o2

testx3:
	eos_nuclei_nompi -set full_results 1 -select-model $(P_FIDUCIAL) \
		-point-nuclei-aws 0.00289088 0.5 0.9910964 \
		~/data/eos/dsh_eos_fiducial.o2

testx4:
	eos_nuclei_nompi -set full_results 1 -select-model $(P_FIDUCIAL) \
		-point-nuclei-aws 0.03 0.1 3.0 \
		~/data/eos/dsh_eos_fiducial.o2

mn1:
	eos_nuclei_nompi -set full_results 1 \
		-mcarlo-nuclei

# -set function_verbose 12211 

# ----------------------------------------------------------------
# Recompute derivatives
# ----------------------------------------------------------------

fid_deriv:
	eos_nuclei_nompi -set full_results 1 -select-model $(P_FIDUCIAL) \
		-eos-deriv ~/data/eos/dsh_eos_fiducial.o2 \
		~/data/eos/dsh_fid_deriv.o2 

fid_eg:
	eos_nuclei_nompi -set full_results 1 -select-model $(P_FIDUCIAL) \
		-add-eg ~/data/eos/dsh_fid_deriv.o2 \
		~/data/eos/dsh_fid_eg.o2 

fid_mt:
	eos_nuclei_nompi  -set full_results 1 -select-model $(P_FIDUCIAL) \
		-maxwell-test ~/data/eos/dsh_fid_eg.o2 0.4 1.0 mt.o2

fid_mt2:
	eos_nuclei_nompi  -set full_results 1 -select-model $(P_FIDUCIAL) \
		-maxwell-test ~/data/eos/dsh_fid_eg.o2 0.4 0.2 mt2.o2
