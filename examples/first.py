
# # E4MMA - First example

# This example is designed to be used in conjuction with the E4mma
# docker images available at
# https://hub.docker.com/repository/docker/awsteiner/e4mma. See the
# documentation at https://np3m.org/code/e4mma for more information.

# +
import matplotlib.pyplot as plot
import numpy
import sys

plots=True
if 'pytest' in sys.modules:
    plots=False
# -

# Use the code to compute the "fiducial" EOS from Du et al. (2022)
# at nB=0.16, Ye=0.5, T=10 MeV. This is a relatively easy example
# because there are few nuclei at these temperatures. Of course
# the docker image includes the executable, so this can be
# done without Jupyter as well.

os.system('/home/np3m/e4mma/eos_nuclei '+
          '-select-model 470 738 0.5 13.0 62.4 32.8 0.9 '+
          '-point-nuclei 0.16 0.5 10.0')

# A more difficult example near the liquid-gas phase transition

os.system('/home/np3m/e4mma/eos_nuclei '+
          '-select-model 470 738 0.5 13.0 62.4 32.8 0.9 '+
          '-point-nuclei 0.08 0.5 1.0')

# Use the code to compute the NRAPR Skyrme model at nB=0.16, Ye0.5,
# T=5 MeV

os.system('/home/np3m/e4mma/eos_nuclei '+
          '-alt-model Skyrme NRAPR '+
          '-point-nuclei 0.16 0.5 10.0')

# Download the EOS table. The "tables" directory is created
# by the docker file.

os.system('curl https://isospin.roam.utk.edu/public_data'+
          '/eos_tables/du21/fid_3_5_22.o2 '+
          '-o /home/np3m/e4mma/tables/fid_3_5_22.o2')

# Import O$_2$sclpy and link the O$_2$scl library:

import o2sclpy
link=o2sclpy.linker()
link.link_o2scl()

# Get a copy (a pointer to) the O$_2$scl unit conversion object, which
# also allows access to the constant library, then get ħc.

cu=link.o2scl_settings.get_convert_units()
ħc=cu.find_unique('hbarc','MeV*fm')
print('ħc = %7.6e\n' % (ħc))

# Use the cloud_file object to download the EOS

cf=o2sclpy.cloud_file(link)
cf.verbose=1
cf.get_file('dsh.o2','https://isospin.roam.utk.edu/public_data'+
            '/eos_tables/du21/fid_3_5_22.o2')

# Read the tensor which stores the average mass number

hf=o2sclpy.hdf_file(link)
tg_A=o2sclpy.tensor_grid(link)
hf.open('dsh.o2')
o2sclpy.hdf_input_tensor_grid(link,hf,tg_A,'A')
hf.close()

# In order to make a plot at fixed Ye, we first need to construct a
# tensor index object. We want to include all values of nB (index 0 in
# the tensor object) and all values of T (index 2 in the tensor
# object), but for Ye, we select the value in the grid which is
# closest to Ye=0.4.

ix=o2sclpy.std_vector_size_t(link)
ix.resize(3)
ix[1]=tg_A.lookup_grid(1,0.4)

# Create a table3d object

t3d=o2sclpy.table3d(link)
tg_A.copy_table3d_align_setxy(0,2,ix,t3d,'nB','T','A')

# Now plot the results. Raw matplotlib works, but o2sclpy has
# a couple functions which make it easier. 
       
if plots:
    pl=o2sclpy.plot_base()
    pl.colbar=True
    pl.xtitle(r'$ n_B~(\mathrm{fm}^{-3}) $')
    pl.ytitle(r'$ T~(\mathrm{MeV}) $')
    pl.ttext(1.25,0.5,u'$ A $',rotation=90)
    pl.den_plot(t3d,'A')
    plot.show()

# For testing purposes
    
def test_fun():
    assert numpy.allclose(t3d.get(0,0,'A'),81,rtol=0.1)
    return


