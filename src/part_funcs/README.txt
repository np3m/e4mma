These files contain high temperature partition functions calculated with 
input from the FRDM (datafile2.txt) and the ETFSI-Q (datafile3.txt)
mass formula. Nuclear properties (microscopic correction) derived from
these mass formulas enter the computation of the nuclear level density,
which in turn is essential for the partition functions. The relevant
publication is Ap. J. Suppl. 147 (2003) 403.

The files are intended to be used as extensions of the partition
function files published in ADNDT 75 (2000) 1 (to be found in the parent 
directory to the one where this readme file is residing). They are extending 
the temperature grid up to T9=275. The level density approach used is that
of Phys. Rev. C 56 (1997) 1613. They are including high temperature
corrections similar to those of Ap. J. 232 (1979) L59. For more details,
see T. Rauscher, Ap. J. Suppl. 147 (2003) 403 (ms1.ps, ms1.pdf).

The structure of the files is quite similar to the ones published in ADNDT,
except that the partition functions are given in e8.2 format. An entry
for a given isotope looks like this:

isotope name (a5) (slightly different than in ADNDT)
Z, A, ground state spin (2i4,f5.1)
48 partition functions (1x,1p,8e9.2)

The 48 partition functions are given for the temperature grid (T9):
12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0,
28.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,
65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,
105.0,110.0,115.0,120.0,125.0,130.0,135.0,140.0,
145.0,150.0,155.0,160.0,165.0,170.0,175.0,180.0,
190.0,200.0,210.0,220.0,230.0,240.0,250.0,275.0



Thomas Rauscher, December 2002 / January 2003
