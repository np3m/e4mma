Todo List
=========
### Experimental
To compute EoS at a certain $\mathrm{nB~(fm^{-3})}$ , $\mathrm{Ye}$ and $\mathrm{T~(MeV)}$, Generate a configuration, using `point_generator.py` in the `src` folder to create a user specific `point.yaml` like:
```
python3 point_generator.py \
	--select_model "470 738 0.5 13.0 62.4 32.8 0.9" \
	--a_virial 10 --b_virial 10 \
	--load ../data/fid_3_5_22.o2 \
	--point_nuclei "0.16 0.4 30" 
```

#### Possible inputs
##### Options:
Select as many as you like.

`select_model`: Select an EOS model. The possible inputs are 
                P_FIDUCIAL=`"470 738 0.5 13.0 62.4 32.8 0.9"`
                P_LARGE_MMAX=`"783 738 0.5 13.0 62.4 32.8 0.9"`
                P_SMALL_R=`"214 738 0.5 13.0 62.4 32.8 0.9"`
                P_SMALLER_R=`"256 738 0.5 13.0 62.4 32.8 0.9"`
                P_LARGE_R=`"0 738 0.5 13.0 62.4 32.8 0.9"`
                P_SMALL_SL=`"470 738 0.5 13.0 23.7 29.5 0.9"`
                P_LARGE_SL=`"470 738 0.5 13.0 100.0 36.0 0.9"`
                default is `"470 738 0.5 13.0 62.4 32.8 0.9"`
                Select only one

`load`: Loads an EOS table in to memory. default: `"../data/fid_3_5_22.o2"`

`select_high_T`: Select 0 for the original DSH fit, 1 for NRAPR, 2 for Sk chi 414, 3 for
                Skchi450, 4 for Skchi500, 5 for ?, "+ and 6 for Sk chi m* (the default).
              default: 6

`eos_deriv`: Compute derivatives numerically.
              default: `'0'`


`a_virial`: Coefficient for modulation of virial EOS. default: `10.0`

`b_virial`: Coefficient for modulation of virial EOS. default: `10.0`

`include_muons`: If true, include muons. default: `false`

`max_ratio`: The maximum value of N/Z or Z/N. default: `7.0`

`mh_tol_rel`: Relative tolerance for the solver in the `eos_fixed_dist()` function. default: `1.0e-06`

`recompute`: If true, recompute all points, irrespective of the value of the convergence flag. default: `false`

#### commands:
Select only one

`get`: This command gets the value of a parameter. default: `'a_virial' `
        
`point`: Evaluate the EOS at one (nB,Ye,T) point.
              default: `"0.16 0.4 0.1"`

`point_nuclei`: Compute and/or show EOS results at one `(n_B,Y_e,T)` point. default: `"0.16 0.4 0.1"`
              
`random`: Select a random EOS, checking several physical constraints and re-selecting a
                new random EOS until all the constraints are met.
              default: `'0'`

Run `test_conf.sh` script inside the `test` folder using
```
bash test_conf.sh
```
To see a terminal output of the computed quantities. However functionality of this script is currently limited and more options will be added later. To use the eos fully, use makefile and CLI which are avilable.
.. todolist::
