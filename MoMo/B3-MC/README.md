# Monte Carlo Code

General performing MC simulations of the 2D Ising Model with a square lattice and analyze the resulting data using binning.

`practice3.f90` : Fortran code developed to perform MC simulations of the Ising Model with square lattice. The main inputs ( $L$, $T$, $n_{MCS}$, $...$) can be changed within the code. 

`r1279.f90` : Random generator used.

`makefile` : Makefile that compiles the program (the main practice3.f90 and the random functions used)

`stats.py` : Python script that reads output from the MC program and estimates statistical error with binning calculations and plots time series.

`/runs3` : Contains the data and plots for each run performed. The `runs_stats.ipynb` in this folder is a notebook that collects the data (time series and mean values and errors) from the different runs and plots it together.

`results2DIsing.ods` : Contains the results gathered from the different runs: statistical mean, error and autocorrelation time for $E$ and $|M|$.

The general workflow used has been first setting up the desired inputs, then compile the program throught the make file using make: 

```
   $ make 
```
run the executable created:
```
   $ ./main.x 
```
Once finished, analize the data with the python script:
```
   $ python3 stats.py 
```
And move the generated files to its correspondent folder in `/runs3`.
