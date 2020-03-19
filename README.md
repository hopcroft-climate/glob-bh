# glob-bh

 This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
-------------------------------------------------------------------------------------
Reversible Jump Markov chain Monte Carlo sampling algorithm for inferring ground temperature history from the global database of borehole temperature-depth profiles.
March 2020.

Users of this code should cite all of the following 3 papers:
1. P.O. Hopcroft, K. Gallagher and C. Pain (2007). Inference of past climate from borehole temperature data using Bayesian Reversible Jump Markov chain Monte Carlo, Geophysical Journal International, 171(3), 1430-1439, doi:10.1111/j.1365-246X.2007.03596.x.
2. P.O. Hopcroft, K. Gallagher and C. Pain (2009). A Bayesian Partition Modelling approach to resolve spatial variability in climate records from borehole temperature inversion , Geophysical Journal International, 178, 2, 651-666, doi:10.1111/j.1365-246X.2009.04192.x.
3. P.O. Hopcroft and K. Gallagher (submitted). Hemispheric divergence in ground warming from a self-adapting Bayesian analysis of the global temperature-depth database.

-------------------------------------------------------------------------------------
Instructions:

1. Contents:
This code should include the following directories:
DATA
Debug
Samples
TNT
Subsampler
matlab

DATA/ contains the data files you wish to analyse for past climate.
In this folder a wget script allows all of the required borehole geothermal data to be obtained from the public archive.
The build_script.scr converts all of the data to the format required by the C++ code.

BH_get_grid_squares/ contains a separate C++ routine to allocate the borehole sites to a 5x5 degree grid.

Debug/ contains the make file for the code

TNT/ is an open source library of matrix and vectors C++ classes required by the main part of the RJ-MCMC code.

Samples/ is where output files from the RJ_MCMC method are written to.

Subsampler/ contains source code for the RJ-MCMC algorithm as it operates within each gridcell of the globe.

matlab/ contains matlab .m scripts for analysing the resultant output and making plots.

--------------------------------------------------
2. Controlling the code:
Most of the code is controlled with the parameters.h file. This sets several global constants including the setup of the finite element 1D heat transfer model (1D FE) and the parameters of the RJ-MCMC.
iterations, a, dt, dx and e control the 1D FE model.
Parameters that control the RJ-MCMC sampler are set in MCMC_Mark_IV.h.

--------------------------------------------------
3. Compiling the code
The code is compiled by running a make command with Debug and copying the resulting .exe file back to the base directory. The code is setup to use icc with O2 optimisation, in the Debug/subdir.mk file.

--------------------------------------------------
4. Running the code
Currently set up to run 1012 borehole profiles across 173 partitions. On a single core this will take around 48 hours and generate a Samples directory of around 2.0 GB.

--------------------------------------------------
5. Analysing the output
This will require matlab code in the matlab/ directory.



First release 19.3.2020:
[![DOI](https://zenodo.org/badge/248197045.svg)](https://zenodo.org/badge/latestdoi/248197045)
