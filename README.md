# mixing-model
Exploring the emergent behavioral patterns in heterogeneous ant colonies using the fixed threshold model. This is the living repository of this code.

## Components of this Repository
This directory has three main files:
* **scripts**: contains all scripts for simulation, analysis, and plotting of the theoretical model.
* **output**: contains all subsequent derived data from simulations as well as any graphs produced from analysis of the simulation data.
* **numerics**: contains all notebooks for calculations in the SI.

The scripts folder has the following structure:
* Scripts starting with **\_\_Util__** contains relevant utility functions used throughout the other scripts. These are broken down by general functionality. **__Util__MASTER.R** is the script sourced in all other scripts that imports all the other utility functions. 
* Scripts starting with **1** contain the general simulation scripts.
* Scripts starting with **2** allow one to plot some of the simulation results.
* **deprecated** contains older versions of scripts that were useful at one point or another but are no longer central to the main text. 
