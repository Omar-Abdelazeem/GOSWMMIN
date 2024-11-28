# GOSWMMIN

Guided Operation of SWMM for Intermittent Networks: a package for SWMM-based simulations of Intermittent Water Supply Networks. This repository containas the GOSWMMIN package in addition to materials and code used to produce the figures and analysis of the associated publication entitled **"Modeling Intermittent Water Supply in SWMM: A Critical Review with Reproducible Recommendations and a Python Package"** 

## Table of Contents

- Introduction
- Installation
- Contents:
  - [Notebooks](#notebooks)
  - [Networks](#networks)
  - [Replication of Figures](#replication-of-figures)
  - [Tools](#tools)
  - [Package](#package)
- License
- [References](#references)

## Introduction

GOSWMMIN is a package designed for the simulation of Intermittent Water Supply Networks using the Storm Water Management Model (SWMM). This package provides tools and utilities to model and analyze the behavior of water supply networks that do not operate continuously.

## Installation

To install GOSWMMIN, clone the repository, install the required dependencies and install the package using pip:

```sh
git clone <repository-url>
cd GOSWMMIN
pip install ./pkg
```

## Contents  
  
### Notebooks  

 This [folder](./Notebooks/) contains a number of Jupyter notebooks that explain the process undertaken by GOSWMMIN to create, execute and process SWMMIN simulations step-by-step. These procedures correspond to and expand upon the details present in the associated publication.
 [Conversion_to_SWMMIN.ipynb](./Notebooks/Conversion_to_SWMMIN.ipynb): Converts an input EPANET simulation into a SWMMIN simulation
 [Run_SWMMIN.ipynb](./Notebooks/Run_SWMMIN.ipynb): Executes a SWMMIN simulation and post-processes and visualizes its output

### Networks

This directory contains all network models used in this study in the form of EPANET and SWMM input files, in addition to some network-specific results (e.g., Continuity error + computational time results) and helper files. Four network models were used in the analysis presented in the main text and the supporting information:

  1. The [Linear network](./Networks/Linear%20Network/): A network consisting of four 1,000 m long pipes serving four demand nodes of different, zig-zaging elevations. The network was originally presented by [Gupta and Bhave (1996)](#references).  
  2. The [Modena Network](./Networks/Modena/): A model of the WDN in Modena, Italy, originally presented by [Bragalli et al. (2012)](#references)  
  3. The [Farina Network](./Networks/Farina%20et%20al%20(2014)/): A model of a WDN in northern Italy originally presented by [Farina et al. (2014)](#references)  
  4. The [Pescara Network](./Networks/Pescara/): A model of the WDN in Pescara, Italy originally presented by [Bragalli et al. (2012)](#references)  

Note that, by default, all generated input files and post-processed result files will be saved in the same directory as their corresponding original input file, i.e., Modena's SWMMIN input files will be in the Modena directory.

### Replication of Figures  

 This [folder](./Figures/) contains Jupyter notebooks that generate some of the reproducible figures in the associated publication which can be found in the [Figure Notebooks](./Figures/). This includes:

- **Figure 3 - Noise reduction**: Regenerates the data for and reproduces Figure 3.  
- **Figure 4 - Spatiotemporal discretization w Timesteps**:  *Only reprdocues* the figure from data stored in the Networks directory. The data for this figure can be generated using the scripts in the [Tools](./Tools/) directory. (see [Tools](#tools))  
- **Figure 5 - Spatiotemporal Discretization w solution speeds**: *Only reprdocues* the figure from data stored in the Networks directory. This notebook also generates Figure S10, The data for this figure can be generated using the scripts in the [Tools](./Tools/) directory. (see [Tools](#tools))  
  
### Tools
  
This directory contains python scripts that can be used to assess how the computational efficiency and mass balance of SWMMIN simulations change with the spatial and temporal discretization resolutions ($\Delta x_{max}$ and $\Delta t$).  
These scripts were used to generate the data in Figures [4](./Figures/Figure%20Files/Figure%204-Modena.png) and [5](/Figures/Figure%20Files/Figure%205%20Modena.png) as well as supporting figures [S8](./Figures/Figure%20Files/Figure%20S8%20Farina%20et%20al.png), [S9](./Figures/Figure%20Files/Figure%20S9%20Pescara.png) and S10 (Figure S11 to Figure 5).  

- [Continuity_Error_Computational_Time_w_timesteps.py](./Tools/Continuity_Error_Computational_Time_w_timesteps.py): Evaluates the continuity error and the average execution time for a specific network at the input spatial resolutions and timesteps.  
- [Continuity_Error_Computational_Time_w_solspeeds.py](./Tools/Continuity_Error_Computational_Time_w_solspeeds.py): Evaluates the continuity error and the average execution time for a specific network at the input spatial resolutions and solution speeds. 

### Package  

The GOSWMMIN package is structured as a new class ```SWMMIN_sim()``` which is initialized by an EPANET .inp file. The class's main methods include:

- ```Convert_to_SWMMIN()```: this method creates a SWMMIN simulation from the input EPANET file using a combination of built-in recommended settings and input arguments
- ```Run_SWMMIN()```: this method executes the SWMMIN simulation and parses its output

![image](./Figures/Figure%20Files/Figure%206.png)
further details on the methods, inputs and outputs of the GOSWMMIN package can be found in the package's [README](./pkg/README.md)

#### Dependencies

The dependencies of the GOSWMMIN package are listed in [requirements.txt](./pkg/requirements.txt). These include:

- WNTR: the Water Network Tool for Resilience, a package for EPANET in python [Klise et al. (2017)](#references)
- PySWMM: a Python wrapper for the Storm Water Management Model [McDonell et al. (2020)](#references)

## License

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.

## References

- Bragalli, C., D’Ambrosio, C., Lee, J., Lodi, A., & Toth, P. (2012). On the optimal design of water distribution networks: a practical MINLP approach. Optimization and Engineering, 13(2), 219–246. <https://doi.org/10.1007/s11081-011-9141-7>  
- Farina, G., Creaco, E., & Franchini, M. (2014). Using EPANET for modelling water distribution systems with users along the pipes. Civil Engineering and Environmental Systems, 31(1), 36–50. <https://doi.org/10.1080/10286608.2013.820279>  
- Gupta, R., & Bhave, P. R. (1996). Comparison of Methods for Predicting Deficient-Network Performance. Journal of Water Resources Planning and Management, 122(3), 214–217. <https://doi.org/10.1061/(ASCE)0733-9496(1996)122:3(214)>  
- Klise, K. A., Bynum, M., Moriarty, D., & Murray, R. (2017). A software framework for assessing the resilience of drinking water systems to disasters with an example earthquake case study. Environmental Modelling & Software, 95, 420–431. <https://doi.org/10.1016/j.envsoft.2017.06.022>
- McDonnell, B., Ratliff, K., Tryby, M., Wu, J., & Mullapudi, A. (2020). PySWMM: The Python Interface to Stormwater Management Model (SWMM). Journal of Open Source Software, 5(52), 2292. <https://doi.org/10.21105/joss.02292>
