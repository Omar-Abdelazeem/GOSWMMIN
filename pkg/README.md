# GOSWMMIN

#### Guided Operation of SWMM for Intermittent Networks: a package for SWMM-based simulations of Intermittent Water Supply Networks.  

 This repository containas the GOSWMMIN package in addition to materials and code used to produce the figures and analysis of the associated publication entitled **"Modeling Intermittent Water Supply in SWMM: A Critical Review with Reproducible Recommendations and a Python Package"** 

## Table of Contents

- Introduction
- Installation
- Dependencies
- Documentation:
  - [Creating a Simulation](#creating-a-simulation)
  - [Running a Simulation](#running-a-simulation)
  - [Processing Results](#processing-results)
  - [Example](#example)

- License
- [References](#references)

## Introduction

GOSWMMIN is a package designed for the simulation of Intermittent Water Supply Networks using the Storm Water Management Model (SWMM). This package provides tools to model and analyze the behavior of water supply networks that do not operate continuously, i.e., Intermittent Water Supply (IWS).

## Installation

To install GOSWMMIN, clone the repository, install the required dependencies and install the package using pip:

```sh
git clone <repository-url>
cd GOSWMMIN
pip install ./pkg/requirements.txt
pip install ./pkg
```

## Dependencies 

The GOSWMMIN package is built on several key dependencies and uses a variety of common Python data handling and visualization packages. A full list of the required packages and specifications for using GOSWMMIN are provided in the [Requirements.txt](./requirements.txt) file. Additionally, the [conda environment](./GOSWMMIN.yaml) file attached can be used to import a conda environment that contains all requirements and dependencies.  
A few prominent dependencies and requirements are listed below: 

1. **Python** = 3.9.20  
2. **WNTR** = 1.2.0
3. **PySWMM** = 2.0.1

## Documentation  

![image](./Resources/Figure%206.png)  
The GOSWMMIN package introduces a new class representing a SWMMIN simulation: ```SWMMIN_sim()```. This section documents the process of creating, executing and processing such a simulation using the ```SWMMIN_sim()``` object and its methods.

### Creating a Simulation

A SWMMIN simualtion is initialized by calling ```SWMMIN_sim(path)``` where ```path``` refers to an EPANET .inp file that is pre-configured. This **EPANET input file** is later used to generate the network model and read important inputs like pipe lengths and diameters, node elevations and demands ($Q_{demand}$), etc. The EPANET input file must contain demands assigned to each demand node (non-zero demand nodes) and that demand is assumed to be the base demand of the node over 24 hours, similar to a CWS base demand.

```python
#Initialize a SWMMIN simulation for the Linear Network
sim = SWIMMIN_sim(path/Linear_Network.inp)
```

Once the simulation is initialized (sim), the user can use ```Convert_to_SWMMIN()``` method to input key conversion inputs and details. The conversion function has the following inputs (inputs with a default value are optional):

- **Supply Duration** (in hours): input the desired duration of the intermittent supply in hours in float, e.g., 4:30 hours should be input as 4.5.  Currently only one supply duration per day are supported.  
- **PDW Parameters**: This group of inputs are used to define the Pressure Dependent Withdrawal (PDW) function for withdrawal outlets. Currently, this package uses the formulation by [Wagner et al. (1989)](#references) and as such, the following input are required:  
  - **Minimum Pressure** (m): The minimum pressure required for withdrawal at the nodes. Float/integer inputs will be used for all nodes. To define a specifc minimum pressure for each node, a path to a CSV file containing node IDs and their corresponding minimum pressure must be passed instead. See example CSV files below.  
  - **Desired Pressure** (m): pressure at which users can withdraw at the desired flow rate. float/integer input for uniform $H_{des}$ or csv input for node-specific assignment.
  - **Desired Flowrate** (m<sup>3</sup>/s): Flow rate withdrawn when pressure equals desired pressure. If None, $Q_{des}$ is inferred from the assigned $Q_{demand}$ and the input supply duration ($t_{supply}$) such that $Q_{des} = Q_{demand} \times 24 hr/ t_{supply}$. Otherwise, a path to a CSV specifying $Q_{des}$ for each node should be input.
  - **PDW Exponent**: the exponent of the PDW formula. Wagner et al. (1989) used 0.4. Float/integer for uniform exponents or CSV for node-specific assignment.  
  - **PDW Variable**: By default, the PDW formula is applied using the pressure difference between the demand node and the user tank ("PRESSURE"). Changing this option to "DEPTH" means that the freeboard depth or pressure at the demand node will be used instead.
- **Number of Days**: number of simulation days in whole numbers only. Default is 1 day. Additional days repeat the same supply duration once a day.  
- **Spatial  Discretization**: the following collection of inputs set the resolution and method of discretizing the SWMMIN simulation spatially. Two methods and an override option are available.
  - **Adaptive Discretization** (bool): default is ```False```. Determines whether the discretization will be adaptive to each pipe's diameter or uniform across all pipes.
  - **Length to Diameter Ratio**: default is 30. The ratio of the maximum discretized pipe length $\Delta x_{max}$ to the pipe diameter. If Adaptive is True, each pipe will have its own $\Delta x_{max}$ equal to $L/D \times D_i$. If Adaptive is False, all pipes will have the same $\Delta x_{max}$ based on the largest pipe diameter and equal to $L/D \times D_{max}$.  
  - **Maximum $\Delta x$**: Override option, default is None. Overrides previous inputs and applies the input $\Delta x_{max}$ unifromly on all pipes.
- **Temporal Discretization**: Sets the timestep for the SWMMIN simulation:
  - **Solution Speed** (m): The ratio of the spatial resolution to the temporal resolution: $\Delta x_{max} / \Delta t$. Default value is 100 m/s and the recommended range is 50-200 m/s. Used along with computed or input $\Delta x_{max}$ to set the timestep.  
  - **timestep** (s): Overrides solution speed input. Timestep input here is used directly. Default is None.  
- **User Tanks**: Defining the height and capacity of user storage tanks:
  - **Tank Heights** (m): default is 1. The heigh of user storage tanks: float/integer input for uniform assignment (e.g., all tanks 1 m high) or CSV for node specifc assignment.  
  - **Tank Areas** (m<sup>2</sup>): The area of storage tanks. Default is None, when None, areas are inferred from user demands where the area of each tank is equal to $Q_{demand} \times 24$ hr $/ h_{tank}$. Otherwise, pass a path ot a CSV with each tank's area.  
- **Leak Fraction**: The fraction of each $Q_{des}$ to be assigned as leakage in each node, e.g., if fraction = 0.1, 10% of each node's $Q_des$ will be taken as leakage using the same PDW parameters ($H_{des}, H_{min}, n$)  
- **Consumption Pattern**: Optional. Input a path to a CSV containing an **hourly** consumption pattern. Only hourly patterns are supported. The CSV should contain exactly 24 pairs of hour-multiplier.

#### Example  

In the [example](./GOSWMMIN/example.py) provided, an example of creating a SWMMIN simulation with all optional CSV inputs used.  For example, the minimum pressure is assigned by specifying a path to

```python
min_pressure = '../Resources/min_pressure.csv'
```

where ```min_pressure.csv``` is formatted as follows with node IDs and corresponding minimum pressures

```csv
1, 1
2, 2
3, 3
4, 4
```

after specifying all paths to CSVs, the conversion function is called

```python
sim.Convert_to_SWMMIN(supply_duration= 8.0, minimum_pressure=min_pressure, 
                  desired_pressure = des_pressure, pdw_exponent=pdw_exponent,
                  q_des=q_des, tank_heights=tank_heights, 
                  tank_areas=tank_areas, consum_pattern=consum_pattern)

```

Note that certain inputs must be defined as the same type: if minimum_pressure is a CSV, then desired_pressure and pdw_exponent must also be CSVs.  

A SWMM .inp file is generated after running the previous method. The path to the .inp file is saved in the ```sim``` object. 

### Running a Simulation 

Once the SWMMIN simulation (.inp file) has been created using ```Convert_to_SWMMIN()```, the simulation can be run in one of two ways:

1. Opening the .inp file in the standard SWMM 5 GUI and running it
2. Run in Python using GOSWMMIN's methods

When running the SWMMIN simulation in Python, modelers can use the ```Run_SWMMIN()``` method. The method requires no input and runs the SWMM input file generated from an earlier call of the ```Convert_to_SWMMIN()``` method. Alternatively, a path can optionally be specified to execute a SWMMIN simulation without converting it first, for running previously converted SWMMIN input files.  

In the [example](./GOSWMMIN/example.py) provided, the method is simply called to execute the simulation:  

```python
# Running the simulation
sim.Run_SWMMIN()
```

### Processing Results

Currently (in version 0.1.0), only a portion of the planned processing tools have been implemented.  Currently, modelers can use methods to retrieve key information from the simulation resutls, which can then be visualized by the modelers. For instance, methods for retrieving timeseries for node pressures, withdrawal rates, consumption rates, stored volumes and height of water inside user tanks are available in 0.1.0. As an example, we demonstrate the process of getting node pressures in the provided [example](./GOSWMMIN/example.py):  

```python
sim.get_pressures(specific_nodes=['DN1','DN2'])
```

This example also specifies certain nodes to retrieve the pressure results for. Similar procedures for retrieving other key inputs are available.  

In upcoming versions, further processing and visualization tools will be added 


## License

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.

## References

- Bragalli, C., D’Ambrosio, C., Lee, J., Lodi, A., & Toth, P. (2012). On the optimal design of water distribution networks: a practical MINLP approach. Optimization and Engineering, 13(2), 219–246. <https://doi.org/10.1007/s11081-011-9141-7>  
- Farina, G., Creaco, E., & Franchini, M. (2014). Using EPANET for modelling water distribution systems with users along the pipes. Civil Engineering and Environmental Systems, 31(1), 36–50. <https://doi.org/10.1080/10286608.2013.820279>  
- Gupta, R., & Bhave, P. R. (1996). Comparison of Methods for Predicting Deficient-Network Performance. Journal of Water Resources Planning and Management, 122(3), 214–217. <https://doi.org/10.1061/(ASCE)0733-9496(1996)122:3(214)>  
- Klise, K. A., Bynum, M., Moriarty, D., & Murray, R. (2017). A software framework for assessing the resilience of drinking water systems to disasters with an example earthquake case study. Environmental Modelling & Software, 95, 420–431. <https://doi.org/10.1016/j.envsoft.2017.06.022>
- McDonnell, B., Ratliff, K., Tryby, M., Wu, J., & Mullapudi, A. (2020). PySWMM: The Python Interface to Stormwater Management Model (SWMM). Journal of Open Source Software, 5(52), 2292. <https://doi.org/10.21105/joss.02292>
