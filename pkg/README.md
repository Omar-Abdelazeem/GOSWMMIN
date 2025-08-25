# GOSWMMIN

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**Guided Operation of SWMM for Intermittent Networks**: A Python package for SWMM-based simulations of Intermittent Water Supply Networks.

## Overview

GOSWMMIN provides tools to model and analyze the behavior of water supply networks that do not operate continuously (Intermittent Water Supply - IWS). The package converts EPANET hydraulic network models to SWMM-compatible formats suitable for intermittent supply simulation with pressure-dependent withdrawal modeling.

## Key Features

- **EPANET to SWMM Conversion**: Seamlessly convert EPANET `.inp` files to SWMM format for IWS modeling
- **Pressure-Dependent Withdrawal (PDW)**: Implement Wagner et al. (1989) PDW formulation
- **Flexible Configuration**: Support for both uniform and node-specific parameters via CSV inputs
- **Adaptive Discretization**: Configurable spatial and temporal discretization schemes
- **User Storage Tanks**: Model household-level storage with customizable tank properties
- **Consumption Patterns**: Apply hourly demand patterns for realistic usage simulation
- **Results Analysis**: Extract and analyze pressure, flow, and tank level time series

## Installation

### Requirements

- Python 3.8 or higher
- SWMM 5.1+ (automatically handled by PySWMM)

### Method 1: Install from Source (Recommended)

```bash
# Clone the repository
git clone https://github.com/omaraliamer/GOSWMMIN.git
cd GOSWMMIN/pkg

# Option A: Direct pip install
pip install -e .

# Option B: Create conda environment first (recommended)
conda env create -f environment.yaml
conda activate goswmmin
pip install -e .
```

### Method 2: Conda Environment Only

If you prefer to use conda for dependency management:

```bash
# Create the environment
conda env create -f environment.yaml
conda activate goswmmin

# Install the package in development mode
pip install -e .
```

### Method 3: Minimal Installation

For a minimal installation with just the core dependencies:

```bash
conda env create -f environment-minimal.yaml
conda activate goswmmin-minimal
pip install -e .
```

### Troubleshooting Installation Issues

#### Common Issues with the Old GOSWMMIN.yaml

The original `GOSWMMIN.yaml` file may fail because:

1. **Platform-specific builds**: Contains Windows-specific packages that won't work on macOS/Linux
2. **Overly restrictive versions**: Exact build strings that may not be available
3. **Outdated packages**: Some versions may no longer be available

#### Solutions

1. **Use the new environment files**: `environment.yaml` or `environment-minimal.yaml`
2. **Manual installation**: Install dependencies individually:
   ```bash
   conda create -n goswmmin python=3.9
   conda activate goswmmin
   conda install numpy pandas matplotlib tqdm -c conda-forge
   pip install wntr pyswmm
   pip install -e .
   ```

3. **Update existing environment**:
   ```bash
   conda env update -f environment.yaml
   ```

### Dependencies

Core dependencies are automatically installed:

- `numpy` - Numerical computations
- `pandas` - Data manipulation and analysis  
- `matplotlib` - Plotting and visualization
- `wntr` - EPANET interface for Python
- `pyswmm` - SWMM interface for Python
- `tqdm` - Progress bars

## Quick Start

### Basic Usage

```python
from goswmmin import SWMMINSimulation

# Initialize simulation with EPANET file
sim = SWMMINSimulation("network.inp")

# Convert to SWMM format with basic configuration
swmm_file = sim.convert_to_swmmin(
    supply_duration=8.0,        # 8-hour supply period
    minimum_pressure=10.0,      # 10 m minimum pressure
    desired_pressure=20.0,      # 20 m desired pressure  
    pdw_exponent=0.5           # PDW exponent
)

# Run simulation
sim.run_swmmin()

# Extract results
pressures = sim.get_pressures()
tank_volumes, tank_heights = sim.get_tank_vols_heights()
```

### Advanced Usage with CSV Inputs

```python
# Use CSV files for node-specific parameters
swmm_file = sim.convert_to_swmmin(
    supply_duration=6.0,
    minimum_pressure="min_pressure.csv",      # Node-specific values
    desired_pressure="des_pressure.csv",      # Node-specific values
    pdw_exponent="pdw_exponent.csv",         # Node-specific values
    q_des="desired_flows.csv",               # Custom desired flows
    tank_areas="tank_areas.csv",             # Tank specifications
    tank_heights="tank_heights.csv",         # Tank heights
    consum_pattern="consumption_pattern.csv", # 24-hour demand pattern
    adaptive_disc=True,                      # Adaptive discretization
    solution_speed=150.0,                    # Solution speed (m/s)
    n_days=2                                 # Multi-day simulation
)
```

## Configuration Parameters

### Required Parameters

- `supply_duration` (float): Supply period duration in hours
- `minimum_pressure` (float/str): Minimum pressure for PDW (m)
- `desired_pressure` (float/str): Desired pressure for PDW (m)  
- `pdw_exponent` (float/str): PDW equation exponent

### Optional Parameters

- `n_days` (int): Number of simulation days (default: 1)
- `length_to_diameter` (float): L/D ratio for discretization (default: 30)
- `adaptive_disc` (bool): Use adaptive discretization (default: False)
- `maximum_xdelta` (float): Override spatial discretization (default: None)
- `solution_speed` (float): Solution speed in m/s (default: 100)
- `timestep` (float): Override temporal discretization (default: None)
- `leak_fraction` (float): Leakage as fraction of demand (default: 0.1)
- `q_des` (str): Path to desired flow CSV (default: None)
- `tank_areas` (str): Path to tank areas CSV (default: None)
- `tank_heights` (float/str): Tank heights in m (default: 1.0)
- `consum_pattern` (str): Path to consumption pattern CSV (default: None)
- `pdw_variable` (str): 'PRESSURE' or 'DEPTH' (default: 'PRESSURE')

## CSV File Formats

### Pressure and Flow Parameters

```csv
NodeID,Value
J1,10.5
J2,12.0
J3,8.5
```

### Consumption Pattern (24-hour)

```csv
Hour,Multiplier
0,0.5
1,0.4
...
23,0.6
```

## Documentation

### Method Reference

#### `SWMMINSimulation(input_file)`

Initialize simulation object with EPANET input file.

#### `convert_to_swmmin(**kwargs)`

Convert EPANET model to SWMM format with specified parameters.

**Returns**: Path to generated SWMM input file

#### `run_swmmin()`

Execute the SWMM simulation.

#### `get_pressures(specific_nodes=None)`

Extract pressure time series results.

**Returns**: pandas DataFrame with pressure data

#### `get_tank_vols_heights(specific_nodes=None)`

Extract tank volume and height time series.

**Returns**: Tuple of (volumes_df, heights_df)

## Examples

See the `examples/` directory for complete usage examples:

- `basic_usage.py` - Basic simulation setup and execution
- `advanced_config.py` - Advanced configuration with CSV inputs
- `results_analysis.py` - Post-processing and visualization

## Testing

Run the test suite:

```bash
pip install pytest pytest-cov
pytest tests/ -v --cov=goswmmin
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Install development dependencies (`pip install -e .[dev]`)
4. Set up pre-commit hooks (`pre-commit install`)
5. Make your changes and add tests
6. Ensure tests pass (`pytest`)
7. Commit your changes (`git commit -m 'Add amazing feature'`)
8. Push to the branch (`git push origin feature/amazing-feature`)
9. Open a Pull Request

### Development Dependencies

Install development tools:

```bash
pip install -e .[dev,docs,examples]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use GOSWMMIN in your research, please cite:

```bibtex
@article{goswmmin2024,
  title={Modeling Intermittent Water Supply in SWMM: A Critical Review with Reproducible Recommendations and a Python Package},
  author={Abdelazeem, Omar and [Co-authors]},
  journal={[Journal Name]},
  year={2024},
  note={In preparation}
}
```

## References

- Wagner, J. M., Shamir, U., & Marks, D. H. (1989). Water distribution reliability: simulation methods. *Journal of Water Resources Planning and Management*, 114(3), 276-294.
- Klise, K. A., et al. (2017). A software framework for assessing the resilience of drinking water systems to disasters with an example earthquake case study. *Environmental Modelling & Software*, 95, 420-431.
- McDonnell, B., et al. (2020). PySWMM: The Python Interface to Stormwater Management Model (SWMM). *Journal of Open Source Software*, 5(52), 2292.

## Support

- **Issues**: [GitHub Issues](https://github.com/omaraliamer/GOSWMMIN/issues)
- **Documentation**: [README](https://github.com/omaraliamer/GOSWMMIN#readme)
- **Email**: omaraliamer98@gmail.com

---

**GOSWMMIN** - Enabling better modeling of intermittent water supply systems

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
