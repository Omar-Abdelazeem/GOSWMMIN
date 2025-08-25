"""Test configuration for pytest."""

import pytest
import pathlib
import tempfile
import shutil
from typing import Generator

import pandas as pd
import numpy as np


@pytest.fixture
def temp_dir() -> Generator[pathlib.Path, None, None]:
    """Create a temporary directory for testing."""
    temp_path = pathlib.Path(tempfile.mkdtemp())
    yield temp_path
    shutil.rmtree(temp_path)


@pytest.fixture  
def sample_epanet_file(temp_dir: pathlib.Path) -> pathlib.Path:
    """Create a sample EPANET input file for testing."""
    epanet_content = """[TITLE]
Sample Network

[JUNCTIONS]
;ID              	Elev        	Demand      	Pattern         
 J1              	100         	10          	                	;
 J2              	95          	15          	                	;
 J3              	90          	20          	                	;

[RESERVOIRS]
;ID              	Head        	Pattern         
 R1              	150         	                	;

[TANKS]
;ID              	Elevation   	InitLevel   	MinLevel    	MaxLevel    	Diameter    	MinVol      	VolCurve
 T1              	120         	10          	5           	15          	50          	0           	                	;

[PIPES]
;ID              	Node1           	Node2           	Length      	Diameter    	Roughness   	MinorLoss   	Status
 P1              	R1              	J1              	1000        	300         	100         	0           	Open  	;
 P2              	J1              	J2              	800         	250         	100         	0           	Open  	;
 P3              	J2              	J3              	600         	200         	100         	0           	Open  	;
 P4              	J2              	T1              	400         	150         	100         	0           	Open  	;

[PUMPS]
;ID              	Node1           	Node2           	Parameters
 PU1             	R1              	J1              	HEAD 1	;

[VALVES]
;ID              	Node1           	Node2           	Diameter    	Type	Setting     	MinorLoss   

[PATTERNS]
;ID              	Multipliers
 1               	1.0 1.2 1.4 1.6 1.4 1.2 1.0 0.8 0.6 0.8 1.0 1.2 1.4 1.6 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.6 0.8

[CURVES]
;ID              	X-Value     	Y-Value
 1               	1000        	10
 1               	1500        	5
 1               	2000        	0

[CONTROLS]

[QUALITY]

[SOURCES]

[REACTIONS]

[MIXING]

[OPTIONS]
 Units           	LPS
 Headloss        	H-W
 Pattern         	1
 Demand Multiplier	1.0
 Emitter Exponent	0.5
 Quality         	None mg/L
 Diffusivity     	1.0
 Viscosity       	1.0
 Trials          	40
 Accuracy        	0.001
 Tolerance       	0.01
 Pressure Exponent	0.5

[REPORT]
 Status          	Yes
 Summary         	No
 Page            	0

[COORDINATES]
;Node            	X-Coord         	Y-Coord
 J1              	1000            	1000            
 J2              	2000            	1000            
 J3              	3000            	1000            
 R1              	0               	1000            
 T1              	2000            	2000            

[END]
"""
    
    epanet_file = temp_dir / "sample_network.inp"
    with open(epanet_file, 'w') as f:
        f.write(epanet_content)
    
    return epanet_file


@pytest.fixture
def sample_csv_data(temp_dir: pathlib.Path) -> dict:
    """Create sample CSV files for testing."""
    csv_files = {}
    
    # Minimum pressure CSV
    min_pressure_data = pd.DataFrame({
        'ID': ['J1', 'J2', 'J3'],
        'Min_Pressure': [5.0, 7.0, 8.0]
    })
    min_pressure_file = temp_dir / "min_pressure.csv"
    min_pressure_data.to_csv(min_pressure_file, header=False, index=False)
    csv_files['min_pressure'] = min_pressure_file
    
    # Desired pressure CSV
    des_pressure_data = pd.DataFrame({
        'ID': ['J1', 'J2', 'J3'],
        'Des_Pressure': [15.0, 17.0, 18.0]
    })
    des_pressure_file = temp_dir / "des_pressure.csv"
    des_pressure_data.to_csv(des_pressure_file, header=False, index=False)
    csv_files['des_pressure'] = des_pressure_file
    
    # PDW exponent CSV
    pdw_exponent_data = pd.DataFrame({
        'ID': ['J1', 'J2', 'J3'],
        'Exponent': [0.5, 0.5, 0.5]
    })
    pdw_exponent_file = temp_dir / "pdw_exponent.csv"
    pdw_exponent_data.to_csv(pdw_exponent_file, header=False, index=False)
    csv_files['pdw_exponent'] = pdw_exponent_file
    
    # Consumption pattern CSV (24 hours)
    pattern_data = pd.DataFrame({
        'Hour': range(24),
        'Multiplier': np.random.uniform(0.5, 1.5, 24)
    })
    pattern_file = temp_dir / "consumption_pattern.csv"
    pattern_data.to_csv(pattern_file, header=False, index=False)
    csv_files['consumption_pattern'] = pattern_file
    
    return csv_files


@pytest.fixture
def simulation_config() -> dict:
    """Standard simulation configuration for testing."""
    return {
        'supply_duration': 8.0,
        'minimum_pressure': 10.0,
        'desired_pressure': 20.0,
        'pdw_exponent': 0.5,
        'n_days': 1,
        'length_to_diameter': 30.0,
        'adaptive_disc': False,
        'solution_speed': 100.0,
        'leak_fraction': 0.1,
        'tank_heights': 1.0,
        'pdw_variable': 'PRESSURE'
    }
