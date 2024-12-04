#%%
from SWMMIN_sim import SWMMIN_sim
import pathlib


# Creating a SWMMIN simulation object with optional CSV inputs for PDW variables
# tank sizes and consumption patterns
filepath = pathlib.Path('../../Networks/Linear Network/Linear_Network.inp')
sim = SWMMIN_sim(filepath)
min_pressure = '../../Networks/Linear Network/min_pressure.csv'
des_pressure = '../../Networks/Linear Network/des_pressure.csv'
pdw_exponent = '../../Networks/Linear Network/pdw_exponent.csv'
q_des = '../../Networks/Linear Network/q_des.csv'
tank_areas = '../../Networks/Linear Network/tank_areas.csv'
tank_heights = '../../Networks/Linear Network/tank_heights.csv'
consum_pattern = '../../Networks/Linear Network/consum_pattern.csv'

# Create the SWMMIN simulation
# Creates an 8-hr supply duration with a non-adaptive
sim.Convert_to_SWMMIN(supply_duration= 8.0, 
                  minimum_pressure=min_pressure, desired_pressure = des_pressure, pdw_exponent=pdw_exponent,
                  n_days= 1,length_to_diameter=30, solution_speed=100, q_des=q_des, tank_heights=tank_heights, 
                  tank_areas=tank_areas, consum_pattern=consum_pattern)

#%%
# Running the simulation
sim.Run_SWMMIN()

#%%
# Getting Pressure Results
sim.get_pressures()