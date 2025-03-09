import wntr
from pathlib import Path
import matplotlib.pyplot as plt

wn = wntr.network.WaterNetworkModel()


# Match Linear Network example from GOSWMMIN documentation
wn.options.hydraulic.headloss = 'H-W'  # Hazen-Williams headloss formula
wn.options.hydraulic.inpfile_units = 'LPS'  # Liters per second
wn.options.hydraulic.trials = 40
wn.options.hydraulic.specific_gravity = 1.0
wn.options.hydraulic.viscosity = 1.0
wn.options.hydraulic.accuracy = 0.001
wn.options.hydraulic.unbalanced = 'CONTINUE'  # Continue for 10 trials if unbalanced
wn.options.hydraulic.demand_multiplier = 1.0
wn.options.hydraulic.emitter_exponent = 0.5


wn.options.time.duration = 24 * 3600
wn.options.time.hydraulic_timestep = 3600
wn.options.time.quality_timestep = 300
wn.options.time.pattern_timestep = 3600
wn.options.time.pattern_start = 0 
wn.options.time.report_timestep = 60
wn.options.time.report_start = 0
wn.options.time.start_clocktime = 0
wn.options.time.statistic = 'NONE'

wn.options.energy.global_efficiency = 75
wn.options.energy.global_price = 0.0
wn.options.energy.demand_charge = 0.0

wn.options.quality.parameter = 'NONE' 
wn.options.quality.diffusivity = 1.0
wn.options.quality.tolerance = 0.01


wn.add_reservoir('1', base_head=200, head_pattern=None)

# TODO: verify that the diameter units are correct
BASE_DEMAND = 0.000666667

# For some reason, the whole simulation errors out if the diameter is too small
PINK_DIAM = 0.3
RED_DIAM = 0.3
GREEN_DIAM = 0.3


wn.add_junction('DN1', elevation=99.35, base_demand=BASE_DEMAND)  # Node 2 in diagram
wn.add_junction('DN2', elevation=99.65, base_demand=BASE_DEMAND)  # Node 7 in diagram
wn.add_junction('DN3', elevation=96.39, base_demand=BASE_DEMAND)  # Node 3 in diagram
wn.add_junction('DN4', elevation=96.30, base_demand=BASE_DEMAND)  # Node 4 in diagram
wn.add_junction('DN5', elevation=96.31, base_demand=BASE_DEMAND)  # Node 5 in diagram
wn.add_junction('DN6', elevation=95.88, base_demand=BASE_DEMAND)  # Node 6 in diagram
wn.add_junction('DN7', elevation=96.44, base_demand=BASE_DEMAND)  # Node 8 in diagram
wn.add_junction('DN8', elevation=96.14, base_demand=BASE_DEMAND)  # Node 9 in diagram
wn.add_junction('DN9', elevation=96.08, base_demand=BASE_DEMAND)  # Node 10 in diagram


wn.add_pipe("P1", "1", "DN1", length=30, diameter=PINK_DIAM) 
wn.add_pipe("P2", "DN1", "DN2", length=150, diameter=GREEN_DIAM)
wn.add_pipe("P3", "DN1", "DN3", length=30, diameter=PINK_DIAM)
wn.add_pipe("P4", "DN3", "DN4", length=40, diameter=PINK_DIAM)
wn.add_pipe("P5", "DN3", "DN7", length=110, diameter=GREEN_DIAM)
wn.add_pipe("P6", "DN4", "DN5", length=50, diameter=RED_DIAM)
wn.add_pipe("P7", "DN4", "DN8", length=50, diameter=GREEN_DIAM)
wn.add_pipe("P8", "DN5", "DN6", length=170, diameter=GREEN_DIAM)
wn.add_pipe("P9", "DN5", "DN9", length=70, diameter=GREEN_DIAM)


wn.get_node("1").coordinates = (0, 0)
wn.get_node("DN1").coordinates = (-30, 0)
wn.get_node("DN2").coordinates = (-30, -150)
wn.get_node("DN3").coordinates = (-30, 30)
wn.get_node("DN4").coordinates = (-30, 70)
wn.get_node("DN5").coordinates = (-30, 120)
wn.get_node("DN6").coordinates = (-30, 290)
wn.get_node("DN7").coordinates = (80, 30)
wn.get_node("DN8").coordinates = (-80, 70)
wn.get_node("DN9").coordinates = (-100, 120)


wntr.graphics.plot_network(wn, 
                           title="Florence Network",
                           show_plot=False,
                           node_labels=True,
                           filename="florence_network.png")


print(f"Network summary: {wn.describe()}")
Path("./Networks/Ismail").mkdir(parents=True, exist_ok=True)
# TODO: verify version of epa-net, 2.0 or 2.2
wntr.epanet.io.InpFile().write(filename="./Networks/Ismail/ismail.inp", 
                               wn=wn, 
                               units="LPS",
                               version=2.2)


# import wntr
# import matplotlib.pyplot as plt

# # Create a new water network model
# wn = wntr.network.WaterNetworkModel()

# # Add a reservoir with base head in meters
# wn.add_reservoir('1', base_head=100.0)  # Changed ID from 'R1' to '1' for consistency

# # Add junctions with elevation in meters and base demand in cubic meters per second (m³/s)
# wn.add_junction('DN1', elevation=90.0, base_demand=0.0666667)  # 66.6667 LPS converted to m³/s
# wn.add_junction('DN2', elevation=88.0, base_demand=0.0666667)  # 66.6667 LPS converted to m³/s
# wn.add_junction('DN3', elevation=90.0, base_demand=0.1)        # 100 LPS converted to m³/s
# wn.add_junction('DN4', elevation=85.0, base_demand=0.0333333)  # 33.3333 LPS converted to m³/s

# # Add pipes with length in meters and diameter in meters
# wn.add_pipe('P1', '1', 'DN1', length=1000.0, diameter=0.4, roughness=130.0, minor_loss=0.0, initial_status='OPEN')
# wn.add_pipe('P2', 'DN1', 'DN2', length=1000.0, diameter=0.35, roughness=130.0, minor_loss=0.0, initial_status='OPEN')
# wn.add_pipe('P3', 'DN2', 'DN3', length=1000.0, diameter=0.3, roughness=130.0, minor_loss=0.0, initial_status='OPEN')
# wn.add_pipe('P4', 'DN3', 'DN4', length=1000.0, diameter=0.3, roughness=150.0, minor_loss=0.0, initial_status='OPEN')

# # Set node coordinates (in meters)
# wn.get_node('1').coordinates = (359.955, 6636.670)
# wn.get_node('DN1').coordinates = (1034.871, 6602.925)
# wn.get_node('DN2').coordinates = (1507.312, 6490.439)
# wn.get_node('DN3').coordinates = (2047.244, 6569.179)
# wn.get_node('DN4').coordinates = (2407.199, 6355.456)

# # Set simulation options
# wn.options.time.duration = 24 * 3600  # 24 hours in seconds
# wn.options.time.hydraulic_timestep = 3600  # 1 hour in seconds
# wn.options.time.quality_timestep = 300  # 5 minutes in seconds
# wn.options.time.pattern_timestep = 3600  # 1 hour in seconds
# wn.options.time.pattern_start = 0  # Start at the beginning
# wn.options.time.report_timestep = 60  # 1 minute in seconds
# wn.options.time.report_start = 0  # Start reporting at the beginning
# wn.options.time.start_clocktime = 0  # Midnight (00:00:00)
# wn.options.time.statistic = 'NONE'  # No statistical reporting

# # Set hydraulic options
# wn.options.hydraulic.headloss = 'H-W'  # Hazen-Williams headloss formula
# wn.options.hydraulic.inpfile_units = 'LPS'  # Liters per second
# wn.options.hydraulic.specific_gravity = 1.0
# wn.options.hydraulic.viscosity = 1.0
# wn.options.hydraulic.trials = 40
# wn.options.hydraulic.accuracy = 0.001
# wn.options.hydraulic.unbalanced = 'CONTINUE'  # Continue for 10 trials if unbalanced
# wn.options.hydraulic.demand_multiplier = 1.0
# wn.options.hydraulic.emitter_exponent = 0.5

# # Set energy options
# wn.options.energy.global_efficiency = 75  # Percent
# wn.options.energy.global_price = 0.0
# wn.options.energy.demand_charge = 0.0

# # Set quality options
# wn.options.quality.parameter = 'NONE'  # No quality analysis
# wn.options.quality.diffusivity = 1.0
# wn.options.quality.tolerance = 0.01

# # Plot the network
# wntr.graphics.plot_network(wn, title="Florence Network", node_labels=True, filename="florence_network.png")

# # Print network summary
# print(f"Network summary: {wn.describe()}")

# # Write the INP file
# # wn.write_inpfile('./Networks/Ismail/ismail.inp', units="LPS")
# wntr.epanet.io.InpFile().write(filename="./Networks/Ismail/ismail.inp", 
#                                wn=wn, 
#                                units="LPS",
#                                version=2.2)