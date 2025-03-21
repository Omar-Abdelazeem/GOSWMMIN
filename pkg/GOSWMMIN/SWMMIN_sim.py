import wntr             # For EPANET file reading
import numpy as np      
import pandas as pd
import re
import math
import pathlib 
import pyswmm
import tqdm
import matplotlib.pyplot as plt
import matplotlib
import scipy.integrate as sci

class SWMMIN_sim:
    def __init__(self, input_file):
        """
        Initialize a SWMMIN simulation object based on an EPANET input file

        Arguments:
        input_file (str): the path to the EPANET input file  
        """

        self.filepath = pathlib.Path(input_file)
        name_only=self.filepath.name
        print("Selected File: ",name_only)
        self.name_only = name_only


    def Convert_to_SWMMIN(self, supply_duration, minimum_pressure, 
                        desired_pressure, pdw_exponent, n_days = 1, length_to_diameter=30, 
                        adaptive_disc=False, maximum_xdelta=None, solution_speed = 100, timestep =None, leak_fraction=0.1,
                        q_des=None, tank_areas=None, tank_heights=1, consum_pattern=None, pdw_variable='PRESSURE'):
        """
        This function reads in an EPANET input file and converts it to a SWMMIN input file.
        The function reads in the input file line by line and converts the EPANET input file to a SWMMIN input file. 
        The function writes the SWMMIN input file to the same directory as the EPANET input fole returns the path to the SWMMIN input file.

        Arguments:
        supply_duration (float or str): the duration of the supply period in hours)
        minimum_pressure (float or str): the minimum pressure for PDW in meters. Either a float for a uniform value for all nodes or a string path to a CSV.
        desired_pressure (float or str): the desired pressure for PDW in meters
        pdw_exponent (float or str): the exponent for the PDW equation
        n_days (int): the number of days to simulate. default is 1.
        length_to_diameter (float): the ratio of maximum pipe length to pipe diameter. Default and recommended value is 30.
        adaptive_disc (Boolean): whether to use adaptive discretization. Default is False.
        delta_x_max (float or str): the maximum allowable distance between nodes in meters. Overrides the length_to_diameter  and the adaptive_disc parameters.
        solution_speed (float): the solution speed in m/s. Default is 100, recommended range is 50 to 200 m/s.
        timestep (float): the time step in seconds. Overrides solution speed based timestep. Default is None.
        leak_fraction (float): fraction of nodal demand to be assigned as leakage. Default is 0.1.
        q_des (float or str): the desired flow rate for PDW in m^3/s. Either a float for a uniform value for all nodes or a string path to a CSV.
        tank_areas (None or str): if None, tank areas are calculated from Q_demand in the EPANET input file. If a string, the path to a CSV with tank areas. Default is None.
        tank_heights (float or str): height of user storage tanks. Either a float for a uniform value for all tanks or a string path to a CSV. Default is 0.
        consum_pattern (None or str): the path to a CSV file with the consumption pattern. Default is None for built-in generic demand pattern.
        pdw_variable (str): Either DEPTH or PRESSURE. Default is PRESSURE.


        Returns: path to the SWMMIN input file
        """

        assert isinstance(supply_duration, (float, int)), "Supply duration must be a float or an integer"
        if isinstance(minimum_pressure, str):
            assert pathlib.Path(minimum_pressure).exists(), "Minimum pressure file does not exist"
            assert isinstance(desired_pressure, str), "Minimum and desired pressure + pdw exponent must match: either all uniform or all node-specific with corresponding CSVs"
            assert isinstance(pdw_exponent, str), "Minimum and desired pressure + pdw exponent must match: either all uniform or all node-specific with corresponding CSVs"
            self.pdw_flag = True
            minimum_pressure = pd.read_csv(minimum_pressure, header = None, names = ['ID', 'Minimum Pressure'])['Minimum Pressure'].to_list()
            desired_pressure = pd.read_csv(desired_pressure, header = None, names = ['ID', 'Desired Pressure'])['Desired Pressure'].to_list()
            pdw_exponent = pd.read_csv(pdw_exponent, header = None, names = ['ID', 'Exponent'])['Exponent'].to_list()
        else: self.pdw_flag = False
        
        # Flags to indicate input types:
        if q_des:
            self.q_des_flag = True
            

        else: self.q_des_flag = False

        if tank_areas:
            self.tank_areas_flag = True
            tank_areas = pd.read_csv(tank_areas,header=None, names = ['ID', 'Area'])['Area'].to_list()
        else: self.tank_areas_flag = False

        if isinstance(tank_heights, str):
            self.tank_heights_flag = True
            tank_heights = pd.read_csv(tank_heights,header=None, names = ['ID', 'Height'])['Height'].to_list()
        else: self.tank_heights_flag = False

        if consum_pattern:
            self.consum_pattern_flag = True
            consum_pattern = pd.read_csv(consum_pattern,header=None, names = ['Time','Pattern'])['Pattern'].to_list()
            assert len(consum_pattern) == 24, "Consumption pattern must be 24 hours long"
        else: consum_pattern = [0.8,0.7,0.6,0.5,0.5,0.5,0.6,0.8,1.2,1.3,1.2,1.2,1.2,1.2,1.2,1.2,1.1,1.1,1.1,1.2,1.3,1.3,1.1,1]



        self.supply_duration = supply_duration
        self.minimum_pressure = minimum_pressure
        self.desired_pressure = desired_pressure
        self.pdw_exponent = pdw_exponent
        self.n_days = n_days - 1
        self.length_to_diameter = length_to_diameter
        self.adaptive_disc = adaptive_disc
        self.maximum_xdelta = maximum_xdelta
        self.solution_speed = solution_speed
        self.timestep = timestep
        self.leak_fraction = leak_fraction
        self.q_des = q_des
        self.tank_areas = tank_areas
        self.tank_heights = tank_heights
        self.consum_pattern = consum_pattern
        self.pdw_variable = pdw_variable

        


        # Read the EPANET input file and extract the network model
        self.__parse_epanet_input__()

        # Match pipes concentrically
        self.__match_pipes_concentric__()
        
        # Discretize the pipes
        self.__discretize__()

        self.__write_junctions__()

        self.__write_outfalls__()

        self.__write_storage_nodes__()

        self.__write_pipes__()

        self.__write_outlets__()

        self.__write_xsections__()

        self.__write_curves__()

        self.__write_controls__()

        self.__write_coords__()

        self.__write_file__()

        return self.swmm_file


    def __parse_epanet_input__(self):
        '''
        This function reads in the EPANET input file and extracts the network model and relevant information.
        '''
        filepath=self.filepath
        supply_duration = self.supply_duration
       
        demand_nodes=[]       # For storing list of nodes that have non-zero demands
        base_demands=[]    # For storing demand rates desired by each node for desired volume calculations
        elevations=[]         # For storing elevations of demand nodes
        coords=dict()         # For storing coordinates corresponding to each node as a tuple with the id as key
        all_nodes=[]          # For storing list of node ids of all nodes
        all_elevations=[]     # For storing elevations of all nodes

        # Creates a network model object using EPANET .inp file
        network=wntr.network.WaterNetworkModel(filepath)

        # Iterates over the junction list in the Network object
        for node in network.junctions():
            all_nodes.append(node[1].name)
            all_elevations.append(node[1].elevation)
            coords[node[1].name]=node[1].coordinates
            # For all nodes that have non-zero demands
            if node[1].base_demand != 0:
                # Record node ID (name), desired demand (base_demand) in CMS, elevations, x and y coordinates
                demand_nodes.append(node[1].name)
                base_demands.append(node[1].base_demand)
                elevations.append(node[1].elevation)

        # If a desired flow rate CSV is provided
        if self.q_des_flag == True:
            #replace desired demands with input Q_des values
            q_des = pd.read_csv(self.q_des,header=None)
            q_des.columns = ['ID','Desired Flow']
            desired_demands = q_des['Desired Flow'].to_list()  
        else: desired_demands = [demand * 24 / self.supply_duration for demand in base_demands]     

        conduit_ids= []       # To store IDs of the original pipes in the EPANET file
        conduit_from= []      # To store the origin node for each pipe
        conduit_to= []        # To store the destination node for each pipe
        conduit_lengths= []   # To store pipe lengths
        conduit_diameters= [] # To store pipe diameters

        # Loop over each link in the EPANET model
        for link in network.pipes():

            # Extract and store each of the aforementioned properties
            conduit_ids.append(link[1].name)
            conduit_from.append(link[1].start_node_name)
            conduit_to.append(link[1].end_node_name)
            conduit_lengths.append(link[1].length)
            conduit_diameters.append(link[1].diameter)

        pump_ids= []       # To store IDs of the original pumps in the EPANET file
        pump_from= []      # To store the origin node for each pump
        pump_to= []        # To store the destination node for each pump
        pump_curves = []   # To store the name of the pump's curve


        # Loop over each pump in the EPANET model
        for pump in network.pumps():

            # Extract and store each of the aforementioned properties
            pump_ids.append(pump[1].name)
            pump_from.append(pump[1].start_node_name)
            pump_to.append(pump[1].end_node_name)
            pump_curves.append(pump[1]._pump_curve_name)


        reservoir_ids=[]      # To store the source reservoirs' IDs
        reservoir_heads={}    # To store the total head of each reservoir indexed by ID
        reservoir_coords={}   # To store the coordinates as tuple (x,y) indexed by ID

        # Loops over each reservoir
        for reservoir in network.reservoirs():
            reservoir_ids.append(reservoir[1].name)
            reservoir_heads[reservoir_ids[-1]]=reservoir[1].base_head
            reservoir_coords[reservoir_ids[-1]]=reservoir[1].coordinates
        reservoir_elevations={reservoir:reservoir_heads[reservoir]-30 for reservoir in reservoir_heads}
        pcurves = {}
        for curve in network.curves():
            pcurves[curve[0]]= curve[1].points

        # Get the supply duration in minutes (/60) as an integer
        supply_duration_EPANET=int(network.options.time.duration/60)
        supply_hh=str(supply_duration//1)     # The hour value of the supply duration (quotient of total supply in minutes/ 60)
        supply_mm=str(supply_duration%1 * 60)      # The minute value of the supply duration (remainder)

        # Corrects the formatting of the HH:MM by adding a 0 if it is a single digit: if minutes =4 -> 04
        if len(supply_mm)<2:
            supply_mm='0'+supply_mm
        if len(supply_hh)<2:
            supply_hh='0'+supply_hh
        
        # Dataframe aggregating all node information gathered from the EPANET file
        junctions = pd.DataFrame(zip(all_nodes,all_elevations,coords.values()),columns=["ID","Elevation","Coordinates"])
        # Set the junction ID as the index of the Dataframe
        junctions.set_index("ID",inplace=True)

        # Dataframe aggregating all conduit information gathered from the EPANET file
        conduits=pd.DataFrame(zip(conduit_ids,conduit_from,conduit_to,conduit_lengths,conduit_diameters),columns=["ID","from node","to node","Length","diameter"])
        # Set the conduit ID as the index
        conduits.set_index("ID",inplace=True)

        reservoirs = pd.DataFrame(zip(reservoir_ids, reservoir_coords.values(), reservoir_elevations.values(), reservoir_heads.values()), columns=["ID","coords","Elevation","Head"])
        reservoirs.set_index("ID",inplace=True)

        self.junctions = junctions
        self.conduits = conduits
        self.reservoirs = reservoirs
        self.demand_nodes = demand_nodes
        self.elevations = elevations
        self.base_demands = base_demands
        self.desired_demands = desired_demands
        self.network = network

        return


    def __match_pipes_concentric__(self):
        '''
        This function matches pipes concentrically based on the diameter of the biggest pipe connected to each node.
        '''
        junctions = self.junctions
        conduits= self.conduits

        concentric=True #If this is False, Pipes will be matched using inverts

        connectivity=pd.DataFrame(index=junctions.index, columns=["US","DS"])

        if concentric:
            # For each node in the network
            for junction in junctions.index:
                upst_connected=[] #List of pipes whose upstream end connects to this node
                downst_connected=[] #List of pipes whose downstream end connects to this node
                # Loop over all conduits
                for conduit in conduits.index:
                    # if conduit starts at this node
                    if conduits.at[conduit,"from node"]==junction:
                        # append to upstream connected list
                        upst_connected.append(conduit)
                    # if conduit ends at this node
                    elif conduits.at[conduit,"to node"]==junction:
                        # append to downstream connected list
                        downst_connected.append(conduit)
                # create a list of all connected conduits
                connected_conduits=upst_connected+downst_connected
                # Get the diameter of all connected conduits
                connected_diameters=[conduits.at[conduit,"diameter"] for conduit in connected_conduits]
                # Find the diameter of the biggest pipe connected to this node
                max_diameter=max(connected_diameters)

                for conduit in upst_connected:
                    # If pipe is not the biggest pipe
                    if conduits.at[conduit,"diameter"]<max_diameter:
                        # offset this pipe by half the difference in their diameters
                        conduits.at[conduit,"InOffset"]=(max_diameter-conduits.at[conduit,"diameter"])/2
                    # if it is the biggest pipe do nothing
                    else: conduits.at[conduit,"InOffset"]=0

                for conduit in downst_connected:
                    # if not the biggest connected pipe to this node
                    if conduits.at[conduit,"diameter"]<max_diameter:
                        # offset this pipe by half the difference in their diameters
                        conduits.at[conduit,"OutOffset"]=(max_diameter-conduits.at[conduit,"diameter"])/2
                    else: conduits.at[conduit,"OutOffset"]=0
            
                connectivity.at[junction,"US"]=upst_connected
                connectivity.at[junction,"DS"]=downst_connected
                connectivity.at[junction,"Max D"]=max_diameter
        
        self.conduits = conduits
        self.junctions = junctions
        self.connectivity = connectivity
        return self

    def __discretize__(self):
        '''
        This function discretizes pipes based on the input to Convert_to_SWMMIN. Adaptive discretization is used if the adaptive_disc flag is set to True.
        If maximum_xdelta is set, it overrides the length_to_diameter parameter. Otherwise, the length_to_diameter parameter is used based on the largest pipe diameter.
        '''
        conduits = self.conduits
        reservoirs = self.reservoirs
        junctions = self.junctions

        reservoir_pipes=[]

        # If the maximum allowable distance between nodes is not set, calculate it based on the largest pipe diameter and len to diameter ratio
        if   self.maximum_xdelta is None:
            self.maximum_xdelta=conduits["diameter"].max() * self.length_to_diameter

        # Loop over each conduit in the original file
        for conduit in conduits.index:

            length=conduits["Length"][conduit]  #Stores the length of the current conduit for shorthand
            diameter = conduits["diameter"][conduit]
            # if adaptive discretization is set to True, set the number of parts adaptively for each conduits
            if self.adaptive_disc:
                n_parts = np.max(round(length / diameter / self.length_to_diameter) - 1, 0)
            # otherwise, use the maximum_xdelta parameter to set the number of parts
            else: n_parts = math.ceil(length / self.maximum_xdelta)

            # If the conduit is bigger than the maximum allowable length (delta x), we will break it down into smaller pipes
            if n_parts>1:
                # Calculate the length of each part 
                part_length=length/n_parts
                # Start node ID (for shorthand)
                start_node=conduits["from node"][conduit]
                # Inlet Offset for the Current Conduit
                in_offset=conduits["InOffset"][conduit]
                # End node ID (for shorthand)
                end_node=conduits["to node"][conduit]
                #Outlet Offset for the Current Conduit
                out_offset=conduits["OutOffset"][conduit]
                # If the start node is a reservoir
                if start_node in reservoirs.index:
                    # MAke the start elevation the same as the end but add 1 (since reservoirs don't have ground elevation in EPANET)
                    start_elevation=junctions.at[end_node,"Elevation"]+1
                    reservoirs.at[start_node,"Elevation"]=start_elevation+1
                # Otherwise make the start elevation equal to the elevation of the start node + the offset of the pipe
                else: start_elevation=junctions.at[start_node,"Elevation"]+in_offset
                
                # If the end node is a reservoir
                if end_node in reservoirs.index:
                    # Make the end elevation the same as the start but subtract 1 (since reservoirs don't have ground elevation in EPANET)
                    end_elevation=start_elevation-1
                # Make the end elevation equal to the elevation of the end node + the offset of the pipe
                else: end_elevation=junctions.at[end_node,"Elevation"]+out_offset
                # Calculate the uniform drop (or rise) in elevation for all the intermediate nodes about to be created when this pipe is broken into several smaller ones
                unit_elev_diff=(end_elevation-start_elevation)/n_parts

                # if the starting node is a reservoir
                if start_node in reservoirs.index:
                    # Get coordinates from reservoir data
                    start_x=reservoirs.at[start_node,"coords"][0]
                    start_y=reservoirs.at[start_node, "coords"][1]
                else:
                    # Get the coordinates from the junction data
                    start_x=junctions.at[start_node,"Coordinates"][0]
                    start_y=junctions.at[start_node,"Coordinates"][1]
                
                # If the end node is a reservoir
                if end_node in reservoirs.index:
                    # Get the coordinates from the reservoir data
                    end_x=reservoirs.at[end_node, "coords"][0]
                    end_y=reservoirs.at[end_node, "coords"][1]
                else:
                    # Get them from the junctions data
                    end_x=junctions.at[end_node,"Coordinates"][0]
                    end_y=junctions.at[end_node,"Coordinates"][1]
                    
                # Calculate the unit difference in x and y coordinates for this pipe and its segments
                unit_x_diff=(end_x-start_x)/n_parts
                unit_y_diff=(end_y-start_y)/n_parts


        # THIS LOOP GENERATES THE SMALLER PIPES TO REPLACE THE ORIGINAL LONG PIPE
                # For each part to be created
                for part in np.arange(1,n_parts+1):

                    # CREATING THE LINKS
                    # Create the ID for the new smaller pipe as OriginPipeID-PartNumber
                    new_id=conduit+"-"+str(part)
                    # Set the new pipe's diameter equal to the original one
                    conduits.at[new_id,"diameter"]=conduits["diameter"][conduit]
                    # Set the start node as OriginStartNode-NewNodeNumber-OriginEndNode  as in the first intermediate nodes between node 13 and 14 will be named 13-1-14
                    conduits.at[new_id,"from node"]=start_node+"-"+str(part-1)+"-"+end_node

                    # Set Inlet and Outlet Offsets to Zero by default
                    conduits.at[new_id,"InOffset"]=0
                    conduits.at[new_id,"OutOffset"]=0

                    # if this is the first part, use the original start node 
                    if part==1:
                        conduits.at[new_id,"from node"]=start_node
                        # Set the first inlet offset to concentrically match the upstream pipe
                        conduits.at[new_id,"InOffset"]=in_offset
                        
                    # Set the end node as OriginStartNode-NewNodeNumber+1-OriginEndNode  as in the second intermediate nodes between node 13 and 14 will be named 13-2-14
                    conduits.at[new_id,"to node"]=start_node+"-"+str(part)+"-"+end_node
                    # If this is the last part, use the original end node as the end node
                    if part==n_parts:
                        conduits.at[new_id,"to node"]=end_node
                        # Set the last outlet offset to concentrically match the downstream pipe
                        conduits.at[new_id,"OutOffset"]=out_offset
                    # Set the new pipe's length to the length of each part
                    conduits.at[new_id,"Length"]=part_length

                    # if this is NOT the last part (as the last pipe segment joins a pre-existing node and does not need a node to be created)
                    if part<n_parts:
                        # Create a new node at the end of this pipe segment whose elevation is translated from the start elevation using the unit slope and the part number
                        junctions.at[conduits.at[new_id,"to node"],"Elevation"]=start_elevation+part*unit_elev_diff
                        # Calculate the coordinates for the new node using the unit difference in x and y coordinates
                        junctions.at[conduits.at[new_id,"to node"],"Coordinates"]=(start_x+part*unit_x_diff,start_y+part*unit_y_diff)

                    if conduits.at[new_id,"from node"] in reservoirs.index or conduits.at[new_id,"to node"] in reservoirs.index:
                        reservoir_pipes.append(new_id)
                    
                # After writing the new smaller pipes, delete the original pipe (since it is now redundant)
                conduits.drop(conduit,inplace=True)

        conduits[["InOffset","OutOffset"]]=conduits[["InOffset","OutOffset"]].fillna(0)
        self.mean_xdelta = conduits['Length'].mean()
        self.conduits = conduits
        self.junctions=junctions
        self.reservoir_pipes = reservoir_pipes

        return self
    
    def __write_junctions__(self):
        '''
        This function writes the junctions section of the SWMM input file
        '''
        junctions = self.junctions

        MaxDepth=[0]*len(junctions)
        InitDepth=MaxDepth
        SurDepth=[100] * len(junctions)  # High value to prevent surcharging
        Aponded=InitDepth

        # Creates dataframe with each row representing one line from the junctions section
        junctions_section=pd.DataFrame(list(zip(junctions.index,junctions["Elevation"],MaxDepth,InitDepth,SurDepth,Aponded)))
        # Converts the dataframe into a list of lines in the junctions section
        junctions_section=junctions_section.to_string(header=False,index=False,col_space=10).splitlines()
        # adds a new line character to the end of each line in the section
        junctions_section=[line+'\n' for line in junctions_section]

        self.junctions_section = junctions_section
        return self

    def __write_outfalls__(self):
        '''
        This function writes the outfalls section of the SWMM input file
        '''
        demand_nodes = self.demand_nodes
        elevations = self.elevations

        # Creates the IDs for the outfalls
        outfall_ids=["Outfall"+ str(node) for node in demand_nodes]
        # Set the elevations of the outfalls to be the same as the elevations of the demand nodes
        outfall_elevations=elevations
        # Free outfall types for all outfalls
        outfall_types=['FREE' for i in demand_nodes]
        # No stage data for the outfalls
        stage_data=["    " for i in demand_nodes]
        # No gated outfalls
        outfall_gated=["NO" for i in demand_nodes]
        # Creates dataframe with each row representing one line from the outfalls section
        outfall_section=pd.DataFrame(zip(outfall_ids,outfall_elevations,outfall_types,stage_data,outfall_gated))

        leak_outfall_id=["L_Outfall"+str(node) for node in demand_nodes]

        leak_outfalls=pd.DataFrame(zip(leak_outfall_id,outfall_elevations,outfall_types,stage_data,outfall_gated))
        outfall_section=pd.concat([outfall_section,leak_outfalls])
        # Converts the dataframe into a list of lines in the outfalls section
        outfall_section=outfall_section.to_string(header=False,index=False,col_space=10).splitlines()
        # adds a new line character to the end of each line in the section
        outfall_section=[line+'\n' for line in outfall_section]
        self.outfall_section = outfall_section
        self.outfall_ids = outfall_ids
        self.leak_outfall_id = leak_outfall_id
        return self

    def __write_storage_nodes__(self):
        '''
        This function writes the storage nodes section of the SWMM input file. Including reservoirs and user tanks'''

        demand_nodes = self.demand_nodes
        tank_heights = self.tank_heights
        tank_areas = self.tank_areas
        reservoirs = self.reservoirs
        conduits = self.conduits

        # Set Storage IDS (correspond to original DN ID)
        storage_ids=["StorageforNode"+id for id in demand_nodes]
        # demand_volumes = pd.DataFrame({'ID':demand_nodes,'Volume':storage_areas})
        # demand_volumes.to_csv(filepath.parent/pathlib.Path(filepath.name+'_DemandVolumes.csv'))
        self.demand_volumes = [demand * 60 * 24 for demand in self.base_demands]

        # OVERRIDE OPTION: Input Tank Heights using the CSV
        if self.tank_heights_flag == True:
            storage_MaxDepth=tank_heights
        else: storage_MaxDepth=[tank_heights]*len(storage_ids)
        
        # OVERRIDE OPTION: Input Tank Areas using the CSV
        if self.tank_areas_flag == True:
            storage_areas = self.tank_areas
        else: storage_areas=[demand * 60 * 24 /height for demand,height in zip(self.base_demands,storage_MaxDepth)]            

        # Set storage invert elevations = elevation of DN + Hmin
        if self.pdw_flag == True:
            storage_elevations=[elevation+minimum for elevation, minimum in zip(self.elevations,self.minimum_pressure)]
        else: storage_elevations=[elevation+self.minimum_pressure for elevation in self.elevations]

        # All tanks start empty
        storage_InitDepth=[0]*len(storage_ids)
        # All tanks are defined using the functional input
        storage_shape=["FUNCTIONAL"]*len(storage_ids)
        zeroes=['0']*len(storage_ids)    # for the other curve parameters which are not required for circular tanks
        # No surcharge depth for the tanks
        storage_SurDepth=[0]*len(storage_ids)
        # No evaporation from storage tanks
        storage_fevap=[0]*len(storage_ids)

        storage_units=pd.DataFrame(zip(storage_ids,storage_elevations,storage_MaxDepth,storage_InitDepth,storage_shape,zeroes,zeroes,storage_areas,storage_SurDepth,storage_fevap))

        # Creating the Storage nodes that represent the supply reservoirs
        # Rename reservoirs
        reservoir_ids_new=["Reservoir-"+str(id) for id in self.reservoirs.index]
        # update references to reservoirs in conduits
        for i in range(0,len(reservoirs.index)):
            conduits.replace(to_replace=reservoirs.index[i],value=reservoir_ids_new[i])
        # set reservoir elevations
        reservoir_elevations=reservoirs["Elevation"]
        # set max depth of reservoirs
        MaxDepth=[max(100,max(reservoirs["Elevation"])+40)]*len(reservoirs.index)
        # Set initial depth of reservoirs
        InitDepth=[head-elevation for head, elevation in zip(reservoirs["Head"],reservoirs["Elevation"])]
        reservoir_shape=["TABULAR"]*len(reservoirs.index)
        reservoir_curves=["Source"+id for id in reservoirs.index]
        blanks=['    ']*len(storage_ids)  
        reservoir_SurDepth=[0]*len(reservoirs.index)
        reservoir_psi=[0]*len(reservoirs.index)

        storage_section=pd.DataFrame(zip(reservoir_ids_new,reservoir_elevations,MaxDepth,InitDepth,reservoir_shape,reservoir_curves,blanks,blanks,reservoir_SurDepth,reservoir_psi))
        storage_section= pd.concat([storage_section,storage_units])
        storage_section=storage_section.to_string(header=False,index=False,col_space=10).splitlines()
        storage_section=[line+'\n' for line in storage_section]

        self.storage_section = storage_section
        self.reservoir_ids_new = reservoir_ids_new
        self.storage_ids = storage_ids
        self.storage_areas = storage_areas
        self.reservoir_curves = reservoir_curves
        return self
    
    def __write_pipes__(self):
        '''
        This function writes the conduits section of the SWMM input file'''

        conduits = self.conduits
        reservoirs = self.reservoirs

        # Create the IDs for the conduits
        conduit_ids=["P-"+str(conduit) for conduit in conduits.index]
        conduits.index = conduit_ids
        # Set the Manning's roughness coefficient for a typical plastic pipe
        roughness=[0.011]*len(conduits)
        # Other inputs that are set to 0
        conduit_zeros=[0]*len(conduits)
        # Replace the reservoir IDs with the new reservoir IDs
        conduits["from node"] = conduits["from node"].replace(reservoirs.index,self.reservoir_ids_new)
        conduits["to node"] = conduits["to node"].replace(reservoirs.index,self.reservoir_ids_new)


        conduits_section=pd.DataFrame(zip(conduit_ids,conduits["from node"],conduits["to node"],conduits["Length"],roughness,conduits["InOffset"],conduits["OutOffset"],conduit_zeros,conduit_zeros))
        conduits_section=conduits_section.to_string(header=False,index=False,col_space=10).splitlines()
        conduits_section=[line+'\n' for line in conduits_section]

        self.conduits_section = conduits_section
        return self
    
    def __write_outlets__(self):
        '''
        This function writes the outlets section of the SWMM input file
        '''
        demand_nodes = self.demand_nodes
        storage_ids = self.storage_ids
        connectivity = self.connectivity
        outfall_ids = self.outfall_ids

        # Withdrawal Outlets
        # Set IDs for the outlets
        outlet_ids = ["Outlet"+id for id in demand_nodes]
        # Set the origin nodes for the outlets
        outlet_from = demand_nodes[:]
        # Set the destination nodes for the outlets
        outlet_to = storage_ids [:]
        # Set the offset for the outlets
        outlet_offset=[connectivity.at[node,"Max D"]/2 for node in demand_nodes]
        # Funcitonal input for all withdrawal outlets
        outlet_type=["FUNCTIONAL/"+self.pdw_variable]*len(outlet_ids)

        # If the PDW flag is set, use the input CSVs to set the coefficients and exponents
        if self.pdw_flag==True:
            outlet_coeff = [demand * 1000 / np.power((desired - minimum), exponent) for demand, desired, minimum, exponent 
                            in zip(self.desired_demands, self.desired_pressure, self.minimum_pressure, self.pdw_exponent)]
            outlet_expon = self.pdw_exponent
        # otherwise use the constant values
        else:
            outlet_coeff=[demand*1000/np.sqrt(self.desired_pressure - self.minimum_pressure) for demand in self.desired_demands]  
            outlet_expon=[str(self.pdw_exponent)]*len(outlet_ids)     
        outlet_gated=["YES"]*len(outlet_ids)

        # Consumption Outlets
        outdemand_ids = ["ConsumptionOutlet"+id for id in demand_nodes]
        outdemand_from = storage_ids [:]
        outdemand_to = outfall_ids [:]
        outdemand_offset=[0]*len(outlet_ids)
        outdemand_type=["TABULAR/DEPTH"]*len(outlet_ids)
        outdemand_coeff=["Demand"+ id for id in demand_nodes]
        outdemand_expon=["     "]*len(outlet_ids)
        outdemand_gated=["YES"]*len(outlet_ids)

        outlets=pd.DataFrame(list(zip(outlet_ids,outlet_from,outlet_to,outlet_offset,outlet_type,outlet_coeff,outlet_expon,outlet_gated)))
        outletdemand=pd.DataFrame(list(zip(outdemand_ids,outdemand_from,outdemand_to,outdemand_offset,outdemand_type,outdemand_coeff,outdemand_expon,outdemand_gated)))
        self.outletwithdraw = outlets
        self.outletdemand = outletdemand

        # Leakage Outlets
        leak_ids=["LeakforNode"+str(node) for node in demand_nodes]
        leak_from=demand_nodes[:]
        leak_to=self.leak_outfall_id[:]
        leak_offset=[connectivity.at[node,"Max D"]/2 for node in demand_nodes]
        leak_type=["FUNCTIONAL/DEPTH"]*len(leak_ids)
        if self.pdw_flag:
            leak_coeff = [self.leak_fraction * demand * 1000 / np.power((desired - minimum), exponent) for demand, desired, minimum, exponent 
                          in zip(self.desired_demands, self.desired_pressure, self.minimum_pressure, self.pdw_exponent)]
            leak_expon = self.pdw_exponent
        else:
            leak_coeff=[demand*1000/np.sqrt(self.desired_pressure-self.minimum_pressure)*self.leak_fraction for demand in self.desired_demands]
            leak_expon=[str(self.pdw_exponent)]*len(outlet_ids)     
        leak_gated=["YES"]*len(leak_ids)

        leak_outlets=pd.DataFrame(zip(leak_ids,leak_from,leak_to,leak_offset,leak_type,leak_coeff,leak_expon,leak_gated))

        outlet_section=pd.concat([outlets,outletdemand,leak_outlets])
        outlet_section=outlet_section.to_string(header=False,index=False,col_space=10).splitlines()
        outlet_section=[line+'\n' for line in outlet_section]

        self.outlet_section = outlet_section
        return self
    
    def __write_xsections__(self):
        conduits = self.conduits
        shape=["FORCE_MAIN"]*len(conduits.index)
        hwcoeffs=[130]*len(shape)
        geom3=[0]*len(shape)
        geom4=geom3
        nbarrels=[1]*len(shape)

        xsections_section=pd.DataFrame(zip(conduits.index,shape,conduits["diameter"],hwcoeffs,geom3,geom4,nbarrels))
        xsections_section=xsections_section.to_string(header=False,index=False, col_space=10).splitlines()
        xsections_section=[line+'\n' for line in xsections_section]

        self.xsections_section = xsections_section
        return self
    
    def __write_curves__(self):

        curves_name=[]
        curves_type=[]
        curves_x=[]
        curves_y=[]


        constant_demands=[demand*float(self.supply_duration//1)/24 for demand in self.desired_demands]
        for i,j in zip(self.demand_nodes,constant_demands):
            curves_name+=["Demand"+str(i),"Demand"+str(i),';']
            curves_type+=["Rating"," ","  "]
            curves_x+=[0,0.01," "]
            curves_y+=[0,j*1000," "]

        # Reservoir Storage Curves
        # reservoir_volume=sum(storage_areas)
        pipe_vols = self.conduits['Length'] * self.conduits['diameter']**2 * math.pi / 4
        pipe_vols = pipe_vols.sum()
        reservoir_volume = (sum(self.storage_areas) + pipe_vols) * (self.n_days + 2) * 2
        for curve,head,elevation in zip(self.reservoir_curves,self.reservoirs["Head"],self.reservoirs["Elevation"]):
            j=float(head)-float(elevation)
            for depth in [0,j-2,j-1,j]:
                curves_name.append(curve)
                curves_x.append(depth)
                if depth<=j-2:
                    curves_y.append(1)
                else: curves_y.append(reservoir_volume)
                if depth==0:
                    curves_type.append("Storage")
                else: 
                    curves_type.append(" ")
            curves_name.append(";")
            curves_type.append(" ")
            curves_x.append(" ")
            curves_y.append(" ")

        # for curve in pcurves:
        #     if len(pcurves[curve]) == 1:
        #         design_flow = pcurves[curve][0][0] * 1000
        #         design_head = pcurves[curve][0][1]

        #         shutoff_head = 1.33 * design_head 
        #         max_flow = 2 * design_flow
        #     # Fit an H(Q) = H_max - B * Q^2 Function
        #         B = (shutoff_head - design_head) / design_flow ** 2
        #         flows = []
        #         heads = []
        #         for flow in np.arange(0, max_flow + max_flow/20, max_flow/20):
        #             flows.append(flow)
        #             heads.append(shutoff_head - B * flow **2)
        #     for flow, head in zip(flows, heads):
        #         curves_name.append("Pump_Curve_"+ curve)
        #         if flow  == flows[0]:
        #             curves_type.append("Pump3")
        #         else: curves_type.append("  ")
        #         curves_x.append(flow)
        #         curves_y.append(head)
        curves=pd.DataFrame(list(zip(curves_name,curves_type,curves_x,curves_y)))
        curves_section=curves.to_string(header=False,index=False,col_space=10)

        self.curves_section = curves_section
        return self
    
    def __write_controls__(self):
        # Step for control curve generation. The smaller the step, the smoother the curve
        step=0.002
        # Parameters for Mauro de Marchis et al. (2015)'s float valve emitter law
        m, n= 5, 10
        if isinstance(self.tank_heights, float) or isinstance(self.tank_heights, int):
            h_min, h_max = 0, self.tank_heights

        # Intializing string for control curves and rules
        control_curves=""
        control_rules=""

        # Diurnal demand pattern defined as a list of multipliers for the base demand
        pattern = self.consum_pattern
        times=np.arange(0,int(self.supply_duration//1)*50,1)
        timesrs_pat=""
        for time in np.arange(0,len(times)):
            if time >=len(pattern):
                days=math.floor(time/len(pattern))
                time_24=time-days*len(pattern)
            else: time_24=time
            timesrs_pat+="Pattern\t"+str(time)+"\t"+str(pattern[time_24])+"\n"

        i = 0
        for outlet in self.outletwithdraw.iloc[:,0]:
            out_name=outlet
            if self.tank_heights_flag == True:
                h_max = self.tank_heights[i]
            h_min = 0
            storage_name=self.outletwithdraw.loc[self.outletwithdraw[0]==outlet,2].iloc[0]
            control_curves+="Control"+out_name+"\tControl"
            for height in np.arange(h_min,h_max+step,step):
                red_coeff=np.tanh(m*(h_max-height)/(h_max-h_min))*np.tanh(n*(h_max-height)/(h_max-h_min))
                if height==h_min:
                    control_curves+="\t\t"+str(round(height,6))+"\t"+str(round(red_coeff,6))+'\n'
                else:
                    control_curves+="Control"+out_name+"\t\t"+"\t"+str(round(height,6))+"\t"+str(round(red_coeff,6))+'\n'
            control_curves+=";\n"

            control_rules+="Rule "+out_name+"\n"
            control_rules+="IF NODE "+storage_name+" DEPTH >= "+str(h_min)+"\n"
            control_rules+="THEN OUTLET "+out_name+" SETTING = CURVE "+"Control"+out_name+"\n\n"
            i+=1

        supply_hh = str(math.floor(self.supply_duration//1))
        supply_mm = str(math.floor(self.supply_duration%1 * 60))
        if len(supply_mm)<2:
            supply_mm='0'+supply_mm
        if len(supply_hh)<2:
            supply_hh='0'+supply_hh

        # Add rule to stop supply by iterating over reservoir adjacent pipes
        control_rules+="Rule STOPSUPPLY\n"
        control_rules+="IF SIMULATION CLOCKTIME > "+supply_hh+":"+supply_mm+"\n"
        reservoirs_list = self.network.reservoir_name_list
        flag = True
        for reservoir in reservoirs_list:
            pipes = self.network.get_links_for_node(reservoir)
            for pipe in pipes:
                if self.network.get_link(pipe).length > self.maximum_xdelta:
                    suffix = '-1'
                else: suffix = ''
                if flag:
                    control_rules+="THEN CONDUIT P-"+pipe+suffix+" STATUS = CLOSED\n"
                    flag = False
                else: control_rules+="AND CONDUIT P-"+pipe+suffix+" STATUS = CLOSED\n"
        control_rules+="\n"

        
        # Add rule to start supply by iterating over reservoir adjacent pipes
        control_rules+="Rule STARTSUPPLY\n"
        control_rules+="IF SIMULATION CLOCKTIME >= 0:00\n"
        control_rules += "AND SIMULATION CLOCKTIME < "+supply_hh+":"+supply_mm+"\n"
        flag = True
        for reservoir in reservoirs_list:
            pipes = self.network.get_links_for_node(reservoir)
            for pipe in pipes:
                if self.network.get_link(pipe).length > self.maximum_xdelta:
                    suffix = '-1'
                else: suffix = ''
                if flag:
                    control_rules+="THEN CONDUIT P-"+pipe+suffix+" STATUS = OPEN\n"
                    flag = False
                else: control_rules+="AND CONDUIT P-"+pipe+suffix+" STATUS = OPEN\n"
        control_rules+="\n"


        control_rules+="Rule Patterns\n"
        control_rules+="IF SIMULATION TIME > 0\n"
        flag=0
        for outlet in self.outletdemand.iloc[:,0]:
            if flag==0:
                control_rules+="THEN OUTLET "+outlet+" SETTING = TIMESERIES Pattern\n"
                flag=1
            else: control_rules+="AND OUTLET "+outlet+" SETTING = TIMESERIES Pattern\n"

        self.curves_section+="\n"+control_curves
        self.control_rules = control_rules
        self.timesrs_pat = timesrs_pat
        return self
    
    def __write_coords__(self):
        junctions = self.junctions

        coords_demand= { node: junctions.at[node, "Coordinates"] for node in self.demand_nodes}
        coords_ids=list(junctions.index)+self.reservoir_ids_new+self.storage_ids+self.outfall_ids+self.leak_outfall_id
        offset=2

        coords_x1=[coord[0] for coord in junctions["Coordinates"]]
        coords_x2=[coord[0] for coord in self.reservoirs["coords"]]
        coords_x3=[coord[0] +2*offset for coord in coords_demand.values()]
        coords_x4=[coord[0] + 3 * offset for coord in coords_demand.values()]
        coords_x5=[coord[0] + offset for coord in coords_demand.values()]
        coords_x=coords_x1+coords_x2+coords_x3+coords_x4+coords_x5

        coords_y1=[coord[1] for coord in junctions["Coordinates"]]
        coords_y2=[coord[1] for coord in self.reservoirs["coords"]]
        coords_y3=[coord[1] for coord in coords_demand.values()]
        coords_y4=[coord[1]+ offset for coord in coords_demand.values()]
        coords_y5=[coord[1] - 1.5 * offset for coord in coords_demand.values()]
        coords_y=coords_y1+coords_y2+coords_y3+coords_y4+coords_y5

        coordinate_section=pd.DataFrame(zip(coords_ids,coords_x,coords_y))
        coordinate_section=coordinate_section.to_string(header=False,index=False,col_space=10).splitlines()
        coordinate_section=[line+'\n' for line in coordinate_section]

        #Setting View Dimensions
        x_left=min(coords_x)
        x_right=max(coords_x)
        y_down=min(coords_y)
        y_up=max(coords_y)
        dimensions_line=str(x_left)+" "+str(y_down)+" "+str(x_right)+" "+str(y_up)+"\n"

        self.coordinate_section = coordinate_section
        self.dimensions_line = dimensions_line

        return self
    
    def __write_file__(self):
        # opens .inp file to read
        template=pathlib.Path("../Resources/Empty_SWMM_Template.inp")
        file=open(template,'r')
        lines=[]              # list to store all lines in the .inp file
        linecount=0           # Counter for the number of lines

        # Loops over each line in the input file 
        for line in file:
            if re.search("^START_DATE",line):
                start_date = line[11:].lstrip().split('/')
            if re.search("^END_TIME",line):
                end_time=linecount
            if re.search("^ROUTING_STEP",line):
                routing_step=linecount
            if re.search("^END_DATE", line):
                end_date = linecount
            if re.search("^DIMENSIONS",line):
                dimensions=linecount
            # Record the position of the phrase [JUNCTIONS] and add 3 to skip the header lines
            if re.search('\[JUNCTIONS\]',line):
                junctions_marker=linecount+3
            # Record the position of the phrase [OUTFALLS] and add 3 to skip the header lines
            if re.search('\[OUTFALLS\]',line):
                outfalls_marker=linecount+3
            # Record the position of the phrase [STORAGE] and add 3 to skip the header lines
            if re.search('\[STORAGE\]',line):
                storage_marker=linecount+3
            # Record the position of the phrase [CONDUITS] and add 3 to skip the header lines
            if re.search('\[CONDUITS\]',line):
                conduits_marker=linecount+3
            # Record the position of the phrase [OUTLETS] and add 3 to skip the header lines
            if re.search('\[PUMPS\]',line):
                pumps_marker=linecount+3
            if re.search('\[OUTLETS\]',line):
                outlets_marker=linecount+3
            # Record the position of the phrase [XSECTIONS] and add 3 to skip the header lines
            if re.search('\[XSECTIONS\]',line):
                xsections_marker=linecount+3
            if re.search('\[CONTROLS\]',line):
                controls_marker=linecount+1
            # Record the position of the phrase [CURVES] and add 3 to skip the header lines
            if re.search('\[CURVES\]',line):
                curves_marker=linecount+3
            if re.search('\[TIMESERIES\]',line):
                timesrs_marker=linecount+3
            # Record the position of the phrase [COORDINATES] and add 3 to skip the header lines
            if re.search('\[COORDINATES\]',line):
                coords_marker=linecount+3
            # Store all lines in a list
            lines.append(line)
            linecount+=1
        file.close()

        try:
            controls_marker
        except NameError:
            controls_marker=curves_marker-4
            control_rules="[CONTROLS]\n"+control_rules

        swmm_file=self.filepath.stem+'_SWMMIN_'+str(self.maximum_xdelta)+'m.inp'
        self.swmm_file = self.filepath.parent / swmm_file
        file=open(self.filepath.parent/swmm_file,'w')
        lines[end_time]="END_TIME             "+str(23)+":"+str(59)+":00\n"
        lines[end_date]="END_DATE             "+start_date[0]+"/"+str(min((int(start_date[1])+self.n_days),30))+"/"+start_date[2]
        if not self.timestep:
            self.timestep = self.mean_xdelta / self.solution_speed
        lines[routing_step]="ROUTING_STEP         "+str(self.timestep)+"\n"
        lines[dimensions]="DIMENSIONS "+self.dimensions_line
        lines[coords_marker:coords_marker]=self.coordinate_section
        lines[timesrs_marker:timesrs_marker]=self.timesrs_pat
        lines[curves_marker:curves_marker]=self.curves_section
        lines[controls_marker:controls_marker]=self.control_rules
        lines[xsections_marker:xsections_marker]=self.xsections_section
        lines[outlets_marker:outlets_marker]=self.outlet_section
        # if len(pump_section)>3:
        #     lines[pumps_marker:pumps_marker]=pump_section
        lines[conduits_marker:conduits_marker]=self.conduits_section
        lines[storage_marker:storage_marker]=self.storage_section
        lines[outfalls_marker:outfalls_marker]=self.outfall_section
        lines[junctions_marker:junctions_marker]=self.junctions_section


        # All lines added by this script are missing a new line character at the end, the conditional statements below add the new line character for these lines only and writes all lines to the file
        for line in lines:
            file.write(line)    
        file.close()
        return self
    
    def Run_SWMMIN(self, path = None):
        '''
        This function runs the SWMMIN simulation. 
        If no SWMM input file was generated using Convert_to_SWMMIN method, a path must be specified for the SWMMIN input file

        path: str: path to the SWMMIN input file to override the built-in path
        '''
        # Verify file path
        if path:
            self.swmm_file = path

        # Name the output file
        self.output_file = str(pathlib.Path(self.swmm_file).with_suffix('.out'))
        # Open the SWMM simulation
        sim=pyswmm.Simulation(str(self.swmm_file))

        # Retrieve demand node IDs
        demand_nodes = self.demand_nodes
        storage_ids = self.storage_ids
        w_outlet_ids = ["Outlet"+str(node) for node in demand_nodes]
        c_outlet_ids = ["ConsumptionOutlet"+str(node) for node in storage_ids]


        number_of_steps = int((sim.end_time - sim.start_time).seconds / self.timestep)
        # runs the simulation step by step
        with tqdm.tqdm(range(number_of_steps), desc = 'Running Simulation') as pbar:
            with sim as sim:
                system_routing = pyswmm.SystemStats(sim)
                for step in sim:
                    pbar.update(1)
                        # print('Current Simulation Time is >> ',sim.current_time,", ",round(sim.percent_complete*100,1),"% Complete")
                    pass
                sim._model.swmm_end()
                print("Continuity Error: ",sim.flow_routing_error,"%\n")

        sim.close()
        return self

    def get_pressures(self, specific_nodes = None):
        '''
        This function retrieves the pressures at the demand nodes and returns a dataframe with the pressures and the time

        specific_nodes (list-like): list of specific nodes to retrieve the pressures from
        '''
        demand_nodes = self.demand_nodes
        if specific_nodes:
            demand_nodes = specific_nodes
        pressures = pd.DataFrame()
        swtch = True
        # Open the SWMM Output
        with pyswmm.Output(self.output_file) as out:
            for node in out.nodes:
                if swtch:
                    pressures['Time'] = pd.Series(out.node_series(node, 'INVERT_DEPTH').keys())
                    swtch = False
                
                if node in demand_nodes:
                    pressures[node] = pd.Series(out.node_series(node, 'INVERT_DEPTH').values())
        self.pressures  = pressures
        return pressures
    
    def get_tank_vols_heights(self, specific_nodes = None):
        '''
        This function retrieves the volumes stored at each user tank and returns a dataframe with the tank volumes and the time

        specific_nodes (list-like): list of specific nodes to retrieve the volumes for
        '''
        storage_ids = self.storage_ids
        if specific_nodes:
            storage_ids = ['StorageforNode' +str(node) for node in specific_nodes]
        tank_vols = pd.DataFrame()
        tank_heights = pd.DataFrame()
        swtch = True
        # Open the SWMM Output
        with pyswmm.Output(self.output_file) as out:
            for node in out.nodes:
                if swtch:
                    tank_vols['Time'] = pd.Series(out.node_series(node, 'PONDED_VOLUME').keys())
                    tank_heights['Time'] = pd.Series(out.node_series(node, 'INVERT_DEPTH').keys())
                    swtch = False
                
                if node in storage_ids:
                    tank_vols[node.split('Node')[1]] = pd.Series(out.node_series(node, 'PONDED_VOLUME').values())
                    tank_heights[node.split('Node')[1]] = pd.Series(out.node_series(node, 'INVERT_DEPTH').values())

        self.tank_vols_result = tank_vols
        self.tank_heights_result = tank_heights

        return self.tank_vols_result, self.tank_heights_result

    def get_withdrawals(self, specific_nodes = None):
        '''
        This function retrieves the withdrawal rates at the demand nodes and returns a dataframe with the withdrawal rates and the time

        specific_nodes (list-like): list of specific nodes to retrieve the withdrawal rates from
        '''
        demand_nodes = self.demand_nodes
        if specific_nodes:
            demand_nodes = specific_nodes 
        with_ids = ["Outlet"+str(node) for node in demand_nodes]
        withdrawals = pd.DataFrame()
        swtch = True
        # Open the SWMM Output
        with pyswmm.Output(self.output_file) as out:
            for link in out.links:
                if swtch:
                    withdrawals['Time'] = pd.Series(out.link_series(link, 'FLOW_RATE').keys())
                    swtch = False
                
                if link in with_ids:
                    withdrawals[link.split('Outlet')[1]] = pd.Series(out.link_series(link, 'FLOW_RATE').values())

        self.withdrawals = withdrawals

        return withdrawals

    def get_consumptions(self, specific_nodes = None):
        '''
        This function retrieves the consumption rates at the demand nodes and returns a dataframe with the consumption rates and the time

        specific_nodes (list-like): list of specific nodes to retrieve the consumption rates from
        '''
        demand_nodes = self.demand_nodes
        if specific_nodes:
            demand_nodes = specific_nodes
        cons_ids = ["ConsumptionOutlet"+str(node) for node in demand_nodes]
        consumptions = pd.DataFrame()
        swtch = True
        # Open the SWMM Output
        with pyswmm.Output(self.output_file) as out:
            for link in out.links:
                if swtch:
                    consumptions['Time'] = pd.Series(out.link_series(link, 'FLOW_RATE').keys())
                    swtch = False
                
                if link in cons_ids:
                    consumptions[link.split("Outlet")[1]] = pd.Series(out.link_series(link, 'FLOW_RATE').values())

        self.consumptions = consumptions
        return consumptions
    
    def get_leaks(self, specific_nodes = None):
        demand_nodes = self.demand_nodes
        if specific_nodes:
            demand_nodes = specific_nodes
        leak_ids = ["LeakforNode"+str(node) for node in demand_nodes]
        leakage = pd.DataFrame()
        swtch = True
        # Open the SWMM Output
        with pyswmm.Output(self.output_file) as out:
            for link in out.links:
                if swtch:
                    leakage['Time'] = pd.Series(out.link_series(link, 'FLOW_RATE').keys())
                    swtch = False
                
                if link in leak_ids:
                    leakage[link.split("Node")[1]] = pd.Series(out.link_series(link, 'FLOW_RATE').values())

        self.leakage = leakage
        return leakage

    def __get_source_volumes__(self):
        '''
        This function retrieves the volumes stored at the source reservoirs and returns a dataframe with the volumes and the time
        '''
        reservoir_ids = self.reservoir_ids_new
        source_volumes = pd.DataFrame()
        swtch = True
        # Open the SWMM Output
        with pyswmm.Output(self.output_file) as out:
            for node in out.nodes:
                if swtch:
                    source_volumes['Time'] = pd.Series((pyswmm.NodeSeries(out)[node].ponded_volume).keys())
                    swtch = False
                
                if node in reservoir_ids:
                    source_volumes[node] = pd.Series((pyswmm.NodeSeries(out)[node].ponded_volume).values())

        self.source_volumes = source_volumes
        return source_volumes
    
    def get_pipe_network_vols(self):
        '''
        Calculates the volume stored in network pipes as a timeseries
        '''

        pipe_ids = self.conduits.index
        pipe_volumes = pd.DataFrame()
        swtch = True
        # Open the SWMM Output
        with pyswmm.Output(self.output_file) as out:
            for link in out.links:
                if swtch:
                    pipe_volumes['Time'] = pd.Series((pyswmm.LinkSeries(out)[link].flow_volume).keys())
                    swtch = False
                if re.search("^P",link):
                    pipe_volumes[link] = pd.Series((pyswmm.LinkSeries(out)[link].flow_volume).values())

        pipe_volumes_total = pipe_volumes.iloc[:,1:].sum(axis = 1)
        self.pipe_volumes = pipe_volumes_total
        return pipe_volumes_total


    def track_water_balances(self, plot = True, all_source_vol = False):
        '''
        Computes the amounts of water in the different components of the network: 
        Source reservoirs, network pipes, user tanks, consumed by users, and leaked water.
        

        returns: water balance dataframe
        '''
        attribs = { 'source_volumes': self.__get_source_volumes__,
                    'tank_vols_result': self.get_tank_vols_heights,
                    'pipe_volumes_total': self.get_pipe_network_vols,
                    'consumptions': self.get_consumptions,
                    'leakage': self.get_leaks}
        # Get the volumes in the reservoirs
        for attr in attribs.keys():
            if not hasattr(self, attr):
                attribs[attr]()

        # self.__get_source_volumes__()
        # self.get_tank_vols_heights()
        # self.get_pipe_network_vols()
        # self.get_consumptions()
        # self.get_leaks()

        #reporting step
        reporting_step = (self.consumptions.at[1,'Time'] - self.consumptions.at[0,'Time']).total_seconds()

        # Get the volumes in the reservoirs
        self.__get_source_volumes__()

        # get cumulative consumption totals
        consumptions = self.consumptions
        consumptions.set_index('Time', inplace = True)
        total_consumptions = consumptions.sum(axis = 1)
        cum_consumptions = sci.cumulative_trapezoid(total_consumptions, dx = reporting_step) / 1000 # in m3
        cum_consumptions = np.append(0,cum_consumptions)

        # get cumulative leakage totals
        leakage = self.leakage
        leakage.set_index('Time', inplace=True)
        total_leakage = leakage.sum(axis = 1)
        cum_leakage = sci.cumulative_trapezoid(total_leakage, dx = reporting_step) / 1000 # in m3
        cum_leakage = np.append(0,cum_leakage)

        # get total pipe volumes
        pipe_vols = self.pipe_volumes
        pipe_vols.index = consumptions.index

        # get total tank volumes
        tank_vols = self.tank_vols_result
        tank_vols.set_index('Time', inplace = True)
        tank_vols_total = tank_vols.sum(axis = 1)

        # source volumes 
        source_volumes = self.source_volumes.iloc[:,1:].sum(axis = 1)
        source_volumes.index = consumptions.index
        if not all_source_vol:
            source_volumes = source_volumes - source_volumes.min()

        self.cum_consumptions = cum_consumptions
        self.cum_leakage = cum_leakage
        self.tank_vols_total = tank_vols_total
        self.source_volumes = source_volumes
        # get total system volumes
        total_sys_volume = source_volumes + cum_consumptions + cum_leakage + tank_vols_total + pipe_vols
        base_volume = total_sys_volume.iloc[0]
        # get time column 
        time = consumptions.index
        time = [(x-time[0]).total_seconds()/3600 for x in time]
        water_balance=pd.DataFrame(zip(time,self.source_volumes,pipe_vols,tank_vols_total,cum_consumptions,cum_leakage,total_sys_volume),
                           columns=["Time","Reservoirs","In Pipes","Stored","Consumed","Leaked","Total"])
        self.water_balance = water_balance
        
        if plot:
            fig, ax = plt.subplots(1,2)
            fig.set_figheight(9)
            fig.set_figwidth(12)

            yaxis = np.vstack([water_balance["Leaked"],water_balance["Consumed"],water_balance["Stored"],water_balance["In Pipes"],water_balance["Reservoirs"]])
            ax[0].stackplot(water_balance["Time"],yaxis,labels=["Leaked","Consumed","Stored","In Pipes","In Reservoirs"],colors=["#e66101","#fdb863","#fee0b6","#b2abd2","#5e3c99"])
            ax[0].set_xlim(0,max(water_balance["Time"]))
            ax[0].set_xlabel("Time (hr)")
            ax[0].set_ylabel("Total Volume")
            ax[0].set_xticks(np.arange(0,max(water_balance["Time"])*1.01,round((max(water_balance["Time"]))/4,0)))
            ax[0].set_xticks(np.arange(0,max(water_balance["Time"])*1.01,round((max(water_balance["Time"]))/12,2)),minor=True)
            ax[0].set_yticks(np.arange(0,int(max(water_balance["Total"])+1),int(max(water_balance["Total"])/10)))
            ax[0].set_ylim(0,round(max(water_balance["Total"])+1))
            ax[0].legend(loc="lower right")
            ax[0].spines[['right', 'top']].set_visible(False)

            line1,=ax[1].plot(water_balance["Time"],water_balance["Reservoirs"],label="In Reservoirs",color="#5e3c99")
            line2,=ax[1].plot(water_balance["Time"],water_balance["In Pipes"],label="In Pipes",color="#b2abd2")
            line3,=ax[1].plot(water_balance["Time"],water_balance["Stored"],label="Stored",color="#fee0b6")
            line4,=ax[1].plot(water_balance["Time"],water_balance["Consumed"],label="Consumed",color="#fdb863")
            line6,=ax[1].plot(water_balance["Time"],water_balance["Leaked"],label="Leaked",color="#e66101")
            line5,=ax[1].plot(water_balance["Time"],water_balance["Total"],label="Total",color='black')


            ax[1].set_xlim(0,max(water_balance["Time"]))
            ax[1].set_xlabel("Time (hr)")
            ax[1].set_ylabel(r"Volume (1000 m3)")
            ax[1].set_xticks(np.arange(0,max(water_balance["Time"])*1.01,round((max(water_balance["Time"]))/4,0)))
            ax[1].set_xticks(np.arange(0,max(water_balance["Time"])*1.01,round((max(water_balance["Time"]))/12,2)),minor=True)
            ax[1].set_yticks(np.arange(0,int(max(water_balance["Total"])+1),int(max(water_balance["Total"])/10)))
            ax[1].set_ylim(0,int(max(water_balance["Total"])+1))
            # ax[0,1].legend(loc="upper right")
            ax[1].spines[['right', 'top']].set_visible(False)


            plt.show()

        return self.water_balance
    
    def mass_balance(self, plot=True):
        '''
        Computes and visualizes the mass balance of the system

        plot: bool: whether to plot the mass balance or not

        returns: mass balance and instantaneous mass balance dataframes
        '''
        # Get the volumes in the reservoirs
        if not hasattr(self, 'water_balance'):
            water_balance = self.track_water_balances(plot = False)
        else: water_balance = self.water_balance


        base_volume = water_balance['Total'].iloc[0]
        cont_error = (water_balance['Total'] - base_volume) / base_volume * 100
        inst_mass_balance = water_balance['Total'].diff().fillna(0)

        fig, ax = plt.subplots(1,2)
        fig.set_figheight(5)
        fig.set_figwidth(12)

        ax[0].plot(water_balance['Time'], cont_error, label = 'Continuity Error (%)', color = 'black')
        ax[0].set_xlabel('Time (hr)')
        ax[0].set_ylabel('Continuity Error (%)')
        ax[0].set_xlim(0, max(water_balance['Time']))
        ax[0].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
        ax[0].set_ylim(-1, 1)

        ax[1].plot(water_balance['Time'], inst_mass_balance, label = 'Instantaneous Mass Balance (m3)', color = 'black')
        ax[1].set_xlabel('Time (hr)')
        ax[1].set_ylabel('Instantaneous Mass Balance (m3)')
        ax[1].set_xlim(0, max(water_balance['Time']))
        # ax[1].set_ylim(-1, 1)

        return cont_error, inst_mass_balance

    def effective_supply_duration(self, pressure_threshold = None, plot = True):
        '''
        Computes the effective supply duration for each node and visualizes it
        If a pressure threshold is defined, the effective supply duration is computed based on the time each node spends with pressure > pressure threshold
        If no pressure threshold is defined, the effective supply duration is the time with non-zero withdrawals
        Parameters:
        pressure_threshold: float: the pressure threshold below which the supply is considered ineffective
        '''
 
        if not pressure_threshold:
            withdrawals = self.get_withdrawals()
            withdrawals.set_index('Time', inplace = True)
            supply_on = withdrawals.apply(lambda x: x > 0)
            reporting_step = (withdrawals.index[1] - withdrawals.index[0]).total_seconds()
            eff_supply_duration = supply_on.sum() * reporting_step / 3600
        else:
            pressures = self.get_pressures()
            pressures.set_index('Time', inplace = True)
            supply_on = pressures.apply(lambda x: x > pressure_threshold)
            reporting_step = (withdrawals.index[1] - withdrawals.index[0]).total_seconds()
            eff_supply_duration = supply_on.sum() * reporting_step / 3600
        
        fig, ax = plt.subplots()
        fig.set_figheight(5)
        fig.set_figwidth(8)
        hist, bins, patches = ax.hist(eff_supply_duration, bins=20, color='black')
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        col = bin_centers - min(bin_centers)
        col /= max(col)
        for c, p in zip(col, patches):
            plt.setp(p, 'facecolor', plt.cm.viridis(c))
        ax.set_xlabel('Effective Supply Duration (hr)')
        ax.set_ylabel('Frequency')

        return eff_supply_duration


    

        



