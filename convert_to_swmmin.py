#%%
import wntr
import numpy as np      
import pandas as pd
import re
import math
from pathlib import Path
import logging

logger=logging.getLogger()
logger.setLevel(logging.DEBUG)
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s')

def load_epanet(
        inp_file: Path, 
        supply_duration_inp: float, 
        concentric: bool,
        len_to_diameter_ratio: float,
        adaptive: bool,
        maximum_xdelta: float,
        tank_height: float,
        minimum_pressure: float,
        pressure_diff: float,
        pdw_variable: str,
        h_min: float,
        h_max: float,
        step: float,
        m: float,
        n: float,
        pattern: list,
        leak_fraction: float,
        n_days: int,
        storage_initial_fullness_factor: float,
        reservoir_height: float,
        tank_areas: Path | None = None
        ):
  
  demand_nodes=[]       #  Modes that have non-zero demands
  base_demands=[]       #  Demand rates desired by each node for desired volume calculations
  elevations=[]         #  Elevations of demand nodes
  coords=dict()         #  Coordinates corresponding to each node as a tuple with the id as key
  all_nodes=[]          #  List of node ids of all nodes
  all_elevations=[]     #  Elevations of all nodes

  network = wntr.network.WaterNetworkModel(inp_file)
  for node in network.junctions():
    all_nodes.append(node[1].name)
    all_elevations.append(node[1].elevation)
    coords[node[1].name]=node[1].coordinates
    if node[1].base_demand != 0:
        demand_nodes.append(node[1].name)
        base_demands.append(node[1].base_demand)
        elevations.append(node[1].elevation)

  conduit_ids= []       # IDs of the original pipes in the EPANET file
  conduit_from= []      # origin node for each pipe
  conduit_to= []        # destination node for each pipe
  conduit_lengths= []   # pipe lengths
  conduit_diameters= [] # pipe diameters

  for link in network.pipes():
    # Extract and store each of the aforementioned properties
    conduit_ids.append(link[1].name)
    conduit_from.append(link[1].start_node_name)
    conduit_to.append(link[1].end_node_name)
    conduit_lengths.append(link[1].length)
    conduit_diameters.append(link[1].diameter)

  pump_ids= []       # IDs of the original pumps in the EPANET file
  pump_from= []      # origin node for each pump
  pump_to= []        # destination node for each pump
  pump_curves = []   # name of the pump's curve


  for pump in network.pumps():
      pump_ids.append(pump[1].name)
      pump_from.append(pump[1].start_node_name)
      pump_to.append(pump[1].end_node_name)
      pump_curves.append(pump[1]._pump_curve_name)


  reservoir_ids=[]      # source reservoirs' IDs
  reservoir_heads={}    # total head of each reservoir indexed by ID
  reservoir_coords={}   # coordinates as tuple (x,y) indexed by ID

  for reservoir in network.reservoirs():
      reservoir_ids.append(reservoir[1].name)
      reservoir_heads[reservoir_ids[-1]]=reservoir[1].base_head
      reservoir_coords[reservoir_ids[-1]]=reservoir[1].coordinates
  reservoir_elevations={reservoir:reservoir_heads[reservoir]-reservoir_height for reservoir in reservoir_heads}

  pcurves = {}
  for curve in network.curves():
      pcurves[curve[0]]= curve[1].points

  supply_duration=int(network.options.time.duration/60)
  if supply_duration_inp:
      supply_duration=int(supply_duration_inp*60)
  supply_hh=str(supply_duration//60)     
  supply_mm=str(supply_duration%60)      

  # Format to HH:MM if necessary
  if len(supply_mm)<2:
      supply_mm='0'+supply_mm
  if len(supply_hh)<2:
      supply_hh='0'+supply_hh

  desired_demands = [demand * 24 / supply_duration_inp for demand in base_demands] # Desired demand rates in CMS for each node
  ## TODO: Input Q_desired.csv file for desired demands   

  junctions=pd.DataFrame(zip(all_nodes,all_elevations,coords.values()),columns=["ID","Elevation","Coordinates"])
  junctions.set_index("ID",inplace=True)
  conduits=pd.DataFrame(zip(conduit_ids,conduit_from,conduit_to,conduit_lengths,conduit_diameters),columns=["ID","from node","to node","Length","diameter"])
  conduits.set_index("ID",inplace=True)
  connectivity=pd.DataFrame(index=junctions.index, columns=["US","DS"])

  if concentric:
      for junction in junctions.index:
          upst_connected=[] #Pipes whose upstream end connects to this node
          downst_connected=[] #Pipes whose downstream end connects to this node
          for conduit in conduits.index:
              # if conduit starts at this node
              if conduits.at[conduit,"from node"]==junction:
                  upst_connected.append(conduit)
              # if conduit ends at this node
              elif conduits.at[conduit,"to node"]==junction:
                  downst_connected.append(conduit)
          # create a list of all connected conduits
          connected_conduits=upst_connected+downst_connected
          # Get the diameter of all connected conduits
          connected_diameters=[conduits.at[conduit,"diameter"] for conduit in connected_conduits]
          max_diameter=max(connected_diameters)

          for conduit in upst_connected:
              if conduits.at[conduit,"diameter"]<max_diameter:
                  conduits.at[conduit,"InOffset"]=(max_diameter-conduits.at[conduit,"diameter"])/2
              else: conduits.at[conduit,"InOffset"]=0

          for conduit in downst_connected:
              if conduits.at[conduit,"diameter"]<max_diameter:
                  conduits.at[conduit,"OutOffset"]=(max_diameter-conduits.at[conduit,"diameter"])/2
              else: conduits.at[conduit,"OutOffset"]=0
      
          connectivity.at[junction,"US"]=upst_connected
          connectivity.at[junction,"DS"]=downst_connected
          connectivity.at[junction,"Max D"]=max_diameter

  reservoir_pipes=[]

  if maximum_xdelta==0:
      maximum_xdelta=conduits["diameter"].max() * len_to_diameter_ratio

  for conduit in conduits.index:
      length=conduits["Length"][conduit] 
      diameter = conduits["diameter"][conduit]
      if adaptive:
          n_parts = np.max(round(length / diameter / len_to_diameter_ratio) - 1, 0)
      else: n_parts = math.ceil(length / maximum_xdelta)

      # If the conduit is bigger than the maximum allowable length (delta x), we will break it down into smaller pipes
      if n_parts>0:
          part_length=length/n_parts
          start_node=conduits["from node"][conduit]
          in_offset=conduits["InOffset"][conduit]
          end_node=conduits["to node"][conduit]
          out_offset=conduits["OutOffset"][conduit]

          if start_node in reservoir_ids:
              # MAke the start elevation the same as the end but add 1 (since reservoirs don't have ground elevation in EPANET)
              start_elevation=junctions.at[end_node,"Elevation"]+1
              reservoir_elevations[start_node]=start_elevation+1
          # Otherwise make the start elevation equal to the elevation of the start node + the offset of the pipe
          else: start_elevation=junctions.at[start_node,"Elevation"]+in_offset
          
          if end_node in reservoir_ids:
              # Make the end elevation the same as the start but subtract 1 (since reservoirs don't have ground elevation in EPANET)
              end_elevation=start_elevation-1
          # Make the end elevation equal to the elevation of the end node + the offset of the pipe
          else: end_elevation=junctions.at[end_node,"Elevation"]+out_offset
          # Calculate the uniform drop (or rise) in elevation for all the intermediate nodes about to be created when this pipe is broken into several smaller ones
          unit_elev_diff=(end_elevation-start_elevation)/n_parts

          # if the starting node is a reservoir
          if start_node in reservoir_ids:
              start_x=reservoir_coords[start_node][0]
              start_y=reservoir_coords[start_node][1]
          else:
              start_x=junctions.at[start_node,"Coordinates"][0]
              start_y=junctions.at[start_node,"Coordinates"][1]
          
          # If the end node is a reservoir
          if end_node in reservoir_ids:
              end_x=reservoir_coords[end_node][0]
              end_y=reservoir_coords[end_node][1]
          else:
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
              conduits.at[new_id,"diameter"]=conduits["diameter"][conduit]
              # Set the start node as OriginStartNode-NewNodeNumber-OriginEndNode  as in the first intermediate nodes between node 13 and 14 will be named 13-1-14
              conduits.at[new_id,"from node"]=start_node+"-"+str(part-1)+"-"+end_node
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

              if conduits.at[new_id,"from node"] in reservoir_ids or conduits.at[new_id,"to node"] in reservoir_ids:
                  reservoir_pipes.append(new_id)
              
          # After writing the new smaller pipes, delete the original pipe (since it is now redundant)
          conduits.drop(conduit,inplace=True)

  conduits[["InOffset","OutOffset"]]=conduits[["InOffset","OutOffset"]].fillna(0)

  # Formatting Junctions
  MaxDepth=[0]*len(junctions)
  InitDepth=MaxDepth
  SurDepth=[100] * len(junctions)  # High value to prevent surcharging
  Aponded=InitDepth
  junctions_section=pd.DataFrame(list(zip(junctions.index,junctions["Elevation"],MaxDepth,InitDepth,SurDepth,Aponded)))
  junctions_section=junctions_section.to_string(header=False,index=False,col_space=10).splitlines()
  junctions_section=[line+'\n' for line in junctions_section]

  # Formatting Outfalls
  outfall_ids=["Outfall"+ str(node) for node in demand_nodes]
  outfall_elevations=elevations
  outfall_types=['FREE' for i in demand_nodes]
  stage_data=["    " for i in demand_nodes]
  outfall_gated=["NO" for i in demand_nodes]
  outfall_section=pd.DataFrame(zip(outfall_ids,outfall_elevations,outfall_types,stage_data,outfall_gated))
  leak_outfall_id=["L_Outfall"+str(node) for node in demand_nodes]
  leak_outfalls=pd.DataFrame(zip(leak_outfall_id,outfall_elevations,outfall_types,stage_data,outfall_gated))
  outfall_section=pd.concat([outfall_section,leak_outfalls])
  outfall_section=outfall_section.to_string(header=False,index=False,col_space=10).splitlines()
  outfall_section=[line+'\n' for line in outfall_section]

  # Writing Storage Nodes
  storage_ids=["StorageforNode"+id for id in demand_nodes]
  # Calculate the area of each storage based on their demand
  demands_total=[demand*60* supply_duration/tank_height for demand in base_demands]
  demand_volumes = pd.DataFrame({'ID':demand_nodes,'Volume':demands_total})
  demand_volumes.to_csv(inp_file.parent/(inp_file.name+'_DemandVolumes.csv'))
  storage_areas=demands_total
  # OVERRIDE OPTION: Input Tank Areas Manually
  if tank_areas is not None:
      storage_areas = pd.read_csv(tank_areas,header=None)
      storage_areas.columns = ['Area']
      storage_areas = storage_areas['Area'].to_list()
      storage_areas = [area*60* supply_duration/tank_height for area in storage_areas]
  # Set storage invert elevations = elevation of DN + Hmin
  storage_elevations=[elevation+minimum_pressure for elevation in elevations]
  if type(tank_height) == str:
      storage_MaxDepth = pd.read_csv(tank_height,header=None)
      storage_MaxDepth.columns = ['Height']
      storage_MaxDepth = storage_MaxDepth['Height'].to_list()
  else: storage_MaxDepth=[tank_height]*len(storage_ids)
  storage_InitDepth = [storage_initial_fullness_factor*depth for depth in storage_MaxDepth] # Depth measured from the bottom of the tank up
  logger.debug(f"Storage Elevations: {storage_elevations}")
  logger.debug(f"Storage Init Depths: {storage_InitDepth}")
  logger.debug(f"Storage Max Depths: {storage_MaxDepth}")
  storage_shape=["FUNCTIONAL"]*len(storage_ids)
  zeroes=['0']*len(storage_ids)    # for the other curve parameters which are not required for circular tanks
  storage_SurDepth=[0]*len(storage_ids)
  storage_fevap=[0]*len(storage_ids)
  storage_units=pd.DataFrame(zip(storage_ids,storage_elevations,storage_MaxDepth,storage_InitDepth,storage_shape,zeroes,zeroes,storage_areas,storage_SurDepth,storage_fevap))

  # Formatting Reservoirs
  reservoir_ids_new=["Reservoir-"+str(id) for id in reservoir_ids]
  # update references to reservoirs in conduits
  for i in range(0,len(reservoir_ids)):
      conduits.replace(to_replace=reservoir_ids[i],value=reservoir_ids_new[i])
  reservoir_elevations=reservoir_elevations.values()
  MaxDepth=[max(100,max(reservoir_heads.values())+10)]*len(reservoir_ids)
  InitDepth=[head-elevation for head, elevation in zip(reservoir_heads.values(),reservoir_elevations)]
  reservoir_shape=["TABULAR"]*len(reservoir_ids)
  reservoir_curves=["Source"+id for id in reservoir_ids]
  blanks=['    ']*len(storage_ids)  
  reservoir_SurDepth=[0]*len(reservoir_ids)
  reservoir_psi=[0]*len(reservoir_ids)
  storage_section=pd.DataFrame(zip(reservoir_ids_new,reservoir_elevations,MaxDepth,InitDepth,reservoir_shape,reservoir_curves,blanks,blanks,reservoir_SurDepth,reservoir_psi))
  storage_section= pd.concat([storage_section,storage_units])
  storage_section=storage_section.to_string(header=False,index=False,col_space=10).splitlines()
  storage_section=[line+'\n' for line in storage_section]

  # Formatting Conduits
  conduit_ids=["P-"+str(conduit) if str(conduit)[0]!="P" else conduit for conduit in conduits.index]
  roughness=[0.011]*len(conduits)
  conduit_zeros=[0]*len(conduits)
  conduits["from node"].replace(reservoir_ids,reservoir_ids_new,inplace=True)
  conduits["to node"].replace(reservoir_ids,reservoir_ids_new,inplace=True)
  conduits_section=pd.DataFrame(zip(conduit_ids,conduits["from node"],conduits["to node"],conduits["Length"],roughness,conduits["InOffset"],conduits["OutOffset"],conduit_zeros,conduit_zeros))
  conduits_section=conduits_section.to_string(header=False,index=False,col_space=10).splitlines()
  conduits_section=[line+'\n' for line in conduits_section]

  # Formatting Outlets
  outlet_ids = ["Outlet"+id for id in demand_nodes]
  outlet_from = demand_nodes[:]
  outlet_to = storage_ids [:]
  outlet_offset=[connectivity.at[node,"Max D"]/2 for node in demand_nodes]
  # Change to FUNCTIONAL/PRESSURE for head difference PDW, rather than upstream depth
  outlet_type=["FUNCTIONAL/"+pdw_variable]*len(outlet_ids)
  # Desired demands are taken as base demands of EPANET input file by default. To change this, either update the EPANET input file or input a new set of desired demands to desired_demands
  outlet_coeff=[demand*1000/np.sqrt(pressure_diff) for demand in desired_demands]  
  # Change to change PDW exponent
  outlet_expon=["0.5"]*len(outlet_ids)
  outlet_gated=["YES"]*len(outlet_ids)
  outdemand_ids = ["DemandOutlet"+id for id in demand_nodes]
  outdemand_from = storage_ids [:]
  outdemand_to = outfall_ids [:]
  outdemand_offset=[0]*len(outlet_ids)
  outdemand_type=["TABULAR/DEPTH"]*len(outlet_ids)
  outdemand_coeff=["Demand"+ id for id in demand_nodes]
  outdemand_expon=["     "]*len(outlet_ids)
  outdemand_gated=["YES"]*len(outlet_ids)
  outlets=pd.DataFrame(list(zip(outlet_ids,outlet_from,outlet_to,outlet_offset,outlet_type,outlet_coeff,outlet_expon,outlet_gated)))
  outletdemand=pd.DataFrame(list(zip(outdemand_ids,outdemand_from,outdemand_to,outdemand_offset,outdemand_type,outdemand_coeff,outdemand_expon,outdemand_gated)))

  # Formatting Leakage Outlets
  leak_ids=["LeakforNode"+str(node) for node in demand_nodes]
  leak_from=demand_nodes[:]
  leak_to=leak_outfall_id[:]
  leak_offset=[connectivity.at[node,"Max D"]/2 for node in demand_nodes]
  leak_type=["FUNCTIONAL/DEPTH"]*len(leak_ids)
  leak_coeff=[demand*1000/np.sqrt(pressure_diff)*leak_fraction for demand in desired_demands]
  leak_expon=["0.5"]*len(leak_ids)
  leak_gated=["YES"]*len(leak_ids)
  leak_outlets=pd.DataFrame(zip(leak_ids,leak_from,leak_to,leak_offset,leak_type,leak_coeff,leak_expon,leak_gated))
  outlet_section=pd.concat([outlets,outletdemand,leak_outlets])
  outlet_section=outlet_section.to_string(header=False,index=False,col_space=10).splitlines()
  outlet_section=[line+'\n' for line in outlet_section]

  # Formatting Pumps
  pump_curves = ["Pump_Curve_"+curve for curve in pump_curves]
  pump_status = ["ON" for pump in pump_ids]
  pump_startup = [0 for pump in pump_ids]
  pump_shutoff = [0 for pump in pump_ids]
  pump_section=pd.DataFrame(zip(pump_ids,pump_from,pump_to,pump_curves,pump_status,pump_startup,pump_shutoff))
  pump_section = pump_section.to_string(header=False,index=False,col_space=10).splitlines()
  pump_section = [line+'\n' for line in pump_section]

  # Formatting Xsections
  shape=["FORCE_MAIN"]*len(conduits.index)
  hwcoeffs=[130]*len(shape)
  geom3=[0]*len(shape)
  geom4=geom3
  nbarrels=[1]*len(shape)
  xsections_section=pd.DataFrame(zip(conduit_ids,shape,conduits["diameter"],hwcoeffs,geom3,geom4,nbarrels))
  xsections_section=xsections_section.to_string(header=False,index=False, col_space=10).splitlines()
  xsections_section=[line+'\n' for line in xsections_section]

  # Formatting Curves
  curves_name=[]
  curves_type=[]
  curves_x=[]
  curves_y=[]
  constant_demands=[demand*float(supply_hh)/24 for demand in desired_demands]
  for i,j in zip(demand_nodes,constant_demands):
      curves_name+=["Demand"+str(i),"Demand"+str(i),';']
      curves_type+=["Rating"," ","  "]
      curves_x+=[0,0.01," "]
      curves_y+=[0,j*1000," "]
  # Reservoir Storage Curves
  # reservoir_volume = sum(storage_areas) * (n_days - 1 + 2) 
  reservoir_volume = 5
  logger.debug(f"Reservoir Volume: {reservoir_volume} m^3")
  for curve,head,elevation in zip(reservoir_curves,reservoir_heads.values(),reservoir_elevations):
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

  for curve in pcurves:
      if len(pcurves[curve]) == 1:
          design_flow = pcurves[curve][0][0] * 1000
          design_head = pcurves[curve][0][1]

          shutoff_head = 1.33 * design_head 
          max_flow = 2 * design_flow
      # Fit an H(Q) = H_max - B * Q^2 Function
          B = (shutoff_head - design_head) / design_flow ** 2
          flows = []
          heads = []
          for flow in np.arange(0, max_flow + max_flow/20, max_flow/20):
              flows.append(flow)
              heads.append(shutoff_head - B * flow **2)
      for flow, head in zip(flows, heads):
          curves_name.append("Pump_Curve_"+ curve)
          if flow  == flows[0]:
              curves_type.append("Pump3")
          else: curves_type.append("  ")
          curves_x.append(flow)
          curves_y.append(head)
  curves=pd.DataFrame(list(zip(curves_name,curves_type,curves_x,curves_y)))
  curves_section=curves.to_string(header=False,index=False,col_space=10)

  # Formatting Controls
  control_curves=""
  control_rules=""

  times=np.arange(0,int(supply_hh)*50,1)
  timesrs_pat=""
  for time in np.arange(0,len(times)):
      if time >=len(pattern):
          days=math.floor(time/len(pattern))
          time_24=time-days*len(pattern)
      else: time_24=time
      timesrs_pat+="Pattern\t"+str(time)+"\t"+str(pattern[time_24])+"\n"


  for outlet in outletdemand.iloc[:,0]:
      out_name="Outlet"+outlet[12:]
      storage_name=outletdemand.loc[outletdemand[0]==outlet,1].iloc[0]
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

  # Add rule to stop supply by iterating over reservoir adjacent pipes
  control_rules+="Rule STOPSUPPLY\n"
  control_rules+="IF SIMULATION CLOCKTIME > "+str(supply_hh)+":"+str(supply_mm)+"\n"
  reservoirs_list = network.reservoir_name_list
  for reservoir in reservoirs_list:
      pipes = network.get_links_for_node(reservoir)
      for pipe in pipes:
          if network.get_link(pipe).length > maximum_xdelta:
              suffix = '-1'
          else: suffix = ''
          if reservoir == reservoirs_list[0]:
              print(pipes[0])
              control_rules+="THEN CONDUIT "+pipe+suffix+" STATUS = CLOSED\n"
          else: control_rules+="AND CONDUIT "+pipe+suffix+" STATUS = CLOSED\n"
  control_rules+="\n"

  # Add rule to start supply by iterating over reservoir adjacent pipes
  control_rules+="Rule STARTSUPPLY\n"
  control_rules+="IF SIMULATION CLOCKTIME >= 0:00\n"
  control_rules += "AND SIMULATION CLOCKTIME < "+str(supply_hh)+":"+str(supply_mm)+"\n"
  for reservoir in reservoirs_list:
      pipes = network.get_links_for_node(reservoir)
      for pipe in pipes:
          if network.get_link(pipe).length > maximum_xdelta:
              suffix = '-1'
          else: suffix = ''
          if reservoir == reservoirs_list[0]:
              control_rules+="THEN CONDUIT "+pipe+suffix+" STATUS = OPEN\n"
          else: control_rules+="AND CONDUIT "+pipe+suffix+" STATUS = OPEN\n"
  control_rules+="\n"
  control_rules+="Rule Patterns\n"
  control_rules+="IF SIMULATION TIME > 0\n"
  flag=0
  for outlet in outletdemand.iloc[:,0]:
      if flag==0:
          control_rules+="THEN OUTLET "+outlet+" SETTING = TIMESERIES Pattern\n"
          flag=1
      else: control_rules+="AND OUTLET "+outlet+" SETTING = TIMESERIES Pattern\n"
  curves_section+="\n"+control_curves

  # Formatting Coordinates
  coords_demand= { node: coords[node] for node in demand_nodes}
  coords_ids=list(junctions.index)+reservoir_ids_new+storage_ids+outfall_ids+leak_outfall_id
  offset=2
  coords_x1=[coord[0] for coord in junctions["Coordinates"]]
  coords_x2=[coord[0] for coord in reservoir_coords.values()]
  coords_x3=[coord[0] +2*offset for coord in coords_demand.values()]
  coords_x4=[coord[0] + 3 * offset for coord in coords_demand.values()]
  coords_x5=[coord[0] + offset for coord in coords_demand.values()]
  coords_x=coords_x1+coords_x2+coords_x3+coords_x4+coords_x5
  coords_y1=[coord[1] for coord in junctions["Coordinates"]]
  coords_y2=[coord[1] for coord in reservoir_coords.values()]
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

  # Return necessary stuff for write_swmmin
  return (
      n_days,
      maximum_xdelta,
      dimensions_line,
      coordinate_section,
      timesrs_pat,
      curves_section,
      control_rules,
      xsections_section,
      outlet_section,
      pump_section,
      conduits_section,
      storage_section,
      outfall_section,
      junctions_section
  )



def write_swmmin(
        template: Path, 
        output: Path,
        n_days: int,
        maximum_xdelta: float,
        timestep: float | None,
        dimensions_line: str,
        coordinate_section: list,
        timesrs_pat: list,
        curves_section: list,
        control_rules: str,
        xsections_section: list,
        outlet_section: list,
        pump_section: list,
        conduits_section: list,
        storage_section: list,
        outfall_section: list,
        junctions_section: list
        ):
  file = open(template, 'r')
  lines=[]
  linecount=0 

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
    if re.search('\\[JUNCTIONS\\]',line):
        junctions_marker=linecount+3
    # Record the position of the phrase [OUTFALLS] and add 3 to skip the header lines
    if re.search('\\[OUTFALLS\\]',line):
        outfalls_marker=linecount+3
    # Record the position of the phrase [STORAGE] and add 3 to skip the header lines
    if re.search('\\[STORAGE\\]',line):
        storage_marker=linecount+3
    # Record the position of the phrase [CONDUITS] and add 3 to skip the header lines
    if re.search('\\[CONDUITS\\]',line):
        conduits_marker=linecount+3
     # Record the position of the phrase [OUTLETS] and add 3 to skip the header lines
    if re.search('\\[PUMPS\\]',line):
        pumps_marker=linecount+3
    if re.search('\\[OUTLETS\\]',line):
        outlets_marker=linecount+3
     # Record the position of the phrase [XSECTIONS] and add 3 to skip the header lines
    if re.search('\\[XSECTIONS\\]',line):
        xsections_marker=linecount+3
    if re.search('\\[CONTROLS\\]',line):
        controls_marker=linecount+1
    # Record the position of the phrase [CURVES] and add 3 to skip the header lines
    if re.search('\\[CURVES\\]',line):
        curves_marker=linecount+3
    if re.search('\\[TIMESERIES\\]',line):
        timesrs_marker=linecount+3
    # Record the position of the phrase [COORDINATES] and add 3 to skip the header lines
    if re.search('\\[COORDINATES\\]',line):
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

  lines[end_time]="END_TIME             "+str(23)+":"+str(59)+":00\n"
  lines[end_date]="END_DATE             "+start_date[0]+"/"+str(min((int(start_date[1])+n_days-1),30))+"/"+start_date[2]
  if not timestep:
      timestep = maximum_xdelta / 100
  lines[routing_step]="ROUTING_STEP         "+str(timestep)+"\n"
  lines[dimensions]="DIMENSIONS "+dimensions_line
  lines[coords_marker:coords_marker]=coordinate_section
  lines[timesrs_marker:timesrs_marker]=timesrs_pat
  lines[curves_marker:curves_marker]=curves_section
  lines[controls_marker:controls_marker]=control_rules
  lines[xsections_marker:xsections_marker]=xsections_section
  lines[outlets_marker:outlets_marker]=outlet_section
  if len(pump_section)>3:
      lines[pumps_marker:pumps_marker]=pump_section
  lines[conduits_marker:conduits_marker]=conduits_section
  lines[storage_marker:storage_marker]=storage_section
  lines[outfalls_marker:outfalls_marker]=outfall_section
  lines[junctions_marker:junctions_marker]=junctions_section

  if output.exists():
      output.unlink()
  output.write_text(''.join(lines))


def define_tank_inflow(filename:str, times:list, inflows:list):
    '''
    Define the inflow pattern for the supply tank

    Parameters
    ----------
    filename : str
        Name of the file to write the inflow pattern to
    times : list
        List of times in hours
    inflows : list
        List of inflows in LPS
    '''
    assert len(times) == len(inflows), "Times and inflows must have the same length"

    file = open(filename, 'r')
    filelines = file.readlines()
    lines=[]
    linecount=0 
    for line in filelines:
        if re.search(r"\[STORAGE\]",line):
            tank_id_place = (linecount + 3)
        if re.search(r"\[INFLOWS\]",line):
            inflows_marker=linecount+3
        if re.search(r"\[REPORT\]",line):
            report_marker=linecount-1
        linecount+=1
        lines.append(line)

    tank_id = lines[tank_id_place].split()[0] 
    inflow_line = [f"{str(tank_id)}\tFLOW\tINFLOW\tFLOW\t1.0\t1.0\t1\n"]

    inflow_pattern = [";\n"]
    for time, inflow in zip(times, inflows):
        inflow_pattern.append(f"INFLOW\t{str(time)}\t{str(inflow)}\n")

    lines[report_marker:report_marker] = inflow_pattern
    lines[inflows_marker:inflows_marker] = inflow_line

    file.close()

    with open(filename, 'w') as file:
        file.writelines(lines)
    
    return filelines

if __name__ == "__main__":
  network_file = Path("Networks/Ismail/ismail.inp")
  logging.info(f"Reading network file: {network_file}")

  desired_pressure = 10
  minimum_pressure = 0 

  # Step for control curve generation. The smaller the step, the smoother the curve
  step = 0.002 # Reduce if continuity error is too high
  # Parameters for Mauro de Marchis et al. (2015)'s float valve emitter law
  m, n= 2.5, 4
  tank_height = 1
  # Diurnal demand pattern defined as a list of multipliers for the base demand
  pattern=[0.8,0.7,0.6,0.5,0.5,0.5,0.6,0.8,1.2,1.3,1.2,1.2,1.2,1.2,1.2,1.2,1.1,1.1,1.1,1.2,1.3,1.3,1.1,1]

  # Disable adaptative discretization to reduce computational load
  maximum_xdelta = 20
  adaptative = False             

  # Initial fullness factor for storage tanks (households)
  storage_initial_fullness_factor  = 0.2

  (
      n_days,
      maximum_xdelta,
      dimensions_line,
      coordinate_section,
      timesrs_pat,
      curves_section,
      control_rules,
      xsections_section,
      outlet_section,
      pump_section,
      conduits_section,
      storage_section,
      outfall_section,
      junctions_section
  ) = load_epanet(
      inp_file=network_file,
      supply_duration_inp=8.0,
      concentric=True,
      len_to_diameter_ratio=30,
      adaptive=adaptative,
      maximum_xdelta=maximum_xdelta,
      tank_height=tank_height,
      minimum_pressure=0,
      pressure_diff=desired_pressure - minimum_pressure,
      pdw_variable="PRESSURE",
      h_min=0.9,
      h_max=tank_height,
      step=step,
      m=m,
      n=n,
      pattern=pattern,
      leak_fraction=0.1,
      n_days=1,
      reservoir_height=5,
      storage_initial_fullness_factor=storage_initial_fullness_factor,
  )

  output_file = Path("Networks/Ismail/ismail_SWMMIN.inp")
  logging.info(f"Writing SWMM file: {output_file}")

  write_swmmin(
      template=Path("Networks/Empty_SWMM_Template.inp"),
      output=output_file,
      n_days=n_days,
      maximum_xdelta=maximum_xdelta,
      timestep=None,
      dimensions_line=dimensions_line,
      coordinate_section=coordinate_section,
      timesrs_pat=timesrs_pat,
      curves_section=curves_section,
      control_rules=control_rules,
      xsections_section=xsections_section,
      outlet_section=outlet_section,
      pump_section=pump_section,
      conduits_section=conduits_section,
      storage_section=storage_section,
      outfall_section=outfall_section,
      junctions_section=junctions_section
  )

    # Define the inflow pattern for the supply tank
  times = np.arange(0, 24, 1)
  inflows = [0.0 for time in times]  ## REPLACE WITH DESIRED INFLOW PATTERN IN LPS (OR LINK IT TO YOUR SOLAR PANEL MODEL)
  
  logging.info(f"Assigning Tank Inflows: {output_file}")
  # a = define_tank_inflow(output_file, times, inflows)
  
    
