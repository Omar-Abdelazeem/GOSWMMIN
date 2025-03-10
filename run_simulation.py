import pyswmm
import re
import pandas as pd
import numpy as np
import datetime
import matplotlib
import matplotlib.pyplot as plt
import logging
from pathlib import Path

logger=logging.getLogger()
logger.setLevel(logging.DEBUG)
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s')

def run_simulation(input_file: Path,output_file: Path):
  sim=pyswmm.Simulation(inputfile=input_file.as_posix(), outputfile=output_file.as_posix())
  nodes=pyswmm.nodes.Nodes(sim)

  tankids=[]
  for node in nodes:
      if re.search('StorageforNode',node.nodeid):
          tankids.append(node.nodeid)
  demand_node_ids=[x[14:] for x in tankids]
  Withdraw_ids = ["Outlet" + x for x in demand_node_ids]
  Consume_ids = ["DemandOutlet" + x for x in demand_node_ids]

  step_count=0       
  with sim as sim:
      logger.info(f"Simulation started")
      system_routing = pyswmm.SystemStats(sim)
      for step in sim:
          if step_count%1000==0:
              logging.debug(f"Current Simulation Time is >> {sim.current_time}, {round(sim.percent_complete*100,1)}% Complete")
          step_count+=1
          pass
      sim._model.swmm_end()

      # Diagnostic metric for how good sim was in terms of mass balance in percentage
      # Want < 5%
      # TODO: reduce time step for more accurate results
      logging.error(f"Continuity Error: {sim.flow_routing_error}") 
      

  return Withdraw_ids,Consume_ids


def analyze_demand_satisfaction(
    output_file: Path,
    demand_volumes_file: Path,
    Withdraw_ids: list,
    Consume_ids: list,
):
  WithdrawRates=pd.DataFrame()   # water depth in tanks
  ConsumeRates=pd.DataFrame()   #water depth in nodes
  swtch=True                   # switch variable for upcoming condition

  with pyswmm.Output(output_file.as_posix()) as out:
      for link in out.links:
          # One time only. Gets the timesteps (the keys in the output series dictionary) and stores them to be used as index
          if swtch:
          # node_series produces a dictionary with the keys corresponding to timestamps and values contain the value of the selected variable (FLOW_RATE) at each timestamp
              index = pd.Series(out.link_series(link, 'Flow_rate').keys())
              WithdrawRates.loc[:,"time"]=index
              ConsumeRates.loc[:,"time"]=index
              swtch=False
          # If node id is in the prepared list of demand nodes (tanks)
          if link in Withdraw_ids:
              WithdrawRates.loc[:,link[6:]] = out.link_series(link, 'Flow_rate').values()  
          elif link in Consume_ids:
              ConsumeRates.loc[:,link[12:]]= out.link_series(link, 'Flow_rate').values()

  reporting_timestep = index.diff().dt.seconds.mode()[0] #reporting time step in seconds
  Withdraw_by_Day = WithdrawRates.groupby(pd.Grouper(key = 'time', freq = 'D'))
  Consumer_by_Day = ConsumeRates.groupby(pd.Grouper(key = 'time', freq = 'D'))
  Daily_Totals_W = WithdrawRates.groupby(pd.Grouper(key = 'time', freq = 'D')).sum().T*reporting_timestep / 1000
  Daily_Totals_C = ConsumeRates.groupby(pd.Grouper(key = 'time', freq = 'D')).sum().T*reporting_timestep / 1000

  demand_volumes = pd.read_csv(demand_volumes_file).set_index(Daily_Totals_W.index)
  Daily_Satis_W = pd.DataFrame(columns=Daily_Totals_W.columns, index = Daily_Totals_W.index)
  Daily_Satis_C = pd.DataFrame(columns = Daily_Totals_C.columns, index = Daily_Totals_C.index)
  for column in range(len(Daily_Totals_W.columns)):
          Daily_Satis_W.iloc[:,column] = Daily_Totals_W.iloc[:,column] / demand_volumes['Volume']*100
  for column in range(len(Daily_Totals_C.columns)):
          Daily_Satis_C.iloc[:,column] = Daily_Totals_C.iloc[:,column] / demand_volumes['Volume']*100

  fig, ax  = plt.subplots()
  ax.bar(Daily_Satis_W.index, Daily_Satis_W.mean(axis=1), label = 'Withdrawal Satisfaction')
  ax.set_xlabel('Demand Node')
  ax.set_ylabel('Satisfaction Consumption-based (%)')
  ax.set_ylim(0,100)
  plt.savefig(output_file.parent / "Withdrawal_Satisfaction.png")



if __name__ == "__main__":
  swmmin_file = Path("Networks/Ismail/ismail_SWMMIN.inp")
  output_file = Path("Networks/Ismail/ismail_SWMMIN.out")

  logging.info(f"Inp file: {swmmin_file}")
  logging.info(f"Output file: {output_file}")

  Withdraw_ids,Consume_ids = run_simulation(swmmin_file,output_file)
  analyze_demand_satisfaction(output_file, Path("Networks/Ismail/ismail.inp_DemandVolumes.csv"),Withdraw_ids,Consume_ids)

  

