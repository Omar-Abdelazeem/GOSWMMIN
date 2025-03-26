import matplotlib.pyplot as plt
from pathlib import Path
from run_simulation import OUTPUT_DIR
import pandas as pd
import pyswmm
from swmm.toolkit.shared_enum import NodeAttribute

import logging


logger=logging.getLogger()
logger.setLevel(logging.DEBUG)
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s')


def plot_flows(sim_output_file: Path):
    Withdraw_ids = ['OutletDN1', 'OutletDN2', 'OutletDN3', 'OutletDN4', 'OutletDN5', 'OutletDN6', 'OutletDN7', 'OutletDN8', 'OutletDN9']
    Consume_ids = ['DemandOutletDN1', 'DemandOutletDN2', 'DemandOutletDN3', 'DemandOutletDN4', 'DemandOutletDN5', 'DemandOutletDN6', 'DemandOutletDN7', 'DemandOutletDN8', 'DemandOutletDN9']

    WithdrawRates = pd.DataFrame()  # water depth in tanks
    ConsumeRates = pd.DataFrame()  # water depth in nodes
    TankVolumes = pd.DataFrame()  # water volume in tanks

    with pyswmm.Output(sim_output_file.as_posix()) as out:
        # Extract the index (timestamps) only once
        index = None
        for link in out.links:
            if index is None:  # One-time operation to get the timestamps
                index = pd.Series(out.link_series(link, 'Flow_rate').keys())
                WithdrawRates.loc[:, "time"] = index
                ConsumeRates.loc[:, "time"] = index

            # If node id is in the prepared list of demand nodes (tanks)
            if link in Withdraw_ids:
                WithdrawRates.loc[:, link[6:]] = out.link_series(link, 'Flow_rate').values()  
            elif link in Consume_ids:
                ConsumeRates.loc[:, link[12:]] = out.link_series(link, 'Flow_rate').values()

        index = None
        for node in out.nodes:
            if "Storage" in node:
                if index is None:
                    index = pd.Series(out.node_series(node, NodeAttribute.PONDED_VOLUME).keys())
                    TankVolumes.loc[:, "time"] = index
                TankVolumes.loc[:, node[14:]] = out.node_series(node, NodeAttribute.PONDED_VOLUME).values()

    fig, ax = plt.subplots(len(WithdrawRates.columns) - 1, 1, figsize=(12, 18))
    for i, column in enumerate(WithdrawRates.columns[1:]):
        ax[i].plot(WithdrawRates['time'], WithdrawRates[column], label=f"Withdraw {column}", color='blue')
        if column in ConsumeRates.columns[1:]:
            ax[i].plot(ConsumeRates['time'], ConsumeRates[column], label=f"Consume {column}", color='red')
        ax[i].set_title(f"Flow Rates - {column}")
        ax[i].set_xlabel("Time")
        ax[i].set_ylabel("Flow Rate")
        ax[i].legend(loc='upper left')

        if column in TankVolumes.columns:
            ax2 = ax[i].twinx()  # Create a twin Axes sharing the x-axis
            ax2.plot(TankVolumes['time'], TankVolumes[column], label=f"Volume {column}", color='green')
            ax2.set_ylabel("Tank Volume")
            ax2.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "flow_rates.png")



if __name__ == "__main__":
    sim_output_file = OUTPUT_DIR / "ismail_SWMMIN.out"
    logger.info(f"Output file: {sim_output_file}")

    plot_flows(sim_output_file)