"""
Basic usage example for GOSWMMIN package.

This example demonstrates how to create and run a SWMMIN simulation
with various configuration options.
"""
#%%
import pathlib
from goswmmin import SWMMINSimulation


def main():
    """Run the basic GOSWMMIN example."""
    
    # Path to the example EPANET network file
    # Note: Update this path to point to your actual network file
    network_file = pathlib.Path('../goswmmin/resources/Linear_Network.inp')
    
    if not network_file.exists():
        print(f"Network file not found: {network_file}")
        print("Please update the path to point to a valid EPANET .inp file")
        return
    
    # Create a SWMMIN simulation object
    print("Creating SWMMIN simulation...")
    sim = SWMMINSimulation(network_file)
    
    # Example 1: Basic configuration with uniform parameters
    print("\n=== Example 1: Basic Configuration ===")
    try:
        swmm_file = sim.convert_to_swmmin(
            supply_duration=8.0,          # 8-hour supply period
            minimum_pressure=10.0,        # 10 m minimum pressure for all nodes
            desired_pressure=20.0,        # 20 m desired pressure for all nodes
            pdw_exponent=0.5,             # PDW exponent (Wagner et al., 1989)
            n_days=1,                     # Simulate for 1 day
            length_to_diameter=30.0,      # L/D ratio for discretization
            solution_speed=100.0,         # Solution speed in m/s
            leak_fraction=0.1,            # 10% leakage
            tank_heights=1.0,             # 1 m tank height for all nodes
            pdw_variable='PRESSURE'       # Use pressure for PDW
        )
        print(f"SWMM input file created: {swmm_file}")
        
        # Run the simulation
        print("Running simulation...")
        sim.run_swmmin()
        print("Simulation completed successfully!")
        
    except Exception as e:
        print(f"Error in basic configuration: {e}")
    
    # Example 2: Advanced configuration with CSV inputs
    print("\n=== Example 2: Advanced Configuration with CSV Files ===")
    
    # Define paths to CSV files (relative to the resources directory)
    resource_dir = pathlib.Path('../goswmmin/resources')
    csv_files = {
        'min_pressure': resource_dir / 'min_pressure.csv',
        'des_pressure': resource_dir / 'des_pressure.csv', 
        'pdw_exponent': resource_dir / 'pdw_exponent.csv',
        'q_des': resource_dir / 'q_des.csv',
        'tank_areas': resource_dir / 'tank_areas.csv',
        'tank_heights': resource_dir / 'tank_heights.csv',
        'consum_pattern': resource_dir / 'consum_pattern.csv'
    }
    
    # Check if CSV files exist
    missing_files = [name for name, path in csv_files.items() if not path.exists()]
    
    if missing_files:
        print(f"CSV files not found: {missing_files}")
        print("Skipping advanced example. Please ensure all CSV files are available.")
    else:
        try:
            # Create simulation with CSV inputs
            sim2 = SWMMINSimulation(network_file)
            
            swmm_file = sim2.convert_to_swmmin(
                supply_duration=6.0,                                    # 6-hour supply
                minimum_pressure=str(csv_files['min_pressure']),        # Node-specific min pressure
                desired_pressure=str(csv_files['des_pressure']),        # Node-specific desired pressure
                pdw_exponent=str(csv_files['pdw_exponent']),            # Node-specific PDW exponent
                q_des=str(csv_files['q_des']),                          # Node-specific desired flow
                tank_areas=str(csv_files['tank_areas']),                # Node-specific tank areas
                tank_heights=str(csv_files['tank_heights']),            # Node-specific tank heights
                consum_pattern=str(csv_files['consum_pattern']),        # 24-hour consumption pattern
                adaptive_disc=True,                                     # Use adaptive discretization
                solution_speed=150.0,                                   # Higher solution speed
                n_days=2                                                # 2-day simulation
            )
            print(f"Advanced SWMM input file created: {swmm_file}")
            
            # Run the simulation
            print("Running advanced simulation...")
            sim2.run_swmmin()
            print("Advanced simulation completed successfully!")
            
            # Get results
            print("\nExtracting results...")
            try:
                # Get pressure results for specific nodes
                pressures = sim2.get_pressures(specific_nodes=['DN1', 'DN2'])
                print(f"Pressure results shape: {pressures.shape}")
                
                # Get tank volume and height results
                tank_vols, tank_heights = sim2.get_tank_vols_heights(specific_nodes=['DN1', 'DN2'])
                print(f"Tank volumes shape: {tank_vols.shape}")
                print(f"Tank heights shape: {tank_heights.shape}")
                
            except Exception as e:
                print(f"Note: Results extraction not fully implemented: {e}")
            
        except Exception as e:
            print(f"Error in advanced configuration: {e}")
    
    print("\n=== Example Complete ===")


def create_sample_csv_files():
    """
    Create sample CSV files for the advanced example.
    
    This function creates example CSV files that can be used with the package.
    In practice, these would be prepared based on your specific network requirements.
    """
    import pandas as pd
    import numpy as np
    
    resource_dir = pathlib.Path('../goswmmin/resources')
    resource_dir.mkdir(exist_ok=True)
    
    # Sample node IDs (adjust based on your network)
    nodes = ['DN1', 'DN2', 'DN3', 'DN4']
    
    # Minimum pressure (m)
    min_pressure_df = pd.DataFrame({
        'ID': nodes,
        'Min_Pressure': [5.0, 7.0, 8.0, 6.0]
    })
    min_pressure_df.to_csv(resource_dir / 'min_pressure.csv', header=False, index=False)
    
    # Desired pressure (m)
    des_pressure_df = pd.DataFrame({
        'ID': nodes,
        'Des_Pressure': [15.0, 17.0, 18.0, 16.0]
    })
    des_pressure_df.to_csv(resource_dir / 'des_pressure.csv', header=False, index=False)
    
    # PDW exponent
    pdw_exponent_df = pd.DataFrame({
        'ID': nodes,
        'Exponent': [0.5, 0.5, 0.4, 0.5]
    })
    pdw_exponent_df.to_csv(resource_dir / 'pdw_exponent.csv', header=False, index=False)
    
    # Desired flow rates (m³/s)
    q_des_df = pd.DataFrame({
        'ID': nodes,
        'Q_des': [0.001, 0.0015, 0.002, 0.0012]
    })
    q_des_df.to_csv(resource_dir / 'q_des.csv', header=False, index=False)
    
    # Tank areas (m²)
    tank_areas_df = pd.DataFrame({
        'ID': nodes,
        'Area': [10.0, 15.0, 20.0, 12.0]
    })
    tank_areas_df.to_csv(resource_dir / 'tank_areas.csv', header=False, index=False)
    
    # Tank heights (m)
    tank_heights_df = pd.DataFrame({
        'ID': nodes,
        'Height': [1.5, 1.2, 1.8, 1.0]
    })
    tank_heights_df.to_csv(resource_dir / 'tank_heights.csv', header=False, index=False)
    
    # 24-hour consumption pattern
    pattern_df = pd.DataFrame({
        'Hour': range(24),
        'Multiplier': [
            0.5, 0.4, 0.3, 0.3, 0.4, 0.6, 0.8, 1.0,
            1.2, 1.0, 0.9, 1.0, 1.1, 1.0, 0.9, 1.0,
            1.2, 1.4, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6
        ]
    })
    pattern_df.to_csv(resource_dir / 'consum_pattern.csv', header=False, index=False)
    
    print(f"Sample CSV files created in: {resource_dir}")


if __name__ == "__main__":
    # Uncomment the line below to create sample CSV files
    # create_sample_csv_files()
    
    main()