#!/usr/bin/env python3
"""
Comprehensive test of the GOSWMMIN package workflow.
Tests all major functionality including result extraction methods.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '.'))

from goswmmin import SWMMINSimulation
import pandas as pd
import numpy as np

def test_full_workflow():
    """Test the complete GOSWMMIN simulation workflow."""
    
    print("=== GOSWMMIN Full Workflow Test ===\n")
    
    # 1. Get network path
    print("1. Locating network file...")
    network_path = os.path.join(os.path.dirname(__file__), 'goswmmin', 'resources', 'Linear_Network.inp')
    if not os.path.exists(network_path):
        print(f"✗ Network file not found: {network_path}")
        return False
    print(f"✓ Network file found: {network_path}")
    
    # 2. Create simulation instance with network file
    print("\n2. Creating SWMMIN simulation instance...")
    try:
        sim = SWMMINSimulation(network_path)
        print("✓ Simulation instance created successfully")
    except Exception as e:
        print(f"✗ Failed to create simulation: {e}")
        return False
    
    # 3. Convert to SWMMIN format and run simulation
    print("\n3. Converting to SWMMIN and running simulation...")
    try:
        # First convert the network to SWMMIN format
        swmmin_file = sim.convert_to_swmmin(
            supply_duration=8.0,
            minimum_pressure=10.0,
            desired_pressure=20.0,
            pdw_exponent=0.5,
            n_days=1,
            timestep=600  # 10 minutes - larger timestep to avoid SWMM error
        )
        print(f"✓ Network converted to SWMMIN format: {swmmin_file}")
        
        # Run the SWMMIN simulation
        sim.run_swmmin()
        print("✓ Simulation completed successfully")
    except Exception as e:
        print(f"✗ Simulation failed: {e}")
        # Try with default timestep if custom one fails
        if "routing time step" in str(e):
            print("   Trying with default timestep...")
            try:
                swmmin_file = sim.convert_to_swmmin(
                    supply_duration=8.0,
                    minimum_pressure=10.0,
                    desired_pressure=20.0,
                    pdw_exponent=0.5,
                    n_days=1
                    # Use default timestep
                )
                sim.run_swmmin()
                print("✓ Simulation completed successfully with default timestep")
            except Exception as e2:
                print(f"✗ Simulation still failed: {e2}")
                return False
        else:
            import traceback
            traceback.print_exc()
            return False
    
    # 4. Test result extraction methods
    print("\n4. Testing result extraction methods...")
    
    # Test pressures
    print("\n4a. Testing pressure extraction...")
    try:
        pressures = sim.get_pressures()
        print(f"✓ Pressures extracted: {type(pressures)} with shape {pressures.shape if hasattr(pressures, 'shape') else len(pressures)}")
        if hasattr(pressures, 'head'):
            print(f"   First few columns: {list(pressures.columns[:5])}")
        else:
            print(f"   Sample data: {pressures[:3] if len(pressures) > 3 else pressures}")
    except Exception as e:
        print(f"✗ Failed to extract pressures: {e}")
    
    # Test tank volumes
    print("\n4b. Testing tank volume extraction...")
    try:
        tank_volumes = sim.get_tank_volumes()
        print(f"✓ Tank volumes extracted: {type(tank_volumes)} with shape {tank_volumes.shape if hasattr(tank_volumes, 'shape') else len(tank_volumes)}")
        if hasattr(tank_volumes, 'head'):
            print(f"   First few columns: {list(tank_volumes.columns[:5])}")
        else:
            print(f"   Sample data: {tank_volumes[:3] if len(tank_volumes) > 3 else tank_volumes}")
    except Exception as e:
        print(f"✗ Failed to extract tank volumes: {e}")
    
    # Test tank heights
    print("\n4c. Testing tank height extraction...")
    try:
        tank_heights = sim.get_tank_heights()
        print(f"✓ Tank heights extracted: {type(tank_heights)} with shape {tank_heights.shape if hasattr(tank_heights, 'shape') else len(tank_heights)}")
        if hasattr(tank_heights, 'head'):
            print(f"   First few columns: {list(tank_heights.columns[:5])}")
        else:
            print(f"   Sample data: {tank_heights[:3] if len(tank_heights) > 3 else tank_heights}")
    except Exception as e:
        print(f"✗ Failed to extract tank heights: {e}")
    
    # Test withdrawals
    print("\n4d. Testing withdrawal extraction...")
    try:
        withdrawals = sim.get_withdrawals()
        print(f"✓ Withdrawals extracted: {type(withdrawals)} with shape {withdrawals.shape if hasattr(withdrawals, 'shape') else len(withdrawals)}")
        if hasattr(withdrawals, 'head'):
            print(f"   First few columns: {list(withdrawals.columns[:5])}")
        else:
            print(f"   Sample data: {withdrawals[:3] if len(withdrawals) > 3 else withdrawals}")
    except Exception as e:
        print(f"✗ Failed to extract withdrawals: {e}")
    
    # Test consumptions
    print("\n4e. Testing consumption extraction...")
    try:
        consumptions = sim.get_consumptions()
        print(f"✓ Consumptions extracted: {type(consumptions)} with shape {consumptions.shape if hasattr(consumptions, 'shape') else len(consumptions)}")
        if hasattr(consumptions, 'head'):
            print(f"   First few columns: {list(consumptions.columns[:5])}")
        else:
            print(f"   Sample data: {consumptions[:3] if len(consumptions) > 3 else consumptions}")
    except Exception as e:
        print(f"✗ Failed to extract consumptions: {e}")
    
    # 5. Test simulation info
    print("\n5. Testing simulation information...")
    try:
        info = sim.get_simulation_info()
        print(f"✓ Simulation info: {info}")
    except Exception as e:
        print(f"✗ Failed to get simulation info: {e}")
    
    # 6. Test backward compatibility
    print("\n6. Testing backward compatibility...")
    try:
        from goswmmin import SWMMIN_sim
        legacy_sim = SWMMIN_sim(network_path)
        print("✓ SWMMIN_sim (legacy) import works")
        
        # Test that it's the same class
        if type(legacy_sim) == type(sim):
            print("✓ SWMMIN_sim is properly aliased to SWMMINSimulation")
        else:
            print(f"⚠ SWMMIN_sim type: {type(legacy_sim)}, SWMMINSimulation type: {type(sim)}")
    except Exception as e:
        print(f"✗ Legacy import failed: {e}")
    
    print("\n=== Test Complete ===")
    return True

if __name__ == "__main__":
    test_full_workflow()
