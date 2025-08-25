"""Tests for the core GOSWMMIN functionality."""

import pytest
import pathlib
from unittest.mock import patch, MagicMock

from goswmmin.core import SWMMINSimulation


class TestSWMMINSimulation:
    """Test cases for SWMMINSimulation class."""
    
    def test_init_with_valid_file(self, sample_epanet_file):
        """Test initialization with a valid EPANET file."""
        sim = SWMMINSimulation(sample_epanet_file)
        assert sim.filepath == sample_epanet_file
        assert sim.name_only == sample_epanet_file.name
    
    def test_init_with_nonexistent_file(self, temp_dir):
        """Test initialization with non-existent file raises error."""
        nonexistent_file = temp_dir / "nonexistent.inp"
        with pytest.raises(FileNotFoundError):
            SWMMINSimulation(nonexistent_file)
    
    def test_init_with_string_path(self, sample_epanet_file):
        """Test initialization with string path."""
        sim = SWMMINSimulation(str(sample_epanet_file))
        assert sim.filepath == sample_epanet_file
    
    def test_init_warns_for_wrong_extension(self, temp_dir):
        """Test that warning is issued for non-.inp files."""
        wrong_ext_file = temp_dir / "test.txt"
        wrong_ext_file.write_text("dummy content")
        
        with pytest.warns(UserWarning):
            SWMMINSimulation(wrong_ext_file)
    
    def test_convert_to_swmmin_basic_config(self, sample_epanet_file, simulation_config):
        """Test basic conversion with uniform parameters."""
        sim = SWMMINSimulation(sample_epanet_file)
        
        with patch.object(sim, '_parse_epanet_input'), \
             patch.object(sim, '_match_pipes_concentric'), \
             patch.object(sim, '_discretize'), \
             patch.object(sim, '_write_junctions'), \
             patch.object(sim, '_write_outfalls'), \
             patch.object(sim, '_write_storage_nodes'), \
             patch.object(sim, '_write_pipes'), \
             patch.object(sim, '_write_outlets'), \
             patch.object(sim, '_write_xsections'), \
             patch.object(sim, '_write_curves'), \
             patch.object(sim, '_write_controls'), \
             patch.object(sim, '_write_coords'), \
             patch.object(sim, '_write_file'):
            
            result = sim.convert_to_swmmin(**simulation_config)
            assert isinstance(result, pathlib.Path)
    
    def test_convert_to_swmmin_with_csv_inputs(self, sample_epanet_file, sample_csv_data):
        """Test conversion with CSV input files."""
        sim = SWMMINSimulation(sample_epanet_file)
        
        config = {
            'supply_duration': 8.0,
            'minimum_pressure': str(sample_csv_data['min_pressure']),
            'desired_pressure': str(sample_csv_data['des_pressure']),
            'pdw_exponent': str(sample_csv_data['pdw_exponent']),
            'consum_pattern': str(sample_csv_data['consumption_pattern'])
        }
        
        with patch.object(sim, '_parse_epanet_input'), \
             patch.object(sim, '_match_pipes_concentric'), \
             patch.object(sim, '_discretize'), \
             patch.object(sim, '_write_junctions'), \
             patch.object(sim, '_write_outfalls'), \
             patch.object(sim, '_write_storage_nodes'), \
             patch.object(sim, '_write_pipes'), \
             patch.object(sim, '_write_outlets'), \
             patch.object(sim, '_write_xsections'), \
             patch.object(sim, '_write_curves'), \
             patch.object(sim, '_write_controls'), \
             patch.object(sim, '_write_coords'), \
             patch.object(sim, '_write_file'):
            
            result = sim.convert_to_swmmin(**config)
            assert isinstance(result, pathlib.Path)
            assert sim.pdw_flag is True
            assert sim.consum_pattern_flag is True
    
    def test_convert_to_swmmin_invalid_supply_duration(self, sample_epanet_file):
        """Test that invalid supply duration raises error."""
        sim = SWMMINSimulation(sample_epanet_file)
        
        with pytest.raises(ValueError, match="Supply duration must be a positive number"):
            sim.convert_to_swmmin(
                supply_duration=-1.0,
                minimum_pressure=10.0,
                desired_pressure=20.0,
                pdw_exponent=0.5
            )
    
    def test_convert_to_swmmin_invalid_pdw_variable(self, sample_epanet_file):
        """Test that invalid PDW variable raises error."""
        sim = SWMMINSimulation(sample_epanet_file)
        
        with pytest.raises(ValueError, match="pdw_variable must be either 'PRESSURE' or 'DEPTH'"):
            sim.convert_to_swmmin(
                supply_duration=8.0,
                minimum_pressure=10.0,
                desired_pressure=20.0,
                pdw_exponent=0.5,
                pdw_variable='INVALID'
            )
    
    def test_convert_to_swmmin_mismatched_csv_inputs(self, sample_epanet_file, sample_csv_data):
        """Test that mismatched CSV inputs raise error."""
        sim = SWMMINSimulation(sample_epanet_file)
        
        with pytest.raises(ValueError, match="must all be either uniform values"):
            sim.convert_to_swmmin(
                supply_duration=8.0,
                minimum_pressure=str(sample_csv_data['min_pressure']),  # CSV
                desired_pressure=20.0,  # Float
                pdw_exponent=0.5  # Float
            )
    
    def test_run_swmmin_without_conversion(self, sample_epanet_file):
        """Test that running simulation without conversion raises error."""
        sim = SWMMINSimulation(sample_epanet_file)
        
        with pytest.raises(RuntimeError, match="Must call convert_to_swmmin"):
            sim.run_swmmin()
    
    def test_run_swmmin_with_missing_file(self, sample_epanet_file, temp_dir):
        """Test that running simulation with missing SWMM file raises error."""
        sim = SWMMINSimulation(sample_epanet_file)
        sim.swmm_file = temp_dir / "nonexistent.inp"
        
        with pytest.raises(FileNotFoundError, match="SWMM input file not found"):
            sim.run_swmmin()
    
    @patch('goswmmin.core.pyswmm.Simulation')
    def test_run_swmmin_success(self, mock_simulation, sample_epanet_file, temp_dir):
        """Test successful simulation run."""
        sim = SWMMINSimulation(sample_epanet_file)
        sim.swmm_file = temp_dir / "test.inp"
        sim.swmm_file.write_text("dummy swmm file")
        
        # Mock the simulation context manager
        mock_sim_instance = MagicMock()
        mock_sim_instance.__enter__ = MagicMock(return_value=mock_sim_instance)
        mock_sim_instance.__exit__ = MagicMock(return_value=None)
        mock_sim_instance.__iter__ = MagicMock(return_value=iter([1, 2, 3]))
        mock_simulation.return_value = mock_sim_instance
        
        # Should not raise an exception
        sim.run_swmmin()
        mock_simulation.assert_called_once()
    
    def test_get_pressures_without_simulation(self, sample_epanet_file):
        """Test getting pressures without running simulation raises error."""
        sim = SWMMINSimulation(sample_epanet_file)
        
        with pytest.raises(RuntimeError, match="Must run simulation before getting results"):
            sim.get_pressures()
    
    def test_get_tank_vols_heights_without_simulation(self, sample_epanet_file):
        """Test getting tank data without running simulation raises error."""
        sim = SWMMINSimulation(sample_epanet_file)
        
        with pytest.raises(RuntimeError, match="Must run simulation before getting results"):
            sim.get_tank_vols_heights()
    
    def test_backward_compatibility(self, sample_epanet_file):
        """Test that SWMMIN_sim alias works for backward compatibility."""
        from goswmmin.core import SWMMIN_sim
        
        sim = SWMMIN_sim(sample_epanet_file)
        assert isinstance(sim, SWMMINSimulation)
        assert sim.filepath == sample_epanet_file


class TestParameterValidation:
    """Test parameter validation methods."""
    
    def test_validate_csv_file_missing(self, sample_epanet_file, temp_dir):
        """Test validation of missing CSV file."""
        sim = SWMMINSimulation(sample_epanet_file)
        
        with pytest.raises(FileNotFoundError):
            sim._validate_csv_file(str(temp_dir / "missing.csv"), "Test file")
    
    def test_read_csv_column(self, sample_epanet_file, sample_csv_data):
        """Test reading CSV column data."""
        sim = SWMMINSimulation(sample_epanet_file)
        
        result = sim._read_csv_column(str(sample_csv_data['min_pressure']), 'Min_Pressure')
        assert isinstance(result, list)
        assert len(result) == 3
        assert all(isinstance(x, (int, float)) for x in result)
