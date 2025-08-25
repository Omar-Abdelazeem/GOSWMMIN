"""Tests for utility functions."""

import pytest
import pathlib
import pandas as pd
from unittest.mock import patch

from goswmmin.utils import (
    validate_input_file,
    load_csv_data,
    create_output_filename,
    format_duration,
    validate_numeric_range,
    ConfigValidator
)


class TestUtilityFunctions:
    """Test utility functions."""
    
    def test_validate_input_file_success(self, sample_epanet_file):
        """Test successful file validation."""
        result = validate_input_file(sample_epanet_file)
        assert result == sample_epanet_file
        assert isinstance(result, pathlib.Path)
    
    def test_validate_input_file_missing(self, temp_dir):
        """Test validation of missing file."""
        missing_file = temp_dir / "missing.txt"
        with pytest.raises(FileNotFoundError):
            validate_input_file(missing_file)
    
    def test_validate_input_file_directory(self, temp_dir):
        """Test validation fails for directory."""
        with pytest.raises(ValueError, match="Path is not a file"):
            validate_input_file(temp_dir)
    
    def test_load_csv_data_success(self, sample_csv_data):
        """Test successful CSV loading."""
        df = load_csv_data(sample_csv_data['min_pressure'])
        assert isinstance(df, pd.DataFrame)
        assert len(df.columns) == 2
        assert len(df) == 3
    
    def test_load_csv_data_with_expected_columns(self, sample_csv_data):
        """Test CSV loading with column count validation."""
        df = load_csv_data(sample_csv_data['min_pressure'], expected_columns=2)
        assert len(df.columns) == 2
    
    def test_load_csv_data_wrong_column_count(self, sample_csv_data):
        """Test CSV loading fails with wrong column count."""
        with pytest.raises(ValueError, match="expected 3"):
            load_csv_data(sample_csv_data['min_pressure'], expected_columns=3)
    
    def test_create_output_filename_default(self, temp_dir):
        """Test output filename creation with default suffix."""
        input_file = temp_dir / "network.inp"
        result = create_output_filename(input_file)
        expected = temp_dir / "network_SWMMIN.inp"
        assert result == expected
    
    def test_create_output_filename_custom_suffix(self, temp_dir):
        """Test output filename creation with custom suffix."""
        input_file = temp_dir / "network.inp"
        result = create_output_filename(input_file, suffix="_custom")
        expected = temp_dir / "network_custom.inp"
        assert result == expected
    
    def test_format_duration_minutes(self):
        """Test duration formatting for minutes."""
        assert format_duration(0.5) == "30 minutes"
        assert format_duration(0.25) == "15 minutes"
    
    def test_format_duration_whole_hours(self):
        """Test duration formatting for whole hours."""
        assert format_duration(1.0) == "1 hours"
        assert format_duration(24.0) == "24 hours"
    
    def test_format_duration_mixed(self):
        """Test duration formatting for mixed hours and minutes."""
        assert format_duration(1.5) == "1 hours 30 minutes"
        assert format_duration(2.25) == "2 hours 15 minutes"
    
    def test_validate_numeric_range_success(self):
        """Test successful numeric range validation."""
        # Should not raise any exception
        validate_numeric_range(5.0, 1.0, 10.0, "test_param")
    
    def test_validate_numeric_range_below_min(self):
        """Test numeric range validation fails below minimum."""
        with pytest.raises(ValueError, match="test_param must be between"):
            validate_numeric_range(0.5, 1.0, 10.0, "test_param")
    
    def test_validate_numeric_range_above_max(self):
        """Test numeric range validation fails above maximum."""
        with pytest.raises(ValueError, match="test_param must be between"):
            validate_numeric_range(15.0, 1.0, 10.0, "test_param")


class TestConfigValidator:
    """Test configuration validation class."""
    
    def test_validate_supply_duration_success(self):
        """Test successful supply duration validation."""
        ConfigValidator.validate_supply_duration(8.0)
        ConfigValidator.validate_supply_duration(12)
    
    def test_validate_supply_duration_wrong_type(self):
        """Test supply duration validation fails for wrong type."""
        with pytest.raises(TypeError, match="Supply duration must be a number"):
            ConfigValidator.validate_supply_duration("8.0")
    
    def test_validate_supply_duration_negative(self):
        """Test supply duration validation fails for negative value."""
        with pytest.raises(ValueError, match="Supply duration must be positive"):
            ConfigValidator.validate_supply_duration(-1.0)
    
    def test_validate_supply_duration_too_large(self):
        """Test supply duration validation fails for value > 24 hours."""
        with pytest.raises(ValueError, match="Supply duration cannot exceed 24 hours"):
            ConfigValidator.validate_supply_duration(25.0)
    
    def test_validate_pressure_values_success(self):
        """Test successful pressure values validation."""
        ConfigValidator.validate_pressure_values(10.0, 20.0)
    
    def test_validate_pressure_values_negative_min(self):
        """Test pressure validation fails for negative minimum."""
        with pytest.raises(ValueError, match="Minimum pressure cannot be negative"):
            ConfigValidator.validate_pressure_values(-5.0, 20.0)
    
    def test_validate_pressure_values_desired_not_greater(self):
        """Test pressure validation fails when desired <= minimum."""
        with pytest.raises(ValueError, match="Desired pressure must be greater"):
            ConfigValidator.validate_pressure_values(20.0, 15.0)
    
    def test_validate_pdw_exponent_success(self):
        """Test successful PDW exponent validation."""
        ConfigValidator.validate_pdw_exponent(0.5)
        ConfigValidator.validate_pdw_exponent(1.0)
    
    def test_validate_pdw_exponent_out_of_range(self):
        """Test PDW exponent validation fails for out of range values."""
        with pytest.raises(ValueError, match="PDW exponent must be between"):
            ConfigValidator.validate_pdw_exponent(0.05)
        
        with pytest.raises(ValueError, match="PDW exponent must be between"):
            ConfigValidator.validate_pdw_exponent(3.0)
    
    def test_validate_solution_speed_success(self):
        """Test successful solution speed validation."""
        ConfigValidator.validate_solution_speed(100.0)
        ConfigValidator.validate_solution_speed(50.0)
    
    def test_validate_solution_speed_out_of_range(self):
        """Test solution speed validation fails for out of range values."""
        with pytest.raises(ValueError, match="Solution speed must be between"):
            ConfigValidator.validate_solution_speed(5.0)
        
        with pytest.raises(ValueError, match="Solution speed must be between"):
            ConfigValidator.validate_solution_speed(2000.0)
