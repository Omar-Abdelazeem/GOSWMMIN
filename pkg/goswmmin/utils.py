"""Utility functions for GOSWMMIN package."""

import pathlib
from typing import Union, Optional, Dict, Any
import pandas as pd


def validate_input_file(filepath: Union[str, pathlib.Path]) -> pathlib.Path:
    """
    Validate that an input file exists and is readable.
    
    Args:
        filepath: Path to the input file
        
    Returns:
        Validated pathlib.Path object
        
    Raises:
        FileNotFoundError: If file doesn't exist
        PermissionError: If file is not readable
    """
    path = pathlib.Path(filepath)
    
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
        
    if not path.is_file():
        raise ValueError(f"Path is not a file: {filepath}")
        
    try:
        # Test readability
        with open(path, 'r') as f:
            f.read(1)
    except PermissionError:
        raise PermissionError(f"Cannot read file: {filepath}")
        
    return path


def load_csv_data(filepath: Union[str, pathlib.Path], 
                  expected_columns: Optional[int] = None) -> pd.DataFrame:
    """
    Load and validate CSV data.
    
    Args:
        filepath: Path to CSV file
        expected_columns: Expected number of columns (optional)
        
    Returns:
        Loaded DataFrame
        
    Raises:
        ValueError: If CSV format is invalid
    """
    path = validate_input_file(filepath)
    
    try:
        df = pd.read_csv(path, header=None)
    except Exception as e:
        raise ValueError(f"Failed to read CSV file {filepath}: {str(e)}")
        
    if expected_columns and len(df.columns) != expected_columns:
        raise ValueError(
            f"CSV file {filepath} has {len(df.columns)} columns, "
            f"expected {expected_columns}"
        )
        
    return df


def get_package_resource_path(resource_name: str) -> pathlib.Path:
    """
    Get the path to a package resource file.
    
    Args:
        resource_name: Name of the resource file
        
    Returns:
        Path to the resource file
        
    Raises:
        FileNotFoundError: If resource doesn't exist
    """
    # Get the package directory
    package_dir = pathlib.Path(__file__).parent
    resource_path = package_dir / "resources" / resource_name
    
    if not resource_path.exists():
        raise FileNotFoundError(f"Resource not found: {resource_name}")
        
    return resource_path


def create_output_filename(input_path: Union[str, pathlib.Path], 
                          suffix: str = "_SWMMIN") -> pathlib.Path:
    """
    Create an output filename based on input filename.
    
    Args:
        input_path: Path to input file
        suffix: Suffix to add to filename
        
    Returns:
        Path for output file
    """
    input_path = pathlib.Path(input_path)
    stem = input_path.stem
    extension = input_path.suffix
    output_name = f"{stem}{suffix}{extension}"
    return input_path.parent / output_name


def format_duration(hours: float) -> str:
    """
    Format duration in hours to human-readable string.
    
    Args:
        hours: Duration in hours
        
    Returns:
        Formatted duration string
    """
    if hours < 1:
        minutes = int(hours * 60)
        return f"{minutes} minutes"
    elif hours == int(hours):
        return f"{int(hours)} hours"
    else:
        whole_hours = int(hours)
        minutes = int((hours - whole_hours) * 60)
        return f"{whole_hours} hours {minutes} minutes"


def validate_numeric_range(value: float, min_val: float, max_val: float, 
                          name: str) -> None:
    """
    Validate that a numeric value is within a specified range.
    
    Args:
        value: Value to validate
        min_val: Minimum allowed value
        max_val: Maximum allowed value  
        name: Name of the parameter for error messages
        
    Raises:
        ValueError: If value is outside the valid range
    """
    if not min_val <= value <= max_val:
        raise ValueError(
            f"{name} must be between {min_val} and {max_val}, got {value}"
        )


def merge_default_config(user_config: Dict[str, Any], 
                        default_config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Merge user configuration with default values.
    
    Args:
        user_config: User-provided configuration
        default_config: Default configuration values
        
    Returns:
        Merged configuration
    """
    merged = default_config.copy()
    merged.update(user_config)
    return merged


class ConfigValidator:
    """Utility class for validating simulation configuration."""
    
    @staticmethod
    def validate_supply_duration(duration: float) -> None:
        """Validate supply duration parameter."""
        if not isinstance(duration, (int, float)):
            raise TypeError("Supply duration must be a number")
        if duration <= 0:
            raise ValueError("Supply duration must be positive")
        if duration > 24:
            raise ValueError("Supply duration cannot exceed 24 hours")
    
    @staticmethod  
    def validate_pressure_values(min_pressure: float, desired_pressure: float) -> None:
        """Validate pressure parameters."""
        if min_pressure < 0:
            raise ValueError("Minimum pressure cannot be negative")
        if desired_pressure <= min_pressure:
            raise ValueError("Desired pressure must be greater than minimum pressure")
    
    @staticmethod
    def validate_pdw_exponent(exponent: float) -> None:
        """Validate PDW exponent parameter."""
        validate_numeric_range(exponent, 0.1, 2.0, "PDW exponent")
    
    @staticmethod
    def validate_solution_speed(speed: float) -> None:
        """Validate solution speed parameter."""
        validate_numeric_range(speed, 10, 1000, "Solution speed")
