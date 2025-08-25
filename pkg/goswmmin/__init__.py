"""
GOSWMMIN: Guided Operation of SWMM for Intermittent Networks

A Python package for SWMM-based simulations of Intermittent Water Supply Networks.
"""

__version__ = "0.2.0"
__author__ = "Omar Abdelazeem"
__email__ = "omaraliamer98@gmail.com"

from .core import SWMMINSimulation, SWMMIN_sim

# Define what gets imported with "from goswmmin import *"
__all__ = [
    "SWMMINSimulation",
    "SWMMIN_sim",  # For backward compatibility
]

# Package metadata
__title__ = "goswmmin"
__description__ = "Guided Operation of SWMM for Intermittent Networks"
__url__ = "https://github.com/omaraliamer/GOSWMMIN"
__license__ = "MIT"