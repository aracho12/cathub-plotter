"""
Cathub-Plotter: Free Energy Diagram Plotting for Catalysis Research

A Python package for plotting free energy diagrams from catalysis-hub.org data
and MKM (Microkinetic Model) files.
"""

__version__ = "0.1.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

# Import main classes for easy access
from .core.calculator import FreeEnergyCalculator
from .core.mkm_calculator import MKMEnergyCalculator
from .plotters.diagram import FreeEnergyDiagramPlotter
from .parsers.cathub import CatalysisHubParser
from .parsers.mkm import MKMFileParser
from .parsers.input import InputFileParser

# Import utility functions
from .utils.thermodynamics import calculate_thermo_correction

__all__ = [
    "FreeEnergyCalculator",
    "MKMEnergyCalculator",
    "FreeEnergyDiagramPlotter", 
    "CatalysisHubParser",
    "MKMFileParser",
    "InputFileParser",
    "calculate_thermo_correction",
]
