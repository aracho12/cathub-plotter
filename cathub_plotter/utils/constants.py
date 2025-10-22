"""
Physical constants and default parameters
"""

class PhysicalConstants:
    """Physical constants in SI units"""
    h = 6.62607015e-34  # Planck constant (J⋅s)
    c = 299792458       # Speed of light (m/s)
    e = 1.602176634e-19 # Elementary charge (C)
    kb = 1.380649e-23   # Boltzmann constant (J/K)
    Na = 6.02214076e23  # Avogadro constant (mol^-1)
    R = 8.314462618     # Gas constant (J/(mol⋅K))
    
    # Conversion factors
    eV_to_J = 1.602176634e-19  # eV to Joules
    J_to_eV = 1 / eV_to_J      # Joules to eV
    cm1_to_eV = 1.239841984e-4 # cm^-1 to eV
    eV_to_cm1 = 1 / cm1_to_eV  # eV to cm^-1

# Default fugacity values in Pa
DEFAULT_FUGACITY_DICT = {
    'H2_g': 101325,
    'CO2_g': 101325,
    'CO_g': 5562,
    'HCOOH_g': 2,
    'CH3OH_g': 6079,
    'H2O_g': 3534,
    'CH4_g': 20467,
}

# Default solvation energy corrections (in eV)
DEFAULT_SOLVATION_DICT = {
    'COH*': -0.25,
    'COOH*': -0.25,
    'CO2_g': +0.4,  
    'H2_g': +0.09,  
}

PHYSICAL_CONSTANTS = PhysicalConstants()
