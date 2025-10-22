"""
Thermodynamic calculations for species and reactions
"""
from typing import List, Dict, Any
import numpy as np

try:
    from ase.thermochemistry import IdealGasThermo, HarmonicThermo
    from ase.build import molecule
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False
    print("Warning: ASE not available. Some functionality may be limited.")

# Import our own frequency database
from .frequency_db import FrequencyDB

# Initialize frequency database
freq_db = FrequencyDB()

# Physical constants
CM1_TO_EV = 1.239841984e-4  # cm^-1 to eV conversion
EV_TO_CM1 = 1 / CM1_TO_EV   # eV to cm^-1 conversion

def cm1_to_ev(frequency_cm1: float) -> float:
    """Convert frequency from cm^-1 to eV"""
    return frequency_cm1 * CM1_TO_EV

def clean_species_name(species: str) -> tuple:
    """
    Clean species name and determine status
    
    Args:
        species: Species name (e.g., 'CO_g', 'CO*')
        
    Returns:
        tuple: (cleaned_name, status)
    """
    if species.endswith('_g'):
        return species.replace('_g', ''), 'gas'
    elif species.endswith('*'):
        return species.replace('*', ''), 'ads'
    else:
        return species, 'unknown'

# Default ideal gas parameters
IDEAL_GAS_PARAMS = {
    'H2': [2, 'linear', 0],
    'CO': [1, 'linear', 0],
    'CO2': [2, 'linear', 0],
    'H2O': [2, 'nonlinear', 0],
    'CH4': [12, 'nonlinear', 0],
    'C2H4': [4, 'nonlinear', 0],
    'CH3OH': [1, 'nonlinear', 0],
    'HCOOH': [1, 'nonlinear', 0],
    'C2H6': [6, 'nonlinear', 0],
    'C3H8': [3, 'nonlinear', 0],
    'HCOO': [1, 'nonlinear', 0],
}

# Default formation energies for gas species (in eV)
E_FORM_GAS_DICT = {
    'H2': 0.0,
    'CO': -1.23,
    'CO2': -4.46,
    'H2O': -2.46,
    'CH4': -0.68,
    'C2H4': 1.33,
    'CH3OH': -1.68,
    'HCOOH': -2.88,
    'C2H6': -0.28,
    'C3H8': 0.20,
    'HCOO': -2.88,
}

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



def calculate_thermo_correction(species: str, temperature: float = 298.15, 
    fugacity: float = None, frequencies: List[float] = None) -> Dict[str, float]:
    """
    Calculate thermodynamic correction for a given species

    Parameters:
        species: str
            Name of species (e.g. CO_g, CO*, COOH*, ...)
        temperature: float
            Temperature in K
        fugacity: float
            Fugacity in Pa
        frequencies: List[float]
            Frequencies in cm-1
        voltage: float
            Voltage in V
    Returns:
        Dict[str, float]
            Dictionary containing thermodynamic corrections:
            {'ZPE', 'Cp', 'H', 'S', 'TS', 'G'}
    """
    # Convert frequencies from cm-1 to eV
    
    if not ASE_AVAILABLE:
        raise ImportError("ASE is required for thermodynamic calculations")
    
    if frequencies is None:
        if species in ['ele_g', 'H_g']:
            frequencies = []
        else:
            frequencies = freq_db.get_frequencies(species)
        #print(f"frequencies: {frequencies}")
    
    if fugacity is None:
        if species in DEFAULT_FUGACITY_DICT:
            fugacity = DEFAULT_FUGACITY_DICT[species]
        else:
            fugacity = 1e5

    vib_energies = [cm1_to_ev(f) for f in frequencies]

    species, status = clean_species_name(species)

    # Gas
    if status == 'gas':
        if species in ['ele']:
            return {'ZPE': 0, 'Cp': 0, 'H': 0, 'S': 0, 'TS': 0, 'G': 0}

        # get gas parameters
        gas_params = IDEAL_GAS_PARAMS.get(species, [1, 'linear', 0])
        symmetry_number, geometry, spin = gas_params[:3]
        atoms = None
        try:
            atoms = molecule(species)
        except:
            print(f"Warning: Could not create molecule for {species}")

        thermo = IdealGasThermo(
        vib_energies=vib_energies,
        geometry=geometry,
        atoms=atoms,
        symmetrynumber=symmetry_number,
        spin=spin
        )

        ZPE = thermo.get_ZPE_correction()
        H = thermo.get_enthalpy(temperature, verbose=False)
        Cp = H - ZPE
        S = thermo.get_entropy(temperature, fugacity, verbose=False)
        TS = temperature * S
        F = H - TS
    else:  # status == 'ads'
        # Calculate surface phase corrections
        thermo = HarmonicThermo(vib_energies)
        
        ZPE = thermo.get_ZPE_correction()
        H = thermo.get_internal_energy(temperature=temperature, verbose=False)
        Cp = H - ZPE
        S = thermo.get_entropy(temperature=temperature, verbose=False)
        TS = temperature * S
        F = thermo.get_helmholtz_energy(temperature=temperature, verbose=False)
    
    return {
        'ZPE': ZPE,
        'Cp': Cp,
        'H': H,
        'S': S,
        'TS': TS,
        'F': F
    }




def calculate_gibbs_free_energy(species: str, E_form: float = None, temperature: float = 298.15, 
    fugacity: float = None, frequencies: List[float] = None, voltage: float = 0.0,
    verbose: bool = False, facet: str = None, surface_name: str = None):
    """
    Calculate Gibbs free energy for a given species
    """
    if voltage is None:
        voltage = 0.0

    default_manual_solv_dict = {
    'COH*':-0.3,
    'COOH*':-0.3,
    #'CO*':+0.0,
    'CO2_g':+0.4, #RPBE
    'H2_g':+0.09, #RPBE
    'CO_g':-0.18, #RPBE
    'H2O_g':-0.21, #RPBE
    }
    
    species_after_clean, status = clean_species_name(species)


    if E_form is None:
        if status == 'gas':
            E_form = E_FORM_GAS_DICT[species]
        else:
            if facet is None:
                facet = '100'
            if surface_name is None:
                surface_name = 'Cu'
            E_form = DB_instance.get_E_form_energy(species, facet=facet, surface_name=surface_name).formation_energy

    if 'H_g' in species:
        thermo_correction = calculate_thermo_correction(species='H2_g', temperature=temperature, fugacity=1e5, frequencies=frequencies)
        thermo_correction['F'] = thermo_correction['F']/2
        solvation_energy = 0
    elif 'ele_g' in species:
        thermo_correction = calculate_thermo_correction(species='ele_g', temperature=temperature, fugacity=1e5, frequencies=frequencies)
        thermo_correction['F'] = 0 - voltage
        solvation_energy = 0
    else:
        #print(f"species: {species}, temperature: {temperature}, fugacity: {fugacity}, frequencies: {frequencies}")
        thermo_correction = calculate_thermo_correction(species=species, temperature=temperature, fugacity=fugacity, frequencies=frequencies)
        if species in default_manual_solv_dict:
            solvation_energy = default_manual_solv_dict[species]
        else:
            solvation_energy = DB_instance.get_solvation_energy(species)
    
    if verbose:
        if status == 'ads':
            print(f"species: {species} on {surface_name} ({facet})")
        else:
            print(f"species: {species}")
        print(f"E_form: {E_form:.2f}, thermo_correction: {thermo_correction['F']:.2f}, solvation_energy: {solvation_energy:.2f}, G: {E_form + thermo_correction['F'] + solvation_energy:.2f}")
    return E_form + thermo_correction['F'] + solvation_energy

def calculate_species_E_form_from_scaling_relation(species: str, dE_CO_ads: float = 0.0, 
facet: str = '100', scaling_relation: dict = None):
    """
    Calculate species E_form from scaling relation
    """
    species_after_clean, status = clean_species_name(species)
    if status == 'gas':
        if species in ['H_g', 'ele_g']:
            return 0.0
        else:
            return E_FORM_GAS_DICT[species]
        #E_form = DB_instance.get_E_form_energy(species, facet=None, surface_name=None).formation_energy
    if scaling_relation is None:
        scaling_relation = DATASET_DICT[facet]['scaling_relation']
    if species not in scaling_relation:
        print(f"species {species} not in scaling relation")
        return None
    return scaling_relation[species][0] * dE_CO_ads + scaling_relation[species][1]


def calculate_free_energy_one_state_from_scaling_relation(state_species_list: List[str], dE_CO_ads: float = 0.0, 
    temperature: float = 298.15, fugacity: float = None, frequencies=None, voltage: float = 0.0, verbose: bool = False):
    """
    Calculate free energy for a given species
    """
    G_state = 0
    for species in state_species_list:
        #print(f"species: {species}")
        if species in ['*']:
            continue
        E_form_species = calculate_species_E_form_from_scaling_relation(species, dE_CO_ads, '100')
        G_state += calculate_gibbs_free_energy(species=species, E_form=E_form_species, temperature=temperature, 
        fugacity=fugacity, frequencies=frequencies, voltage=voltage, verbose=verbose)
        if verbose:
            print(f"species: {species}, E_form: {E_form_species:.2f} -> G: {G_state:.2f}")
    return G_state


if __name__ == "__main__":
    #print(calculate_species_E_form_from_scaling_relation('H*', 0.0, '100'))
    #print(calculate_free_energy_one_state_from_scaling_relation(['H_g', 'ele_g', 'COH*'], 0.0, 298.15, verbose=True))
    a=calculate_thermo_correction('H2_g', temperature=400, fugacity=1e5)
    print(a)