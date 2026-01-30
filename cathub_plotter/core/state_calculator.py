"""
State Energy Calculator
Read formation energies and frequencies from input file and calculate free energies for states
"""
from typing import List, Dict, Optional
import pandas as pd
import os
from ..utils.thermodynamics import calculate_thermo_correction

class InputDataParser:
    """
    Parse input file to extract species data
    """
    
    def __init__(self, input_file: str = 'input.txt'):
        """
        Initialize parser with input file path
        
        Args:
            input_file: Path to input file (.txt, .xlsx, .xls, .csv)
        """
        self.input_file = input_file
        self.data = self._parse_input_file()
    
    def _parse_input_file(self) -> pd.DataFrame:
        """
        Parse input file into pandas DataFrame
        
        Returns:
            DataFrame with species data
        """
        # Determine file type and read accordingly
        _, ext = os.path.splitext(self.input_file)
        
        if ext.lower() in ['.xlsx', '.xls']:
            # Read Excel file
            df = pd.read_excel(self.input_file)
        elif ext.lower() == '.csv':
            # Read CSV file
            df = pd.read_csv(self.input_file)
        else:
            # Read tab-separated text file
            df = pd.read_csv(self.input_file, sep='\t')
        
        return df
    
    def get_species_data(self, species: str) -> Optional[Dict]:
        """
        Get data for a specific species
        
        Args:
            species: Species name (e.g., 'CO_g', 'H*', 'ele_g', 'H_g', 'H-ele*')
        
        Returns:
            Dictionary with species data or None if not found
        """
        # Convert species name to match input.txt format
        # 'CO_g' -> species_name='CO', status='gas'
        # 'H*' -> species_name='H', status='ads'
        # 'ele_g' -> species_name='ele', status='gas'
        # 'H_g' -> species_name='H', status='gas'
        # 'H-ele*' -> species_name='H-ele', status='ts'
        
        if species == '*':
            # Empty site - no energy contribution
            return None
        
        if species.endswith('_g'):
            # Gas phase
            species_name = species[:-2]  # Remove '_g'
            status = 'gas'
        elif species.endswith('*'):
            # Could be adsorbed or transition state
            species_name = species[:-1]  # Remove '*'
            
            # First try to find as TS (for species like 'H-ele*', 'COH-H-ele*')
            mask_ts = (self.data['species_name'] == species_name) & (self.data['status'] == 'ts')
            matches_ts = self.data[mask_ts]
            
            if len(matches_ts) > 0:
                status = 'ts'
            else:
                # If not found as TS, treat as adsorbate
                status = 'ads'
        else:
            # Default to gas phase
            species_name = species
            status = 'gas'
        
        # Search in dataframe
        mask = (self.data['species_name'] == species_name) & (self.data['status'] == status)
        matches = self.data[mask]
        
        if len(matches) == 0:
            print(f"Warning: Species '{species}' (name={species_name}, status={status}) not found in input file")
            return None
        
        # Use first match
        row = matches.iloc[0]
        
        # Parse frequencies from string to list
        freq_str = row['frequencies']
        if pd.isna(freq_str) or freq_str == '[]':
            frequencies = []
        else:
            # Remove brackets and split by comma
            freq_str = freq_str.strip('[]')
            if freq_str:
                frequencies = [float(f.strip()) for f in freq_str.split(',')]
            else:
                frequencies = []
        
        return {
            'species': species,
            'species_name': row['species_name'],
            'status': row['status'],
            'surface_name': row['surface_name'],
            'site_name': row['site_name'],
            'formation_energy': row['formation_energy'],
            'frequencies': frequencies,
            'solvation_energy': row.get('solvation_energy', 0.0)
        }


class StateEnergyCalculator:
    """
    Calculate free energy for states (collections of species)
    """
    
    def __init__(self, input_file: str = 'input.txt', 
                 temperature: float = 298.15, voltage: float = 0.0,
                 fugacity_dict: Dict = None,
                 use_gibbs_energy: bool = False):
        """
        Initialize calculator
        
        Args:
            input_file: Path to input.txt file
            temperature: Temperature in K
            voltage: Voltage in V (for electrochemical reactions)
            fugacity_dict: Dictionary of fugacities for gas species {species: Pa}
            use_gibbs_energy: If True, formation_energy is treated as Gibbs free energy at U=0V
                            (no thermodynamic corrections applied, only voltage correction)
        """
        self.parser = InputDataParser(input_file)
        self.temperature = temperature
        self.voltage = voltage
        self.fugacity_dict = fugacity_dict or {}
        self.use_gibbs_energy = use_gibbs_energy
        
        # Default fugacity values
        self.default_fugacity_dict = {
            'H2_g': 101325,
            'CO2_g': 101325,
            'CO_g': 5562,
            'HCOOH_g': 2,
            'CH3OH_g': 6079,
            'H2O_g': 3534,
            'CH4_g': 20467,
        }
    
    def calculate_species_free_energy(self, species: str, verbose: bool = False, 
                                      is_transition_state: bool = False,
                                      beta: float = 0.5,
                                      reactants: List[str] = None) -> Dict:
        """
        Calculate free energy for a single species
        
        If use_gibbs_energy=False (default):
            G = E_form + F_thermo + E_solv
        
        If use_gibbs_energy=True:
            G = G_0V + voltage_correction
            (formation_energy is treated as Gibbs energy at U=0V, no thermo corrections)
        
        For transition states in electrochemical reactions:
        G_TS includes voltage-dependent correction based on beta
        
        Args:
            species: Species name (e.g., 'CO_g', 'H*', 'ele_g')
            verbose: Print detailed information
            is_transition_state: Whether this species is a TS
            beta: Symmetry factor for electrochemical TS (default 0.5)
            reactants: List of reactant species (for TS correction)
        
        Returns:
            Dictionary with energy components
        """
        # Special case: empty site
        if species == '*':
            return {
                'species': '*',
                'E_form': 0.0,
                'F_thermo': 0.0,
                'E_solv': 0.0,
                'G': 0.0
            }
        
        # Get species data from input file
        species_data = self.parser.get_species_data(species)
        if species_data is None:
            raise ValueError(f"Could not find species '{species}' in input file")
        
        E_form = species_data['formation_energy']
        frequencies = species_data['frequencies']
        
        # If using Gibbs energy mode, skip thermodynamic corrections
        if self.use_gibbs_energy:
            # formation_energy is already Gibbs free energy at U=0V
            # Only apply voltage correction for electrons
            if species == 'ele_g':
                # Electron with electrochemical potential
                F_thermo = -self.voltage
                E_solv = 0.0
            elif species == 'H_g':
                # H_g: half of H2_g Gibbs energy + voltage correction
                # Since H_g comes with ele_g in reactions, voltage is already in ele_g
                F_thermo = 0.0
                E_solv = 0.0
            else:
                # No thermodynamic correction needed
                F_thermo = 0.0
                E_solv = 0.0
            
            G = E_form + F_thermo + E_solv
        else:
            # Original mode: calculate thermodynamic corrections
            # Get fugacity
            fugacity = self.fugacity_dict.get(species, None)
            if fugacity is None:
                fugacity = self.default_fugacity_dict.get(species, 1e5)
            
            # Special cases for H_g and ele_g
            if species == 'H_g':
                # H_g is half of H2_g
                # Important: Use H2_g frequencies, not H_g frequencies!
                h2_species_data = self.parser.get_species_data('H2_g')
                if h2_species_data is not None:
                    h2_frequencies = h2_species_data['frequencies']
                else:
                    h2_frequencies = []
                
                thermo = calculate_thermo_correction(
                    species='H2_g',
                    temperature=self.temperature,
                    fugacity=fugacity,
                    frequencies=h2_frequencies
                )
                F_thermo = thermo['F'] / 2.0
                E_solv = 0.0
            elif species == 'ele_g':
                # Electron with electrochemical potential
                thermo = calculate_thermo_correction(
                    species='ele_g',
                    temperature=self.temperature,
                    fugacity=1e5,
                    frequencies=[]
                )
                F_thermo = -self.voltage
                E_solv = 0.0
            else:
                # Normal species
                thermo = calculate_thermo_correction(
                    species=species,
                    temperature=self.temperature,
                    fugacity=fugacity,
                    frequencies=frequencies
                )
                F_thermo = thermo['F']
                
                # Solvation energy from input file
                E_solv_raw = species_data.get('solvation_energy', 0.0)
                if pd.isna(E_solv_raw) or E_solv_raw == '':
                    E_solv = 0.0
                else:
                    # Parse solvation energy (might be "0.00 eV" format)
                    if isinstance(E_solv_raw, str):
                        E_solv = float(E_solv_raw.split()[0])
                    else:
                        E_solv = float(E_solv_raw)
            
            G = E_form + F_thermo + E_solv
        
        # Apply voltage-dependent correction for electrochemical transition states
        # Following CatMAP's simple_electrochemical mode:
        # thermo_dict[TS] = -voltage + beta * (voltage - voltage_ref)
        #                 = -voltage + beta * voltage  (when voltage_ref = 0)
        #                 = -(1 - beta) * voltage
        # This accounts for the partial charge transfer in the TS
        voltage_correction = 0.0
        if is_transition_state and reactants:
            # Count electrons in reactants (IS)
            n_ele = sum(1 for r in reactants if r in ['ele_g', 'pe_g'])
            if n_ele > 0:
                # Apply electrochemical TS correction
                # TS gets -(1-beta)*voltage correction
                # IS gets -voltage correction from ele_g
                # So relative: Ga changes by beta*voltage per electron
                voltage_correction = -(1 - beta) * self.voltage
                G += voltage_correction
                
                if verbose:
                    print(f"  [Electrochemical TS correction]")
                    print(f"  n_electrons: {n_ele}")
                    print(f"  beta: {beta}")
                    print(f"  voltage_correction: {voltage_correction:.4f} eV")
                    print(f"  Effect on Ga: {-beta * self.voltage:.4f} eV per electron")
        
        if verbose:
            print(f"Species: {species}")
            if self.use_gibbs_energy:
                print(f"  G(U=0V): {E_form:.4f} eV")
                if F_thermo != 0:
                    print(f"  V_corr: {F_thermo:.4f} eV")
            else:
                print(f"  E_form: {E_form:.4f} eV")
                print(f"  F_thermo: {F_thermo:.4f} eV")
                print(f"  E_solv: {E_solv:.4f} eV")
            if voltage_correction != 0:
                print(f"  V_corr (TS): {voltage_correction:.4f} eV")
            print(f"  G: {G:.4f} eV")
        
        return {
            'species': species,
            'E_form': E_form,
            'F_thermo': F_thermo,
            'E_solv': E_solv,
            'voltage_correction': voltage_correction,
            'G': G,
            'frequencies': frequencies,
            'temperature': self.temperature,
            'voltage': self.voltage,
            'beta': beta if is_transition_state else None
        }
    
    def calculate_state_free_energy(self, state: List[str], verbose: bool = False,
                                   is_transition_state: bool = False,
                                   beta: float = 0.5,
                                   reactants: List[str] = None) -> Dict:
        """
        Calculate total free energy for a state
        
        A state is a collection of species, e.g.:
        - ['CO_g', '*'] : CO in gas phase + empty site
        - ['H_g', 'ele_g', 'H*'] : H + electron + H adsorbed
        - ['H-ele*'] : Transition state (with is_transition_state=True)
        
        Args:
            state: List of species in the state
            verbose: Print detailed information
            is_transition_state: Whether this state is a TS (single species)
            beta: Symmetry factor for electrochemical TS (default 0.5)
            reactants: Reactant species list (for TS correction)
        
        Returns:
            Dictionary with state energy information
        """
        if verbose:
            print(f"\n{'='*60}")
            print(f"Calculating free energy for state: {state}")
            print(f"Temperature: {self.temperature} K, Voltage: {self.voltage} V")
            print(f"{'='*60}")
        
        species_energies = []
        G_total = 0.0
        E_form_total = 0.0
        F_thermo_total = 0.0
        E_solv_total = 0.0
        voltage_correction_total = 0.0
        
        for species in state:
            # For TS, pass the parameters
            if is_transition_state and len(state) == 1:
                species_result = self.calculate_species_free_energy(
                    species, verbose=verbose, 
                    is_transition_state=True,
                    beta=beta,
                    reactants=reactants
                )
            else:
                species_result = self.calculate_species_free_energy(species, verbose=verbose)
            
            species_energies.append(species_result)
            
            G_total += species_result['G']
            E_form_total += species_result['E_form']
            F_thermo_total += species_result['F_thermo']
            E_solv_total += species_result['E_solv']
            voltage_correction_total += species_result.get('voltage_correction', 0.0)
        
        if verbose:
            print(f"\n{'='*60}")
            print(f"Total State Energy:")
            print(f"  E_form_total: {E_form_total:.4f} eV")
            print(f"  F_thermo_total: {F_thermo_total:.4f} eV")
            print(f"  E_solv_total: {E_solv_total:.4f} eV")
            if voltage_correction_total != 0:
                print(f"  V_corr_total (TS): {voltage_correction_total:.4f} eV")
            print(f"  G_total: {G_total:.4f} eV")
            print(f"{'='*60}\n")
        
        return {
            'state': state,
            'species_energies': species_energies,
            'G_total': G_total,
            'E_form_total': E_form_total,
            'F_thermo_total': F_thermo_total,
            'E_solv_total': E_solv_total,
            'voltage_correction_total': voltage_correction_total,
            'temperature': self.temperature,
            'voltage': self.voltage,
            'beta': beta if is_transition_state else None
        }
    
    def calculate_reaction_free_energy(self, reactants: List[str], products: List[str], 
                                      verbose: bool = False) -> Dict:
        """
        Calculate free energy change for a reaction
        
        ΔG = G(products) - G(reactants)
        
        Args:
            reactants: List of reactant species
            products: List of product species
            verbose: Print detailed information
        
        Returns:
            Dictionary with reaction energy information
        """
        if verbose:
            print(f"\n{'#'*60}")
            print(f"REACTION: {' + '.join(reactants)} -> {' + '.join(products)}")
            print(f"{'#'*60}")
        
        reactants_result = self.calculate_state_free_energy(reactants, verbose=verbose)
        products_result = self.calculate_state_free_energy(products, verbose=verbose)
        
        delta_G = products_result['G_total'] - reactants_result['G_total']
        delta_E_form = products_result['E_form_total'] - reactants_result['E_form_total']
        delta_F_thermo = products_result['F_thermo_total'] - reactants_result['F_thermo_total']
        delta_E_solv = products_result['E_solv_total'] - reactants_result['E_solv_total']
        
        if verbose:
            print(f"\n{'#'*60}")
            print(f"REACTION ENERGY CHANGE:")
            print(f"  ΔE_form: {delta_E_form:.4f} eV")
            print(f"  ΔF_thermo: {delta_F_thermo:.4f} eV")
            print(f"  ΔE_solv: {delta_E_solv:.4f} eV")
            print(f"  ΔG: {delta_G:.4f} eV")
            print(f"{'#'*60}\n")
        
        return {
            'reactants': reactants,
            'products': products,
            'reactants_result': reactants_result,
            'products_result': products_result,
            'delta_G': delta_G,
            'delta_E_form': delta_E_form,
            'delta_F_thermo': delta_F_thermo,
            'delta_E_solv': delta_E_solv,
            'temperature': self.temperature,
            'voltage': self.voltage
        }


if __name__ == "__main__":
    # Example usage
    print("="*60)
    print("State Energy Calculator - Example Usage")
    print("="*60)
    
    # Initialize calculator
    calculator = StateEnergyCalculator(
        input_file='input.txt',
        temperature=298.15,
        voltage=0.0
    )
    
    # Example 1: Calculate free energy for state ['CO_g', '*']
    print("\n" + "="*60)
    print("Example 1: State ['CO_g', '*']")
    print("="*60)
    state1 = ['CO_g', '*']
    result1 = calculator.calculate_state_free_energy(state1, verbose=True)
    
    # Example 2: Calculate free energy for state ['H_g', 'ele_g', 'H*']
    print("\n" + "="*60)
    print("Example 2: State ['H_g', 'ele_g', 'H*']")
    print("="*60)
    state2 = ['H_g', 'ele_g', 'H*']
    result2 = calculator.calculate_state_free_energy(state2, verbose=True)
    
    # Example 3: Calculate reaction free energy
    # CO_g + * -> CO*
    print("\n" + "="*60)
    print("Example 3: Reaction CO_g + * -> CO*")
    print("="*60)
    reaction_result = calculator.calculate_reaction_free_energy(
        reactants=['CO_g', '*'],
        products=['CO*'],
        verbose=True
    )
    
    # Example 4: CO2 reduction reaction
    # CO2_g + H_g + ele_g + * -> COOH*
    print("\n" + "="*60)
    print("Example 4: CO2 reduction - CO2_g + H_g + ele_g + * -> COOH*")
    print("="*60)
    reaction_result2 = calculator.calculate_reaction_free_energy(
        reactants=['CO2_g', 'H_g', 'ele_g', '*'],
        products=['COOH*'],
        verbose=True
    )

