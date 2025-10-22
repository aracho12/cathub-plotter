"""
MKM Energy Calculator
Calculate reaction energies from parsed .mkm files using input.txt data
"""
from ..parsers.mkm import MKMFileParser
from .state_calculator import StateEnergyCalculator
from typing import Dict, List, Optional
import pandas as pd


class MKMEnergyCalculator:
    """
    Calculate reaction energies for reactions parsed from .mkm files
    """
    
    def __init__(self, 
                 input_file: str = 'input.txt',
                 temperature: float = 298.15,
                 voltage: float = 0.0,
                 fugacity_dict: Dict = None):
        """
        Initialize calculator
        
        Args:
            input_file: Path to input.txt file with species data
            temperature: Temperature in K
            voltage: Voltage in V
            fugacity_dict: Dictionary of fugacities for gas species
        """
        self.state_calculator = StateEnergyCalculator(
            input_file=input_file,
            temperature=temperature,
            voltage=voltage,
            fugacity_dict=fugacity_dict
        )
        self.temperature = temperature
        self.voltage = voltage
    
    def calculate_reaction_energy(self, 
                                  reactants: List[str],
                                  products: List[str],
                                  ts: Optional[List[str]] = None,
                                  beta: float = 0.5,
                                  verbose: bool = False) -> Dict:
        """
        Calculate reaction energies (ΔG and Ga if TS exists)
        
        For electrochemical reactions (with ele_g or pe_g):
        - IS and FS energies include -voltage correction for each electron
        - TS energy includes voltage-dependent correction: -voltage + beta*voltage
        - This makes Ga voltage-dependent for electrochemical steps
        
        Args:
            reactants: List of reactant species (Initial State)
            products: List of product species (Final State)
            ts: List of transition state species (optional)
            beta: Symmetry factor for electrochemical TS (default 0.5)
            verbose: Print detailed information
        
        Returns:
            Dictionary with:
                - G_IS: Free energy of initial state
                - G_FS: Free energy of final state
                - G_TS: Free energy of transition state (if exists)
                - delta_G: Reaction free energy (G_FS - G_IS)
                - Ga: Activation energy (G_TS - G_IS) (if TS exists)
                - Ga_reverse: Reverse activation energy (G_TS - G_FS) (if TS exists)
                - beta: Symmetry factor used (if TS exists)
        """
        if verbose:
            print(f"\n{'='*70}")
            print(f"Calculating reaction energies")
            print(f"Reactants (IS): {reactants}")
            print(f"Products (FS): {products}")
            if ts:
                print(f"Transition State (TS): {ts}")
            print(f"Temperature: {self.temperature} K, Voltage: {self.voltage} V")
            print(f"{'='*70}")
        
        # Calculate Initial State (IS) energy
        IS_result = self.state_calculator.calculate_state_free_energy(
            reactants, verbose=verbose
        )
        G_IS = IS_result['G_total']
        
        # Calculate Final State (FS) energy
        FS_result = self.state_calculator.calculate_state_free_energy(
            products, verbose=verbose
        )
        G_FS = FS_result['G_total']
        
        # Calculate reaction free energy
        delta_G = G_FS - G_IS
        
        result = {
            'reactants': reactants,
            'products': products,
            'ts': ts,
            'G_IS': G_IS,
            'G_FS': G_FS,
            'delta_G': delta_G,
            'IS_result': IS_result,
            'FS_result': FS_result,
            'temperature': self.temperature,
            'voltage': self.voltage
        }
        
        # Calculate Transition State (TS) energy if provided
        if ts:
            try:
                TS_result = self.state_calculator.calculate_state_free_energy(
                    ts, verbose=verbose,
                    is_transition_state=True,
                    beta=beta,
                    reactants=reactants
                )
                G_TS = TS_result['G_total']
                
                # Activation energy (forward)
                # For electrochemical reactions, Ga is voltage-dependent due to beta
                Ga = G_TS - G_IS
                
                # Reverse activation energy
                Ga_reverse = G_TS - G_FS
                
                result.update({
                    'G_TS': G_TS,
                    'Ga': Ga,
                    'Ga_reverse': Ga_reverse,
                    'TS_result': TS_result,
                    'beta': beta
                })
            except (ValueError, KeyError) as e:
                # TS data not found - continue without barrier
                if verbose:
                    print(f"Warning: TS data not available ({e}). Continuing without activation barrier.")
                result.update({
                    'G_TS': None,
                    'Ga': None,
                    'Ga_reverse': None,
                    'TS_result': None,
                    'beta': beta
                })
            
            if verbose:
                # Check if this is an electrochemical reaction
                n_ele = sum(1 for r in reactants if r in ['ele_g', 'pe_g'])
                
                print(f"\n{'='*70}")
                print(f"ENERGETICS SUMMARY:")
                print(f"  G(IS):  {G_IS:8.4f} eV")
                print(f"  G(TS):  {G_TS:8.4f} eV")
                print(f"  G(FS):  {G_FS:8.4f} eV")
                print(f"  ΔG:     {delta_G:8.4f} eV  (G_FS - G_IS)")
                print(f"  Ga:     {Ga:8.4f} eV  (G_TS - G_IS, forward)")
                print(f"  Ga_rev: {Ga_reverse:8.4f} eV  (G_TS - G_FS, reverse)")
                if n_ele > 0:
                    print(f"  [Electrochemical reaction: {n_ele} electron(s), beta={beta}]")
                    print(f"  Note: Ga is voltage-dependent (V={self.voltage:.3f} V)")
                print(f"{'='*70}\n")
        else:
            if verbose:
                print(f"\n{'='*70}")
                print(f"ENERGETICS SUMMARY:")
                print(f"  G(IS):  {G_IS:8.4f} eV")
                print(f"  G(FS):  {G_FS:8.4f} eV")
                print(f"  ΔG:     {delta_G:8.4f} eV  (G_FS - G_IS)")
                print(f"  (No transition state)")
                print(f"{'='*70}\n")
        
        return result
    
    def calculate_mkm_reactions(self, 
                                mkm_file: str,
                                verbose: bool = False,
                                reaction_numbers: List[int] = None) -> pd.DataFrame:
        """
        Calculate energies for all reactions in a .mkm file
        
        Args:
            mkm_file: Path to .mkm file
            verbose: Print detailed information for each reaction
            reaction_numbers: List of specific reaction numbers to calculate (1-indexed)
                            If None, calculate all reactions
        
        Returns:
            DataFrame with reaction energies
        """
        # Parse mkm file
        mkm_data = MKMFileParser.parse_file(mkm_file)
        parsed_reactions = mkm_data['rxn_expressions']
        
        results = []
        
        for step_num, rxn_data in parsed_reactions.items():
            # Skip if specific reactions requested and this isn't one
            if reaction_numbers and step_num not in reaction_numbers:
                continue
            
            try:
                if verbose:
                    print(f"\n{'#'*70}")
                    print(f"REACTION {step_num}: {rxn_data['original']}")
                    print(f"{'#'*70}")
                
                # Calculate energy
                # Use beta from mkm file if available, otherwise default to 0.5
                beta = rxn_data.get('beta', 0.5)
                if beta is None:
                    beta = 0.5
                
                energy_result = self.calculate_reaction_energy(
                    reactants=rxn_data['reactants'],
                    products=rxn_data['products'],
                    ts=rxn_data['ts'],
                    beta=beta,
                    verbose=verbose
                )
                
                # Compile result
                result_row = {
                    'step_num': step_num,
                    'equation': rxn_data['equation'],
                    'original': rxn_data['original'],
                    'reversible': rxn_data['reversible'],
                    'beta': rxn_data['beta'],
                    'G_IS': energy_result['G_IS'],
                    'G_FS': energy_result['G_FS'],
                    'delta_G': energy_result['delta_G'],
                }
                
                # Add TS info if exists
                if energy_result['ts']:
                    result_row['has_TS'] = True
                    result_row['G_TS'] = energy_result['G_TS']
                    result_row['Ga'] = energy_result['Ga']
                    result_row['Ga_reverse'] = energy_result['Ga_reverse']
                else:
                    result_row['has_TS'] = False
                    result_row['G_TS'] = None
                    result_row['Ga'] = None
                    result_row['Ga_reverse'] = None
                
                results.append(result_row)
                
            except Exception as e:
                print(f"Error calculating reaction {step_num}: {e}")
                import traceback
                traceback.print_exc()
                continue
        
        df = pd.DataFrame(results)
        return df
    
    def export_results(self, df: pd.DataFrame, output_file: str):
        """
        Export results to CSV file
        
        Args:
            df: DataFrame with reaction energies
            output_file: Path to output CSV file
        """
        df.to_csv(output_file, index=False)
        print(f"Results exported to {output_file}")


if __name__ == "__main__":
    print("="*70)
    print("MKM Energy Calculator - Test")
    print("="*70)
    
    # Initialize calculator with specific conditions
    calculator = MKMEnergyCalculator(
        input_file='input.txt',
        temperature=298.15,
        voltage=-0.5  # -0.5 V vs reference
    )
    
    # Test 1: Single reaction calculation
    print("\n" + "="*70)
    print("Test 1: Single reaction with TS")
    print("="*70)
    
    # Volmer reaction: H_g + ele_g + * <-> H-ele* <-> H*
    result1 = calculator.calculate_reaction_energy(
        reactants=['H_g', 'ele_g', '*'],
        products=['H*'],
        ts=['H-ele*'],
        verbose=True
    )
    
    # Test 2: Reaction without TS
    print("\n" + "="*70)
    print("Test 2: CO adsorption without TS")
    print("="*70)
    
    result2 = calculator.calculate_reaction_energy(
        reactants=['CO_g', '*'],
        products=['CO*'],
        ts=None,
        verbose=True
    )
    
    # Test 3: Calculate first 5 reactions from mkm file
    print("\n" + "="*70)
    print("Test 3: Calculate reactions from COR_template.mkm")
    print("="*70)
    
    df = calculator.calculate_mkm_reactions(
        mkm_file='COR_template.mkm',
        verbose=False,  # Set to True for detailed output
        reaction_numbers=[1, 2, 5, 11, 14]  # Calculate specific reactions
    )
    
    print("\n" + "="*70)
    print("Summary Table:")
    print("="*70)
    
    # Print summary table
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', 50)
    
    summary_df = df[['step_num', 'equation', 'delta_G', 'Ga', 'has_TS']]
    print(summary_df.to_string(index=False))
    
    # Export full results
    calculator.export_results(df, 'reaction_energies.csv')
    
    print("\n" + "="*70)
    print("Done!")
    print("="*70)

