"""
Free Energy Calculator for Catalysis-Hub data
"""
from ..parsers.cathub import CatalysisHubParser
from ..utils.thermodynamics import calculate_thermo_correction
from ..utils.frequency_db import FrequencyDB
from typing import List, Dict, Optional
import pandas as pd

class FreeEnergyCalculator:
    """
    Calculate Gibbs free energies for reactions from Catalysis-Hub
    """
    
    # Default manual solvation energy corrections (in eV)
    DEFAULT_SOLVATION_DICT = {
        'COH*': -0.25,
        'COOH*': -0.25,
        'CO2_g': +0.4,  
        'H2_g': +0.09,  
    }
    
    def __init__(self, temperature=298.15, voltage=0.0, fugacity_dict=None, 
                 use_solvation=True, solvation_dict=None):
        """
        Initialize calculator with thermodynamic conditions
        
        Args:
            temperature: Temperature in Kelvin
            voltage: Electrochemical potential in Volts
            fugacity_dict: Dictionary of gas fugacities {species: fugacity_in_Pa}
            use_solvation: Whether to apply solvation energy corrections (default: True)
            solvation_dict: Custom solvation energies {species: energy_in_eV}
                           If None, uses DEFAULT_SOLVATION_DICT + database values
        """
        self.parser = CatalysisHubParser()
        self.temperature = temperature
        self.voltage = voltage
        self.fugacity_dict = fugacity_dict or {}
        self.use_solvation = use_solvation
        self.solvation_dict = solvation_dict or self.DEFAULT_SOLVATION_DICT.copy()
        
        # Initialize frequency database
        self.freq_db = FrequencyDB()
    
    def search_and_calculate(self,
                            reactants: List[str] = None,
                            products: List[str] = None,
                            surface: str = None,
                            facet: str = None,
                            dft_code: str = None,
                            dft_functional: str = None,
                            limit: int = 50,
                            sort_by: str = 'composition',
                            ascending: bool = True) -> pd.DataFrame:
        """
        Search Catalysis-Hub and calculate free energies
        
        Args:
            reactants: List of reactant species (e.g., ['CO_g', '*'])
            products: List of product species (e.g., ['CO*'])
            surface: Surface composition filter (e.g., 'Cu', 'Pt')
            facet: Facet filter (e.g., '111', '100')
            dft_code: DFT code filter with partial matching
            dft_functional: DFT functional filter with partial matching
            limit: Maximum number of results
            sort_by: Column to sort by ('composition', 'delta_G', 'delta_E_DFT', 'delta_F_thermo')
            ascending: Sort order (True for ascending, False for descending)
        
        Returns:
            DataFrame with reaction data and calculated free energies
        """
        # Search reactions from Catalysis-Hub
        reactions = self.parser.search_reactions(
            reactants=reactants,
            products=products,
            chemical_composition=surface,
            facet=facet,
            dft_code=dft_code,
            dft_functional=dft_functional,
            limit=limit
        )
        
        if not reactions:
            return pd.DataFrame()
        
        # Calculate free energies for each reaction
        results = []
        for reaction in reactions:
            converted = self.parser.convert_reaction_to_your_format(reaction)
            
            try:
                result = self.calculate_reaction_free_energy(converted)
                results.append(result)
            except Exception as e:
                print(f"Error calculating reaction {converted['id']}: {e}")
                continue
        
        df = pd.DataFrame(results)
        
        # Sort results
        if not df.empty and sort_by in df.columns:
            # Custom sorting for composition to prioritize exact matches
            if sort_by == 'composition' and surface:
                # Put exact matches first, then alphabetical
                df['_sort_priority'] = df['composition'].apply(
                    lambda x: 0 if x.startswith(surface) else 1
                )
                df = df.sort_values(
                    by=['_sort_priority', 'composition', 'delta_G'], 
                    ascending=[True, ascending, True]
                ).drop('_sort_priority', axis=1).reset_index(drop=True)
            else:
                df = df.sort_values(by=sort_by, ascending=ascending).reset_index(drop=True)
        
        return df
    
    def get_solvation_energy(self, species: str) -> float:
        """
        Get solvation energy for a species
        
        Priority:
        1. Custom solvation_dict provided at initialization
        2. Database lookup
        3. Return 0 if not found
        
        Args:
            species: Species name
            
        Returns:
            Solvation energy in eV
        """
        if not self.use_solvation:
            return 0.0
        
        # Check custom solvation dict first
        if species in self.solvation_dict:
            return self.solvation_dict[species]
        
        # Try frequency database lookup
        try:
            return self.freq_db.get_solvation_energy(species)
        except:
            return 0.0
    
    def calculate_reaction_free_energy(self, reaction: Dict) -> Dict:
        """
        Calculate free energy for a single reaction
        
        ΔG_reaction = ΔE_DFT + ΔF_thermo + ΔE_solv
        where:
            ΔE_DFT: DFT energy from Catalysis-Hub
            ΔF_thermo: Thermodynamic correction (ZPE - TS)
            ΔE_solv: Solvation energy correction
        
        Args:
            reaction: Reaction dictionary from parser
        
        Returns:
            Dictionary with calculated energies
        """
        # Calculate thermodynamic and solvation corrections for reactants
        F_reactants = 0
        E_solv_reactants = 0
        F_reactants_detail = {}
        E_solv_reactants_detail = {}
        
        for species, coeff in reaction['reactants_dict'].items():
            if species == '*':  # Skip empty sites
                continue
            
            try:
                # Use custom fugacity if provided
                fugacity = self.fugacity_dict.get(species, None)
                
                # Thermodynamic correction
                thermo = calculate_thermo_correction(
                    species=species,
                    temperature=self.temperature,
                    fugacity=fugacity,
                    frequencies=None
                )
                # Multiply by stoichiometric coefficient
                F_reactants += thermo['F'] * coeff
                F_reactants_detail[f"{coeff} {species}"] = thermo['F'] * coeff
                
                # Solvation energy correction
                solv_energy = self.get_solvation_energy(species)
                E_solv_reactants += solv_energy * coeff
                E_solv_reactants_detail[f"{coeff} {species}"] = solv_energy * coeff
                
            except Exception as e:
                print(f"Warning: Could not calculate thermo for {species}: {e}")
                F_reactants_detail[f"{coeff} {species}"] = 0.0
                E_solv_reactants_detail[f"{coeff} {species}"] = 0.0
        
        # Calculate thermodynamic and solvation corrections for products
        F_products = 0
        E_solv_products = 0
        F_products_detail = {}
        E_solv_products_detail = {}
        
        for species, coeff in reaction['products_dict'].items():
            if species == '*':
                continue
            
            try:
                fugacity = self.fugacity_dict.get(species, None)
                
                # Thermodynamic correction
                thermo = calculate_thermo_correction(
                    species=species,
                    temperature=self.temperature,
                    fugacity=fugacity,
                    frequencies=None
                )
                F_products += thermo['F'] * coeff
                F_products_detail[f"{coeff} {species}"] = thermo['F'] * coeff
                
                # Solvation energy correction
                solv_energy = self.get_solvation_energy(species)
                E_solv_products += solv_energy * coeff
                E_solv_products_detail[f"{coeff} {species}"] = solv_energy * coeff
                
            except Exception as e:
                print(f"Warning: Could not calculate thermo for {species}: {e}")
                F_products_detail[f"{coeff} {species}"] = 0.0
                E_solv_products_detail[f"{coeff} {species}"] = 0.0
        
        # Calculate energy changes
        delta_F_thermo = F_products - F_reactants
        delta_E_solv = E_solv_products - E_solv_reactants
        delta_E_DFT = reaction['reaction_energy']
        delta_G = delta_E_DFT + delta_F_thermo + delta_E_solv
        
        return {
            'id': reaction['id'],
            'equation': reaction['equation'],
            'surface': reaction['surface'],
            'composition': reaction['composition'],
            'facet': reaction['facet'],
            'sites_str': reaction['sites_str'],
            'dft_code': reaction.get('dft_code', 'Unknown'),
            'dft_functional': reaction.get('dft_functional', 'Unknown'),
            'pub_id': reaction.get('pub_id', 'Unknown'),
            'pub_authors': reaction.get('pub_authors', 'Unknown'),
            'pub_year': reaction.get('pub_year', 'Unknown'),
            'delta_E_DFT': delta_E_DFT,
            'delta_F_thermo': delta_F_thermo,
            'delta_E_solv': delta_E_solv,
            'delta_G': delta_G,
            'F_reactants': F_reactants,
            'F_products': F_products,
            'E_solv_reactants': E_solv_reactants,
            'E_solv_products': E_solv_products,
            'F_reactants_detail': F_reactants_detail,
            'F_products_detail': F_products_detail,
            'E_solv_reactants_detail': E_solv_reactants_detail,
            'E_solv_products_detail': E_solv_products_detail,
            'activation_energy': reaction.get('activation_energy'),
        }