"""
Simple frequency database using frequencies.csv
"""
import pandas as pd
import ast
import os
from typing import List, Optional


class FrequencyDB:
    """Simple frequency database using frequencies.csv file"""
    
    def __init__(self, csv_file: str = None):
        """
        Initialize frequency database
        
        Args:
            csv_file: Path to frequencies.csv file
        """
        if csv_file is None:
            # Try to find frequencies.csv in data directory
            current_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
            csv_file = os.path.join(current_dir, 'data', 'frequencies.csv')
        
        self.csv_file = csv_file
        self._load_frequencies()
    
    def _load_frequencies(self):
        """Load frequencies from CSV file"""
        try:
            self.df = pd.read_csv(self.csv_file)
            print(f"Loaded {len(self.df)} frequency entries from {self.csv_file}")
        except FileNotFoundError:
            print(f"Warning: Could not find {self.csv_file}")
            self.df = pd.DataFrame()
    
    def get_frequencies(self, species: str) -> List[float]:
        """
        Get frequencies for a species
        
        Args:
            species: Species name (e.g., 'CO_g', 'CO*')
            
        Returns:
            List of frequencies in cm^-1
        """
        # Clean species name
        clean_species = self._clean_species_name(species)
        
        # Search for exact match first
        exact_match = self.df[self.df['species_name'] == clean_species]
        if not exact_match.empty:
            return self._parse_frequencies(exact_match.iloc[0]['frequencies'])
        
        # Search for partial match (without status suffix)
        if '_g' in clean_species:
            base_name = clean_species.replace('_g', '')
            partial_match = self.df[
                (self.df['species_name'] == base_name) & 
                (self.df['status'] == 'gas')
            ]
            if not partial_match.empty:
                return self._parse_frequencies(partial_match.iloc[0]['frequencies'])
        
        if '*' in clean_species:
            base_name = clean_species.replace('*', '')
            partial_match = self.df[
                (self.df['species_name'] == base_name) & 
                (self.df['status'] == 'ads')
            ]
            if not partial_match.empty:
                return self._parse_frequencies(partial_match.iloc[0]['frequencies'])
        
        # If no match found, return empty list
        print(f"Warning: No frequencies found for species '{species}'")
        return []
    
    def get_solvation_energy(self, species: str) -> float:
        """
        Get solvation energy for a species
        For now, return 0.0 as we don't have solvation data in frequencies.csv
        
        Args:
            species: Species name
            
        Returns:
            Solvation energy in eV (default: 0.0)
        """
        # Default solvation energies for common species
        default_solvation = {
            'COH*': -0.25,
            'COOH*': -0.25,
            'CO2_g': +0.4,
            'H2_g': +0.09,
            'CO_g': -0.18,
            'H2O_g': -0.21,
        }
        
        return default_solvation.get(species, 0.0)
    
    def _clean_species_name(self, species: str) -> str:
        """Clean species name for matching"""
        return species.strip()
    
    def _parse_frequencies(self, freq_str: str) -> List[float]:
        """Parse frequency string to list of floats"""
        try:
            if isinstance(freq_str, str):
                # Remove quotes and parse as list
                freq_str = freq_str.strip('"\'')
                frequencies = ast.literal_eval(freq_str)
                if isinstance(frequencies, list):
                    return [float(f) for f in frequencies]
            return []
        except (ValueError, SyntaxError):
            print(f"Warning: Could not parse frequencies: {freq_str}")
            return []
    
    def list_species(self) -> List[str]:
        """List all available species"""
        return self.df['species_name'].tolist()
    
    def get_species_info(self, species: str) -> dict:
        """Get detailed information about a species"""
        clean_species = self._clean_species_name(species)
        match = self.df[self.df['species_name'] == clean_species]
        
        if not match.empty:
            row = match.iloc[0]
            return {
                'species_name': row['species_name'],
                'status': row['status'],
                'reference': row['reference'],
                'frequencies': self._parse_frequencies(row['frequencies'])
            }
        return {}
