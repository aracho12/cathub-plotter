"""
Input File Parser for species data
"""
import pandas as pd
import ast
import os
from typing import Dict, Any

class InputFileParser:
    """Parser for input files containing species data"""
    
    @staticmethod
    def parse_file(filename: str) -> Dict[str, Dict[str, Any]]:
        """
        Reads input file (text or Excel) and converts it into a dictionary keyed by species.

        Args:
            filename: The path to input file (.txt, .xlsx, .xls, .csv).

        Returns:
            dict: {species_name: {status, formation_energy, frequencies, ...}}
        """
        # Determine file type and read accordingly
        _, ext = os.path.splitext(filename)
        
        if ext.lower() in ['.xlsx', '.xls']:
            # Read Excel file
            df = pd.read_excel(filename)
        elif ext.lower() == '.csv':
            # Read CSV file
            df = pd.read_csv(filename)
        else:
            # Read tab-separated text file
            df = pd.read_csv(filename, sep='\t')

        species_data = {}

        for _, row in df.iterrows():
            status = row['status']
            species_name = row['species_name']

            # Determine full species name
            if status == 'gas':
                # Add _g for gas species
                full_name = f"{species_name}_g"
            elif status == 'ads':
                # Add * for adsorbed species
                full_name = f"{species_name}*"
            elif status == 'ts':
                # Add * for transition states
                full_name = f"{species_name}*"
            elif status == 'slab':
                # Keep slab species as is
                full_name = species_name
            else:
                full_name = species_name

            # Parse frequencies (string to list)
            freq_str = row['frequencies']
            try:
                frequencies = ast.literal_eval(freq_str)
                if not isinstance(frequencies, list):
                    frequencies = []
            except Exception:
                frequencies = []

            # Parse formation_energy
            try:
                formation_energy = float(row['formation_energy'])
            except Exception:
                formation_energy = None

            # Build species data dict
            species_data[full_name] = {
                'status': status,
                'species_name': species_name,
                'surface_name': row['surface_name'],
                'site_name': row['site_name'],
                'formation_energy': formation_energy,
                'frequencies': frequencies,
                'reference': row['reference'],
                'solvation_energy': row.get('solvation_energy', None)
            }

        return species_data


# Backward compatibility function
def parse_input_file(filename: str) -> Dict[str, Dict[str, Any]]:
    """Backward compatibility function"""
    return InputFileParser.parse_file(filename)


if __name__ == "__main__":
    # Test parsing
    species_data = parse_input_file('input.txt')

    # Print a few gas species
    print("=== GAS SPECIES ===")
    for name in ['CO_g', 'H2_g', 'H2O_g']:
        if name in species_data:
            print(f"\n{name}:")
            print(f"  Formation energy: {species_data[name]['formation_energy']}")
            print(f"  Frequencies: {species_data[name]['frequencies'][:3]}...")  # Show first 3

    # Print a few adsorbed species
    print("\n=== ADSORBED SPECIES ===")
    for name in ['CO*', 'H*', 'COH*']:
        if name in species_data:
            print(f"\n{name}:")
            print(f"  Formation energy: {species_data[name]['formation_energy']}")
            print(f"  Surface: {species_data[name]['surface_name']}")
            print(f"  Site: {species_data[name]['site_name']}")
            print(f"  Frequencies: {species_data[name]['frequencies'][:3]}...")

    print(f"\nTotal {len(species_data)} species loaded")