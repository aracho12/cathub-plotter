"""
MKM File Parser for .mkm and .yaml configuration files
"""
import re
import os
import yaml
from typing import Dict, List, Any


class MKMFileParser:
    """Parser for MKM (.mkm) and YAML configuration files"""
    
    @staticmethod
    def parse_file(filename: str) -> Dict[str, Any]:
        """
        Reads a .mkm or .yaml/.yml file and extracts relevant data.

        Args:
            filename: Path to the .mkm or .yaml/.yml file

        Returns:
            dict: {
                'rxn_expressions': parsed reactions,
                'rxn_mechanisms': mechanism dictionary,
                'species_definitions': species definitions,
                'surface_names': list of surfaces,
                ...
            }
        """
        # Check file extension
        _, ext = os.path.splitext(filename)
        
        if ext.lower() in ['.yaml', '.yml']:
            return MKMFileParser._parse_yaml_file(filename)
        else:
            return MKMFileParser._parse_mkm_python_file(filename)

    @staticmethod
    def _parse_yaml_file(filename: str) -> Dict[str, Any]:
        """Parse YAML configuration file"""
        with open(filename, 'r') as f:
            data = yaml.safe_load(f)
        
        rxn_expressions_raw = data.get('rxn_expressions', {})
        rxn_mechanisms_nested = data.get('rxn_mechanisms', {})
        
        # Parse reactions from YAML format
        parsed_reactions = MKMFileParser._parse_rxn_expressions_yaml(rxn_expressions_raw)
        
        # Flatten nested mechanism structure
        flattened_mechanisms = {}
        mechanism_colors = {}
        
        for product, mechanisms in rxn_mechanisms_nested.items():
            if isinstance(mechanisms, dict):
                for mech_name, mech_data in mechanisms.items():
                    if isinstance(mech_data, dict):
                        steps = mech_data.get('steps', [])
                        color = mech_data.get('color', None)
                    else:
                        steps = mech_data
                        color = None
                    
                    # Create combined name: product_mechanism
                    combined_name = f"{product}_{mech_name}"
                    flattened_mechanisms[combined_name] = steps
                    if color:
                        mechanism_colors[combined_name] = color
                    
                    # Also add without product prefix for backward compatibility
                    if mech_name not in flattened_mechanisms:
                        flattened_mechanisms[mech_name] = steps
                        if color:
                            mechanism_colors[mech_name] = color
            else:
                flattened_mechanisms[product] = mechanisms
        
        return {
            'rxn_expressions': parsed_reactions,
            'rxn_expressions_raw': rxn_expressions_raw,
            'rxn_mechanisms': flattened_mechanisms,
            'mechanism_colors': mechanism_colors,
            'species_definitions': data.get('species_definitions', {}),
            'surface_names': data.get('surface_names', []),
            'prefactor_list': data.get('prefactor_list', []),
            'descriptor_names': data.get('descriptor_names', []),
            'descriptor_ranges': data.get('descriptor_ranges', []),
            'products': list(rxn_mechanisms_nested.keys()),
            'mechanisms_by_product': rxn_mechanisms_nested,
        }

    @staticmethod
    def _parse_mkm_python_file(filename: str) -> Dict[str, Any]:
        """Parse .mkm Python file"""
        with open(filename, 'r') as f:
            content = f.read()

        # Execute the .mkm Python code to extract variables
        namespace = {}
        exec(content, namespace)

        rxn_expressions_raw = namespace.get('rxn_expressions', [])
        parsed_reactions = MKMFileParser._parse_rxn_expressions(rxn_expressions_raw)

        return {
            'rxn_expressions': parsed_reactions,
            'rxn_expressions_raw': rxn_expressions_raw,
            'rxn_mechanisms': namespace.get('rxn_mechanisms', {}),
            'species_definitions': namespace.get('species_definitions', {}),
            'surface_names': namespace.get('surface_names', []),
            'prefactor_list': namespace.get('prefactor_list', []),
            'descriptor_names': namespace.get('descriptor_names', []),
            'descriptor_ranges': namespace.get('descriptor_ranges', []),
        }

    @staticmethod
    def _parse_rxn_expressions_yaml(rxn_expressions_dict: Dict[str, str]) -> Dict[int, Dict[str, Any]]:
        """Parse reaction expressions from YAML format"""
        parsed = {}
        
        for step_num, rxn_str in rxn_expressions_dict.items():
            # Remove inline comment
            rxn_str = str(rxn_str).split('#')[0].strip()
            
            # Extract parameters after semicolon
            beta = None
            prefactor = None
            
            if ';' in rxn_str:
                parts = rxn_str.split(';')
                rxn_str = parts[0].strip()
                
                param_str = ';'.join(parts[1:])
                
                # Extract beta
                beta_match = re.search(r'beta\s*=\s*([0-9.]+)', param_str, re.IGNORECASE)
                if beta_match:
                    beta = float(beta_match.group(1))
                
                # Extract prefactor
                prefactor_match = re.search(r'prefactor\s*=\s*([0-9.eE+-]+)', param_str, re.IGNORECASE)
                if prefactor_match:
                    prefactor = float(prefactor_match.group(1))
            
            parsed_rxn = MKMFileParser._parse_single_reaction(rxn_str)
            parsed_rxn['beta'] = beta
            parsed_rxn['prefactor'] = prefactor
            parsed_rxn['step_num'] = step_num
            parsed_rxn['original'] = rxn_expressions_dict[step_num]
            
            parsed[step_num] = parsed_rxn
        
        return parsed

    @staticmethod
    def _parse_rxn_expressions(rxn_expressions_raw: List[str]) -> Dict[int, Dict[str, Any]]:
        """Parse reaction expressions from list format"""
        parsed = {}

        for idx, rxn_str in enumerate(rxn_expressions_raw):
            step_num = idx + 1  # 1-indexed

            # Remove inline comment
            rxn_str = rxn_str.split('#')[0].strip()

            # Extract parameters after semicolon
            beta = None
            prefactor = None
            
            if ';' in rxn_str:
                parts = rxn_str.split(';')
                rxn_str = parts[0].strip()
                
                param_str = ';'.join(parts[1:])
                
                # Extract beta
                beta_match = re.search(r'beta\s*=\s*([0-9.]+)', param_str, re.IGNORECASE)
                if beta_match:
                    beta = float(beta_match.group(1))
                
                # Extract prefactor
                prefactor_match = re.search(r'prefactor\s*=\s*([0-9.eE+-]+)', param_str, re.IGNORECASE)
                if prefactor_match:
                    prefactor = float(prefactor_match.group(1))

            parsed_rxn = MKMFileParser._parse_single_reaction(rxn_str)
            parsed_rxn['beta'] = beta
            parsed_rxn['prefactor'] = prefactor
            parsed_rxn['step_num'] = step_num
            parsed_rxn['original'] = rxn_expressions_raw[idx]

            parsed[step_num] = parsed_rxn

        return parsed

    @staticmethod
    def _parse_single_reaction(rxn_str: str) -> Dict[str, Any]:
        """Parse a single reaction string into its components"""
        if '<->' in rxn_str:
            parts = rxn_str.split('<->')
            reversible = True
        elif '->' in rxn_str:
            parts = rxn_str.split('->')
            reversible = False
        else:
            return {
                'reactants': [],
                'products': [],
                'ts': None,
                'reversible': False,
                'equation': rxn_str
            }

        parts = [p.strip() for p in parts]

        ts = None
        if len(parts) == 3:
            # 'reactants <-> TS <-> products'
            reactants_str = parts[0]
            ts = parts[1]
            products_str = parts[2]
        elif len(parts) == 2:
            # 'reactants <-> products' or 'reactants -> products'
            reactants_str = parts[0]
            products_str = parts[1]
        else:
            reactants_str = parts[0]
            products_str = parts[-1]
            if len(parts) > 2:
                ts = parts[1]

        reactants = [MKMFileParser._normalize_species_name(s) for s in MKMFileParser._parse_species_list(reactants_str)]
        products = [MKMFileParser._normalize_species_name(s) for s in MKMFileParser._parse_species_list(products_str)]
        
        # Filter out None values
        reactants = [r for r in reactants if r is not None]
        products = [p for p in products if p is not None]
        
        if ts:
            ts = [MKMFileParser._normalize_species_name(s) for s in MKMFileParser._parse_species_list(ts)]
            ts = [t for t in ts if t is not None]
            if not ts:
                ts = None

        reactants_eq = ' + '.join(reactants)
        products_eq = ' + '.join(products)
        arrow = '<->' if reversible else '->'
        equation = f"{reactants_eq} {arrow} {products_eq}"

        return {
            'reactants': reactants,
            'products': products,
            'ts': ts,
            'reversible': reversible,
            'equation': equation
        }

    @staticmethod
    def _parse_species_list(species_str: str) -> List[str]:
        """Split a string of species into a list, handling stoichiometric coefficients"""
        species_raw = [s.strip() for s in species_str.split('+')]
        species_expanded = []
        
        for sp in species_raw:
            # Extract stoichiometric coefficient if present
            match = re.match(r'^(\d+\.?\d*)\s*(.+)$', sp)
            if match:
                coeff = float(match.group(1))
                species_name = match.group(2)
                # Expand: add species multiple times based on coefficient
                for _ in range(int(coeff)):
                    species_expanded.append(species_name)
            else:
                species_expanded.append(sp)
        
        return species_expanded

    @staticmethod
    def _normalize_species_name(species: str) -> str:
        """Normalize a species name"""
        # Check if this is an activation energy notation (^X.XXeV)
        if species.startswith('^') and 'eV' in species:
            return None  # This will be filtered out as a special TS marker
        
        # Convert *_t to *
        species = species.replace('*_t', '*')
        
        # Convert X_t to X* (e.g., CO_t -> CO*, H_t -> H*)
        if species.endswith('_t') and not species.endswith('*_t'):
            species = species.replace('_t', '*')
        
        return species


# Backward compatibility function
def parse_mkm_file(filename: str) -> Dict[str, Any]:
    """Backward compatibility function"""
    return MKMFileParser.parse_file(filename)