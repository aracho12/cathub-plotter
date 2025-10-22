"""
Catalysis-Hub API Parser
"""
import json
import requests
from typing import List, Dict, Optional, Tuple

class CatalysisHubParser:
    """Class for parsing data from the Catalysis-hub API"""
    
    BASE_URL = "https://api.catalysis-hub.org/graphql"
    
    def __init__(self):
        pass
    
    def convert_to_cathub_format(self, species: str) -> str:
        """Convert our format to Catalysis-hub format"""
        if species == '*':
            return 'star'
        elif species.endswith('*'):
            return species.replace('*', 'star')
        elif species.endswith('_g'):
            return species.replace('_g', '')
        else:
            return species
    
    def convert_from_cathub_format(self, species: str) -> str:
        """Convert Catalysis-hub format to our format"""
        species_lower = species.lower()
        
        if species_lower == 'star':
            return '*'
        elif 'star' in species_lower:
            return species.replace('star', '').replace('Star', '') + '*'
        elif 'gas' in species_lower:
            return species.replace('gas', '').replace('Gas', '') + '_g'
        else:
            gas_molecules = ['H2', 'O2', 'N2', 'CO', 'CO2', 'H2O', 'CH4', 
                             'NH3', 'NO', 'NO2', 'SO2', 'CH3OH', 'HCOOH']
            if species in gas_molecules:
                return species + '_g'
            return species + '_g'
    
    def parse_sites(self, sites_dict: Dict) -> Dict[str, str]:
        """
        Parse the sites dictionary
        
        Example: {"COstar": "bridge"} -> {"CO*": "bridge"}
        """
        if not sites_dict or not isinstance(sites_dict, dict):
            return {}
        
        parsed_sites = {}
        for cathub_species, site in sites_dict.items():
            our_species = self.convert_from_cathub_format(cathub_species)
            parsed_sites[our_species] = site
        
        return parsed_sites
    
    def parse_species_dict(self, species_str: str) -> Dict[str, float]:
        """Convert a species string to a dictionary"""
        try:
            species_dict = json.loads(species_str)
            
            converted = {}
            for species, coeff in species_dict.items():
                if coeff == 0:
                    continue
                    
                converted_name = self.convert_from_cathub_format(species)
                
                if converted_name in converted:
                    converted[converted_name] += coeff
                else:
                    converted[converted_name] = coeff
            
            return converted
        except (json.JSONDecodeError, TypeError):
            print(f"Failed to parse: {species_str}")
            return {}
    
    def species_dict_to_string(self, species_dict: Dict[str, float]) -> str:
        """Convert a dictionary to a human-readable string"""
        parts = []
        for species, coeff in sorted(species_dict.items()):
            if coeff == 1:
                parts.append(species)
            elif coeff == int(coeff):
                parts.append(f"{int(coeff)} {species}")
            else:
                parts.append(f"{coeff} {species}")
        
        return ' + '.join(parts) if parts else '(empty)'
    
    def get_reaction_equation(self, reaction: Dict) -> str:
        """Return the reaction equation as a string"""
        reactants_dict = self.parse_species_dict(reaction['reactants'])
        products_dict = self.parse_species_dict(reaction['products'])
        
        reactants_str = self.species_dict_to_string(reactants_dict)
        products_str = self.species_dict_to_string(products_dict)
        
        return f"{reactants_str} -> {products_str}"
    
    def validate_reaction_stoichiometry(self, reaction: Dict) -> bool:
        """Check if the reaction stoichiometry is valid"""
        reactants_dict = self.parse_species_dict(reaction['reactants'])
        products_dict = self.parse_species_dict(reaction['products'])
        
        reactants_species = set(reactants_dict.keys()) - {'*'}
        products_species = set(products_dict.keys()) - {'*'}
        
        if len(reactants_species) == 0 or len(products_species) == 0:
            return False
        
        return True
    
    def search_reactions(self, 
                        reactants: List[str] = None, 
                        products: List[str] = None,
                        surface: str = None,
                        facet: str = None,
                        chemical_composition: str = None,
                        dft_code: str = None,
                        dft_functional: str = None,
                        limit: int = 50,
                        validate_stoichiometry: bool = True) -> List[Dict]:
        """
        Search for reactions in Catalysis-hub
        
        Args:
            dft_code: DFT code filter (e.g., 'VASP', 'GPAW')
            dft_functional: DFT functional filter (e.g., 'RPBE', 'PBE')
                            Supports partial matching (e.g., 'RPBE' matches 'BEEF-vdW-RPBE')
        """
        cathub_reactants = None
        cathub_products = None
        
        if reactants:
            converted = [self.convert_to_cathub_format(r) for r in reactants]
            cathub_reactants = '+'.join(converted)
        
        if products:
            converted = [self.convert_to_cathub_format(p) for p in products]
            cathub_products = '+'.join(converted)
        
        filters = []
        
        if cathub_reactants:
            filters.append(f'reactants: "{cathub_reactants}"')
        if cathub_products:
            filters.append(f'products: "{cathub_products}"')
        if chemical_composition:
            filters.append(f'chemicalComposition: "~{chemical_composition}"')
        if facet:
            filters.append(f'facet: "{facet}"')
        if dft_code:
            filters.append(f'dftCode: "~{dft_code}"')
        if dft_functional:
            filters.append(f'dftFunctional: "~{dft_functional}"')
        
        filter_str = ", ".join(filters) if filters else ""
        
        query = f"""
        {{
          reactions(first: {limit}, {filter_str}) {{
            totalCount
            edges {{
              node {{
                id
                reactants
                products
                reactionEnergy
                chemicalComposition
                surfaceComposition
                facet
                sites
                activationEnergy
                dftCode
                dftFunctional
                pubId
                publication {{
                  title
                  authors
                  year
                  doi
                }}
              }}
            }}
          }}
        }}
        """
        
        response = requests.post(self.BASE_URL, json={'query': query})
        data = response.json()
        
        if 'errors' in data:
            print(f"API Error: {data['errors']}")
            return []
        
        total_count = data.get('data', {}).get('reactions', {}).get('totalCount', 0)
        
        reactions = []
        filtered_count = 0
        for edge in data.get('data', {}).get('reactions', {}).get('edges', []):
            reaction = edge['node']
            
            if validate_stoichiometry and not self.validate_reaction_stoichiometry(reaction):
                filtered_count += 1
                continue
            
            reactions.append(reaction)
        
        return reactions
    
    def convert_reaction_to_your_format(self, reaction: Dict) -> Dict:
        """Convert a Catalysis-hub reaction to our format"""
        reactants_dict = self.parse_species_dict(reaction['reactants'])
        products_dict = self.parse_species_dict(reaction['products'])
        
        equation = self.get_reaction_equation(reaction)
        
        # Parse sites
        sites_raw = reaction.get('sites', {})
        if isinstance(sites_raw, str):
            try:
                sites_raw = json.loads(sites_raw)
            except:
                sites_raw = {}
        sites = self.parse_sites(sites_raw)
        
        # Parse publication information
        publication = reaction.get('publication', {}) or {}
        
        return {
            'id': reaction['id'],
            'reactants_dict': reactants_dict,
            'products_dict': products_dict,
            'equation': equation,
            'reaction_energy': reaction['reactionEnergy'],
            'activation_energy': reaction.get('activationEnergy'),
            'surface': reaction.get('surfaceComposition', 'Unknown'),
            'composition': reaction.get('chemicalComposition', 'Unknown'),
            'facet': reaction.get('facet', 'Unknown'),
            'sites': sites,
            'sites_str': ', '.join([f"{sp}: {site}" for sp, site in sites.items()]) if sites else 'N/A',
            'dft_code': reaction.get('dftCode', 'Unknown'),
            'dft_functional': reaction.get('dftFunctional', 'Unknown'),
            'pub_id': reaction.get('pubId', 'Unknown'),
            'pub_title': publication.get('title', 'Unknown'),
            'pub_authors': publication.get('authors', 'Unknown'),
            'pub_year': publication.get('year', 'Unknown'),
            'pub_doi': publication.get('doi', 'Unknown'),
        }