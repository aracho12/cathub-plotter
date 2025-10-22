"""
Tests for parsers
"""
import unittest
import tempfile
import os
from pathlib import Path

from cathub_plotter.parsers.mkm import MKMFileParser
from cathub_plotter.parsers.input import InputFileParser


class TestMKMFileParser(unittest.TestCase):
    """Test MKM file parser"""
    
    def setUp(self):
        """Set up test data"""
        self.test_yaml_content = """
rxn_expressions:
  1: 'H_g + ele_g + *_t <-> H-ele*_t <-> H*_t; beta=0.65'
  2: 'H_g + ele_g + H*_t -> H2-ele*_t -> H2_g + *_t; beta=0.65'

rxn_mechanisms:
  H2:
    HER_Heyrovsky:
      steps: [1, 2]
      color: '#7F7F7F'
"""
    
    def test_parse_yaml_file(self):
        """Test parsing YAML file"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(self.test_yaml_content)
            temp_file = f.name
        
        try:
            result = MKMFileParser.parse_file(temp_file)
            
            # Check basic structure
            self.assertIn('rxn_expressions', result)
            self.assertIn('rxn_mechanisms', result)
            self.assertIn('mechanism_colors', result)
            
            # Check reactions
            self.assertEqual(len(result['rxn_expressions']), 2)
            self.assertIn(1, result['rxn_expressions'])
            self.assertIn(2, result['rxn_expressions'])
            
            # Check mechanisms
            self.assertIn('HER_Heyrovsky', result['rxn_mechanisms'])
            self.assertEqual(result['rxn_mechanisms']['HER_Heyrovsky'], [1, 2])
            
            # Check colors
            self.assertIn('HER_Heyrovsky', result['mechanism_colors'])
            self.assertEqual(result['mechanism_colors']['HER_Heyrovsky'], '#7F7F7F')
            
        finally:
            os.unlink(temp_file)


class TestInputFileParser(unittest.TestCase):
    """Test input file parser"""
    
    def setUp(self):
        """Set up test data"""
        self.test_input_content = """species_name	status	formation_energy	frequencies	surface_name	site_name	reference
CO_g	gas	-1.23	[32.9, 48.7, 2150.5]		bridge	Ara
CO*	ads	-0.85	[2040, 307, 268, 261, 100, 69]	Cu	bridge	Ara
"""
    
    def test_parse_input_file(self):
        """Test parsing input file"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(self.test_input_content)
            temp_file = f.name
        
        try:
            result = InputFileParser.parse_file(temp_file)
            
            # Check basic structure
            self.assertIn('CO_g', result)
            self.assertIn('CO*', result)
            
            # Check gas species
            co_g = result['CO_g']
            self.assertEqual(co_g['status'], 'gas')
            self.assertEqual(co_g['formation_energy'], -1.23)
            self.assertEqual(co_g['frequencies'], [32.9, 48.7, 2150.5])
            
            # Check adsorbed species
            co_ads = result['CO*']
            self.assertEqual(co_ads['status'], 'ads')
            self.assertEqual(co_ads['formation_energy'], -0.85)
            self.assertEqual(co_ads['surface_name'], 'Cu')
            self.assertEqual(co_ads['site_name'], 'bridge')
            
        finally:
            os.unlink(temp_file)


if __name__ == '__main__':
    unittest.main()
