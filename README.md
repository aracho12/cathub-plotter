# Cathub-Plotter

A Python package for plotting free energy diagrams from catalysis-hub.org data and MKM (Microkinetic Model) files.

## Features

- 🔍 **Database Integration**: Access reaction energies from catalysis-hub.org
- 🌡️ **Thermodynamic Calculations**: Convert DFT energies to Gibbs free energies using frequency data
- 📊 **Free Energy Diagrams**: Generate publication-ready free energy diagrams
- ⚙️ **MKM Support**: Parse and visualize reaction mechanisms from .mkm and .yaml files
- 🎨 **Interactive Plots**: Create interactive plots with Plotly
- 🌐 **Web Interface**: Streamlit-based web application for easy use

## Installation

### From Source

```bash
git clone https://github.com/yourusername/cathub-plotter.git
cd cathub-plotter
pip install -e .
```

### Development Installation

```bash
pip install -e .[dev]
```

## Quick Start

### Command Line Interface

```bash
# Plot free energy diagram from MKM file
cathub-plotter plot --mkm-file config.yaml --input-file input.txt

# Search catalysis-hub and calculate free energies
cathub-plotter search --reactants "CO_g,*" --products "CO*" --surface "Cu"
```

### Python API

```python
from cathub_plotter import FreeEnergyCalculator, FreeEnergyDiagramPlotter

# Search and calculate free energies
calc = FreeEnergyCalculator(temperature=298.15)
results = calc.search_and_calculate(
    reactants=['CO_g', '*'],
    products=['CO*'],
    surface='Cu'
)

# Plot free energy diagram
plotter = FreeEnergyDiagramPlotter('config.yaml', 'input.txt')
plotter.plot_mechanism('HER_Heyrovsky')
```

### Web Interface

```bash
streamlit run cathub_plotter/app.py
```

## Usage Examples

### 1. Free Energy Diagram from MKM File

```python
from cathub_plotter import FreeEnergyDiagramPlotter

# Initialize plotter with MKM file and input data
plotter = FreeEnergyDiagramPlotter(
    mkm_file='examples/config.yaml',
    input_file='examples/input.txt',
    temperature=298.15
)

# Plot a specific mechanism
plotter.plot_mechanism('HER_Heyrovsky', save_path='her_diagram.png')
```

### 2. Search Catalysis-Hub Database

```python
from cathub_plotter import FreeEnergyCalculator

# Initialize calculator
calc = FreeEnergyCalculator(temperature=400.0, voltage=0.0)

# Search for CO adsorption reactions
results = calc.search_and_calculate(
    reactants=['CO_g', '*'],
    products=['CO*'],
    surface='Cu',
    facet='111',
    limit=10
)

print(results[['equation', 'delta_G', 'composition', 'facet']])
```

### 3. Compare Multiple Mechanisms

```python
from cathub_plotter import compare_mechanisms

# Compare different HER mechanisms
mechanisms = ['HER_Heyrovsky', 'HER_Tafel']
compare_mechanisms(
    mkm_file='examples/config.yaml',
    input_file='examples/input.txt',
    mechanisms=mechanisms,
    temperature=298.15
)
```

## File Formats

### MKM/YAML Configuration

```yaml
rxn_expressions:
  1: 'H_g + ele_g + *_t <-> H-ele*_t <-> H*_t; beta=0.65'
  2: 'H_g + ele_g + H*_t -> H2-ele*_t -> H2_g + *_t; beta=0.65'

rxn_mechanisms:
  H2:
    HER_Heyrovsky:
      steps: [1, 2]
      color: '#7F7F7F'
```

### Input Data Format

```csv
species_name,status,formation_energy,frequencies,surface_name,site_name
CO_g,gas,-1.23,"[32.9, 48.7, 2150.5]",,,
CO*,ads,-0.85,"[2040, 307, 268, 261, 100, 69]",Cu,bridge
```

## API Reference

### Core Classes

- `FreeEnergyCalculator`: Calculate Gibbs free energies from catalysis-hub data
- `FreeEnergyDiagramPlotter`: Generate free energy diagrams from MKM files
- `CatalysisHubParser`: Parse data from catalysis-hub.org API
- `MKMFileParser`: Parse .mkm and .yaml configuration files

### Key Functions

- `calculate_thermo_correction()`: Calculate thermodynamic corrections
- `parse_mkm_file()`: Parse MKM configuration files
- `plot_mechanism()`: Plot free energy diagram for a mechanism

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this package in your research, please cite:

```bibtex
@software{cathub_plotter,
  title={Cathub-Plotter: Free Energy Diagram Plotting for Catalysis Research},
  author={Your Name},
  year={2024},
  url={https://github.com/yourusername/cathub-plotter}
}
```

## Acknowledgments

- [Catalysis-Hub](https://www.catalysis-hub.org/) for providing the reaction database
- The ASE (Atomic Simulation Environment) project for thermodynamic calculations
- The scientific Python community for excellent tools and libraries
