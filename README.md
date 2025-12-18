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
git clone https://github.com/aracho12/cathub-plotter.git
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
cathub-plotter plot --mkm-file config.yaml --input-file input.txt --mechanism HER_Heyrovsky

# Plot for specific surface and facet (e.g., Cu(100))
cathub-plotter plot --mkm-file config.yaml --input-file input.txt --surface Cu --facet 100 --mechanism CO_via_COOH

# Plot at different temperature and voltage
cathub-plotter plot --mkm-file config.yaml --input-file input.txt -T 400 -U -0.5 --mechanism CO_via_COOH

# Compare mechanisms
cathub-plotter compare --mkm-file config.yaml --input-file input.txt \
  --vary mechanism --mechanisms HER_Heyrovsky,HER_Tafel \
  --surface Cu --facet 100 -U -0.5

# Compare voltages (voltage sweep)
# Note: Use = sign and quotes for negative values
cathub-plotter compare --mkm-file config.yaml --input-file input.txt \
  --vary voltage --voltages="-1.0,-0.8,-0.6,-0.4,-0.2,0.0" \
  --mechanism HER_Heyrovsky --surface Cu --facet 100

# 2D comparison: mechanism vs voltage
cathub-plotter compare --mkm-file config.yaml --input-file input.txt \
  --vary mechanism,voltage \
  --mechanisms HER_Heyrovsky,HER_Tafel --voltages="-1.0,-0.5,0.0" \
  --surface Cu --facet 100 --layout subplots

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

## CLI Options

### Plot Command

```bash
cathub-plotter plot [OPTIONS]

Options:
  --mkm-file PATH          Path to MKM/YAML file (required)
  --input-file PATH        Path to input file with species data (required)
  --mechanism TEXT         Specific mechanism to plot (default: all)
  --surface TEXT           Surface name to filter (e.g., Cu)
  --facet TEXT             Facet name to filter (e.g., 100)
  -T, --temperature FLOAT  Temperature in K (default: 298.15)
  -U, --voltage FLOAT      Voltage in V vs RHE (default: 0.0)
  --y-margin FLOAT         Y-axis margin fraction (default: 0.15)
  --ylim YMIN YMAX         Y-axis limits (e.g., --ylim -0.1 0.3)
  -o, --output PATH        Output file path
  --save-data              Save raw data to CSV file
```

**Examples:**

```bash
# Plot HER mechanism at standard conditions
cathub-plotter plot --mkm-file config.yaml --input-file input.txt --mechanism HER_Heyrovsky

# Plot CO2R mechanism for Cu(100) at -0.5 V vs RHE and 400 K
cathub-plotter plot --mkm-file config.yaml --input-file input.txt \
  --surface Cu --facet 100 --mechanism CO2R_CO -U -0.5 -T 400 -o output.png

# Save raw data along with the plot
cathub-plotter plot --mkm-file config.yaml --input-file input.txt \
  --mechanism HER_Heyrovsky --save-data -o her_diagram.png
```

### Compare Command

Compare free energy diagrams across different conditions (mechanisms, temperatures, voltages, surfaces, facets).

```bash
cathub-plotter compare [OPTIONS]

Options:
  --mkm-file PATH          Path to MKM/YAML file (required)
  --input-file PATH        Path to input file with species data (required)
  --vary TEXT              Parameter(s) to vary (required)
                          Single: mechanism, temperature, voltage, surface, facet
                          2D: "mechanism,voltage", "surface,facet", etc.
  
  # Values for varying parameters
  --mechanisms TEXT        Comma-separated mechanism names
  --temperatures TEXT      Comma-separated temperatures in K
  --voltages TEXT          Comma-separated voltages in V
  --surfaces TEXT          Comma-separated surface names
  --facets TEXT            Comma-separated facet names
  
  # Fixed conditions
  --mechanism TEXT         Mechanism (when not varying)
  --surface TEXT           Surface (when not varying)
  --facet TEXT             Facet (when not varying)
  -T, --temperature FLOAT  Temperature in K (when not varying, default: 298.15)
  -U, --voltage FLOAT      Voltage in V (when not varying, default: 0.0)
  
  # Plot options
  --layout {subplots,overlay}  Layout for 2D comparison (default: subplots)
  --colors TEXT            Comma-separated colors
  --labels TEXT            Comma-separated custom labels
  --show-barriers          Show activation barriers
  --show-labels            Show state labels
  --legend-position TEXT   Legend position (default: best)
  -o, --output PATH        Output file path
  --save-data              Save raw comparison data to CSV file
```

**1D Comparison Examples:**

```bash
# Compare different mechanisms
cathub-plotter compare --mkm-file config.yaml --input-file input.txt \
  --vary mechanism \
  --mechanisms HER_Heyrovsky,HER_Tafel \
  --surface Cu --facet 100 -T 298.15 -U -0.5 \
  --show-barriers -o mechanism_comparison.png

# Voltage sweep (very useful for electrochemistry!)
# Note: Use quotes for negative values to avoid parsing issues
cathub-plotter compare --mkm-file config.yaml --input-file input.txt \
  --vary voltage \
  --voltages="-1.0,-0.8,-0.6,-0.4,-0.2,0.0" \
  --mechanism HER_Heyrovsky --surface Cu --facet 100 -T 298.15 \
  --show-barriers -o voltage_sweep.png --save-data

# Temperature dependence
cathub-plotter compare --mkm-file config.yaml --input-file input.txt \
  --vary temperature \
  --temperatures 298.15,400,500,600 \
  --mechanism CO2R_CO --surface Cu --facet 100 -U -0.5 \
  -o temperature_comparison.png

# Metal screening
cathub-plotter compare --mkm-file config.yaml --input-file input.txt \
  --vary surface \
  --surfaces Cu,Ag,Au,Pt \
  --mechanism CO2R_CO --facet 100 -T 298.15 -U -0.5 \
  -o metal_screening.png

# Facet comparison
cathub-plotter compare --mkm-file config.yaml --input-file input.txt \
  --vary facet \
  --facets 100,111,211 \
  --mechanism CO2R_CO --surface Cu -T 298.15 -U -0.5 \
  -o facet_comparison.png
```

**2D Comparison Examples:**

```bash
# Mechanism vs Voltage (subplots layout)
# Note: Use = sign and quotes for negative values
cathub-plotter compare --mkm-file config.yaml --input-file input.txt \
  --vary mechanism,voltage \
  --mechanisms HER_Heyrovsky,HER_Tafel \
  --voltages="-1.0,-0.5,0.0" \
  --surface Cu --facet 100 -T 298.15 \
  --layout subplots --show-barriers -o mech_vs_voltage.png --save-data

# Surface vs Temperature (overlay layout)
cathub-plotter compare --mkm-file config.yaml --input-file input.txt \
  --vary surface,temperature \
  --surfaces Cu,Ag,Au \
  --temperatures 298.15,400,500 \
  --mechanism CO2R_CO --facet 100 -U -0.5 \
  --layout overlay -o surface_vs_temp.png

# Voltage vs Temperature (subplots)
cathub-plotter compare --mkm-file config.yaml --input-file input.txt \
  --vary voltage,temperature \
  --voltages="-1.0,-0.5,0.0" \
  --temperatures 298.15,400,500 \
  --mechanism HER_Heyrovsky --surface Cu --facet 100 \
  --layout subplots -o voltage_vs_temp.png
```

## Usage Examples

### 1. Free Energy Diagram from MKM File

```python
from cathub_plotter import FreeEnergyDiagramPlotter

# Initialize plotter with MKM file and input data
plotter = FreeEnergyDiagramPlotter(
    mkm_file='examples/config.yaml',
    input_file='examples/input.txt',
    temperature=298.15,
    voltage=-0.5  # Voltage in V vs RHE
)

# Plot a specific mechanism
plotter.plot_mechanism('HER_Heyrovsky', save_path='her_diagram.png')

# Compare different voltages
plotter.compare_voltages(
    'HER_Heyrovsky',
    voltages=[-1.0, -0.8, -0.6, -0.4, -0.2, 0.0],
    show_barriers=True
)
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
