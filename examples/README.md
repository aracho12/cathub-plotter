# Examples

This directory contains example files and tutorials for using cathub-plotter.

## Files

- `config.yaml` - Example MKM configuration file with reaction mechanisms
- `input.txt` - Example input file with species data
- `frequencies.csv` - Frequency data for thermodynamic calculations
- `tutorial.ipynb` - Jupyter notebook tutorial (coming soon)

## Quick Start

### 1. Plot Free Energy Diagram

```python
from cathub_plotter import FreeEnergyDiagramPlotter

# Initialize plotter
plotter = FreeEnergyDiagramPlotter(
    mkm_file='config.yaml',
    input_file='input.txt',
    temperature=298.15
)

# Plot a specific mechanism
plotter.plot_mechanism('HER_Heyrovsky')
```

### 2. Search Catalysis-Hub

```python
from cathub_plotter import FreeEnergyCalculator

# Initialize calculator
calc = FreeEnergyCalculator(temperature=298.15)

# Search for reactions
results = calc.search_and_calculate(
    reactants=['CO_g', '*'],
    products=['CO*'],
    surface='Cu',
    facet='111'
)

print(results[['equation', 'delta_G', 'composition', 'facet']])
```

### 3. Command Line Usage

```bash
# Plot free energy diagram
cathub-plotter plot --mkm-file config.yaml --input-file input.txt --mechanism HER_Heyrovsky

# Search catalysis-hub
cathub-plotter search --reactants "CO_g,*" --products "CO*" --surface "Cu"
```
