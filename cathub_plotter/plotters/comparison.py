"""
Comparison plotting functionality
"""
from typing import List, Dict, Optional
import matplotlib.pyplot as plt
import numpy as np
from .diagram import FreeEnergyDiagramPlotter


def compare_mechanisms(mkm_file: str, 
                      input_file: str,
                      mechanisms: List[str],
                      temperature: float = 298.15,
                      voltage: float = 0.0,
                      save_path: Optional[str] = None) -> None:
    """
    Compare multiple mechanisms on the same plot
    
    Args:
        mkm_file: Path to MKM/YAML file
        input_file: Path to input file with species data
        mechanisms: List of mechanism names to compare
        temperature: Temperature in K
        voltage: Voltage in V
        save_path: Optional path to save the plot
    """
    plotter = FreeEnergyDiagramPlotter(
        mkm_file=mkm_file,
        input_file=input_file,
        temperature=temperature,
        voltage=voltage
    )
    
    # Create comparison plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    colors = plt.cm.Set1(np.linspace(0, 1, len(mechanisms)))
    
    for i, mechanism in enumerate(mechanisms):
        if mechanism not in plotter.rxn_mechanisms:
            print(f"Warning: Mechanism '{mechanism}' not found")
            continue
            
        # Calculate mechanism energies
        energy_data = plotter.calculate_mechanism_energies(mechanism)
        
        # Plot the mechanism
        x_positions = np.arange(len(energy_data['states']))
        ax.plot(x_positions, energy_data['states'], 
                'o-', color=colors[i], linewidth=2, markersize=8,
                label=mechanism)
        
        # Add barriers if they exist
        if energy_data['barriers']:
            barrier_positions = np.arange(len(energy_data['barriers']))
            ax.plot(barrier_positions, energy_data['barriers'], 
                    '^', color=colors[i], markersize=6, alpha=0.7)
    
    ax.set_xlabel('Reaction Coordinate', fontsize=12)
    ax.set_ylabel('Free Energy (eV)', fontsize=12)
    ax.set_title(f'Mechanism Comparison at {temperature}K', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Comparison plot saved to {save_path}")
    
    plt.show()
