"""
Free Energy Diagram Plotter
Plot potential energy diagrams from mkm files using calculated energies
"""
from ..parsers.mkm import MKMFileParser
from ..core.mkm_calculator import MKMEnergyCalculator
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
from scipy.interpolate import make_interp_spline
from typing import List, Dict, Optional
import pandas as pd
import os

# Configure Helvetica font
def setup_helvetica_font():
    """Setup Helvetica font for matplotlib"""
    try:
        # Try to use system Helvetica
        plt.rcParams['font.family'] = 'Helvetica'
        plt.rcParams['font.sans-serif'] = ['Helvetica']
    except:
        # Try to load custom Helvetica font
        font_path = 'resource/Helvetica.ttf'
        if os.path.exists(font_path):
            font_prop = fm.FontProperties(fname=font_path)
            plt.rcParams['font.family'] = font_prop.get_name()
        else:
            # Fallback to Arial or sans-serif
            plt.rcParams['font.family'] = 'sans-serif'
            plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']

# Setup font at module import
setup_helvetica_font()


class FreeEnergyDiagramPlotter:
    """
    Generate free energy diagrams from reaction mechanisms
    """
    
    def __init__(self, 
                 mkm_file: str,
                 input_file: str = 'input.txt',
                 temperature: float = 298.15,
                 voltage: float = 0.0,
                 use_gibbs_energy: bool = False):
        """
        Initialize plotter
        
        Args:
            mkm_file: Path to .mkm file
            input_file: Path to input.txt with species data
            temperature: Temperature in K
            voltage: Voltage in V
            use_gibbs_energy: If True, formation_energy in input.txt is treated as Gibbs free energy at U=0V
                            (no thermodynamic corrections applied, only voltage correction)
        """
        self.mkm_file = mkm_file
        self.mkm_data = MKMFileParser.parse_file(mkm_file)
        self.calculator = MKMEnergyCalculator(
            input_file=input_file,
            temperature=temperature,
            voltage=voltage,
            use_gibbs_energy=use_gibbs_energy
        )
        self.temperature = temperature
        self.voltage = voltage
        self.use_gibbs_energy = use_gibbs_energy
        
        # Parse data
        self.rxn_expressions = self.mkm_data['rxn_expressions']
        self.rxn_mechanisms = self.mkm_data['rxn_mechanisms']
    
    def calculate_mechanism_energies(self, mechanism_name: str, verbose: bool = False) -> Dict:
        """
        Calculate energies for all steps in a mechanism
        
        Args:
            mechanism_name: Name of mechanism (e.g., 'C1_via_CO-H-ele')
            verbose: Print detailed information
        
        Returns:
            Dictionary with energy data for the mechanism
        """
        if mechanism_name not in self.rxn_mechanisms:
            raise ValueError(f"Mechanism '{mechanism_name}' not found in mkm file")
        
        step_indices = self.rxn_mechanisms[mechanism_name]
        
        if verbose:
            print(f"\n{'='*70}")
            print(f"Calculating mechanism: {mechanism_name}")
            print(f"Steps: {step_indices}")
            print(f"Temperature: {self.temperature} K, Voltage: {self.voltage} V")
            print(f"{'='*70}\n")
        
        # Store results
        states = []  # List of state energies
        barriers = []  # List of barrier energies
        labels = []  # List of state labels
        state_species = []  # List of species at each state
        is_electrochemical = []  # List indicating if each step is electrochemical
        delta_G_values = []  # List of reaction free energies
        Ga_values = []  # List of activation energies
        n_electrons = []  # List of cumulative electron transfers (for x-axis)
        
        # Initial state (reference = 0)
        cumulative_energy = 0.0
        cumulative_electrons = 0  # Start with n(H+ + e-) = 0
        states.append(cumulative_energy)
        labels.append('Initial')
        n_electrons.append(cumulative_electrons)
        
        # Calculate energy for each step
        for idx, step_num in enumerate(step_indices):
            rxn_data = self.rxn_expressions[step_num]
            
            if verbose:
                print(f"\nStep {idx+1} (Reaction {step_num}): {rxn_data['equation']}")
            
            # Get beta value
            beta = rxn_data.get('beta', 0.5)
            if beta is None:
                beta = 0.5
            
            # Calculate reaction energy
            result = self.calculator.calculate_reaction_energy(
                reactants=rxn_data['reactants'],
                products=rxn_data['products'],
                ts=rxn_data['ts'],
                beta=beta,
                verbose=verbose
            )
            
            # Update cumulative energy
            delta_G = result['delta_G']
            cumulative_energy += delta_G
            
            # Check if electrochemical and count electrons
            # Count all electrons (ele_g or pe_g) in reactants OR products
            electron_count_reactants = sum(1 for r in rxn_data['reactants'] if r in ['ele_g', 'pe_g'])
            electron_count_products = sum(1 for p in rxn_data['products'] if p in ['ele_g', 'pe_g'])
            electron_count = electron_count_reactants  # For cumulative count, use reactants
            has_electron = (electron_count_reactants > 0) or (electron_count_products > 0)
            is_electrochemical.append(has_electron)
            
            # Update cumulative electron count with actual number of electrons
            # This handles stoichiometric coefficients like 5ele_g correctly
            cumulative_electrons += electron_count
            
            # Store energies
            delta_G_values.append(delta_G)
            
            # Add barrier if TS exists and has valid energy
            if result['ts'] and result.get('Ga') is not None:
                Ga = result['Ga']
                barrier_energy = states[-1] + Ga  # Absolute barrier height
                barriers.append(barrier_energy)
                Ga_values.append(Ga)
            else:
                barriers.append(None)
                Ga_values.append(None)
            
            # Add final state
            states.append(cumulative_energy)
            n_electrons.append(cumulative_electrons)
            
            # Create label
            products = rxn_data['products']
            label = self._make_label(products)
            labels.append(label)
            state_species.append(products)
        
        return {
            'mechanism_name': mechanism_name,
            'states': states,
            'barriers': barriers,
            'labels': labels,
            'state_species': state_species,
            'step_indices': step_indices,
            'is_electrochemical': is_electrochemical,
            'delta_G_values': delta_G_values,
            'Ga_values': Ga_values,
            'n_electrons': n_electrons,  # Number of (H+ + e-) transfers at each state
            'temperature': self.temperature,
            'voltage': self.voltage
        }
    
    def _calculate_x_offsets(self, x_base: np.ndarray, states: List[float], 
                             offset_val: float = 0.5) -> np.ndarray:
        """
        Calculate x-axis offsets for consecutive states with the same n_electrons value
        
        For A->B transition where both have same n:
        - A gets -offset_val (left)
        - B stays at 0 (center)
        - Next state C starts at its original position (not affected)
        
        Args:
            x_base: Base x positions (n_electrons)
            states: Energy values
            offset_val: Offset value for the initial state
        
        Returns:
            Array of x positions with offsets applied
        """
        x_positions = x_base.copy().astype(float)
        
        # Check consecutive pairs
        for i in range(len(x_base) - 1):
            if x_base[i] == x_base[i+1]:
                # A and B have same n
                # A (initial state) gets -offset_val
                x_positions[i] = x_base[i] - offset_val
                # B (final state) stays at original position (0 offset)
                x_positions[i+1] = x_base[i+1]
        
        return x_positions
    
    def _format_subscripts(self, text: str) -> str:
        """
        Convert numbers in chemical formulas to subscripts
        e.g., 'CO2(g)' -> 'CO$_2$(g)', 'CH4(g)' -> 'CH$_4$(g)'
        But keep leading numbers: '2CO2(g)' -> '2CO$_2$(g)'
        
        Args:
            text: Text to format
        
        Returns:
            Formatted text with subscripts
        """
        import re
        
        # Pattern: non-digit followed by digit(s) not at the start
        # This will match O2, H4, etc. but not leading numbers like "2CO2"
        def replace_subscript(match):
            before = match.group(1)  # Character before number
            number = match.group(2)  # The number
            return f'{before}$_{{{number}}}$'
        
        # Match: (letter or *)+(digit+) but not at start of string
        # Use lookahead/lookbehind to preserve context
        result = re.sub(r'([A-Za-z\*])(\d+)', replace_subscript, text)
        
        return result
    
    def _normalize_label(self, label: str) -> str:
        """
        Normalize label by removing gas phase products and format subscripts
        e.g., 'CH* + H2O(g)' -> 'CH*'
        e.g., 'CO2(g)' -> 'CO$_2$(g)'
        e.g., 'C2H4(g) + H2O(g)' -> 'C2H4(g)' (removes H2O(g))
        
        Args:
            label: Original label string
        
        Returns:
            Normalized label string with subscripts
        """
        if not label or label == 'Initial':
            return label
        
        # Split by ' + ' and filter out gas phase species
        parts = [p.strip() for p in label.split('+')]
        
        # Filter out H2O(g) specifically
        parts = [p for p in parts if p != 'H2O(g)']
        
        # Then filter out other gas phase species if there are surface species
        surface_species = [p for p in parts if not p.endswith('(g)')]
        
        if not surface_species:
            # If only gas species, return with subscripts but no H2O(g)
            result = ' + '.join(parts)
            return self._format_subscripts(result)
        
        result = ' + '.join(surface_species)
        return self._format_subscripts(result)
    
    def _make_label(self, species_list: List[str]) -> str:
        """
        Create a label for a state from species list
        
        Args:
            species_list: List of species
        
        Returns:
            Formatted label string
        """
        # Filter out empty sites
        species = [s for s in species_list if s != '*']
        
        if not species:
            return '*'
        
        # Format species names
        formatted = []
        for sp in species:
            if sp.endswith('_g'):
                name = sp[:-2] + '(g)'
            else:
                name = sp.replace('*', '') + '*'
            formatted.append(name)
        
        return ' + '.join(formatted)
    
    def plot_mechanism(self, 
                      mechanism_name: str,
                      ax: Optional[plt.Axes] = None,
                      color: str = 'blue',
                      label: Optional[str] = None,
                      show_barriers: bool = True,
                      show_labels: bool = True,
                      show_legend: bool = True,
                      verbose: bool = False,
                      bar_width_small: float = 0.25,
                      bar_width_large: float = 0.3,
                      same_n_offset: float = 0.4,
                      initial_label: Optional[str] = None,
                      label_offset_crowded: float = 0,
                      use_annotations: bool = False,
                      show_rds: bool = False,
                      rds_linewidth: float = 4.0,
                      legend_position: str = 'right',
                      width_per_n: float = 1.5,
                      base_width: float = 4.0,
                      height: float = 6.0,
                      fontsize_label: int = 11,
                      fontsize_axis: int = 14,
                      fontsize_title: int = 15,
                      fontsize_legend: int = 11,
                      fontsize_ticks: int = 12,
                      y_margin_fraction: float = 0.15,
                      show_voltage_text: bool = False,
                      zorder: int = 2,
                      line_alpha: float = 1.0,
                      ylim: Optional[tuple] = None,
                      x_axis_mode: str = 'step',
                      save_path: Optional[str] = None,
                      save_data: bool = True) -> plt.Figure:
        """
        Plot free energy diagram for a mechanism
        
        Args:
            mechanism_name: Name of mechanism to plot
            ax: Matplotlib axes object (creates new if None)
            color: Color for the plot
            label: Label for the legend
            show_barriers: Whether to show activation barriers
            show_labels: Whether to show state labels
            show_legend: Whether to show legend
            verbose: Print detailed information
            bar_width_small: Width for bars in same-n transitions (default: 0.25)
            bar_width_large: Width for regular bars (default: 0.3)
            same_n_offset: X-offset for initial state in same-n transitions (default: 0.4)
            initial_label: Custom label for initial state (e.g., 'CO(g)') instead of 'Initial'
            label_offset_crowded: (Deprecated, not used)
            use_annotations: (Deprecated, not used)
            show_rds: Whether to highlight the Rate Determining Step (largest ΔG)
            rds_linewidth: Line width for the RDS step (default: 4.0)
            legend_position: Legend position ('right', 'upper right', 'best', etc.) or None to disable
            width_per_n: Width per n(H+ + e-) unit (default: 1.5)
            base_width: Minimum figure width (default: 4.0)
            height: Figure height (default: 6.0)
            fontsize_label: Font size for state labels (default: 11)
            fontsize_axis: Font size for axis labels (default: 14)
            fontsize_title: Font size for title (default: 15)
            fontsize_legend: Font size for legend (default: 11)
            fontsize_ticks: Font size for tick labels (default: 12)
            y_margin_fraction: Fraction of y-range to add as margin (default: 0.15)
            show_voltage_text: Whether to show voltage text in plot (default: False)
            zorder: Drawing order (higher values drawn on top, default: 2)
            line_alpha: Alpha transparency for lines (0-1, default: 1.0)
            ylim: Y-axis limits as tuple (ymin, ymax). If None, automatically calculated.
            x_axis_mode: X-axis mode - 'step' for reaction coordinate (default) or 'electron' for n(H+ + e-)
            save_data: Whether to save raw data to CSV file (default: True)
        
        Returns:
            Matplotlib figure object
        """
        # Calculate energies
        data = self.calculate_mechanism_energies(mechanism_name, verbose=verbose)
        
        states = data['states']
        barriers = data['barriers']
        labels = data['labels'].copy()  # Make a copy to modify
        is_electrochemical = data.get('is_electrochemical', [])
        n_electrons = data['n_electrons']  # Get cumulative electron transfer counts
        delta_G_values = data.get('delta_G_values', [])
        
        # Find RDS (Rate Determining Step) - step with largest ΔG
        rds_step = None
        if show_rds and delta_G_values:
            rds_step = np.argmax(delta_G_values)
        
        # Replace 'Initial' label if custom label provided
        if initial_label and labels[0] == 'Initial':
            labels[0] = initial_label
        
        # Plot energy levels
        n_states = len(states)
        
        # Determine x-axis mode
        if x_axis_mode == 'electron':
            # Use n(H+ + e-) as x-axis
            x_base = np.array(n_electrons)
            # Calculate x offsets for states with same n_electrons
            x_positions = self._calculate_x_offsets(x_base, states, offset_val=same_n_offset)
            
            # Calculate line widths for each state
            # States in same-n transitions get bar_width_small, others get bar_width_large
            line_widths = []
            for i in range(n_states):
                # Check if this state is part of a same-n transition
                is_same_n_transition = False
                if i < n_states - 1 and x_base[i] == x_base[i+1]:
                    is_same_n_transition = True  # A in A->B
                elif i > 0 and x_base[i-1] == x_base[i]:
                    is_same_n_transition = True  # B in A->B
                
                line_widths.append(bar_width_small if is_same_n_transition else bar_width_large)
            
            # Calculate dynamic width based on electron range
            x_range = max(n_electrons) - min(n_electrons)
            dynamic_width = base_width + width_per_n * x_range
        else:  # x_axis_mode == 'step' (default)
            # Use mechanism step order as x-axis (0, 1, 2, 3, ...)
            x_base = np.arange(n_states)
            x_positions = x_base.astype(float)
            # All states use the same line width
            line_widths = [bar_width_large] * n_states
            
            # Calculate dynamic width based on number of steps
            n_steps = len(states) - 1
            dynamic_width = base_width + width_per_n * n_steps
        
        # Create figure if needed
        if ax is None:
            fig, ax = plt.subplots(figsize=(dynamic_width, height))
        else:
            fig = ax.get_figure()
        
        # Draw horizontal lines for each state
        for i in range(n_states):
            line_width = line_widths[i]
            x_start = x_positions[i] - line_width/2
            x_end = x_positions[i] + line_width/2
            ax.plot([x_start, x_end], [states[i], states[i]], 
                   color=color, linewidth=2, solid_capstyle='round', alpha=line_alpha, zorder=zorder)
        
        # Connect states with lines
        for i in range(n_states - 1):
            x_start = x_positions[i] + line_widths[i]/2
            x_end = x_positions[i+1] - line_widths[i+1]/2
            
            # Determine line style (electrochemical: solid, thermal: dashed)
            is_echem = is_electrochemical[i] if i < len(is_electrochemical) else False
            linestyle = '-' if is_echem else '--'
            
            # Determine line width (thicker for RDS)
            connection_linewidth = rds_linewidth if (show_rds and i == rds_step) else 1.5
            
            if show_barriers and barriers[i] is not None:
                # Draw smooth barrier using spline
                x_barrier = (x_start + x_end) / 2
                barrier_height = barriers[i]
                
                # Create smooth curve through IS -> TS -> FS
                x_points = np.array([x_start, x_barrier, x_end])
                y_points = np.array([states[i], barrier_height, states[i+1]])
                
                # Generate smooth spline
                x_smooth = np.linspace(x_start, x_end, 50)
                try:
                    spl = make_interp_spline(x_points, y_points, k=2)
                    y_smooth = spl(x_smooth)
                    ax.plot(x_smooth, y_smooth, color=color, linewidth=connection_linewidth, 
                           linestyle=linestyle, alpha=line_alpha, zorder=zorder)
                except:
                    # Fallback to straight lines if spline fails
                    ax.plot([x_start, x_barrier], [states[i], barrier_height],
                           color=color, linewidth=connection_linewidth, linestyle=linestyle, alpha=line_alpha, zorder=zorder)
                    ax.plot([x_barrier, x_end], [barrier_height, states[i+1]],
                           color=color, linewidth=connection_linewidth, linestyle=linestyle, alpha=line_alpha, zorder=zorder)
            else:
                # Direct connection (no barrier)
                ax.plot([x_start, x_end], [states[i], states[i+1]],
                       color=color, linewidth=connection_linewidth, linestyle=linestyle, alpha=line_alpha, zorder=zorder)
        
        # Add labels
        if show_labels:
            # Calculate y_offset based on data range (before setting ylim)
            y_min_data = min(states)
            y_max_data = max(states)
            y_range_data = y_max_data - y_min_data
            y_offset_up = y_range_data * 0.015   # 1.5% above
            y_offset_down = -y_range_data * 0.02  # 2% below
            
            # Track shown labels at each x position to avoid duplicates
            if not hasattr(ax, '_shown_labels_at_position'):
                ax._shown_labels_at_position = set()
            
            for i, (x, energy, lbl) in enumerate(zip(x_positions, states, labels)):
                # Normalize label (remove gas phase products and format subscripts)
                normalized_lbl = self._normalize_label(lbl)
                
                # Round x position to avoid floating point precision issues
                x_rounded = round(x, 2)
                label_key = (x_rounded, normalized_lbl)
                
                # Only show if this label hasn't been shown at this x position
                if label_key not in ax._shown_labels_at_position:
                    ax._shown_labels_at_position.add(label_key)
                    
                    # Special positioning for specific species at n=2
                    # HCOOH(g) at n=2 → above bar
                    # CO(g), H2(g) at n=2 → below bar
                    if abs(x_rounded - 2.0) < 0.1:  # at n=2
                        if 'HCOOH' in normalized_lbl or 'HCOO' in normalized_lbl:
                            # HCOOH above
                            y_offset = y_offset_up
                            va = 'bottom'
                        elif 'CO' in normalized_lbl and 'HCOO' not in normalized_lbl and 'COH' not in normalized_lbl:
                            # CO (not HCOOH, not COH) below
                            y_offset = y_offset_down
                            va = 'top'
                        elif 'H$_{2}$' in normalized_lbl or 'H2' in normalized_lbl:
                            # H2 below
                            y_offset = y_offset_down
                            va = 'top'
                        else:
                            # Default alternating
                            y_offset = y_offset_up if i % 2 == 0 else y_offset_down
                            va = 'bottom' if i % 2 == 0 else 'top'
                    else:
                        # Default positioning for other positions
                        y_offset = y_offset_up if i % 2 == 0 else y_offset_down
                        va = 'bottom' if i % 2 == 0 else 'top'
                    
                    ax.text(x, energy + y_offset, normalized_lbl, 
                           ha='center', va=va,
                           fontsize=fontsize_label, rotation=0, color=color, alpha=line_alpha)
        
        # Formatting
        if x_axis_mode == 'electron':
            ax.set_xlabel('n(H⁺ + e⁻)', fontsize=fontsize_axis)
            # Use unique n_electrons for ticks (without offsets)
            unique_n = sorted(set(x_base))
            ax.set_xticks(unique_n)
            ax.set_xticklabels([f'{int(n)}' for n in unique_n])
        else:  # x_axis_mode == 'step'
            ax.set_xlabel('Reaction Coordinate', fontsize=fontsize_axis)
            # Use step numbers for ticks
            ax.set_xticks(x_base)
            ax.set_xticklabels([f'{int(n)}' for n in x_base])
        
        ax.set_ylabel('G (eV)', fontsize=fontsize_axis)
        
        # Set tick label font sizes
        ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
        
        # Add horizontal line at y=0
        ax.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.3, zorder=0)
        
        # Add vertical dotted grid lines between steps (at 0.5, 1.5, 2.5, ...)
        if n_states > 1:
            x_min = min(x_base)
            x_max = max(x_base)
            grid_positions = np.arange(0.5, x_max + 0.5, 1.0)
            for x_grid in grid_positions:
                if x_min <= x_grid <= x_max:
                    ax.axvline(x=x_grid, color='gray', linestyle=':', linewidth=1, alpha=0.5, zorder=0)
        
        # Set y-axis limits with margin to prevent label clipping
        if ylim is not None:
            # Use user-specified ylim
            ax.set_ylim(ylim[0], ylim[1])
        else:
            # Calculate the actual min/max including label positions and barriers (TS)
            y_min_data = min(states)
            y_max_data = max(states)
            
            # Include barriers (transition states) in max calculation
            valid_barriers = [b for b in barriers if b is not None]
            if valid_barriers:
                y_max_data = max(y_max_data, max(valid_barriers))
            
            y_range = y_max_data - y_min_data
            
            if show_labels:
                # Calculate where labels actually appear
                y_offset_up = y_range * 0.015
                y_offset_down = -y_range * 0.02  # negative because it's below
                text_height = y_range * 0.03  # approximate text height
                
                # Find the actual min and max y positions including labels
                y_min_with_labels = y_min_data + y_offset_down - text_height/2
                y_max_with_labels = y_max_data + y_offset_up + text_height/2
                
                # Add base margin
                y_margin = y_range * y_margin_fraction
                ax.set_ylim(y_min_with_labels - y_margin, y_max_with_labels + y_margin)
            else:
                # No labels, just use data range with margin
                y_margin = y_range * y_margin_fraction
                ax.set_ylim(y_min_data - y_margin, y_max_data + y_margin)
        
        # Add label to legend if provided (even if show_legend=False)
        # This allows external functions to collect labels and show legend later
        if label:
            ax.plot([], [], color=color, linewidth=2, label=label)
        
        # Add legend if requested
        if show_legend and legend_position is not None:
            # Add electrochemical/thermal legend
            ax.plot([], [], 'k-', linewidth=1.5, label='electrochemical', alpha=0.7)
            ax.plot([], [], 'k--', linewidth=1.5, label='chemical', alpha=0.7)
            
            # Determine legend location based on legend_position
            if legend_position == 'right':
                # Place legend outside to the right
                ax.legend(fontsize=fontsize_legend, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
            elif legend_position == 'upper':
                # Place legend above the plot
                ax.legend(fontsize=fontsize_legend, loc='lower left', bbox_to_anchor=(0, 1.02), 
                         ncol=3, frameon=False)
            else:
                # Use standard matplotlib location
                ax.legend(fontsize=fontsize_legend, loc=legend_position, frameon=False)
        
        # Add voltage text inside the plot
        if show_voltage_text:
            voltage_text = f'U = {self.voltage:.2f} V vs. RHE'
            ax.text(0.02, 0.98, voltage_text, transform=ax.transAxes,
                   fontsize=fontsize_ticks, verticalalignment='top',
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='gray', alpha=0.8))
        
        # Title
        title = f'{mechanism_name} at U={self.voltage:.2f} V vs RHE'
        ax.set_title(title, fontsize=fontsize_title, fontweight='bold')
        
        plt.tight_layout()
        
        # Save figure if path is provided
        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")
            
            # Save raw data if requested
            if save_data:
                # Determine CSV filename from save_path
                csv_path = save_path.rsplit('.', 1)[0] + '_data.csv'
                self.save_mechanism_data(data, csv_path, x_axis_mode=x_axis_mode)
                print(f"Raw data saved to {csv_path}")
        
        return fig
    
    def compare_voltages(self,
                        mechanism_name: str,
                        voltages: List[float],
                        colors: Optional[List[str]] = None,
                        figsize: Optional[tuple] = None,
                        show_barriers: bool = True,
                        show_labels: bool = False,
                        legend_position: str = 'right',
                        width_per_n: float = 1.5,
                        base_width: float = 4.0,
                        height: float = 6.0,
                        fontsize_label: int = 11,
                        fontsize_axis: int = 14,
                        fontsize_title: int = 15,
                        fontsize_legend: int = 11,
                        fontsize_ticks: int = 12,
                        y_margin_fraction: float = 0.15,
                        x_axis_mode: str = 'step') -> plt.Figure:
        """
        Compare free energy diagrams at different voltages
        
        Args:
            mechanism_name: Name of mechanism to plot
            voltages: List of voltages to compare
            colors: List of colors for each voltage
            figsize: Figure size (if None, calculated dynamically)
            show_barriers: Whether to show activation barriers
            show_labels: Whether to show state labels
            legend_position: Legend position ('right', 'upper', 'best', etc.)
            width_per_n: Width per n(H+ + e-) unit (default: 1.5)
            base_width: Minimum figure width (default: 4.0)
            height: Figure height (default: 6.0)
            fontsize_label: Font size for state labels
            fontsize_axis: Font size for axis labels
            fontsize_title: Font size for title
            fontsize_legend: Font size for legend
            fontsize_ticks: Font size for tick labels
            y_margin_fraction: Fraction of y-range to add as margin
            x_axis_mode: X-axis mode - 'step' for reaction coordinate (default) or 'electron' for n(H+ + e-)
        
        Returns:
            Matplotlib figure object
        """
        if colors is None:
            # Generate colors
            try:
                cmap = plt.colormaps['coolwarm']
            except:
                cmap = plt.cm.get_cmap('coolwarm')
            v_normalized = (np.array(voltages) - min(voltages)) / (max(voltages) - min(voltages) + 1e-10)
            colors = [cmap(v) for v in v_normalized]
        
        # Calculate dynamic figure size if not provided
        if figsize is None:
            data = self.calculate_mechanism_energies(mechanism_name, verbose=False)
            if x_axis_mode == 'electron':
                n_electrons = data['n_electrons']
                x_range = max(n_electrons) - min(n_electrons)
                dynamic_width = base_width + width_per_n * x_range
            else:  # x_axis_mode == 'step'
                n_steps = len(data['states']) - 1
                dynamic_width = base_width + width_per_n * n_steps
            figsize = (dynamic_width, height)
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Plot each voltage
        for idx, (voltage, color) in enumerate(zip(voltages, colors)):
            # Update calculator voltage
            self.calculator.voltage = voltage
            self.calculator.state_calculator.voltage = voltage
            self.voltage = voltage
            
            # Show labels only for the first voltage (to avoid clutter)
            show_labels_here = show_labels and (idx == 0)
            
            # Plot
            self.plot_mechanism(
                mechanism_name,
                ax=ax,
                color=color,
                label=f'U = {voltage:.2f} V',
                show_barriers=show_barriers,
                show_labels=show_labels_here,  # Labels only once
                show_legend=False,  # Will add custom legend later
                verbose=False,
                legend_position=None,  # Don't show legend per mechanism
                fontsize_label=fontsize_label,
                fontsize_axis=fontsize_axis,
                fontsize_title=fontsize_title,
                fontsize_legend=fontsize_legend,
                fontsize_ticks=fontsize_ticks,
                y_margin_fraction=y_margin_fraction,
                show_voltage_text=False,  # Don't show voltage text when comparing multiple voltages
                x_axis_mode=x_axis_mode
            )
        
        # Add legend
        handles, labels_list = ax.get_legend_handles_labels()
        if handles:
            if legend_position == 'right':
                ax.legend(fontsize=fontsize_legend, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
            elif legend_position == 'upper':
                ax.legend(fontsize=fontsize_legend, loc='lower left', bbox_to_anchor=(0, 1.02), 
                         ncol=min(len(voltages), 4), frameon=False)
            else:
                ax.legend(fontsize=fontsize_legend, loc=legend_position, frameon=False)
        ax.set_title(f'{mechanism_name} - Voltage Comparison\nT = {self.temperature:.2f} K',
                    fontsize=fontsize_title, fontweight='bold')
        
        plt.tight_layout()
        return fig
    
    def _align_mechanisms(self, mechanism_names: List[str]) -> Dict:
        """
        Align mechanism steps to show common intermediates at same x-position
        
        Args:
            mechanism_names: List of mechanism names to align
        
        Returns:
            Dict with aligned steps for each mechanism
        """
        if len(mechanism_names) < 2:
            # No alignment needed for single mechanism
            steps_dict = {}
            for name in mechanism_names:
                steps_dict[name] = self.rxn_mechanisms[name]
            return steps_dict
        
        # Get all mechanism steps
        all_steps = {name: self.rxn_mechanisms[name] for name in mechanism_names}
        
        # Find common suffix (steps at the end that are identical)
        min_len = min(len(steps) for steps in all_steps.values())
        common_suffix_len = 0
        
        for i in range(1, min_len + 1):
            # Check if last i steps are the same in all mechanisms
            last_steps = [steps[-i] for steps in all_steps.values()]
            if len(set(last_steps)) == 1:  # All same
                common_suffix_len = i
            else:
                break
        
        # Build aligned sequences
        aligned = {}
        max_prefix_len = max(len(steps) - common_suffix_len for steps in all_steps.values())
        
        for name, steps in all_steps.items():
            prefix = steps[:-common_suffix_len] if common_suffix_len > 0 else steps
            suffix = steps[-common_suffix_len:] if common_suffix_len > 0 else []
            
            # Pad prefix to align suffix
            padding_needed = max_prefix_len - len(prefix)
            aligned_steps = prefix + [None] * padding_needed + suffix
            aligned[name] = aligned_steps
        
        return aligned
    
    def _plot_aligned_mechanisms(self, ax, mechanism_names, aligned_steps, colors,
                                 show_barriers, show_labels,
                                 bar_width_small=0.25, bar_width_large=0.5, same_n_offset=0.5,
                                 initial_label=None, fontsize_label=11, fontsize_axis=14,
                                 fontsize_ticks=12, y_margin_fraction=0.15, x_axis_mode='step'):
        """
        Plot aligned mechanisms with user-defined placeholders (None or 'x')
        Placeholders are skipped and states are connected naturally
        Each mechanism uses n(H+ + e-) as x-axis positions
        """
        # Store all labels at each n_electron position for deduplication
        labels_at_position = {}  # n_electron -> {label: [(energy, mech_idx)]}
        
        # Process each mechanism
        all_mechanism_data = []
        for mech_idx, (mech_name, color) in enumerate(zip(mechanism_names, colors)):
            steps = aligned_steps[mech_name]
            
            # Calculate energies for each step
            states = [0.0]  # Start at 0
            labels = [initial_label if initial_label else 'Initial']
            barriers = []
            is_electrochemical = []
            n_electrons_list = [0]  # Start with n(H+ + e-) = 0
            cumulative_energy = 0.0
            cumulative_electrons = 0
            
            for step_num in steps:
                # Check if placeholder (None or 'x')
                if step_num is None or step_num == 'x':
                    # Placeholder - skip this position (no horizontal line here)
                    states.append(None)  # Mark as skipped
                    labels.append('')
                    barriers.append(None)
                    is_electrochemical.append(False)
                    n_electrons_list.append(cumulative_electrons)  # Keep same n
                else:
                    # Calculate actual step
                    rxn_data = self.rxn_expressions[step_num]
                    beta = rxn_data.get('beta', 0.5) or 0.5
                    
                    result = self.calculator.calculate_reaction_energy(
                        reactants=rxn_data['reactants'],
                        products=rxn_data['products'],
                        ts=rxn_data['ts'],
                        beta=beta,
                        verbose=False
                    )
                    
                    delta_G = result['delta_G']
                    cumulative_energy += delta_G
                    
                    # Check if electrochemical and count electrons
                    # Count all electrons (ele_g or pe_g) in reactants OR products
                    electron_count_reactants = sum(1 for r in rxn_data['reactants'] if r in ['ele_g', 'pe_g'])
                    electron_count_products = sum(1 for p in rxn_data['products'] if p in ['ele_g', 'pe_g'])
                    electron_count = electron_count_reactants  # For cumulative count, use reactants
                    has_electron = (electron_count_reactants > 0) or (electron_count_products > 0)
                    is_electrochemical.append(has_electron)
                    
                    # Update cumulative electron count with actual number of electrons
                    # This handles stoichiometric coefficients like 5ele_g correctly
                    cumulative_electrons += electron_count
                    
                    if result['ts'] and result.get('Ga') is not None:
                        Ga = result['Ga']
                        barrier_energy = cumulative_energy - delta_G + Ga  # Relative to previous real state
                        barriers.append(barrier_energy)
                    else:
                        barriers.append(None)
                    
                    states.append(cumulative_energy)
                    n_electrons_list.append(cumulative_electrons)
                    products = rxn_data['products']
                    label = self._make_label(products)
                    labels.append(label)
            
            all_mechanism_data.append({
                'name': mech_name,
                'color': color,
                'states': states,
                'labels': labels,
                'barriers': barriers,
                'is_electrochemical': is_electrochemical,
                'steps': steps,
                'n_electrons': n_electrons_list
            })
            
            # Labels will be added after plotting, with actual x positions
        
        # Store label information for each mechanism
        label_info = []  # List of (x_pos, energy, label, color)
        
        # Plot each mechanism
        for mech_idx, mech_data in enumerate(all_mechanism_data):
            mech_name = mech_data['name']
            color = mech_data['color']
            states = mech_data['states']
            labels = mech_data['labels']
            barriers = mech_data['barriers']
            is_electrochemical = mech_data['is_electrochemical']
            n_electrons = mech_data['n_electrons']
            
            # Determine x positions based on x_axis_mode
            if x_axis_mode == 'electron':
                # Use n(H+ + e-) as x positions
                x_base = np.array(n_electrons, dtype=float)
                x_pos_mech = x_base.copy()
                
                # Find real (non-None) state indices for consecutive checking
                real_indices = [i for i, s in enumerate(states) if s is not None]
                
                # Track which states are in same-n transitions
                same_n_states = set()
                
                # Check consecutive real states
                for idx in range(len(real_indices) - 1):
                    i = real_indices[idx]
                    j = real_indices[idx + 1]
                    
                    if x_base[i] == x_base[j]:
                        # Consecutive states with same n
                        # First state (i) gets -same_n_offset
                        x_pos_mech[i] = x_base[i] - same_n_offset
                        # Second state (j) stays at original position
                        x_pos_mech[j] = x_base[j]
                        # Mark both as part of same-n transition
                        same_n_states.add(i)
                        same_n_states.add(j)
            else:  # x_axis_mode == 'step'
                # Use step order as x positions (0, 1, 2, 3, ...)
                x_pos_mech = np.arange(len(states), dtype=float)
                # No same-n transitions when using step order
                same_n_states = set()
            
            # Draw horizontal lines only for non-None states
            for i in range(len(states)):
                if states[i] is not None:
                    # Use smaller width for same-n transitions
                    width = bar_width_small if i in same_n_states else bar_width_large
                    x_start = x_pos_mech[i] - width/2
                    x_end = x_pos_mech[i] + width/2
                    ax.plot([x_start, x_end], [states[i], states[i]], 
                           color=color, linewidth=2, solid_capstyle='round', 
                           label=mech_name if i == 0 else '')
                    
                    # Store label info with actual x position and color
                    if labels[i]:
                        label_info.append((x_pos_mech[i], states[i], labels[i], color))
            
            # Connect states, skipping None positions
            i = 0
            while i < len(states) - 1:
                # Find next non-None state
                j = i + 1
                while j < len(states) and states[j] is None:
                    j += 1
                
                if j < len(states):
                    # Get widths for states i and j
                    width_i = bar_width_small if i in same_n_states else bar_width_large
                    width_j = bar_width_small if j in same_n_states else bar_width_large
                    
                    # Connect state i to state j
                    x_start = x_pos_mech[i] + width_i/2
                    x_end = x_pos_mech[j] - width_j/2
                    
                    # Use electrochemical info from step i
                    # electrochemical: solid, thermal: dashed
                    is_echem = is_electrochemical[i] if i < len(is_electrochemical) else False
                    linestyle = '-' if is_echem else '--'
                    
                    # Check if there's a barrier
                    if show_barriers and i < len(barriers) and barriers[i] is not None:
                        # Draw smooth barrier
                        x_barrier = (x_start + x_end) / 2
                        barrier_height = barriers[i]
                        
                        x_points = np.array([x_start, x_barrier, x_end])
                        y_points = np.array([states[i], barrier_height, states[j]])
                        
                        x_smooth = np.linspace(x_start, x_end, 50)
                        try:
                            spl = make_interp_spline(x_points, y_points, k=2)
                            y_smooth = spl(x_smooth)
                            ax.plot(x_smooth, y_smooth, color=color, linewidth=1.5, 
                                   linestyle=linestyle, alpha=0.8)
                        except:
                            ax.plot([x_start, x_barrier], [states[i], barrier_height],
                                   color=color, linewidth=1.5, linestyle=linestyle, alpha=0.8)
                            ax.plot([x_barrier, x_end], [barrier_height, states[j]],
                                   color=color, linewidth=1.5, linestyle=linestyle, alpha=0.8)
                    else:
                        # Direct connection
                        ax.plot([x_start, x_end], [states[i], states[j]],
                               color=color, linewidth=1.5, linestyle=linestyle, alpha=0.8)
                
                i = j if j < len(states) else len(states)
        
        # Add labels with actual x positions and mechanism colors
        if show_labels:
            y_min, y_max = ax.get_ylim()
            y_range = y_max - y_min
            y_offset_up = y_range * 0.015   # 1.5% above
            y_offset_down = -y_range * 0.02  # 2% below
            
            # Track shown labels to avoid duplicates (normalize labels)
            shown_labels = set()
            
            # Draw labels for each state
            for idx, (x_pos, energy, lbl, label_color) in enumerate(label_info):
                # Normalize label (remove gas phase products)
                normalized_lbl = self._normalize_label(lbl)
                
                # Only show if not already shown
                if normalized_lbl not in shown_labels:
                    shown_labels.add(normalized_lbl)
                    # Alternate label position based on index
                    y_offset = y_offset_up if idx % 2 == 0 else y_offset_down
                    ax.text(x_pos, energy + y_offset, normalized_lbl, 
                           ha='center', va='bottom' if idx % 2 == 0 else 'top',
                           fontsize=fontsize_label, rotation=0, color=label_color)  # Use mechanism color
        
        # Formatting
        ax.set_ylabel('G (eV)', fontsize=fontsize_axis)
        
        if x_axis_mode == 'electron':
            # Get all unique n_electron values for x-axis ticks
            all_n_electrons = set()
            for mech_data in all_mechanism_data:
                all_n_electrons.update([n for n in mech_data['n_electrons'] if n is not None])
            all_n_electrons = sorted(list(all_n_electrons))
            
            ax.set_xlabel('n(H⁺ + e⁻)', fontsize=fontsize_axis)
            if all_n_electrons:
                ax.set_xticks(all_n_electrons)
                ax.set_xticklabels([f'{int(n)}' for n in all_n_electrons])
        else:  # x_axis_mode == 'step'
            # Get maximum number of states for x-axis ticks
            max_states = max(len(mech_data['states']) for mech_data in all_mechanism_data)
            all_steps = list(range(max_states))
            
            ax.set_xlabel('Reaction Coordinate', fontsize=fontsize_axis)
            if all_steps:
                ax.set_xticks(all_steps)
                ax.set_xticklabels([f'{int(n)}' for n in all_steps])
        
        # Set tick label font sizes
        ax.tick_params(axis='both', which='major', labelsize=fontsize_ticks)
        
        # Set y-axis limits with margin to prevent label clipping
        all_energies = []
        for mech_data in all_mechanism_data:
            all_energies.extend([s for s in mech_data['states'] if s is not None])
        if all_energies:
            y_min_data = min(all_energies)
            y_max_data = max(all_energies)
            y_range = y_max_data - y_min_data
            
            if show_labels:
                # Calculate where labels actually appear
                y_offset_up = y_range * 0.015
                y_offset_down = -y_range * 0.02  # negative because it's below
                text_height = y_range * 0.03  # approximate text height
                
                # Find the actual min and max y positions including labels
                y_min_with_labels = y_min_data + y_offset_down - text_height/2
                y_max_with_labels = y_max_data + y_offset_up + text_height/2
                
                # Add base margin
                y_margin = y_range * y_margin_fraction
                ax.set_ylim(y_min_with_labels - y_margin, y_max_with_labels + y_margin)
            else:
                # No labels, just use data range with margin
                y_margin = y_range * y_margin_fraction
                ax.set_ylim(y_min_data - y_margin, y_max_data + y_margin)
    
    def compare_mechanisms(self,
                          mechanism_names: Optional[List[str]] = None,
                          mechanism_steps: Optional[Dict[str, List]] = None,
                          colors: Optional[List[str]] = None,
                          figsize: Optional[tuple] = None,
                          show_barriers: bool = True,
                          show_labels: bool = False,
                          bar_width_small: float = 0.3,
                          bar_width_large: float = 0.6,
                          same_n_offset: float = 0.5,
                          initial_label: Optional[str] = None,
                          legend_position: str = 'right',
                          width_per_n: float = 1.5,
                          base_width: float = 4.0,
                          height: float = 6.0,
                          fontsize_label: int = 11,
                          fontsize_axis: int = 14,
                          fontsize_title: int = 15,
                          fontsize_legend: int = 11,
                          fontsize_ticks: int = 12,
                          y_margin_fraction: float = 0.15,
                          show_voltage_text: bool = False,
                          x_axis_mode: str = 'step') -> plt.Figure:
        """
        Compare different reaction mechanisms
        
        Args:
            mechanism_names: List of mechanism names (if using predefined mechanisms)
            mechanism_steps: Dict of mechanism_name -> step_list (can include None/'x' to skip positions)
                           Example: {'C1_via_CO-H-ele': [5, None, 11, 12, 13, 8, 9, 10],
                                    'C1_via_H-CO': [5, 1, 14, 6, 7, 8, 9, 10]}
            colors: List of colors for each mechanism
            figsize: Figure size (if None, calculated dynamically)
            show_barriers: Whether to show activation barriers
            show_labels: Whether to show state labels
            bar_width_small: Width for bars in same-n transitions (default: 0.3)
            bar_width_large: Width for regular bars (default: 0.6)
            same_n_offset: X-offset for initial state in same-n transitions (default: 0.5)
            initial_label: Custom label for initial state
            legend_position: Legend position ('right', 'upper', 'best', etc.)
            width_per_n: Width per n(H+ + e-) unit (default: 1.5)
            base_width: Minimum figure width (default: 4.0)
            height: Figure height (default: 6.0)
            fontsize_label: Font size for state labels
            fontsize_axis: Font size for axis labels
            fontsize_title: Font size for title
            fontsize_legend: Font size for legend
            fontsize_ticks: Font size for tick labels
            y_margin_fraction: Fraction of y-range to add as margin
            show_voltage_text: Whether to show voltage text in plot
            x_axis_mode: X-axis mode - 'step' for reaction coordinate (default) or 'electron' for n(H+ + e-)
        
        Returns:
            Matplotlib figure object
        """
        # Determine which mechanisms to plot
        if mechanism_steps is not None:
            # User provided custom steps with None placeholders
            mech_names = list(mechanism_steps.keys())
            steps_to_use = mechanism_steps
        elif mechanism_names is not None:
            # Use predefined mechanisms
            mech_names = mechanism_names
            steps_to_use = {name: self.rxn_mechanisms[name] for name in mechanism_names}
        else:
            raise ValueError("Must provide either mechanism_names or mechanism_steps")
        
        if colors is None:
            colors = plt.cm.tab10(np.linspace(0, 1, len(mech_names)))
        
        # Calculate dynamic figure size if not provided
        if figsize is None:
            if x_axis_mode == 'electron':
                # Find maximum x range across all mechanisms
                max_x_range = 0
                for name in mech_names:
                    try:
                        data = self.calculate_mechanism_energies(name, verbose=False)
                        n_electrons = data['n_electrons']
                        x_range = max(n_electrons) - min(n_electrons)
                        max_x_range = max(max_x_range, x_range)
                    except:
                        pass
                dynamic_width = base_width + width_per_n * max_x_range
            else:  # x_axis_mode == 'step'
                # Find maximum number of steps across all mechanisms
                max_steps = 0
                for name in mech_names:
                    try:
                        data = self.calculate_mechanism_energies(name, verbose=False)
                        n_steps = len(data['states']) - 1
                        max_steps = max(max_steps, n_steps)
                    except:
                        pass
                dynamic_width = base_width + width_per_n * max_steps
            figsize = (dynamic_width, height)
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Use aligned plotting with user-controlled steps
        self._plot_aligned_mechanisms(ax, mech_names, steps_to_use, colors, 
                                     show_barriers, show_labels,
                                     bar_width_small, bar_width_large, same_n_offset,
                                     initial_label, fontsize_label, fontsize_axis,
                                     fontsize_ticks, y_margin_fraction, x_axis_mode)
        
        # Add horizontal line at y=0
        ax.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.3, zorder=0)
        
        # Add vertical dotted grid lines at 0.5, 1.5, 2.5, ...
        x_lim = ax.get_xlim()
        grid_positions = np.arange(0.5, x_lim[1] + 0.5, 1.0)
        for x_grid in grid_positions:
            if x_lim[0] <= x_grid <= x_lim[1]:
                ax.axvline(x=x_grid, color='gray', linestyle=':', linewidth=1, alpha=0.5, zorder=0)
        
        # Add voltage text inside the plot
        if show_voltage_text:
            voltage_text = f'U = {self.voltage:.2f} V vs. RHE'
            ax.text(0.02, 0.98, voltage_text, transform=ax.transAxes,
                   fontsize=fontsize_ticks, verticalalignment='top',
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='gray', alpha=0.8))
        
        # Add legend
        handles, labels_list = ax.get_legend_handles_labels()
        if handles:
            if legend_position == 'right':
                ax.legend(fontsize=fontsize_legend, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
            elif legend_position == 'upper':
                ax.legend(fontsize=fontsize_legend, loc='lower left', bbox_to_anchor=(0, 1.02), 
                         ncol=min(len(mech_names), 4), frameon=False)
            else:
                ax.legend(fontsize=fontsize_legend, loc=legend_position, frameon=False)
        ax.set_title(f'Mechanism Comparison at U={self.voltage:.2f} V vs RHE',
                    fontsize=fontsize_title, fontweight='bold')
        
        plt.tight_layout()
        return fig
    
    def save_mechanism_data(self, data: Dict, csv_path: str, x_axis_mode: str = 'step') -> None:
        """
        Save mechanism raw data to CSV file
        
        Args:
            data: Dictionary returned by calculate_mechanism_energies
            csv_path: Path to save CSV file
            x_axis_mode: X-axis mode - 'step' or 'electron' (default: 'step')
        """
        # Extract data
        mechanism_name = data['mechanism_name']
        states = data['states']
        labels = data['labels']
        delta_G_values = data.get('delta_G_values', [])
        Ga_values = data.get('Ga_values', [])
        n_electrons = data['n_electrons']
        is_electrochemical = data.get('is_electrochemical', [])
        step_indices = data['step_indices']
        
        # Calculate x coordinates based on mode
        n_states = len(states)
        if x_axis_mode == 'electron':
            x_coords = n_electrons
        else:  # x_axis_mode == 'step'
            x_coords = list(range(n_states))
        
        # Prepare data for CSV
        rows = []
        
        # Initial state
        rows.append({
            'Step': 0,
            'Label': labels[0],
            'X': x_coords[0],
            'G (eV)': round(states[0], 4),
            'ΔG (eV)': 0.0,
            'Ga (eV)': '',
            'Reaction_Step': 'Initial',
            'n(H+ + e-)': n_electrons[0],
            'Electrochemical': '',
            'Temperature (K)': data['temperature'],
            'Voltage (V)': data['voltage']
        })
        
        # Each reaction step
        for i in range(len(delta_G_values)):
            step_num = step_indices[i]
            rxn_data = self.rxn_expressions[step_num]
            
            rows.append({
                'Step': i + 1,
                'Label': labels[i + 1],
                'X': x_coords[i + 1],
                'G (eV)': round(states[i + 1], 4),
                'ΔG (eV)': round(delta_G_values[i], 4),
                'Ga (eV)': round(Ga_values[i], 4) if Ga_values[i] is not None else '',
                'Reaction_Step': rxn_data['equation'],
                'n(H+ + e-)': n_electrons[i + 1],
                'Electrochemical': 'Yes' if is_electrochemical[i] else 'No',
                'Temperature (K)': data['temperature'],
                'Voltage (V)': data['voltage']
            })
        
        # Create DataFrame and save
        df = pd.DataFrame(rows)
        df.to_csv(csv_path, index=False)
    
    def create_energy_table(self, 
                           mechanism_names: Optional[List[str]] = None,
                           save_to_file: Optional[str] = None) -> pd.DataFrame:
        """
        Create a summary table of energies for mechanisms
        
        Args:
            mechanism_names: List of mechanism names (None = all)
            save_to_file: Filename to save table (CSV format)
        
        Returns:
            DataFrame with energy values
        """
        if mechanism_names is None:
            mechanism_names = list(self.rxn_mechanisms.keys())
        
        table_data = []
        
        for mech_name in mechanism_names:
            try:
                data = self.calculate_mechanism_energies(mech_name, verbose=False)
                
                states = data['states']
                delta_G_values = data['delta_G_values']
                Ga_values = data['Ga_values']
                labels = data['labels'][1:]  # Skip 'Initial'
                is_echem = data['is_electrochemical']
                
                # Overall mechanism energetics
                overall_deltaG = states[-1] - states[0]
                max_barrier = max([g for g in Ga_values if g is not None], default=None)
                
                # Add overall row
                table_data.append({
                    'Mechanism': mech_name,
                    'Step': 'Overall',
                    'State': labels[-1] if labels else '',
                    'G (eV)': states[-1],
                    'ΔG (eV)': overall_deltaG,
                    'Ga (eV)': max_barrier,
                    'Electrochemical': 'Mixed' if any(is_echem) and not all(is_echem) else ('Yes' if all(is_echem) else 'No'),
                    'V (V)': self.voltage,
                    'T (K)': self.temperature
                })
                
                # Add individual steps
                for i, (lbl, dG, Ga, echem) in enumerate(zip(labels, delta_G_values, Ga_values, is_echem)):
                    table_data.append({
                        'Mechanism': mech_name,
                        'Step': f'Step {i+1}',
                        'State': lbl,
                        'G (eV)': states[i+1],
                        'ΔG (eV)': dG,
                        'Ga (eV)': Ga if Ga is not None else '-',
                        'Electrochemical': 'Yes' if echem else 'No',
                        'V (V)': self.voltage,
                        'T (K)': self.temperature
                    })
                
            except Exception as e:
                print(f"Error calculating {mech_name}: {e}")
                continue
        
        df = pd.DataFrame(table_data)
        
        # Save to file if requested
        if save_to_file:
            df.to_csv(save_to_file, index=False)
            print(f"Energy table saved to: {save_to_file}")
        
        return df


def compare_catalysts_mechanisms(
    plotters: Dict[str, 'FreeEnergyDiagramPlotter'],
    mechanisms: List[str],
    subplot_by: str = 'mechanism',
    figsize_per_plot: Optional[tuple] = None,
    show_barriers: bool = False,
    show_labels: bool = True,
    show_legend: bool = True,
    bar_width_small: float = 0.3,
    bar_width_large: float = 0.6,
    same_n_offset: float = 0.5,
    initial_label: Optional[str] = None,
    colors: Optional[Dict[str, str]] = None,
    legend_position: str = 'right',
    width_per_n: float = 1.5,
    base_width: float = 4.0,
    height: float = 6.0,
    fontsize_label: int = 11,
    fontsize_axis: int = 14,
    fontsize_title: int = 15,
    fontsize_legend: int = 11,
    fontsize_ticks: int = 12,
    y_margin_fraction: float = 0.15,
    show_voltage_text: bool = False,
    x_axis_mode: str = 'step') -> plt.Figure:
    """
    Compare multiple catalysts and mechanisms with flexible subplot arrangement
    
    Args:
        plotters: Dict of {catalyst_label: FreeEnergyDiagramPlotter}
                 e.g., {'Cu(100)': plotter1, 'Cu(211)': plotter2}
        mechanisms: List of mechanism names to compare
        subplot_by: 'mechanism' or 'catalyst'
                   'mechanism': each subplot shows one mechanism with multiple catalysts
                   'catalyst': each subplot shows one catalyst with multiple mechanisms
        figsize_per_plot: Size of each subplot (width, height) - if None, calculated dynamically
        show_barriers: Whether to show activation barriers
        show_labels: Whether to show state labels
        show_legend: Whether to show legend (default: True)
        bar_width_small: Width for bars in same-n transitions
        bar_width_large: Width for regular bars
        same_n_offset: X-offset for initial state in same-n transitions
        initial_label: Custom label for initial state
        colors: Dict mapping catalyst_label or mechanism to color
               If None, automatic colors will be assigned
        legend_position: Legend position ('right', 'upper', 'best', etc.)
        width_per_n: Width per n(H+ + e-) unit (default: 1.5)
        base_width: Minimum figure width (default: 4.0)
        height: Figure height (default: 6.0)
        fontsize_label: Font size for state labels
        fontsize_axis: Font size for axis labels
        fontsize_title: Font size for title
        fontsize_legend: Font size for legend
        fontsize_ticks: Font size for tick labels
        y_margin_fraction: Fraction of y-range to add as margin
        show_voltage_text: Whether to show voltage text in plot
        x_axis_mode: X-axis mode - 'step' for reaction coordinate (default) or 'electron' for n(H+ + e-)
    
    Returns:
        Matplotlib figure object with subplots
    
    Examples:
        # Compare Cu(100) and Cu(211) for CO2R_CO mechanism (single plot)
        plotters = {'Cu(100)': plotter1, 'Cu(211)': plotter2}
        fig = compare_catalysts_mechanisms(plotters, ['CO2R_CO'], subplot_by='mechanism')
        
        # Compare Cu(100) and Cu(211) for multiple mechanisms (subplots by mechanism)
        fig = compare_catalysts_mechanisms(
            plotters, 
            ['CO2R_CO', 'CO2R_HCOOH'],
            subplot_by='mechanism'
        )
        
        # Compare Cu(100) and Cu(211) for multiple mechanisms (subplots by catalyst)
        fig = compare_catalysts_mechanisms(
            plotters,
            ['CO2R_CO', 'CO2R_HCOOH'],
            subplot_by='catalyst'
        )
    """
    import matplotlib.pyplot as plt
    
    catalyst_labels = list(plotters.keys())
    n_catalysts = len(catalyst_labels)
    n_mechanisms = len(mechanisms)
    
    # Determine subplot layout
    if subplot_by == 'mechanism':
        # Each subplot = one mechanism, multiple catalysts
        n_subplots = n_mechanisms
        subplot_titles = mechanisms
    elif subplot_by == 'catalyst':
        # Each subplot = one catalyst, multiple mechanisms
        n_subplots = n_catalysts
        subplot_titles = catalyst_labels
    else:
        raise ValueError("subplot_by must be 'mechanism' or 'catalyst'")
    
    # Calculate dynamic figure size if not provided
    if figsize_per_plot is None:
        if x_axis_mode == 'electron':
            # Find maximum x range across all mechanisms and catalysts
            max_x_range = 0
            for plotter in plotters.values():
                for mech in mechanisms:
                    try:
                        data = plotter.calculate_mechanism_energies(mech, verbose=False)
                        n_electrons = data['n_electrons']
                        x_range = max(n_electrons) - min(n_electrons)
                        max_x_range = max(max_x_range, x_range)
                    except:
                        pass
            dynamic_width = base_width + width_per_n * max_x_range
        else:  # x_axis_mode == 'step'
            # Find maximum number of steps across all mechanisms and catalysts
            max_steps = 0
            for plotter in plotters.values():
                for mech in mechanisms:
                    try:
                        data = plotter.calculate_mechanism_energies(mech, verbose=False)
                        n_steps = len(data['states']) - 1
                        max_steps = max(max_steps, n_steps)
                    except:
                        pass
            dynamic_width = base_width + width_per_n * max_steps
        figsize_per_plot = (dynamic_width, height)
    
    # Create figure with subplots
    if n_subplots == 1:
        fig, axes = plt.subplots(1, 1, figsize=figsize_per_plot)
        axes = [axes]
    else:
        n_cols = min(2, n_subplots)
        n_rows = (n_subplots + n_cols - 1) // n_cols
        figsize = (figsize_per_plot[0] * n_cols, figsize_per_plot[1] * n_rows)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
        axes = axes.flatten() if n_subplots > 1 else [axes]
    
    # Assign colors
    if colors is None:
        if subplot_by == 'mechanism':
            # Color by catalyst
            color_map = dict(zip(catalyst_labels, 
                               plt.cm.tab10(np.linspace(0, 1, n_catalysts))))
        else:
            # Color by mechanism
            color_map = dict(zip(mechanisms,
                               plt.cm.tab10(np.linspace(0, 1, n_mechanisms))))
    else:
        color_map = colors
    
    # Track global y limits for consistent scaling
    all_y_values = []
    all_barrier_values = []
    
    # Plot each subplot
    for idx in range(n_subplots):
        ax = axes[idx]
        
        if subplot_by == 'mechanism':
            # One mechanism, multiple catalysts
            mech = mechanisms[idx]
            
            # Track shown labels per subplot to avoid duplicates within same plot
            subplot_shown_labels = set()
            
            for cat_idx, cat_label in enumerate(catalyst_labels):
                plotter = plotters[cat_label]
                color = color_map.get(cat_label, 'blue')
                
                # Show labels for first catalyst only, but track which labels are shown
                should_show_labels = show_labels and (cat_idx == 0)
                
                # Plot this catalyst's data
                plotter.plot_mechanism(
                    mech,
                    ax=ax,
                    color=color,
                    label=cat_label,
                    show_barriers=show_barriers,
                    show_labels=should_show_labels,
                    show_legend=False,  # Will add legend later
                    bar_width_small=bar_width_small,
                    bar_width_large=bar_width_large,
                    same_n_offset=same_n_offset,
                    initial_label=initial_label,
                    legend_position=None,
                    fontsize_label=fontsize_label,
                    fontsize_axis=fontsize_axis,
                    fontsize_title=fontsize_title,
                    fontsize_legend=fontsize_legend,
                    fontsize_ticks=fontsize_ticks,
                    y_margin_fraction=y_margin_fraction,
                    show_voltage_text=show_voltage_text,
                    x_axis_mode=x_axis_mode
                )
                
                # Collect y values for ylim
                data = plotter.calculate_mechanism_energies(mech)
                all_y_values.extend(data['states'])
                
                # Include barriers if showing them
                if show_barriers and data.get('barriers'):
                    valid_barriers = [b for b in data['barriers'] if b is not None]
                    all_barrier_values.extend(valid_barriers)
            
            ax.set_title(f'{mech}', fontsize=fontsize_title, fontweight='bold')
            
        else:  # subplot_by == 'catalyst'
            # One catalyst, multiple mechanisms
            cat_label = catalyst_labels[idx]
            plotter = plotters[cat_label]
            
            for mech_idx, mech in enumerate(mechanisms):
                color = color_map.get(mech, 'blue')
                
                # Show labels for all mechanisms (each gets its own color and should be distinguishable)
                # For single subplot with multiple mechanisms, show all labels
                should_show_labels = show_labels
                
                # Plot this mechanism's data
                plotter.plot_mechanism(
                    mech,
                    ax=ax,
                    color=color,
                    label=mech,
                    show_barriers=show_barriers,
                    show_labels=should_show_labels,
                    show_legend=False,  # Will add legend later
                    bar_width_small=bar_width_small,
                    bar_width_large=bar_width_large,
                    same_n_offset=same_n_offset,
                    initial_label=initial_label,
                    legend_position=None,
                    fontsize_label=fontsize_label,
                    fontsize_axis=fontsize_axis,
                    fontsize_title=fontsize_title,
                    fontsize_legend=fontsize_legend,
                    fontsize_ticks=fontsize_ticks,
                    y_margin_fraction=y_margin_fraction,
                    show_voltage_text=show_voltage_text,
                    x_axis_mode=x_axis_mode
                )
                
                # Collect y values for ylim
                data = plotter.calculate_mechanism_energies(mech)
                all_y_values.extend(data['states'])
                
                # Include barriers if showing them
                if show_barriers and data.get('barriers'):
                    valid_barriers = [b for b in data['barriers'] if b is not None]
                    all_barrier_values.extend(valid_barriers)
            
            ax.set_title(f'{cat_label}', fontsize=fontsize_title, fontweight='bold')
        
        # Add legend for this subplot if requested
        if show_legend:
            handles, labels_list = ax.get_legend_handles_labels()
            if handles:
                if legend_position == 'right':
                    ax.legend(fontsize=fontsize_legend, loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
                elif legend_position == 'upper':
                    ax.legend(fontsize=fontsize_legend, loc='lower left', bbox_to_anchor=(0, 1.02), 
                             ncol=min(len(handles), 4), frameon=False)
                else:
                    ax.legend(fontsize=fontsize_legend, loc=legend_position, frameon=False)
    
    # Set consistent y limits across all subplots
    if all_y_values:
        y_min = min(all_y_values)
        y_max = max(all_y_values)
        
        # Include barriers in max calculation
        if all_barrier_values:
            y_max = max(y_max, max(all_barrier_values))
        
        y_range = y_max - y_min
        
        if show_labels:
            # Calculate where labels actually appear
            y_offset_up = y_range * 0.015
            y_offset_down = -y_range * 0.02  # negative because it's below
            text_height = y_range * 0.03  # approximate text height
            
            # Find the actual min and max y positions including labels
            y_min_with_labels = y_min + y_offset_down - text_height/2
            y_max_with_labels = y_max + y_offset_up + text_height/2
            
            # Add base margin
            y_margin = y_range * 0.1
            for ax in axes[:n_subplots]:
                ax.set_ylim(y_min_with_labels - y_margin, y_max_with_labels + y_margin)
        else:
            # No labels, just use data range with margin
            y_margin = y_range * 0.1
            for ax in axes[:n_subplots]:
                ax.set_ylim(y_min - y_margin, y_max + y_margin)
    
    # Hide extra subplots if any
    for idx in range(n_subplots, len(axes)):
        axes[idx].axis('off')
    
    # Add overall title
    voltage = list(plotters.values())[0].voltage
    fig.suptitle(f'Catalyst & Mechanism Comparison at U={voltage:.2f} V vs RHE',
                fontsize=fontsize_title, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


if __name__ == "__main__":
    print("="*70)
    print("Free Energy Diagram Plotter - Test")
    print("="*70)
    
    # Initialize plotter
    plotter = FreeEnergyDiagramPlotter(
        mkm_file='COR_template.mkm',
        input_file='input.txt',
        temperature=298.15,
        voltage=-0.5
    )
    
    # List available mechanisms
    print("\nAvailable mechanisms:")
    for name in plotter.rxn_mechanisms.keys():
        steps = plotter.rxn_mechanisms[name]
        print(f"  {name}: steps {steps}")
    
    # Test 1: Plot single mechanism
    print("\n" + "="*70)
    print("Test 1: Plotting HER_Heyrovsky mechanism")
    print("="*70)
    
    fig1 = plotter.plot_mechanism(
        'HER_Heyrovsky',
        show_barriers=True,
        show_labels=True,
        verbose=True
    )
    fig1.savefig('FED_HER_Heyrovsky.png', dpi=150, bbox_inches='tight')
    print("\nSaved: FED_HER_Heyrovsky.png")
    
    # Test 2: Compare voltages
    print("\n" + "="*70)
    print("Test 2: Comparing HER at different voltages")
    print("="*70)
    
    voltages = [-1.0, -0.8, -0.6, -0.4, -0.2, 0.0]
    fig2 = plotter.compare_voltages(
        'HER_Heyrovsky',
        voltages=voltages,
        show_barriers=True,
        show_labels=False
    )
    fig2.savefig('FED_HER_voltage_comparison.png', dpi=150, bbox_inches='tight')
    print("\nSaved: FED_HER_voltage_comparison.png")
    
    # Test 3: Plot C1 mechanism
    print("\n" + "="*70)
    print("Test 3: Plotting C1_via_CO-H-ele mechanism")
    print("="*70)
    
    plotter.voltage = -0.5
    plotter.calculator.voltage = -0.5
    plotter.calculator.state_calculator.voltage = -0.5
    
    fig3 = plotter.plot_mechanism(
        'C1_via_CO-H-ele',
        show_barriers=True,
        show_labels=True,
        verbose=False
    )
    fig3.savefig('FED_C1_via_CO-H-ele.png', dpi=150, bbox_inches='tight')
    print("\nSaved: FED_C1_via_CO-H-ele.png")
    
    # Test 4: Compare different C1 mechanisms
    print("\n" + "="*70)
    print("Test 4: Comparing C1 mechanisms")
    print("="*70)
    
    fig4 = plotter.compare_mechanisms(
        ['C1_via_CO-H-ele', 'C1_via_H-CO'],
        show_barriers=True,
        show_labels=True  # Show labels in comparison
    )
    fig4.savefig('FED_C1_comparison.png', dpi=150, bbox_inches='tight')
    print("\nSaved: FED_C1_comparison.png")
    
    # Test 5: Create energy table
    print("\n" + "="*70)
    print("Test 5: Creating energy tables")
    print("="*70)
    
    # Table for HER
    df_her = plotter.create_energy_table(
        mechanism_names=['HER_Heyrovsky'],
        save_to_file='energy_table_HER.csv'
    )
    print("\nHER Mechanism Energy Table:")
    print(df_her.to_string(index=False))
    
    # Table for C1 mechanisms
    df_c1 = plotter.create_energy_table(
        mechanism_names=['C1_via_CO-H-ele', 'C1_via_H-CO'],
        save_to_file='energy_table_C1.csv'
    )
    print("\nC1 Mechanisms Energy Table:")
    print(df_c1[['Mechanism', 'Step', 'State', 'ΔG (eV)', 'Ga (eV)', 'Electrochemical']].to_string(index=False))
    
    # Table for all mechanisms
    df_all = plotter.create_energy_table(save_to_file='energy_table_all.csv')
    print(f"\nAll mechanisms table saved to: energy_table_all.csv")
    print(f"Total mechanisms: {len(plotter.rxn_mechanisms)}")
    
    print("\n" + "="*70)
    print("All tests completed!")
    print("="*70)

