"""
Command-line interface for cathub-plotter
"""
import argparse
import sys
import os
import tempfile
import pandas as pd
from pathlib import Path
from typing import List, Optional

from ..core.calculator import FreeEnergyCalculator
from ..plotters.diagram import FreeEnergyDiagramPlotter


def filter_input_file(input_file: str, surface: Optional[str] = None, facet: Optional[str] = None) -> str:
    """
    Filter input file based on surface_name and site_name
    
    Args:
        input_file: Path to original input file
        surface: Surface name to filter (e.g., 'Cu')
        facet: Facet name to filter (e.g., '100')
    
    Returns:
        Path to filtered input file (temporary file)
    """
    # Determine file type and read accordingly
    _, ext = os.path.splitext(input_file)
    
    if ext.lower() in ['.xlsx', '.xls']:
        df = pd.read_excel(input_file)
    elif ext.lower() == '.csv':
        df = pd.read_csv(input_file)
    else:
        df = pd.read_csv(input_file, sep='\t')
    
    # Filter based on surface_name and site_name
    # Always keep gas phase and slab species
    gas_slab_mask = df['status'].isin(['gas', 'slab'])
    
    # Build filter for adsorbed species and transition states
    ads_ts_mask = df['status'].isin(['ads', 'ts'])
    
    if surface:
        ads_ts_mask = ads_ts_mask & (df['surface_name'] == surface)
    
    if facet:
        # Convert facet to string for comparison (handles both '100' and 100)
        facet_str = str(facet)
        ads_ts_mask = ads_ts_mask & (df['site_name'].astype(str) == facet_str)
    
    # Combine: keep gas/slab OR matching ads/ts
    mask = gas_slab_mask | ads_ts_mask
    
    # Apply filter
    filtered_df = df[mask].copy()
    
    # Create temporary file with same extension
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix=ext, delete=False)
    temp_path = temp_file.name
    temp_file.close()
    
    # Save filtered data
    if ext.lower() in ['.xlsx', '.xls']:
        filtered_df.to_excel(temp_path, index=False)
    elif ext.lower() == '.csv':
        filtered_df.to_csv(temp_path, index=False)
    else:
        filtered_df.to_csv(temp_path, sep='\t', index=False)
    
    return temp_path


def plot_command(args):
    """Plot free energy diagram from MKM file"""
    try:
        # Check if surface/facet filtering is requested
        if args.surface or args.facet:
            # Filter input file based on surface_name and site_name
            input_file = filter_input_file(
                args.input_file, 
                surface=args.surface, 
                facet=args.facet
            )
            print(f"Filtering species for surface={args.surface or 'all'}, facet={args.facet or 'all'}")
        else:
            input_file = args.input_file
        
        plotter = FreeEnergyDiagramPlotter(
            mkm_file=args.mkm_file,
            input_file=input_file,
            temperature=args.temperature,
            voltage=args.voltage,
            use_gibbs_energy=args.use_gibbs if hasattr(args, 'use_gibbs') else False
        )
        
        # Prepare plot_mechanism kwargs
        plot_kwargs = {}
        if hasattr(args, 'y_margin') and args.y_margin is not None:
            plot_kwargs['y_margin_fraction'] = args.y_margin
        if hasattr(args, 'ylim') and args.ylim is not None:
            plot_kwargs['ylim'] = args.ylim
        if hasattr(args, 'x_axis') and args.x_axis is not None:
            plot_kwargs['x_axis_mode'] = args.x_axis
        if hasattr(args, 'figsize') and args.figsize is not None:
            plot_kwargs['figsize'] = tuple(args.figsize)
        if hasattr(args, 'bar_width') and args.bar_width is not None:
            plot_kwargs['bar_width_small'] = args.bar_width
            plot_kwargs['bar_width_large'] = args.bar_width
        
        if args.mechanism:
            plotter.plot_mechanism(args.mechanism, save_path=args.output, **plot_kwargs)
            
            # Save raw data if requested
            if args.save_data:
                # Determine CSV filename
                if args.output:
                    csv_path = args.output.rsplit('.', 1)[0] + '_data.csv'
                else:
                    csv_path = f"{args.mechanism}_data.csv"
                
                # Calculate and save data
                data = plotter.calculate_mechanism_energies(args.mechanism, verbose=False)
                x_axis_mode = args.x_axis if hasattr(args, 'x_axis') else 'step'
                plotter.save_mechanism_data(data, csv_path, x_axis_mode=x_axis_mode)
                print(f"Raw data saved to {csv_path}")
        else:
            # Plot all mechanisms
            mechanisms = list(plotter.rxn_mechanisms.keys())
            for mech in mechanisms:
                output_path = f"{args.output}_{mech}.png" if args.output else None
                plotter.plot_mechanism(mech, save_path=output_path, **plot_kwargs)
                
                # Save raw data if requested
                if args.save_data:
                    if args.output:
                        csv_path = f"{args.output}_{mech}_data.csv"
                    else:
                        csv_path = f"{mech}_data.csv"
                    
                    data = plotter.calculate_mechanism_energies(mech, verbose=False)
                    x_axis_mode = args.x_axis if hasattr(args, 'x_axis') else 'step'
                    plotter.save_mechanism_data(data, csv_path, x_axis_mode=x_axis_mode)
                    print(f"Raw data saved to {csv_path}")
        
        print("Plotting completed successfully!")
        
    except Exception as e:
        print(f"Error plotting diagram: {e}")
        sys.exit(1)


def compare_command(args):
    """Compare free energy diagrams across different conditions"""
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        from ..plotters.diagram import compare_catalysts_mechanisms
        
        # Parse what to vary
        vary_params = [v.strip() for v in args.vary.split(',')]
        
        if len(vary_params) == 1:
            # 1D comparison
            _compare_1d(args, vary_params[0])
        elif len(vary_params) == 2:
            # 2D comparison
            _compare_2d(args, vary_params[0], vary_params[1])
        else:
            print("Error: Can only vary 1 or 2 parameters at a time")
            sys.exit(1)
        
        print("Comparison completed successfully!")
        
    except Exception as e:
        print(f"Error creating comparison: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def _save_comparison_data(args, vary_param: str, values: list, plotter_data: list):
    """
    Save comparison data to CSV
    
    Args:
        args: Command line arguments
        vary_param: Parameter being varied
        values: List of values for the varying parameter
        plotter_data: List of tuples (value, plotter, mechanism)
    """
    import pandas as pd
    
    # Determine CSV filename
    if args.output:
        csv_path = args.output.rsplit('.', 1)[0] + '_data.csv'
    else:
        csv_path = f"comparison_{vary_param}_data.csv"
    
    # Collect all data
    all_rows = []
    
    for value, plotter, mechanism in plotter_data:
        # Calculate energies
        data = plotter.calculate_mechanism_energies(mechanism, verbose=False)
        
        # Create label for this condition
        label = _create_label(vary_param, value, args)
        
        # Extract data
        states = data['states']
        labels_list = data['labels']
        delta_G_values = data.get('delta_G_values', [])
        Ga_values = data.get('Ga_values', [])
        n_electrons = data['n_electrons']
        is_electrochemical = data.get('is_electrochemical', [])
        step_indices = data['step_indices']
        
        # Initial state
        all_rows.append({
            'Condition': label,
            'Step': 0,
            'Reaction_Step': 'Initial',
            'State': labels_list[0],
            'n(H+ + e-)': n_electrons[0],
            'G (eV)': round(states[0], 4),
            'ΔG (eV)': 0.0,
            'Ga (eV)': '',
            'Electrochemical': '',
            'Temperature (K)': data['temperature'],
            'Voltage (V)': data['voltage'],
            'Mechanism': mechanism
        })
        
        # Each reaction step
        for i in range(len(delta_G_values)):
            step_num = step_indices[i]
            rxn_data = plotter.rxn_expressions[step_num]
            
            all_rows.append({
                'Condition': label,
                'Step': i + 1,
                'Reaction_Step': rxn_data['equation'],
                'State': labels_list[i + 1],
                'n(H+ + e-)': n_electrons[i + 1],
                'G (eV)': round(states[i + 1], 4),
                'ΔG (eV)': round(delta_G_values[i], 4),
                'Ga (eV)': round(Ga_values[i], 4) if Ga_values[i] is not None else '',
                'Electrochemical': 'Yes' if is_electrochemical[i] else 'No',
                'Temperature (K)': data['temperature'],
                'Voltage (V)': data['voltage'],
                'Mechanism': mechanism
            })
    
    # Create DataFrame and save
    df = pd.DataFrame(all_rows)
    df.to_csv(csv_path, index=False)
    print(f"Comparison data saved to {csv_path}")


def _save_comparison_data_2d(args, vary_param1: str, values1: list, vary_param2: str, values2: list, plotter_data: list):
    """
    Save 2D comparison data to CSV
    
    Args:
        args: Command line arguments
        vary_param1: First parameter being varied
        values1: List of values for first parameter
        vary_param2: Second parameter being varied
        values2: List of values for second parameter
        plotter_data: List of tuples (val1, val2, plotter, mechanism)
    """
    import pandas as pd
    
    # Determine CSV filename
    if args.output:
        csv_path = args.output.rsplit('.', 1)[0] + '_data.csv'
    else:
        csv_path = f"comparison_{vary_param1}_vs_{vary_param2}_data.csv"
    
    # Collect all data
    all_rows = []
    
    for val1, val2, plotter, mechanism in plotter_data:
        # Calculate energies
        data = plotter.calculate_mechanism_energies(mechanism, verbose=False)
        
        # Create labels for both conditions
        label1 = _create_label(vary_param1, val1, args)
        label2 = _create_label(vary_param2, val2, args)
        combined_label = f"{label1}, {label2}"
        
        # Extract data
        states = data['states']
        labels_list = data['labels']
        delta_G_values = data.get('delta_G_values', [])
        Ga_values = data.get('Ga_values', [])
        n_electrons = data['n_electrons']
        is_electrochemical = data.get('is_electrochemical', [])
        step_indices = data['step_indices']
        
        # Initial state
        all_rows.append({
            'Condition': combined_label,
            f'{vary_param1}': val1,
            f'{vary_param2}': val2,
            'Step': 0,
            'Reaction_Step': 'Initial',
            'State': labels_list[0],
            'n(H+ + e-)': n_electrons[0],
            'G (eV)': round(states[0], 4),
            'ΔG (eV)': 0.0,
            'Ga (eV)': '',
            'Electrochemical': '',
            'Temperature (K)': data['temperature'],
            'Voltage (V)': data['voltage'],
            'Mechanism': mechanism
        })
        
        # Each reaction step
        for i in range(len(delta_G_values)):
            step_num = step_indices[i]
            rxn_data = plotter.rxn_expressions[step_num]
            
            all_rows.append({
                'Condition': combined_label,
                f'{vary_param1}': val1,
                f'{vary_param2}': val2,
                'Step': i + 1,
                'Reaction_Step': rxn_data['equation'],
                'State': labels_list[i + 1],
                'n(H+ + e-)': n_electrons[i + 1],
                'G (eV)': round(states[i + 1], 4),
                'ΔG (eV)': round(delta_G_values[i], 4),
                'Ga (eV)': round(Ga_values[i], 4) if Ga_values[i] is not None else '',
                'Electrochemical': 'Yes' if is_electrochemical[i] else 'No',
                'Temperature (K)': data['temperature'],
                'Voltage (V)': data['voltage'],
                'Mechanism': mechanism
            })
    
    # Create DataFrame and save
    df = pd.DataFrame(all_rows)
    df.to_csv(csv_path, index=False)
    print(f"2D comparison data saved to {csv_path}")


def _compare_1d(args, vary_param: str):
    """Compare along a single parameter"""
    import matplotlib.pyplot as plt
    
    # Get the values to compare
    if vary_param == 'mechanism':
        if not args.mechanisms:
            print("Error: --mechanisms required when varying mechanism")
            sys.exit(1)
        values = [m.strip() for m in args.mechanisms.split(',')]
        base_mechanism = None
    elif vary_param == 'temperature':
        if not args.temperatures:
            print("Error: --temperatures required when varying temperature")
            sys.exit(1)
        values = [float(t.strip()) for t in args.temperatures.split(',')]
        base_mechanism = args.mechanism
    elif vary_param == 'voltage':
        if not args.voltages:
            print("Error: --voltages required when varying voltage")
            sys.exit(1)
        values = [float(v.strip()) for v in args.voltages.split(',')]
        base_mechanism = args.mechanism
    elif vary_param == 'surface':
        if not args.surfaces:
            print("Error: --surfaces required when varying surface")
            sys.exit(1)
        values = [s.strip() for s in args.surfaces.split(',')]
        base_mechanism = args.mechanism
    elif vary_param == 'facet':
        if not args.facets:
            print("Error: --facets required when varying facet")
            sys.exit(1)
        values = [f.strip() for f in args.facets.split(',')]
        base_mechanism = args.mechanism
    else:
        print(f"Error: Unknown parameter to vary: {vary_param}")
        sys.exit(1)
    
    # Parse colors if provided
    colors = None
    if args.colors:
        colors = [c.strip() for c in args.colors.split(',')]
        if len(colors) != len(values):
            print(f"Warning: Number of colors ({len(colors)}) doesn't match number of values ({len(values)})")
            colors = None
    
    # Parse labels if provided
    labels = None
    if args.labels:
        labels = [l.strip() for l in args.labels.split(',')]
        if len(labels) != len(values):
            print(f"Warning: Number of labels ({len(labels)}) doesn't match number of values ({len(values)})")
            labels = None
    
    # Create plotters for each case
    plotter_data = []  # For saving data: list of (value, plotter, mechanism)
    
    if vary_param == 'mechanism':
        # Varying mechanism: single plotter, compare mechanisms
        input_file = _get_filtered_input_file(args)
        plotter = FreeEnergyDiagramPlotter(
            mkm_file=args.mkm_file,
            input_file=input_file,
            temperature=args.temperature,
            voltage=args.voltage,
            use_gibbs_energy=args.use_gibbs if hasattr(args, 'use_gibbs') else False
        )
        
        # Collect data for saving
        for mech in values:
            plotter_data.append((mech, plotter, mech))
        
        # Prepare figsize
        figsize = tuple(args.figsize) if hasattr(args, 'figsize') and args.figsize else None
        
        # Prepare bar_width kwargs
        bar_width_kwargs = {}
        if hasattr(args, 'bar_width') and args.bar_width is not None:
            bar_width_kwargs['bar_width_small'] = args.bar_width
            bar_width_kwargs['bar_width_large'] = args.bar_width
        
        fig = plotter.compare_mechanisms(
            mechanism_names=values,
            colors=colors,
            figsize=figsize,
            show_barriers=args.show_barriers,
            show_labels=args.show_labels,
            legend_position=args.legend_position,
            x_axis_mode=args.x_axis if hasattr(args, 'x_axis') else 'step',
            **bar_width_kwargs
        )
    
    elif vary_param == 'voltage':
        # Varying voltage: single plotter, compare voltages
        input_file = _get_filtered_input_file(args)
        
        # Collect data for saving (need separate plotters for each voltage)
        for voltage in values:
            plotter_v = FreeEnergyDiagramPlotter(
                mkm_file=args.mkm_file,
                input_file=input_file,
                temperature=args.temperature,
                voltage=voltage,
                use_gibbs_energy=args.use_gibbs if hasattr(args, 'use_gibbs') else False
            )
            plotter_data.append((voltage, plotter_v, base_mechanism))
        
        # Create main plotter for plotting
        plotter = FreeEnergyDiagramPlotter(
            mkm_file=args.mkm_file,
            input_file=input_file,
            temperature=args.temperature,
            voltage=values[0],  # Will be updated in compare_voltages
            use_gibbs_energy=args.use_gibbs if hasattr(args, 'use_gibbs') else False
        )
        
        # Prepare figsize
        figsize = tuple(args.figsize) if hasattr(args, 'figsize') and args.figsize else None
        
        # Prepare bar_width kwargs
        bar_width_kwargs = {}
        if hasattr(args, 'bar_width') and args.bar_width is not None:
            bar_width_kwargs['bar_width_small'] = args.bar_width
            bar_width_kwargs['bar_width_large'] = args.bar_width
        
        fig = plotter.compare_voltages(
            mechanism_name=base_mechanism,
            voltages=values,
            colors=colors,
            figsize=figsize,
            show_barriers=args.show_barriers,
            show_labels=args.show_labels,
            legend_position=args.legend_position,
            x_axis_mode=args.x_axis if hasattr(args, 'x_axis') else 'step',
            **bar_width_kwargs
        )
    
    elif vary_param in ['temperature', 'surface', 'facet']:
        # Varying temperature, surface, or facet: create multiple plotters
        plotters = {}
        
        for idx, value in enumerate(values):
            # Determine label
            if labels and idx < len(labels):
                label = labels[idx]
            elif vary_param == 'temperature':
                label = f'T={value:.1f} K'
            elif vary_param == 'surface':
                label = f'{value}({args.facet})' if args.facet else value
            elif vary_param == 'facet':
                label = f'{args.surface}({value})' if args.surface else value
            
            # Get filtered input file
            if vary_param == 'surface':
                input_file = filter_input_file(args.input_file, surface=value, facet=args.facet)
            elif vary_param == 'facet':
                input_file = filter_input_file(args.input_file, surface=args.surface, facet=value)
            else:
                input_file = _get_filtered_input_file(args)
            
            # Create plotter
            temp = value if vary_param == 'temperature' else args.temperature
            plotter = FreeEnergyDiagramPlotter(
                mkm_file=args.mkm_file,
                input_file=input_file,
                temperature=temp,
                voltage=args.voltage,
                use_gibbs_energy=args.use_gibbs if hasattr(args, 'use_gibbs') else False
            )
            plotters[label] = plotter
            
            # Collect data for saving
            plotter_data.append((value, plotter, base_mechanism))
        
        # Use compare_catalysts_mechanisms for plotting
        from ..plotters.diagram import compare_catalysts_mechanisms
        
        # Prepare colors
        color_map = None
        if colors:
            color_map = {label: colors[idx] for idx, label in enumerate(plotters.keys()) if idx < len(colors)}
        
        # Prepare bar_width kwargs
        bar_width_kwargs = {}
        if hasattr(args, 'bar_width') and args.bar_width is not None:
            bar_width_kwargs['bar_width_small'] = args.bar_width
            bar_width_kwargs['bar_width_large'] = args.bar_width
        
        fig = compare_catalysts_mechanisms(
            plotters=plotters,
            mechanisms=[base_mechanism],
            subplot_by='catalyst',
            colors=color_map,
            show_barriers=args.show_barriers,
            show_labels=args.show_labels,
            show_legend=True,
            legend_position=args.legend_position,
            x_axis_mode=args.x_axis if hasattr(args, 'x_axis') else 'step',
            **bar_width_kwargs
        )
    
    # Save figure
    if args.output:
        fig.savefig(args.output, dpi=300, bbox_inches='tight')
        print(f"Comparison plot saved to {args.output}")
    else:
        plt.show()
    
    # Save data if requested
    if args.save_data and plotter_data:
        _save_comparison_data(args, vary_param, values, plotter_data)


def _compare_2d(args, vary_param1: str, vary_param2: str):
    """Compare along two parameters (2D comparison with subplots or overlay)"""
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Get values for both parameters
    values1 = _get_vary_values(args, vary_param1)
    values2 = _get_vary_values(args, vary_param2)
    
    print(f"Creating 2D comparison: {vary_param1} vs {vary_param2}")
    print(f"  {vary_param1}: {values1}")
    print(f"  {vary_param2}: {values2}")
    
    layout = args.layout if hasattr(args, 'layout') and args.layout else 'subplots'
    
    if layout == 'subplots':
        _compare_2d_subplots(args, vary_param1, values1, vary_param2, values2)
    else:  # overlay
        _compare_2d_overlay(args, vary_param1, values1, vary_param2, values2)


def _compare_2d_subplots(args, vary_param1: str, values1: list, vary_param2: str, values2: list):
    """2D comparison with subplots layout"""
    import matplotlib.pyplot as plt
    from ..plotters.diagram import compare_catalysts_mechanisms
    
    # Determine which parameter goes to subplots and which to overlay
    # Strategy: param1 -> subplots (rows/cols), param2 -> overlay (colors in each subplot)
    
    n_subplots = len(values1)
    n_cols = min(2, n_subplots)
    n_rows = (n_subplots + n_cols - 1) // n_cols
    
    # Determine figsize
    if hasattr(args, 'figsize') and args.figsize:
        figsize = tuple(args.figsize)
    else:
        figsize = (8*n_cols, 6*n_rows)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_subplots == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if n_subplots > 1 else [axes]
    
    # Generate colors for param2
    cmap = plt.cm.tab10
    colors2 = [cmap(i % 10) for i in range(len(values2))]
    
    # Collect all y-values across all subplots for consistent ylim
    all_y_values = []
    all_barrier_values = []
    
    # Collect plotter data for saving
    plotter_data = []  # List of (val1, val2, plotter, mechanism)
    
    # For each value of param1, create a subplot
    for idx1, val1 in enumerate(values1):
        ax = axes[idx1]
        subplot_y_values = []
        subplot_barrier_values = []
        
        # For each value of param2, plot on the same subplot
        for idx2, val2 in enumerate(values2):
            # Create plotter with specific conditions
            plotter = _create_plotter_for_conditions(
                args, vary_param1, val1, vary_param2, val2
            )
            
            # Get mechanism to plot
            mechanism = _get_mechanism_for_conditions(args, vary_param1, val1, vary_param2, val2)
            
            # Collect plotter data for saving
            plotter_data.append((val1, val2, plotter, mechanism))
            
            # Calculate energies to get y-values
            data = plotter.calculate_mechanism_energies(mechanism, verbose=False)
            subplot_y_values.extend(data['states'])
            all_y_values.extend(data['states'])
            
            # Include barriers if showing them
            if args.show_barriers and data.get('barriers'):
                valid_barriers = [b for b in data['barriers'] if b is not None]
                subplot_barrier_values.extend(valid_barriers)
                all_barrier_values.extend(valid_barriers)
            
            # Create label
            label = _create_label(vary_param2, val2, args)
            
            # Plot on this axis
            plotter.plot_mechanism(
                mechanism,
                ax=ax,
                color=colors2[idx2],
                label=label,
                show_barriers=args.show_barriers,
                show_labels=(idx2 == 0) and args.show_labels,  # Only first overlay shows labels
                show_legend=False,
                legend_position=None
            )
        
        # Set ylim for this subplot based on its data
        if subplot_y_values:
            y_min = min(subplot_y_values)
            y_max = max(subplot_y_values)
            
            # Include barriers in max calculation
            if subplot_barrier_values:
                y_max = max(y_max, max(subplot_barrier_values))
            
            y_range = y_max - y_min
            y_margin = y_range * 0.15  # 15% margin
            
            # Add extra margin for labels if showing them
            if args.show_labels:
                y_offset_up = y_range * 0.015
                y_offset_down = -y_range * 0.02
                text_height = y_range * 0.03
                y_min_with_labels = y_min + y_offset_down - text_height/2
                y_max_with_labels = y_max + y_offset_up + text_height/2
                ax.set_ylim(y_min_with_labels - y_margin, y_max_with_labels + y_margin)
            else:
                ax.set_ylim(y_min - y_margin, y_max + y_margin)
        
        # Set subplot title
        title = _create_label(vary_param1, val1, args)
        ax.set_title(title, fontsize=14, fontweight='bold')
        
        # Add legend to each subplot
        ax.legend(loc='best', fontsize=10, frameon=False)
    
    # Hide extra subplots
    for idx in range(n_subplots, len(axes)):
        axes[idx].axis('off')
    
    # Overall title
    fig.suptitle(f'Comparison: {vary_param1} vs {vary_param2}', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # Save figure
    if args.output:
        fig.savefig(args.output, dpi=300, bbox_inches='tight')
        print(f"2D comparison plot saved to {args.output}")
    else:
        plt.show()
    
    # Save data if requested
    if args.save_data and plotter_data:
        _save_comparison_data_2d(args, vary_param1, values1, vary_param2, values2, plotter_data)


def _compare_2d_overlay(args, vary_param1: str, values1: list, vary_param2: str, values2: list):
    """2D comparison with overlay layout (all on one plot)"""
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Determine figsize
    if hasattr(args, 'figsize') and args.figsize:
        figsize = tuple(args.figsize)
    else:
        figsize = (10, 8)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Generate colors: use 2D colormap or distinct colors
    n_total = len(values1) * len(values2)
    colors = plt.cm.tab20(np.linspace(0, 1, n_total))
    
    # Collect all y-values for ylim calculation
    all_y_values = []
    all_barrier_values = []
    
    # Collect plotter data for saving
    plotter_data = []  # List of (val1, val2, plotter, mechanism)
    
    color_idx = 0
    for idx1, val1 in enumerate(values1):
        for idx2, val2 in enumerate(values2):
            # Create plotter
            plotter = _create_plotter_for_conditions(
                args, vary_param1, val1, vary_param2, val2
            )
            
            # Get mechanism
            mechanism = _get_mechanism_for_conditions(args, vary_param1, val1, vary_param2, val2)
            
            # Collect plotter data for saving
            plotter_data.append((val1, val2, plotter, mechanism))
            
            # Calculate energies to get y-values
            data = plotter.calculate_mechanism_energies(mechanism, verbose=False)
            all_y_values.extend(data['states'])
            
            # Include barriers if showing them
            if args.show_barriers and data.get('barriers'):
                valid_barriers = [b for b in data['barriers'] if b is not None]
                all_barrier_values.extend(valid_barriers)
            
            # Create combined label
            label1 = _create_label(vary_param1, val1, args)
            label2 = _create_label(vary_param2, val2, args)
            label = f"{label1}, {label2}"
            
            # Plot
            plotter.plot_mechanism(
                mechanism,
                ax=ax,
                color=colors[color_idx],
                label=label,
                show_barriers=args.show_barriers,
                show_labels=(color_idx == 0) and args.show_labels,
                show_legend=False,
                legend_position=None,
                line_alpha=0.7
            )
            
            color_idx += 1
    
    # Set ylim based on all data
    if all_y_values:
        y_min = min(all_y_values)
        y_max = max(all_y_values)
        
        # Include barriers in max calculation
        if all_barrier_values:
            y_max = max(y_max, max(all_barrier_values))
        
        y_range = y_max - y_min
        y_margin = y_range * 0.15  # 15% margin
        
        # Add extra margin for labels if showing them
        if args.show_labels:
            y_offset_up = y_range * 0.015
            y_offset_down = -y_range * 0.02
            text_height = y_range * 0.03
            y_min_with_labels = y_min + y_offset_down - text_height/2
            y_max_with_labels = y_max + y_offset_up + text_height/2
            ax.set_ylim(y_min_with_labels - y_margin, y_max_with_labels + y_margin)
        else:
            ax.set_ylim(y_min - y_margin, y_max + y_margin)
    
    # Add legend
    ax.legend(loc='best', fontsize=9, frameon=False, ncol=2)
    ax.set_title(f'Comparison: {vary_param1} vs {vary_param2}', fontsize=15, fontweight='bold')
    
    plt.tight_layout()
    
    # Save figure
    if args.output:
        fig.savefig(args.output, dpi=300, bbox_inches='tight')
        print(f"2D overlay comparison plot saved to {args.output}")
    else:
        plt.show()
    
    # Save data if requested
    if args.save_data and plotter_data:
        _save_comparison_data_2d(args, vary_param1, values1, vary_param2, values2, plotter_data)


def _get_vary_values(args, vary_param: str) -> list:
    """Get values for a parameter to vary"""
    if vary_param == 'mechanism':
        if not args.mechanisms:
            print(f"Error: --mechanisms required when varying mechanism")
            sys.exit(1)
        return [m.strip() for m in args.mechanisms.split(',')]
    elif vary_param == 'temperature':
        if not args.temperatures:
            print(f"Error: --temperatures required when varying temperature")
            sys.exit(1)
        return [float(t.strip()) for t in args.temperatures.split(',')]
    elif vary_param == 'voltage':
        if not args.voltages:
            print(f"Error: --voltages required when varying voltage")
            sys.exit(1)
        return [float(v.strip()) for v in args.voltages.split(',')]
    elif vary_param == 'surface':
        if not args.surfaces:
            print(f"Error: --surfaces required when varying surface")
            sys.exit(1)
        return [s.strip() for s in args.surfaces.split(',')]
    elif vary_param == 'facet':
        if not args.facets:
            print(f"Error: --facets required when varying facet")
            sys.exit(1)
        return [f.strip() for f in args.facets.split(',')]
    else:
        print(f"Error: Unknown parameter to vary: {vary_param}")
        sys.exit(1)


def _create_plotter_for_conditions(args, vary_param1: str, val1, vary_param2: str, val2):
    """Create a plotter with specific conditions"""
    # Determine conditions
    conditions = {
        'temperature': args.temperature,
        'voltage': args.voltage,
        'surface': args.surface,
        'facet': args.facet
    }
    
    # Update with varied parameters
    if vary_param1 == 'temperature':
        conditions['temperature'] = val1
    elif vary_param1 == 'voltage':
        conditions['voltage'] = val1
    elif vary_param1 == 'surface':
        conditions['surface'] = val1
    elif vary_param1 == 'facet':
        conditions['facet'] = val1
    
    if vary_param2 == 'temperature':
        conditions['temperature'] = val2
    elif vary_param2 == 'voltage':
        conditions['voltage'] = val2
    elif vary_param2 == 'surface':
        conditions['surface'] = val2
    elif vary_param2 == 'facet':
        conditions['facet'] = val2
    
    # Get filtered input file
    input_file = filter_input_file(
        args.input_file,
        surface=conditions['surface'],
        facet=conditions['facet']
    )
    
    # Create plotter
    plotter = FreeEnergyDiagramPlotter(
        mkm_file=args.mkm_file,
        input_file=input_file,
        temperature=conditions['temperature'],
        voltage=conditions['voltage'],
        use_gibbs_energy=args.use_gibbs if hasattr(args, 'use_gibbs') else False
    )
    
    return plotter


def _get_mechanism_for_conditions(args, vary_param1: str, val1, vary_param2: str, val2) -> str:
    """Get mechanism name for specific conditions"""
    if vary_param1 == 'mechanism':
        return val1
    elif vary_param2 == 'mechanism':
        return val2
    else:
        return args.mechanism


def _create_label(vary_param: str, value, args) -> str:
    """Create a label for a specific parameter value"""
    if vary_param == 'mechanism':
        return str(value)
    elif vary_param == 'temperature':
        return f'T={value:.1f}K'
    elif vary_param == 'voltage':
        return f'U={value:.2f}V'
    elif vary_param == 'surface':
        facet = args.facet if hasattr(args, 'facet') and args.facet else ''
        return f'{value}({facet})' if facet else value
    elif vary_param == 'facet':
        surface = args.surface if hasattr(args, 'surface') and args.surface else ''
        return f'{surface}({value})' if surface else value
    else:
        return str(value)


def _get_filtered_input_file(args) -> str:
    """Get filtered input file based on surface/facet"""
    if args.surface or args.facet:
        return filter_input_file(args.input_file, surface=args.surface, facet=args.facet)
    else:
        return args.input_file


def search_command(args):
    """Search catalysis-hub and calculate free energies"""
    try:
        calc = FreeEnergyCalculator(
            temperature=args.temperature,
            voltage=args.voltage,
            use_solvation=args.solvation
        )
        
        results = calc.search_and_calculate(
            reactants=args.reactants.split(',') if args.reactants else None,
            products=args.products.split(',') if args.products else None,
            surface=args.surface,
            facet=args.facet,
            dft_code=args.dft_code,
            dft_functional=args.dft_functional,
            limit=args.limit
        )
        
        if results.empty:
            print("No reactions found!")
            return
        
        print(f"Found {len(results)} reactions:")
        # Select and rename columns for display
        # ∆F_thermo = ∆H - T∆S = ∆ZPE + ∆Cp - T∆S
        # ∆G = ∆E + ∆F_thermo + ∆E_solv
        if args.solvation:
            display_cols = ['equation', 'composition', 'facet', 'delta_E_DFT', 'delta_ZPE', 'delta_Cp', 'temperature', 'delta_S', 'minus_TdS', 'delta_F_thermo', 'delta_E_solv', 'delta_G']
            col_names = ['equation', 'composition', 'facet', '∆E', '∆ZPE', '∆Cp', 'T (K)', '∆S', '-T∆S', '∆F_thermo', '∆E_solv', '∆G']
            numeric_cols = ['∆E', '∆ZPE', '∆Cp', 'T (K)', '∆S', '-T∆S', '∆F_thermo', '∆E_solv', '∆G']
        else:
            display_cols = ['equation', 'composition', 'facet', 'delta_E_DFT', 'delta_ZPE', 'delta_Cp', 'temperature', 'delta_S', 'minus_TdS', 'delta_F_thermo', 'delta_G']
            col_names = ['equation', 'composition', 'facet', '∆E', '∆ZPE', '∆Cp', 'T (K)', '∆S', '-T∆S', '∆F_thermo', '∆G']
            numeric_cols = ['∆E', '∆ZPE', '∆Cp', 'T (K)', '∆S', '-T∆S', '∆F_thermo', '∆G']
        display_df = results[display_cols].copy()
        display_df.columns = col_names
        # Format numeric columns to 2 decimal places
        display_df[numeric_cols] = display_df[numeric_cols].round(2)
        print(display_df.to_string())
        
        if args.output:
            results.to_csv(args.output, index=False)
            print(f"Results saved to {args.output}")
        
    except Exception as e:
        print(f"Error searching reactions: {e}")
        sys.exit(1)


def main():
    """Main CLI entry point"""
    # Allow negative numbers as arguments (not as options)
    parser = argparse.ArgumentParser(
        description="Cathub-Plotter: Free Energy Diagram Plotting for Catalysis Research",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        allow_abbrev=False,  # Disable abbreviations to avoid conflicts
        epilog="""
Examples:
  # Plot free energy diagram from MKM file
  cathub-plotter plot --mkm-file config.yaml --input-file input.txt --mechanism HER_Heyrovsky
  
  # Plot for specific surface and facet (e.g., Cu(100))
  cathub-plotter plot --mkm-file config.yaml --input-file input.txt --surface Cu --facet 100 --mechanism CO2R_CO
  
  # Plot at different temperature and voltage
  cathub-plotter plot --mkm-file config.yaml --input-file input.txt -T 400 -U -0.5 --mechanism CO2R_CO
  
  # Compare mechanisms (1D)
  cathub-plotter compare --mkm-file config.yaml --input-file input.txt --vary mechanism \\
    --mechanisms HER_Heyrovsky,HER_Tafel --surface Cu --facet 100 -T 298.15 -U -0.5
  
  # Compare voltages (1D) - Note: use = sign and quotes for negative values
  cathub-plotter compare --mkm-file config.yaml --input-file input.txt --vary voltage \\
    --voltages="-1.0,-0.8,-0.6,-0.4,-0.2,0.0" --mechanism HER_Heyrovsky --surface Cu --facet 100
  
  # Compare mechanism vs voltage (2D subplots)
  cathub-plotter compare --mkm-file config.yaml --input-file input.txt --vary mechanism,voltage \\
    --mechanisms HER_Heyrovsky,HER_Tafel --voltages="-1.0,-0.5,0.0" --surface Cu --facet 100 \\
    --layout subplots -o comparison.png
  
  # Compare surfaces (1D)
  cathub-plotter compare --mkm-file config.yaml --input-file input.txt --vary surface \\
    --surfaces Cu,Ag,Au --mechanism CO2R_CO --facet 100 -T 298.15 -U -0.5
  
  # Compare with custom figure size
  cathub-plotter compare --mkm-file config.yaml --input-file input.txt --vary surface,temperature \\
    --surfaces Cu,Ni --temperatures 298.15,400,500 --mechanism CO_via_COOH --facet 100 -U -0.5 \\
    --layout overlay --figsize 12 8 -o surface_vs_temp.png --show-labels
  
  # Search catalysis-hub for CO adsorption reactions
  cathub-plotter search --reactants "CO_g,*" --products "CO*" --surface "Cu" --facet "111"
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Plot command
    plot_parser = subparsers.add_parser('plot', help='Plot free energy diagrams')
    plot_parser.add_argument('--mkm-file', required=True, help='Path to MKM/YAML file')
    plot_parser.add_argument('--input-file', required=True, help='Path to input file with species data')
    plot_parser.add_argument('--mechanism', help='Specific mechanism to plot (default: all)')
    plot_parser.add_argument('--surface', help='Surface name to filter (e.g., Cu)')
    plot_parser.add_argument('--facet', help='Facet name to filter (e.g., 100)')
    plot_parser.add_argument('--temperature', '-T', type=float, default=298.15, help='Temperature in K (default: 298.15)')
    plot_parser.add_argument('--voltage', '-U', type=float, default=0.0, help='Voltage in V vs RHE (default: 0.0)')
    plot_parser.add_argument('--y-margin', type=float, default=None, help='Y-axis margin fraction (default: 0.15, increase if labels are clipped)')
    plot_parser.add_argument('--ylim', nargs=2, type=float, metavar=('YMIN', 'YMAX'), help='Y-axis limits (e.g., --ylim -0.1 0.3)')
    plot_parser.add_argument('--x-axis', choices=['step', 'electron'], default='step', help='X-axis mode: "step" for reaction coordinate (default) or "electron" for n(H+ + e-)')
    plot_parser.add_argument('--use-gibbs', action='store_true', help='Treat formation_energy as Gibbs free energy at U=0V (no thermo corrections, only voltage correction)')
    plot_parser.add_argument('--figsize', nargs=2, type=float, metavar=('WIDTH', 'HEIGHT'), help='Figure size in inches (e.g., --figsize 8 6)')
    plot_parser.add_argument('--bar-width', type=float, default=0.6, help='Width of bars in the diagram (default: 0.6)')
    plot_parser.add_argument('--output', '-o', help='Output file path')
    plot_parser.add_argument('--save-data', action='store_true', help='Save raw data (reaction steps, energies) to CSV file')
    plot_parser.set_defaults(func=plot_command)
    
    # Compare command
    compare_parser = subparsers.add_parser('compare', help='Compare free energy diagrams across conditions')
    compare_parser.add_argument('--mkm-file', required=True, help='Path to MKM/YAML file')
    compare_parser.add_argument('--input-file', required=True, help='Path to input file with species data')
    
    # What to vary (required)
    compare_parser.add_argument('--vary', required=True, 
                               help='Parameter(s) to vary: mechanism, temperature, voltage, surface, facet. '
                                    'For 2D comparison, use comma-separated (e.g., "mechanism,voltage")')
    
    # Values for each parameter (provide based on --vary)
    compare_parser.add_argument('--mechanisms', help='Comma-separated mechanism names (e.g., "HER_Heyrovsky,HER_Tafel")')
    compare_parser.add_argument('--temperatures', help='Comma-separated temperatures in K (e.g., "298.15,400,500")')
    compare_parser.add_argument('--voltages', help='Comma-separated voltages in V (e.g., "-1.0,-0.5,0.0")')
    compare_parser.add_argument('--surfaces', help='Comma-separated surface names (e.g., "Cu,Ag,Au")')
    compare_parser.add_argument('--facets', help='Comma-separated facet names (e.g., "100,111,211")')
    
    # Fixed conditions (for parameters not being varied)
    compare_parser.add_argument('--mechanism', help='Mechanism name (when not varying mechanism)')
    compare_parser.add_argument('--surface', help='Surface name (when not varying surface)')
    compare_parser.add_argument('--facet', help='Facet name (when not varying facet)')
    compare_parser.add_argument('-T', '--temperature', type=float, default=298.15, 
                               help='Temperature in K (when not varying temperature, default: 298.15)')
    compare_parser.add_argument('-U', '--voltage', type=float, default=0.0, 
                               help='Voltage in V vs RHE (when not varying voltage, default: 0.0)')
    
    # Plot options
    compare_parser.add_argument('--layout', choices=['subplots', 'overlay'], default='subplots',
                               help='Layout for 2D comparison (default: subplots)')
    compare_parser.add_argument('--figsize', nargs=2, type=float, metavar=('WIDTH', 'HEIGHT'),
                               help='Figure size in inches (e.g., --figsize 12 8). If not specified, size is calculated automatically.')
    compare_parser.add_argument('--colors', help='Comma-separated colors (hex or named colors)')
    compare_parser.add_argument('--labels', help='Comma-separated custom labels')
    compare_parser.add_argument('--show-barriers', action='store_true', help='Show activation barriers')
    compare_parser.add_argument('--show-labels', action='store_true', help='Show state labels')
    compare_parser.add_argument('--legend-position', default='best', 
                               help='Legend position (default: best)')
    compare_parser.add_argument('--x-axis', choices=['step', 'electron'], default='step', 
                               help='X-axis mode: "step" for reaction coordinate (default) or "electron" for n(H+ + e-)')
    compare_parser.add_argument('--use-gibbs', action='store_true', 
                               help='Treat formation_energy as Gibbs free energy at U=0V (no thermo corrections, only voltage correction)')
    compare_parser.add_argument('--bar-width', type=float, default=0.6, help='Width of bars in the diagram (default: 0.6)')
    compare_parser.add_argument('--output', '-o', help='Output file path')
    compare_parser.add_argument('--save-data', action='store_true', 
                               help='Save raw comparison data (energies for all conditions) to CSV file')
    compare_parser.set_defaults(func=compare_command)
    
    # Search command
    search_parser = subparsers.add_parser('search', help='Search catalysis-hub database')
    search_parser.add_argument('--reactants', help='Reactants (comma-separated)')
    search_parser.add_argument('--products', help='Products (comma-separated)')
    search_parser.add_argument('--surface', help='Surface composition (e.g., Cu, Pt)')
    search_parser.add_argument('--facet', help='Surface facet (e.g., 111, 100)')
    search_parser.add_argument('--dft-code', help='DFT code (e.g., VASP, GPAW)')
    search_parser.add_argument('--dft-functional', help='DFT functional (e.g., RPBE, PBE)')
    search_parser.add_argument('--temperature', type=float, default=298.15, help='Temperature in K (default: 298.15)')
    search_parser.add_argument('--voltage', type=float, default=0.0, help='Voltage in V (default: 0.0)')
    search_parser.add_argument('--solvation', action='store_true', help='Apply solvation energy corrections')
    search_parser.add_argument('--limit', type=int, default=50, help='Maximum number of results (default: 50)')
    search_parser.add_argument('--output', '-o', help='Output CSV file path')
    search_parser.set_defaults(func=search_command)
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    args.func(args)


if __name__ == '__main__':
    main()
