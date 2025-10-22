"""
Command-line interface for cathub-plotter
"""
import argparse
import sys
from pathlib import Path
from typing import List, Optional

from ..core.calculator import FreeEnergyCalculator
from ..plotters.diagram import FreeEnergyDiagramPlotter


def plot_command(args):
    """Plot free energy diagram from MKM file"""
    try:
        plotter = FreeEnergyDiagramPlotter(
            mkm_file=args.mkm_file,
            input_file=args.input_file,
            temperature=args.temperature,
            voltage=args.voltage
        )
        
        if args.mechanism:
            plotter.plot_mechanism(args.mechanism, save_path=args.output)
        else:
            # Plot all mechanisms
            mechanisms = list(plotter.rxn_mechanisms.keys())
            for mech in mechanisms:
                output_path = f"{args.output}_{mech}.png" if args.output else None
                plotter.plot_mechanism(mech, save_path=output_path)
        
        print("Plotting completed successfully!")
        
    except Exception as e:
        print(f"Error plotting diagram: {e}")
        sys.exit(1)


def search_command(args):
    """Search catalysis-hub and calculate free energies"""
    try:
        calc = FreeEnergyCalculator(
            temperature=args.temperature,
            voltage=args.voltage
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
        print(results[['equation', 'delta_G', 'composition', 'facet']].to_string())
        
        if args.output:
            results.to_csv(args.output, index=False)
            print(f"Results saved to {args.output}")
        
    except Exception as e:
        print(f"Error searching reactions: {e}")
        sys.exit(1)


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        description="Cathub-Plotter: Free Energy Diagram Plotting for Catalysis Research",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Plot free energy diagram from MKM file
  cathub-plotter plot --mkm-file config.yaml --input-file input.txt --mechanism HER_Heyrovsky
  
  # Search catalysis-hub for CO adsorption reactions
  cathub-plotter search --reactants "CO_g,*" --products "CO*" --surface "Cu" --facet "111"
  
  # Calculate at different temperature
  cathub-plotter plot --mkm-file config.yaml --input-file input.txt --temperature 400
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Plot command
    plot_parser = subparsers.add_parser('plot', help='Plot free energy diagrams')
    plot_parser.add_argument('--mkm-file', required=True, help='Path to MKM/YAML file')
    plot_parser.add_argument('--input-file', required=True, help='Path to input file with species data')
    plot_parser.add_argument('--mechanism', help='Specific mechanism to plot (default: all)')
    plot_parser.add_argument('--temperature', type=float, default=298.15, help='Temperature in K (default: 298.15)')
    plot_parser.add_argument('--voltage', type=float, default=0.0, help='Voltage in V (default: 0.0)')
    plot_parser.add_argument('--output', '-o', help='Output file path')
    plot_parser.set_defaults(func=plot_command)
    
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
