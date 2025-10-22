# app.py - Streamlit Application (with Free Energy Diagram)
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from cathub_plotter import FreeEnergyCalculator

st.set_page_config(page_title="Catalysis-Hub Free Energy Calculator", layout="wide")

def create_free_energy_diagram(equation, delta_G, temperature, composition, facet):
    """
    Create a smooth free energy diagram for a single reaction step
    
    Args:
        equation: Reaction equation string (e.g., "CO_g + * -> CO*")
        delta_G: Gibbs free energy change in eV
        temperature: Temperature in K
        composition: Surface composition (e.g., "Cu36")
        facet: Surface facet (e.g., "100")
    
    Returns:
        Plotly figure object
    """
    # Parse equation to get reactants and products
    parts = equation.split(' -> ')
    if len(parts) != 2:
        return None
    
    reactants = parts[0].strip()
    products = parts[1].strip()
    
    # Create smooth curve for the diagram
    import numpy as np
    
    # Define the x positions for flat regions and transition
    x_initial_start = 0
    x_initial_end = 0.35
    x_transition_start = 0.35
    x_transition_end = 0.65
    x_final_start = 0.65
    x_final_end = 1.0
    
    # Create smooth transition using a sigmoid-like curve
    n_transition_points = 50
    x_transition = np.linspace(x_transition_start, x_transition_end, n_transition_points)
    
    # Use a smooth S-curve (cubic interpolation)
    t = np.linspace(0, 1, n_transition_points)
    # Smooth cubic curve: 3t^2 - 2t^3
    smooth_t = 3 * t**2 - 2 * t**3
    y_transition = smooth_t * delta_G
    
    # Combine all segments
    x_positions = np.concatenate([
        [x_initial_start, x_initial_end],  # Initial flat
        x_transition,                       # Smooth transition
        [x_final_start, x_final_end]       # Final flat
    ])
    
    y_positions = np.concatenate([
        [0, 0],                            # Initial flat
        y_transition,                       # Smooth transition
        [delta_G, delta_G]                 # Final flat
    ])
    
    fig = go.Figure()
    
    # Add the smooth energy profile line
    fig.add_trace(go.Scatter(
        x=x_positions,
        y=y_positions,
        mode='lines',
        line=dict(color='black', width=4),
        showlegend=False,
        hoverinfo='skip'
    ))
    
    # Add labels for initial state
    fig.add_annotation(
        x=(x_initial_start + x_initial_end) / 2,
        y=0,
        text=f"<b>{reactants}</b>",
        showarrow=False,
        yshift=-25,
        font=dict(size=14, color='black')
    )
    
    # Add labels for final state
    fig.add_annotation(
        x=(x_final_start + x_final_end) / 2,
        y=delta_G,
        text=f"<b>{products}</b>",
        showarrow=False,
        yshift=-25,
        font=dict(size=14, color='black')
    )
    
    # Determine color based on sign of delta_G
    arrow_color = 'red' if delta_G > 0 else 'green'
    reaction_type = "Endergonic" if delta_G > 0 else "Exergonic"
    
    # Update layout
    y_range_buffer = max(abs(delta_G) * 0.4, 0.3)
    y_min = min(0, delta_G) - y_range_buffer
    y_max = max(0, delta_G) + y_range_buffer
    
    # Create title with reaction details
    title_text = (
        f"<b style='font-size:22px;'>{equation} on {composition}({facet}) at {temperature} K</b><br>"
        f"<span style='font-size:18px; color:{arrow_color};'>ΔG = {delta_G:.3f} eV ({reaction_type})</span>"
    )
    
    fig.update_layout(
        title=dict(
            text=title_text,
            x=0.5,
            xanchor='center',
            y=0.95,
            yanchor='top'
        ),
        xaxis=dict(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            range=[-0.05, 1.05]
        ),
        yaxis=dict(
            title="<b>Relative Free Energy (eV)</b>",
            title_font=dict(size=14),
            showgrid=True,
            gridcolor='lightgray',
            gridwidth=1,
            zeroline=True,
            zerolinewidth=2,
            zerolinecolor='black',
            range=[y_min, y_max],
            tickfont=dict(size=12)
        ),
        plot_bgcolor='white',
        width=800,
        height=500,
        margin=dict(l=80, r=80, t=120, b=60)
    )
    
    return fig

st.title("🧪 Catalysis-Hub Free Energy Calculator")
st.markdown("Search reactions from catalysis-hub.org and calculate Gibbs free energies")

# Sidebar - Search Parameters
st.sidebar.header("🔍 Search Parameters")

with st.sidebar.expander("Reaction Species", expanded=True):
    reactants_input = st.text_input("Reactants (comma-separated)", 
                                     placeholder="e.g., CO_g, *",
                                     help="Use * for empty site, _g for gas")
    products_input = st.text_input("Products (comma-separated)", 
                                    placeholder="e.g., CO*",
                                    help="Use * for empty site, _g for gas")

with st.sidebar.expander("Surface Parameters", expanded=True):
    surface = st.text_input("Surface composition", 
                           placeholder="e.g., Cu, Pt, Pd",
                           help="Leave empty for all surfaces")
    facet = st.selectbox("Facet", 
                        ["All", "111", "100", "110", "211"],
                        help="Crystal facet")

with st.sidebar.expander("DFT Parameters"):
    dft_code = st.text_input("DFT Code", 
                             placeholder="e.g., VASP, GPAW",
                             help="Partial match supported")
    dft_functional = st.text_input("DFT Functional", 
                                   placeholder="e.g., RPBE, PBE",
                                   help="Partial match (e.g., RPBE matches BEEF-vdW-RPBE)")

limit = st.sidebar.slider("Max results", 10, 200, 50)

# Sorting options
st.sidebar.header("📊 Display Options")
sort_by = st.sidebar.selectbox("Sort by", 
                               ["composition", "delta_G", "delta_E_DFT", "delta_F_thermo"])
ascending = st.sidebar.radio("Order", ["Ascending", "Descending"]) == "Ascending"

# Search button
search_button = st.sidebar.button("🔍 Search & Calculate", type="primary")

# Initialize session state for storing results
if 'df' not in st.session_state:
    st.session_state.df = None
if 'selected_reaction' not in st.session_state:
    st.session_state.selected_reaction = None

# Main content
if search_button:
    # Parse inputs
    reactants = [r.strip() for r in reactants_input.split(",")] if reactants_input else None
    products = [p.strip() for p in products_input.split(",")] if products_input else None
    facet_val = None if facet == "All" else facet
    surface_val = surface if surface else None
    dft_code_val = dft_code if dft_code else None
    dft_functional_val = dft_functional if dft_functional else None
    
    # Calculate with default conditions first
    with st.spinner("Searching catalysis-hub.org and calculating free energies..."):
        calc = FreeEnergyCalculator(
            temperature=298.15,  # Default temperature
            voltage=0.0,
            fugacity_dict=None
        )
        
        df = calc.search_and_calculate(
            reactants=reactants,
            products=products,
            surface=surface_val,
            facet=facet_val,
            dft_code=dft_code_val,
            dft_functional=dft_functional_val,
            limit=limit,
            sort_by=sort_by,
            ascending=ascending
        )
        
        st.session_state.df = df
        st.session_state.selected_reaction = None

# Display results if available
if st.session_state.df is not None and not st.session_state.df.empty:
    df = st.session_state.df
    
    st.success(f"Found {len(df)} reactions!")
    
    # Thermodynamic parameter controls (above the table)
    st.markdown("### 🌡️ Thermodynamic Conditions")
    
    col1, col2, col3 = st.columns([2, 2, 1])
    
    with col1:
        temperature = st.number_input("Temperature (K)", 
                                     min_value=200.0, 
                                     max_value=1000.0, 
                                     value=298.15, 
                                     step=10.0,
                                     key='temp_input')
    
    with col2:
        use_custom_fugacity = st.checkbox("Use custom fugacity (Pa)", key='fugacity_checkbox')
    
    with col3:
        recalculate_button = st.button("🔄 Recalculate", type="secondary")
    
    # Fugacity inputs
    fugacity_dict = {}
    if use_custom_fugacity:
        st.markdown("##### Gas Fugacity (Pa)")
        col1, col2, col3 = st.columns(3)
        with col1:
            fugacity_dict['H2_g'] = st.number_input("H₂", value=101325, format="%d", key='h2_fug')
            fugacity_dict['CO_g'] = st.number_input("CO", value=5562, format="%d", key='co_fug')
        with col2:
            fugacity_dict['CO2_g'] = st.number_input("CO₂", value=101325, format="%d", key='co2_fug')
            fugacity_dict['H2O_g'] = st.number_input("H₂O", value=3534, format="%d", key='h2o_fug')
        with col3:
            fugacity_dict['CH4_g'] = st.number_input("CH₄", value=20467, format="%d", key='ch4_fug')
            fugacity_dict['CH3OH_g'] = st.number_input("CH₃OH", value=6079, format="%d", key='ch3oh_fug')
    
    # Recalculate with new conditions
    if recalculate_button:
        with st.spinner("Recalculating free energies..."):
            # Get original search parameters
            reactants = [r.strip() for r in reactants_input.split(",")] if reactants_input else None
            products = [p.strip() for p in products_input.split(",")] if products_input else None
            facet_val = None if facet == "All" else facet
            surface_val = surface if surface else None
            dft_code_val = dft_code if dft_code else None
            dft_functional_val = dft_functional if dft_functional else None
            
            calc = FreeEnergyCalculator(
                temperature=temperature,
                voltage=0.0,
                fugacity_dict=fugacity_dict if use_custom_fugacity else None
            )
            
            df = calc.search_and_calculate(
                reactants=reactants,
                products=products,
                surface=surface_val,
                facet=facet_val,
                dft_code=dft_code_val,
                dft_functional=dft_functional_val,
                limit=limit,
                sort_by=sort_by,
                ascending=ascending
            )
            
            st.session_state.df = df
            st.rerun()
    
    st.markdown("---")
    
    # Tabs
    tab1, tab2, tab3 = st.tabs(["📋 Data Table", "📊 Visualizations", "📥 Export"])
    
    with tab1:
        st.subheader("Reaction Data")
        
        # Display columns selector
        all_columns = ['equation', 'composition', 'facet', 'sites_str', 
                      'dft_code', 'dft_functional', 'delta_E_DFT', 
                      'delta_F_thermo', 'delta_G', 'pub_authors', 'pub_year']
        
        display_cols = st.multiselect(
            "Select columns to display",
            [col for col in all_columns if col in df.columns],
            default=['equation', 'composition', 'facet', 'sites_str', 
                    'dft_functional', 'delta_E_DFT', 'delta_G']
        )
        
        # Create clickable dataframe
        st.markdown("*Click on a row to view detailed information*")
        
        # Display dataframe with selection
        event = st.dataframe(
            df[display_cols], 
            use_container_width=True, 
            height=400,
            on_select="rerun",
            selection_mode="single-row"
        )
        
        # Handle row selection
        if event.selection and event.selection.rows:
            selected_idx = event.selection.rows[0]
            st.session_state.selected_reaction = selected_idx
        
        # Show detailed view if a reaction is selected
        if st.session_state.selected_reaction is not None:
            st.markdown("---")
            st.markdown("### 🔬 Detailed Reaction Information")
            
            row = df.iloc[st.session_state.selected_reaction]
            
            # Free Energy Diagram (full width at top)
            st.markdown("#### Free Energy Diagram")
            fig = create_free_energy_diagram(
                equation=row['equation'],
                delta_G=row['delta_G'],
                temperature=temperature,
                composition=row['composition'],  # 추가
                facet=row['facet']              # 추가
            )
            if fig:
                st.plotly_chart(fig, use_container_width=True)
            
            st.markdown("---")
            
            # Two column layout for details
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("#### Reaction Details")
                st.markdown(f"**Equation:** {row['equation']}")
                st.markdown(f"**Surface:** {row['composition']} ({row['facet']})")
                st.markdown(f"**Active Sites:** {row['sites_str']}")
                st.markdown(f"**DFT Code:** {row['dft_code']}")
                st.markdown(f"**DFT Functional:** {row['dft_functional']}")
                
                st.markdown("#### Publication")
                st.markdown(f"**Authors:** {row['pub_authors']}")
                st.markdown(f"**Year:** {row['pub_year']}")
                st.markdown(f"**Pub ID:** {row['pub_id']}")
            
            with col2:
                st.markdown("#### Energetics")
                st.markdown(f"**ΔE_DFT:** {row['delta_E_DFT']:.4f} eV")
                st.markdown(f"**ΔF_thermo:** {row['delta_F_thermo']:.4f} eV")
                st.markdown(f"**ΔG:** {row['delta_G']:.4f} eV")
                
                if pd.notna(row['activation_energy']):
                    st.markdown(f"**E_activation:** {row['activation_energy']:.4f} eV")
                
                st.markdown("#### Free Energy Correction Breakdown")
                st.markdown(f"**F(reactants):** {row['F_reactants']:.4f} eV")
                st.markdown(f"**F(products):** {row['F_products']:.4f} eV")
                st.markdown(f"**ΔF_thermo = F(products) - F(reactants)**")
                st.markdown(f"**ΔF_thermo = {row['F_products']:.4f} - {row['F_reactants']:.4f} = {row['delta_F_thermo']:.4f} eV**")
                
                st.markdown("#### Final Calculation")
                st.code(f"""
ΔG_reaction = ΔE_DFT + ΔF_thermo
ΔG_reaction = {row['delta_E_DFT']:.4f} + {row['delta_F_thermo']:.4f}
ΔG_reaction = {row['delta_G']:.4f} eV

where:
  ΔE_DFT: DFT energy from catalysis-hub
  ΔF_thermo: Thermodynamic correction (ZPE - TS)
  Temperature: {temperature} K
                """, language="text")
            
            # Add a close button
            if st.button("✖ Close Details"):
                st.session_state.selected_reaction = None
                st.rerun()
    
    with tab2:
        st.subheader("Energy Comparison")
        
        # Bar chart: ΔE_DFT vs ΔG
        fig1 = go.Figure()
        fig1.add_trace(go.Bar(
            name='ΔE_DFT',
            x=df.index,
            y=df['delta_E_DFT'],
            marker_color='lightblue'
        ))
        fig1.add_trace(go.Bar(
            name='ΔG',
            x=df.index,
            y=df['delta_G'],
            marker_color='darkblue'
        ))
        fig1.update_layout(
            title=f"Energy Comparison at {temperature}K",
            xaxis_title="Reaction Index",
            yaxis_title="Energy (eV)",
            barmode='group',
            height=400
        )
        st.plotly_chart(fig1, use_container_width=True)
        
        # Scatter plot: ΔE_DFT vs ΔG
        fig2 = px.scatter(
            df, 
            x='delta_E_DFT', 
            y='delta_G',
            color='dft_functional',
            hover_data=['equation', 'composition', 'facet'],
            title="ΔE_DFT vs ΔG",
            labels={'delta_E_DFT': 'ΔE_DFT (eV)', 'delta_G': 'ΔG (eV)'}
        )
        fig2.add_shape(
            type='line',
            x0=df['delta_E_DFT'].min(),
            y0=df['delta_E_DFT'].min(),
            x1=df['delta_E_DFT'].max(),
            y1=df['delta_E_DFT'].max(),
            line=dict(color='red', dash='dash', width=1),
            name='y=x'
        )
        st.plotly_chart(fig2, use_container_width=True)
        
        # Thermodynamic correction distribution
        fig3 = px.histogram(
            df,
            x='delta_F_thermo',
            nbins=30,
            title="Distribution of Thermodynamic Corrections",
            labels={'delta_F_thermo': 'ΔF_thermo (eV)', 'count': 'Frequency'}
        )
        st.plotly_chart(fig3, use_container_width=True)
    
    with tab3:
        st.subheader("Export Data")
        
        # CSV export
        csv = df.to_csv(index=False)
        st.download_button(
            label="📥 Download as CSV",
            data=csv,
            file_name=f"cathub_reactions_{temperature}K.csv",
            mime="text/csv"
        )
        
        # Show summary statistics
        st.markdown("### Data Summary")
        st.dataframe(df[['delta_E_DFT', 'delta_F_thermo', 'delta_G']].describe())

elif search_button and (st.session_state.df is None or st.session_state.df.empty):
    st.warning("No reactions found! Try different search parameters.")

else:
    st.info("👈 Configure search parameters in the sidebar and click 'Search & Calculate'")
    
    # Show example
    with st.expander("💡 Example Searches"):
        st.markdown("""
        **Example 1: CO adsorption on Cu(111)**
        - Reactants: `CO_g, *`
        - Products: `CO*`
        - Surface: `Cu`
        - Facet: `111`
        
        **Example 2: CO oxidation with RPBE functional**
        - Reactants: `CO*, O*`
        - Products: `CO2_g, *, *`
        - DFT Functional: `RPBE`
        
        **Example 3: H2 dissociation on Pt**
        - Reactants: `H2_g, *, *`
        - Products: `H*, H*`
        - Surface: `Pt`
        """)
    
    st.markdown("---")
    st.markdown("""
    ### Features
    - 🔍 Search catalysis-hub.org database
    - 🌡️ Calculate free energies at different temperatures
    - 💨 Adjust gas fugacities/concentrations
    - 📈 Free energy diagrams for each reaction
    - 📊 Interactive visualizations
    - 🔬 Detailed thermodynamic breakdown
    - 📥 Export results to CSV
    """)