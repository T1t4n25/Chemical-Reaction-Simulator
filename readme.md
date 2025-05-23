# Chemical Reaction Simulation System

A comprehensive Python-based simulation system for modeling multi-reactant chemical reactions using object-oriented programming principles and numerical integration methods.

## 🔬 Overview

This simulation system allows you to model complex chemical reactions involving multiple reactants and products. It uses Euler's method for numerical integration to solve the system of differential equations that govern concentration changes over time, providing accurate visualization of reaction kinetics.

### Key Features

- **Object-Oriented Design**: Clean, modular architecture with separate classes for chemical species and reactions
- **Multi-Reactant Support**: Handle complex reactions with multiple reactants and products
- **Numerical Integration**: Robust Euler's method implementation with physical constraints
- **Real-time Visualization**: Generate publication-quality concentration vs. time plots
- **Physical Validation**: Built-in constraints to ensure concentrations remain non-negative
- **Flexible Rate Laws**: Support for custom reaction orders and stoichiometric coefficients

## 📋 Requirements

```bash
numpy>=1.20.0
matplotlib>=3.3.0
```

## 🚀 Quick Start

### Installation

1. Ensure you have Python 3.7+ installed
2. Install required dependencies:
```bash
pip install numpy matplotlib
```
3. Download the `chemical_simulation.py` file

### Basic Usage

```python
from chemical_simulation import Chemical, MultiReactantReaction

# Create chemical species
chemicals = {
    'A': Chemical('Reactant A', 'A', 1.0),  # name, formula, initial_concentration
    'B': Chemical('Reactant B', 'B', 0.8),
    'C': Chemical('Product C', 'C', 0.0)
}

# Define reaction: A + B → C with rate law r = k[A]¹[B]¹
reaction = MultiReactantReaction(
    k=0.5,                                    # rate constant
    reaction_orders={'A': 1, 'B': 1},        # reaction orders
    coeffs={'A': -1, 'B': -1, 'C': 1}       # stoichiometric coefficients
)

# Run simulation and plot results
reaction.plot(t_max=10.0, concentrations=chemicals)
```

## 🏗️ Architecture

### Chemical Class

Represents individual chemical species with the following features:

**Attributes:**
- `_name`: Species name (e.g., "Glucose")
- `_formula`: Chemical formula (e.g., "C6H12O6")  
- `_concentration`: Current concentration in mol/L

**Methods:**
- `get_name()`, `get_formula()`, `get_concentration()`: Accessor methods
- `set_concentration(value)`: Updates concentration with validation
- `__str__()`: Returns formatted string representation

### MultiReactantReaction Class

Manages reaction parameters and simulation logic:

**Attributes:**
- `_k`: Rate constant (units depend on overall reaction order)
- `_reaction_orders`: Dictionary mapping species to their reaction orders
- `_coeffs`: Dictionary mapping species to stoichiometric coefficients

**Methods:**
- `rate(concentrations)`: Calculates instantaneous reaction rate
- `get_concentration_profile(t_max, concentrations)`: Returns concentration histories
- `plot(t_max, concentrations)`: Simulates and visualizes the reaction
- `__str__()`: Returns formatted rate law and reaction equation

## 📊 Examples

### Example 1: Hydrogen Peroxide Decomposition
```python
# H2O2 + 3I- + 2H+ → I3- + 2H2O
chemicals = {
    'H2O2': Chemical('Hydrogen Peroxide', 'H2O2', 0.1),
    'I-': Chemical('Iodide Ion', 'I-', 0.3),
    'H+': Chemical('Hydrogen Ion', 'H+', 0.2),
    'I3-': Chemical('Triiodide Ion', 'I3-', 0.0),
    'H2O': Chemical('Water', 'H2O', 0.0)
}

reaction = MultiReactantReaction(
    k=10.0,
    reaction_orders={'H2O2': 1, 'I-': 1, 'H+': 0},
    coeffs={'H2O2': -1, 'I-': -3, 'H+': -2, 'I3-': 1, 'H2O': 2}
)

reaction.plot(t_max=5.0, concentrations=chemicals)
```

### Example 2: Second-Order Reaction
```python
# 2A → B with rate law r = k[A]²
chemicals = {
    'A': Chemical('Reactant A', 'A', 2.0),
    'B': Chemical('Product B', 'B', 0.0)
}

reaction = MultiReactantReaction(
    k=0.1,
    reaction_orders={'A': 2},
    coeffs={'A': -2, 'B': 1}
)

# Get numerical data
profiles = reaction.get_concentration_profile(t_max=20.0, concentrations=chemicals)
print("Final [A]:", profiles['A'][-1])
print("Final [B]:", profiles['B'][-1])
```

## ⚙️ Technical Details

### Numerical Integration

The system uses **Euler's method** for solving the system of ordinary differential equations:

```
[Species]ₜ₊Δₜ = max(0, [Species]ₑ + (d[Species]/dt) × Δt)
```

- **Time step**: Δt = t_max / 1000 (1000 integration steps)
- **Non-negativity constraint**: `max(0, ...)` ensures physical realism
- **Accuracy**: First-order method suitable for most chemical kinetics applications

### Rate Law Calculation

For a general reaction with rate law `r = k[A]ᵃ[B]ᵇ[C]ᶜ...`:

```python
rate = k × ∏(concentration[species]^order[species])
```

### Concentration Updates

Each species concentration changes according to:

```
d[Species]/dt = stoichiometric_coefficient × reaction_rate
```

## 🎨 Visualization Features

The plotting system generates professional-quality graphs with:

- **Time series plots** for all species concentrations
- **Automatic legend** with species identification  
- **Grid lines** for easy reading
- **Axis labels** with proper units
- **Non-negative y-axis** constraint
- **Customizable time range**

## ⚠️ Important Notes

### Physical Constraints
- Concentrations cannot be negative (enforced by both Euler's method and Chemical class)
- Rate constants must be positive
- Stoichiometric coefficients: negative for reactants, positive for products

### Numerical Considerations
- **Step size**: Fixed at t_max/1000 for smooth curves
- **Stability**: Generally stable for most chemical systems
- **Accuracy**: Suitable for visualization and educational purposes
- **Limitations**: First-order method may accumulate errors for very long simulations

### Performance
- **Simulation time**: ~1000 time steps per run
- **Memory usage**: Stores complete concentration history
- **Scalability**: Handles dozens of species efficiently

## 🔧 Customization

### Custom Rate Laws
```python
# Zero-order kinetics: r = k
reaction_orders = {}  # No concentration dependence

# Mixed orders: r = k[A]^0.5[B]^1.5  
reaction_orders = {'A': 0.5, 'B': 1.5}
```

### Complex Stoichiometry
```python
# 2A + 3B → 4C + D
coeffs = {'A': -2, 'B': -3, 'C': 4, 'D': 1}
```

## 📈 Output Data Structure

The `get_concentration_profile()` method returns a dictionary:

```python
{
    'species_name': [conc_at_t0, conc_at_t1, ..., conc_at_tfinal],
    # ... for all species
}
```

Each list contains 1001 concentration values (including t=0).

## 🧪 Validation and Testing

### Built-in Validation
- Negative concentration detection and prevention
- Input parameter type checking
- Physical constraint enforcement

### Testing Your Simulations
```python
# Mass balance check for A → B reaction
initial_total = chemicals['A'].get_concentration()
profiles = reaction.get_concentration_profile(10.0, chemicals)
final_total = profiles['A'][-1] + profiles['B'][-1]
print(f"Mass balance error: {abs(initial_total - final_total):.6f}")
```

### Potential Extensions
- Higher-order numerical methods (Runge-Kutta)
- Reversible reactions
- Temperature-dependent rate constants
- Reaction networks and parallel pathways
- Sensitivity analysis tools

## 📚 Scientific Background

The simulation is based on fundamental principles of chemical kinetics:

- **Rate laws** describe how reaction rates depend on concentrations
- **Stoichiometry** determines the proportional consumption/formation of species
- **Differential equations** govern concentration changes over time
- **Numerical integration** approximates solutions to these equations

This approach is widely used in chemical engineering, pharmaceutical research, and environmental modeling.

