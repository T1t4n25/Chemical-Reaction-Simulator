import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Union


class Chemical:
    """
    Represents a chemical species with name, formula, and concentration.
    
    This class encapsulates the properties and state of a single chemical species,
    including validation to ensure physical constraints are maintained.
    """
    
    def __init__(self, name: str, formula: str, concentration: float):
        """
        Initialize a Chemical object.
        
        Args:
            name (str): Name of the chemical species (e.g., "Glucose")
            formula (str): Chemical formula (e.g., "C6H12O6")
            concentration (float): Initial concentration in mol/L
            
        Raises:
            ValueError: If concentration is negative
        """
        self._name = name
        self._formula = formula
        if concentration < 0:
            raise ValueError("Concentration cannot be negative")
        self._concentration = concentration
    
    def get_name(self) -> str:
        """Return the name of the chemical species."""
        return self._name
    
    def get_formula(self) -> str:
        """Return the chemical formula."""
        return self._formula
    
    def get_concentration(self) -> float:
        """Return the current concentration in mol/L."""
        return self._concentration
    
    def set_concentration(self, value: float) -> None:
        """
        Set the concentration of the chemical species.
        
        Args:
            value (float): New concentration in mol/L
            
        Raises:
            ValueError: If value is negative
        """
        if value < 0:
            raise ValueError("Concentration cannot be negative")
        self._concentration = value
    
    def __str__(self) -> str:
        """Return a string representation of the chemical species."""
        return f"{self._name} ({self._formula}), {self._concentration:.3f} mol/L"


class MultiReactantReaction:
    """
    Represents a multi-reactant chemical reaction with kinetics simulation capabilities.
    
    This class encapsulates the reaction parameters and provides methods for calculating
    reaction rates and simulating concentration profiles over time using Euler's method.
    """
    
    def __init__(self, k: float, reaction_orders: Dict[str, float], coeffs: Dict[str, float]):
        """
        Initialize a MultiReactantReaction object.
        
        Args:
            k (float): Rate constant (units depend on overall reaction order)
            reaction_orders (dict): Mapping of species names to their reaction orders
            coeffs (dict): Mapping of species names to stoichiometric coefficients
                          (negative for reactants, positive for products)
        """
        self._k = k
        self._reaction_orders = reaction_orders
        self._coeffs = coeffs
    
    def rate(self, concentrations: Dict[str, float]) -> float:
        """
        Calculate the instantaneous reaction rate.
        
        Args:
            concentrations (dict): Current concentrations of all species
            
        Returns:
            float: Instantaneous reaction rate r
        """
        rate_value = self._k
        for species, order in self._reaction_orders.items():
            if species in concentrations:
                rate_value *= concentrations[species] ** order
        return rate_value
    
    def get_concentration_profile(self, t_max: float, concentrations: Dict[str, Chemical]) -> Dict[str, List[float]]:
        """
        Simulate the reaction and return concentration profiles over time.
        
        Uses Euler's method with 1000 time steps to integrate the differential equations
        governing the concentration changes.
        
        Args:
            t_max (float): Maximum simulation time in seconds
            concentrations (dict): Dictionary mapping species names to Chemical objects
            
        Returns:
            dict: Mapping of species names to lists of concentrations over time
        """
        # Initialize simulation parameters
        delta_t = t_max / 1000
        num_steps = 1000
        
        # Initialize concentration profiles
        profiles = {}
        for species_name in concentrations:
            profiles[species_name] = [concentrations[species_name].get_concentration()]
        
        # Simulation loop using Euler's method
        for step in range(num_steps):
            # Get current concentrations
            current_concs = {name: chem.get_concentration() 
                           for name, chem in concentrations.items()}
            
            # Calculate instantaneous reaction rate
            r = self.rate(current_concs)
            
            # Update concentrations for each species
            for species_name, chemical in concentrations.items():
                if species_name in self._coeffs:
                    # Calculate rate of change: d[species]/dt = coeff * r
                    rate_of_change = self._coeffs[species_name] * r
                    
                    # Apply Euler's method with non-negativity constraint
                    current_conc = chemical.get_concentration()
                    new_conc = max(0.0, current_conc + rate_of_change * delta_t)
                    
                    # Update the chemical object
                    chemical.set_concentration(new_conc)
                    
                    # Store in profile
                    profiles[species_name].append(new_conc)
        
        return profiles
    
    def plot(self, t_max: float, concentrations: Dict[str, Chemical]) -> None:
        """
        Simulate the reaction and plot concentration profiles.
        
        Args:
            t_max (float): Maximum simulation time in seconds
            concentrations (dict): Dictionary mapping species names to Chemical objects
        """
        # Get concentration profiles
        profiles = self.get_concentration_profile(t_max, concentrations)
        
        # Create time array
        delta_t = t_max / 1000
        time_array = np.arange(0, t_max + delta_t, delta_t)
        
        # Create the plot
        plt.figure(figsize=(10, 6))
        
        for species_name, conc_profile in profiles.items():
            # Ensure time and concentration arrays have the same length
            min_length = min(len(time_array), len(conc_profile))
            plt.plot(time_array[:min_length], conc_profile[:min_length], 
                    label=species_name, linewidth=2, marker='o', markersize=2)
        
        # Customize the plot
        plt.xlabel('Time (s)', fontsize=12)
        plt.ylabel('Concentration (mol/L)', fontsize=12)
        plt.title('Concentration vs Time for Multi-Reactant Chemical Reaction', fontsize=14)
        plt.legend(fontsize=10)
        plt.grid(True, alpha=0.3)
        plt.xlim(0, t_max)
        plt.ylim(bottom=0)
        
        # Adjust layout and display
        plt.tight_layout()
        plt.show()
    
    def __str__(self) -> str:
        """
        Return a string representation of the reaction rate law and equation.
        
        Returns:
            str: Formatted string showing rate law and reaction equation
        """
        # Build rate law string
        rate_law = f"r = {self._k}"
        for species, order in self._reaction_orders.items():
            rate_law += f" × [{species}]^{order}"
        
        # Build reaction equation
        reactants = []
        products = []
        
        for species, coeff in self._coeffs.items():
            if coeff < 0:
                reactants.append(f"{abs(coeff)} {species}")
            elif coeff > 0:
                products.append(f"{coeff} {species}")
        
        reaction_eq = " + ".join(reactants) + " → " + " + ".join(products)
        
        return f"{rate_law} (Reaction: {reaction_eq})"


# Example usage and demonstration
def demo_simulation():
    """
    Demonstrate the chemical reaction simulation with a sample reaction.
    
    Example: H2O2 + 3I- + 2H+ → I3- + 2H2O
    Rate law: r = k[H2O2]^1[I-]^1[H+]^0
    """
    print("=== Chemical Reaction Simulation Demo ===\n")
    
    # Create chemical species
    chemicals = {
        'H2O2': Chemical('Hydrogen Peroxide', 'H2O2', 0.1),
        'I-': Chemical('Iodide Ion', 'I-', 0.3),
        'H+': Chemical('Hydrogen Ion', 'H+', 0.2),
        'I3-': Chemical('Triiodide Ion', 'I3-', 0.0),
        'H2O': Chemical('Water', 'H2O', 0.0)
    }
    
    # Display initial concentrations
    print("Initial Concentrations:")
    for name, chemical in chemicals.items():
        print(f"  {chemical}")
    print()
    
    # Define reaction parameters
    rate_constant = 10.0  # L²/(mol²·s)
    reaction_orders = {'H2O2': 1, 'I-': 1, 'H+': 0}
    stoichiometric_coeffs = {
        'H2O2': -1,  # reactant
        'I-': -3,    # reactant
        'H+': -2,    # reactant
        'I3-': 1,    # product
        'H2O': 2     # product
    }
    
    # Create reaction object
    reaction = MultiReactantReaction(rate_constant, reaction_orders, stoichiometric_coeffs)
    
    # Display reaction information
    print("Reaction Details:")
    print(f"  {reaction}")
    print()
    
    # Run simulation
    t_max = 5.0  # seconds
    print(f"Running simulation for {t_max} seconds...")
    
    # Make copies for simulation (to preserve original objects)
    sim_chemicals = {}
    for name, chem in chemicals.items():
        sim_chemicals[name] = Chemical(chem.get_name(), chem.get_formula(), chem.get_concentration())
    
    # Get concentration profiles
    profiles = reaction.get_concentration_profile(t_max, sim_chemicals)
    
    # Display final concentrations
    print("\nFinal Concentrations:")
    for name, chemical in sim_chemicals.items():
        print(f"  {chemical}")
    print()
    
    # Calculate and display reaction progress
    initial_h2o2 = chemicals['H2O2'].get_concentration()
    final_h2o2 = sim_chemicals['H2O2'].get_concentration()
    conversion = (initial_h2o2 - final_h2o2) / initial_h2o2 * 100
    print(f"H2O2 Conversion: {conversion:.1f}%")
    
    # Plot results
    print("\nGenerating concentration vs time plot...")
    reaction.plot(t_max, chemicals)  # Use fresh chemicals for plotting
    
    return profiles


if __name__ == "__main__":
    # Run the demonstration
    concentration_profiles = demo_simulation()
    
    # Additional example: Simple A + B → C reaction
    print("\n" + "="*60)
    print("=== Simple A + B → C Reaction Example ===\n")
    
    # Create simple reaction system
    simple_chemicals = {
        'A': Chemical('Reactant A', 'A', 1.0),
        'B': Chemical('Reactant B', 'B', 0.8),
        'C': Chemical('Product C', 'C', 0.0)
    }
    
    simple_reaction = MultiReactantReaction(
        k=0.5,
        reaction_orders={'A': 1, 'B': 1},
        coeffs={'A': -1, 'B': -1, 'C': 1}
    )
    
    print(f"Reaction: {simple_reaction}")
    print("Plotting simple reaction...")
    simple_reaction.plot(10.0, simple_chemicals)