import sys
import importlib
from Environment import Environment

# Store the current environment globally
current_env = None

def get_environment_from_user():
    print("\nSelect the external factors you want to modify (press Enter to use default):")
    try:
        temperature = float(input("Temperature [default: 39]: ") or 39)
        kinase_level = float(input("Kinase level [default: 1.5]: ") or 1.5)
        phosphatase_level = float(input("Phosphatase level [default: 1.0]: ") or 1.0)
        protease_level = float(input("Protease level [default: 1.0]: ") or 1.0)
        oxidative_stress = float(input("Oxidative stress [default: 0.2]: ") or 0.2)
    except ValueError:
        print("Invalid input. Using default values.")
        temperature, kinase_level, phosphatase_level, protease_level, oxidative_stress = 39, 1.5, 1.0, 1.0, 0.2

    return Environment(
        temperature=temperature,
        kinase_level=kinase_level,
        phosphatase_level=phosphatase_level,
        protease_level=protease_level,
        oxidative_stress=oxidative_stress
    )

def run_tau_simulation():
    import tau_simulation
    global current_env
    if current_env is None:
        current_env = Environment(temperature=39, kinase_level=1.5, oxidative_stress=0.2)
    tau_simulation.run_and_plot_simulation(current_env)

def show_tau_visualization():
    print("\n[Showing tau protein phosphorylation & aggregation visualization...]")
    import tau_simulation
    tau_simulation.run_and_plot_simulation(plot_sites=True, plot_heatmap=True)
    print("\n[Visualization complete.]")

def show_disease_visualization():
    import Disease_sim
    Disease_sim.run_and_plot_disease_simulation()

def show_help():
    print("""
Tau Protein Simulator Chatbot Menu:
1. Run tau protein simulation: Simulates tau phosphorylation and aggregation, and shows the main plot.
2. Show tau protein visualization: Displays the main tau simulation plot (phosphorylation & aggregation).
3. Show disease simulation visualization: Displays a plot of phosphorylation percentage over time under disease-like conditions.
4. Exit: Quit the chatbot.
""")

def set_environment():
    global current_env
    current_env = get_environment_from_user()
    print("\n[Environment updated!]\n")

def main_menu():
    while True:
        try:
            print("\n=== Tau Protein Simulator Menu ===")
            print("1. Set external environment factors")
            print("2. Run tau protein simulation (with current environment)")
            print("3. Show tau protein phosphorylation & aggregation visualization")
            print("4. Show disease simulation visualization")
            print("5. Help/About")
            print("6. Exit")
            choice = input("Select an option (1-6): ").strip()
            if choice == "1":
                set_environment()
            elif choice == "2":
                run_tau_simulation()
            elif choice == "3":
                show_tau_visualization()
            elif choice == "4":
                show_disease_visualization()
            elif choice == "5":
                show_help()
            elif choice == "6":
                print("Goodbye!")
                break
            else:
                print("Invalid choice. Please select 1-6.")
        except Exception as e:
            print(f"\n[Error]: {e}\nReturning to menu...")

if __name__ == "__main__":
    main_menu()
