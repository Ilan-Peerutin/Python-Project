"""
tau_simulation.py
Runs a simulation of tau protein phosphorylation and aggregation over time under specified environmental conditions, and visualizes the results.
"""

import numpy as np
from TauProtein import TauProtein
from Environment import Environment
import matplotlib.pyplot as plt

# --- Main Simulation and Visualization Functions ---

def plot_tau_summary(history):
    """
    Plots phosphorylation count (P > 0.5), average phosphorylation probability, and aggregation state over time.
    """
    # Extract time series and relevant metrics from simulation history
    times = [entry['minute'] for entry in history]
    phospho_counts = [entry['phospho_count'] for entry in history]
    avg_probs = [entry.get('avg_prob', 0) for entry in history]
    aggregation_states = [entry['aggregation_state'] for entry in history]
    agg_state_numeric = {'monomer': 0, 'oligomer': 1, 'fibril': 2}
    aggregation_numeric = [agg_state_numeric[state] for state in aggregation_states]

    # Create the main plot with dual y-axes
    fig, ax1 = plt.subplots(figsize=(10,6))
    ax1.set_xlabel('Time Step')
    ax1.set_ylabel('Phosphorylation Count (P > 0.5)', color='tab:blue')
    ax1.plot(times, phospho_counts, color='tab:blue', label='Phosphorylation Count')
    ax1.plot(times, avg_probs, color='tab:orange', label='Average Phosphorylation Probability')
    ax1.tick_params(axis='y', labelcolor='tab:blue')

    # Add aggregation state as a step plot on a secondary y-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel('Aggregation State', color='tab:green')
    ax2.step(times, aggregation_numeric, color='tab:green', where='post', label='Aggregation State')
    ax2.set_yticks([0, 1, 2])
    ax2.set_yticklabels(['Monomer', 'Oligomer', 'Fibril'])
    ax2.tick_params(axis='y', labelcolor='tab:green')

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left')

    plt.title('Tau Phosphorylation and Aggregation Over Time')
    fig.tight_layout()
    plt.show()

    # Print summary statistics for insight
    print(f"Final average phosphorylation probability: {avg_probs[-1]:.3f}")
    print(f"Final aggregation state: {aggregation_states[-1]}")
    print(f"Maximum phosphorylation count: {max(phospho_counts)} out of {max(phospho_counts) if phospho_counts else 0} sites")


def plot_site_probabilities(site_probabilities, timepoints):
    """
    Plots the probability trajectory for each phosphorylation site as a line plot.
    """
    plt.figure(figsize=(12, 6))
    # Plot each site's probability trajectory
    for site, probs in site_probabilities.items():
        plt.plot(timepoints, probs, label=f'Site {site}')
    plt.xlabel('Time Step')
    plt.ylabel('Phosphorylation Probability')
    plt.title('Site-specific Phosphorylation Probability Over Time')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small', ncol=2)
    plt.tight_layout()
    plt.show()


def plot_phosphorylation_heatmap(site_probabilities, timepoints):
    """
    Plots a heatmap of phosphorylation probabilities for all sites over time.
    """
    import seaborn as sns
    import pandas as pd
    # Convert the site probabilities to a DataFrame for heatmap plotting
    df = pd.DataFrame(site_probabilities, index=timepoints)
    plt.figure(figsize=(14, 6))
    sns.heatmap(df.T, cmap='viridis', cbar_kws={'label': 'Phosphorylation Probability'})
    plt.xlabel('Time Step')
    plt.ylabel('Site')
    plt.title('Phosphorylation Probability Heatmap (Sites x Time)')
    plt.tight_layout()
    plt.show()


def run_and_plot_simulation(env=None, plot_sites=False, plot_heatmap=False):
    """
    Runs the tau protein simulation and generates summary and (optionally) site-specific visualizations.
    Args:
        env: Environment object (optional). If None, uses default parameters.
        plot_sites: If True, plot all site probability trajectories.
        plot_heatmap: If True, plot a heatmap of all site probabilities.
    """
    # Set up environment if not provided
    if env is None:
        env = Environment(temperature=39, kinase_level=1.5, oxidative_stress=0.2)
    # Initialize tau protein and run simulation
    tau = TauProtein(isoform='4R')
    timepoints = np.arange(100)
    site_probabilities = tau.update_state(env, timepoints)
    # Main summary plot
    plot_tau_summary(tau.history)
    # Optional: plot all site probability trajectories
    if plot_sites:
        plot_site_probabilities(site_probabilities, timepoints)
    # Optional: plot heatmap of all site probabilities
    if plot_heatmap:
        plot_phosphorylation_heatmap(site_probabilities, timepoints)
    # Print tau protein object summary
    print(tau)

# If run as a script, show the main summary plot and advanced visualizations
if __name__ == "__main__":
    # By default, show both site-specific and heatmap visualizations
    run_and_plot_simulation(plot_sites=True, plot_heatmap=True) 