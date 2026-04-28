from os import path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def visualize_transverse_sweep(csv_path):
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found.")
        return

    # Load the exported pathlengths
    df = pd.read_csv(csv_path)
    
    # Ensure wavelength is treated as a categorical variable for clean legends
    df['wavelength_nm'] = df['wavelength_nm'].astype(str) + ' nm'

    # Set up a 3-panel plot
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    
    # --- 1. Detected Intensity ---
    # Using raw string (r) for the title
    sns.lineplot(data=df, x='x_location_cm', y='i_detected', 
                 hue='wavelength_nm', marker='o', ax=axes[0], 
                 palette=['#d62728', '#8c564b', '#7f7f7f'])
    axes[0].set_title(r'1. Total Detected Intensity ($I_{det}$) vs. Array Position', fontsize=14, fontweight='bold')
    axes[0].set_ylabel('Intensity (per launched photon)')
    
    # --- 2. Vessel Mean Pathlength (<Lv>) ---
    # Using raw strings (r) for title and ylabel
    sns.lineplot(data=df, x='x_location_cm', y='L_vessel_cm', 
                 hue='wavelength_nm', marker='s', ax=axes[1],
                 palette=['#d62728', '#8c564b', '#7f7f7f'])
    axes[1].set_title(r'2. Vessel Mean Pathlength ($\langle L_v \rangle$) vs. Array Position', fontsize=14, fontweight='bold')
    axes[1].set_ylabel(r'$\langle L_v \rangle$ (cm)')
    
    # --- 3. Muscle Mean Pathlength (<Lm>) ---
    # Using raw strings (r) for title and ylabel
    sns.lineplot(data=df, x='x_location_cm', y='L_muscle_cm', 
                 hue='wavelength_nm', marker='^', ax=axes[2],
                 palette=['#d62728', '#8c564b', '#7f7f7f'])
    axes[2].set_title(r'3. Bulk Muscle Mean Pathlength ($\langle L_m \rangle$) vs. Array Position', fontsize=14, fontweight='bold')
    axes[2].set_ylabel(r'$\langle L_m \rangle$ (cm)')
    axes[2].set_xlabel('Sensor X-Location (cm)  [0 = Directly over Vessel]', fontsize=12)

    # Clean up legends
    for ax in axes:
        ax.legend(title="Wavelength")

    plt.tight_layout()
    
    # Save the figure
    config = "s75_r0p0575mm"
    save_path = f"reports/results/{config}_transverse_sweep_visualization.png"
    plt.savefig(save_path, dpi=300)
    print(f"Visualization saved to {save_path}")
    plt.show()

# Run the visualizer
config = "s75_r0p0575mm"
path = f"reports/results/{config}_transverse_sweep_pathlengths.csv"
visualize_transverse_sweep(path)