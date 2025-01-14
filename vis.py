# @author Facundo Costa with additions by Gianluca Seaford

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm  # For color mapping
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from os import path, makedirs

file_name = 'Gauss-Seidel_Electrostatics.nc'

plots_dir = 'plots'
output_plot_file = path.join(plots_dir, "Particle_Positions.png")  # Specify the name for the saved plot file
output_ex_plot_file = path.join(plots_dir, "E_x_Field_Intensity.png")  # Specify the name for the saved E_x color plot file
output_ey_plot_file = path.join(plots_dir, "E_y_Field_Intensity.png")
output_rho_plot_file = path.join(plots_dir, "Charge_Density_Pseudocolour.png")
output_pot_plot_file = path.join(plots_dir, "Potential_Pseudocolour.png")
output_animation_file = path.join(plots_dir, "Particle_Motion_Simulation.mp4")  # Output animation file

# Initialize variables to store the arrays
position = E_x = E_y = potential = velocity = acceleration = rho = None

# Open the NetCDF file
try:
    dataset = nc.Dataset(file_name, mode='r')
    print(f"File '{file_name}' opened successfully.\n")

    # Extract variables individually
    if 'positions' in dataset.variables:
        position = dataset.variables['positions'][:]
        print(f"Extracted 'positions': shape {position.shape}")
    else:
        print("Variable 'positions' not found.")

    if 'Ex' in dataset.variables:
        E_x = dataset.variables['Ex'][:]
        print(f"Extracted 'Ex': shape {E_x.shape}")
    else:
        print("Variable 'Ex' not found.")

    if 'Ey' in dataset.variables:
        E_y = dataset.variables['Ey'][:]
        print(f"Extracted 'Ey': shape {E_y.shape}")
    else:
        print("Variable 'Ey' not found.")

    if 'phi' in dataset.variables:
        potential = dataset.variables['phi'][:]
        print(f"Extracted 'phi': shape {potential.shape}")
    else:
        print("Variable 'phi' not found.")
    
    if 'velocities' in dataset.variables:
        velocity = dataset.variables['velocities'][:]
        print(f"Extracted 'velocities': shape {velocity.shape}")
    else:
        print("Variable 'velocities' not found.")

    if 'accelerations' in dataset.variables:
        acceleration = dataset.variables['accelerations'][:]
        print(f"Extracted 'accelerations': shape {acceleration.shape}")
    else:
        print("Variable 'accelerations' not found.")
    
    if 'rho' in dataset.variables:
        rho = dataset.variables['rho'][:]
        print(f"Extracted 'rho': shape {rho.shape}")
    else:
        print("Variable 'rho' not found.")
    


    # Close the dataset
    dataset.close()
    print("\nFile closed successfully.")

    # Access and print details of extracted variables
    print("\nAccessing extracted variables:")
    if position is not None:
        print(f"Position array: shape {position.shape}")
    if E_x is not None:
        print(f"E_x array: shape {E_x.shape}")
    if E_y is not None:
        print(f"E_y array: shape {E_y.shape}")
    if potential is not None:
        print(f"Potential array: shape {potential.shape}")
    if velocity is not None:
        print(f"Velocity array: shape {velocity.shape}")
    if acceleration is not None:
        print(f"acceleration array: shape {acceleration.shape}")
    if rho is not None:
        print(f"rho array: shape {rho.shape}")
        
    makedirs(plots_dir, exist_ok=True)
    
    # All the variables are here now for our use
    if position is not None:
        x = position[:, 0]
        y = position[:, 1]
        num_points = len(x)

        fig, ax = plt.subplots(figsize=(8, 6))
        scatter_main = ax.scatter(x, y, c=range(num_points), cmap='viridis', s=20)
        fig.colorbar(scatter_main, ax=ax, label="Iteration Number")

        ax.set_title("Particle Positions across Grid")
        ax.set_xlabel("X Position (arb.)")
        ax.set_ylabel("Y Position (arb.)")
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.grid(True)

        # Create the inset axis
        inset_ax = inset_axes(
            ax,
            width="30%",   
            height="30%",
            loc="lower left",
            borderpad=4.0
        )

        # Plot the same data in the inset
        scatter_inset = inset_ax.scatter(x, y, c=range(num_points), cmap='viridis', s=20)

        # Zoom in on the min/max region
        inset_ax.set_xlim(x.min(), x.max())
        inset_ax.set_ylim(y.min(), y.max())

        # Give the inset its own axis labels
        inset_ax.set_xlabel("X (arb.)", fontsize=8)
        inset_ax.set_ylabel("Y (arb.)", fontsize=8)
        inset_ax.tick_params(axis='both', which='major', labelsize=8)

        plt.savefig(output_plot_file, dpi=300)
        print(f"Position scatter plot saved as '{output_plot_file}'.")
        
    if E_x is not None:
        plt.figure(figsize=(8, 6))
        plt.imshow(E_x, cmap='viridis', origin='lower', extent=(-1, 1, -1, 1))
        plt.colorbar(label="Field Strength (arb.)")
        plt.title("E_x Field Intensity across Grid.")
        plt.xlabel("X Position (arb.)")
        plt.ylabel("Y Position (arb.)")
        plt.grid(False)

        # Save the color plot
        plt.savefig(output_ex_plot_file)
        print(f"E_x field color plot saved as '{output_ex_plot_file}'.")
    else:
        print("E_x data not available for plotting.")
        
    if E_y is not None:
        plt.figure(figsize=(8, 6))
        plt.imshow(E_y, cmap='viridis', origin='lower', extent=(-1, 1, -1, 1))
        plt.colorbar(label="Field Strength (arb.)")
        plt.title("E_y Field Intensity across Grid.")
        plt.xlabel("X Position (arb.)")
        plt.ylabel("Y Position (arb.)")
        plt.grid(False)

        # Save the color plot
        plt.savefig(output_ey_plot_file)
        print(f"E_x field color plot saved as '{output_ey_plot_file}'.")
    else:
        print("E_x data not available for plotting.")
        
    if rho is not None:
        plt.figure(figsize=(8, 6))
        plt.imshow(rho, cmap='viridis', origin='lower', extent=(-1, 1, -1, 1))
        plt.colorbar(label="Charge Density (arb.)")
        plt.title("Distribution of Charge Density across Grid.")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.grid(False)

        # Save the color plot
        plt.savefig(output_rho_plot_file)
        print(f"rho field color plot saved as '{output_rho_plot_file}'.")
    else:
        print("rho data not available for plotting.")

    if potential is not None:
        plt.figure(figsize=(8, 6))
        plt.imshow(potential, cmap='viridis', origin='lower', extent=(-1, 1, -1, 1))
        plt.colorbar(label="rho")
        plt.title(" Gauss-Seidel Potential across Grid.")
        plt.xlabel("X Position (arb.)")
        plt.ylabel("Y Position (arb.)")
        plt.grid(False)

        # Save the color plot
        plt.savefig(output_pot_plot_file)
        print(f"potential color plot saved as '{output_pot_plot_file}'.")
    else:
        print("potential data not available for plotting.")

except FileNotFoundError:
    print(f"Error: File '{file_name}' not found.")
except Exception as e:
    print(f"An error occurred: {e}")
