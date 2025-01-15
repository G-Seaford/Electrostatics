import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm  # For color mapping
# By facu

file_name = 'Gauss-Seidel_Electrostatics.nc'
output_plot_file = "position_plot.png"  # Specify the name for the saved plot file
output_ex_plot_file = "Ex_color_plot.png"  # Specify the name for the saved E_x color plot file
output_ey_plot_file = "Ey_color_plot.png"
output_rho_plot_file = "rho_color_plot.png"
output_pot_plot_file = "pot_color_plot.png"
output_animation_file = "particle_motion_simulation.mp4"  # Output animation file
# Initialize variables to store the arrays
position = E_x = E_y = potential = velocity = accelerations = rho = None

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
    
    # All the variables are here now for our use
    if position is not None:
        x = position[:, 0]
        y = position[:, 1]

        # Color mapping by index
        num_points = len(x)
        colors = cm.jet(np.linspace(0, 1, num_points))  # Gradient color map

        plt.figure(figsize=(8, 6))
        plt.scatter(x, y, c=range(num_points), cmap='viridis', s=20)
        plt.colorbar(label="# iteration")
        plt.title("Particle Positions")
        plt.xlabel("Position (Axis 1)")
        plt.ylabel("Position (Axis 2)")
        plt.xlim(-1, 1)
        plt.ylim(-1, 1)
        plt.grid(True)

        # Save the plot
        plt.savefig(output_plot_file)
        print(f"Position scatter plot saved as '{output_plot_file}'.")
    if E_x is not None:
        plt.figure(figsize=(8, 6))
        plt.imshow(E_x, cmap='viridis', origin='lower', extent=(-1, 1, -1, 1))
        plt.colorbar(label="E_x Field")
        plt.title("E_x Field Color Plot")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.grid(False)

        # Save the color plot
        plt.savefig(output_ex_plot_file)
        print(f"E_x field color plot saved as '{output_ex_plot_file}'.")
    else:
        print("E_x data not available for plotting.")
    if E_y is not None:
        plt.figure(figsize=(8, 6))
        plt.imshow(E_y, cmap='viridis', origin='lower', extent=(-1, 1, -1, 1))
        plt.colorbar(label="E_y Field")
        plt.title("E_y Field Color Plot")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.grid(False)

        # Save the color plot
        plt.savefig(output_ey_plot_file)
        print(f"E_x field color plot saved as '{output_ey_plot_file}'.")
    else:
        print("E_x data not available for plotting.")
    if rho is not None:
        plt.figure(figsize=(8, 6))
        plt.imshow(rho, cmap='viridis', origin='lower', extent=(-1, 1, -1, 1))
        plt.colorbar(label="rho")
        plt.title("rho Plot")
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
        plt.title("potential Plot")
        plt.xlabel("X")
        plt.ylabel("Y")
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
