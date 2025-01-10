import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Specify the NetCDF file name
file_name = "saved.nc"
output_plot_file = "position_plot.png"  # Specify the name for the saved plot file

# Initialize variables to store the arrays
position = E_x = E_y = potential = velocity =acceleration=rho= None

# Open the NetCDF file
try:
    dataset = nc.Dataset(file_name, mode='r')
    print(f"File '{file_name}' opened successfully.\n")

    # Extract variables individually
    if 'position' in dataset.variables:
        position = dataset.variables['position'][:]
        print(f"Extracted 'position': shape {position.shape}")
    else:
        print("Variable 'position' not found.")

    if 'E_x' in dataset.variables:
        E_x = dataset.variables['E_x'][:]
        print(f"Extracted 'E_x': shape {E_x.shape}")
    else:
        print("Variable 'E_x' not found.")

    if 'E_y' in dataset.variables:
        E_y = dataset.variables['E_y'][:]
        print(f"Extracted 'E_y': shape {E_y.shape}")
    else:
        print("Variable 'E_y' not found.")

    if 'potential' in dataset.variables:
        potential = dataset.variables['potential'][:]
        print(f"Extracted 'potential': shape {potential.shape}")
    else:
        print("Variable 'potential' not found.")
    
    if 'velocity' in dataset.variables:
        velocity = dataset.variables['velocity'][:]
        print(f"Extracted 'velocity': shape {velocity.shape}")
    else:
        print("Variable 'velocity' not found.")

    if 'acceleration' in dataset.variables:
        acceleration = dataset.variables['acceleration'][:]
        print(f"Extracted 'acceleration': shape {acceleration.shape}")
    else:
        print("Variable 'acceleration' not found.")
    
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
    
    # All the variables are here now
    # Plot position using the first axis vs second axis
    if position is not None:
        # Assuming 'position' is a 2D array and we are plotting the first axis vs the second axis.
        x = position[:, 0]  # Taking the first column (axis 1)
        y = position[:, 1]  # Taking the second column (axis 2)

        plt.figure(figsize=(8, 6))
        plt.plot(x, y, marker='o', linestyle='-', color='b')
        plt.title("Position Plot")
        plt.xlabel("Position (Axis 1)")
        plt.ylabel("Position (Axis 2)")
        plt.grid(True)

        # Save the plot to a file instead of showing it
        plt.savefig(output_plot_file)
        print(f"Plot saved as '{output_plot_file}'.")

except FileNotFoundError:
    print(f"Error: File '{file_name}' not found.")
except Exception as e:
    print(f"An error occurred: {e}")
