import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from scipy.optimize import curve_fit  # Importer la fonction pour l'ajustement

# Fonction gaussienne
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def main(argv):
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Mixing layer analysis with velocity profile.')
    
    # Adding arguments for the calculation method
    parser.add_argument('--du', action='store_true', help='Use the velocity gradient method (du/dy) to calculate the mixing layer.')
    parser.add_argument('--u', action='store_true', help='Use the method based on max and min values of u to calculate the mixing layer.')
    
    # Parse the arguments
    args = parser.parse_args(argv)

    # Determine which method to use based on the command-line arguments
    if args.du:
        arg = 2  # Gradient-based method
    elif args.u:
        arg = 1  # Max/min-based method
    else:
        print("Error: you must specify a method with --du or --u")
        return

    # Load data from the specified file
    data_file = 'outputs/profil_y.dat'
    if not os.path.exists(data_file):
        raise FileNotFoundError(f"The file {data_file} was not found.")
    
    data = np.loadtxt(data_file)
    
    # Check the dimensions of the data
    if data.shape[1] < 5:
        raise ValueError("The data must contain at least 5 columns.")
    
    # Extract columns
    y = data[:, 0]
    u = data[:, 1]
    v = data[:, 2]
    w = data[:, 3]
    p = data[:, 4]
    
    if (arg == 1):
        # Calculate the mixing layer thickness based on u max/min
        u_max = np.max(u)
        u_min = np.min(u)
        
        # Define thresholds
        lower_threshold = 0.99 * u_min
        upper_threshold = 0.99 * u_max
        
        mix_layer_indices = np.where((u >= lower_threshold) & (u <= upper_threshold))[0]
    elif (arg == 2):
        # Calculate the velocity gradient du/dy
        du_dy = np.gradient(u, y)
        
        # Use a threshold value to identify regions with significant velocity variation
        gradient_threshold = 0.01 * np.max(np.abs(du_dy))  # Threshold based on 1% of the maximum gradient
        mix_layer_indices = np.where(np.abs(du_dy) > gradient_threshold)[0]

    if len(mix_layer_indices) > 0:
        mix_layer_start = y[mix_layer_indices[0]]
        mix_layer_end = y[mix_layer_indices[-1]]
        mix_layer_thickness = mix_layer_end - mix_layer_start
    else:
        mix_layer_thickness = 0
        mix_layer_start = mix_layer_end = None  # No mixing layer identified

    # Calculate the shear stress (tau) based on the velocity gradient
    # Assuming constant viscosity mu; thus, tau = mu * (du/dy)
    du_dy = np.gradient(u, y)  # Calculate du/dy
    mu = 1  # Assumed constant viscosity (adjust as necessary)
    shear_stress = mu * du_dy  # Calculate shear stress

    # Fit a Gaussian curve to the shear stress data
    initial_guess = [np.max(shear_stress), y[np.argmax(shear_stress)], 1]  # Initial guess: [amplitude, mean, stddev]
    popt, pcov = curve_fit(gaussian, y, shear_stress, p0=initial_guess)
    
    # Display the mixing layer thickness
    if mix_layer_thickness > 0:
        shear_stress_max = np.max(shear_stress)
        print(f"Mixing layer thickness: {mix_layer_thickness:.2f} units")
        print(f"Shear Stress max: {popt[0]:.2f}")
        print(f"Gaussian fit parameters: amplitude={popt[0]:.2f}, mean={popt[1]:.2f}, stddev={popt[2]:.2f}")
    else:
        print("No identifiable mixing layer.")

    # Create a figure with 3 subplots side by side
    plt.figure(figsize=(20, 9))

    colors = {
        'u': '#0072B2',
        'v': '#009E73',
        'w': '#D55E00',
        'shear_stress': '#FF6347',
        'gaussian_fit': '#999999',
        'mix_layer_start': '#E69F00',
        'mix_layer_end': '#CC79A7'
    }

    # Plot for u
    ax1 = plt.subplot(1, 3, 1)
    ax1.plot(u, y, label='u', color=colors['u'])
    if mix_layer_start is not None:
        ax1.axhline(mix_layer_start, color=colors['mix_layer_start'], linestyle='--', label='Mixing Layer Start')
        ax1.axhline(mix_layer_end, color=colors['mix_layer_end'], linestyle='--', label='Mixing Layer End')
    ax1.set_title('Velocity Profile u')
    ax1.set_xlabel('Velocity (u) or Shear Stress (tau)')
    ax1.set_ylabel('Height (y)')
    ax1.legend(loc='upper left')
    ax1.grid()

    # Add shear stress and gaussian fit on the same plot
    ax2 = ax1.twinx()  # Create a second y-axis for shear stress
    ax2.plot(shear_stress, y, label='Shear Stress', color=colors['shear_stress'], linestyle='--')
    ax2.plot(gaussian(y, *popt), y, label='Gaussian Fit', color=colors['gaussian_fit'], linestyle='--')
    ax2.set_ylabel('')  # Remove the y-label
    ax2.set_yticks([])  # Remove y-ticks
    ax2.legend(loc='upper right')  # Legend for shear stress and Gaussian fit
    ax2.grid(False)  # Turn off the grid for the second axis

    # Plot for v
    ax3 = plt.subplot(1, 3, 2)
    ax3.plot(v, y, label='v', color=colors['v'])
    if mix_layer_start is not None:
        ax3.axhline(mix_layer_start, color=colors['mix_layer_start'], linestyle='--')
        ax3.axhline(mix_layer_end, color=colors['mix_layer_end'], linestyle='--')
    ax3.set_title('Velocity Profile v')
    ax3.set_xlabel('Velocity (v)')
    ax3.set_ylabel('')  # Remove the y-label
    ax3.legend(loc='upper left')
    ax3.grid()

    # Plot for w
    ax4 = plt.subplot(1, 3, 3)
    ax4.plot(w, y, label='w', color=colors['w'])
    if mix_layer_start is not None:
        ax4.axhline(mix_layer_start, color=colors['mix_layer_start'], linestyle='--')
        ax4.axhline(mix_layer_end, color=colors['mix_layer_end'], linestyle='--')
    ax4.set_title('Velocity Profile w')
    ax4.set_xlabel('Velocity (w)')
    ax4.set_ylabel('')  # Remove the y-label
    ax4.legend(loc='upper left')
    ax4.grid()

    plt.tight_layout()
    plt.show()

# Entry point of the script
if __name__ == "__main__":
    import sys
    main(sys.argv[1:])


