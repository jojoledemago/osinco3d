import numpy as np
import argparse
import matplotlib.pyplot as plt

def read_fields(filename):
    """Read fields from a binary file."""
    with open(filename, 'rb') as file:
        # Read time
        time = np.fromfile(file, dtype=np.float64, count=1)[0]

        # Read dimensions nx, ny, nz
        nx, ny, nz = np.fromfile(file, dtype=np.int32, count=3)

        # Read coordinates x, y, z
        x = np.fromfile(file, dtype=np.float64, count=nx)
        y = np.fromfile(file, dtype=np.float64, count=ny)
        z = np.fromfile(file, dtype=np.float64, count=nz)

        # Read velocity components ux, uy, uz
        ux = np.fromfile(file, dtype=np.float64, count=nx * ny * nz).reshape((nx, ny, nz), order='F')
        uy = np.fromfile(file, dtype=np.float64, count=nx * ny * nz).reshape((nx, ny, nz), order='F')
        uz = np.fromfile(file, dtype=np.float64, count=nx * ny * nz).reshape((nx, ny, nz), order='F')

        # Read pressure component pp
        pp = np.fromfile(file, dtype=np.float64, count=nx * ny * nz).reshape((nx, ny, nz), order='F')

    return time, nx, ny, nz, x, y, z, ux, uy, uz, pp

def compute_average_velocity(ux, uy, uz, x_index):
    """Compute the average velocity profiles along y for a given x index."""
    avg_ux = np.mean(ux[x_index, :, :], axis=1)  # Average over y
    avg_uy = np.mean(uy[x_index, :, :], axis=1)
    avg_uz = np.mean(uz[x_index, :, :], axis=1)
    return avg_ux, avg_uy, avg_uz

def compute_mix_layer_limits(avg_ux, y):
    """Compute the lower and upper limits of the mixing layer using min-max and gradient methods."""
    u_min = np.min(avg_ux)
    u_max = np.max(avg_ux)
    
    # Min-Max Method
    target_low = 0.99 * u_min
    target_high = 0.99 * u_max
    
    # Find indices for min-max method
    idx_low = np.where(avg_ux >= target_low)[0][0]  # First index where avg_ux >= target_low
    idx_high = np.where(avg_ux <= target_high)[0][-1]  # Last index where avg_ux <= target_high
    lower_limit_min_max = y[idx_low]
    upper_limit_min_max = y[idx_high]
    
    # Gradient Method
    du_dy = np.gradient(avg_ux, y)
    max_gradient = np.max(np.abs(du_dy))
    gradient_threshold = 0.01 * max_gradient

    # Find indices for gradient method
    idx_gradient = np.where(np.abs(du_dy) > gradient_threshold)[0]
    if len(idx_gradient) > 0:
        lower_limit_gradient = y[idx_gradient[0]]
        upper_limit_gradient = y[idx_gradient[-1]]
    else:
        lower_limit_gradient = upper_limit_gradient = None  # No layer detected

    return (lower_limit_min_max, upper_limit_min_max), (lower_limit_gradient, upper_limit_gradient)

def plot_velocity_profiles(y, avg_ux, avg_uy, avg_uz, viscosity, plot_shear_stress):
    """Plot the average velocity profiles and their derivatives."""
    # Calculate the derivatives using np.gradient
    d_avg_ux = np.gradient(avg_ux, y)
    d_avg_uy = np.gradient(avg_uy, y)
    d_avg_uz = np.gradient(avg_uz, y)

    # Calculate shear stress
    shear_stress = viscosity * d_avg_ux  # Assuming d_avg_ux is du/dy

    # Calculate mixing layer limits
    (lower_limit_min_max, upper_limit_min_max), (lower_limit_gradient, upper_limit_gradient) = compute_mix_layer_limits(avg_ux, y)

    # Define colors suitable for color-blind users
    color_ux = '#E69F00'  # orange
    color_uy = '#56B4E9'  # sky blue
    color_uz = '#009E73'  # teal
    color_shear = '#D55E00'  # reddish orange for shear stress
    color_thickness = '#CC79A7'  # pink for mixing layer limits

    plt.figure(figsize=(12, 10))

    # Plotting the average velocity profile for u_x
    plt.subplot(1, 3, 1)
    plt.plot(avg_ux, y, label='$\\langle u_x \\rangle$', color=color_ux)
    if plot_shear_stress:
        plt.plot(shear_stress, y, label='Shear Stress $\\tau$', color=color_shear, linestyle='--')

    # Plotting the limits of the mixing layer
    plt.axhline(y=lower_limit_min_max, color=color_thickness, linestyle=':', label='Mixing Layer Lower Limit (Min-Max)')
    plt.axhline(y=upper_limit_min_max, color=color_thickness, linestyle=':')
    if lower_limit_gradient is not None and upper_limit_gradient is not None:
        plt.axhline(y=lower_limit_gradient, color=color_thickness, linestyle='-.', label='Mixing Layer Lower Limit (Gradient)')
        plt.axhline(y=upper_limit_gradient, color=color_thickness, linestyle='-.')

    plt.title('Velocity Profile $\\langle u_x \\rangle$ and Shear Stress')
    plt.xlabel('Velocity ($u_x$) / Shear Stress ($\\tau$)')
    plt.ylabel('Height (y)')
    plt.ylim([min(y), max(y)])  # Set y-axis limits
    plt.grid()
    plt.legend()

    # Plotting the average velocity profiles for u_y and u_z
    plt.subplot(1, 3, 2)
    plt.plot(avg_uy, y, label='$\\langle u_y \\rangle$', color=color_uy, linestyle='--', linewidth=1.5)
    plt.plot(avg_uz, y, label='$\\langle u_z \\rangle$', color=color_uz, linestyle='-.', linewidth=1.5)
    # Plotting the limits of the mixing layer
    plt.axhline(y=lower_limit_min_max, color=color_thickness, linestyle=':')
    plt.axhline(y=upper_limit_min_max, color=color_thickness, linestyle=':')
    if lower_limit_gradient is not None and upper_limit_gradient is not None:
        plt.axhline(y=lower_limit_gradient, color=color_thickness, linestyle='-.')
        plt.axhline(y=upper_limit_gradient, color=color_thickness, linestyle='-.')
    plt.title('Velocity Profiles $\\langle u_y \\rangle$ and $\\langle u_z \\rangle$')
    plt.xlabel('Velocity')
    plt.ylabel('Height (y)')
    plt.ylim([min(y), max(y)])  # Set y-axis limits
    plt.grid()
    plt.legend()

    # Plotting the derivatives of the velocity profiles
    plt.subplot(1, 3, 3)
    plt.plot(d_avg_ux, y, label='$\\frac{d(\\langle u_x \\rangle)}{dy}$', color=color_ux, linestyle='--', linewidth=1.5)
    plt.plot(d_avg_uy, y, label='$\\frac{d(\\langle u_y \\rangle)}{dy}$', color=color_uy, linestyle='--', linewidth=1.5)
    plt.plot(d_avg_uz, y, label='$\\frac{d(\\langle u_z \\rangle)}{dy}$', color=color_uz, linestyle='--', linewidth=1.5)
    plt.axhline(y=lower_limit_min_max, color=color_thickness, linestyle=':')
    plt.axhline(y=upper_limit_min_max, color=color_thickness, linestyle=':')
    if lower_limit_gradient is not None and upper_limit_gradient is not None:
        plt.axhline(y=lower_limit_gradient, color=color_thickness, linestyle='-.')
        plt.axhline(y=upper_limit_gradient, color=color_thickness, linestyle='-.')
    plt.title('Derivatives of Velocity Profiles')
    plt.xlabel('Derivative')
    plt.ylabel('Height (y)')
    plt.ylim([min(y), max(y)])  # Set y-axis limits
    plt.grid()
    plt.legend()

    plt.tight_layout()
    plt.show()

# Example usage
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Read simulation fields from a binary file.')
    parser.add_argument('--filename', '-f', type=str, default='fields.bin',
                        help='The path to the binary file to read')
    parser.add_argument('--x_index', '-x', type=int, default=0,
                        help='The x index for which to compute velocity profiles (default is 0)')
    parser.add_argument('--viscosity', '-v', type=float, default=1.0,
                        help='Dynamic viscosity of the fluid (default is 1.0)')
    parser.add_argument('--shearstress', '-s', action='store_true',
                        help='Plot the shear stress on the velocity profile')
    args = parser.parse_args()

    # Read the fields
    time, nx, ny, nz, x, y, z, ux, uy, uz, pp = read_fields(args.filename)

    # Check index validity
    if args.x_index < 0 or args.x_index >= nx:
        raise ValueError(f"x_index ({args.x_index}) is out of bounds for nx ({nx}).")

    # Compute average velocity profiles
    avg_ux, avg_uy, avg_uz = compute_average_velocity(ux, uy, uz, args.x_index)

    # Plot the profiles
    plot_velocity_profiles(y, avg_ux, avg_uy, avg_uz, args.viscosity, args.shearstress)

