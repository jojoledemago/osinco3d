import numpy as np
import argparse
import matplotlib.pyplot as plt
import sys

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

def compute_fluctuations(ux, uy, uz, avg_ux, avg_uy, avg_uz, x_index):
    """Compute the fluctuation components of velocity at a given x index."""
    ux_fluct = ux[x_index, :, :] - avg_ux[:, np.newaxis]
    uy_fluct = uy[x_index, :, :] - avg_uy[:, np.newaxis]
    uz_fluct = uz[x_index, :, :] - avg_uz[:, np.newaxis]
    return ux_fluct, uy_fluct, uz_fluct

def compute_mix_layer_limits(avg_ux, y, min_max_threshold, gradient_threshold_factor):
    """Compute the lower and upper limits of the mixing layer using min-max and gradient methods."""
    u_min = avg_ux[ 0]
    u_max = avg_ux[-1]
    
    # Min-Max Method
    target_low = min_max_threshold * u_min
    target_high = min_max_threshold * u_max
    
    # Find indices for min-max method
    idx_low = np.where(avg_ux >= target_low)[0][0]  # First index where avg_ux >= target_low
    idx_high = np.where(avg_ux <= target_high)[0][-1]  # Last index where avg_ux <= target_high
    lower_limit_min_max = y[idx_low]
    upper_limit_min_max = y[idx_high]
    
    # Gradient Method
    du_dy = np.gradient(avg_ux, y)
    max_gradient = np.max(np.abs(du_dy))
    gradient_threshold = gradient_threshold_factor * max_gradient

    # Find indices for gradient method
    idx_gradient = np.where(np.abs(du_dy) > gradient_threshold)[0]
    if len(idx_gradient) > 0:
        lower_limit_gradient = y[idx_gradient[0]]
        upper_limit_gradient = y[idx_gradient[-1]]
    else:
        lower_limit_gradient = upper_limit_gradient = None  # No layer detected

    return (lower_limit_min_max, upper_limit_min_max), (lower_limit_gradient, upper_limit_gradient)

def plot_velocity_profiles(y, avg_ux, avg_uy, avg_uz, viscosity, 
        plot_shear_stress, min_max_limits, gradient_limits, x, x_index):
    """Plot the average velocity profiles and their derivatives."""
    # Calculate the derivatives using np.gradient
    d_avg_ux = np.gradient(avg_ux, y)
    d_avg_uy = np.gradient(avg_uy, y)
    d_avg_uz = np.gradient(avg_uz, y)

    # Calculate shear stress
    shear_stress = viscosity * d_avg_ux  # Assuming d_avg_ux is du/dy

    # Define colors suitable for color-blind users
    color_ux = '#E69F00'  # orange
    color_uy = '#56B4E9'  # sky blue
    color_uz = '#009E73'  # teal
    color_shear = '#D55E00'  # reddish orange for shear stress
    color_thickness = '#CC79A7'  # pink for mixing layer limits

    plt.figure(figsize=(12, 9))
    # Set the overall title for the figure, including x[x_index]
    plt.suptitle(f'Velocity Profiles and Shear Stress for x = {x[x_index]:.2f}', fontsize=14)

    # Plotting the average velocity profile for u_x
    plt.subplot(1, 3, 1)
    plt.plot(avg_ux, y, label='$\\langle u_x \\rangle$', color=color_ux)
    if plot_shear_stress:
        plt.plot(shear_stress, y, label='Shear Stress $\\tau$', color=color_shear, linestyle='--')

    # Plotting the limits of the mixing layer
    lower_limit_min_max, upper_limit_min_max = min_max_limits
    plt.axhline(y=lower_limit_min_max, color=color_thickness, linestyle=':', label='Mixing Layer Limit (Min-Max)')
    plt.axhline(y=upper_limit_min_max, color=color_thickness, linestyle=':')

    lower_limit_gradient, upper_limit_gradient = gradient_limits
    if lower_limit_gradient is not None and upper_limit_gradient is not None:
        plt.axhline(y=lower_limit_gradient, color=color_thickness, linestyle='-.', label='Mixing Layer Limit (Gradient)')
        plt.axhline(y=upper_limit_gradient, color=color_thickness, linestyle='-.')

    if plot_shear_stress:
        plt.title('Velocity Profile $\\langle u_x \\rangle$ and Shear Stress')
    else:
        plt.title('Velocity Profile $\\langle u_x \\rangle$')
    plt.xlabel('Velocity ($u_x$) / Shear Stress ($\\tau$)')
    plt.ylabel('Height (y)')
    plt.ylim([min(y), max(y)])  # Set y-axis limits
    plt.grid()
    plt.legend()

    # Plotting the average velocity profiles for u_y and u_z
    plt.subplot(1, 3, 2)
    plt.plot(avg_uy, y, label='$\\langle u_y \\rangle$', color=color_uy, linestyle='--', linewidth=1.5)
    plt.plot(avg_uz, y, label='$\\langle u_z \\rangle$', color=color_uz, linestyle='-.', linewidth=1.5)
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
    plt.plot(d_avg_ux, y, label='$d\\langle u_x \\rangle/dy$', color=color_ux)
    plt.plot(d_avg_uy, y, label='$d\\langle u_y \\rangle/dy$', color=color_uy, linestyle='--', linewidth=1.5)
    plt.plot(d_avg_uz, y, label='$d\\langle u_z \\rangle/dy$', color=color_uz, linestyle='-.', linewidth=1.5)
    plt.axhline(y=lower_limit_min_max, color=color_thickness, linestyle=':')
    plt.axhline(y=upper_limit_min_max, color=color_thickness, linestyle=':')
    if lower_limit_gradient is not None and upper_limit_gradient is not None:
        plt.axhline(y=lower_limit_gradient, color=color_thickness, linestyle='-.')
        plt.axhline(y=upper_limit_gradient, color=color_thickness, linestyle='-.')

    plt.title('Velocity Gradients')
    plt.xlabel('Velocity Gradient')
    plt.ylabel('Height (y)')
    plt.ylim([min(y), max(y)])  # Set y-axis limits
    plt.grid()
    plt.legend()

    plt.tight_layout()
    plt.show()

def plot_fluctuations(ux_fluct, uy_fluct, uz_fluct, x, y, z, x_index):
    """Plot the velocity fluctuations in 3D with a single colormap and shared colorbar."""
    fig = plt.figure(figsize=(18, 6))
    # Set the overall title for the figure, including x[x_index]
    fig.suptitle(f'Velocity Fluctuations for x = {x[x_index]:.2f}', fontsize=14)

    # Create a meshgrid for y and z
    Y, Z = np.meshgrid(y, z, indexing='ij')

    # Find the global minimum and maximum values for the color scale
    vmin = min(np.min(ux_fluct), np.min(uy_fluct), np.min(uz_fluct))
    vmax = max(np.max(ux_fluct), np.max(uy_fluct), np.max(uz_fluct))

    # Fluctuations for u_x
    ax1 = fig.add_subplot(131, projection='3d')
    surf_ux = ax1.plot_surface(Y, Z, ux_fluct, cmap='viridis', vmin=vmin, vmax=vmax)
    ax1.set_title("Fluctuations $u'_x$")
    ax1.set_xlabel("Height (y)")
    ax1.set_ylabel("Span (z)")
    ax1.set_zlabel("Fluctuations $u'_x$")

    # Fluctuations for u_y
    ax2 = fig.add_subplot(132, projection='3d')
    surf_uy = ax2.plot_surface(Y, Z, uy_fluct, cmap='viridis', vmin=vmin, vmax=vmax)
    ax2.set_title("Fluctuations $u'_y$")
    ax2.set_xlabel("Height (y)")
    ax2.set_ylabel("Span (z)")
    ax2.set_zlabel("Fluctuations $u'_y$")

    # Fluctuations for u_z
    ax3 = fig.add_subplot(133, projection='3d')
    surf_uz = ax3.plot_surface(Y, Z, uz_fluct, cmap='viridis', vmin=vmin, vmax=vmax)
    ax3.set_title("Fluctuations $u'_z$")
    ax3.set_xlabel("Height (y)")
    ax3.set_ylabel("Span (z)")
    ax3.set_zlabel("Fluctuations $u'_z$")

    # Adjust layout manually instead of using tight_layout
    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1, wspace=0.3)

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CFD post-processing script")
    parser.add_argument("--filename", '-f', type=str, required=True, help="Input binary file containing the CFD data")
    parser.add_argument("--x_index", '-x', type=int, default=0, help="Index in the x direction to calculate profiles")
    parser.add_argument("--viscosity", '-v', type=float, default=0.89, help="Dynamic viscosity")
    parser.add_argument("--shearstress", '-s', action='store_true', help="Include shear stress in the plots")
    parser.add_argument("--fluctuations", '-fl', action='store_true', help="Include fluctuations in the plots")
    parser.add_argument("--min_max_threshold", type=float, default=0.99, help="Threshold for min-max method (default=0.1)")
    parser.add_argument("--gradient_threshold_factor", type=float, default=0.05, help="Threshold factor for gradient method (default=0.5)")
    args = parser.parse_args()

    # Read the data
    time, nx, ny, nz, x, y, z, ux, uy, uz, pp = read_fields(args.filename)

    # Compute average velocity profiles
    avg_ux, avg_uy, avg_uz = compute_average_velocity(ux, uy, uz, args.x_index)

    # Compute mixing layer limits
    min_max_limits, gradient_limits = compute_mix_layer_limits(avg_ux, y, args.min_max_threshold, args.gradient_threshold_factor)

    # Plot velocity profiles
    plot_velocity_profiles(y, avg_ux, avg_uy, avg_uz, args.viscosity, args.shearstress, 
            min_max_limits, gradient_limits, x, args.x_index)

    # Plot fluctuations if required
    if args.fluctuations:
        ux_fluct, uy_fluct, uz_fluct = compute_fluctuations(ux, uy, uz, avg_ux, avg_uy, avg_uz, args.x_index)
        plot_fluctuations(ux_fluct, uy_fluct, uz_fluct, x, y, z, args.x_index)
    
    # End of execution
    sys.exit(0)
