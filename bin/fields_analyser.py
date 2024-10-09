import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import glob
import sys
import os

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

def compute_dissipation_rate(ux, uy, uz, dx, dy, dz):
    """
    Estimer le taux de dissipation d'énergie à partir des gradients des vitesses.
    
    ε = 2ν * (Sij * Sij), où Sij est le tenseur des taux de déformation
    """
    # Calcul des dérivées des composantes de vitesse dans chaque direction
    dux_dx = np.gradient(ux, dx, axis=0, edge_order = 2)
    duy_dy = np.gradient(uy, dy, axis=1, edge_order = 2)
    duz_dz = np.gradient(uz, dz, axis=2, edge_order = 2)
    
    dux_dy = np.gradient(ux, dy, axis=1, edge_order = 2)
    dux_dz = np.gradient(ux, dz, axis=2, edge_order = 2)
    duy_dx = np.gradient(uy, dx, axis=0, edge_order = 2)
    duy_dz = np.gradient(uy, dz, axis=2, edge_order = 2)
    duz_dx = np.gradient(uz, dx, axis=0, edge_order = 2)
    duz_dy = np.gradient(uz, dy, axis=1, edge_order = 2)
    
    # Tenseur des taux de déformation Sij
    Sij = 0.5 * (dux_dx**2 + duy_dy**2 + duz_dz**2 +
                 dux_dy**2 + dux_dz**2 +
                 duy_dx**2 + duy_dz**2 +
                 duz_dx**2 + duz_dy**2)

    # Taux de dissipation (en moyenne sur le volume)
    epsilon = np.mean(Sij)
    
    return epsilon

def compute_kolmogorov_scales(epsilon, viscosity):
    """
    Calculer l'échelle de Kolmogorov en longueur, temps et vitesse.
    """
    eta = (viscosity**3 / epsilon)**0.25
    tau_eta = (viscosity / epsilon)**0.5
    v_eta = (viscosity * epsilon)**0.25
    return eta, tau_eta, v_eta

def plot_velocity_profiles(y, avg_ux, avg_uy, avg_uz, viscosity, 
        min_max_limits, gradient_limits, x, x_index, time,
        filename, png):
    """Plot the average velocity profiles and their derivatives."""
    # Calculate the derivatives using np.gradient
    d_avg_ux = np.gradient(avg_ux, y, edge_order = 2)
    d_avg_uy = np.gradient(avg_uy, y, edge_order = 2)
    d_avg_uz = np.gradient(avg_uz, y, edge_order = 2)

    # Define colors suitable for color-blind users
    color_ux = '#E69F00'  # orange
    color_uy = '#56B4E9'  # sky blue
    color_uz = '#009E73'  # teal
    color_thickness = '#CC79A7'  # pink for mixing layer limits

    plt.figure(figsize=(12, 9))
    # Set the overall title for the figure, including x[x_index]
    plt.suptitle(f'Velocity Profiles for x = {x[x_index]:.2f} at time = {time:.0f}', fontsize=14)

    # Plotting the average velocity profile for u_x
    plt.subplot(1, 3, 1)
    plt.plot(avg_ux, y, label='$\\langle u_x \\rangle$', color=color_ux)

    if (min_max_limits):
        # Plotting the limits of the mixing layer
        lower_limit_min_max, upper_limit_min_max = min_max_limits
        plt.axhline(y=lower_limit_min_max, color=color_thickness, linestyle=':', label='Mixing Layer Limit (Min-Max)')
        plt.axhline(y=upper_limit_min_max, color=color_thickness, linestyle=':')

        lower_limit_gradient, upper_limit_gradient = gradient_limits
        if lower_limit_gradient is not None and upper_limit_gradient is not None:
            plt.axhline(y=lower_limit_gradient, color=color_thickness, linestyle='-.', label='Mixing Layer Limit (Gradient)')
            plt.axhline(y=upper_limit_gradient, color=color_thickness, linestyle='-.')

    plt.title('Velocity Profile $\\langle u_x \\rangle$')
    plt.xlabel('Velocity ($u_x$)')
    plt.ylabel('Height (y)')
    plt.ylim([min(y), max(y)])  # Set y-axis limits
    plt.grid()
    plt.legend()

    # Plotting the average velocity profiles for u_y and u_z
    plt.subplot(1, 3, 2)
    plt.plot(avg_uy, y, label='$\\langle u_y \\rangle$', color=color_uy, linestyle='--', linewidth=1.5)
    plt.plot(avg_uz, y, label='$\\langle u_z \\rangle$', color=color_uz, linestyle='-.', linewidth=1.5)
    if (min_max_limits):
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
    if (min_max_limits):
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

    output_filename = f'./images/velocity_profile_{os.path.basename(filename).replace(".bin", "")}.png'
    plt.tight_layout()
    if (png):
        plt.savefig(output_filename)
        print(f"Velocity profile saved as {output_filename}")
    else:
        plt.show()
    plt.close()

def plot_fluctuations(ux_fluct, uy_fluct, uz_fluct, x, y, z, x_index, time, filename, png):
    """Plot the velocity fluctuations in 3D with a single colormap and shared colorbar."""
    fig = plt.figure(figsize=(18, 6))
    # Set the overall title for the figure, including x[x_index]
    fig.suptitle(f'Velocity Fluctuations for x = {x[x_index]:.2f} at time = {time:.0f}', fontsize=14)

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
    set_axis_ticks(ax1, ux_fluct, axis='z')

    # Fluctuations for u_y
    ax2 = fig.add_subplot(132, projection='3d')
    surf_uy = ax2.plot_surface(Y, Z, uy_fluct, cmap='viridis', vmin=vmin, vmax=vmax)
    ax2.set_title("Fluctuations $u'_y$")
    ax2.set_xlabel("Height (y)")
    ax2.set_ylabel("Span (z)")
    ax2.set_zlabel("Fluctuations $u'_y$")
    set_axis_ticks(ax2, uy_fluct, axis='z')

    # Fluctuations for u_z
    ax3 = fig.add_subplot(133, projection='3d')
    surf_uz = ax3.plot_surface(Y, Z, uz_fluct, cmap='viridis', vmin=vmin, vmax=vmax)
    ax3.set_title("Fluctuations $u'_z$")
    ax3.set_xlabel("Height (y)")
    ax3.set_ylabel("Span (z)")
    ax3.set_zlabel("Fluctuations $u'_z$")
    set_axis_ticks(ax3, uz_fluct, axis='z')

    # Adjust layout manually instead of using tight_layout
    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1, wspace=0.3)

    output_filename = f'./images/fluctuations_{os.path.basename(filename).replace(".bin", "")}.png'
    if png:
        plt.savefig(output_filename)
        print(f"Fluctuations plot saved as {output_filename}")
    else:
        plt.show()
    plt.close()

def plot_kolmogorov_spectrum(ux, uy, uz, x, y, z, epsilon, time, filename, png, num_bins=50, k_min_cutoff=1e-3):
    """
    Plot the energy spectrum of the velocity field, with band-averaged smoothing and correction for edge effects.
    
    Parameters:
        ux, uy, uz : np.ndarray
            Velocity components in 3D space.
        x, y, z : np.ndarray
            Spatial coordinates in the x, y, and z directions.
        epsilon : float
            Dissipation rate 
        time : float
            Time at which the spectrum is computed.
        filename : str
            Filename used for saving the plot.
        png : bool
            If True, saves the plot as a PNG file. Otherwise, shows the plot.
        num_bins : int
            Number of bands for the smoothing process (default is 50).
        k_min_cutoff : float
            Minimum value of k to include in the plot to avoid edge effects (default is 1e-3).
    """
    # Compute the 3D Fourier Transform of velocity components
    ux_fft = np.fft.fftn(ux) / ((x[-1]-x[0])**2 * (y[-1]-y[0])**2 * (z[-1]-z[0])**2)
    uy_fft = np.fft.fftn(uy) / ((x[-1]-x[0])**2 * (y[-1]-y[0])**2 * (z[-1]-z[0])**2) 
    uz_fft = np.fft.fftn(uz) / ((x[-1]-x[0])**2 * (y[-1]-y[0])**2 * (z[-1]-z[0])**2) 
    
    # Calculate the energy spectrum (sum of squares of the FFTs)
    energy_spectrum = 0.5 * (np.abs(ux_fft)**2 + np.abs(uy_fft)**2 + np.abs(uz_fft)**2)

    # Shift zero frequency to the center
    energy_spectrum = np.fft.fftshift(energy_spectrum)

    # Define wave numbers in each direction
    dx = x[1] - x[0]  # Spacing in x-direction
    dy = y[1] - y[0]  # Spacing in y-direction
    dz = z[1] - z[0]  # Spacing in z-direction
    
    kx = np.fft.fftfreq(len(x), d=dx)
    ky = np.fft.fftfreq(len(y), d=dy)
    kz = np.fft.fftfreq(len(z), d=dz)

    kx = np.fft.fftshift(kx)
    ky = np.fft.fftshift(ky)
    kz = np.fft.fftshift(kz)

    # Create a meshgrid of wave numbers
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    k_magnitude = np.sqrt(KX**2 + KY**2 + KZ**2)

    # Flatten the arrays for easy manipulation
    k_magnitude = k_magnitude.flatten()
    energy_spectrum = energy_spectrum.flatten()

    # Sort by wave number
    sort_indices = np.argsort(k_magnitude)
    k_magnitude = k_magnitude[sort_indices]
    energy_spectrum = energy_spectrum[sort_indices]

    # Filter out k = 0 and apply cutoff for small k values to avoid edge effects
    non_zero_indices = (k_magnitude > k_min_cutoff)
    k_magnitude = k_magnitude[non_zero_indices]
    energy_spectrum = energy_spectrum[non_zero_indices]

    # Bin the data to average the energy spectrum in bands of k
    k_min = k_magnitude.min()
    k_max = k_magnitude.max()
    
    # Define logarithmic bins for the wave numbers
    bins = np.logspace(np.log10(k_min), np.log10(k_max), num=num_bins)
    
    # Initialize arrays for binned wave numbers and energy
    binned_k = 0.5 * (bins[1:] + bins[:-1])  # Midpoints of the bins
    binned_energy = np.zeros(len(binned_k))
    
    # Average energy spectrum in each bin
    for i in range(len(bins) - 1):
        in_bin = (k_magnitude >= bins[i]) & (k_magnitude < bins[i + 1])
        if np.sum(in_bin) > 0:
            binned_energy[i] = np.mean(energy_spectrum[in_bin])
        else:
            # If no points are in this bin, assign NaN to handle later
            binned_energy[i] = np.nan

    # Remove NaN values (empty bins)
    valid_indices = ~np.isnan(binned_energy)
    binned_k = binned_k[valid_indices]
    binned_energy = binned_energy[valid_indices]

    # Plotting the energy spectrum
    plt.figure(figsize=(10, 6))
    plt.loglog(binned_k, binned_energy, label='Binned Energy Spectrum', color='b')

    # Kolmogorov -5/3 scaling line for reference
    kolmogorov_constant = 1.5  # Adjust as appropriate
    plt.loglog(binned_k, kolmogorov_constant * binned_k**(-5/3) * epsilon ** (2/3), linestyle='--', color='r', label='$k^{-5/3}$ slope')
    #plt.loglog(binned_k, np.exp(-2. * binned_k), linestyle='-', color='g', label='$e^{-2k}$ slope')

    plt.xlabel('Wave number (k)')
    plt.ylabel('Energy')
    plt.title(f'Kolmogorov Energy Spectrum at time = {time:.0f}')
    plt.grid(True, which="both", ls="--")
    plt.legend()

    output_filename = f'./images/kolmogorov_spectrum_{os.path.basename(filename).replace(".bin", "")}.png'
    if png:
        plt.savefig(output_filename)
        print(f"Kolmogorov spectrum saved as {output_filename}")
    else:
        plt.show()
    plt.close()

def set_axis_ticks(ax, fluct_data, axis='z', n=5, round_value=1):
    """
    Set axis ticks and formatting for a 3D axis.

    Parameters:
    - ax: The axis to modify (3D axes like x, y, or z).
    - fluct_data: The fluctuation data used to determine min and max values for the axis.
    - axis: The axis to which the ticks are applied ('x', 'y', or 'z'). Default is 'z'.
    - n: Number of ticks. Default is 5.
    - round_value: The decimal place for rounding the step. Default is 1 (i.e., round to one decimal).
    """
    min_val, max_val = np.min(fluct_data), np.max(fluct_data)
    step = np.round((max_val - min_val) / n, round_value)
    if step == 0.0:
        step = 0.1  # Default to 0.1 if the step is too small

    if axis == 'z':
        ax.set_zticks(np.arange(min_val, max_val + step, step))
        ax.zaxis.set_major_formatter(FormatStrFormatter(f'%.{round_value}f'))
    elif axis == 'x':
        ax.set_xticks(np.arange(min_val, max_val + step, step))
        ax.xaxis.set_major_formatter(FormatStrFormatter(f'%.{round_value}f'))
    elif axis == 'y':
        ax.set_yticks(np.arange(min_val, max_val + step, step))
        ax.yaxis.set_major_formatter(FormatStrFormatter(f'%.{round_value}f'))


def remove_images(filenames):
    """Remove existing image files."""
    for filename in filenames:
        if os.path.exists(filename):
            os.remove(filename)
            print(f"Removed {filename}")

def main(args):

    # Remove existing images if the option is set
    if args.remove:
        velocity_files = glob.glob('./images/velocity_profile_fields_??????.png')
        fluctuation_files = glob.glob('./images/fluctuations_fields_??????.png')
        kolmogorov_files = glob.glob('./images/kolmogorov_spectrum_fields_??????.png')
        remove_images(velocity_files + fluctuation_files + kolmogorov_files)

    nu = 1. / args.reynolds
    print(f"Dynamic viscosity nu = {nu:.4e}")

    for filename in args.filenames:
        print(f"Processing of the file {filename}")
        # Read the data
        time, nx, ny, nz, x, y, z, ux, uy, uz, pp = read_fields(filename)

        # Compute average velocity profiles
        avg_ux, avg_uy, avg_uz = compute_average_velocity(ux, uy, uz, args.x_index)
        
        if (args.mixing_layer):
            # Compute mixing layer limits
            min_max_limits, gradient_limits = compute_mix_layer_limits(
                avg_ux, y, args.min_max_threshold, args.gradient_threshold_factor)
        else:
            min_max_limits = None
            gradient_limits = None

        # Plot velocity profiles
        plot_velocity_profiles(y, avg_ux, avg_uy, avg_uz, nu, 
                min_max_limits, gradient_limits, x, args.x_index, time,
                filename, args.png)

        # Plot fluctuations if required
        if args.fluctuations:
            ux_fluct, uy_fluct, uz_fluct = compute_fluctuations(ux, uy, 
                uz, avg_ux, avg_uy, avg_uz, args.x_index)
            plot_fluctuations(ux_fluct, uy_fluct, uz_fluct, x, y, z, args.x_index, time, filename, args.png)

        if args.kolmogorov:
            dx = x[1] - x[0]
            dy = y[1] - y[0]
            dz = z[1] - z[0]
            epsilon = compute_dissipation_rate(ux, uy, uz, dx, dy, dz)
            eta, tau_eta, v_eta = compute_kolmogorov_scales(epsilon, nu)
            print(f"Dissipation rate epsilon = {epsilon:.4e}")
            print(f"Kolmogorov scale eta = {eta:.4e}")
            # Plot Kolmogorov energy spectrum
            plot_kolmogorov_spectrum(ux, uy, uz, x, y, z, epsilon, time, filename, args.png)

    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CFD post-processing script", 
            epilog="Example command: python3 fields_analyser.py fields_000000.bin fields_001000.bin fields_002000.bin -re 2500 -x 128 -fl -r")
    default_filenames = glob.glob("fields_??????.bin")
    parser.add_argument("--filenames", '-f', nargs='+', type=str, default=default_filenames, help="List of input binary files containing the CFD data")
    parser.add_argument("--x_index", '-x', type=int, default=0, help="Index in the x direction to calculate profiles")
    parser.add_argument("--reynolds", '-re', type=float, default=1000., help="Reynolds number (nu=1/re)")
    parser.add_argument("--fluctuations", '-fl', action='store_true', help="Include fluctuations in the plots")
    parser.add_argument("--kolmogorov", '-k', action='store_true', help="Include Kolmogorov scale process and plot")
    parser.add_argument("--mixing_layer", '-ml', action='store_true', help="Compute and plot mixing layer imits")
    parser.add_argument("--min_max_threshold", type=float, default=0.99, help="Threshold for min-max method (default=0.1)")
    parser.add_argument("--gradient_threshold_factor", type=float, default=0.05, help="Threshold factor for gradient method (default=0.5)")
    parser.add_argument("--remove", '-r', action='store_true', help='Remove existing images')
    parser.add_argument('-png', action='store_true', help='Generate PNG plots in ./images directory')
    args = parser.parse_args()

    err = main(args)
    
    # End of execution
    sys.exit(err)
