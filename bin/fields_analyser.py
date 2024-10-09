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

def plot_kolmogorov_spectrum(ux, uy, uz, x, y, z, nx, ny, nz, 
        epsilon, time, filename, png):
    """
    Plots the Kolmogorov energy spectrum for a 3D turbulent flow.
    
    Parameters:
    - ux, uy, uz: 3D arrays representing the velocity components of the flow.
    - x, y, z: 1D arrays representing the spatial grid points in each direction.
    - nx, ny, nz: Integers representing the number of grid points in x, y, and z directions.
    - epsilon: Scalar, the turbulent dissipation rate.
    - time: Scalar, the time at which the spectrum is plotted.
    - filename: String, the name of the file associated with the data.
    - png: Boolean, whether to save the figure as a PNG or display it.
    
    This function calculates the 3D Fourier transform of the velocity fields, computes the energy spectrum, 
    reduces it to 1D based on the magnitude of the wave vector |k|, and then compares it to the theoretical 
    Kolmogorov spectrum (k^(-5/3)).
    """

    # Compute the 3D Fourier transform of the velocity fields and normalize by grid size
    uxf = np.fft.fftn(ux) / (nx * ny * nz)
    uyf = np.fft.fftn(uy) / (nx * ny * nz)
    uzf = np.fft.fftn(uz) / (nx * ny * nz)
    
    # Spectral energy: E(k) = 0.5 * (|ûx(k)|² + |ûy(k)|² + |ûz(k)|²)
    energy_spectrum = 0.5 * (np.abs(uxf)**2 + np.abs(uyf)**2 + np.abs(uzf)**2)
    
    # Wave numbers in 3D
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    kx = np.fft.fftfreq(nx, d=dx) * 2 * np.pi
    ky = np.fft.fftfreq(ny, d=dy) * 2 * np.pi
    kz = np.fft.fftfreq(nz, d=dz) * 2 * np.pi
    
    # 3D wave number grids
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    
    # Magnitude of the wave vector |k|
    k = np.sqrt(KX**2 + KY**2 + KZ**2)
    
    # Convert 3D energy spectrum to 1D by averaging in spherical shells
    k_bins = np.linspace(0, np.max(k), nx//2)
    energy_bins = np.zeros(len(k_bins) - 1)
    
    for i in range(len(k_bins) - 1):
        mask = (k >= k_bins[i]) & (k < k_bins[i+1])
        energy_bins[i] = np.sum(energy_spectrum[mask])
    
    # Center the k values for display
    k_values = 0.5 * (k_bins[:-1] + k_bins[1:])
    
    # Kolmogorov constant (typically around 1.5)
    kolmogorov_constant = 1.5
    
    # Plotting the energy spectrum
    plt.figure(figsize=(10, 6))
    # Plot the Kolmogorov -5/3 slope
    plt.loglog(k_values, kolmogorov_constant * k_values**(-5/3) * epsilon ** (2/3), 
               linestyle='--', color='r', label='$k^{-5/3}$ slope')
    
    # Plot the calculated energy spectrum
    plt.loglog(k_values, energy_bins, label='Energy spectrum')
    
    # Set limits for the y-axis: [1e-7, 1]
    plt.ylim(1e-7, 1)
    # Set specific ticks for the y-axis
    plt.yticks([1e-6, 1e-4, 1e-2, 1e0])
    
    plt.xlabel('Wave number $k$')
    plt.ylabel('Energy $E(k)$')
    plt.title(f'Kolmogorov Energy Spectrum at time = {time:.0f}')
    plt.grid(True, which="both", ls="--")
    plt.legend()

    output_dir = './images/'
    output_filename = os.path.join(output_dir, f'kolmogorov_spectrum_{os.path.basename(filename).replace(".bin", "")}.png')
    
    if png:
        plt.savefig(output_filename)
        print(f"Kolmogorov spectrum saved as {output_filename}")
    else:
        plt.show()
    
    # Close the plot to free memory
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
            plot_kolmogorov_spectrum(ux, uy, uz, x, y, z, nx, ny, nz, 
                    epsilon, time, filename, args.png)

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
