#!/bin/python3
# encoding: utf-8

import sys
import glob
import argparse
import itertools

import matplotlib.pyplot as plt

import numpy as np
from scipy.interpolate import RegularGridInterpolator

def print_block(title, entries):
    """Print a well-formatted block of labeled entries."""
    print(f"\n--- {title} ---")
    for label, val in entries:
        print(f"  {label:<28}: {val}")
    print("-" * 48)


def read_fields(filename):
    with open(filename, 'rb') as f:
        time = np.fromfile(f, dtype=np.float64, count=1)[0]
        nx, ny, nz = np.fromfile(f, dtype=np.int32, count=3)
        x = np.fromfile(f, dtype=np.float64, count=nx)
        y = np.fromfile(f, dtype=np.float64, count=ny)
        z = np.fromfile(f, dtype=np.float64, count=nz)
        n = nx * ny * nz
        ux = np.fromfile(f, dtype=np.float64, count=n).reshape((nx, ny, nz), order='F')
        uy = np.fromfile(f, dtype=np.float64, count=n).reshape((nx, ny, nz), order='F')
        uz = np.fromfile(f, dtype=np.float64, count=n).reshape((nx, ny, nz), order='F')
        pp = np.fromfile(f, dtype=np.float64, count=n).reshape((nx, ny, nz), order='F')
    return time, x, y, z, ux, uy, uz, pp

def compute_momentum_thickness(y, u_profile, U1, U2):
    """Compute momentum thickness delta_theta for a given velocity profile."""
    delta_u = U1 - U2
    integrand = ((U1 - u_profile) / delta_u) * ((u_profile - U2) / delta_u)
    return np.trapz(integrand, y)

def plot_velocity_profiles(time, x, y, z, ux, uy, uz, x_plot, y_plot, z_plot, direction='x'):
    """
    Plot 1D velocity profiles along a given direction, with proper axis orientation.
    If direction is 'y', the profile is plotted horizontally (u vs y).
    """

    # === 1D profile extraction ===
    interp_ux = RegularGridInterpolator((x, y, z), ux, bounds_error=False, fill_value=None)
    interp_uy = RegularGridInterpolator((x, y, z), uy, bounds_error=False, fill_value=None)
    interp_uz = RegularGridInterpolator((x, y, z), uz, bounds_error=False, fill_value=None)

    if direction == 'x':
        coord_range = x
        pts = np.array([[xi, y_plot, z_plot] for xi in x])
    elif direction == 'y':
        coord_range = y
        pts = np.array([[x_plot, yi, z_plot] for yi in y])
    elif direction == 'z':
        coord_range = z
        pts = np.array([[x_plot, y_plot, zi] for zi in z])
    else:
        raise ValueError("Direction must be 'x', 'y', or 'z'")

    ux_vals = interp_ux(pts)
    uy_vals = interp_uy(pts)
    uz_vals = interp_uz(pts)

    plt.figure(figsize=(16, 9))
    # Daltonian-friendly colors
    if direction == 'y':
        # Horizontal plot: velocity vs y
        plt.plot(ux_vals, coord_range, label='u_x', color='#0072B2')
        plt.plot(uy_vals, coord_range, label='u_y', color='#E69F00')
        plt.plot(uz_vals, coord_range, label='u_z', color='#009E73')
        plt.xlabel('Velocity')
        plt.ylabel('y')
        plt.title(f'Velocity profiles along y at (x={x_plot}, z={z_plot})\n$t = {time:.1f}$')
    else:
        # Default: direction on x-axis
        plt.plot(coord_range, ux_vals, label='u_x', color='#0072B2')
        plt.plot(coord_range, uy_vals, label='u_y', color='#E69F00')
        plt.plot(coord_range, uz_vals, label='u_z', color='#009E73')
        plt.xlabel(f'{direction}-coordinate')
        plt.ylabel('Velocity')
        plt.title(f'Velocity profiles at (x={x_plot}, y={y_plot}, z={z_plot})\n$t = {time:.1f}$')

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_mixing_layer_profiles(x, y, z, ux, uy, uz, x_plot, z_plot, time, re_init, threshold_ratio=0.05):
    """
    Analyze and plot the mixing layer at a given (x,z) position:
    - Mean velocity profiles (averaged over x and z)
    - Local velocity line at (x_plot, z_plot)
    - Fluctuations relative to the mean profiles
    "
    """
    # Compute mean profiles along y (averaged over x and z)
    ux_profile = np.mean(ux, axis=(0, 2))
    uy_profile = np.mean(uy, axis=(0, 2))
    uz_profile = np.mean(uz, axis=(0, 2))

    # Extract profile at (x_plot, z_plot)
    ix = np.argmin(np.abs(x - x_plot))
    iz = np.argmin(np.abs(z - z_plot))
    ux_line = ux[ix, :, iz]
    uy_line = uy[ix, :, iz]
    uz_line = uz[ix, :, iz]

    # Compute fluctuations
    ux_fluct = ux_line - ux_profile
    uy_fluct = uy_line - uy_profile
    uz_fluct = uz_line - uz_profile

    # Estimate mixing layer bounds using ux_profile
    dux_dy = np.gradient(ux_profile, y, edge_order=2)
    grad_max_idx = np.argmax(np.abs(dux_dy))
    grad_max_val = np.abs(dux_dy[grad_max_idx])
    y_center = y[grad_max_idx]

    ux_min, ux_max = ux_profile[0], ux_profile[-1]
    ux_01 = ux_min + 0.01 * (ux_max - ux_min)
    ux_99 = ux_min + 0.99 * (ux_max - ux_min)

    y_01 = np.interp(ux_01, ux_profile, y)
    y_99 = np.interp(ux_99, ux_profile, y)

    lower_idx = grad_max_idx
    while lower_idx > 0 and np.abs(dux_dy[lower_idx]) > threshold_ratio * grad_max_val:
        lower_idx -= 1
    upper_idx = grad_max_idx
    while upper_idx < len(y) - 1 and np.abs(dux_dy[upper_idx]) > threshold_ratio * grad_max_val:
        upper_idx += 1

    y_lower = y[lower_idx]
    y_upper = y[upper_idx]
    delta = abs(y_upper - y_lower)

    # Estimate effective Reynolds number
    du = abs(ux_max - ux_min)
    delta_init = 1.0
    nu = du * delta_init / re_init
    re_new = du * delta / nu

    # Plotting
    fig, axs = plt.subplots(1, 3, figsize=(16, 9), sharey=True)
    y_min, y_max = y[0], y[-1]

    def draw_mixing_lines(ax, show_legend=False, delta=None, delta_theta=None):
        ax.axhline(y_01, color='y', linestyle='--', linewidth=1, label='99% bounds' if show_legend else None)
        ax.axhline(y_99, color='y', linestyle='--', linewidth=1)
        ax.axhline(y_lower, color='k', linestyle='--', linewidth=1, label='Mixing layer bounds' if show_legend else None)
        ax.axhline(y_upper, color='k', linestyle='--', linewidth=1)
        ax.axhline(y_center, color='k', linestyle=':', linewidth=1, label='Center' if show_legend else None)
        x_limits = ax.get_xlim()
        ax.fill_betweenx(y=[y_lower, y_upper], x1=x_limits[0], x2=x_limits[1],
                         color='gray', alpha=0.1, label='Mixing layer' if show_legend else None)
                
        y_low_t = y_center - 0.5 * delta_theta
        y_up_t = y_center + 0.5 * delta_theta
        ax.fill_betweenx(
            y=[y_low_t, y_up_t],
            x1=x_limits[0],
            x2=x_limits[1],
            color='blue',
            alpha=0.08,
            label='Momentum thickness' if show_legend else None
        )

        if show_legend:
            # Position under legend (bottom-left of axis)
            ax.annotate(
                rf'$\delta ~= {delta:.3f}$',
                xy=(0.75, 0.35),
                xycoords='axes fraction',
                fontsize=12,
                backgroundcolor='white',
                bbox=dict(boxstyle='round,pad=0.3', edgecolor='gray', facecolor='white')
            )
            ax.annotate(
                rf'$\delta_\theta = {delta_theta:.3f}$',
                xy=(0.75, 0.30),
                xycoords='axes fraction',
                fontsize=12,
                backgroundcolor='white',
                bbox=dict(boxstyle='round,pad=0.3', edgecolor='gray', facecolor='white')
            )

    # Compute momentum thickness δ_theta
    u1 = ux_profile[-1]
    u2 = ux_profile[ 0]
    delta_theta = compute_momentum_thickness(y, ux_profile, u1, u2)
    # Compute Reynolds based on momentum thickness
    re_theta = du * delta_theta / nu

    print_block("Mixing layer Reynolds analysis", [
        ("du", f"{du:.4f}"),
        ("nu (from Re_init)", f"{nu:.2e}"),
        ("delta (gradient-based)", f"{delta:.4f}"),
        ("delta_theta (momentum)", f"{delta_theta:.4f}"),
        ("Re_delta", f"{re_new:.1f}"),
        ("Re_delta_theta", f"{re_theta:.1f}")
    ])

    # Plot <u_x>
    axs[0].plot(ux_profile, y, label=r'$\langle u_x \rangle$', color='tab:blue')
    draw_mixing_lines(axs[0], show_legend=True, delta=delta, delta_theta=delta_theta)
    axs[0].set_title(r'Profile $\langle u_x(y) \rangle$', fontsize=14)
    axs[0].set_xlabel(r'$u_x$', fontsize=12)
    axs[0].set_ylabel(r'$y$', fontsize=12)
    axs[0].set_ylim(y_min, y_max)
    axs[0].grid(True, linestyle='--', alpha=0.7)
    axs[0].legend()

    # Plot <u_y> and <u_z>
    axs[1].plot(uy_profile, y, label=r'$\langle u_y \rangle$', color='tab:green')
    axs[1].plot(uz_profile, y, label=r'$\langle u_z \rangle$', color='tab:orange')
    draw_mixing_lines(axs[1], show_legend=False, delta=delta, delta_theta=delta_theta)
    axs[1].set_title(r'Profiles $\langle u_y(y) \rangle$, $\langle u_z(y) \rangle$', fontsize=14)
    axs[1].set_xlabel(r'$u_y$, $u_z$', fontsize=12)
    axs[1].set_ylim(y_min, y_max)
    axs[1].grid(True, linestyle='--', alpha=0.7)
    axs[1].legend()

    # Plot fluctuations u'_i = u_i(x_plot, z_plot) - <u_i>
    axs[2].plot(ux_fluct, y, label=r"$u'_x$", color='tab:red')
    axs[2].plot(uy_fluct, y, label=r"$u'_y$", color='tab:purple')
    axs[2].plot(uz_fluct, y, label=r"$u'_z$", color='tab:brown')
    draw_mixing_lines(axs[2], show_legend=False, delta=delta, delta_theta=delta_theta)
    axs[2].set_title(r"Fluctuations $u'_i(y)$ at $x,z$", fontsize=14)
    axs[2].set_xlabel(r"$u'_i$", fontsize=12)
    axs[2].set_ylim(y_min, y_max)
    axs[2].grid(True, linestyle='--', alpha=0.7)
    axs[2].legend()

    fig.suptitle(
        rf"Mixing layer analysis at $x = {x_plot}$, $z = {z_plot}$, $t = {time:.1f}$",
        fontsize=16
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

def plot_velocity_fluctuations_plan(x, y, z, ux, uy, uz, x_plot, y_plot, z_plot, time, direction='x'):
    """
    Plot velocity fluctuations in a plane orthogonal to the specified direction.
    If direction = 'x', plot in (y,z)
    If direction = 'y', plot in (x,z)
    If direction = 'z', plot in (x,y)
    """

    if direction == 'x':
        ix = np.argmin(np.abs(x - x_plot))
        ux_plane = ux[ix, :, :]
        uy_plane = uy[ix, :, :]
        uz_plane = uz[ix, :, :]
        ux_mean = np.mean(ux, axis=0)
        uy_mean = np.mean(uy, axis=0)
        uz_mean = np.mean(uz, axis=0)
        coord1, coord2 = y, z
        label1, label2 = 'y', 'z'
    elif direction == 'y':
        iy = np.argmin(np.abs(y - y_plot))
        ux_plane = ux[:, iy, :]
        uy_plane = uy[:, iy, :]
        uz_plane = uz[:, iy, :]
        ux_mean = np.mean(ux, axis=1)
        uy_mean = np.mean(uy, axis=1)
        uz_mean = np.mean(uz, axis=1)
        coord1, coord2 = x, z
        label1, label2 = 'x', 'z'
    elif direction == 'z':
        iz = np.argmin(np.abs(z - z_plot))
        ux_plane = ux[:, :, iz]
        uy_plane = uy[:, :, iz]
        uz_plane = uz[:, :, iz]
        ux_mean = np.mean(ux, axis=2)
        uy_mean = np.mean(uy, axis=2)
        uz_mean = np.mean(uz, axis=2)
        coord1, coord2 = x, y
        label1, label2 = 'x', 'y'
    else:
        raise ValueError("Direction must be 'x', 'y', or 'z'")

    # Compute fluctuations
    fluct_x = ux_plane - ux_mean
    fluct_y = uy_plane - uy_mean
    fluct_z = uz_plane - uz_mean

    # Color scale
    vabs = max(np.max(np.abs(fluct_x)), np.max(np.abs(fluct_y)), np.max(np.abs(fluct_z)))
    vmin, vmax = -vabs, vabs

    fig, axs = plt.subplots(1, 3, figsize=(16, 9))
    cmap = 'seismic'

    extent = [coord2.min(), coord2.max(), coord1.min(), coord1.max()]
    im0 = axs[0].imshow(fluct_x, extent=extent, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
    axs[0].set_title(rf"Fluctuation $u'_x({label1},{label2})$")
    axs[0].set_xlabel(f"{label2}")
    axs[0].set_ylabel(f"{label1}")
    fig.colorbar(im0, ax=axs[0])

    im1 = axs[1].imshow(fluct_y, extent=extent, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
    axs[1].set_title(rf"Fluctuation $u'_y({label1},{label2})$")
    axs[1].set_xlabel(f"{label2}")
    axs[1].set_ylabel(f"{label1}")
    fig.colorbar(im1, ax=axs[1])

    im2 = axs[2].imshow(fluct_z, extent=extent, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
    axs[2].set_title(rf"Fluctuation $u'_z({label1},{label2})$")
    axs[2].set_xlabel(f"{label2}")
    axs[2].set_ylabel(f"{label1}")
    fig.colorbar(im2, ax=axs[2])

    fig.suptitle(
        f"Velocity fluctuations in the {label1}-{label2} plane at time t = {time:.2f}",
        fontsize=16
    )
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.show()

def plot_energy_spectrum(ux, uy, uz, x, y, z, time, LES=False):
    """
    Compute and plot the isotropic energy spectrum E(k) from a 3D FFT
    of the velocity field. The energy is binned by the magnitude of 
    the wavevector |k| and plotted on a log-log scale.

    Parameters:
    - ux, uy, uz: 3D velocity field components (numpy arrays)
    - x, y, z: grid coordinates (1D arrays)
    - time: simulation time (float)
    """

    nx, ny, nz = len(x), len(y), len(z)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]

    # Perform 3D FFT of each velocity component
    uxf = np.fft.fftn(ux)
    uyf = np.fft.fftn(uy)
    uzf = np.fft.fftn(uz)

    # Compute spectral energy density
    energy_spectral = 0.5 * (np.abs(uxf)**2 + np.abs(uyf)**2 + np.abs(uzf)**2)

    # Create frequency grids
    kx = np.fft.fftfreq(nx, dx) * 2 * np.pi
    ky = np.fft.fftfreq(ny, dy) * 2 * np.pi
    kz = np.fft.fftfreq(nz, dz) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')

    # Compute |k| magnitude
    k_mag = np.sqrt(KX**2 + KY**2 + KZ**2).flatten()
    energy_spectral = energy_spectral.flatten()

    # Bin energy values by |k|
    k_bins = np.linspace(0, k_mag.max(), 100)
    k_centers = 0.5 * (k_bins[:-1] + k_bins[1:])
    E_k = np.zeros_like(k_centers)

    for i in range(len(k_bins) - 1):
        mask = (k_mag >= k_bins[i]) & (k_mag < k_bins[i+1])
        E_k[i] = np.sum(energy_spectral[mask])

    # Plotting
    plt.figure(figsize=(16, 9))
    plt.loglog(k_centers, E_k, '-.', markersize=4, label=r"$E(k)$")
    plt.title("Isotropic energy spectrum", fontsize=14)
    plt.xlabel(r"Wavenumber $k = |\mathbf{k}|$", fontsize=12)
    plt.ylabel(r"$E(k)$", fontsize=12)
    plt.grid(True, which='both', linestyle=':', alpha=0.7)

    # Add reference k^-5/3 slope
    ref_k = np.array([k_centers[10], k_centers[40]])
    ref_E = ref_k**(-5/3) * E_k[10] / (k_centers[10]**(-5/3))
    plt.loglog(ref_k, ref_E, 'k--', label=r"$k^{-5/3}$")

    if LES:
        delta = (dx * dy * dz)**(1/3)
        k_cutoff = np.pi / delta
        lambda_cutoff = 2 * delta  # corresponding physical scale
        plt.axvline(k_cutoff, color='red', linestyle='--', linewidth=1.5, label=r"$k_c = \pi / \Delta$")
        plt.text(k_cutoff * 1.1, max(E_k) * 1e-2,
                 f"LES cutoff: λ ≈ {lambda_cutoff:.3f}", color='red', fontsize=10)
    
    plt.legend(loc='lower left', frameon=True)
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.tight_layout()
    plt.show()

def analyze_turbulent_development(x, y, z, ux, uy, uz, time, LES=False):
    """
    Analyze turbulence development using mean TKE and enstrophy.
    Compare to empirical thresholds to automatically assess turbulent state.
    """

    # Remove mean in x-direction (streamwise)
    ux_fluct = ux - np.mean(ux, axis=0)
    uy_fluct = uy - np.mean(uy, axis=0)
    uz_fluct = uz - np.mean(uz, axis=0)

    # Turbulent kinetic energy
    tke = 0.5 * (ux_fluct**2 + uy_fluct**2 + uz_fluct**2)
    mean_tke = np.mean(tke)

    # Vorticity components (central gradients)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]

    dudy = np.gradient(ux, axis=1, edge_order=2) / dy
    dudz = np.gradient(ux, axis=2, edge_order=2) / dz
    dvdx = np.gradient(uy, axis=0, edge_order=2) / dx
    dvdz = np.gradient(uy, axis=2, edge_order=2) / dz
    dwdx = np.gradient(uz, axis=0, edge_order=2) / dx
    dwdy = np.gradient(uz, axis=1, edge_order=2) / dy

    wx = dwdy - dvdz
    wy = dudz - dwdx
    wz = dvdx - dudy

    enstrophy = 0.5 * (wx**2 + wy**2 + wz**2)
    mean_enstrophy = np.mean(enstrophy)

    # Thresholds (adjustable)
    TKE_THRESHOLD = 1e-2
    ENSTROPHY_THRESHOLD = 1.0

    block = [
        ("Time", f"{time:.3f}"),
        ("Mean TKE", f"{mean_tke:.4e}"),
        ("Mean enstrophy", f"{mean_enstrophy:.4e}")
    ]

    if LES:
        nut, delta = compute_smagorinsky_viscosity(x, y, z, ux, uy, uz)
        nut_min = np.min(nut)
        nut_max = np.max(nut)
        nut_mean = np.mean(nut)
        lambda_cutoff = 2 * delta

        block += [
            ("Smagorinsky nu_t min", f"{nut_min:.3e}"),
            ("Smagorinsky nu_t max", f"{nut_max:.3e}"),
            ("Smagorinsky nu_t mean", f"{nut_mean:.3e}"),
            ("LES cutoff length scale", f"lambda = {lambda_cutoff:.4f}")
        ]

    print_block("Turbulence development analysis", block)

    # Regime assessment remains outside for clarity
    print(f"\n  Regime assessment         : ", end="")
    if mean_tke > TKE_THRESHOLD and mean_enstrophy > ENSTROPHY_THRESHOLD:
        print("TURBULENT")
    elif mean_tke > 0.1 * TKE_THRESHOLD:
        print("TRANSITIONAL")
    else:
        print("LAMINAR or UNDERDEVELOPED")
    print("-" * 48)


def compute_smagorinsky_viscosity(x, y, z, ux, uy, uz, Cs=0.17):
    """
    Compute the Smagorinsky eddy viscosity nu_t using local strain rates.
    
    Parameters:
        x, y, z : 1D grid coordinates
        ux, uy, uz : 3D velocity fields
        Cs : Smagorinsky constant (default: 0.17)
    
    Returns:
        nut : 3D array of turbulent viscosity (same shape as ux)
        delta : filter scale (constant)
    """
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    delta = (dx * dy * dz) ** (1 / 3)

    # Compute velocity gradients
    dux_dx = np.gradient(ux, axis=0, edge_order=2) / dx
    dux_dy = np.gradient(ux, axis=1, edge_order=2) / dy
    dux_dz = np.gradient(ux, axis=2, edge_order=2) / dz
    duy_dx = np.gradient(uy, axis=0, edge_order=2) / dx
    duy_dy = np.gradient(uy, axis=1, edge_order=2) / dy
    duy_dz = np.gradient(uy, axis=2, edge_order=2) / dz
    duz_dx = np.gradient(uz, axis=0, edge_order=2) / dx
    duz_dy = np.gradient(uz, axis=1, edge_order=2) / dy
    duz_dz = np.gradient(uz, axis=2, edge_order=2) / dz

    # Compute strain-rate magnitude |S|
    S2 = (
        2 * dux_dx**2 +
        2 * duy_dy**2 +
        2 * duz_dz**2 +
        (dux_dy + duy_dx)**2 +
        (dux_dz + duz_dx)**2 +
        (duy_dz + duz_dy)**2
    )
    nut = (Cs * delta) ** 2 * np.sqrt(0.5 * S2)
    return nut, delta

def main():

    parser = argparse.ArgumentParser(
        description="Plot velocity profiles and turbulence statistics from CFD binary output."
    )
    parser.add_argument(
        '--filename', '-fi',
        nargs='*',
        default=glob.glob("fields_??????.bin"),
        help="Fichiers binaires à analyser (wildcards autorisés). Par défaut : fields_??????.bin"
    )
    parser.add_argument('--x', type=float, default=0.0, help="x coordinate for profile extraction")
    parser.add_argument('--y', type=float, default=0.0, help="y coordinate for profile extraction")
    parser.add_argument('--z', type=float, default=0.0, help="z coordinate for profile extraction")

    parser.add_argument('--summary', '-s', action='store_true',
    help="Print a summary of the field (filename, time, min/max of velocity and pressure)")
    
    parser.add_argument('--plot_velocity_profiles', '-p',  action='store_true',
            help="Display velocity profiles in i-direction (default x). You can change the direction with --direction option")

    parser.add_argument('--direction', choices=['x', 'y', 'z'], default='x',
            help="Direction along which to extract the velocity profile (default: z)")

    parser.add_argument('--mixing_layer', '-ml', type=float, metavar='RE_INIT',
            help="Activate mixing layer mode (3-panel plot in the y-direction)")

    parser.add_argument('--fluctuations', '-f', action='store_true',
            help="Display velocity fluctuations u'_i in the (y,z) plane at given x")

    parser.add_argument('--turbulence', '-t', action='store_true',
            help="Display turbulence statistics and plot the energy spectrum")
    parser.add_argument('--LES', action='store_true',
            help="Enable LES mode to show expected spectral cutoff due to Smagorinsky model")

    args = parser.parse_args()
 
    file_list = list(itertools.chain.from_iterable(glob.glob(f) for f in args.filename))
    file_list = sorted(set(file_list))  # dédoublonner + trier

    if not file_list:
        print("Aucun fichier trouvé pour les motifs donnés.")
        return 1

    for fname in file_list:
        print(f"\n========= Analyse de {fname} =========")
        time, x, y, z, ux, uy, uz, pp = read_fields(fname)
        if args.summary:
            block = [
                ("File", fname),
                ("Time", f"{time:.1f}"),
                ("nx, ny, nz", f"{len(x)}, {len(y)}, {len(z)}"),
                ("ux min/max", f"{ux.min():.3e} / {ux.max():.3e} "),
                ("uy min/max", f"{uy.min():.3e} / {uy.max():.3e} "),
                ("uz min/max", f"{uz.min():.3e} / {uz.max():.3e} "),
                ("Pressure min", f"{pp.min():.3e}"),
                ("Pressure max", f"{pp.max():.3e}")
            ]
            print_block("Field summary", block)

    
        if args.plot_velocity_profiles:
            plot_velocity_profiles(time, x, y, z, ux, uy, uz, args.x, args.y, args.z, args.direction)
    
        if args.mixing_layer:
            plot_mixing_layer_profiles(x, y, z, ux, uy, uz, args.x, args.z, time, re_init=args.mixing_layer)
    
        if args.fluctuations:
            plot_velocity_fluctuations_plan(
                    x, y, z, ux, uy, uz,
                    args.x, args.y, args.z,
                    time, direction=args.direction
                    )

        if args.turbulence:
            plot_energy_spectrum(ux, uy, uz, x, y, z, time, LES=args.LES)
            analyze_turbulent_development(x, y, z, ux, uy, uz, time, args.LES)

    return 0

if __name__ == "__main__":

    err = main()
    sys.exit(0)
