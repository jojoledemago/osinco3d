#!/usr/bin/env python3
# encoding: utf-8

import os
import sys
import glob
import argparse
import itertools
import numpy as np
import matplotlib.pyplot as plt

def print_block(title, entries):
    """
    Print a formatted block of labeled entries.

    Parameters
    ----------
    title : str
        Title of the block.
    entries : list of tuple
        List of (label, value) tuples to display.
    """
    print(f"\n--- {title} ---")
    for label, value in entries:
        print(f"  {label:<30}: {value}")
    print("-" * 48)

def read_fields(filename):
    """
    Read OSINCO3D binary field file.

    Parameters
    ----------
    filename : str

    Returns
    -------
    time : float
    x, y, z : ndarray
    ux, uy, uz, pp : ndarray
    """
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

def compute_rms_profiles(u, axis_mean, axis_profile):
    """
    Compute RMS and fluctuations for a velocity component.

    Parameters
    ----------
    u : ndarray
        Velocity component (3D).
    axis_mean : tuple of int
        Axes over which to average to define mean.
    axis_profile : int
        Axis along which the profile is returned.

    Returns
    -------
    u_rms : 1D array (profile along axis_profile)
    u_fluct : ndarray (same shape as u)
    """
    u_mean = np.mean(u, axis=axis_mean, keepdims=True)
    u_fluct = u - u_mean
    u_rms = np.sqrt(np.mean(u_fluct**2, axis=axis_mean))
    return u_rms, u_fluct

def plot_rms_profiles(coord, ux_rms, uy_rms, uz_rms, tke, intensity, direction, time):
    """
    Plot RMS profiles and TKE.

    Parameters
    ----------
    coord : ndarray
        Coordinates in the direction of profile.
    ux_rms, uy_rms, uz_rms : ndarray
        RMS profiles for velocity components.
    tke : ndarray
        Turbulent kinetic energy profile.
    intensity : ndarray
        Turbulent intensity profile.
    direction : str
        Direction of profile (x/y/z)
    time : float
        Simulation time
    """
    # Determine orthogonal plane name
    planes = {'x': 'y-z', 'y': 'x-z', 'z': 'x-y'}
    plane_name = planes[direction]

    plt.figure(figsize=(16, 9))
    plt.plot(ux_rms, coord, label="u'_x RMS")
    plt.plot(uy_rms, coord, label="u'_y RMS")
    plt.plot(uz_rms, coord, label="u'_z RMS")
    plt.plot(np.sqrt(2 * tke), coord, label="Turbulent intensity √(2k)", linestyle='--')

    plt.xlabel("RMS Value")
    plt.ylabel(direction)
    plt.title(f"RMS Velocity Profiles in {plane_name} plane at time: {time:.2f}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def  analyze_rms_and_tke(ux, uy, uz, coords, direction, time):
    """
    Analyse turbulence quantities and display summary and plots.

    Parameters
    ----------
    ux, uy, uz : ndarray
        Velocity components.
    coords : dict
        Dictionary with keys 'x', 'y', 'z' → coordinate arrays.
    direction : str
        Direction to compute profiles along: 'x', 'y', or 'z'.
    title : str
        Label to use for figure title and output.
    """
    axis_dict = {'x': 0, 'y': 1, 'z': 2}
    profile_axis = axis_dict[direction]
    mean_axes = tuple(ax for ax in (0, 1, 2) if ax != profile_axis)
    coord = coords[direction]

    ux_rms, _ = compute_rms_profiles(ux, mean_axes, profile_axis)
    uy_rms, _ = compute_rms_profiles(uy, mean_axes, profile_axis)
    uz_rms, _ = compute_rms_profiles(uz, mean_axes, profile_axis)

    tke_profile = 0.5 * (ux_rms**2 + uy_rms**2 + uz_rms**2)
    intensity_profile = np.sqrt(2 * tke_profile)

    ux_rms_avg = np.mean(ux_rms)
    uy_rms_avg = np.mean(uy_rms)
    uz_rms_avg = np.mean(uz_rms)

    anisotropy_ratios = [
        ("ux'/uy'", f"{ux_rms_avg / uy_rms_avg if uy_rms_avg else np.nan:.3f}"),
        ("ux'/uz'", f"{ux_rms_avg / uz_rms_avg if uz_rms_avg else np.nan:.3f}"),
        ("uy'/uz'", f"{uy_rms_avg / uz_rms_avg if uz_rms_avg else np.nan:.3f}"),
    ]

    block_rms = [
        ("ux RMS avg", f"{ux_rms_avg:.3e}"),
        ("uy RMS avg", f"{uy_rms_avg:.3e}"),
        ("uz RMS avg", f"{uz_rms_avg:.3e}"),
        ("TKE avg", f"{np.mean(tke_profile):.3e}"),
        ("Turbulent intensity avg", f"{np.mean(intensity_profile):.3e}")
    ] + anisotropy_ratios

    print_block("Turbulence RMS Statistics", block_rms)

    plot_rms_profiles(coord, ux_rms, uy_rms, uz_rms, tke_profile, intensity_profile, direction, time)

def plot_fluctuations(ux_fluct, uy_fluct, uz_fluct, x, y, z, direction, coord_value, time):
    """
    Plot 2D slices of velocity fluctuations in the specified direction.

    Parameters
    ----------
    ux_fluct, uy_fluct, uz_fluct : ndarray
        Fluctuation fields for each velocity component.
    x, y, z : ndarray
        Coordinates in each direction.
    direction : str
        Direction orthogonal to the plane to extract (x, y, or z).
    coord_value : float
        Coordinate value at which to extract the slice.
    time : float
        Simulation time.
    """

    dir_map = {
        'x': (x, 'x', ('z', 'y')),
        'y': (y, 'y', ('z', 'x')),
        'z': (z, 'z', ('y', 'x')),
    }

    coords, label_dir, (label1, label2) = dir_map[direction]

    # Index of the closest coordinate value
    index = np.argmin(np.abs(coords - coord_value))

    if direction == 'x':
        uxp = ux_fluct[index, :, :]
        uyp = uy_fluct[index, :, :]
        uzp = uz_fluct[index, :, :]
        X, Y = np.meshgrid(z, y, indexing='ij')
    elif direction == 'y':
        uxp = ux_fluct[:, index, :]
        uyp = uy_fluct[:, index, :]
        uzp = uz_fluct[:, index, :]
        X, Y = np.meshgrid(z, x, indexing='ij')
    elif direction == 'z':
        uxp = ux_fluct[:, :, index]
        uyp = uy_fluct[:, :, index]
        uzp = uz_fluct[:, :, index]
        X, Y = np.meshgrid(x, y, indexing='ij')

    fig, axs = plt.subplots(1, 3, figsize=(16, 9), constrained_layout=True)

    fields = [(uxp, "u'_x"), (uyp, "u'_y"), (uzp, "u'_z")]
    for ax, (field, label) in zip(axs, fields):
        if (direction == 'z'):
            im = ax.pcolormesh(X, Y, field, shading='auto', cmap='coolwarm')
            ax.set_xlabel(label2)
            ax.set_ylabel(label1)
        else:
            im = ax.pcolormesh(X, Y, field.T, shading='auto', cmap='coolwarm')
            ax.set_xlabel(label1)
            ax.set_ylabel(label2)
        ax.set_title(label)
        fig.colorbar(im, ax=ax)

    title = f"Velocity fluctuations in the {label1}-{label2} plane at time t = {time:.2f} and {label_dir} = {coord_value:.2f}"
    fig.suptitle(title, fontsize=14)
    plt.show()

def compute_isotropic_spectrum(ux, uy, uz, dx, dy, dz, nbins=150):
    """
    Compute the isotropic 3D energy spectrum E(k) and 1D spectra Ex(kx), Ey(ky), Ez(kz).

    Returns
    -------
    k_bin_centers : ndarray
        Midpoints of the k-bins.
    E_k : ndarray
        Isotropic energy spectrum E(k).
    kx, Ex : ndarray
        Wavenumbers in x and corresponding Ex(kx)
    ky, Ey : ndarray
        Wavenumbers in y and corresponding Ey(ky)
    kz, Ez : ndarray
        Wavenumbers in z and corresponding Ez(kz)
    """
    nx, ny, nz = ux.shape

    # FFTs
    uxf = np.fft.fftn(ux)
    uyf = np.fft.fftn(uy)
    uzf = np.fft.fftn(uz)

    # Grid of wave numbers
    kx = np.fft.fftfreq(nx, d=dx) * 2 * np.pi
    ky = np.fft.fftfreq(ny, d=dy) * 2 * np.pi
    kz = np.fft.fftfreq(nz, d=dz) * 2 * np.pi

    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    K_mag = np.sqrt(KX**2 + KY**2 + KZ**2)

    # Spectral energy density (normalized)
    E_local = 0.5 * (np.abs(uxf)**2 + np.abs(uyf)**2 + np.abs(uzf)**2) / (nx * ny * nz)

    # Isotropic spectrum E(k)
    k_flat = K_mag.ravel()
    E_flat = E_local.ravel()
    k_max = np.max(k_flat)
    k_bins = np.linspace(0, k_max, nbins + 1)
    E_k, _ = np.histogram(k_flat, bins=k_bins, weights=E_flat)
    N_k, _ = np.histogram(k_flat, bins=k_bins)
    E_k = np.divide(E_k, N_k, out=np.zeros_like(E_k), where=N_k > 0)
    k_bin_centers = 0.5 * (k_bins[:-1] + k_bins[1:])

    # 1D spectra by averaging over other directions
    # Ex(kx) = mean over ky,kz of 0.5*|û_x|^2
    Ex_kx = 0.5 * np.mean(np.abs(uxf)**2, axis=(1,2)) / (nx * ny * nz)
    Ey_ky = 0.5 * np.mean(np.abs(uyf)**2, axis=(0,2)) / (nx * ny * nz)
    Ez_kz = 0.5 * np.mean(np.abs(uzf)**2, axis=(0,1)) / (nx * ny * nz)

    # Return all
    return k_bin_centers, E_k, kx, Ex_kx, ky, Ey_ky, kz, Ez_kz


def plot_isotropic_spectrum(k, E_k, kx, Ex, ky, Ey, kz, Ez, time, les_active=False, delta=None, dx=None):
    plt.figure(figsize=(16, 9))
    
    # Isotropic 3D spectrum
    plt.loglog(k, E_k, label='Isotropic E(k)')
    
    # 1D spectra
    plt.loglog(np.abs(kx), Ex, label='Ex(kx)')
    plt.loglog(np.abs(ky), Ey, label='Ey(ky)')
    plt.loglog(np.abs(kz), Ez, label='Ez(kz)')

    # Référence -5/3
    if np.count_nonzero(E_k) > 5:
        valid = E_k > 0
        k_ref = k[valid][1:-1]
        E_ref = E_k[valid][1] * (k_ref / k_ref[0])**(-5/3)
        plt.loglog(k_ref, E_ref, '--', label='-5/3 slope')

    # Si LES, tracer lignes verticales
    if les_active and delta is not None and dx is not None:
        k_filter = 2 * np.pi / delta
        k_min_resolved = np.pi / dx  # plus petite taille résolue (échelle de Nyquist)
        plt.axvline(k_filter, color='red', linestyle=':', label='Filter scale (delta)')
        plt.axvline(k_min_resolved, color='green', linestyle=':', label='Smallest resolved scale')

    plt.xlabel("Wavenumber |k|")
    plt.ylabel("E(k)")
    plt.title(f"Isotropic Energy Spectrum at t = {time:.2f}")
    plt.grid(True, which='both', ls=':')
    plt.legend()
    plt.tight_layout()
    plt.show()

def analyze_les(ux, uy, uz, dx, dy, dz, cs, Re):
    """
    Evaluate the Smagorinsky LES model performance.

    Parameters
    ----------
    ux, uy, uz : ndarray
        Velocity components (3D).
    dx, dy, dz : float
        Grid spacings.
    cs : float
        Smagorinsky constant.
    Re : float
        Reynolds number (dimensionless viscosity = 1/Re).
    """
    delta = (dx * dy * dz) ** (1.0 / 3.0)
    nu_molecular = 1.0 / Re
    nu_t = np.zeros_like(ux)
    epsilon_sgs = np.zeros_like(ux)

    for u in [ux, uy, uz]:
        dudx = np.gradient(u, dx, axis=0, edge_order=2)
        dudy = np.gradient(u, dy, axis=1, edge_order=2)
        dudz = np.gradient(u, dz, axis=2, edge_order=2)
        S = np.sqrt(2 * (dudx**2 + dudy**2 + dudz**2))
        nu_t += (cs**2) * delta**2 * S
        epsilon_sgs += (cs**2) * delta**2 * S**3  # nu_t * S^2

    nu_t /= 3.0
    epsilon_sgs /= 3.0

    nu_t_avg = np.mean(nu_t)
    nu_ratio = nu_t.max() / nu_molecular
    epsilon_avg = np.mean(epsilon_sgs)

    print_block("LES Smagorinsky Evaluation", [
        ("Samgoransky constant Cs", f"{cs:.3f}"),
        ("Filter scale Delta", f"{delta:.3e}"),
        ("Molecular viscosity (1/Re)", f"{nu_molecular:.3e}"),
        ("nu_t min/max", f"{nu_t.min():.3e} / {nu_t.max():.3e}"),
        ("nu_t mean", f"{nu_t_avg:.3e}"),
        ("nu_t / nu", f"{nu_ratio:.2f}"),
        ("SGS dissipation <epsilon_sgs>", f"{epsilon_avg:.3e}"),
    ])

    # ASCII diagnostics
    print("-" * 48)
    if nu_ratio < 0.5:
        print("WARNING: nu_t is small compared to molecular viscosity.")
        print("         (Possible under-resolution or low shear.)")
    elif nu_ratio > 100:
        print("WARNING: nu_t is large compared to molecular viscosity.")
        print("         (Possible over-dissipation, high Cs, or coarse grid.)")
    else:
        print("Model check: nu_t / nu in typical LES range.")
    print("-" * 48)

def main():
    """
    Main routine for parsing arguments and running requested analyses.
    """
    parser = argparse.ArgumentParser(
        description="Post-processing and summary for OSINCO3D binary output."
    )

    parser.add_argument('--x', type=float, default=0.0, help="x coordinate for profile extraction")
    parser.add_argument('--y', type=float, default=0.0, help="y coordinate for profile extraction")
    parser.add_argument('--z', type=float, default=0.0, help="z coordinate for profile extraction")

    parser.add_argument('--filename', '-fi', nargs='*',
        default=glob.glob("fields_??????.bin"),
        help="Binary files to analyze (wildcards allowed).")

    parser.add_argument('--summary', '-s', action='store_true',
        help="Print summary of the fields (time, min/max of velocity and pressure)")

    parser.add_argument('--rms', action='store_true',
        help="Compute and plot RMS velocity profiles and turbulent statistics")

    parser.add_argument('--fluctuations', '-f', action='store_true',
        help="Plot 2D fields of velocity fluctuations u' in a plane orthogonal to direction (x/y/z).")

    parser.add_argument('--direction', choices=['x', 'y', 'z'], default='x',
        help="Direction along which RMS profiles are computed (default: x)")

    parser.add_argument('--turbulence', '-t', action='store_true',
        help="Perform turbulence analysis: energy spectra E(k), Kolmogorov scaling.")

    parser.add_argument('--les', action='store_true',
        help="Activate LES diagnostics for the Smagorinsky model.")

    parser.add_argument('-cs', type=float, default=0.17,
        help="Smagorinsky constant Cs (default: 0.17)")

    parser.add_argument('--reynolds', '-re', type=float, default=None,
        help="Reynolds number for LES analysis.")

    parser.add_argument('--nbins', type=int, default=150,
        help="Number of bins for the isotropic energy spectrum E(k). Default: 150")

    args = parser.parse_args()

    file_list = list(itertools.chain.from_iterable(glob.glob(f) for f in args.filename))
    file_list = sorted(set(file_list))

    if not file_list:
        print("No file found matching the given patterns.")
        return 1

    for fname in file_list:
        print(f"\n========= Analyzing {fname} =========")
        time, x, y, z, ux, uy, uz, pp = read_fields(fname)

        if args.summary:
            block = [
                ("File", fname),
                ("Time", f"{time:.1f}"),
                ("nx, ny, nz", f"{len(x)}, {len(y)}, {len(z)}"),
                ("Lx, Ly, Lz", f"{x[-1]-x[0]:.1f}, {y[-1]-y[0]:.1f}, {z[-1]-z[0]:.1f}"),
                ("ux min/max", f"{ux.min():.3e} / {ux.max():.3e}"),
                ("uy min/max", f"{uy.min():.3e} / {uy.max():.3e}"),
                ("uz min/max", f"{uz.min():.3e} / {uz.max():.3e}"),
                ("Pressure min", f"{pp.min():.3e}"),
                ("Pressure max", f"{pp.max():.3e}")
            ]
            print_block("Field Summary", block)

        if args.rms:
            coords = {'x': x, 'y': y, 'z': z}
            analyze_rms_and_tke(ux, uy, uz, coords, args.direction, time)

        if args.fluctuations:
            axis_dict = {'x': 0, 'y': 1, 'z': 2}
            profile_axis = axis_dict[args.direction]
            mean_axes = tuple(ax for ax in (0, 1, 2) if ax != profile_axis)

            ux_rms, ux_fluct = compute_rms_profiles(ux, mean_axes, profile_axis)
            uy_rms, uy_fluct = compute_rms_profiles(uy, mean_axes, profile_axis)
            uz_rms, uz_fluct = compute_rms_profiles(uz, mean_axes, profile_axis)

            coord_value = {'x': args.x, 'y': args.y, 'z': args.z}[args.direction]
            plot_fluctuations(ux_fluct, uy_fluct, uz_fluct, x, y, z, args.direction, coord_value, time)

        if args.turbulence:
            dx_ = x[1] - x[0]
            dy_ = y[1] - y[0]
            dz_ = z[1] - z[0]
            k, E, kx, Ex, ky, Ey, kz, Ez = compute_isotropic_spectrum(ux, uy, uz, dx_, dy_, dz_, nbins=args.nbins)
            plot_isotropic_spectrum(k, E, kx, Ex, ky, Ey, kz, Ez, time, les_active=args.les, delta=(dx_*dy_*dz_)**(1/3), dx=max(dx_,dy_,dz_))

        if args.les:
            if args.reynolds is None:
                print("Please provide --reynolds for the LES analysis.")
            else:
                dx_ = x[1] - x[0]
                dy_ = y[1] - y[0]
                dz_ = z[1] - z[0]
                analyze_les(ux, uy, uz, dx_, dy_, dz_, args.cs, args.reynolds)

    return 0

if __name__ == "__main__":
    sys.exit(main())

