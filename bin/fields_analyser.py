#!/bin/python3
# encoding: utf-8

import os
import sys
import glob
import argparse
import itertools

import numpy as np
import matplotlib.pyplot as plt

def print_block(title, entries):
    """Print a well-formatted block of labeled entries."""
    print(f"\n--- {title} ---")
    for label, value in entries:
        print(f"  {label:<30}: {value}")
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

def plot_mixing_layer_profiles(x, y, z, ux, uy, uz, x_plot, z_plot, time, re_init):
    """
    Analyze and plot the mixing layer at a given (x,z) position:
    - Mean velocity profiles (averaged over x and z)
    - Local velocity line at (x_plot, z_plot)
    - Fluctuations relative to the mean profiles
    "
    """
    threshold_ratio=0.05
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
    # Nombre de points dans la couche de mélange
    mixing_layer_mask = (y >= y_lower) & (y <= y_upper)
    n_points_mixing_layer = np.count_nonzero(mixing_layer_mask)

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
        ("Re_delta_theta", f"{re_theta:.1f}"),
        ("N_points in layer", f"{n_points_mixing_layer}")
    ])

    # Plot <u_x>
    axs[0].plot(ux_profile, y, '--o', label=r'$\langle u_x \rangle$', color='tab:blue', 
            markersize=4, markerfacecolor='none', linewidth=0.8)
    draw_mixing_lines(axs[0], show_legend=True, delta=delta, delta_theta=delta_theta)
    axs[0].set_title(r'Profile $\langle u_x(y) \rangle$', fontsize=14)
    axs[0].set_xlabel(r'$u_x$', fontsize=12,)
    axs[0].set_ylabel(r'$y$', fontsize=12)
    axs[0].set_ylim(y_min, y_max)
    axs[0].grid(True, linestyle='--', alpha=0.7)
    axs[0].legend()

    # Plot <u_y> and <u_z>
    axs[1].plot(uy_profile, y, '--o', label=r'$\langle u_y \rangle$', color='tab:green', 
            markersize=4, markerfacecolor='none', linewidth=0.8)
    axs[1].plot(uz_profile, y, '--o', label=r'$\langle u_z \rangle$', color='tab:orange',
            markersize=4, markerfacecolor='none', linewidth=0.8)
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

    return

def plot_velocity_fluctuations_plan(x, y, z, ux, uy, uz, 
        x_plot, y_plot, z_plot, time, direction='x'):
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
        coord_plot = ['x', x_plot]
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
        coord_plot = ['y', y_plot]
    elif direction == 'z':
        iz = np.argmin(np.abs(z - z_plot))
        ux_plane = ux[:, :, iz]
        uy_plane = uy[:, :, iz]
        uz_plane = uz[:, :, iz]
        ux_mean = np.mean(ux, axis=2)
        uy_mean = np.mean(uy, axis=2)
        uz_mean = np.mean(uz, axis=2)
        coord1, coord2 = y, x
        label1, label2 = 'y', 'x'
        coord_plot = ['z', z_plot]
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
    if (direction == 'z'):
        extent = [coord2.min(), coord2.max(), coord1.min(), coord1.max()]
        im0 = axs[0].imshow(fluct_x.T, extent=extent, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
        axs[0].set_title(rf"Fluctuation $u'_x({label1},{label2})$")
        axs[0].set_xlabel(f"{label2}")
        axs[0].set_ylabel(f"{label1}")
        fig.colorbar(im0, ax=axs[0])

        im1 = axs[1].imshow(fluct_y.T, extent=extent, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
        axs[1].set_title(rf"Fluctuation $u'_y({label1},{label2})$")
        axs[1].set_xlabel(f"{label2}")
        axs[1].set_ylabel(f"{label1}")
        fig.colorbar(im1, ax=axs[1])

        im2 = axs[2].imshow(fluct_z.T, extent=extent, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
        axs[2].set_title(rf"Fluctuation $u'_z({label1},{label2})$")
        axs[2].set_xlabel(f"{label2}")
        axs[2].set_ylabel(f"{label1}")
        fig.colorbar(im2, ax=axs[2])
    else:
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
            rf"Velocity fluctuations in the ${label1}$-${label2}$ plane at time $t = {time:.2f}$ and ${coord_plot[0]:1s} = {coord_plot[1]:.2f}$",
        fontsize=16
    )
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    
    plt.show()

def divergence(ux, uy, uz, x, y, z, time,
               bc_x='periodic', bc_y='free_slip', bc_z='periodic', 
               x_plot=None, y_plot=None, z_plot=None):
    """
    Compute divergence of velocity field and plot div(u) in planes.
    ux, uy, uz : 3D arrays (nx, ny, nz)
    x, y, z : coordinate arrays
    bc_x/y/z : 'periodic' or 'periodic'
    x_plot, y_plot, z_plot : coordinates at which to plot planes (default: center)
    """

    nx, ny, nz = ux.shape
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]

    # Allocate divergence array
    div = np.zeros_like(ux)

    # --- Derivatives --- #
    # x-direction
    if bc_x == 'periodic':
        dux_dx = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) / (2*dx)
    else:  # not periodic
        dux_dx = np.zeros_like(ux)
        dux_dx[1:-1,:,:] = (ux[2:,:,:] - ux[:-2,:,:]) / (2*dx)
        dux_dx[0,:,:] = (ux[1,:,:] - ux[0,:,:]) / dx
        dux_dx[-1,:,:] = (ux[-1,:,:] - ux[-2,:,:]) / dx

    # y-direction
    if bc_y == 'periodic':
        duy_dy = (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1)) / (2*dy)
    else:
        duy_dy = np.zeros_like(uy)
        duy_dy[:,1:-1,:] = (uy[:,2:,:] - uy[:,:-2,:]) / (2*dy)
        duy_dy[:,0,:] = (uy[:,1,:] - uy[:,0,:]) / dy
        duy_dy[:,-1,:] = (uy[:,-1,:] - uy[:,-2,:]) / dy

    # z-direction
    if bc_z == 'periodic':
        duz_dz = (np.roll(uz, -1, axis=2) - np.roll(uz, 1, axis=2)) / (2*dz)
    else:
        duz_dz = np.zeros_like(uz)
        duz_dz[:,:,1:-1] = (uz[:,:,2:] - uz[:,:,:-2]) / (2*dz)
        duz_dz[:,:,0] = (uz[:,:,1] - uz[:,:,0]) / dz
        duz_dz[:,:,-1] = (uz[:,:,-1] - uz[:,:,-2]) / dz

    # --- Divergence --- #
    div = dux_dx + duy_dy + duz_dz

    div_abs = np.abs(div)
    max_div = np.max(div_abs)
    idx_max = np.unravel_index(np.argmax(div_abs), div.shape)
    mean_div = np.mean(div_abs)
    std_div = np.std(div_abs)

    print_block("Divergence analysis", [
        ("max(|div(u)|)", f"{max_div:.4e}"),
        ("Coordinates of max", f"x={x[idx_max[0]]:.3f}, y={y[idx_max[1]]:.3f}, z={z[idx_max[2]]:.3f}"),
        ("mean(|div(u)|)", f"{mean_div:.4e}"),
        ("std(|div(u)|)", f"{std_div:.4e}")
    ])

    # --- Determine plot positions --- #
    ix = np.argmin(np.abs(x - x_plot)) if x_plot is not None else nx//2
    iy = np.argmin(np.abs(y - y_plot)) if y_plot is not None else ny//2
    iz = np.argmin(np.abs(z - z_plot)) if z_plot is not None else nz//2

    planes = {
        'x-y': div[:, :, iz],
        'x-z': div[:, iy, :],
        'y-z': div[ix, :, :]
    }
    coords = {
        'x-y': (x, y),
        'x-z': (x, z),
        'y-z': (y, z)
    }

    # --- Plotting --- #
    fig, axs = plt.subplots(1, 3, figsize=(16, 9))
    cmap = 'seismic'
    vabs = np.max(np.abs(div))
    vmin, vmax = -vabs, vabs

    for i, key in enumerate(['x-y', 'x-z', 'y-z']):
        if (key == 'y-z'):
            c2, c1 = coords[key]
            im = axs[i].imshow(
                planes[key], origin='lower', extent=[c2[0], c2[-1], c1[0], c1[-1]],
                cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto'
            )
            axs[i].set_xlabel(key.split('-')[1])
            axs[i].set_ylabel(key.split('-')[0])
            axs[i].set_title(f"div(u) in {key} plane")
            fig.colorbar(im, ax=axs[i])
        else:
            c1, c2 = coords[key]
            im = axs[i].imshow(
                planes[key].T, origin='lower', extent=[c1[0], c1[-1], c2[0], c2[-1]],
                cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto'
            )
            axs[i].set_xlabel(key.split('-')[0])
            axs[i].set_ylabel(key.split('-')[1])
            axs[i].set_title(f"div(u) in {key} plane")
            fig.colorbar(im, ax=axs[i])

    fig.suptitle(rf"Divergence plots at $t = {time:.2f}$ and $x={x[ix]:.2f}, y={y[iy]:.2f}, z={z[iz]:.2f}$", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.show()

    return div

def turbulence_rms_analysis(x, y, z, ux, uy, uz, time):
    """
    Compute RMS turbulence profiles along y and extract global RMS indicators.
    """
    # Mean profiles (average over x,z)
    ux_mean = np.mean(ux, axis=(0, 2))
    uy_mean = np.mean(uy, axis=(0, 2))
    uz_mean = np.mean(uz, axis=(0, 2))

    # Fluctuations
    ux_fluct = ux - ux_mean[np.newaxis, :, np.newaxis]
    uy_fluct = uy - uy_mean[np.newaxis, :, np.newaxis]
    uz_fluct = uz - uz_mean[np.newaxis, :, np.newaxis]

    # RMS per y
    ux_rms_y = np.sqrt(np.mean(ux_fluct**2, axis=(0, 2)))
    uy_rms_y = np.sqrt(np.mean(uy_fluct**2, axis=(0, 2)))
    uz_rms_y = np.sqrt(np.mean(uz_fluct**2, axis=(0, 2)))

    # Max values and corresponding y
    ix_max = np.argmax(ux_rms_y)
    iy_max = np.argmax(uy_rms_y)
    iz_max = np.argmax(uz_rms_y)
    yx_peak, yy_peak, yz_peak = y[ix_max], y[iy_max], y[iz_max]

    ux_rms_max = ux_rms_y[ix_max]
    uy_rms_max = uy_rms_y[iy_max]
    uz_rms_max = uz_rms_y[iz_max]

    # Turbulent kinetic energy and sqrt(2k)
    sqrt2k = np.sqrt(ux_rms_y**2 + uy_rms_y**2 + uz_rms_y**2)
    sqrt2k_max = np.max(sqrt2k)

    # Anisotropy ratios
    ratio_yz = uy_rms_max / uz_rms_max if uz_rms_max != 0 else np.nan
    ratio_xz = ux_rms_max / uz_rms_max if uz_rms_max != 0 else np.nan

    # Print summary
    print_block("RMS turbulence summary", [
        ("Max RMS u'_x", f"{ux_rms_max:.4f}"),
        ("Max RMS u'_y", f"{uy_rms_max:.4f}"),
        ("Max RMS u'_z", f"{uz_rms_max:.4f}"),
        ("Peak y (u'_x RMS)", f"{yx_peak:.4f}"),
        ("Peak y (u'_y RMS)", f"{yy_peak:.4f}"),
        ("Peak y (u'_z RMS)", f"{yz_peak:.4f}"),
        ("Max sqrt(2k)", f"{sqrt2k_max:.4f}")
    ])

    print_block("Anisotropy indicators", [
        ("u'_y / u'_z", f"{ratio_yz:.2f}"),
        ("u'_x / u'_z", f"{ratio_xz:.2f}")
    ])

    # Plot RMS profiles
    plt.figure(figsize=(16,9))
    plt.plot(ux_rms_y, y, label=r"$u'_x$ RMS")
    plt.plot(uy_rms_y, y, label=r"$u'_y$ RMS")
    plt.plot(uz_rms_y, y, label=r"$u'_z$ RMS")
    plt.plot(sqrt2k, y, 'k--', label=r"$\sqrt{2k}$")
    plt.xlabel("Velocity RMS")
    plt.ylabel("y")
    plt.title(rf"Turbulent RMS profiles at $t = {time:.2f}$", fontsize=16)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.show()

# Corrected turbulence_spectrum function for fields_analyser.py
# Replaces the original turbulence_spectrum function. Fixes:
# - compute 1D energy spectra properly by summing velocity-component spectra
# - interpolate directional spectra onto a common k grid to build E_total
# - optional Hann window to reduce spectral leakage
# - keeps original plotting layout and diagnostics

def turbulence_spectrum(ux, uy, uz, x, y, z, time, Re,
        bc_x='periodic', bc_y='free_slip', bc_z='periodic', window=False):
    """
    Compute 1D energy spectra along each Cartesian direction, plus an isotropic
    estimate E_total(k) obtained by averaging directional spectra interpolated
    onto a common k-grid.

    Parameters
    ----------
    ux, uy, uz : ndarray
        Velocity components with shape (nx, ny, nz).
    x, y, z : 1D arrays
        Coordinates (assumed uniform spacing in each direction).
    time : float
        Physical time (for titles).
    Re : float
        Reynolds number (used to compute nu and scales).
    bc_* : str
        Boundary condition tags (kept for API compatibility; not used here).
    window : bool
        If True apply a Hann window along the FFT axis to reduce leakage.

    Returns
    -------
    k_center, E_total : ndarray
        Wavenumber grid and isotropic energy estimate on that grid.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d

    nx, ny, nz = ux.shape
    dx, dy, dz = x[1]-x[0], y[1]-y[0], z[1]-z[0]
    nu = 1.0 / Re

    # Helper: optional Hann window applied along given axis
    def apply_window(f, axis):
        if not window:
            return f
        n = f.shape[axis]
        w = np.hanning(n)
        # reshape w for broadcasting
        shape = [1]*f.ndim
        shape[axis] = n
        w = w.reshape(shape)
        return f * w

    # === 1D energy spectra along each direction ===
    def fft1d_energy_component(field, axis):
        # optionally window
        fwin = apply_window(field, axis)
        # rFFT along axis
        u_hat = np.fft.rfft(fwin, axis=axis)
        # normalization: use Parseval-consistent normalization
        # we'll divide by number of points along axis so that energy has correct units
        N = field.shape[axis]
        u_hat = u_hat / np.sqrt(N)
        axes_avg = tuple(i for i in range(field.ndim) if i != axis)
        # mean over the other axes
        return 0.5 * np.mean(np.abs(u_hat)**2, axis=axes_avg)

    # Compute spectra for each velocity component along each axis
    E_kx_x = fft1d_energy_component(ux, axis=0)
    E_kx_y = fft1d_energy_component(uy, axis=0)
    E_kx_z = fft1d_energy_component(uz, axis=0)
    E_kx = E_kx_x + E_kx_y + E_kx_z

    E_ky_x = fft1d_energy_component(ux, axis=1)
    E_ky_y = fft1d_energy_component(uy, axis=1)
    E_ky_z = fft1d_energy_component(uz, axis=1)
    E_ky = E_ky_x + E_ky_y + E_ky_z

    E_kz_x = fft1d_energy_component(ux, axis=2)
    E_kz_y = fft1d_energy_component(uy, axis=2)
    E_kz_z = fft1d_energy_component(uz, axis=2)
    E_kz = E_kz_x + E_kz_y + E_kz_z

    # Wavenumber axes
    kx = 2*np.pi*np.fft.rfftfreq(nx, dx)
    ky = 2*np.pi*np.fft.rfftfreq(ny, dy)
    kz = 2*np.pi*np.fft.rfftfreq(nz, dz)

    # Build a common k grid (use kx by default, but avoid zero-frequency issues)
    k_center = kx

    # Interpolate ky,kz spectra onto kx grid (extrapolate with constant fill)
    interp_Eky = interp1d(ky, E_ky, bounds_error=False, fill_value=(E_ky[0], E_ky[-1]))
    interp_Ekz = interp1d(kz, E_kz, bounds_error=False, fill_value=(E_kz[0], E_kz[-1]))
    E_total = (E_kx + interp_Eky(k_center) + interp_Ekz(k_center)) / 3.0

    # === Derivatives with BCs (kept from original, used to compute enstrophy/dissipation) ===
    def grad(f, axis, d, bc):
        df = np.zeros_like(f)
        if bc == 'periodic':
            df = (np.roll(f, -1, axis=axis) - np.roll(f, 1, axis=axis)) / (2*d)
        else:
            sl = [slice(None)] * f.ndim
            sl[axis] = slice(1, -1)
            df[tuple(sl)] = (np.take(f, range(2, f.shape[axis]), axis=axis)
                            - np.take(f, range(0, f.shape[axis]-2), axis=axis)) / (2*d)
            if axis == 0:
                df[0,:,:] = (f[1,:,:] - f[0,:,:]) / d
                df[-1,:,:] = (f[-1,:,:] - f[-2,:,:]) / d
            elif axis == 1:
                df[:,0,:] = (f[:,1,:] - f[:,0,:]) / d
                df[:,-1,:] = (f[:,-1,:] - f[:,-2,:]) / d
            elif axis == 2:
                df[:,:,0] = (f[:,:,1] - f[:,:,0]) / d
                df[:,:,-1] = (f[:,:,-1] - f[:,:,-2]) / d
        return df

    duxdx = grad(ux, 0, dx, bc_x); duxdy = grad(ux, 1, dy, bc_y); duxdz = grad(ux, 2, dz, bc_z)
    duydx = grad(uy, 0, dx, bc_x); duydy = grad(uy, 1, dy, bc_y); duydz = grad(uy, 2, dz, bc_z)
    duzdx = grad(uz, 0, dx, bc_x); duzdy = grad(uz, 1, dy, bc_y); duzdz = grad(uz, 2, dz, bc_z)

    # Vorticity components (correct expressions)
    omega_x = duzdy - duydz
    omega_y = duxdz - duzdx
    omega_z = duydx - duxdy

    enstrophy = 0.5*np.mean(omega_x**2 + omega_y**2 + omega_z**2)
    epsilon = 2 * nu * enstrophy
    # avoid division by zero
    if epsilon <= 0:
        eta = np.nan
    else:
        eta = (nu**3 / epsilon)**0.25
    delta = min(dx, dy, dz)
    delta_over_eta = delta / eta if eta and (not np.isnan(eta)) else np.nan

    # Integral scale (approx)
    Lx = x[-1]-x[0]; Ly = y[-1]-y[0]; Lz = z[-1]-z[0]
    L_int = (Lx * Ly * Lz)**(1/3)

    # Fit inertial range slope (choose valid mask)
    nk1, nk2 = 2, 10
    mask = (k_center > nk1) & (k_center < nk2) & (~np.isnan(E_total)) & (E_total>0)
    slope_total, intercept_total = (np.nan, np.nan)
    if np.count_nonzero(mask) >= 2:
        slope_total, intercept_total = np.polyfit(np.log10(k_center[mask]), np.log10(E_total[mask]), 1)

    # Diagnostics print
    print_block("Turbulence spectrum diagnostics", [
        ("nu (viscosity)", f"{nu:.2e}"),
        ("Mean enstrophy", f"{enstrophy:.5f}"),
        ("Dissipation epsilon", f"{epsilon:.5e}"),
        ("Kolmogorov scale eta", f"{eta:.5f}"),
        ("Delta (grid scale)", f"{delta:.5f}"),
        ("Delta / eta", f"{delta_over_eta:.2f}"),
        (f"Inertial range slope (k={nk1}-{nk2})", f"{slope_total:.2f}" if not np.isnan(slope_total) else "N/A"),
        ("", "(expected: -1.67)")
    ])

    # === Plot spectra (preserve original layout) ===
    plt.figure(figsize=(16,9))

    # (1) Raw spectra
    plt.subplot(1,2,1)
    plt.loglog(kx, E_kx, label=r'$E(k_x)$', linestyle='dotted')
    plt.loglog(ky, E_ky, label=r'$E(k_y)$', linestyle='dotted')
    plt.loglog(kz, E_kz, label=r'$E(k_z)$', linestyle='dotted')
    plt.loglog(k_center, E_total, 'k', linewidth=2, label=r'$E_{total}(k)$')

    # Reference line k^-5/3 (use magnitude close to E_total at nk1 if available)
    k_ref = np.linspace(nk1, nk2+12, 10)
    if np.any((~np.isnan(E_total)) & (k_center>0)):
        idx_ref = np.argmin(np.abs(k_center - nk1))
        C = max(E_total[idx_ref], 1e-16) * (nk1**(5/3))
        plt.loglog(k_ref, C*k_ref**(-5/3), 'k--', label=r'$k^{-5/3}$')

    plt.axvline(1/L_int, color='green', linestyle='--', label=r'$k_L$ (large scale)')
    plt.axvline(np.pi/delta, color='red', linestyle='--', label=r'$k_\Delta$ (LES cutoff)')
    if (not np.isnan(eta)):
        plt.axvline(1/eta, color='orange', linestyle='--', label=r'$k_\eta$ (Kolmogorov)')

    if not np.isnan(slope_total):
        k_fit_line = np.linspace(nk1, nk2, 100)
        E_fit_line = 10**intercept_total * k_fit_line**slope_total
        plt.loglog(k_fit_line, E_fit_line, 'r--', linewidth=2, label=rf"Measured slope = {slope_total:.2f}")

    plt.xlabel(r"Wavenumber $k$")
    plt.ylabel("Energy")
    plt.title("1D Energy Spectra")
    plt.grid(True, which='both', ls='--', alpha=0.5)
    plt.legend()

    # (2) Compensated spectra
    plt.subplot(1,2,2)
    # avoid k=0 when compensating
    mask_k = kx > 0
    if np.any(mask_k):
        plt.loglog(kx[mask_k], E_kx[mask_k]*kx[mask_k]**(5/3), label=r'$k_x^{5/3} E(k_x)$')
    mask_k = ky > 0
    if np.any(mask_k):
        plt.loglog(ky[mask_k], E_ky[mask_k]*ky[mask_k]**(5/3), label=r'$k_y^{5/3} E(k_y)$')
    mask_k = kz > 0
    if np.any(mask_k):
        plt.loglog(kz[mask_k], E_kz[mask_k]*kz[mask_k]**(5/3), label=r'$k_z^{5/3} E(k_z)$')
    mask_kc = k_center > 0
    if np.any(mask_kc):
        plt.loglog(k_center[mask_kc], E_total[mask_kc]*k_center[mask_kc]**(5/3), 'k', linewidth=2, label=r'$k^{5/3} E_{total}(k)$')

    plt.xlabel(r"Wavenumber $k$")
    plt.ylabel("Compensated Energy")
    plt.title("Compensated Spectra")
    plt.grid(True, which='both', ls='--', alpha=0.5)
    plt.legend()

    plt.suptitle(rf"1D Turbulence Energy Spectrum at $t = {time:.2f}$", fontsize=16)
    plt.tight_layout(rect=[0,0,1,0.95])
    plt.show()

    return k_center, E_total

# ============================================================
# --- Turbulence scales analysis ---
# ============================================================

def turbulence_scales_analysis(ux, uy, uz, x, y, z, re, bcx='periodic',
                               bcy='free_slip', bcz='periodic'):
    """
    Compute key turbulence scales: Taylor microscale, Re_lambda, Kolmogorov scale, etc.
    Handles periodic and free-slip boundary conditions.
    """
    nu = 1. / re
    # === Helper for derivative with boundary conditions ===
    def ddx(f, d, axis, bc):
        if bc == 'periodic':
            fp = np.roll(f, -1, axis=axis)
            fm = np.roll(f,  1, axis=axis)
            df = (fp - fm) / (2*d)
        else:
            # zero gradient at walls → Neumann condition
            df = np.zeros_like(f)
            if axis == 0:
                df[1:-1,:,:] = (f[2:,:,:] - f[:-2,:,:]) / (2*d)
                df[0,:,:] = (f[1,:,:] - f[0,:,:]) / d
                df[-1,:,:] = (f[-1,:,:] - f[-2,:,:]) / d
            elif axis == 1:
                df[:,1:-1,:] = (f[:,2:,:] - f[:,:-2,:]) / (2*d)
                df[:,0,:] = (f[:,1,:] - f[:,0,:]) / d
                df[:,-1,:] = (f[:,-1,:] - f[:,-2,:]) / d
            elif axis == 2:
                df[:,:,1:-1] = (f[:,:,2:] - f[:,:,:-2]) / (2*d)
                df[:,:,0] = (f[:,:,1] - f[:,:,0]) / d
                df[:,:,-1] = (f[:,:,-1] - f[:,:,-2]) / d
        return df

    # === Remove mean velocity ===
    ux_p = ux - np.mean(ux)
    uy_p = uy - np.mean(uy)
    uz_p = uz - np.mean(uz)

    # === RMS velocity ===
    u_rms = np.sqrt(np.mean(ux_p**2 + uy_p**2 + uz_p**2) / 3)
    
    # === Velocity gradients ===
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    dudx = ddx(ux_p, dx, 0, bcx)
    dudy = ddx(ux_p, dy, 1, bcy)
    dudz = ddx(ux_p, dz, 2, bcz)

    dvdx = ddx(uy_p, dx, 0, bcx)
    dvdy = ddx(uy_p, dy, 1, bcy)
    dvdz = ddx(uy_p, dz, 2, bcz)

    dwdx = ddx(uz_p, dx, 0, bcx)
    dwdy = ddx(uz_p, dy, 1, bcy)
    dwdz = ddx(uz_p, dz, 2, bcz)

    # === Strain rate tensor ===
    s11 = dudx
    s22 = dvdy
    s33 = dwdz
    s12 = 0.5 * (dudy + dvdx)
    s13 = 0.5 * (dudz + dwdx)
    s23 = 0.5 * (dvdz + dwdy)

    # === Dissipation ε = 2 ν <s_ij s_ij> ===
    eps = 2 * nu * np.mean(
        s11**2 + s22**2 + s33**2 + 2*(s12**2 + s13**2 + s23**2)
    )

    # === Taylor microscale ===
    lambda_taylor = np.sqrt(15 * nu * u_rms**2 / eps)

    # === Reynolds number based on lambda ===
    Re_lambda = u_rms * lambda_taylor / nu

    # === Kolmogorov scale ===
    eta = (nu**3 / eps)**0.25

    # === Integral length scale (approximation by box size) ===
    Lx = dx * ux.shape[0]
    Ly = dy * ux.shape[1]
    Lz = dz * ux.shape[2]
    L_int = (Lx * Ly * Lz)**(1/3)
    L_over_eta = L_int / eta

    print_block("Turbulence scales analysis", [
        ("u_rms", f"{u_rms:.5f}"),
        ("Dissipation", f"{eps:.5e}"),
        ("Taylor microscale", f"{lambda_taylor:.5f}"),
        ("Re_lambda", f"{Re_lambda:.2f}"),
        ("L_over_eta", f"{L_over_eta:.1f}")
    ])

def les_analysis(ux, uy, uz, x, y, z, time, cs=0.17, x_plot=None, y_plot=None, z_plot=None,
                 Re=1000.0, bc_x='periodic', bc_y='free_slip', bc_z='periodic'):
    """
    Compute Smagorinsky LES viscosity nu_t with proper boundary conditions
    and plot it in x-y, x-z, y-z planes at specified coordinates.
    Also print LES diagnostics.
    """
    nx, ny, nz = ux.shape
    dx, dy, dz = x[1]-x[0], y[1]-y[0], z[1]-z[0]
    nu = 1.0 / Re
    delta = (dx*dy*dz)**(1/3)

    def grad(f, axis, dx, bc):
        df = np.zeros_like(f)
        if bc == 'periodic':
            df = (np.roll(f, -1, axis=axis) - np.roll(f, 1, axis=axis)) / (2*dx)
        else:  # periodic
            sl = [slice(None)]*f.ndim
            sl[axis] = slice(1,-1)
            slp = [slice(None)]*f.ndim
            slp[axis] = slice(2,None)
            slm = [slice(None)]*f.ndim
            slm[axis] = slice(None,-2)
            df[tuple(sl)] = (f[tuple(slp)] - f[tuple(slm)]) / (2*dx)
            # boundaries: forward/backward difference
            sl0 = [slice(None)]*f.ndim; sl0[axis]=0
            sl1 = [slice(None)]*f.ndim; sl1[axis]=-1
            sl01 = [slice(None)]*f.ndim; sl01[axis]=1
            slm1 = [slice(None)]*f.ndim; slm1[axis]=-2
            df[tuple(sl0)] = (f[tuple(sl01)] - f[tuple(sl0)]) / dx
            df[tuple(sl1)] = (f[tuple(sl1)] - f[tuple(slm1)]) / dx
        return df

    # Compute gradients with BC
    dudx = grad(ux, 0, dx, bc_x); dudy = grad(ux, 1, dy, bc_y); dudz = grad(ux, 2, dz, bc_z)
    dvdx = grad(uy, 0, dx, bc_x); dvdy = grad(uy, 1, dy, bc_y); dvdz = grad(uy, 2, dz, bc_z)
    dwdx = grad(uz, 0, dx, bc_x); dwdy = grad(uz, 1, dy, bc_y); dwdz = grad(uz, 2, dz, bc_z)

    # Rate-of-strain tensor S_ij
    Sxx = dudx; Syy = dvdy; Szz = dwdz
    Sxy = 0.5*(dudy + dvdx); Sxz = 0.5*(dudz + dwdx); Syz = 0.5*(dvdz + dwdy)

    # Magnitude of strain tensor
    S_mag = np.sqrt(2*(Sxx**2 + Syy**2 + Szz**2 + 2*(Sxy**2 + Sxz**2 + Syz**2)))

    # Smagorinsky viscosity
    nu_t = (cs*delta)**2 * S_mag

    # Mean TKE
    tke = 0.5*(ux**2 + uy**2 + uz**2)
    mean_tke = np.mean(tke)

    # Enstrophy (using the same BCs for derivatives)
    omega_x = dwdy - dvdz
    omega_y = dwdz - dudx
    omega_z = dvdx - dudy
    enstrophy = 0.5*np.mean(omega_x**2 + omega_y**2 + omega_z**2)

    print_block("Turbulence development analysis", [
        ("Mean TKE", f"{mean_tke:.5f}"),
        ("Mean nu_t", f"{nu_t.mean():.3e}"),
        ("Max nu_t", f"{nu_t.max():.3e}"),
        ("Filter size delta", f"{delta:.5f}"),
        ("Smagorinsky constant C_s", f"{cs}")
    ])

    # Determine plot indices
    ix = np.argmin(np.abs(x - x_plot)) if x_plot is not None else nx//2
    iy = np.argmin(np.abs(y - y_plot)) if y_plot is not None else ny//2
    iz = np.argmin(np.abs(z - z_plot)) if z_plot is not None else nz//2

    planes = {
        'x-y': nu_t[:, :, iz],
        'x-z': nu_t[:, iy, :],
        'y-z': nu_t[ix, :, :]
    }
    coords = {
        'x-y': (x, y),
        'x-z': (x, z),
        'y-z': (y, z)
    }

    fig, axs = plt.subplots(1,3, figsize=(16,9))
    cmap='viridis'
    vmin, vmax = nu_t.min(), nu_t.max()

    for i, key in enumerate(['x-y','x-z','y-z']):
        c1, c2 = coords[key]
        if key == 'y-z':
            c2, c1 = c1, c2  # swap for imshow
            im = axs[i].imshow(planes[key], origin='lower', 
                extent=[c2[0], c2[-1], c1[0], c1[-1]],
                cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
            axs[i].set_xlabel(key.split('-')[1])
            axs[i].set_ylabel(key.split('-')[0])
        else:
            im = axs[i].imshow(planes[key].T, origin='lower', 
                extent=[c1[0], c1[-1], c2[0], c2[-1]],
                cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
            axs[i].set_xlabel(key.split('-')[0])
            axs[i].set_ylabel(key.split('-')[1])
        axs[i].set_title(rf"$\nu_t$ in {key} plane")
        fig.colorbar(im, ax=axs[i])

    fig.suptitle(rf"Smagorinsky viscosity $\nu_t$ at $t = {time:.2f}$", fontsize=16)
    plt.tight_layout(rect=[0,0,1,0.92])
    plt.show()

    return nu_t

def check_cfl(x, y, z, ux, uy, uz, dt):
    """
    Compute CFL numbers in each direction and print a summary.
    Warning if any CFL > 0.1
    """
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]

    cfl_x = np.max(np.abs(ux)) * dt / dx
    cfl_y = np.max(np.abs(uy)) * dt / dy
    cfl_z = np.max(np.abs(uz)) * dt / dz

    print_block("CFL Analysis", [
        ("dt", f"{dt:.5f}"),
        ("dx, dy, dz", f"{dx:.5f}, {dy:.5f}, {dz:.5f}"),
        ("CFL_x", f"{cfl_x:.4f}"),
        ("CFL_y", f"{cfl_y:.4f}"),
        ("CFL_z", f"{cfl_z:.4f}")
    ])

    cfl_max = max(cfl_x, cfl_y, cfl_z)
    if cfl_max > 0.3:
        print(f"WARNING: CFL number {cfl_max:.3f} > 0.3, consider reducing dt!")

def argsparser():
    parser = argparse.ArgumentParser(
        description="Plot velocity profiles and turbulence statistics from CFD binary output."
    )
    parser.add_argument(
        '--filename', '-fi',
        nargs='*',
        default=glob.glob("fields_??????.bin"),
        help="Fichiers binaires à analyser (wildcards autorisés). Par défaut : fields_??????.bin"
    )

    parser.add_argument('--x', type=float, default=0.0, 
            help="x coordinate for profile extraction")
    parser.add_argument('--y', type=float, default=0.0, 
            help="y coordinate for profile extraction")
    parser.add_argument('--z', type=float, default=0.0, 
            help="z coordinate for profile extraction")

    parser.add_argument('--direction', choices=['x', 'y', 'z'], default='x',
            help="Direction along which to extract the velocity profile (default: z)")

    parser.add_argument('--fluctuations', '-f', action='store_true',
            help="Display velocity fluctuations u'_i in the (y,z) plane at given x")

    parser.add_argument('--divergence', '-div', action='store_true',
            help="Calculate divergence of the velocity field and print max and mean.")
    
    parser.add_argument('--rms', action='store_true',
            help="Compute RMS turbulence statistics and anisotropy indicators.")

    parser.add_argument('--turbulence', '-turb', action='store_true',
            help="Calculate 3D energy spectra E(k), Ex(k), Ey(k), Ez(k) and plot.")

    parser.add_argument('--les', action='store_true', 
            help="Activate LES Smagorinsky analysis")
    parser.add_argument('--cs', type=float, default=0.17, 
            help="Smagorinsky constant (default 0.17)")

    parser.add_argument('--mixing_layer', '-ml', action='store_true',
            help="Activate mixing layer mode (3-panel plot in the y-direction)")

    parser.add_argument('--reynolds', '-re', type=float, default=1000.0,
            help="Reynolds number initial")

    parser.add_argument('--cfl', type=float, default=None,
            help="Calculate CFL numbers and print a warning if dt > CFL limit. Usage: --cfl 0.0005")

    return parser.parse_args()

def main():

    args = argsparser()
    file_list = list(itertools.chain.from_iterable(glob.glob(f) 
        for f in args.filename))
    file_list = sorted(set(file_list))

    for fname in file_list:
        print(f"\n========= Analyse de {fname} =========")
        basename = os.path.basename(fname)
        file_id = os.path.splitext(basename)[0].split('_')[-1]
        time, x, y, z, ux, uy, uz, pp = read_fields(fname)
        block = [
            ("File", fname),
            ("Time", f"{time:.1f}"),
            ("nx, ny, nz", f"{len(x)}, {len(y)}, {len(z)}"),
            ("Lx, Ly, Lz", f"{x[-1]-x[0]:.1f}, {y[-1]-y[0]:.1f}, {z[-1]-z[0]:.1f}"),
            ("dx, dy, dz", f"{x[1]-x[0]:.3f}, {y[1]-y[0]:.3f}, {z[1]-z[0]:.3f}"),
            ("ux min/max", f"{ux.min():.3e} / {ux.max():.3e} "),
            ("uy min/max", f"{uy.min():.3e} / {uy.max():.3e} "),
            ("uz min/max", f"{uz.min():.3e} / {uz.max():.3e} "),
            ("Pressure min", f"{pp.min():.3e}"),
            ("Pressure max", f"{pp.max():.3e}")
        ]
        print_block("Field summary", block)
        
        if args.cfl is not None:
            check_cfl(x, y, z, ux, uy, uz, args.cfl)
        
        if args.fluctuations:
            plot_velocity_fluctuations_plan(
                x, y, z, ux, uy, uz,
                args.x, args.y, args.z,
                time, direction=args.direction)

        if args.divergence:
            divergence(ux, uy, uz, x, y, z, time,
                    bc_x='periodic', bc_y='free_slip', bc_z='periodic',
                    x_plot=args.x, y_plot=args.y, z_plot=args.z)

        if args.turbulence:
            turbulence_spectrum(ux, uy, uz, x, y, z, time, args.reynolds)
            turbulence_scales_analysis(ux, uy, uz, x, y, z, args.reynolds)
            if args.rms:
                turbulence_rms_analysis(x, y, z, ux, uy, uz, time)


        if args.les:
            les_analysis(ux, uy, uz, x, y, z, time,
                 cs=args.cs, x_plot=args.x, y_plot=args.y, z_plot=args.z,
                 Re=args.reynolds, 
                 bc_x='periodic', bc_y='free_slip', bc_z='periodic')

        if args.mixing_layer:
            plot_mixing_layer_profiles(
                    x, y, z, ux, uy, uz, 
                    args.x, args.z, time, 
                    re_init=args.reynolds)

if __name__ == "__main__":

    err = main()
    sys.exit(0)
