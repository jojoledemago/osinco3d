#!/bin/python3
# encoding: utf-8

import os
import sys
import glob
import argparse
import itertools
from PIL import Image
import subprocess
from collections import defaultdict


import numpy as np
from scipy.interpolate import RegularGridInterpolator
from matplotlib.table import Table
import matplotlib.pyplot as plt

def generate_latex_figure_pdf(image_folder="png", output_dir="pdf", field_time_map=None):
    images = sorted([f for f in os.listdir(image_folder) if f.endswith(".png")])
    if not images:
        print("Warning: No PNG files found in ./png.")
        return

    os.makedirs(output_dir, exist_ok=True)
    output_tex = os.path.join(output_dir, "figures.tex")

    grouped = defaultdict(list)
    for fname in images:
        key = fname.split('_')[0]
        grouped[key].append(fname)

    with open(output_tex, "w") as f:
        f.write(r"""\documentclass[a4paper,10pt]{article}
                    \usepackage{graphicx}
                    \usepackage{caption}
                    \usepackage{subcaption}
                    \usepackage[margin=1.5cm]{geometry}
                    \usepackage{float}
                    \usepackage{parskip}
                    \begin{document}
                    """)
        for key, fnames in grouped.items():
            # Try to extract time from filename if map available
            time_str = ""
            if key in field_time_map:
                time_str = f" at t = {field_time_map[key]:.0f}"
            f.write(f"\\section*{{Figures for field {key}{time_str}}}\n")
            for i, fname in enumerate(fnames):
                if i % 2 == 0:
                    f.write("\\begin{figure}[H]\n")
                escaped_name = fname.replace('_', r'\_')
                f.write("  \\begin{subfigure}[t]{0.48\\linewidth}\n")
                f.write(f"    \\includegraphics[width=\\linewidth]{{../{image_folder}/{fname}}}\n")
                f.write(f"    \\caption{{{escaped_name}}}\n")
                f.write("  \\end{subfigure}\n")
                if i % 2 == 1 or i == len(fnames) - 1:
                    f.write("\\end{figure}\n\n")

            f.write("\\clearpage\n\n")
        f.write("\\end{document}\n")

    print(f"LaTeX file generated: {output_tex}")

    try:
        subprocess.run(
            ["pdflatex", "-output-directory", output_dir, output_tex],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        print(f"PDF generated successfully: {os.path.join(output_dir, 'figures.pdf')}")
    except FileNotFoundError:
        print("Error: pdflatex not found.")
        return
    except subprocess.CalledProcessError:
        print("Error: pdflatex compilation failed.")
        return

    # Clean up auxiliary files
    for ext in [".aux", ".log", ".out"]:
        aux_file = os.path.join(output_dir, "figures" + ext)
        if os.path.exists(aux_file):
            os.remove(aux_file)

    print("Auxiliary files cleaned up.")

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

def periodic_gradient(f, dx, axis):
    return (np.roll(f, -1, axis=axis) - np.roll(f, 1, axis=axis)) / (2 * dx)

def plot_divergence_planes(x, y, z, ux, uy, uz, x_plot, y_plot, z_plot, time, save_fig=None, fig_name=None):
    # Compute divergence of the velocity field
    div = (
        np.gradient(ux, x, axis=0, edge_order=2) +
        np.gradient(uy, y, axis=1, edge_order=2) +
        np.gradient(uz, z, axis=2, edge_order=2)
    )

    # Get global statistics
    div_abs = np.abs(div)
    div_min = np.min(div_abs)
    div_max = np.max(div_abs)
    div_mean = np.mean(div)
    div_std = np.std(div)

    # Find location of min and max divergence
    idx_min = np.unravel_index(np.argmin(div), div.shape)
    idx_max = np.unravel_index(np.argmax(div), div.shape)
    coord_min = (x[idx_min[0]], y[idx_min[1]], z[idx_min[2]])
    coord_max = (x[idx_max[0]], y[idx_max[1]], z[idx_max[2]])

    # Print diagnostic info
    print_block("Divergence diagnostics", [
        ("div(u) min", f"{div_min:.3e}"),
        ("at index", f"{idx_min[0]:3d}, {idx_min[1]:3d}, {idx_min[2]:3d}"),
        ("   coord", f"{coord_min[0]:.1f}, {coord_min[1]:.1f}, {coord_min[2]:.1f}"),
        ("div(u) max", f"{div_max:.3e}"),
        ("at index", f"{idx_max[0]:3d}, {idx_max[1]:3d}, {idx_max[2]:3d}"),
        ("   coord", f"{coord_max[0]:.1f}, {coord_max[1]:.1f}, {coord_max[2]:.1f}"),
        ("div(u) mean", f"{div_mean:.3e}"),
        ("div(u) std", f"{div_std:.3e}")
        ])

    # Select slices for plotting
    ix = np.argmin(np.abs(x - x_plot))
    iy = np.argmin(np.abs(y - y_plot))
    iz = np.argmin(np.abs(z - z_plot))

    yz = div[ix, :, :]
    xz = div[:, iy, :]
    xy = div[:, :, iz]

    # Plotting
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    cmap = 'seismic'
    vlim = np.max(np.abs(div))
    kwargs = dict(origin='lower', cmap=cmap, vmin=-vlim, vmax=vlim)

    im0 = axs[0].imshow(yz, extent=[z.min(), z.max(), y.min(), y.max()], aspect='auto', **kwargs)
    axs[0].set_title(rf"$\nabla \cdot \mathbf{{u}}$ in $z$--$y$ at $x = {x_plot:.2f}$")
    axs[0].set_xlabel(r'$z$'); axs[0].set_ylabel(r'$y$')
    plt.colorbar(im0, ax=axs[0])

    im1 = axs[1].imshow(xz.T, extent=[x.min(), x.max(), z.min(), z.max()], aspect='auto', **kwargs)
    axs[1].set_title(rf"$\nabla \cdot \mathbf{{u}}$ in $x$--$z$ at $y = {y_plot:.2f}$")
    axs[1].set_xlabel(r'$x$'); axs[1].set_ylabel(r'$z$')
    plt.colorbar(im1, ax=axs[1])

    im2 = axs[2].imshow(xy.T, extent=[x.min(), x.max(), y.min(), y.max()], aspect='auto', **kwargs)
    axs[2].set_title(rf"$\nabla \cdot \mathbf{{u}}$ in $x$--$y$ at $z = {z_plot:.2f}$")
    axs[2].set_xlabel(r'$x$'); axs[2].set_ylabel(r'$y$')
    plt.colorbar(im2, ax=axs[2])

    fig.suptitle(f"Divergence field at t = {time:.2f}", fontsize=16)
    plt.tight_layout()
    if save_fig:
        if not os.path.exists("png"):
            os.makedirs("png")
        plt.savefig(f"png/{fig_name}", dpi=300, bbox_inches='tight', transparent=False)
    else:
        plt.show()

def compute_momentum_thickness(y, u_profile, U1, U2):
    """Compute momentum thickness delta_theta for a given velocity profile."""
    delta_u = U1 - U2
    integrand = ((U1 - u_profile) / delta_u) * ((u_profile - U2) / delta_u)
    return np.trapz(integrand, y)

def plot_velocity_profiles(time, x, y, z, ux, uy, uz, x_plot, y_plot, z_plot, direction='x', save_fig=None, fig_name=None):
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

    if save_fig:
        if not os.path.exists("png"):
            os.makedirs("png")
        plt.savefig(f"png/{fig_name}", dpi=300, bbox_inches='tight', transparent=False)
    else:
        plt.show()

def plot_mixing_layer_profiles(x, y, z, ux, uy, uz, x_plot, z_plot, time, re_init, threshold_ratio=0.05, save_fig=None, fig_name=None):
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
    if save_fig:
        if not os.path.exists("png"):
            os.makedirs("png")
        plt.savefig(f"png/{fig_name}", dpi=300, bbox_inches='tight', transparent=False)
    else:
        plt.show()

    return delta, delta_theta, re_new, re_theta 

def plot_planar_jet_profiles(x, y, z, ux, uy, uz, x_plot, z_plot, time, u_jet, re_init, threshold_ratio=0.05, save_fig=None, fig_name=None):
    """
    Analyze a planar jet: compute mean and fluctuating velocity profiles,
    jet thickness (delta_j), and shear-layer momentum thickness (theta).
    """
    nu = u_jet / re_init
    # Mean velocity profiles across y (averaged over x and z)
    ux_profile = np.mean(ux, axis=(0, 2))
    uy_profile = np.mean(uy, axis=(0, 2))
    uz_profile = np.mean(uz, axis=(0, 2))

    # Extract line at (x_plot, z_plot)
    ix = np.argmin(np.abs(x - x_plot))
    iz = np.argmin(np.abs(z - z_plot))
    ux_line = ux[ix, :, iz]
    uy_line = uy[ix, :, iz]
    uz_line = uz[ix, :, iz]

    # Velocity fluctuations
    ux_fluct = ux_line - ux_profile
    uy_fluct = uy_line - uy_profile
    uz_fluct = uz_line - uz_profile

    # Estimate jet thickness using threshold (default 5% of u_jet)
    u_min = ux_profile[0]
    u_max = max(ux_profile)
    u_thresh = u_max * threshold_ratio

    y_down = y_up = None
    for j in range(len(y)):
        if ux_profile[j] > u_thresh:
            y_down = y[j]
            break
    for j in range(len(y)-1, 0, -1):
        if ux_profile[j] > u_thresh:
            y_up = y[j]
            break

    if y_up is None or y_down is None:
        raise ValueError("Jet bounds could not be determined. Adjust `threshold_ratio` or check profile.")

    delta_j = abs(y_up - y_down)

    # Compute shear-layer momentum thickness theta
    u_norm = (ux_profile - u_min) / (u_max - u_min)
    integrand = u_norm * (1 - u_norm)
    theta = np.trapz(integrand, y)

    re_delta_j = delta_j / nu
    re_theta = theta / nu

    print_block("Planar Jet Analysis", [
        ("U_jet (reference)", f"{u_jet:.3f}"),
        ("u_min / u_max", f"{u_min:.3f} / {u_max:.3f}"),
        ("y_up", f"{y_up:.3f}"),
        ("y_down", f"{y_down:.3f}"),
        ("Jet thickness delta_j", f"{delta_j:.4f}"),
        ("Momentum thickness theta", f"{theta:.4f}"),
        ("Re_delta_j", f"{re_delta_j:.1f}"),
        ("Re_theta", f"{re_theta:.1f}")
    ])

    # Plotting
    fig, axs = plt.subplots(1, 3, figsize=(16, 9), sharey=True)
    y_min, y_max = y[0], y[-1]

    def draw_jet_bounds(ax, x_=False, d_=False):
        if (x_):
            ax.axhline(y_down, color='y', linestyle='--', linewidth=1, label='Threshold bounds')
            ax.fill_betweenx(y=[y_down, y_up], x1=ax.get_xlim()[0], x2=ax.get_xlim()[1],
                            color='gray', alpha=0.1, label='Jet region')
        else:

            ax.axhline(y_down, color='y', linestyle='--', linewidth=1, label=None)
            ax.fill_betweenx(y=[y_down, y_up], x1=ax.get_xlim()[0], x2=ax.get_xlim()[1],
                            color='gray', alpha=0.1, label=None)
        ax.axhline(y_up, color='y', linestyle='--', linewidth=1)
        ax.set_ylim(y_min, y_max)
        ax.grid(True, linestyle='--', alpha=0.7)

        if d_:
            # Position under legend (bottom-left of axis)
            ax.annotate(
                rf'$\delta_j ~= {delta_j:.3f}$',
                xy=(0.75, 0.4),
                xycoords='axes fraction',
                fontsize=12,
                backgroundcolor='white',
                bbox=dict(boxstyle='round,pad=0.3', edgecolor='gray', facecolor='white')
            )


    axs[0].plot(ux_profile, y, label=r'$\langle u_x \rangle$', color='tab:blue')
    draw_jet_bounds(axs[0], True, True)
    axs[0].set_title(r'Mean profile $\langle u_x(y) \rangle$', fontsize=14)
    axs[0].set_xlabel(r'$u_x$', fontsize=12)
    axs[0].set_ylabel(r'$y$', fontsize=12)
    axs[0].legend()

    axs[1].plot(uy_profile, y, label=r'$\langle u_y \rangle$', color='tab:green')
    axs[1].plot(uz_profile, y, label=r'$\langle u_z \rangle$', color='tab:orange')
    draw_jet_bounds(axs[1])
    axs[1].set_title(r'Mean profiles $\langle u_y \rangle$, $\langle u_z \rangle$', fontsize=14)
    axs[1].set_xlabel(r'$u_i$', fontsize=12)
    axs[1].legend()

    axs[2].plot(ux_fluct, y, label=r"$u'_x$", color='tab:red')
    axs[2].plot(uy_fluct, y, label=r"$u'_y$", color='tab:purple')
    axs[2].plot(uz_fluct, y, label=r"$u'_z$", color='tab:brown')
    draw_jet_bounds(axs[2])
    axs[2].set_title(r"Velocity fluctuations $u'_i$", fontsize=14)
    axs[2].set_xlabel(r"$u'_i$", fontsize=12)
    axs[2].legend()

    fig.suptitle(
        rf"Planar jet analysis at $x = {x_plot}$, $z = {z_plot}$, $t = {time:.1f}$",
        fontsize=16
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    if save_fig:
        if not os.path.exists("png"):
            os.makedirs("png")
        plt.savefig(f"png/{fig_name}", dpi=300, bbox_inches='tight', transparent=False)
    else:
        plt.show()

    return delta_j, theta, re_delta_j, re_theta

def compute_Re_lambda(ux, uy, uz, x, y, z, nu):
    """
    Compute the Taylor-scale Reynolds number from a 3D velocity field.

    Parameters
    ----------
    ux, uy, uz : ndarray (nx, ny, nz)
        Components of the velocity field.
    x, y, z :  ndarray (nx), (ny), (nz)
        Grid coordinates in x, y, z directions.
    nu : float
        Kinematic viscosity.

    Returns
    -------
    Re_lambda : float
        Taylor-scale Reynolds number.
    """
    # Compute RMS velocity
    u2 = ux**2 + uy**2 + uz**2
    urms = np.sqrt(np.mean(u2))

    # Compute velocity gradients
    dudx = np.gradient(ux, x, axis=0, edge_order=2)
    dudy = np.gradient(ux, y, axis=1, edge_order=2)
    dudz = np.gradient(ux, z, axis=2, edge_order=2)

    dvdx = np.gradient(uy, x, axis=0, edge_order=2)
    dvdy = np.gradient(uy, y, axis=1, edge_order=2)
    dvdz = np.gradient(uy, z, axis=2, edge_order=2)

    dwdx = np.gradient(uz, x, axis=0, edge_order=2)
    dwdy = np.gradient(uz, y, axis=1, edge_order=2)
    dwdz = np.gradient(uz, z, axis=2, edge_order=2)

    # Compute strain-rate tensor squared S_ij^2
    sij2 = (
        dudx**2 + dvdy**2 + dwdz**2 +
        0.5 * ((dudy + dvdx)**2 + (dudz + dwdx)**2 + (dvdz + dwdy)**2)
    ) * 0.5

    # Compute dissipation rate
    eps = 2 * nu * np.mean(sij2)

    # Taylor microscale
    lambda_t = np.sqrt(15 * nu * urms**2 / eps)

    # Reynolds number based on Taylor microscale
    Re_lambda = urms * lambda_t / nu

    print_block("Reynolds de Taylor", [
        ("eps:", f"{eps:.3e}"),
        ("lambda:", f"{lambda_t:.3e}"),
        ("Re_lambda:", f"{Re_lambda:.1f}")
        ])

    return Re_lambda

def plot_velocity_fluctuations_plan(x, y, z, ux, uy, uz, x_plot, y_plot, z_plot, time, direction='x', save_fig=None, fig_name=None):
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
            f"Velocity fluctuations in the {label1}-{label2} plane at time t = {time:.2f} and {coord_plot[0]:1s} = {coord_plot[1]:.2f}",
        fontsize=16
    )
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    if save_fig:
        if not os.path.exists("png"):
            os.makedirs("png")
        plt.savefig(f"png/{fig_name}", dpi=300, bbox_inches='tight', transparent=False)
    else:
        plt.show()

def plot_energy_spectrum(ux, uy, uz, x, y, z, time, LES=False, normalize=False, return_data=False, re_init=1000.0, save_fig=None, fig_name=None):
    """
    Compute and plot the isotropic energy spectrum and diagnose spectral resolution.

    Parameters:
        ux, uy, uz : 3D velocity field components
        x, y, z : grid vectors
        time : simulation time
        LES : show LES cutoff lines
        normalize : normalize spectrum
        return_data : return dict of spectra
        re_init : initial Reynolds number (for eta calculation)

    Returns:
        dict of spectra if return_data is True
    """

    # Grid info
    nx, ny, nz = len(x), len(y), len(z)
    dx, dy, dz = x[1] - x[0], y[1] - y[0], z[1] - z[0]
    delta = (dx * dy * dz) ** (1 / 3)
    nu = 1.0 / re_init

    # FFT
    uxf = np.fft.fftn(ux)
    uyf = np.fft.fftn(uy)
    uzf = np.fft.fftn(uz)

    ex = 0.5 * np.abs(uxf)**2
    ey = 0.5 * np.abs(uyf)**2
    ez = 0.5 * np.abs(uzf)**2
    etot = ex + ey + ez

    # Frequencies
    kx = np.fft.fftfreq(nx, d=dx) * 2 * np.pi
    ky = np.fft.fftfreq(ny, d=dy) * 2 * np.pi
    kz = np.fft.fftfreq(nz, d=dz) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    k_mag = np.sqrt(KX**2 + KY**2 + KZ**2).flatten()

    ex, ey, ez, etot = ex.flatten(), ey.flatten(), ez.flatten(), etot.flatten()

    # Energy binning
    k_bins = np.linspace(0, k_mag.max(), 100)
    k_centers = 0.5 * (k_bins[:-1] + k_bins[1:])
    E_k = np.zeros_like(k_centers)
    E_x = np.zeros_like(k_centers)
    E_y = np.zeros_like(k_centers)
    E_z = np.zeros_like(k_centers)

    for i in range(len(k_bins) - 1):
        mask = (k_mag >= k_bins[i]) & (k_mag < k_bins[i + 1])
        if np.any(mask):
            E_k[i] = np.sum(etot[mask]) / np.sum(mask)
            E_x[i] = np.sum(ex[mask]) / np.sum(mask)
            E_y[i] = np.sum(ey[mask]) / np.sum(mask)
            E_z[i] = np.sum(ez[mask]) / np.sum(mask)

    if normalize:
        E_k /= np.trapz(E_k, k_centers)
        E_x /= np.trapz(E_x, k_centers)
        E_y /= np.trapz(E_y, k_centers)
        E_z /= np.trapz(E_z, k_centers)

    # Compute enstrophy-based Kolmogorov scale
    grad = np.gradient
    omega_x = grad(uz, y, axis=1) - grad(uy, z, axis=2)
    omega_y = grad(ux, z, axis=2) - grad(uz, x, axis=0)
    omega_z = grad(uy, x, axis=0) - grad(ux, y, axis=1)
    enstrophy = np.mean(omega_x**2 + omega_y**2 + omega_z**2)
    epsilon = 2 * nu * enstrophy
    eta = (nu**3 / epsilon)**0.25 if epsilon > 0 else np.inf
    k_eta = 1.0 / eta if eta < np.inf else 0.0

    # Nyquist frequency
    k_nyquist = np.pi / min(dx, dy, dz)

    # Resolution ratio
    delta_over_eta = delta / eta if eta < np.inf else np.inf

    # Formatted terminal output using print_block
    res_quality = (
        "OK (DNS)" if delta_over_eta <= 1.5 else
        "Marginal" if delta_over_eta <= 3.0 else
        "Under-resolved"
    )

    # Find peak energy position (excluding k ≈ 0)
    valid = (k_centers > 1e-8)  # exclude k=0 or near-zero
    k_peak = k_centers[valid][np.argmax(E_k[valid])]

    # Compute local slope in log-log space
    log_k = np.log10(k_centers)
    log_E = np.log10(E_k)
    slope_local = np.gradient(log_E, log_k)

    # Define inertial range as slope ≈ -5/3 ± 0.3
    mask_inertial = (slope_local < -1.3) & (slope_local > -2.0)

    print_block("Spectral resolution diagnostics", [
        ("nu (viscosity)",        f"{nu:.2e}"),
        ("Mean enstrophy",        f"{enstrophy:.4e}"),
        ("Dissipation epsilon",   f"{epsilon:.4e}"),
        ("Kolmogorov scale eta",  f"{eta:.4e}"),
        ("Delta (grid scale)",    f"{delta:.4f}"),
        ("Delta / eta",           f"{delta_over_eta:.2f}"),
        ("k_max (Nyquist)",       f"{k_nyquist:.2f}"),
        ("k_eta = 1 / eta",       f"{k_eta:.2f}"),
        ("Resolution check",      res_quality),
        ("k_peak (eddy scale)"  , f"{k_peak:.2f}")
    ])
    if np.any(mask_inertial):
        k_inertial = k_centers[mask_inertial]
        inertial_width = np.log10(k_inertial[-1]) - np.log10(k_inertial[0])
        print(f"  Inertial range width log10(k): {inertial_width:.2f}")
        print(f"  From k = {k_inertial[0]:.2f} to {k_inertial[-1]:.2f}")
    else:
        print("   No clear inertial range detected.")

    print("-" * 48)

    # Plot
    plt.figure(figsize=(16, 9))
    plt.loglog(k_centers, E_k, '-', color='black', lw=2, label=r"$E(k)$")
    plt.loglog(k_centers, E_x, '-', color='tab:red', label=r"$E_x(k)$")
    plt.loglog(k_centers, E_y, '--', color='tab:green', label=r"$E_y(k)$")
    plt.loglog(k_centers, E_z, '-.', color='tab:blue', label=r"$E_z(k)$")

    # Reference slope
    i1 = int(0.05 * len(k_centers))
    i2 = int(0.40 * len(k_centers))
    slope, _ = np.polyfit(np.log10(k_centers[i1:i2]), np.log10(E_k[i1:i2]), 1)
    ref_k = k_centers[i1:i2]
    ref_E = ref_k**(-5/3) * E_k[i1] / (ref_k[0]**(-5/3))
    plt.loglog(ref_k, ref_E, 'k--', label=r"$k^{-5/3}$")
    plt.text(ref_k[0], ref_E[0]*1.2, f"slope ≈ {slope:.2f}", fontsize=10, color='gray')

    # k_max and k_eta
    plt.axvline(k_nyquist, color='gray', linestyle=':', label=r"$k_{\max}$ (Nyquist)")
    if eta < np.inf:
        plt.axvline(k_eta, color='orange', linestyle='--', label=r"$k_\eta = 1/\eta$")

    if LES:
        k_cutoff = np.pi / delta
        plt.axvline(k_cutoff, color='red', linestyle='--', label=fr"$k_c = \pi/\Delta$")

    plt.title(f"Energy spectrum at t = {time:.2f}", fontsize=14)
    plt.xlabel(r"Wavenumber $k$", fontsize=12)
    plt.ylabel(r"Spectral energy $E(k)$", fontsize=12)
    plt.grid(True, which='both', linestyle=':', alpha=0.6)
    plt.legend(loc='lower left')
    plt.tight_layout()
    if save_fig:
        if not os.path.exists("png"):
            os.makedirs("png")
        plt.savefig(f"png/{fig_name}", dpi=300, bbox_inches='tight', transparent=False)
    else:
        plt.show()   

    if return_data:
        return {
            'k': k_centers,
            'E_k': E_k,
            'E_x': E_x,
            'E_y': E_y,
            'E_z': E_z,
            'eta': eta,
            'delta_over_eta': delta_over_eta
        }

def compute_smagorinsky_viscosity(x, y, z, ux, uy, uz, cs=0.17):
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
    dux_dx = np.gradient(ux, x, axis=0, edge_order=2)
    dux_dy = np.gradient(ux, y, axis=1, edge_order=2)
    dux_dz = np.gradient(ux, z, axis=2, edge_order=2)
    duy_dx = np.gradient(uy, x, axis=0, edge_order=2)
    duy_dy = np.gradient(uy, y, axis=1, edge_order=2)
    duy_dz = np.gradient(uy, z, axis=2, edge_order=2)
    duz_dx = np.gradient(uz, x, axis=0, edge_order=2)
    duz_dy = np.gradient(uz, y, axis=1, edge_order=2)
    duz_dz = np.gradient(uz, z, axis=2, edge_order=2)

    # Compute strain-rate magnitude |S|
    S2 = (
        2 * dux_dx**2 +
        2 * duy_dy**2 +
        2 * duz_dz**2 +
        (dux_dy + duy_dx)**2 +
        (dux_dz + duz_dx)**2 +
        (duy_dz + duz_dy)**2
    )
    nut = (cs * delta) ** 2 * np.sqrt(S2)
    return nut, delta

def plot_nut_slices_combined(nut, x, y, z, time, x_plot=None, y_plot=None, z_plot=None, vmin=None, vmax=None, 
        save_fig=None, fig_name=None):
    """
    Affiche les trois coupes (x-y, x-z, y-z) de log10(nu_t) dans une seule fenêtre matplotlib.
    """

    # Trouve les indices proches des coordonnées demandées
    i_plot = np.argmin(np.abs(x - x_plot)) if x_plot is not None else len(x) // 2
    j_plot = np.argmin(np.abs(y - y_plot)) if y_plot is not None else len(y) // 2
    k_plot = np.argmin(np.abs(z - z_plot)) if z_plot is not None else len(z) // 2

    # Évite les valeurs nulles ou négatives avant log10
    nut_safe = np.where(nut > 1e-20, nut, 1e-20)
    lognut = np.log10(nut_safe)

    # Détermine les niveaux de couleur
    vmin = vmin if vmin is not None else np.min(lognut)
    vmax = vmax if vmax is not None else np.max(lognut)

    fig, axs = plt.subplots(1, 3, figsize=(16, 9), constrained_layout=True)
    cmap = "viridis"

    # x-y @ z = z_plot
    im0 = axs[0].contourf(x, y, lognut[:, :, k_plot].T, levels=50, cmap=cmap, vmin=vmin, vmax=vmax)
    axs[0].set_title(rf"$\log(\nu_t)$ in x-y at $z =${z[k_plot]:.2f}")
    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")

    # x-z @ y = y_plot
    im1 = axs[1].contourf(z, x, lognut[:, j_plot, :], levels=50, cmap=cmap, vmin=vmin, vmax=vmax)
    axs[1].set_title(rf"$\log(\nu_t)$ in x-z at $y =${y[j_plot]:.2f}")
    axs[1].set_xlabel("z")
    axs[1].set_ylabel("x")

    # y-z @ x = x_plot
    im2 = axs[2].contourf(z, y, lognut[i_plot, :, :], levels=50, cmap=cmap, vmin=vmin, vmax=vmax)
    axs[2].set_title(rf"$\log(\nu_t)$ in y-z at $x =${x[i_plot]:.2f}")
    axs[2].set_xlabel("z")
    axs[2].set_ylabel("y")

    # Barre de couleur partagée
    fig.colorbar(im2, ax=axs.ravel().tolist(), shrink=0.9, label=r"$\log(\nu_t)$")

    plt.suptitle(rf"Subgrid-scale viscosity slices $\log(\nu_t)$ at $t =${time:.1f}", fontsize=14)
    if save_fig:
        if not os.path.exists("png"):
            os.makedirs("png")
        plt.savefig(f"png/{fig_name}", dpi=300, bbox_inches='tight', transparent=False)
    else:
        plt.show()   

def analyze_turbulent_development(x, y, z, ux, uy, uz, time, re_init=1000.0, cs=0.17, LES=False):
    u_mag_sq = ux**2 + uy**2 + uz**2
    tke = 0.5 * np.mean(u_mag_sq)
    enstrophy = np.mean(
        (np.gradient(uz, axis=1) - np.gradient(uy, axis=2))**2 +
        (np.gradient(ux, axis=2) - np.gradient(uz, axis=0))**2 +
        (np.gradient(uy, axis=0) - np.gradient(ux, axis=1))**2
    )

    block = [
        ("Mean TKE", f"{tke:.4e}"),
        ("Mean enstrophy", f"{enstrophy:.4e}"),
    ]

    TKE_THRESHOLD = 1e-3
    ENSTROPHY_THRESHOLD = 1e-1

    if LES:
        nut, delta = compute_smagorinsky_viscosity(x, y, z, ux, uy, uz, cs=cs)
        nut_min = np.min(nut)
        nut_max = np.max(nut)
        nut_mean = np.mean(nut)
        lambda_cutoff = 2 * delta
        nu = 1.0 / re_init
        resolved_ratio = nut_max / (nut_max + nu)
        eta = 1.0 / (re_init ** 0.75)
        delta_over_eta = delta / eta

        if resolved_ratio < 0.2:
            if delta_over_eta <= 2.0:
                quality = "DNS-like (well resolved)"
            else:
                quality = "Under-resolved (model inactive, resolution too coarse)"
        elif resolved_ratio < 0.6:
            quality = "Good LES"
        else:
            quality = "Under-resolved LES"

        block += [
            ("Smagorinsky constant Cs", f"{cs:.2f}"),
            ("Smagorinsky nu_t min", f"{nut_min:.3e}"),
            ("Smagorinsky nu_t max", f"{nut_max:.3e}"),
            ("Smagorinsky nu_t mean", f"{nut_mean:.3e}"),
            ("LES cutoff length scale", f"{lambda_cutoff:.4f}"),
            ("nu", f"{nu:.3e}"),
            ("nu_t max", f"{nut_max:.3e}"),
            ("nu_t mean", f"{nut_mean:.3e}"),
            ("nu_t / (nu + nu_t)", f"{resolved_ratio:.2f}"),
            ("Filter size delta", f"{delta:.4e}"),
            ("Kolmogorov scale (eta)", f"{eta:.4e}"),
            ("Delta / eta (resolution)", f"{delta_over_eta:.2f}")
        ]

        print_block("Turbulence development analysis", block)
        print(f"\n  LES quality assessment     : {quality}")
        print(f"  Filter resolution check    : delta / eta = {delta_over_eta:.2f}")
        if delta_over_eta <= 1.0:
            print("  LES nearly DNS-resolved (resolves Kolmogorov scale)")
        elif delta_over_eta <= 2.0:
            print("  LES resolves inertial range adequately")
        else:
            print("  LES under-resolved (filter too coarse)")
        print("-" * 48)
    else:
        print_block("Turbulence development analysis", block)
        print("\n  Regime assessment         : ", end="")
        if tke > TKE_THRESHOLD and enstrophy > ENSTROPHY_THRESHOLD:
            print("TURBULENT")
        elif tke > 0.1 * TKE_THRESHOLD:
            print("TRANSITIONAL")
        else:
            print("LAMINAR or UNDERDEVELOPED")
        print("-" * 48)
        nut, delta=None, None

    return nut, delta

def plot_rms_profiles(x, y, z, ux, uy, uz, time, save_fig=None, fig_name=None):
    """
    Compute and plot the RMS profiles of velocity fluctuations.
    Fluctuations are defined as u' = u - <u>, where <u> is the mean over (x, z).
    """

    # Compute mean velocity profiles along y (averaged over x and z)
    ux_mean = np.mean(ux, axis=(0, 2))  # shape (ny,)
    uy_mean = np.mean(uy, axis=(0, 2))
    uz_mean = np.mean(uz, axis=(0, 2))

    # Compute fluctuations: u' = u - <u>
    ux_fluct = ux - ux_mean[np.newaxis, :, np.newaxis]
    uy_fluct = uy - uy_mean[np.newaxis, :, np.newaxis]
    uz_fluct = uz - uz_mean[np.newaxis, :, np.newaxis]

    # Compute RMS of fluctuations: sqrt(mean(u'^2)) over x and z
    ux_rms = np.sqrt(np.mean(ux_fluct**2, axis=(0, 2)))
    uy_rms = np.sqrt(np.mean(uy_fluct**2, axis=(0, 2)))
    uz_rms = np.sqrt(np.mean(uz_fluct**2, axis=(0, 2)))

    # Local turbulent kinetic energy (TKE) per y-slice
    tke_profile = 0.5 * (ux_rms**2 + uy_rms**2 + uz_rms**2)

    # Plotting RMS profiles
    plt.figure(figsize=(16, 9))
    plt.plot(ux_rms, y, label=r"$u'_x$ RMS", color='tab:red')
    plt.plot(uy_rms, y, label=r"$u'_y$ RMS", color='tab:green')
    plt.plot(uz_rms, y, label=r"$u'_z$ RMS", color='tab:blue')

    # Plot turbulent intensity: sqrt(2k)
    plt.plot(np.sqrt(2 * tke_profile), y,
             label=r"Turbulent intensity $\sqrt{2\,k}$",
             color='black', linestyle='--')

    print_block("RMS turbulence summary", [
        ("Max RMS u'_x", f"{np.max(ux_rms):.4f}"),
        ("Max RMS u'_y", f"{np.max(uy_rms):.4f}"),
        ("Max RMS u'_z", f"{np.max(uz_rms):.4f}"),
        ("Peak y (u'_x RMS)", f"{y[np.argmax(ux_rms)]:.4f}"),
        ("Peak y (u'_y RMS)", f"{y[np.argmax(uy_rms)]:.4f}"),
        ("Peak y (u'_z RMS)", f"{y[np.argmax(uz_rms)]:.4f}"),
        ("Max sqrt(2k)", f"{np.max(np.sqrt(2 * tke_profile)):.4f}")
    ])

    anisotropy_yz = np.max(uy_rms) / np.max(uz_rms)
    anisotropy_xz = np.max(ux_rms) / np.max(uz_rms)

    print_block("Anisotropy indicators", [
        ("u'_y / u'_z", f"{anisotropy_yz:.2f}"),
        ("u'_x / u'_z", f"{anisotropy_xz:.2f}")
    ])

    plt.xlabel("RMS velocity fluctuations", fontsize=12)
    plt.ylabel("y", fontsize=12)
    plt.title(rf"Turbulent RMS profiles at $t = {time:.2f}$", fontsize=14)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    if save_fig:
        if not os.path.exists("png"):
            os.makedirs("png")
        plt.savefig(f"png/{fig_name}", dpi=300, bbox_inches='tight', transparent=False)
    else:
        plt.show()

    return ux_rms, uy_rms, uz_rms, tke_profile

def plot_time_evolution(data_log, save_fig=None, fig_name=None):
    """
    Plot temporal evolution of mixing layer and turbulence indicators.
    """
    times = np.array([d['time'] for d in data_log])
    Re_delta = np.array([d['Re_delta'] for d in data_log])
    Re_theta = np.array([d['Re_theta'] for d in data_log])
    sqrt2k = np.array([d['sqrt2k'] for d in data_log])
    rms_ux = np.array([d['rms_ux'] for d in data_log])
    rms_uy = np.array([d['rms_uy'] for d in data_log])
    rms_uz = np.array([d['rms_uz'] for d in data_log])

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))

    axs[0, 0].plot(times, Re_delta, 'o-', label=r"$Re_\delta$")
    axs[0, 0].plot(times, Re_theta, 's-', label=r"$Re_\theta$")
    axs[0, 0].set_title("Reynolds numbers")
    axs[0, 0].set_xlabel("Time")
    axs[0, 0].legend()
    axs[0, 0].grid()

    axs[0, 1].plot(times, sqrt2k, 'k--', label=r"$\sqrt{2k}$")
    axs[0, 1].set_title("Max turbulent intensity")
    axs[0, 1].set_xlabel("Time")
    axs[0, 1].grid()

    axs[1, 0].plot(times, rms_ux, 'r-', label=r"$u'_x$")
    axs[1, 0].plot(times, rms_uy, 'g-', label=r"$u'_y$")
    axs[1, 0].plot(times, rms_uz, 'b-', label=r"$u'_z$")
    axs[1, 0].set_title("Max RMS fluctuations")
    axs[1, 0].set_xlabel("Time")
    axs[1, 0].legend()
    axs[1, 0].grid()

    axs[1, 1].axis('off')
    # Extract columns
    times = [d['time'] for d in data_log]
    Re_delta = [d['Re_delta'] for d in data_log]
    sqrt2k = [d['sqrt2k'] for d in data_log]
    rms_ux = [d['rms_ux'] for d in data_log]
    rms_uy = [d['rms_uy'] for d in data_log]
    rms_uz = [d['rms_uz'] for d in data_log]
    table_data = [["t", "Re_δ", "√2k", "RMS u'_x", "u'_y", "u'_z"]]
    for d in data_log:
        table_data.append([
            f"{d['time']:.1f}",
            f"{d['Re_delta']:.0f}",
            f"{d['sqrt2k']:.4f}",
            f"{d['rms_ux']:.4f}",
            f"{d['rms_uy']:.4f}",
            f"{d['rms_uz']:.4f}"
        ])
    
    table = axs[1, 1].table(cellText=table_data,
                            cellLoc='center',
                            loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    fig.suptitle("Time evolution of turbulence indicators", fontsize=16)
    plt.tight_layout()
    if save_fig:
        if not os.path.exists("png"):
            os.makedirs("png")
        plt.savefig(f"png/{fig_name}", dpi=300, bbox_inches='tight', transparent=False)
    else:
        plt.show()

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

    parser.add_argument('--reynolds', '-re', type=float, default=1000.0,
            help="Reynolds number initial")

    parser.add_argument('--smagorinsky_constant', '-cs', type=float, default=0.17,
            help="Smagorinsky constant (default=0.17)")

    parser.add_argument('--mixing_layer', '-ml', action='store_true',
            help="Activate mixing layer mode (3-panel plot in the y-direction)")

    parser.add_argument('--planar_jet', '-pj', action='store_true',
            help="Analyse du jet plan : profils + épaisseur + theta. Fournir U_j (vitesse du jet)")

    parser.add_argument('--fluctuations', '-f', action='store_true',
            help="Display velocity fluctuations u'_i in the (y,z) plane at given x")

    parser.add_argument('-rms', action='store_true', help='Plot RMS velocity fluctuations')

    parser.add_argument('--turbulence', '-t', action='store_true',
            help="Display turbulence statistics and plot the energy spectrum")
    parser.add_argument('-les', action='store_true',
            help="Enable les mode to show expected spectral cutoff due to Smagorinsky model")
    parser.add_argument('-d', '--divergence', action='store_true',
                        help="Plot divergence div(u) in planes yz, xz, and xy at given positions")

    parser.add_argument('--hit', action='store_true',
                        help="Compute Re_lambda for Homogeneous Isotrpic Turbulence")

    parser.add_argument('--evolution', '-e', action='store_true',
            help="Plot time evolution of key turbulence indicators from all files")

    parser.add_argument('--pdf', action='store_true',
        help="Save all figures to ./png and compile them into a single PDF at the end.")

    args = parser.parse_args()
 
    file_list = list(itertools.chain.from_iterable(glob.glob(f) for f in args.filename))
    file_list = sorted(set(file_list))  # dédoublonner + trier

    field_time_map = {}

    if not file_list:
        print("Aucun fichier trouvé pour les motifs donnés.")
        return 1

    if args.pdf:
        import shutil

        os.makedirs("png", exist_ok=True)
        # Nettoyer ./png
        for f in os.listdir("png"):
            os.remove(os.path.join("png", f))

        os.makedirs("pdf", exist_ok=True)

    data_log = []  # store time series data
    for fname in file_list:
        print(f"\n========= Analyse de {fname} =========")
        basename = os.path.basename(fname)
        file_id = os.path.splitext(basename)[0].split('_')[-1]  # extract XXXXXX
        time, x, y, z, ux, uy, uz, pp = read_fields(fname)
        field_time_map[file_id] = time
        block = [
            ("File", fname),
            ("Time", f"{time:.1f}"),
            ("nx, ny, nz", f"{len(x)}, {len(y)}, {len(z)}"),
            ("Lx, Ly, Lz", f"{x[-1]-x[0]:.1f}, {y[-1]-y[0]:.1f}, {z[-1]-z[0]:.1f}"),
            ("ux min/max", f"{ux.min():.3e} / {ux.max():.3e} "),
            ("uy min/max", f"{uy.min():.3e} / {uy.max():.3e} "),
            ("uz min/max", f"{uz.min():.3e} / {uz.max():.3e} "),
            ("Pressure min", f"{pp.min():.3e}"),
            ("Pressure max", f"{pp.max():.3e}")
        ]
        print_block("Field summary", block)

        if args.plot_velocity_profiles:
            plot_velocity_profiles(time, x, y, z, ux, uy, uz,
                args.x, args.y, args.z, args.direction,
                save_fig=args.pdf,
                fig_name=f"{file_id}_velocity_profile_{args.direction}.png")
        
        if args.mixing_layer:
            delta_value, theta_value, re_delta_value, re_theta_value = plot_mixing_layer_profiles(
                x, y, z, ux, uy, uz,
                args.x, args.z, time, re_init=args.reynolds,
                save_fig=args.pdf,
                fig_name=f"{file_id}_mixing_layer.png")
        
        if args.planar_jet:
            delta_value, theta_value, re_delta_value, re_theta_value = plot_planar_jet_profiles(
                x, y, z, ux, uy, uz,
                args.x, args.z, time, u_jet=1.,
                re_init=args.reynolds,
                save_fig=args.pdf,
                fig_name=f"{file_id}_planar_jet.png")
        
        if args.fluctuations:
            plot_velocity_fluctuations_plan(
                x, y, z, ux, uy, uz,
                args.x, args.y, args.z,
                time, direction=args.direction,
                save_fig=args.pdf,
                fig_name=f"{file_id}_fluctuations_{args.direction}.png")
        
        if args.turbulence:
            plot_energy_spectrum(ux, uy, uz, x, y, z, time,
                LES=args.les, normalize=True,
                save_fig=args.pdf, re_init=args.reynolds, 
                fig_name=f"{file_id}_kolmogorov_spectrum.png")
            nut, delta_filter = analyze_turbulent_development(x, y, z, ux, uy, uz, time,
                re_init=args.reynolds, cs=args.smagorinsky_constant,
                LES=args.les)
            if args.les:
                plot_nut_slices_combined(nut, x, y, z, time, x_plot=args.x, y_plot=args.y, z_plot=args.z,
                        save_fig=args.pdf, fig_name=f"{file_id}_nut_planes.png")
        
        if args.rms:
            ux_rms, uy_rms, uz_rms, tke_profile = plot_rms_profiles(
                x, y, z, ux, uy, uz, time,
                save_fig=args.pdf,
                fig_name=f"{file_id}_rms_profiles.png")
        
        if args.divergence:
            plot_divergence_planes(x, y, z, ux, uy, uz,
                args.x, args.y, args.z, time,
                save_fig=args.pdf,
                fig_name=f"{file_id}_divergence_planes.png")
        if args.hit:
            compute_Re_lambda(ux, uy, uz, x, y, z, 1./args.reynolds)
        

        # Collect data for time evolution plots
        if args.evolution:
            rms_vals = {
                    'time': time,
                    'delta': delta_value,  # à récupérer dans plot_mixing_layer_profiles
                    'theta': theta_value,  # idem
                    'Re_delta': re_delta_value,  # à récupérer dans plot_mixing_layer_profiles
                    'Re_theta': re_theta_value,  # idem
                    'rms_ux': np.max(ux_rms),    # depuis plot_rms_profiles
                    'rms_uy': np.max(uy_rms),
                    'rms_uz': np.max(uz_rms),
                    'sqrt2k': np.max(np.sqrt(2 * tke_profile)),
                    }
            data_log.append(rms_vals)

    if args.evolution:
        plot_time_evolution(data_log,
                save_fig=args.pdf,
                fig_name=f"evolution.png")

    if args.pdf:
        generate_latex_figure_pdf(field_time_map=field_time_map)
    
    return 0

if __name__ == "__main__":

    err = main()
    sys.exit(0)
