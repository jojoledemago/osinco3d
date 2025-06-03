#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os

def read_stats(filepath):
    """
    Read CFD statistical data from the given path.
    Expected columns:
        0: time
        1: kinetic energy E_k
        2: temporal dissipation epsilon_t = -dE_k/dt
        3: viscous dissipation epsilon
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")
    data = np.loadtxt(filepath)
    return data

def plot_energy_dissipation(time, E_k, epsilon_t, epsilon):
    """
    Plot kinetic energy on the left y-axis, and dissipation terms on the right y-axis.
    Includes -dE_k/dt computed numerically for comparison with epsilon_t.
    """

    # Compute -d(E_k)/dt numerically
    dEk_dt_numeric = -np.gradient(E_k, time)

    fig, ax1 = plt.subplots(figsize=(16, 9))

    # Left Y-axis: E_k
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Kinetic Energy $E_k$", color='tab:blue')
    l1 = ax1.plot(time, E_k, color='tab:blue', label="Kinetic Energy $E_k$")
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax1.grid(True, linestyle='--', alpha=0.5)

    # Right Y-axis: dissipation
    ax2 = ax1.twinx()
    ax2.set_ylabel("Dissipation", color='tab:red')
    l2 = ax2.plot(time, epsilon_t, color='tab:red', linestyle='-', label=r"Temporal dissipation $\varepsilon_t$")
    l3 = ax2.plot(time, dEk_dt_numeric, color='gray', linestyle=':', label=r"Numerical $-dE_k/dt$")
    l4 = ax2.plot(time, epsilon, color='tab:green', linestyle='--', label=r"Viscous dissipation $\varepsilon$")
    ax2.tick_params(axis='y', labelcolor='tab:red')

    # Combine legends
    lines = l1 + l2 + l3 + l4
    labels = [line.get_label() for line in lines]
    ax1.legend(lines, labels, loc='upper right')

    plt.title("Kinetic energy and dissipation over time")
    fig.tight_layout()
    plt.show()

def assess_flow_regime(time, E_k, epsilon_t, epsilon, Re_init=None):
    """
    Print a qualitative assessment of the flow regime based on energy trends.
    """
    print("\n--- Flow Regime Assessment ---")
    print(f"Total time interval: [{time[0]:.2f}, {time[-1]:.2f}]")

    # Energy analysis
    dE = E_k[-1] - E_k[0]
    mean_eps_t = np.mean(epsilon_t)
    mean_eps = np.mean(epsilon)

    print(f"Kinetic energy change: delta E_k = {dE:+.4e}")
    print(f"Mean epsilon_t        : {mean_eps_t:.4e}")
    print(f"Mean viscous epsilon  : {mean_eps:.4e}")

    # Heuristic thresholds
    eps_threshold = 1e-3
    ek_threshold = 1e-4

    if mean_eps < eps_threshold and abs(dE) < ek_threshold:
        print("=> Likely a **laminar or steady** regime.")
    elif abs(dE) > ek_threshold and mean_eps < eps_threshold:
        print("=> Likely a **developing transient** flow.")
    elif mean_eps > eps_threshold and mean_eps_t < 1e-4:
        print("=> Likely a **turbulent but statistically stationary** regime.")
    elif mean_eps > eps_threshold and abs(mean_eps_t) > 1e-4:
        print("=> Likely a **fully turbulent** and **non-stationary** regime.")
    else:
        print("=> Regime unclear. Further diagnostics required.")

    if Re_init is not None:
        print(f"Initial Reynolds number: Re = {Re_init:.1f}")
    print("-------------------------------\n")

def main():

    # Inside main(), before reading file:
    parser = argparse.ArgumentParser(description="Analyze CFD global statistics over time.")
    parser.add_argument(
        "-re", "--reynolds",
        type=float,
        default=None,
        help="Initial Reynolds number (optional, used in diagnostic)"
    )
    args = parser.parse_args()

    stats_file = "outputs/stats.dat"
    try:
        data = read_stats(stats_file)
    except Exception as e:
        print(e)
        return

    # Extract relevant columns
    time = data[:, 0]
    E_k = data[:, 1]
    epsilon_t = data[:, 2]
    epsilon = data[:, 3]

    # Plot and assess
    plot_energy_dissipation(time, E_k, epsilon_t, epsilon)
    assess_flow_regime(time, E_k, epsilon_t, epsilon, Re_init=args.reynolds)

if __name__ == "__main__":
    main()

