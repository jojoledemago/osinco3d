#!/usr/bin/env python3
"""
case_hit.py -- decaying homogeneous isotropic turbulence (fully periodic cube).

Same diagnostics as before; only the data access changed: the velocity
gradients are now produced ON DEMAND, one at a time, and reduced immediately
(means of squares, divergence, skewness...), so the nine components of the
gradient tensor are never stored simultaneously.  The numerical operators
(FD6 and spectral) are unchanged.
"""

import numpy as np
import matplotlib.pyplot as plt

from fields_core import (Case, VELOCITY, print_block,
                         fd6_first_derivative, spectral_first_derivative)

# ----------------------------------------------------------------------------
# Spectrum
# ----------------------------------------------------------------------------
def energy_spectrum(ux, uy, uz, L):
    """
    Shell-averaged energy spectrum on a cubic periodic box.

    Returns
    -------
    k_shell : wavenumbers (2*pi/L * integer shell index), shell 0 excluded
    E_shell : energy contained in each shell  (sum = kinetic energy, minus
              the corner modes beyond n//2 and the mean flow)
    e3d     : 3D spectral energy density (for exact Parseval sums)
    k2      : 3D array of |k|^2 (physical units)
    """
    n = ux.shape[0]
    n_tot = ux.size
    e3d = 0.5 * sum(np.abs(np.fft.fftn(u) / n_tot) ** 2 for u in (ux, uy, uz))

    k1 = np.fft.fftfreq(n, d=1.0 / n)               # integer wavenumbers
    kx, ky, kz = np.meshgrid(k1, k1, k1, indexing="ij")
    kmag_int = np.sqrt(kx**2 + ky**2 + kz**2)
    shells = np.rint(kmag_int).astype(int)

    kmax = n // 2
    E_shell = np.bincount(shells.ravel(), weights=e3d.ravel(),
                          minlength=kmax + 1)[: kmax + 1]
    dk = 2.0 * np.pi / L
    k_shell = np.arange(kmax + 1) * dk
    k2 = (kmag_int * dk) ** 2
    return k_shell[1:], E_shell[1:], e3d, k2


# ----------------------------------------------------------------------------
# HIT analysis
# ----------------------------------------------------------------------------
def gradient_diagnostics(grad, moments=True):
    """
    Reduce the velocity-gradient tensor without ever storing all of it.

    `grad(i, j)` must return d u_i / d x_j (computed on demand).

    Returns enstrophy 1/2<w^2>, <S_ij S_ij>, divergence (max and rms) and,
    if `moments`, the mean skewness / flatness of the diagonal components.
    """
    div = None
    sij2 = 0.0
    skew, flat = [], []

    for i in range(3):                       # diagonal components
        g = grad(i, i)
        div = g.copy() if div is None else div + g
        m2 = np.mean(g**2)
        sij2 += m2                           # S_ii S_ii = (du_i/dx_i)^2
        if moments:
            skew.append(np.mean(g**3) / m2**1.5)
            flat.append(np.mean(g**4) / m2**2)

    w2 = 0.0
    for i in range(3):                       # off-diagonal pairs (i < j)
        for j in range(i + 1, 3):
            gij, gji = grad(i, j), grad(j, i)
            w2 += np.mean((gji - gij) ** 2)              # one vorticity comp.
            sij2 += 2.0 * np.mean((0.5 * (gij + gji)) ** 2)   # S_ij, S_ji

    out = {
        "enstrophy": 0.5 * w2,
        "sij2": sij2,
        "div_max": np.abs(div).max(),
        "div_rms": np.sqrt(np.mean(div**2)),
    }
    if moments:
        out["skew"] = np.mean(skew)
        out["flat"] = np.mean(flat)
    return out


def hit_analysis(ff, nu):
    """Compute all HIT diagnostics for one snapshot (FieldFile `ff`)."""
    g = ff.grid
    if any(b != "periodic" for b in g.bc):
        raise ValueError("HIT analysis requires periodic BC in x, y and z.")
    nx, ny, nz = g.n
    if not (nx == ny == nz):
        raise ValueError("HIT analysis requires nx = ny = nz.")
    dx, dy, dz = g.d
    if not (np.isclose(dx, dy) and np.isclose(dx, dz)):
        raise ValueError("HIT analysis requires dx = dy = dz.")
    n, L = nx, nx * dx

    vel = tuple(ff.load(name) for name in VELOCITY)
    ux, uy, uz = vel

    # --- gradients: FD6 (same operator as the solver), on demand ------------
    gd_fd = gradient_diagnostics(
        lambda i, j: fd6_first_derivative(vel[i], dx, j))

    # --- same quantities, spectral (reference) ------------------------------
    gd_sp = gradient_diagnostics(
        lambda i, j: spectral_first_derivative(vel[i], L, j), moments=False)

    # --- global quantities ---------------------------------------------------
    mean_u = [u.mean() for u in vel]
    ke = 0.5 * sum(np.mean(u**2) for u in vel)
    urms_comp = [np.sqrt(np.mean((u - m) ** 2)) for u, m in zip(vel, mean_u)]
    urms = np.sqrt(2.0 * ke / 3.0)

    enstrophy_fd = gd_fd["enstrophy"]
    enstrophy_sp = gd_sp["enstrophy"]
    eps_fd = 2.0 * nu * enstrophy_fd           # homogeneous: eps = nu <w^2>
    eps_sp = 2.0 * nu * enstrophy_sp

    # strain-based dissipation: eps = 2 nu <S_ij S_ij> (FD6)
    eps_strain = 2.0 * nu * gd_fd["sij2"]

    # --- spectrum ---------------------------------------------------------------
    k, E_shell, e3d, k2 = energy_spectrum(ux, uy, uz, L)
    dk = k[0]
    Ek = E_shell / dk                                   # E(k) density
    eps_spec = 2.0 * nu * np.sum(k2 * e3d)              # exact Parseval sum
    ke_spec = np.sum(e3d)
    captured = E_shell.sum() / ke_spec

    # --- scales -------------------------------------------------------------------
    eps_ref = eps_sp                                     # spectral = reference
    eta = (nu**3 / eps_ref) ** 0.25
    tau_eta = np.sqrt(nu / eps_ref)
    lam = np.sqrt(15.0 * nu * urms**2 / eps_ref)         # Taylor microscale
    re_lambda = urms * lam / nu
    l_int = (3.0 * np.pi / 4.0) * np.sum(E_shell / k) / np.sum(E_shell)
    tau_e = l_int / urms                                 # eddy turnover time
    kmax_brut = np.pi / dx
    kmax_23 = (2.0 / 3.0) * kmax_brut

    return {
        "iteration": ff.iteration, "time": ff.time,
        "n": n, "L": L, "dx": dx,
        "mean_u": mean_u, "urms_comp": urms_comp,
        "E": ke, "E_spec": ke_spec, "captured": captured,
        "enstrophy": enstrophy_fd,
        "eps_fd": eps_fd, "eps_sp": eps_sp,
        "eps_strain": eps_strain, "eps_spec": eps_spec,
        "div_fd_max": gd_fd["div_max"],
        "div_fd_rms": gd_fd["div_rms"],
        "div_sp_max": gd_sp["div_max"],
        "urms": urms, "eta": eta, "tau_eta": tau_eta,
        "lambda": lam, "Re_lambda": re_lambda,
        "L_int": l_int, "tau_e": tau_e,
        "L_over_Lint": L / l_int,
        "kmax_eta": kmax_brut * eta, "kmax23_eta": kmax_23 * eta,
        "dx_over_eta": dx / eta,
        "skew": gd_fd["skew"], "flat": gd_fd["flat"],
        "k": k, "Ek": Ek,
    }

def print_hit_block(s):
    print_block(
        f"HIT diagnostics  it={s['iteration']:06d}  t={s['time']:.4g}",
        [
            ("Kinetic energy E", s["E"]),
            ("E (spectrum, Parseval)", s["E_spec"]),
            ("Energy captured by E(k)", s["captured"]),
            ("<u>, <v>, <w> (should be ~0)",
             "  ".join(f"{m:.2e}" for m in s["mean_u"])),
            ("u'x, u'y, u'z (isotropy)",
             "  ".join(f"{u:.4g}" for u in s["urms_comp"])),
            ("Enstrophy (FD6)", s["enstrophy"]),
            ("eps = 2 nu <w^2>/2 (FD6)", s["eps_fd"]),
            ("eps = 2 nu <S:S>   (FD6)", s["eps_strain"]),
            ("eps from vorticity (spectral)", s["eps_sp"]),
            ("eps from spectrum (Parseval)", s["eps_spec"]),
            ("max|div u| (FD6)", f"{s['div_fd_max']:.3e}"),
            ("rms|div u| (FD6)", f"{s['div_fd_rms']:.3e}"),
            ("max|div u| (spectral)", f"{s['div_sp_max']:.3e}"),
            ("u_rms", s["urms"]),
            ("Kolmogorov eta", s["eta"]),
            ("Taylor lambda", s["lambda"]),
            ("Integral scale L_int", s["L_int"]),
            ("Re_lambda", s["Re_lambda"]),
            ("Eddy turnover L_int/u'", s["tau_e"]),
            ("L_box / L_int", s["L_over_Lint"]),
            ("dx / eta", s["dx_over_eta"]),
            ("kmax*eta (kmax=pi/dx)", s["kmax_eta"]),
            ("kmax*eta (kmax=2/3 pi/dx)", s["kmax23_eta"]),
            ("Skewness du_i/dx_i", s["skew"]),
            ("Flatness du_i/dx_i", s["flat"]),
        ],
    )

def hit_warnings(s):
    """Print simple pass/fail hints."""
    msgs = []
    if s["kmax_eta"] < 1.5:
        msgs.append(f"kmax*eta = {s['kmax_eta']:.2f} < 1.5 : DNS under-resolved.")
    if s["L_over_Lint"] < 6.0:
        msgs.append(f"L/L_int = {s['L_over_Lint']:.1f} < 6 : "
                    "periodicity may contaminate large scales.")
    rel = abs(s["eps_fd"] - s["eps_sp"]) / max(s["eps_sp"], 1e-300)
    if rel > 0.05:
        msgs.append(f"eps(FD6) differs from eps(spectral) by {100*rel:.1f} % : "
                    "numerical error in gradients/under-resolution.")
    if s["captured"] < 0.99:
        msgs.append(f"Only {100*s['captured']:.1f} % of energy in shells "
                    "k <= n/2 (rest is in the cube corners).")
    if abs(s["mean_u"][0]) + abs(s["mean_u"][1]) + abs(s["mean_u"][2]) > 1e-8:
        msgs.append("Non-zero mean velocity.")
    if msgs:
        print("  Warnings:")
        for m in msgs:
            print(f"   * {m}")

def hit_time_analysis(history):
    """Energy budget dE/dt = -eps over the series of snapshots."""
    t = np.array([h["time"] for h in history])
    E = np.array([h["E"] for h in history])
    eps = np.array([h["eps_sp"] for h in history])
    if len(t) < 3:
        return None
    dEdt = np.gradient(E, t)
    rel = np.abs(dEdt + eps) / np.maximum(eps, 1e-300)
    return dEdt, rel

def save_history(history, filename="hit_history.dat"):
    """Write time series for gnuplot."""
    cols = ["it", "t", "E", "eps_fd", "eps_sp", "eps_spec", "div_fd_max",
            "urms", "eta", "lambda", "Re_lambda", "L_int", "kmax_eta",
            "skew", "flat"]
    keys = ["iteration", "time", "E", "eps_fd", "eps_sp", "eps_spec",
            "div_fd_max", "urms", "eta", "lambda", "Re_lambda", "L_int",
            "kmax_eta", "skew", "flat"]
    with open(filename, "w") as f:
        f.write("# " + " ".join(f"{c:>14}" for c in cols) + "\n")
        for h in history:
            f.write("  " + " ".join(f"{h[k]:14.6e}" for k in keys) + "\n")
    print(f"\nTime series written to {filename}")

def plot_hit(history, show=True, outfile="hit_analysis.png"):
    t = np.array([h["time"] for h in history])
    get = lambda key: np.array([h[key] for h in history])

    fig, axs = plt.subplots(2, 3, figsize=(16, 9))

    # (0,0) energy and dissipation
    ax = axs[0, 0]
    ax.plot(t, get("E"), "o-", label="E")
    ax.set_xlabel("t")
    ax.set_ylabel("E")
    ax2 = ax.twinx()
    ax2.plot(t, get("eps_sp"), "s-", color="C1", label="eps (spectral)")
    ax2.plot(t, get("eps_fd"), "^--", color="C3", label="eps (FD6)")
    if len(t) >= 3:
        dEdt, _ = hit_time_analysis(history)
        ax2.plot(t, -dEdt, "k:", label="-dE/dt")
    ax2.set_ylabel("eps")
    ax.set_title("Energy and dissipation")
    ax2.legend(loc="upper right")

    # (0,1) spectra
    ax = axs[0, 1]
    idx = np.unique(np.linspace(0, len(history) - 1,
                                min(len(history), 8)).astype(int))
    colors = plt.cm.viridis(np.linspace(0, 1, len(idx)))
    for c, i in zip(colors, idx):
        h = history[i]
        ax.loglog(h["k"], h["Ek"], color=c, label=f"t={h['time']:.3g}")
    h = history[idx[0]]
    kref = h["k"][(h["k"] > h["k"][2]) & (h["k"] < h["k"][len(h["k"]) // 2])]
    if kref.size:
        ref = kref ** (-5.0 / 3.0)
        ref *= h["Ek"][2] / ref[0] * 0.5
        ax.loglog(kref, ref, "k--", label=r"$k^{-5/3}$")
    ax.set_ylim(bottom=max(1e-16, np.min([hh["Ek"][hh["Ek"] > 0].min()
                                          for hh in history])))
    ax.set_xlabel("k")
    ax.set_ylabel("E(k)")
    ax.set_title("Energy spectrum")
    ax.legend(fontsize=7)

    # (0,2) compensated spectrum (Kolmogorov)
    ax = axs[0, 2]
    for c, i in zip(colors, idx):
        h = history[i]
        comp = h["k"] ** (5.0 / 3.0) * h["Ek"] / h["eps_sp"] ** (2.0 / 3.0)
        ax.semilogx(h["k"], comp, color=c)
    ax.axhline(1.5, color="k", ls="--", lw=0.8, label="C_K ~ 1.5")
    ax.set_xlabel("k")
    ax.set_ylabel(r"$k^{5/3}E(k)/\epsilon^{2/3}$")
    ax.set_title("Compensated spectrum")
    ax.legend()

    # (1,0) divergence
    ax = axs[1, 0]
    ax.semilogy(t, get("div_fd_max"), "o-", label="max|div u| FD6")
    ax.semilogy(t, get("div_sp_max"), "s--", label="max|div u| spectral")
    ax.set_xlabel("t")
    ax.set_title("Divergence")
    ax.legend()

    # (1,1) Re_lambda and resolution
    ax = axs[1, 1]
    ax.plot(t, get("Re_lambda"), "o-", label=r"$Re_\lambda$")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$Re_\lambda$")
    ax2 = ax.twinx()
    ax2.plot(t, get("kmax_eta"), "s-", color="C1", label=r"$k_{max}\eta$")
    ax2.axhline(1.5, color="C1", ls=":", lw=0.8)
    ax2.set_ylabel(r"$k_{max}\eta$")
    ax.set_title("Reynolds and resolution")
    ax.legend(loc="upper left")
    ax2.legend(loc="upper right")

    # (1,2) gradient statistics
    ax = axs[1, 2]
    ax.plot(t, get("skew"), "o-", label="skewness")
    ax.axhline(-0.5, color="k", ls="--", lw=0.8)
    ax.set_xlabel("t")
    ax.set_ylabel("skewness")
    ax2 = ax.twinx()
    ax2.plot(t, get("flat"), "s-", color="C1", label="flatness")
    ax2.set_ylabel("flatness")
    ax.set_title("Velocity-gradient statistics")
    ax.legend(loc="upper left")
    ax2.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    print(f"Figure saved to {outfile}")


class HitCase(Case):
    """Decaying homogeneous isotropic turbulence, fully periodic cube."""
    name = "hit"
    bc = ("periodic", "periodic", "periodic")

    def analyse(self, ff, args):
        return hit_analysis(ff, args.nu)

    def report(self, s, args):
        if args.quiet:
            print(f"  t={s['time']:.4g} E={s['E']:.5g} "
                  f"eps={s['eps_sp']:.5g} Re_l={s['Re_lambda']:.4g} "
                  f"kmax*eta={s['kmax_eta']:.3g} "
                  f"div={s['div_fd_max']:.2e}")
        else:
            print_hit_block(s)
        hit_warnings(s)

    def finalize(self, history, args):
        res = hit_time_analysis(history)
        if res is not None:
            dEdt, rel = res
            print_block(
                "Energy budget  dE/dt = -eps",
                [("mean |dE/dt+eps|/eps", float(np.mean(rel))),
                 ("max  |dE/dt+eps|/eps", float(np.max(rel)))],
            )
        save_history(history, args.outdir / "hit_history.dat")
        if args.plot:
            plot_hit(history, show=True, outfile=args.outdir / "hit_analysis.png")
