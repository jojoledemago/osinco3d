#!/usr/bin/env python3
"""
les_sgs.py -- Smagorinsky sub-grid-scale (SGS) diagnostics, independent of the
flow case (hit, mixing_layer, tgv, ...).  Activated with --les [--cs X].

Model assumed (reconstructed from the resolved velocity, since nu_t is not
stored in the field files)
-----------------------------------------------------------------------------
    tau_ij - 1/3 tau_kk delta_ij = -2 nu_t S_ij ,    S_ij = (du_i/dx_j + du_j/dx_i)/2
    nu_t = (Cs Delta)^2 |S| ,                        |S| = sqrt(2 S_ij S_ij)
    Delta = (dx dy dz)^(1/3)    (see --les-delta)
    Volume averages use uniform weights on periodic directions and trapezoid
    weights on free-slip ones (boundary nodes carry half a cell).
    nu_eff = nu + nu_t ,  nu = 1/Re (non-dimensional equations)

The gradients are the 6th-order centred ones used everywhere else in this
post-processing.  If your solver defines Delta, |S| or the discrete operators
differently, the numbers below are a close reconstruction, not a bit-for-bit
copy of what the solver does; the budget check (below) tells you how well it
matches.

Quantities
----------
Statistics    <nu_t>/nu, max nu_t/nu, volume fraction where nu_t > nu and
              > 10 nu, PDF of nu_t/nu, plane profiles in y.
Dissipation   eps_nu  = 2 nu <S_ij S_ij>        (resolved, molecular)
              eps_sgs = <nu_t |S|^2> = (Cs Delta)^2 <|S|^3>  (SGS transfer)
              share   = eps_sgs / (eps_nu + eps_sgs)
Quality       Delta/eta with eta = (nu^3/eps_tot)^(1/4);
              LES_IQ_nu = 1/(1 + 0.05 (nu_eff/nu)^0.53)   (Celik et al. 2005;
              >= 0.8 : acceptable LES, >= 0.95 : close to DNS).
Budget        With no forcing and free-slip / periodic boundaries,
              dE/dt = -(eps_nu + eps_sgs)  with E = <u_i u_i>/2 (volume mean),
              exactly, for any case.  Over a series of snapshots this gives
                 Cs_inferred^2 = (-dE/dt - eps_nu) / (Delta^2 <|S|^3>)
              i.e. the Smagorinsky constant that the run actually behaves
              with, to compare with the Cs requested by --cs.  Numerical
              dissipation of the scheme/time integration, explicit filtering
              or forcing appear as a discrepancy.  Cs_inferred subtracts two
              close quantities when the SGS share is small and uses a time
              derivative of E: it needs closely spaced snapshots (checked on a
              1D Smagorinsky shear: error of a few % with dt_snap = 0.5 and an
              SGS share of 6 %, 0.2 % or less with finer spacing).
Reference     Lilly (1967): Cs = (1/pi) (3 C_K / 2)^(-3/4) = 0.173 for C_K = 1.5.
"""

import numpy as np
import matplotlib.pyplot as plt

from fields_core import VELOCITY, print_block

CS_LILLY = (1.0 / np.pi) * (3.0 * 1.5 / 2.0) ** -0.75
IQ_ALPHA, IQ_N = 0.05, 0.53                  # Celik et al. (2005), nu-based
HIST_LO, HIST_HI, HIST_N = -6.0, 4.0, 100    # log10(nu_t/nu) bins


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def filter_width(grid, spec):
    """Delta from 'cbrt' = (dx dy dz)^(1/3), 'max' = max(dx,dy,dz), or a number."""
    dx, dy, dz = grid.d
    if spec == "cbrt":
        return (dx * dy * dz) ** (1.0 / 3.0)
    if spec == "max":
        return max(dx, dy, dz)
    try:
        val = float(spec)
    except ValueError:
        raise ValueError(f"--les-delta must be 'cbrt', 'max' or a number, "
                         f"got '{spec}'.")
    if val <= 0.0:
        raise ValueError("--les-delta must be positive.")
    return val


def celik_iq_nu(nu_eff_over_nu):
    return 1.0 / (1.0 + IQ_ALPHA * np.asarray(nu_eff_over_nu) ** IQ_N)


# ----------------------------------------------------------------------------
# Snapshot analysis (one pass over y-slabs, gradients on demand)
# ----------------------------------------------------------------------------
def sgs_analysis(ff, args):
    g = ff.grid
    nx, ny, nz = g.n
    nu, cs = args.nu, args.cs
    delta = filter_width(g, args.les_delta)
    c2 = (cs * delta) ** 2
    wx, yw, wz = (g.weights(a) for a in range(3))
    pmean = lambda f: np.einsum("i,ijk,k->j", wx, f, wz)   # (x,z) average
    edges = np.linspace(HIST_LO, HIST_HI, HIST_N + 1)
    hist = np.zeros(HIST_N)

    S1 = np.zeros(ny)        # <|S|>
    S2m = np.zeros(ny)       # <S_ij S_ij>
    S3 = np.zeros(ny)        # <|S|^3>
    nut_max = np.zeros(ny)
    f1 = np.zeros(ny)        # fraction nu_t > nu
    f10 = np.zeros(ny)       # fraction nu_t > 10 nu
    ke = np.zeros(ny)        # <u_i u_i>/2

    for slab in ff.iter_slabs(args.slab):
        sl = slice(slab.j0, slab.j1)
        S2 = 0.0                                   # S_ij S_ij, on demand
        for i in range(3):
            gii = slab.vel_grad(i, i)
            S2 = S2 + gii * gii
            slab.clear()
        for i in range(3):
            for j in range(i + 1, 3):
                gs = slab.vel_grad(i, j) + slab.vel_grad(j, i)
                S2 = S2 + 0.5 * gs * gs            # 2 S_ij^2, S_ij = gs/2
                slab.clear()
        smag = np.sqrt(2.0 * S2)
        nut = c2 * smag
        ratio = nut / nu

        S1[sl] = pmean(smag)
        S2m[sl] = pmean(S2)
        S3[sl] = pmean(smag ** 3)
        nut_max[sl] = nut.max(axis=(0, 2))
        f1[sl] = pmean((ratio > 1.0).astype(float))
        f10[sl] = pmean((ratio > 10.0).astype(float))
        ke[sl] = 0.5 * sum(pmean(slab[name] ** 2) for name in VELOCITY)

        logr = np.clip(np.log10(np.maximum(ratio, 1e-300)),
                       HIST_LO, HIST_HI - 1e-9)
        w3 = (wx[:, None, None] * yw[sl][None, :, None]
              * wz[None, None, :])
        h, _ = np.histogram(logr.ravel(), bins=edges, weights=w3.ravel())
        hist += h

    vm = lambda f: float(np.sum(yw * f))            # volume mean of a profile
    nut_mean = c2 * S1                               # plane means
    eps_nu = 2.0 * nu * S2m
    eps_sgs = c2 * S3
    eps_tot = eps_nu + eps_sgs
    iq = celik_iq_nu((nu + nut_mean) / nu)
    with np.errstate(divide="ignore", invalid="ignore"):
        eta_y = np.where(eps_tot > 0.0,
                         (nu ** 3 / np.where(eps_tot > 0.0, eps_tot, 1.0))
                         ** 0.25, np.inf)

    nut_vol = vm(nut_mean)
    eps_nu_vol, eps_sgs_vol = vm(eps_nu), vm(eps_sgs)
    eps_tot_vol = eps_nu_vol + eps_sgs_vol
    eta = (nu ** 3 / eps_tot_vol) ** 0.25 if eps_tot_vol > 0 else float("inf")
    dlog = edges[1] - edges[0]

    return {
        "iteration": ff.iteration, "time": ff.time, "nu": nu, "cs": cs,
        "delta": delta, "delta_spec": str(args.les_delta),
        "delta_over_d": [delta / d for d in g.d],
        "y": g.coords[1], "yw": yw,
        # profiles
        "nut_mean_y": nut_mean / nu, "nut_max_y": nut_max / nu,
        "frac1_y": f1, "frac10_y": f10,
        "eps_nu_y": eps_nu, "eps_sgs_y": eps_sgs,
        "share_y": np.where(eps_tot > 0, eps_sgs / np.where(eps_tot > 0,
                            eps_tot, 1.0), 0.0),
        "iq_y": iq, "eta_y": eta_y,
        # volume quantities
        "nut_mean": nut_vol / nu, "nut_max": float(nut_max.max()) / nu,
        "frac1": vm(f1), "frac10": vm(f10),
        "S_mean": vm(S1), "S3_mean": vm(S3),
        "eps_nu": eps_nu_vol, "eps_sgs": eps_sgs_vol,
        "eps_tot": eps_tot_vol,
        "share": eps_sgs_vol / eps_tot_vol if eps_tot_vol > 0 else 0.0,
        "E": vm(ke), "eta": eta, "delta_over_eta": delta / eta,
        "iq": float(celik_iq_nu(1.0 + nut_vol / nu)), "iq_min": float(iq.min()),
        # PDF of log10(nu_t/nu)
        "pdf_x": 0.5 * (edges[1:] + edges[:-1]), "pdf": hist / dlog,
    }


# ----------------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------------
def print_sgs_block(s):
    print_block(f"LES / Smagorinsky SGS   it={s['iteration']:06d}  "
                f"t={s['time']:.4g}", [
        ("Cs (--cs)", s["cs"]),
        ("Lilly estimate (C_K = 1.5)", CS_LILLY),
        ("Filter width Delta", s["delta"]),
        ("  definition", s["delta_spec"]),
        ("  Delta/dx, /dy, /dz",
         "  ".join(f"{v:.3g}" for v in s["delta_over_d"])),
        ("<nu_t>/nu (volume)", s["nut_mean"]),
        ("max nu_t/nu", s["nut_max"]),
        ("volume fraction nu_t > nu", s["frac1"]),
        ("volume fraction nu_t > 10 nu", s["frac10"]),
        ("<|S|>", s["S_mean"]),
        ("eps_nu  = 2 nu <S_ij S_ij>", s["eps_nu"]),
        ("eps_sgs = <nu_t |S|^2>", s["eps_sgs"]),
        ("eps_sgs / (eps_nu + eps_sgs)", s["share"]),
        ("Delta / eta  (eps_tot)", s["delta_over_eta"]),
        ("LES_IQ_nu (Celik), mean", s["iq"]),
        ("LES_IQ_nu (Celik), min over y", s["iq_min"]),
        ("Energy <u_i u_i>/2 (volume)", s["E"]),
    ])


def sgs_warnings(s):
    msgs = []
    if s["cs"] <= 0.0:
        msgs.append("Cs <= 0: the SGS model is switched off.")
    elif s["nut_mean"] < 1e-3 and s["eps_tot"] > 0:
        msgs.append(f"<nu_t>/nu = {s['nut_mean']:.1e}: the model is almost "
                    "inactive (flow resolved like a DNS, or Cs/Delta very "
                    "small).")
    if s["iq"] < 0.8:
        msgs.append(f"LES_IQ_nu = {s['iq']:.2f} < 0.8 on average: poorly "
                    "resolved LES (nu_t/nu above about 20).")
    elif s["iq_min"] < 0.8:
        msgs.append(f"LES_IQ_nu drops to {s['iq_min']:.2f} < 0.8 somewhere "
                    "in y (see profile).")
    if s["share"] > 0.5:
        msgs.append(f"SGS dissipation is {100 * s['share']:.0f} % of the "
                    "total: the model dominates the dissipation.")
    if msgs:
        print("  Warnings (SGS):")
        for m in msgs:
            print(f"   * {m}")


# ----------------------------------------------------------------------------
# Time analysis: energy budget and inferred Cs
# ----------------------------------------------------------------------------
def sgs_time_analysis(history):
    if len(history) < 3:
        return None
    t = np.array([h["time"] for h in history], dtype=float)
    get = lambda k: np.array([h[k] for h in history], dtype=float)
    E, e_nu, e_sgs = get("E"), get("eps_nu"), get("eps_sgs")
    delta, S3 = get("delta"), get("S3_mean")
    dEdt = np.gradient(E, t, edge_order=2)
    eps_tot = e_nu + e_sgs
    resid = np.abs(dEdt + eps_tot) / np.maximum(eps_tot, 1e-300)
    num = -dEdt - e_nu
    with np.errstate(divide="ignore", invalid="ignore"):
        cs2 = num / (delta ** 2 * S3)
    cs_inf = np.where(cs2 > 0.0, np.sqrt(np.abs(cs2)), np.nan)
    return {"t": t, "dEdt": dEdt, "resid": resid, "cs_inf": cs_inf,
            "cs_in": get("cs"), "eps_tot": eps_tot,
            "share": e_sgs / np.maximum(eps_tot, 1e-300)}


def print_sgs_time_analysis(ta):
    ok = np.isfinite(ta["cs_inf"])
    ent = [("mean |dE/dt+eps_tot|/eps_tot", float(np.mean(ta["resid"]))),
           ("max  |dE/dt+eps_tot|/eps_tot", float(np.max(ta["resid"]))),
           ("Cs requested (--cs)", float(ta["cs_in"][0]))]
    if ok.any():
        ent += [("Cs inferred (mean)", float(np.nanmean(ta["cs_inf"]))),
                ("Cs inferred (min / max)",
                 f"{np.nanmin(ta['cs_inf']):.4g} / {np.nanmax(ta['cs_inf']):.4g}")]
    print_block("SGS time analysis: energy budget", ent)
    print("  dE/dt = -(eps_nu + eps_sgs)  (valid without forcing)")
    if float(np.mean(ta["share"])) < 0.1:
        print(f"  Caution: the SGS share is only {100 * np.mean(ta['share']):.1f} %"
              " of the dissipation, so Cs_inferred (a difference of two close "
              "quantities) amplifies the time-derivative error; use closely "
              "spaced snapshots.")
    print(f"   {'t':>10} {'-dE/dt':>12} {'eps_nu+eps_sgs':>15} "
          f"{'resid':>9} {'Cs_inferred':>12}")
    for n in range(len(ta["t"])):
        print(f"   {ta['t'][n]:10.4g} {-ta['dEdt'][n]:12.4e} "
              f"{ta['eps_tot'][n]:15.4e} {ta['resid'][n]:9.2e} "
              f"{ta['cs_inf'][n]:12.4f}")


# ----------------------------------------------------------------------------
# Files and figures
# ----------------------------------------------------------------------------
def save_sgs_history(history, ta, filename):
    keys = ["iteration", "time", "cs", "delta", "nut_mean", "nut_max",
            "frac1", "frac10", "eps_nu", "eps_sgs", "share", "delta_over_eta",
            "iq", "iq_min", "E"]
    cols = list(keys)
    extra = []
    if ta is not None:
        extra = [("dEdt", ta["dEdt"]), ("resid", ta["resid"]),
                 ("cs_inferred", ta["cs_inf"])]
    with open(filename, "w") as f:
        f.write("# " + " ".join(f"{c:>14}" for c in
                                cols + [e[0] for e in extra]) + "\n")
        for n, h in enumerate(history):
            vals = [h[k] for k in keys] + [e[1][n] for e in extra]
            f.write("  " + " ".join(f"{float(v):14.6e}" for v in vals) + "\n")
    print(f"\nSGS time series written to {filename}")


def save_sgs_profiles(history, outdir):
    for h in history:
        fn = outdir / f"les_profiles_{h['iteration']:06d}.dat"
        names = ["y", "nut_mean/nu", "nut_max/nu", "frac_nut>nu",
                 "frac_nut>10nu", "eps_nu", "eps_sgs", "share", "LES_IQ_nu"]
        cols = [h["y"], h["nut_mean_y"], h["nut_max_y"], h["frac1_y"],
                h["frac10_y"], h["eps_nu_y"], h["eps_sgs_y"], h["share_y"],
                h["iq_y"]]
        with open(fn, "w") as f:
            f.write(f"# t = {float(h['time']):.8e}   Cs = {h['cs']}   "
                    f"Delta = {h['delta']:.8e}\n")
            f.write("# " + " ".join(f"{c:>14}" for c in names) + "\n")
            for row in zip(*cols):
                f.write("  " + " ".join(f"{v:14.6e}" for v in row) + "\n")
        fp = outdir / f"les_pdf_{h['iteration']:06d}.dat"
        with open(fp, "w") as f:
            f.write("# log10(nu_t/nu)   pdf (per unit log10)\n")
            for x, p in zip(h["pdf_x"], h["pdf"]):
                f.write(f"  {x:14.6e} {p:14.6e}\n")
    print(f"SGS profiles / PDFs written to {outdir}/les_profiles_XXXXXX.dat, "
          "les_pdf_XXXXXX.dat")


def plot_sgs_snapshot(s, outfile):
    y = s["y"]
    fig, axs = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle(f"Smagorinsky SGS   it={s['iteration']:06d}   "
                 f"t={s['time']:.4g}   Cs={s['cs']}   "
                 rf"$\Delta$={s['delta']:.3g}", fontsize=11)

    ax = axs[0, 0]
    ok = s["pdf"] > 0
    ax.semilogy(10.0 ** s["pdf_x"][ok], s["pdf"][ok], "C0-")
    ax.set_xscale("log")
    ax.axvline(1.0, color="k", ls="--", lw=0.8, label=r"$\nu_t=\nu$")
    ax.axvline(max(s["nut_mean"], 1e-300), color="C3", ls=":", lw=1.2,
               label=rf"$\langle\nu_t\rangle/\nu$ = {s['nut_mean']:.3g}")
    ax.set_xlabel(r"$\nu_t/\nu$")
    ax.set_ylabel(r"PDF per unit $\log_{10}$")
    ax.legend(fontsize=8)
    ax.set_title(r"Distribution of the eddy viscosity")

    ax = axs[0, 1]
    ax.semilogx(np.maximum(s["nut_mean_y"], 1e-12), y, "C0-",
                label=r"$\langle\nu_t\rangle_{xz}/\nu$")
    ax.semilogx(np.maximum(s["nut_max_y"], 1e-12), y, "C1--",
                label=r"$\max_{xz}\,\nu_t/\nu$")
    ax.axvline(1.0, color="k", ls=":", lw=0.8)
    ax.set_xlabel(r"$\nu_t/\nu$")
    ax.set_ylabel("y")
    ax.legend(fontsize=8)
    ax.set_title("Eddy viscosity profile")

    ax = axs[1, 0]
    ax.plot(s["eps_nu_y"], y, "C0-", label=r"$\epsilon_\nu=2\nu\langle S_{ij}S_{ij}\rangle$")
    ax.plot(s["eps_sgs_y"], y, "C3-", label=r"$\epsilon_{sgs}=\langle\nu_t|S|^2\rangle$")
    ax.set_xlabel("dissipation")
    ax.set_ylabel("y")
    axt = ax.twiny()
    axt.plot(s["share_y"], y, "k--", lw=1)
    axt.set_xlim(0, 1)
    axt.set_xlabel(r"$\epsilon_{sgs}/\epsilon_{tot}$  (dashed)")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("Resolved and SGS dissipation")

    ax = axs[1, 1]
    ax.plot(s["iq_y"], y, "C2-", label=r"$LES\_IQ_\nu(y)$")
    ax.axvline(0.8, color="k", ls="--", lw=0.8, label="0.8 (acceptable LES)")
    ax.axvline(0.95, color="k", ls=":", lw=0.8, label="0.95 (~DNS)")
    ax.set_xlim(min(0.5, s["iq_min"] - 0.02), 1.0)
    ax.set_xlabel(r"$1/(1+0.05\,(\nu_{eff}/\nu)^{0.53})$")
    ax.set_ylabel("y")
    ax.legend(fontsize=8, loc="lower left")
    ax.set_title("LES quality index (Celik et al. 2005)")

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outfile, dpi=130)
    print(f"Figure saved to {outfile}")
    return fig


def plot_sgs_summary(history, ta, outfile):
    t = np.array([h["time"] for h in history], dtype=float)
    get = lambda k: np.array([h[k] for h in history], dtype=float)
    fig, axs = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("Smagorinsky SGS: time evolution", fontsize=11)

    ax = axs[0, 0]
    ax.semilogy(t, np.maximum(get("nut_mean"), 1e-12), "o-",
                label=r"$\langle\nu_t\rangle/\nu$")
    ax.semilogy(t, np.maximum(get("nut_max"), 1e-12), "s--",
                label=r"$\max\,\nu_t/\nu$")
    ax.axhline(1.0, color="k", ls=":", lw=0.8)
    ax.set_xlabel("t")
    ax.legend(fontsize=8)
    ax.set_title("Eddy viscosity")

    ax = axs[0, 1]
    ax.plot(t, get("share"), "o-", label=r"$\epsilon_{sgs}/\epsilon_{tot}$")
    ax.plot(t, get("iq"), "s-", label=r"$LES\_IQ_\nu$ (mean)")
    ax.plot(t, get("iq_min"), "^--", label=r"$LES\_IQ_\nu$ (min over y)")
    ax.axhline(0.8, color="k", ls="--", lw=0.8)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("t")
    ax.legend(fontsize=8)
    ax.set_title("SGS share and resolution quality")

    ax = axs[1, 0]
    ax.plot(t, get("eps_nu"), "o-", label=r"$\epsilon_\nu$")
    ax.plot(t, get("eps_sgs"), "s-", label=r"$\epsilon_{sgs}$")
    ax.plot(t, get("eps_nu") + get("eps_sgs"), "k-", lw=1.2,
            label=r"$\epsilon_\nu+\epsilon_{sgs}$")
    if ta is not None:
        ax.plot(t, -ta["dEdt"], "r--", lw=1.5, label=r"$-dE/dt$")
    ax.set_xlabel("t")
    ax.legend(fontsize=8)
    ax.set_title("Energy budget")

    ax = axs[1, 1]
    if ta is not None:
        ax.plot(t, ta["cs_inf"], "o-", label="Cs inferred from the budget")
    ax.plot(t, get("cs"), "k--", label="Cs requested")
    ax.axhline(CS_LILLY, color="C2", ls=":", label=f"Lilly ({CS_LILLY:.3f})")
    ax.set_ylim(bottom=0)
    ax.set_xlabel("t")
    ax.set_ylabel("Cs")
    ax.legend(fontsize=8)
    ax.set_title("Effective Smagorinsky constant")

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outfile, dpi=130)
    print(f"Figure saved to {outfile}")
    return fig


# ----------------------------------------------------------------------------
# Hooks used by the dispatcher
# ----------------------------------------------------------------------------
def sgs_report(s, args):
    if args.quiet:
        print(f"  [SGS] t={s['time']:.4g} <nu_t>/nu={s['nut_mean']:.3g} "
              f"share={s['share']:.3g} LES_IQ={s['iq']:.3f} "
              f"D/eta={s['delta_over_eta']:.3g}")
    else:
        print_sgs_block(s)
    sgs_warnings(s)
    if args.plot:
        args.outdir.mkdir(parents=True, exist_ok=True)
        fig = plot_sgs_snapshot(s, args.outdir / f"les_{s['iteration']:06d}.png")
        if args.no_show:
            plt.close(fig)


def sgs_finalize(history, args):
    ta = sgs_time_analysis(history)
    if ta is not None:
        print_sgs_time_analysis(ta)
    args.outdir.mkdir(parents=True, exist_ok=True)
    save_sgs_history(history, ta, args.outdir / "les_history.dat")
    save_sgs_profiles(history, args.outdir)
    if args.plot and len(history) >= 2:
        fig = plot_sgs_summary(history, ta, args.outdir / "les_summary.png")
        if args.no_show:
            plt.close(fig)
    if args.plot and not args.no_show:
        plt.show()
