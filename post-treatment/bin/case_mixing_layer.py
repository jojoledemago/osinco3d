#!/usr/bin/env python3
"""
case_mixing_layer.py -- temporal mixing layer (periodic x, z ; free-slip y).

Conventions
-----------
x : streamwise (periodic)     y : cross-stream, free-slip walls at y[0], y[-1]
z : spanwise (periodic)       U(y) = <ux>(x,z-average), dU = U(y_max)-U(y_min)

All statistics are plane averages over (x, z) (the two homogeneous
directions), computed slab by slab: no full 3D gradient tensor is stored.
Primes denote fluctuations about the plane average.

Quantities computed for each snapshot
-------------------------------------
Thicknesses  vorticity   d_w   = |dU| / max|dU/dy|
             momentum    theta = int (U1-U)(U-U2)/dU^2 dy
             visual      d_10-90 : distance between U = 10 % and 90 % of dU
             plus the equivalent widths of the TKE and of |<u'v'>|.
             Reynolds numbers Re_w = |dU| d_w / nu, Re_theta.
Turbulence   Reynolds stresses, k, production P = -<u'_i u'_j> dU_i/dx_j,
             dissipation eps = 2 nu <s'_ij s'_ij> (+ mean-flow part),
             centre-plane Re_lambda and Kolmogorov scale, resolution k_max*eta.
KH           amplitude of the streamwise Fourier modes of v' (2D: kz = 0, and
             all kz), 2D/3D energy partition, inviscid Rayleigh growth rate
             computed on the instantaneous mean profile, and the classical
             hyperbolic-tangent prediction (Michalke 1964) for reference.
Budgets      dE/dt = -int eps dy (exact for free-slip walls) and
             dK/dt = int (P - eps') dy, checked over the series of snapshots.
Box          two-point correlations at the centre plane (periodic box size).

Normalisation: lengths by d_w, velocities by |dU|, so that the usual
literature plots (self-similar profiles, growth rate d(d_w)/dt / |dU|) can be
compared directly.
"""

import numpy as np
import matplotlib.pyplot as plt

from fields_core import Case, VELOCITY, SCALAR_FIELD, print_block

# Inviscid linear theory of the profile U = (dU/2) tanh(y/h) (Michalke 1964):
# most unstable wavenumber alpha*h = 0.4446, growth rate 0.1897*(dU/2)/h.
# With d_w = 2h:   k*d_w = 0.889,  sigma = 0.1897 * dU / d_w,  neutral k*d_w = 2.
KDW_TANH = 2.0 * 0.4446
SIGMA_TANH = 0.1897            # in units of |dU| / d_w
KDW_NEUTRAL = 2.0

PAIRS = {"uu": (0, 0), "vv": (1, 1), "ww": (2, 2),
         "uv": (0, 1), "uw": (0, 2), "vw": (1, 2)}


# ----------------------------------------------------------------------------
# Small numerical helpers
# ----------------------------------------------------------------------------
def trapz(f, x):
    """Trapezoidal integral (local implementation, numpy-version proof)."""
    return float(np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(x)))


def crossing(y, f, level, near=None):
    """
    Position where f crosses `level` (linear interpolation).  If several
    crossings exist, the one closest to `near` is returned (first one if
    `near` is None).  NaN if there is none.
    """
    s = f - level
    idx = np.where(s[:-1] * s[1:] <= 0.0)[0]
    if idx.size == 0:
        return np.nan
    den = s[idx + 1] - s[idx]
    t = np.where(den == 0.0, 0.0, -s[idx] / np.where(den == 0.0, 1.0, den))
    pos = y[idx] + t * (y[idx + 1] - y[idx])
    if near is None or not np.isfinite(near):
        return float(pos[0])
    return float(pos[np.argmin(np.abs(pos - near))])


def refined_max(y, f):
    """Location and value of max(f), refined by a parabola through the
    three points around the discrete maximum (uniform grid)."""
    i = int(np.argmax(f))
    if 0 < i < f.size - 1:
        f0, f1, f2 = f[i - 1], f[i], f[i + 1]
        den = f0 - 2.0 * f1 + f2
        if den != 0.0:
            off = 0.5 * (f0 - f2) / den
            return (float(y[i] + off * (y[1] - y[0])),
                    float(f1 - 0.25 * (f0 - f2) * off))
    return float(y[i]), float(f[i])


# ----------------------------------------------------------------------------
# Plane statistics (one pass over the y-slabs)
# ----------------------------------------------------------------------------
def plane_statistics(ff, nu, nplanes):
    """
    Single pass over the y-slabs.  Returns the (x,z)-averaged profiles and the
    accumulators of the spectral analysis.

    G[i, j, y] = <d u_i / d x_j>  (plane-mean velocity gradient; only the
    j = y column is non-zero up to round-off).
    """
    g = ff.grid
    nx, ny, nz = g.n
    dy = g.d[1]
    w = np.full(ny, dy)                    # trapezoid weights along y
    w[0] *= 0.5
    w[-1] *= 0.5

    U3 = np.zeros((3, ny))                 # plane means of ux, uy, uz
    R = {key: np.zeros(ny) for key in PAIRS}
    G = np.zeros((3, 3, ny))
    s2 = np.zeros(ny)                      # <s'_ij s'_ij> (sum over i, j)
    Ev2d = np.zeros((nx, nz))              # int |v^(kx,kz)|^2 dy
    Ex = np.zeros((3, ny, nx))             # 1D spectra in kx, per plane
    Ez = np.zeros((3, ny, nz))             # 1D spectra in kz, per plane
    phi_m = np.zeros(ny) if ff.has_scalar else None
    phi_r = np.zeros(ny) if ff.has_scalar else None
    div_max = 0.0

    for slab in ff.iter_slabs(nplanes):
        sl = slice(slab.j0, slab.j1)
        vel = [slab[name] for name in VELOCITY]
        mean = [v.mean(axis=(0, 2)) for v in vel]
        fl = [v - m[None, :, None] for v, m in zip(vel, mean)]
        for c in range(3):
            U3[c, sl] = mean[c]
        for key, (a, b) in PAIRS.items():
            R[key][sl] = (fl[a] * fl[b]).mean(axis=(0, 2))

        if ff.has_scalar:
            ph = slab[SCALAR_FIELD]
            pm = ph.mean(axis=(0, 2))
            phi_m[sl] = pm
            phi_r[sl] = np.sqrt(((ph - pm[None, :, None]) ** 2).mean(axis=(0, 2)))

        # --- gradients on demand: diagonal first, then off-diagonal pairs ---
        m_loc = slab.nplanes
        Gl = np.zeros((3, 3, m_loc))
        s2l = np.zeros(m_loc)
        div = None
        for i in range(3):
            gii = slab.vel_grad(i, i)
            Gl[i, i] = gii.mean(axis=(0, 2))
            div = gii if div is None else div + gii
            s2l += ((gii - Gl[i, i][None, :, None]) ** 2).mean(axis=(0, 2))
            slab.clear()
        for i in range(3):
            for j in range(i + 1, 3):
                gij, gji = slab.vel_grad(i, j), slab.vel_grad(j, i)
                Gl[i, j] = gij.mean(axis=(0, 2))
                Gl[j, i] = gji.mean(axis=(0, 2))
                sij = 0.5 * (gij + gji)
                sm = 0.5 * (Gl[i, j] + Gl[j, i])
                s2l += 2.0 * ((sij - sm[None, :, None]) ** 2).mean(axis=(0, 2))
                slab.clear()
        G[:, :, sl] = Gl
        s2[sl] = s2l
        div_max = max(div_max, float(np.abs(div).max()))

        # --- spectral content of the fluctuations (x, z) --------------------
        for c in range(3):
            F = np.fft.fft2(fl[c], axes=(0, 2)) / (nx * nz)
            P = np.abs(F) ** 2                      # (nx, m, nz)
            Ex[c, sl, :] = P.sum(axis=2).T
            Ez[c, sl, :] = P.sum(axis=0)
            if c == 1:                              # v' : KH analysis
                Ev2d += np.einsum("j,ijk->ik", w[sl], P)

    Rfull = np.array([[R["uu"], R["uv"], R["uw"]],
                      [R["uv"], R["vv"], R["vw"]],
                      [R["uw"], R["vw"], R["ww"]]])
    Sbar = 0.5 * (G + G.transpose(1, 0, 2))
    return {
        "y": g.coords[1], "w": w, "U": U3[0], "V": U3[1], "W": U3[2],
        "R": R, "G": G,
        "k": 0.5 * (R["uu"] + R["vv"] + R["ww"]),
        "prod": -np.einsum("ijy,ijy->y", Rfull, G),
        "eps_turb": 2.0 * nu * s2,
        "eps_mean": 2.0 * nu * np.einsum("ijy,ijy->y", Sbar, Sbar),
        "Ev2d": Ev2d, "Ex": Ex, "Ez": Ez,
        "phi_mean": phi_m, "phi_rms": phi_r,
        "div_max": div_max,
    }


# ----------------------------------------------------------------------------
# Thicknesses
# ----------------------------------------------------------------------------
def thickness_analysis(y, U, dUdy):
    """Vorticity, momentum and 10-90 % thicknesses of the mean profile."""
    dU = float(U[-1] - U[0])
    nan = float("nan")
    out = {"dU": dU, "Uc": 0.5 * float(U[0] + U[-1]), "U_lo": float(U[0]),
           "U_hi": float(U[-1])}
    y_shear, gmax = refined_max(y, np.abs(dUdy))
    out.update(y_shear=y_shear, gmax=gmax)
    if abs(dU) < 1e-12 or gmax <= 0.0:
        out.update(delta_w=nan, theta=nan, d1090=nan, yc=nan, y10=nan,
                   y90=nan, shape=nan)
        return out
    s = (U - U[0]) / dU                       # 0 -> 1 across the layer
    delta_w = abs(dU) / gmax
    theta = trapz(s * (1.0 - s), y)           # (U1-U)(U-U2)/dU^2
    yc = crossing(y, U, out["Uc"], near=y_shear)
    y10 = crossing(y, s, 0.1, near=yc)
    y90 = crossing(y, s, 0.9, near=yc)
    out.update(delta_w=delta_w, theta=theta, d1090=abs(y90 - y10), yc=yc,
               y10=y10, y90=y90, shape=delta_w / theta if theta else nan)
    return out


# ----------------------------------------------------------------------------
# Inviscid linear stability (Rayleigh) on the instantaneous mean profile
# ----------------------------------------------------------------------------
def rayleigh_growth(y, U, k_list, npts=161):
    """
    Temporal growth rate sigma(k) = k * max Im(c) of the Rayleigh equation
        (U - c)(phi'' - k^2 phi) - U'' phi = 0 ,   phi = 0 at the walls
    (v = 0 on free-slip walls), second-order finite differences.  The profile
    is interpolated on `npts` uniform points if it has more.

    Returns (sigma, c_r) arrays, or (None, None) if npts <= 0.
    """
    if npts <= 0:
        return None, None
    if y.size > npts:
        yy = np.linspace(y[0], y[-1], npts)
        UU = np.interp(yy, y, U)
    else:
        yy, UU = y, U
    h = yy[1] - yy[0]
    ni = yy.size - 2
    Ui = UU[1:-1]
    Upp = (UU[2:] - 2.0 * UU[1:-1] + UU[:-2]) / h**2
    D2 = (np.diag(-2.0 * np.ones(ni)) + np.diag(np.ones(ni - 1), 1)
          + np.diag(np.ones(ni - 1), -1)) / h**2
    tol = 1e-6 * max(float(np.ptp(UU)), 1e-300)
    sigma = np.zeros(len(k_list))
    c_r = np.zeros(len(k_list))
    eye = np.eye(ni)
    for q, k in enumerate(k_list):
        B = D2 - k * k * eye
        A = Ui[:, None] * B - np.diag(Upp)
        ev = np.linalg.eigvals(np.linalg.solve(B, A))
        c = ev[np.argmax(ev.imag)]
        if c.imag > tol:
            sigma[q] = k * c.imag
            c_r[q] = c.real
    return sigma, c_r


# ----------------------------------------------------------------------------
# Snapshot analysis
# ----------------------------------------------------------------------------
def ml_analysis(ff, args):
    g = ff.grid
    if g.bc != ("periodic", "freeslip", "periodic"):
        raise ValueError("Mixing-layer analysis requires BC "
                         "periodic / freeslip / periodic in x, y, z.")
    nu = args.nu
    nx, ny, nz = g.n
    dxs = g.d
    Lx, Ly, Lz = g.length
    st = plane_statistics(ff, nu, args.slab)
    y, U, V, W = st["y"], st["U"], st["V"], st["W"]
    R, k = st["R"], st["k"]
    dUdy = st["G"][0, 1]
    th = thickness_analysis(y, U, dUdy)
    dw, dU = th["delta_w"], th["dU"]
    nan = float("nan")
    # a vanishing dU (no layer) must give NaN, not a ZeroDivisionError
    adU = abs(dU) if abs(dU) > 1e-12 else nan
    dU2 = dU**2 if abs(dU) > 1e-12 else nan

    # --- integrals over y ---------------------------------------------------
    E_mean = 0.5 * trapz(U**2 + V**2 + W**2, y)
    K_int = trapz(k, y)
    P_int = trapz(st["prod"], y)
    eps_t_int = trapz(st["eps_turb"], y)
    eps_m_int = trapz(st["eps_mean"], y)

    # --- turbulence statistics (normalised) --------------------------------
    def peak(f):
        yp, vp = refined_max(y, f)
        return vp, (yp - th["yc"]) / dw if np.isfinite(dw) else nan
    rms = {c: np.sqrt(np.maximum(R[c], 0.0)) for c in ("uu", "vv", "ww")}
    pk = {c: peak(rms[c]) for c in rms}
    # -<u'v'> sgn(dU): positive for a gradient-diffusion (down-gradient)
    # transfer; the peak is the extremum of largest MAGNITUDE, sign kept.
    uv_s = -R["uv"] * np.sign(dU if dU else 1.0)
    j_uv = int(np.argmax(np.abs(uv_s)))
    uv_pk = peak(np.abs(uv_s))
    uv_pk = (uv_pk[0] * np.sign(uv_s[j_uv] or 1.0), uv_pk[1])
    k_pk = peak(k)
    d_k = K_int / k_pk[0] if k_pk[0] > 0 else nan
    uv_abs = np.abs(R["uv"])
    d_uv = trapz(uv_abs, y) / uv_abs.max() if uv_abs.max() > 0 else nan

    # --- centre-plane quantities, Kolmogorov scale, resolution -------------
    yc = th["yc"]
    jc = int(np.argmin(np.abs(y - yc))) if np.isfinite(yc) else ny // 2
    eps_c = float(np.interp(yc, y, st["eps_turb"])) if np.isfinite(yc) else nan
    k_c = float(np.interp(yc, y, k)) if np.isfinite(yc) else nan
    if np.isfinite(eps_c) and eps_c > 0 and k_c > 0:
        eta_c = (nu**3 / eps_c) ** 0.25
        up = np.sqrt(2.0 * k_c / 3.0)
        lam_c = np.sqrt(15.0 * nu * up**2 / eps_c)
        re_lam_c = up * lam_c / nu
    else:
        eta_c = lam_c = re_lam_c = nan
    eps_pos = st["eps_turb"] > 1e-12 * max(st["eps_turb"].max(), 1e-300)
    with np.errstate(divide="ignore", invalid="ignore"):
        eta_y = np.where(eps_pos, (nu**3 / np.where(eps_pos,
                         st["eps_turb"], 1.0)) ** 0.25, np.inf)
    kmax_eta_y = [np.pi * eta_y / d for d in dxs]            # per direction
    kmax_eta_min = [float(a.min()) if eps_pos.any() else nan
                    for a in kmax_eta_y]

    # --- KH instability: measured modes ------------------------------------
    M = max(1, min(int(args.kh_modes), nx // 2 - 1))
    mi = np.arange(1, M + 1)
    kx = 2.0 * np.pi * mi / Lx
    E2 = st["Ev2d"]
    ekz0 = E2[:, 0]
    A2d = np.sqrt(ekz0[mi] + ekz0[-mi])                     # kz = 0
    Aall = np.sqrt(E2[mi, :].sum(axis=1) + E2[-mi, :].sum(axis=1))
    e_tot = float(E2.sum())
    e_2d = float(ekz0[1:].sum())                            # kx != 0, kz = 0
    e_3d = float(E2[:, 1:].sum())                           # kz != 0
    m_dom = int(mi[np.argmax(A2d)]) if A2d.max() > 0 else 0

    sig_ray, cr_ray = rayleigh_growth(y, U, kx, args.rayleigh_npts)
    if sig_ray is not None and sig_ray.max() > 0:
        m_ray = int(mi[np.argmax(sig_ray)])
        sig_ray_max = float(sig_ray.max())
    else:
        m_ray, sig_ray_max = 0, nan
    k_tanh = KDW_TANH / dw if np.isfinite(dw) else nan
    lam_tanh = 2.0 * np.pi / k_tanh if np.isfinite(k_tanh) else nan

    # --- box-size check: two-point correlations at the centre plane ----------
    def corr_half(E, n):
        out = []
        for c in range(3):
            e = E[c, jc, :]
            var = e.sum()
            if var <= 0.0:
                out.append(nan)
                continue
            r = np.real(np.fft.ifft(e)) * n
            out.append(float(r[n // 2] / r[0]))
        return out
    rx_half = corr_half(st["Ex"], nx)
    rz_half = corr_half(st["Ez"], nz)

    s = {
        "iteration": ff.iteration, "time": ff.time, "has_scalar": ff.has_scalar,
        "nx": nx, "ny": ny, "nz": nz, "Lx": Lx, "Ly": Ly, "Lz": Lz,
        "d": dxs, "nu": nu,
        # profiles
        "y": y, "U": U, "V": V, "W": W, "dUdy": dUdy,
        "uu": R["uu"], "vv": R["vv"], "ww": R["ww"], "uv": R["uv"],
        "uw": R["uw"], "vw": R["vw"], "k": k, "prod": st["prod"],
        "eps_turb": st["eps_turb"], "eps_mean": st["eps_mean"],
        "kmax_eta_y": kmax_eta_y,
        "phi_mean": st["phi_mean"], "phi_rms": st["phi_rms"],
        # thicknesses
        **th,
        "Re_w": adU * dw / nu, "Re_theta": adU * th["theta"] / nu,
        "Lx_dw": Lx / dw, "Ly_dw": Ly / dw, "Lz_dw": Lz / dw,
        "d_k": d_k, "d_uv": d_uv,
        # turbulence
        "urms_pk": pk["uu"][0] / adU, "vrms_pk": pk["vv"][0] / adU,
        "wrms_pk": pk["ww"][0] / adU,
        "urms_pk_y": pk["uu"][1], "vrms_pk_y": pk["vv"][1],
        "wrms_pk_y": pk["ww"][1],
        "uv_pk": uv_pk[0] / dU2, "uv_pk_y": uv_pk[1],
        "k_pk": k_pk[0] / dU2, "k_pk_y": k_pk[1],
        "eps_c": eps_c, "eta_c": eta_c, "lambda_c": lam_c,
        "Re_lambda_c": re_lam_c, "kmax_eta_min": kmax_eta_min,
        # KH
        "kh_m": mi, "kh_kx": kx, "kh_A2d": A2d, "kh_Aall": Aall,
        "kh_m_dom": m_dom, "kh_E2d": e_2d, "kh_E3d": e_3d, "kh_Etot": e_tot,
        "ray_sigma": sig_ray, "ray_cr": cr_ray, "ray_m": m_ray,
        "ray_sigma_max": sig_ray_max,
        "tanh_k": k_tanh, "tanh_lambda": lam_tanh,
        "tanh_m": Lx / lam_tanh if np.isfinite(lam_tanh) else nan,
        "tanh_sigma": SIGMA_TANH * adU / dw,
        # budgets
        "E_mean": E_mean, "K_int": K_int, "E_tot": E_mean + K_int,
        "P_int": P_int, "eps_t_int": eps_t_int, "eps_m_int": eps_m_int,
        "eps_tot_int": eps_t_int + eps_m_int,
        "div_max": st["div_max"],
        "Rx_half": rx_half, "Rz_half": rz_half,
    }
    if ff.has_scalar:
        pm = st["phi_mean"]
        dphi = pm[-1] - pm[0]
        if abs(dphi) > 1e-12:
            sp = (pm - pm[0]) / dphi
            yc_p = crossing(y, sp, 0.5, near=yc)
            s["delta_phi"] = abs(crossing(y, sp, 0.9, near=yc_p)
                                 - crossing(y, sp, 0.1, near=yc_p))
        else:
            s["delta_phi"] = nan
        s["phi_var_int"] = trapz(st["phi_rms"] ** 2, y)
    return s


# ----------------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------------
def print_ml_blocks(s):
    tag = f"it={s['iteration']:06d}  t={s['time']:.4g}"
    print_block(f"Mixing layer: thicknesses  {tag}", [
        ("dU = U(top) - U(bottom)", s["dU"]),
        ("Uc = (U_top + U_bottom)/2", s["Uc"]),
        ("y_c  (U = Uc)", s["yc"]),
        ("Vorticity thickness d_w", s["delta_w"]),
        ("Momentum thickness theta", s["theta"]),
        ("Visual thickness d_10-90", s["d1090"]),
        ("Shape factor d_w/theta", s["shape"]),
        ("d_w from TKE (int k / kmax)", s["d_k"]),
        ("d_w from |<u'v'>| (equiv.)", s["d_uv"]),
        ("Re_w = |dU| d_w / nu", s["Re_w"]),
        ("Re_theta = |dU| theta / nu", s["Re_theta"]),
        ("Lx / d_w", s["Lx_dw"]),
        ("Ly / d_w", s["Ly_dw"]),
        ("Lz / d_w", s["Lz_dw"]),
    ] + ([("Scalar d_10-90 (phi)", s["delta_phi"])] if s["has_scalar"] else []))

    print_block("Turbulence statistics (normalised by dU)", [
        ("max u'_rms / |dU|", s["urms_pk"]),
        ("max v'_rms / |dU|", s["vrms_pk"]),
        ("max w'_rms / |dU|", s["wrms_pk"]),
        ("max -<u'v'> / dU^2", s["uv_pk"]),
        ("  at (y-y_c)/d_w", s["uv_pk_y"]),
        ("max k / dU^2", s["k_pk"]),
        ("  at (y-y_c)/d_w", s["k_pk_y"]),
        ("eps (turb.) at y_c", s["eps_c"]),
        ("Kolmogorov eta at y_c", s["eta_c"]),
        ("Taylor lambda at y_c", s["lambda_c"]),
        ("Re_lambda at y_c", s["Re_lambda_c"]),
        ("min kmax*eta  (x, y, z)",
         "  ".join(f"{v:.3g}" for v in s["kmax_eta_min"])),
    ])

    m_dom = s["kh_m_dom"]
    ent = [
        ("Modes analysed (kx = 2 pi m/Lx)", f"m = 1 .. {s['kh_m'][-1]}"),
        ("Dominant measured mode m", m_dom if m_dom else "none"),
    ]
    if m_dom:
        kdw = s["kh_kx"][m_dom - 1] * s["delta_w"]
        ent += [("  k_x d_w of that mode", kdw),
                ("  wavelength / d_w", 2 * np.pi / kdw),
                ("  amplitude A_m (kz = 0)", s["kh_A2d"][m_dom - 1])]
    e_tot = max(s["kh_Etot"], 1e-300)
    ent += [
        ("v' energy: 2D (kz=0) fraction", s["kh_E2d"] / e_tot),
        ("v' energy: 3D (kz!=0) fraction", s["kh_E3d"] / e_tot),
        ("Tanh theory: k d_w", KDW_TANH),
        ("Tanh theory: wavelength / d_w", 2 * np.pi / KDW_TANH),
        ("Tanh theory: Lx / lambda", s["tanh_m"]),
        ("Tanh theory: sigma d_w/|dU|", SIGMA_TANH),
    ]
    if s["ray_sigma"] is not None:
        ent += [("Rayleigh: most unstable m", s["ray_m"] if s["ray_m"]
                 else "none (stable)"),
                ("Rayleigh: sigma_max d_w/|dU|",
                 s["ray_sigma_max"] * s["delta_w"] / abs(s["dU"]))]
        if s["ray_m"]:
            ent += [("Rayleigh: sigma_max (abs.)", s["ray_sigma_max"]),
                    ("Rayleigh: k d_w of that mode",
                     s["kh_kx"][s["ray_m"] - 1] * s["delta_w"])]
    print_block("Kelvin-Helmholtz instability", ent)

    print_block("Budgets and numerics", [
        ("Mean-flow energy  1/2 int U_i^2", s["E_mean"]),
        ("Turbulent energy  int k", s["K_int"]),
        ("Production  int P", s["P_int"]),
        ("Dissipation int eps (turb.)", s["eps_t_int"]),
        ("Dissipation int eps (mean)", s["eps_m_int"]),
        ("max|div u| (FD6)", f"{s['div_max']:.3e}"),
        ("Corr. at Lx/2 (u,v,w)",
         "  ".join(f"{v:.2f}" for v in s["Rx_half"])),
        ("Corr. at Lz/2 (u,v,w)",
         "  ".join(f"{v:.2f}" for v in s["Rz_half"])),
    ])


def ml_warnings(s, les=False):
    """Simple pass/fail hints.  With --les the DNS resolution criterion (kmax*eta from the
    molecular dissipation only) is not applicable and is skipped."""
    msgs = []
    if not np.isfinite(s["delta_w"]):
        msgs.append("No shear layer detected (dU ~ 0 or flat profile).")
        print("  Warnings:\n   * " + msgs[0])
        return
    gmax = s["gmax"]
    if abs(s["dUdy"][0]) > 1e-2 * gmax or abs(s["dUdy"][-1]) > 1e-2 * gmax:
        msgs.append("dU/dy at a wall exceeds 1 % of its maximum: the free "
                    "streams are not uniform (layer too thick for Ly).")
    if s["Ly_dw"] < 6.0:
        msgs.append(f"Ly/d_w = {s['Ly_dw']:.2f} < 6 : confinement by the "
                    "free-slip walls is likely.")
    kref = max(s["k_pk"], 1e-300) * s["dU"] ** 2
    if kref > 0 and (s["k"][0] > 1e-3 * kref or s["k"][-1] > 1e-3 * kref):
        msgs.append("Turbulent kinetic energy reaches the walls "
                    "(> 0.1 % of its peak).")
    if s["tanh_m"] < 1.0:
        msgs.append(f"Lx = {s['tanh_m']:.2f} KH wavelength(s): the box is "
                    "shorter than the most unstable wavelength.")
    turbulent = s["vrms_pk"] > 1e-3
    if turbulent:
        rz = max(abs(v) for v in s["Rz_half"] if np.isfinite(v)) \
            if any(np.isfinite(s["Rz_half"])) else 0.0
        if rz > 0.1:
            msgs.append(f"Spanwise correlation at Lz/2 = {rz:.2f} > 0.1 : "
                        "box too narrow in z.")
        kmin = [v for v in s["kmax_eta_min"] if np.isfinite(v)]
        if kmin and min(kmin) < 1.5 and not les:
            msgs.append(f"kmax*eta = {min(kmin):.2f} < 1.5 : DNS "
                        "under-resolved somewhere in the layer.")
    if s["ray_sigma"] is not None and s["ray_m"] == 0 \
            and s["kh_m_dom"] == 0:
        msgs.append("No unstable KH mode is resolved by the box "
                    "(Rayleigh stable for all modes analysed).")
    if msgs:
        print("  Warnings:")
        for m in msgs:
            print(f"   * {m}")


# ----------------------------------------------------------------------------
# Time series
# ----------------------------------------------------------------------------
def choose_kh_mode(history, args):
    """Mode followed in time: user choice, else the most unstable Rayleigh
    mode of the first snapshot, else its dominant measured mode."""
    M = len(history[0]["kh_m"])
    if args.kh_mode > 0:
        return min(args.kh_mode, M)
    h0 = history[0]
    if h0["ray_m"]:
        return h0["ray_m"]
    return h0["kh_m_dom"] if h0["kh_m_dom"] else 1


def ml_time_analysis(history, args):
    """Growth rates and budget checks over the series of snapshots."""
    if len(history) < 3:
        return None
    t = np.array([h["time"] for h in history], dtype=float)
    get = lambda key: np.array([h[key] for h in history], dtype=float)
    grad = lambda f: np.gradient(f, t, edge_order=2)

    dw, theta, adU = get("delta_w"), get("theta"), np.abs(get("dU"))
    out = {"t": t, "ddw_dt": grad(dw), "dtheta_dt": grad(theta)}
    out["growth_dw"] = out["ddw_dt"] / adU
    out["growth_theta"] = out["dtheta_dt"] / adU

    # self-similar growth rate: linear fit of d_w(t) over the last half
    half = len(t) // 2
    if len(t) - half >= 3 and np.all(np.isfinite(dw[half:])):
        slope = np.polyfit(t[half:], dw[half:], 1)[0]
        out["fit_growth_dw"] = float(slope / adU[half:].mean())
        slope_t = np.polyfit(t[half:], theta[half:], 1)[0]
        out["fit_growth_theta"] = float(slope_t / adU[half:].mean())
    else:
        out["fit_growth_dw"] = out["fit_growth_theta"] = float("nan")

    # KH: measured growth of one 2D mode vs Rayleigh prediction
    m_star = choose_kh_mode(history, args)
    A = np.array([h["kh_A2d"][m_star - 1] for h in history])
    out["m_star"] = m_star
    out["A_star"] = A
    with np.errstate(divide="ignore"):
        out["sigma_meas"] = grad(np.log(np.maximum(A, 1e-300)))
    out["sigma_ray"] = np.array([
        h["ray_sigma"][m_star - 1] if h["ray_sigma"] is not None
        else float("nan") for h in history])

    # global energy budget  dE/dt = - int eps dy   (exact, free-slip walls)
    E, eps = get("E_tot"), get("eps_tot_int")
    dEdt = grad(E)
    out["dEdt"] = dEdt
    out["rel_E"] = np.abs(dEdt + eps) / np.maximum(eps, 1e-300)
    # turbulent kinetic energy  dK/dt = int (P - eps')
    K, rhs = get("K_int"), get("P_int") - get("eps_t_int")
    dKdt = grad(K)
    out["dKdt"], out["rhsK"] = dKdt, rhs
    scale = np.maximum(np.maximum(get("P_int"), get("eps_t_int")), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        out["rel_K"] = np.where(scale > 1e-12 * max(scale.max(), 1e-300),
                                np.abs(dKdt - rhs) / scale, np.nan)
    return out


def print_ml_time_analysis(ta):
    ent = [("KH mode followed (m)", ta["m_star"]),
           ("mean d(d_w)/dt / |dU| (2nd half)", ta["fit_growth_dw"]),
           ("mean d(theta)/dt / |dU| (2nd half)", ta["fit_growth_theta"]),
           ("mean |dE/dt+eps|/eps", float(np.mean(ta["rel_E"]))),
           ("max  |dE/dt+eps|/eps", float(np.max(ta["rel_E"])))]
    if np.any(np.isfinite(ta["rel_K"])):
        ent += [("mean |dK/dt-(P-eps')|/scale", float(np.nanmean(ta["rel_K"]))),
                ("max  |dK/dt-(P-eps')|/scale", float(np.nanmax(ta["rel_K"])))]
    print_block("Time analysis (mixing layer)", ent)
    print("  Growth rate of the followed KH mode "
          "(valid while the mode is in its linear stage):")
    print(f"   {'t':>10} {'A_m':>12} {'sigma_meas':>12} "
          f"{'sigma_Rayleigh':>15} {'ratio':>8}")
    for t, a, sm, sr in zip(ta["t"], ta["A_star"], ta["sigma_meas"],
                            ta["sigma_ray"]):
        ratio = sm / sr if np.isfinite(sr) and sr > 0 else float("nan")
        print(f"   {t:10.4g} {a:12.4e} {sm:12.4e} {sr:15.4e} {ratio:8.3f}")


# ----------------------------------------------------------------------------
# Files (gnuplot-friendly)
# ----------------------------------------------------------------------------
def save_ml_history(history, ta, filename):
    cols = ["it", "t", "dU", "delta_w", "theta", "d_10-90", "shape", "Re_w",
            "Re_theta", "yc", "urms_pk", "vrms_pk", "wrms_pk", "uv_pk",
            "k_pk", "K_int", "E_mean", "P_int", "eps_t_int", "eps_m_int",
            "div_max", "kh_m_dom", "ray_m", "tanh_sigma"]
    keys = ["iteration", "time", "dU", "delta_w", "theta", "d1090", "shape",
            "Re_w", "Re_theta", "yc", "urms_pk", "vrms_pk", "wrms_pk",
            "uv_pk", "k_pk", "K_int", "E_mean", "P_int", "eps_t_int",
            "eps_m_int", "div_max", "kh_m_dom", "ray_m", "tanh_sigma"]
    extra = []
    if ta is not None:
        extra = [("growth_dw", ta["growth_dw"]),
                 ("growth_theta", ta["growth_theta"]),
                 ("sigma_meas", ta["sigma_meas"]),
                 ("sigma_ray", ta["sigma_ray"])]
    with open(filename, "w") as f:
        f.write("# " + " ".join(f"{c:>14}" for c in
                                cols + [e[0] for e in extra]) + "\n")
        for n, h in enumerate(history):
            vals = [h[k] for k in keys] + [e[1][n] for e in extra]
            f.write("  " + " ".join(f"{float(v):14.6e}" for v in vals) + "\n")
    print(f"\nTime series written to {filename}")


def save_ml_profiles(history, outdir):
    """One text file per snapshot: plane-averaged profiles."""
    names = ["y", "U", "V", "W", "uu", "vv", "ww", "uv", "uw", "vw", "k",
             "prod", "eps_turb", "eps_mean"]
    for h in history:
        cols = [h[n] for n in names]
        labels = list(names)
        if h["has_scalar"]:
            cols += [h["phi_mean"], h["phi_rms"]]
            labels += ["phi_mean", "phi_rms"]
        fn = outdir / f"ml_profiles_{h['iteration']:06d}.dat"
        with open(fn, "w") as f:
            f.write(f"# t = {float(h['time']):.8e}\n")
            f.write("# " + " ".join(f"{c:>14}" for c in labels) + "\n")
            for row in zip(*cols):
                f.write("  " + " ".join(f"{v:14.6e}" for v in row) + "\n")
    print(f"Profiles written to {outdir}/ml_profiles_XXXXXX.dat")


def save_ml_kh_modes(history, filename):
    M = len(history[0]["kh_m"])
    with open(filename, "w") as f:
        head = (["it", "t"] + [f"A2d_m{m}" for m in range(1, M + 1)]
                + [f"Aall_m{m}" for m in range(1, M + 1)]
                + [f"sig_ray_m{m}" for m in range(1, M + 1)])
        f.write("# " + " ".join(f"{c:>14}" for c in head) + "\n")
        for h in history:
            sr = h["ray_sigma"] if h["ray_sigma"] is not None \
                else np.full(M, np.nan)
            vals = [h["iteration"], h["time"], *h["kh_A2d"], *h["kh_Aall"],
                    *sr]
            f.write("  " + " ".join(f"{float(v):14.6e}" for v in vals) + "\n")
    print(f"KH mode amplitudes written to {filename}")


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
def _ylim(s, zoom):
    y = s["y"]
    if zoom > 0 and np.isfinite(s["yc"]) and np.isfinite(s["delta_w"]):
        return (max(y[0], s["yc"] - zoom * s["delta_w"]),
                min(y[-1], s["yc"] + zoom * s["delta_w"]))
    return (y[0], y[-1])


def plot_ml_snapshot(s, outfile, zoom=0.0):
    """Six panels for one snapshot (profiles are plotted against y)."""
    y, dw, yc, theta = s["y"], s["delta_w"], s["yc"], s["theta"]
    adU, dU = abs(s["dU"]), s["dU"]
    ylim = _ylim(s, zoom)
    fig, axs = plt.subplots(2, 3, figsize=(17, 9.5))
    fig.suptitle(f"Mixing layer   it={s['iteration']:06d}   t={s['time']:.4g}"
                 f"    Re_w={s['Re_w']:.4g}", fontsize=11)

    # (0,0) mean streamwise velocity and thicknesses -------------------------
    ax = axs[0, 0]
    ax.plot(s["U"], y, "k-", lw=1.6, label=r"$\langle u_x\rangle$")
    if np.isfinite(dw):
        ax.plot([s["U_lo"], s["U_hi"]], [yc - dw / 2, yc + dw / 2], "C3--",
                lw=1, label=rf"$\delta_\omega$ = {dw:.3g}  (max-slope tangent)")
        ax.axhline(s["y10"], color="C2", ls=":", lw=1,
                   label=rf"$\delta_{{10-90}}$ = {s['d1090']:.3g}")
        ax.axhline(s["y90"], color="C2", ls=":", lw=1)
        ax.axhline(yc - 2 * theta, color="C4", ls="-.", lw=1,
                   label=rf"$\theta$ = {theta:.3g}  (band $y_c\pm2\theta$)")
        ax.axhline(yc + 2 * theta, color="C4", ls="-.", lw=1)
        ax.axvline(s["Uc"], color="0.7", lw=0.8)
        ax.plot([s["Uc"]], [yc], "ko", ms=4)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$\langle u_x\rangle$")
    ax.set_ylabel("y")
    ax.legend(fontsize=7, loc="upper left")
    ax2 = ax.twiny()
    ax2.plot(s["dUdy"], y, color="0.55", lw=1)
    ax2.set_xlabel(r"$d\langle u_x\rangle/dy$  (grey)", color="0.4")
    ax.set_title("Mean velocity and thicknesses")

    # (0,1) cross-stream and spanwise velocities --------------------------------
    ax = axs[0, 1]
    for key, rk, col, lab in (("V", "vv", "C0", "u_y"), ("W", "ww", "C1", "u_z")):
        m, r = s[key], np.sqrt(np.maximum(s[rk], 0.0))
        ax.plot(m, y, color=col, label=rf"$\langle {lab}\rangle$")
        ax.fill_betweenx(y, m - r, m + r, color=col, alpha=0.2,
                         label=rf"$\pm$ rms of ${lab}'$")
    ax.axvline(0.0, color="0.6", lw=0.8)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$\langle u_y\rangle,\ \langle u_z\rangle$")
    ax.set_ylabel("y")
    ax.legend(fontsize=7)
    ax.set_title("Cross-stream / spanwise velocity")

    # (0,2) fluctuation intensities and shear stress ------------------------------
    ax = axs[0, 2]
    for rk, col, lab in (("uu", "C0", "u'"), ("vv", "C1", "v'"),
                         ("ww", "C2", "w'")):
        ax.plot(np.sqrt(np.maximum(s[rk], 0.0)) / adU, y, color=col,
                label=rf"${lab}_{{rms}}/|\Delta U|$")
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$u'_{rms}/|\Delta U|$")
    ax.set_ylabel("y")
    axt = ax.twiny()
    axt.plot(-s["uv"] * np.sign(dU or 1.0) / dU**2, y, "k--",
             label=r"$-\langle u'v'\rangle\,\mathrm{sgn}(\Delta U)/\Delta U^2$")
    axt.set_xlabel(r"$-\langle u'v'\rangle/\Delta U^2$  (dashed)")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = axt.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7)
    ax.set_title("Turbulence intensities and Reynolds shear stress")

    # (1,0) Kelvin-Helmholtz modes ------------------------------------------------
    ax = axs[1, 0]
    kdw = s["kh_kx"] * dw
    pos = s["kh_A2d"] > 0
    if pos.any():
        ax.semilogy(kdw[pos], s["kh_A2d"][pos], "o-", color="C0",
                    label=r"$A_m$  ($k_z=0$, 2D)")
    pos = s["kh_Aall"] > 0
    if pos.any():
        ax.semilogy(kdw[pos], s["kh_Aall"][pos], "s--", color="C1", ms=4,
                    label=r"$A_m$  (all $k_z$)")
    ax.axvline(KDW_TANH, color="0.5", ls=":", lw=1)
    ax.axvline(KDW_NEUTRAL, color="0.5", ls=":", lw=1)
    ax.text(KDW_TANH, 0.97, "tanh max", rotation=90, va="top", ha="right",
            fontsize=7, transform=ax.get_xaxis_transform())
    ax.text(KDW_NEUTRAL, 0.97, "tanh neutral", rotation=90, va="top",
            ha="right", fontsize=7, transform=ax.get_xaxis_transform())
    amax = max(float(s["kh_A2d"].max()), float(s["kh_Aall"].max()))
    if amax > 0:                         # hide round-off floor (10 decades)
        ax.set_ylim(bottom=amax * 1e-10, top=amax * 3.0)
    ax.set_xlabel(r"$k_x\,\delta_\omega$")
    ax.set_ylabel(r"$A_m=\left(\int|\hat v'(k_x)|^2dy\right)^{1/2}$")
    if s["ray_sigma"] is not None:
        axr = ax.twinx()
        axr.plot(kdw, s["ray_sigma"] * dw / adU, "^-", color="C3", ms=5,
                 label=r"Rayleigh $\sigma\,\delta_\omega/|\Delta U|$")
        axr.axhline(SIGMA_TANH, color="C3", ls=":", lw=0.8)
        axr.set_ylabel(r"$\sigma\,\delta_\omega/|\Delta U|$", color="C3")
        axr.set_ylim(bottom=0)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = axr.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="lower left")
    else:
        ax.legend(fontsize=7)
    ax.set_title("KH modes: measured amplitude and linear growth rate")

    # (1,1) TKE production and dissipation -------------------------------------------
    ax = axs[1, 1]
    sc = dw / adU**3 if np.isfinite(dw) else 1.0
    ax.plot(s["prod"] * sc, y, "C3", label=r"$P$")
    ax.plot(s["eps_turb"] * sc, y, "C0", label=r"$\epsilon'$ (turbulent)")
    ax.plot(s["eps_mean"] * sc, y, "C2", ls="--", label=r"$\epsilon$ (mean flow)")
    ax.axvline(0.0, color="0.6", lw=0.8)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$\times\ \delta_\omega/|\Delta U|^3$")
    ax.set_ylabel("y")
    ax.legend(fontsize=7)
    ax.set_title("TKE production and dissipation")

    # (1,2) resolution ----------------------------------------------------------------
    ax = axs[1, 2]
    for a, lab, col in zip(s["kmax_eta_y"], "xyz", ("C0", "C1", "C2")):
        ok = np.isfinite(a)
        if ok.any():
            ax.semilogx(a[ok], y[ok], color=col,
                        label=rf"$k_{{max,{lab}}}\,\eta(y)$")
    ax.axvline(1.5, color="k", ls="--", lw=0.8, label="1.5")
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$\pi\,\eta(y)/\Delta_i$")
    ax.set_ylabel("y")
    ax.legend(fontsize=7)
    ax.set_title("Resolution of the dissipative scales")

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outfile, dpi=130)
    print(f"Figure saved to {outfile}")
    return fig


def plot_ml_summary(history, ta, outfile):
    """Time evolution and self-similarity (needs several snapshots)."""
    t = np.array([h["time"] for h in history], dtype=float)
    get = lambda key: np.array([h[key] for h in history], dtype=float)
    fig, axs = plt.subplots(2, 3, figsize=(17, 9.5))
    fig.suptitle("Mixing layer: time evolution", fontsize=11)

    ax = axs[0, 0]
    ax.plot(t, get("delta_w"), "o-", label=r"$\delta_\omega$")
    ax.plot(t, get("d1090"), "s-", label=r"$\delta_{10-90}$")
    ax.plot(t, 4.0 * get("theta"), "^-", label=r"$4\theta$")
    ax.set_xlabel("t")
    ax.set_ylabel("thickness")
    ax.legend(fontsize=8)
    ax.set_title("Thicknesses")

    ax = axs[0, 1]
    if ta is not None:
        ax.plot(t, ta["growth_dw"], "o-", label=r"$d\delta_\omega/dt\,/|\Delta U|$")
        ax.plot(t, ta["growth_theta"], "s-",
                label=r"$d\theta/dt\,/|\Delta U|$")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "needs >= 3 snapshots", ha="center",
                transform=ax.transAxes)
    ax.set_xlabel("t")
    ax.set_title("Growth rates")

    ax = axs[0, 2]
    M = min(4, len(history[0]["kh_m"]))
    for m in range(1, M + 1):
        a = np.array([h["kh_A2d"][m - 1] for h in history])
        if np.any(a > 0):
            ax.semilogy(t, np.where(a > 0, a, np.nan), "o-", ms=3,
                        label=f"m={m}")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$A_m$ ($k_z=0$)")
    ax.legend(fontsize=8)
    ax.set_title("KH mode amplitudes")

    ax = axs[1, 0]
    if ta is not None:
        ax.plot(t, ta["sigma_meas"], "o-", label=f"measured (m={ta['m_star']})")
        if np.any(np.isfinite(ta["sigma_ray"])):
            ax.plot(t, ta["sigma_ray"], "^--", color="C3",
                    label="Rayleigh (mean profile)")
        ax.axhline(0.0, color="0.6", lw=0.8)
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "needs >= 3 snapshots", ha="center",
                transform=ax.transAxes)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$\sigma$")
    ax.set_title("KH growth rate: measured vs linear theory")

    ax = axs[1, 1]
    if ta is not None:
        ax.plot(t, ta["dEdt"], "o-", label=r"$dE/dt$")
        ax.plot(t, -get("eps_tot_int"), "k:", lw=1.5, label=r"$-\int\epsilon\,dy$")
        ax2 = ax.twinx()
        ax2.plot(t, ta["dKdt"], "s-", color="C1", ms=4, label=r"$dK/dt$")
        ax2.plot(t, ta["rhsK"], "^--", color="C3", ms=4,
                 label=r"$\int(P-\epsilon')dy$")
        ax2.set_ylabel("TKE budget")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=7)
    else:
        ax.text(0.5, 0.5, "needs >= 3 snapshots", ha="center",
                transform=ax.transAxes)
    ax.set_xlabel("t")
    ax.set_ylabel("energy budget")
    ax.set_title("Energy budgets")

    ax = axs[1, 2]
    idx = np.unique(np.linspace(0, len(history) - 1,
                                min(len(history), 6)).astype(int))
    colors = plt.cm.viridis(np.linspace(0, 1, len(idx)))
    axr = ax.twinx()
    for c, i in zip(colors, idx):
        h = history[i]
        if not np.isfinite(h["delta_w"]):
            continue
        eta = (h["y"] - h["yc"]) / h["delta_w"]
        ax.plot(eta, (h["U"] - h["Uc"]) / h["dU"], color=c,
                label=f"t={h['time']:.3g}")
        axr.plot(eta, -h["uv"] * np.sign(h["dU"] or 1.0) / h["dU"] ** 2,
                 "--", color=c, lw=1)
    ax.set_xlabel(r"$(y-y_c)/\delta_\omega$")
    ax.set_ylabel(r"$(\langle u_x\rangle-U_c)/\Delta U$")
    axr.set_ylabel(r"$-\langle u'v'\rangle/\Delta U^2$  (dashed)")
    ax.legend(fontsize=7)
    ax.set_title("Self-similarity")

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(outfile, dpi=130)
    print(f"Figure saved to {outfile}")
    return fig


# ----------------------------------------------------------------------------
# Case
# ----------------------------------------------------------------------------
class MixingLayerCase(Case):
    """Temporal mixing layer: periodic x and z, free-slip y."""
    name = "mixing_layer"
    bc = ("periodic", "freeslip", "periodic")

    def analyse(self, ff, args):
        return ml_analysis(ff, args)

    def report(self, s, args):
        if args.quiet:
            print(f"  t={s['time']:.4g} d_w={s['delta_w']:.4g} "
                  f"theta={s['theta']:.4g} Re_w={s['Re_w']:.4g} "
                  f"m_dom={s['kh_m_dom']} "
                  f"vrms_pk={s['vrms_pk']:.3g} "
                  f"div={s['div_max']:.2e}")
        else:
            print_ml_blocks(s)
        ml_warnings(s, les=args.les)
        if args.plot:
            args.outdir.mkdir(parents=True, exist_ok=True)
            fig = plot_ml_snapshot(
                s, args.outdir / f"ml_{s['iteration']:06d}.png",
                zoom=args.ml_zoom)
            if args.no_show:
                plt.close(fig)

    def finalize(self, history, args):
        ta = ml_time_analysis(history, args)
        if ta is not None:
            print_ml_time_analysis(ta)
        save_ml_history(history, ta, args.outdir / "ml_history.dat")
        save_ml_kh_modes(history, args.outdir / "ml_kh_modes.dat")
        save_ml_profiles(history, args.outdir)
        if args.plot:
            fig = plot_ml_summary(history, ta, args.outdir / "ml_summary.png")
            if args.no_show:
                plt.close(fig)
            else:
                plt.show()
