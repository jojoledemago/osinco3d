#!/usr/bin/env python3
"""
fields_core.py -- common core of the osinco3d post-processing tools.

Contents
--------
* Printing helpers and file discovery.
* Grid          : uniform Cartesian grid + boundary condition per direction.
* FieldFile     : lazy (memory-mapped) reader of fields_XXXXXX.bin, with an
                  optional passive scalar `phi`.
* Slab          : a block of y-planes (+ ghost planes) giving gradients ON
                  DEMAND, so that a whole 3D gradient tensor never has to sit
                  in memory.
* Case          : minimal interface that every case analysis implements
                  (hit, mixing_layer, ...), used by the dispatcher.

File format (Fortran stream, unformatted)
-----------------------------------------
    time (float64) | nx, ny, nz (int32) | x(nx) y(ny) z(nz) (float64)
    ux, uy, uz, p [, phi]   each nx*ny*nz float64, Fortran order

`phi` is optional: its presence is detected from the file size
(4 or 5 arrays after the header).

Boundary conditions (per direction)
-----------------------------------
"periodic" : n distinct points, the stencils wrap modulo n, so the period
             seen by the solver is L = n*dx.
"freeslip" : the grid contains both boundary nodes (dy = Ly/(n-1)).  Ghost
             points are obtained by MIRRORING about the boundary node:
               - the velocity component normal to the boundary is ODD,
               - every other field (tangential velocity, p, phi) is EVEN.
             This is the usual symmetry/antisymmetry closure of free-slip
             boundaries with centred schemes; check that it matches the
             closure used by your Fortran derivative routine.

Slabs are cut along y; any combination of "periodic" / "freeslip" is
supported in the three directions (mixing layer, HIT, free-slip TGV box...).
Derivatives along x and z use the same mirror closure as along y.
"""

import re
from pathlib import Path

import numpy as np

VELOCITY = ("ux", "uy", "uz")
BASE_FIELDS = VELOCITY + ("p",)
SCALAR_FIELD = "phi"

AXES = "xyz"
BC_CHOICES = ("periodic", "freeslip")

# 6th-order centred first derivative: coefficients of f(i+m) - f(i-m)
FD6 = (45.0 / 60.0, -9.0 / 60.0, 1.0 / 60.0)
HALO = len(FD6)                       # ghost planes needed on each side

_FILENAME_RE = re.compile(r"^fields_(\d{6})\.bin$")

# ----------------------------------------------------------------------------
# Printing helpers
# ----------------------------------------------------------------------------
def format_value(value):
    """Format numerical values with 4 significant digits."""
    if isinstance(value, (float, np.floating)):
        return f"{value:.4g}"
    return value


def print_block(title, entries):
    """Print a well-formatted block of labeled entries."""
    print(f"\n--- {title} ---")
    for label, value in entries:
        print(f"  {label:<30}: {format_value(value)}")
    print("-" * 48)


# ----------------------------------------------------------------------------
# Grid and field reader
# ----------------------------------------------------------------------------
def parity(name, axis):
    """
    Parity of field `name` under a free-slip mirror in direction `axis`:
    -1 (odd) for the velocity component normal to the boundary, +1 (even)
    for everything else (tangential velocity, p, phi).
    """
    return -1 if name == VELOCITY[axis] else 1


class Grid:
    """Uniform Cartesian grid with one boundary condition per direction."""

    def __init__(self, coords, bc, drop_last=False):
        bc = tuple(bc)
        if len(bc) != 3 or any(b not in BC_CHOICES for b in bc):
            raise ValueError(f"bc must be 3 values among {BC_CHOICES}, "
                             f"got {bc}")
        self.bc = bc
        self.n_file = tuple(int(c.size) for c in coords)
        # extent exactly as in the file: x[-1] - x[0]
        self.extent = tuple(float(c[-1] - c[0]) for c in coords)

        d, n, xyz = [], [], []
        for ax, c in enumerate(coords):
            if c.size < 2:
                raise ValueError(f"Direction {AXES[ax]}: fewer than 2 points.")
            h = float(c[1] - c[0])
            if np.max(np.abs(np.diff(c) - h)) > 1e-6 * abs(h):
                raise NotImplementedError(
                    f"Direction {AXES[ax]}: non-uniform grid (stretching) is "
                    "not supported yet.")
            # legacy files: last point duplicates the first (periodic only)
            ne = c.size - 1 if (drop_last and bc[ax] == "periodic") else c.size
            if bc[ax] == "freeslip" and ne <= HALO:
                raise ValueError(f"Direction {AXES[ax]}: need more than "
                                 f"{HALO} points for free-slip ghost cells.")
            d.append(h)
            n.append(int(ne))
            xyz.append(c[:ne])
        self.d = tuple(d)
        self.n = tuple(n)
        self.coords = tuple(xyz)

    def weights(self, axis):
        """Quadrature weights (sum = 1) along `axis` for volume averages:
        uniform if periodic, trapezoid if free-slip (both boundary nodes are
        stored, so they carry half a cell)."""
        n = self.n[axis]
        if self.bc[axis] == "freeslip":
            w = np.ones(n)
            w[0] = w[-1] = 0.5
            return w / w.sum()
        return np.full(n, 1.0 / n)

    @property
    def length(self):
        """Length per direction: n*d if periodic (period seen by the
        stencils), file extent otherwise."""
        return tuple(n * h if b == "periodic" else e
                     for n, h, b, e in zip(self.n, self.d, self.bc,
                                           self.extent))


class FieldFile:
    """
    Lazy reader of one fields_XXXXXX.bin file.

    Fields are memory-mapped: nothing is read until it is accessed.

    Attributes
    ----------
    path, iteration, time, grid, fields (tuple of available names)
    """

    def __init__(self, filename, bc=("periodic",) * 3, drop_last=False):
        self.path = Path(filename)
        match = _FILENAME_RE.match(self.path.name)
        if match is None:
            raise ValueError(f"Invalid field filename format: {filename}")
        self.iteration = int(match.group(1))

        with self.path.open("rb") as f:
            t = np.fromfile(f, dtype=np.float64, count=1)
            dims = np.fromfile(f, dtype=np.int32, count=3)
            if t.size != 1 or dims.size != 3:
                raise ValueError(f"Invalid header in file: {filename}")
            nx, ny, nz = (int(d) for d in dims)
            coords = [np.fromfile(f, dtype=np.float64, count=n)
                      for n in (nx, ny, nz)]
            if any(c.size != n for c, n in zip(coords, (nx, ny, nz))):
                raise ValueError(f"Unexpected end of file while reading: "
                                 f"{filename}")
        self.time = t[0]
        self.grid = Grid(coords, bc, drop_last=drop_last)
        self._shape = (nx, ny, nz)

        # --- which fields are stored? (phi is optional) --------------------
        n_values = nx * ny * nz
        header_bytes = 8 + 3 * 4 + 8 * (nx + ny + nz)
        payload = self.path.stat().st_size - header_bytes
        n_arrays, rest = divmod(payload, 8 * n_values)
        if rest != 0 or n_arrays not in (len(BASE_FIELDS),
                                         len(BASE_FIELDS) + 1):
            raise ValueError(
                f"{filename}: {payload} bytes after the header, expected "
                f"{len(BASE_FIELDS)} or {len(BASE_FIELDS) + 1} arrays of "
                f"{nx}x{ny}x{nz} float64 (ux,uy,uz,p[,phi]).")
        self.fields = BASE_FIELDS + ((SCALAR_FIELD,)
                                     if n_arrays > len(BASE_FIELDS) else ())

        self._mm = {
            name: np.memmap(self.path, dtype=np.float64, mode="r",
                            offset=header_bytes + 8 * n_values * k,
                            shape=self._shape, order="F")
            for k, name in enumerate(self.fields)
        }

    # -- access ---------------------------------------------------------------
    @property
    def has_scalar(self):
        return SCALAR_FIELD in self.fields

    def has(self, name):
        return name in self._mm

    def raw(self, name):
        """Full memory-mapped array, exactly as stored in the file."""
        if name not in self._mm:
            raise KeyError(f"Field '{name}' is not present in {self.path.name} "
                           f"(available: {', '.join(self.fields)}).")
        return self._mm[name]

    def array(self, name):
        """Memory-mapped view restricted to the points used by the analysis
        (differs from raw() only with --drop-last)."""
        nx, ny, nz = self.grid.n
        return self.raw(name)[:nx, :ny, :nz]

    def load(self, name):
        """Whole field copied into memory (for global operations such as
        3D FFTs)."""
        return np.array(self.array(name))

    def iter_slabs(self, nplanes=16):
        """Yield Slab objects covering all y-planes, `nplanes` at a time."""
        g = self.grid
        ny = g.n[1]
        nplanes = max(1, int(nplanes))
        for j0 in range(0, ny, nplanes):
            yield Slab(self, j0, min(j0 + nplanes, ny))


# ----------------------------------------------------------------------------
# Slabs: y-planes with ghost cells, gradients on demand
# ----------------------------------------------------------------------------
def ghost_indices(j0, j1, n, bc, halo=HALO):
    """
    Indices (and reflection mask) of planes j0-halo ... j1+halo-1 mapped into
    [0, n-1] according to the boundary condition along the slab direction.
    """
    j = np.arange(j0 - halo, j1 + halo)
    if bc == "periodic":
        return j % n, np.zeros(j.size, dtype=bool)
    # free-slip: mirror about the boundary nodes 0 and n-1
    refl = (j < 0) | (j > n - 1)
    jm = np.where(j < 0, -j, j)
    jm = np.where(jm > n - 1, 2 * (n - 1) - jm, jm)
    return jm, refl


def fd6_from_extended(ext, d, axis, halo=HALO):
    """
    6th-order centred first derivative of an array that already carries
    `halo` ghost planes on each side along `axis`.  Returns the derivative on
    the core (ghost planes removed).
    """
    n_core = ext.shape[axis] - 2 * halo

    def cut(a, b):
        sl = [slice(None)] * ext.ndim
        sl[axis] = slice(a, b)
        return ext[tuple(sl)]

    out = 0.0
    for m, cm in enumerate(FD6, start=1):
        out = out + cm * (cut(halo + m, halo + m + n_core)
                          - cut(halo - m, halo - m + n_core))
    return out / d


def fd6_derivative(f, d, axis, bc, par=1):
    """
    6th-order centred first derivative of a full array along `axis`.
    periodic : wrap-around (same operator as fd6_first_derivative);
    freeslip : mirror ghost points about the boundary nodes, with parity
               `par` (+1 even, -1 odd), see parity().
    """
    if bc == "periodic":
        return fd6_first_derivative(f, d, axis)
    n = f.shape[axis]
    idx, refl = ghost_indices(0, n, n, bc)          # planes -HALO .. n+HALO-1
    ext = np.take(f, idx, axis=axis)
    if par == -1:
        sl = [slice(None)] * f.ndim
        sl[axis] = refl
        ext[tuple(sl)] *= -1.0
    return fd6_from_extended(ext, d, axis)


class Slab:
    """
    Planes j0 <= j < j1 of a FieldFile (all x, all z), plus the ghost planes
    needed by the y-derivative.

        slab["ux"]        -> core values, shape (nx, j1-j0, nz)
        slab.grad("ux",1) -> d(ux)/dy on the core, computed on first request
        slab.vel_grad(i,j)-> d u_i / d x_j

    Everything is lazy: only the fields/gradients you ask for are read or
    computed.  Gradients are cached for the lifetime of the slab (use
    `slab.drop(name, axis)` or just let the slab go out of scope).
    """

    def __init__(self, ff, j0, j1):
        self.ff = ff
        self.j0, self.j1 = j0, j1
        g = ff.grid
        self._jm, self._refl = ghost_indices(j0, j1, g.n[1], g.bc[1])
        self._ext = {}
        self._grad = {}

    @property
    def y(self):
        return self.ff.grid.coords[1][self.j0:self.j1]

    @property
    def nplanes(self):
        return self.j1 - self.j0

    def _extended(self, name):
        if name not in self._ext:
            ext = np.asarray(self.ff.array(name)[:, self._jm, :],
                             dtype=np.float64)
            g = self.ff.grid
            if g.bc[1] == "freeslip" and parity(name, 1) == -1:
                ext[:, self._refl, :] *= -1.0
            self._ext[name] = ext
        return self._ext[name]

    def __getitem__(self, name):
        return self._extended(name)[:, HALO:HALO + self.nplanes, :]

    def grad(self, name, axis):
        """d(name)/d(x_axis) on the core planes (6th-order centred)."""
        key = (name, axis)
        if key not in self._grad:
            d = self.ff.grid.d[axis]
            if axis == 1:
                val = fd6_from_extended(self._extended(name), d, axis=1)
            else:                      # x, z: local to the slab
                val = fd6_derivative(self[name], d, axis,
                                     self.ff.grid.bc[axis],
                                     parity(name, axis))
            self._grad[key] = val
        return self._grad[key]

    def vel_grad(self, i, j):
        """d u_i / d x_j."""
        return self.grad(VELOCITY[i], j)

    def drop(self, name, axis):
        self._grad.pop((name, axis), None)

    def clear(self):
        """Forget all cached gradients (keeps the loaded fields)."""
        self._grad.clear()


# ----------------------------------------------------------------------------
# Diagnostics shared by several cases
# ----------------------------------------------------------------------------
def periodicity_check(ff):
    """max |u(last)-u(first)| per field and direction, with the neighbour
    difference |u(2)-u(1)| as a reference scale.
    wrap << neighbour  -> duplicated end point (use --drop-last)
    wrap ~  neighbour  -> no duplicate, all n points are unique."""
    rows = []
    for name in BASE_FIELDS:
        u = ff.raw(name)
        wrap = (np.abs(u[-1, :, :] - u[0, :, :]).max(),
                np.abs(u[:, -1, :] - u[:, 0, :]).max(),
                np.abs(u[:, :, -1] - u[:, :, 0]).max())
        nb   = (np.abs(u[1, :, :] - u[0, :, :]).max(),
                np.abs(u[:, 1, :] - u[:, 0, :]).max(),
                np.abs(u[:, :, 1] - u[:, :, 0]).max())
        fmt = lambda vals: "  ".join(
            f"{v:.2e}" if b == "periodic" else f"{'n/a':>8}"
            for v, b in zip(vals, ff.grid.bc))
        rows.append((f"|{name}(n)-{name}(1)| x,y,z", fmt(wrap)))
        rows.append((f"|{name}(2)-{name}(1)| x,y,z (ref.)", fmt(nb)))
    return rows

# ----------------------------------------------------------------------------
# Case interface (used by the dispatcher in fields_analyser.py)
# ----------------------------------------------------------------------------
class Case:
    """
    Interface of a case analysis.

    name         : key used on the command line (--<name>)
    bc           : default boundary conditions (x, y, z)
    implemented  : False -> the dispatcher only reports it is pending
    analyse(ff, args) -> dict : diagnostics of ONE snapshot (ff = FieldFile;
                         args.nu is the viscosity 1/Re)
    report(s, args)   : print the diagnostics of one snapshot
    finalize(history, args) : time series, files, plots (after all snapshots)
    """
    name = "case"
    bc = ("periodic", "periodic", "periodic")
    implemented = True

    def analyse(self, ff, args):
        raise NotImplementedError

    def report(self, s, args):
        pass

    def finalize(self, history, args):
        pass


class PendingCase(Case):
    """Placeholder for a case that is registered but not implemented yet."""
    implemented = False

    def __init__(self, name, bc=("periodic", "periodic", "periodic")):
        self.name = name
        self.bc = tuple(bc)


def get_filenames(filenames=None, directory="."):
    """Return the sorted list of field files to process (explicit list, or
    every fields_XXXXXX.bin found in `directory`)."""
    if filenames is not None:
        return sorted(directory / Path(f) for f in filenames)

    pattern = re.compile(r"^fields_(\d{6})\.bin$")
    matched = []
    for path in Path(directory).iterdir():
        m = pattern.match(path.name)
        if m:
            matched.append((int(m.group(1)), path))
    matched.sort(key=lambda item: item[0])
    return [path for _, path in matched]


# ----------------------------------------------------------------------------
# Numerical tools
# ----------------------------------------------------------------------------
def fd6_first_derivative(f, dx, axis):
    """6th-order centered first derivative, periodic (same as derx_00)."""
    c = (45.0 / 60.0, -9.0 / 60.0, 1.0 / 60.0)
    out = np.zeros_like(f)
    for m, cm in enumerate(c, start=1):
        out += cm * (np.roll(f, -m, axis) - np.roll(f, m, axis))
    return out / dx

def spectral_first_derivative(f, L, axis):
    """Exact (spectral) first derivative on a periodic box."""
    n = f.shape[axis]
    k = 2.0 * np.pi * np.fft.fftfreq(n, d=L / n)
    shape = [1, 1, 1]
    shape[axis] = n
    k = k.reshape(shape)
    return np.real(np.fft.ifftn(1j * k * np.fft.fftn(f)))
