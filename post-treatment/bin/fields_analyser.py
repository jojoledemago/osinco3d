#!/usr/bin/env python3
"""
Post-processing of osinco3d flow fields (fields_XXXXXX.bin).

Layout
------
fields_core.py          reader (optional phi), grid/BC, slabs, gradients on
                        demand, Case interface
case_hit.py             --hit            decaying homogeneous isotropic turb.
case_mixing_layer.py    --mixing_layer   temporal mixing layer
les_sgs.py              --les [--cs X]   Smagorinsky SGS diagnostics, any case
fields_analyser.py      command line + dispatch (this file)

Directory layout (example)
--------------------------
post-treatment/
    bin/            <- all the .py files (this script is run from here or
                       from anywhere: imports resolve next to the script)
    hit_case/       <- fields_XXXXXX.bin (or a link to them); results/ is
    ml_case/           created here
    ...

Typical use, from a case directory:

    cd post-treatment/hit_case
    python ../bin/fields_analyser.py --hit --reynolds 1600 --plot
    python ../bin/fields_analyser.py --hit --input-dir /path/to/ns_run/out

To add a case: write a Case subclass (analyse / report / finalize) in its own
module and register it in CASES below.
"""

import argparse
from pathlib import Path

from fields_core import (BC_CHOICES, FieldFile, PendingCase, AXES,
                         format_value, get_filenames, periodicity_check,
                         print_block)
from case_hit import HitCase
from case_tgv import TgvCase
from les_sgs import sgs_analysis, sgs_report, sgs_finalize
from case_mixing_layer import MixingLayerCase

CASES = {
    "mixing_layer": MixingLayerCase(),
    "coplanar": PendingCase("coplanar"),
    "planar": PendingCase("planar"),
    "tgv": TgvCase(),
    "hit": HitCase(),
}


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Post-process instantaneous CFD flow fields."
    )
    parser.add_argument("--filenames", nargs="+", type=Path,
                        help="List of field files to process.")
    parser.add_argument("--input-dir", type=Path, default=Path("."),
                        help="Directory searched for fields_XXXXXX.bin when "
                             "--filenames is not given (default: current "
                             "directory).")
    parser.add_argument("--outdir", type=Path, default=Path("results"),
                        help="Directory for the output files (time series, "
                             "figures); created if needed (default: "
                             "./results).")
    parser.add_argument("--plot", action="store_true",
                        help="Display and save the analysis plots.")
    parser.add_argument("--les", action="store_true",
                        help="Smagorinsky SGS diagnostics (nu_t, SGS "
                             "dissipation, LES quality, energy budget); "
                             "works with any case, or alone.")
    parser.add_argument("--cs", type=float, default=0.17,
                        help="Smagorinsky coefficient used by --les "
                             "(default: 0.17)")
    parser.add_argument("--les-delta", default="cbrt",
                        help="Filter width for --les: 'cbrt' = (dx dy dz)^1/3 "
                             "(default), 'max' = max(dx,dy,dz), or a number.")
    parser.add_argument("--turbulence", action="store_true",
                        help="Perform turbulence-scale analysis")
    parser.add_argument("--reynolds", type=float, default=1000.0,
                        help="Initial Reynolds number, nu = 1/Re "
                             "(default: 1000)")
    parser.add_argument("--drop-last", action="store_true",
                        help="Drop the last point of each periodic direction "
                             "(legacy files with duplicated periodic point).")
    parser.add_argument("--check-periodicity", action="store_true",
                        help="Print |u(n)-u(1)| for each direction.")
    parser.add_argument("--quiet", action="store_true",
                        help="Print only a one-line summary per file.")
    parser.add_argument("--bc", nargs=3, choices=BC_CHOICES,
                        metavar=("BCX", "BCY", "BCZ"),
                        help="Boundary conditions in x, y, z "
                             f"({'|'.join(BC_CHOICES)}). Default: the one of "
                             "the selected case.")
    parser.add_argument("--slab", type=int, default=16,
                        help="Number of y-planes processed at a time by "
                             "slab-based analyses (default: 16).")
    parser.add_argument("--kh-modes", type=int, default=16,
                        help="Mixing layer: number of streamwise Fourier "
                             "modes m (kx = 2 pi m/Lx) analysed (default: 16).")
    parser.add_argument("--kh-mode", type=int, default=0,
                        help="Mixing layer: mode m whose growth rate is "
                             "followed in time (default: 0 = most unstable "
                             "Rayleigh mode of the first snapshot).")
    parser.add_argument("--rayleigh-npts", type=int, default=401,
                        help="Mixing layer: max points used by the Rayleigh "
                             "stability solver (the profile is only "
                             "re-sampled if it has more points); <= 0 "
                             "disables it (default: 401).")
    parser.add_argument("--ml-zoom", type=float, default=5.0,
                        help="Mixing layer plots: show only y_c +/- ZOOM*d_w "
                             "(default: 5; 0 = whole domain).")
    parser.add_argument("--no-show", action="store_true",
                        help="With --plot: save the figures but do not open "
                             "windows.")

    case_group = parser.add_mutually_exclusive_group()
    for name in CASES:
        case_group.add_argument(f"--{name}", action="store_true",
                                help=f"Select the {name} case.")
    return parser.parse_args()


def selected_case(args):
    for name, case in CASES.items():
        if getattr(args, name):
            return case
    return None


def main():
    args = parse_arguments()
    args.nu = 1.0 / args.reynolds
    case = selected_case(args)
    bc = tuple(args.bc) if args.bc else (case.bc if case else
                                         ("periodic",) * 3)

    filenames = get_filenames(args.filenames, args.input_dir)
    if not filenames:
        raise FileNotFoundError("No field files found.")
    print(f"Found {len(filenames)} field file(s).")

    if args.les:
        print(f"SGS diagnostics: Smagorinsky, Cs = {args.cs}, "
              f"Delta = {args.les_delta}.")
        if case is not None and case.implemented:
            print("Note: the energy budgets and the Kolmogorov-scale / "
                  "kmax*eta diagnostics of the case analysis only account "
                  "for the molecular dissipation; for LES use the SGS "
                  "analysis (energy budget, Delta/eta, LES_IQ).")

    history = []
    sgs_history = []
    for filename in filenames:
        print(f"\n Reading: {filename}")
        ff = FieldFile(filename, bc=bc, drop_last=args.drop_last)
        g = ff.grid
        print_block(
            "Field file",
            [
                ("Filename", ff.path),
                ("Iteration", f"{ff.iteration:06d}"),
                ("Simulation time", ff.time),
                ("Grid size", " x ".join(str(n) for n in g.n_file)),
                ("Domain size (x[-1]-x[0])",
                 " x ".join(format_value(e) for e in g.extent)),
                ("Length used (n*d if periodic)",
                 " x ".join(format_value(e) for e in g.length)),
                ("BC (x, y, z)", "  ".join(g.bc)),
                ("Fields", " ".join(ff.fields)),
            ],
        )

        if args.check_periodicity:
            print_block("Periodicity check", periodicity_check(ff))

        if case is not None:
            if case.implemented:
                s = case.analyse(ff, args)
                history.append(s)
                case.report(s, args)
            else:
                print("  This case analysis is not implemented yet.")

        if args.les:
            s_sgs = sgs_analysis(ff, args)
            sgs_history.append(s_sgs)
            sgs_report(s_sgs, args)

    if case is not None and case.implemented and history:
        args.outdir.mkdir(parents=True, exist_ok=True)
        case.finalize(history, args)
    if args.les and sgs_history:
        sgs_finalize(sgs_history, args)


if __name__ == "__main__":
    main()
