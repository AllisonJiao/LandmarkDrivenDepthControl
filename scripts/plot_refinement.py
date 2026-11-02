"""Refinement study (part_b_plan.md B5.2): the same position-based handle set
on base scenes of increasing resolution, metrics against the mean edge
length on log-log axes.

  python make_grid_scenes.py                   # writes data/ellipse_nb*.json
  python plot_refinement.py --scenes data/ellipse_nb*.json \
      --handles configs/ellipse_x1.5.json --out out/refinement
  options: --methods G1,P1,P3  --dashed G2,P2  --R-px 27  (or --R edges)

Use a handle file with "at" handles: vertex ids change with the resolution.
The falloff radius of G2/P2 is held fixed in canvas pixels (--R-px) and
converted to mean edge lengths per mesh, so the disk keeps its size while the
mesh is refined; --R (in edge lengths) instead shrinks it with the mesh.
"""

import argparse
import csv
import glob
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402

from common import length_to_px, load_base  # noqa: E402
from handles import load_handles  # noqa: E402
from metrics import all_metrics  # noqa: E402
from plot_handles import run_methods  # noqa: E402
from results import mean_edge_length  # noqa: E402

# (metric, axis label, how to reduce it)
PANELS = (
    ("peak_dp", r"peak $|dp_{out}|/c$"),
    ("dp_energy", r"$dp^\top M\,dp$"),
    ("area_moved", "area moved (fraction)"),
)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scenes", nargs="+", required=True, help="base scenes (globs allowed)")
    ap.add_argument("--handles", required=True, help="handle file with 'at' handles")
    ap.add_argument("--methods", default="G1,P1,P3")
    ap.add_argument("--dashed", default="G2,P2", help="methods drawn dashed (fixed R)")
    R = ap.add_mutually_exclusive_group()
    R.add_argument("--R-px", type=float, default=27.0,
                   help="falloff radius in canvas pixels, the same on every mesh (default)")
    R.add_argument("--R", type=float, default=None,
                   help="falloff radius in mean edge lengths of each mesh instead")
    ap.add_argument("--out", default="out/refinement")
    args = ap.parse_args()

    # natural order (nb60 before nb120): shorter names first, then by name
    paths = sorted({p for pat in args.scenes for p in glob.glob(pat)},
                   key=lambda p: (len(p), p))
    if not paths:
        raise SystemExit("no scenes matched")
    solid = [m for m in args.methods.split(",") if m]
    dashed = [m for m in args.dashed.split(",") if m]
    methods = solid + [m for m in dashed if m not in solid]
    hs = load_handles(args.handles)
    if any(h.vertex is not None for h in hs.handles):
        print("warning: vertex handles only work on their own mesh; use 'at' handles")

    table = []  # (scene, n, mean_edge_px, method, metrics)
    for path in paths:
        base = load_base(path)
        edge_px = float(length_to_px(base, mean_edge_length(base)))
        R_edges = args.R if args.R is not None else args.R_px / edge_px
        scene, ops, runs = run_methods(path, hs, methods, R_edges, None)
        h_px = float(length_to_px(scene, ops.mean_edge))
        for method, res, f, err in runs:
            m = all_metrics(scene, ops, res, f) if f is not None else {}
            if m:
                m["peak_dp"] = max(abs(m.get("dp_max", np.nan)), abs(m.get("dp_min", np.nan)))
            table.append((path, scene.n, h_px, method, m))
        print(f"{path}: n = {scene.n}, mean edge {h_px:.2f} px, R = {R_edges:.2f} edges")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    keys = sorted({k for *_, m in table for k in m})
    with open(out / "refinement.csv", "w", newline="") as f:
        # TODO(B5): update when Part A lands interior-row P3 (Delta 1) in C++
        R_desc = f"{args.R} edges" if args.R is not None else f"{args.R_px} px"
        f.write(f"# handle_file: {args.handles}\n# R: {R_desc}\n# p3: C++ (every row)\n")
        w = csv.writer(f)
        w.writerow(["scene", "n", "mean_edge_px", "method", *keys])
        for path, n, h_px, method, m in table:
            w.writerow([path, n, h_px, method, *(m.get(k, float("nan")) for k in keys)])

    fig, axes = plt.subplots(1, len(PANELS), figsize=(4.2 * len(PANELS), 3.4),
                             constrained_layout=True)
    for ax, (key, label) in zip(axes, PANELS):
        for method in methods:
            pts = [(h_px, m.get(key, np.nan)) for _, _, h_px, mm, m in table if mm == method]
            pts = [(x, y) for x, y in pts if np.isfinite(y) and y > 0]
            if not pts:
                continue
            x, y = np.array(sorted(pts)).T
            ax.loglog(x, y, "--o" if method in dashed and method not in solid else "-o",
                      ms=4, label=method)
        ax.xaxis.set_minor_formatter(NullFormatter())  # minor labels overlap
        ax.set_xlabel("mean edge length (px)")
        ax.set_title(label, fontsize=10)
        ax.invert_xaxis()  # finer meshes to the right
        if ax.lines:
            ax.legend(fontsize=8)
        else:
            ax.text(0.5, 0.5, "no data yet", ha="center", va="center",
                    transform=ax.transAxes, color="0.4")
    fig.savefig(out / "refinement.png", dpi=200)
    plt.close(fig)
    print(f"wrote {out / 'refinement.csv'} and {out / 'refinement.png'}")


if __name__ == "__main__":
    main()
