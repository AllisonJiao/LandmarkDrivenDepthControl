"""Method comparison on one base scene and one handle set (part_b_plan.md B5.2):
the dh/dp grid, the height profile and the metrics table, for every method
(and every sweep depth, if the handle file has a sweep).

  python plot_handles.py --scene data/ellipse_nb120.json \
      --handles configs/ellipse_x1.5.json --out out/ellipse
  options: --methods G1,G2,P1,P2,P3  --R 4  --weight W  --symlog
           --no-wireframe  --name NAME

The methods are solved by the app's C++ solver (cpp_bridge.py; build it with
./build-web.sh). A method that fails shows as an empty panel and empty cells.
"""

import argparse
from pathlib import Path

import numpy as np

import cpp_bridge
from common import length_to_px, load_base
from handles import load_handles, resolve_vertices
from metrics import METRIC_NAMES, all_metrics
from plotting import plot_method_grid, plot_profile, write_metrics
from results import METHODS, operators


def run_methods(scene_path, hs, methods, R_edges, weight, depth_override=None):
    """Solve one handle set under every method with the C++ solver.

    Returns (scene, ops, runs) with runs = [(method, Result or None,
    Fields or None, error message or None)]. scene and ops come from the C++
    output (its L, Voronoi M and sheet labels), so metrics and plots use the
    same operators as the solve."""
    out_scene, runs = cpp_bridge.run(scene_path, hs, methods, R_edges, weight, depth_override)
    return out_scene, operators(out_scene), runs


def figures(scene, ops, hs, runs, out_dir, name, header, symlog, wireframe):
    vertices = resolve_vertices(scene, hs)
    anchors = [h.anchor for h in hs.handles]
    # the figures show the sheet of the first moving handle
    moving = [v for v, a in zip(vertices, anchors) if not a] or list(vertices)
    focus = int(moving[0])
    s = int(ops.sheet[focus])

    columns = [(m, f, vertices, anchors) for m, _, f, _ in runs]
    title = f"{name}  ({header['handles']}, R = {header['R_edges']} edges)"
    plot_method_grid(scene, ops.sheet, s, columns, out_dir / f"{name}_grid.png",
                     symlog=symlog, wireframe=wireframe, title=title)

    curves = [(m, length_to_px(scene, f.h)) for m, _, f, _ in runs if f is not None]
    plot_profile(scene, ops.sheet, s, length_to_px(scene, scene.h0), curves, focus,
                 out_dir / f"{name}_profile.png", title=name)

    rows = [(m, all_metrics(scene, ops, res, f) if f is not None else {})
            for m, res, f, _ in runs]
    write_metrics(rows, header, out_dir / f"{name}_metrics.csv",
                  out_dir / f"{name}_metrics.md", METRIC_NAMES)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scene", required=True, help="base scene JSON")
    ap.add_argument("--handles", required=True, help="handle file JSON")
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--R", type=float, default=4.0, help="falloff radius in mean edge lengths")
    ap.add_argument("--weight", type=float, default=None, help="P3 soft weight (default hard)")
    ap.add_argument("--out", default="out")
    ap.add_argument("--name", default=None, help="file prefix (default: handle file name)")
    ap.add_argument("--symlog", action="store_true", help="symmetric-log colours for dp")
    ap.add_argument("--no-wireframe", action="store_true")
    args = ap.parse_args()

    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    scene, hs = load_base(args.scene), load_handles(args.handles)
    ops = operators(scene)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    name = args.name or Path(args.handles).stem

    base_header = {
        "scene": args.scene, "mesh_name": scene.mesh_name or "-", "n": scene.n,
        "handle_file": args.handles, "depth_unit": hs.depth_unit,
        "R_edges": args.R, "weight": "hard" if args.weight is None else args.weight,
        # TODO(B5): update when Part A lands interior-row P3 (Delta 1) in C++
        "p3": "C++ (every row)",
        "c": ops.c, "mean_edge_px": round(float(length_to_px(scene, ops.mean_edge)), 3),
    }

    if hs.sweep:
        k = int(hs.sweep["handle"])
        depth_runs = [(f"{name}_d{d:g}", (k, float(d))) for d in hs.sweep["depths"]]
    else:
        depth_runs = [(name, None)]

    for run_name, override in depth_runs:
        run_scene, run_ops, runs = run_methods(args.scene, hs, methods, args.R, args.weight,
                                               override)
        depths = [d for d in (hs.handles[i].depth for i in range(len(hs.handles)))]
        if override is not None:
            depths[override[0]] = override[1]
        header = dict(base_header, handles=", ".join(
            f"{'anchor' if h.anchor else 'handle'} {d if d is not None else 'hold'}"
            for h, d in zip(hs.handles, depths)))
        figures(run_scene, run_ops, hs, runs, out_dir, run_name, header,
                args.symlog, not args.no_wireframe)
        missing = [f"{m} ({err})" for m, res, _, err in runs if res is None]
        print(f"{run_name}: wrote grid, profile, metrics to {out_dir}"
              + (f"  (missing: {', '.join(missing)})" if missing else ""))


if __name__ == "__main__":
    main()
