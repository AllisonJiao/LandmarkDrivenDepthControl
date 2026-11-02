"""Runs the app's C++ solver (src/tools/landmark_cli.cpp, built by
build-web.sh into build/web/cli/) on a base scene and a handle set, and
returns its results as the Result / Fields of results.py, so the
plotting and metrics code works unchanged.

The C++ is the single implementation of the methods: whatever the app does
(including Part A's later changes, e.g. interior-row P3 or P3-R) is what the
figures show. Python only resolves handles, runs the CLI and converts units.

Environment overrides: MM_NODE (node executable), MM_CLI (path to
landmark_cli.js).

  python cpp_bridge.py --scene data/ellipse_nb120.json --handles configs/ellipse_x1.5.json
"""

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from common import ROOT, SOLVE_KEYS, load_base, read_json, write_json
from results import METHODS, Fields, Result
from handles import handle_depths, handle_factors, load_handles, resolve_vertices

CLI = Path(os.environ.get("MM_CLI", ROOT / "build" / "web" / "cli" / "landmark_cli.js"))
NODE = os.environ.get("MM_NODE", "node")


class CppSolverError(RuntimeError):
    pass


def check_cli():
    if not CLI.exists():
        raise CppSolverError(f"{CLI} not found: run ./build-web.sh (it builds landmark_cli)")
    if shutil.which(NODE) is None and not Path(NODE).exists():
        raise CppSolverError(f"node not found ({NODE}); set MM_NODE, e.g. to emsdk's node.exe")


def fields_from_export(out):
    """Fields of an exported solve (app units) in normalized units: lengths
    / R, u-like fields / R^2, pressure unchanged. dp on the outline is NaN,
    as in the app's heat map."""
    f, R = out.fields, out.R
    du = f["du"] / R**2
    dp_rel = f["dp_out"] / out.c
    dp_rel = np.where(out.is_fixed, np.nan, dp_rel)
    max_h0 = np.abs(out.V_input[:, 2]).max()
    return Fields(
        u=out.u0 + du, h=f["h"] / R, du=du, dp=f["dp"],
        du_out=f["du_out"] / R**2, dp_out=f["dp_out"], dh_out=f["dh_out"] / R,
        dh_rel=f["dh_out"] / max_h0, dp_rel=dp_rel)


def run(scene_path, hs, methods=METHODS, R_edges=4.0, weight=None, depth_override=None,
        keep_dir=None):
    """Solve one handle set under several methods with the C++ solver.

    Returns (out_scene, runs): out_scene is the first exported result loaded
    as a Scene; it carries the C++ operators (L, Voronoi M, component), so
    use it for operators() and metrics. runs is
    [(method, Result or None, Fields or None, error or None)].
    keep_dir: keep the CLI's input and outputs there instead of a temp dir."""
    check_cli()
    base = load_base(scene_path)
    vertices = resolve_vertices(base, hs)
    factors = handle_factors(base, vertices, handle_depths(hs, depth_override), hs.depth_unit)

    scene_json = {k: v for k, v in read_json(scene_path).items() if k not in SOLVE_KEYS}
    scene_json["landmarks"] = [
        {"vertex": int(v), "factor": float(f), "anchor": bool(h.anchor)}
        for v, f, h in zip(vertices, factors, hs.handles)]

    work = Path(keep_dir) if keep_dir else Path(tempfile.mkdtemp(prefix="mm_cli_"))
    try:
        write_json(work / "input.json", scene_json)
        cmd = [NODE, str(CLI), "--scene", str(work / "input.json"), "--out-dir", str(work),
               "--prefix", "run", "--methods", ",".join(methods), "--radius", str(R_edges)]
        if weight is not None:
            cmd += ["--weight", str(weight)]
        # input="": give the CLI an empty pipe as stdin. Emscripten's NODERAWFS
        # can't stat the Windows NUL device, so an inherited NUL stdin (e.g.
        # from a non-interactive shell, or subprocess.DEVNULL) makes node die
        # with an unhandled ErrnoError before main() runs.
        proc = subprocess.run(cmd, capture_output=True, text=True, input="")
        if proc.returncode not in (0, 1):  # 1: some method failed, files still written
            raise CppSolverError(f"landmark_cli failed ({proc.returncode}):\n{proc.stderr}")

        missing = [m for m in methods if not (work / f"run_{m}.json").exists()]
        if missing:
            raise CppSolverError(
                f"landmark_cli wrote no result for {', '.join(missing)} in {work} "
                f"(exit {proc.returncode}).\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")

        out_scene, runs = None, []
        for method in methods:
            out = load_base(work / f"run_{method}.json")
            if out_scene is None:
                out_scene = out
            if not out.ok:
                runs.append((method, None, None, "C++ solve failed"))
                continue
            f = fields_from_export(out)
            # targets in u, normalized: (f^2 - 1) u0 at the handles
            r = (factors**2 - 1) * out.u0[vertices]
            settings = {"R_edges": out.radius, "R_norm": (out.radius_px or 0) / out.R,
                        "weight": out.weight}
            res = Result(method=method, du=f.du, vertices=vertices, r=r, settings=settings)
            runs.append((method, res, f, None))
        return out_scene, runs
    finally:
        if not keep_dir:
            shutil.rmtree(work, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description="run the C++ solver on a scene + handle file")
    ap.add_argument("--scene", required=True)
    ap.add_argument("--handles", required=True)
    ap.add_argument("--methods", default=",".join(METHODS))
    ap.add_argument("--R", type=float, default=4.0)
    ap.add_argument("--weight", type=float, default=None)
    ap.add_argument("--keep", default=None, help="keep the CLI files in this folder")
    args = ap.parse_args()

    methods = [m for m in args.methods.split(",") if m]
    out, runs = run(args.scene, load_handles(args.handles), methods, args.R, args.weight,
                    keep_dir=args.keep)
    for method, res, f, err in runs:
        if f is None:
            print(f"{method}: {err}")
            continue
        v = res.vertices
        ratio = np.abs(f.h[v]) / np.abs(out.h0[v])
        print(f"{method}: |h|/|h0| at handles {np.round(ratio, 4)}"
              f"  peak |dp|/c {np.nanmax(np.abs(f.dp_rel)):.3g}"
              f"  max |dh_rel| {np.abs(f.dh_rel).max():.3g}")


if __name__ == "__main__":
    main()
