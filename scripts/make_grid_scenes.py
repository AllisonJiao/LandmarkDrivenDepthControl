"""Base scenes for the grid tests: the synthetic two-sheet ellipse of
sds/make_test_scene.py at several resolutions, without landmarks (handles come
from handle files, e.g. configs/ellipse_x1.5.json).

  python make_grid_scenes.py                     # nb = 60 120 240 480
  python make_grid_scenes.py --nb 60 120 --out data

Each file is data/ellipse_nb<nb>.json. nb is the number of outline vertices;
the interior spacing follows it, so doubling nb halves the edge length.
In these scenes the first sheet has h0 < 0, so its handles use "side": -1.

make_test_scene.py solves the default inflation with a lumped (barycentric)
mass matrix, while the C++ solver (LandmarkSolver::prepare) recomputes a
Voronoi one, so on the raw scene p0 = M^-1 L u0 is not a uniform +-c in the
C++. By default the heights are re-solved here with libigl's Voronoi mass, so
the C++ sees exactly the default pressure +-c (--mass lumped keeps the raw
scene).
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg

from common import ROOT, SCRIPTS, read_json, write_json

MAKE = ROOT / "sds" / "make_test_scene.py"
C = 2.0  # default inflation amount, as in the app and make_test_scene.py


def voronoi_mass(V, F):
    """Diagonal of libigl's MASSMATRIX_TYPE_VORONOI on the flat mesh (x, y):
    the circumcentric (Voronoi) area per corner, and for an obtuse triangle
    half its area to the obtuse corner and a quarter to the others
    (igl::massmatrix_intrinsic)."""
    xy = V[:, :2]
    P = [xy[F[:, k]] for k in range(3)]
    # l[:, k]: length of the edge opposite corner k
    l = np.stack([np.linalg.norm(P[1] - P[2], axis=1), np.linalg.norm(P[2] - P[0], axis=1),
                  np.linalg.norm(P[0] - P[1], axis=1)], 1)
    d1, d2 = P[1] - P[0], P[2] - P[0]
    area = 0.5 * np.abs(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0])
    l2 = l**2
    cos = np.stack([(l2[:, 1] + l2[:, 2] - l2[:, 0]) / (2 * l[:, 1] * l[:, 2]),
                    (l2[:, 2] + l2[:, 0] - l2[:, 1]) / (2 * l[:, 2] * l[:, 0]),
                    (l2[:, 0] + l2[:, 1] - l2[:, 2]) / (2 * l[:, 0] * l[:, 1])], 1)
    bary = cos * l  # circumcenter barycentric coordinates (unnormalized)
    bary /= bary.sum(1, keepdims=True)
    partial = bary * area[:, None]
    quads = 0.5 * np.stack([partial[:, 1] + partial[:, 2], partial[:, 2] + partial[:, 0],
                            partial[:, 0] + partial[:, 1]], 1)
    for k in range(3):  # obtuse at corner k
        obtuse = cos[:, k] < 0
        quads[obtuse] = 0.25 * area[obtuse, None]
        quads[obtuse, k] = 0.5 * area[obtuse]
    m = np.zeros(len(V))
    for k in range(3):
        np.add.at(m, F[:, k], quads[:, k])
    return m


def resolve_with_voronoi(path):
    """Re-solve the default heights of a make_test_scene.py file with the
    Voronoi mass: L_II u = M_II p with p = -C sgn(u) on each sheet (the sign
    make_test_scene.py uses), u = 0 on the outline."""
    d = read_json(path)
    V = np.asarray(d["V"], dtype=np.float64)
    F = np.asarray(d["F"], dtype=np.int64)
    n = len(V)
    L = sp.csr_matrix((d["L"]["v"], (d["L"]["i"], d["L"]["j"])), shape=(n, n))
    fixed = np.zeros(n, dtype=bool)
    fixed[np.asarray(d["fixed"], dtype=np.int64)] = True
    I = np.nonzero(~fixed)[0]

    m = voronoi_mass(V, F)
    p = -C * np.sign(V[I, 2])  # +C on the h0 < 0 sheet, as in make_test_scene.py
    u = np.zeros(n)
    u[I] = scipy.sparse.linalg.spsolve(L[I][:, I].tocsc(), m[I] * p)
    V[:, 2] = np.sign(u) * np.sqrt(np.abs(u))
    d["V"], d["M"] = V.tolist(), m.tolist()
    write_json(path, d)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--nb", type=int, nargs="+", default=[60, 120, 240, 480])
    ap.add_argument("--a", type=float, default=150.0, help="semi-axis x (px)")
    ap.add_argument("--b", type=float, default=100.0, help="semi-axis y (px)")
    ap.add_argument("--mass", default="voronoi", choices=("voronoi", "lumped"),
                    help="mass matrix of the default inflation (voronoi matches the C++)")
    ap.add_argument("--out", default=str(SCRIPTS / "data"))
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for nb in args.nb:
        path = out / f"ellipse_nb{nb}.json"
        subprocess.run([sys.executable, str(MAKE), "--out", str(path), "--nb", str(nb),
                        "--a", str(args.a), "--b", str(args.b)], check=True)
        if args.mass == "voronoi":
            resolve_with_voronoi(path)


if __name__ == "__main__":
    main()
