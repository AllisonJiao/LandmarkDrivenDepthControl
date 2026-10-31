"""Synthetic stand-in for an exported Monster Mash scene: an inflated ellipse
made of a front and a back sheet that share the outline, with the default
inflation solved exactly as in reconstruction.cpp (M^-1 L z = c, z = 0 on
the outline, c = 2 on the front and -2 on the back).

  python make_test_scene.py --out test_ellipse.json [--landmark X Y FACTOR ...]
"""

import argparse
import json

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg
from scipy.spatial import Delaunay


def cotan_laplacian(V, F):
    n = len(V)
    rows, cols, vals = [], [], []
    for a, b, c in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        # angle at c, opposite edge (a, b)
        u = V[F[:, a]] - V[F[:, c]]
        w = V[F[:, b]] - V[F[:, c]]
        cos = (u * w).sum(1)
        sin = np.linalg.norm(np.cross(u, w), axis=1)
        half_cot = 0.5 * cos / sin
        for i, j in ((a, b), (b, a)):
            rows.append(F[:, i]); cols.append(F[:, j]); vals.append(half_cot)
        for i in (a, b):
            rows.append(F[:, i]); cols.append(F[:, i]); vals.append(-half_cot)
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(n, n))


def lumped_mass(V, F):
    area = 0.5 * np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1)
    M = np.zeros(len(V))
    for j in range(3):
        np.add.at(M, F[:, j], area / 3)
    return M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="test_ellipse.json")
    ap.add_argument("--a", type=float, default=150.0)
    ap.add_argument("--b", type=float, default=100.0)
    ap.add_argument("--nb", type=int, default=120, help="outline vertices")
    ap.add_argument("--landmark", type=float, nargs=3, action="append", default=[],
                    metavar=("X", "Y", "FACTOR"), help="front landmark at (X, Y) in [-1, 1]")
    args = ap.parse_args()

    a, b = args.a, args.b
    t = np.linspace(0, 2 * np.pi, args.nb, endpoint=False)
    outline = np.stack([a * np.cos(t), b * np.sin(t)], 1)
    step = 2 * np.pi * np.sqrt((a * a + b * b) / 2) / args.nb
    gx, gy = np.meshgrid(np.arange(-a, a, step), np.arange(-b, b, step * np.sqrt(3) / 2))
    gx = gx + (np.arange(gx.shape[0])[:, None] % 2) * step / 2
    grid = np.stack([gx.ravel(), gy.ravel()], 1)
    inside = (grid[:, 0] / (a - step / 2)) ** 2 + (grid[:, 1] / (b - step / 2)) ** 2 < 1
    interior = grid[inside]

    pts = np.concatenate([outline, interior])
    tri = Delaunay(pts).simplices
    nb, ni = len(outline), len(interior)

    # back sheet: duplicate interior vertices, share the outline, flip winding
    back = np.arange(nb + ni)
    back[nb:] += ni
    F = np.concatenate([tri, back[tri][:, ::-1]])
    V = np.zeros((nb + 2 * ni, 3))
    V[: nb + ni, :2] = pts
    V[nb + ni:, :2] = interior

    L = cotan_laplacian(V, F)
    M = lumped_mass(V, F)
    fixed = np.arange(nb)
    I = np.arange(nb, len(V))
    c = np.where(I < nb + ni, 2.0, -2.0)
    z = np.zeros(len(V))
    z[I] = scipy.sparse.linalg.spsolve(L[I][:, I].tocsc(), M[I] * c)
    V[:, 2] = np.sign(z) * np.sqrt(np.abs(z))

    landmarks = []
    for x, y, factor in args.landmark:
        p = np.array([x * a, y * b])
        for f in range(len(tri)):
            A, B, C = pts[tri[f]]
            T = np.array([B - A, C - A]).T
            w1, w2 = np.linalg.solve(T, p - A)
            if w1 >= 0 and w2 >= 0 and w1 + w2 <= 1:
                landmarks.append({"face": f, "bary": [1 - w1 - w2, w1, w2],
                                  "factor": factor, "side": -1, "anchor": False})
                break

    Lc = L.tocoo()
    scene = {
        "version": 1,
        "V": V.tolist(), "F": F.tolist(), "fixed": fixed.tolist(),
        "L": {"i": Lc.row.tolist(), "j": Lc.col.tolist(), "v": Lc.data.tolist()},
        "M": M.tolist(), "landmarks": landmarks,
    }
    with open(args.out, "w") as f:
        json.dump(scene, f)
    print(f"wrote {args.out}: {len(V)} vertices, {len(F)} faces, {len(landmarks)} landmarks")


if __name__ == "__main__":
    main()
