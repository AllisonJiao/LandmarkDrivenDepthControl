"""Scene exported from Monster Mash ("Export for SDS"): the default rest mesh,
the exact flat operators used by the C++ landmark solver, and the landmarks.

Coordinates are normalized so that the sketch spans about [-1, 1]. Under the
scaling x -> x / R the 2D cotangent Laplacian is unchanged, the mass matrix
scales by 1 / R^2 and so does u = sgn(h) h^2, which leaves the pressure
p = M^-1 L u unchanged.
"""

import json
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp


@dataclass
class Landmark:
    vertex: int  # mesh vertex the handle sits on (never on the outline)
    factor: float  # target height = factor * default height
    side: int  # sign of the rest z of the sheet it lies on
    anchor: bool


class Scene:
    def __init__(self, path):
        with open(path) as f:
            d = json.load(f)

        V = np.asarray(d["V"], dtype=np.float64)
        self.F = np.asarray(d["F"], dtype=np.int64)
        n = len(V)
        self.n = n
        self.V_input = V.copy()

        # normalization (undone by to_original_heights)
        lo, hi = V[:, :2].min(0), V[:, :2].max(0)
        self.center = 0.5 * (lo + hi)
        self.R = 0.5 * float((hi - lo).max())
        self.V = V.copy()
        self.V[:, :2] = (V[:, :2] - self.center) / self.R
        self.V[:, 2] = V[:, 2] / self.R

        Ld = d["L"]
        self.L = sp.csr_matrix((Ld["v"], (Ld["i"], Ld["j"])), shape=(n, n))
        self.M = np.asarray(d["M"], dtype=np.float64) / self.R**2

        self.fixed = np.unique(np.asarray(d["fixed"], dtype=np.int64))
        is_fixed = np.zeros(n, dtype=bool)
        is_fixed[self.fixed] = True
        self.interior = np.nonzero(~is_fixed)[0]

        h0 = self.V[:, 2]
        self.h0 = h0
        self.u0 = np.sign(h0) * h0**2

        self.landmarks = []
        for l in d.get("landmarks", []):
            v = self._landmark_vertex(l, is_fixed)
            if v is None:
                continue
            self.landmarks.append(Landmark(
                vertex=v,
                factor=float(l["factor"]),
                side=int(l["side"]),
                anchor=bool(l.get("anchor", False)),
            ))

    def _landmark_vertex(self, l, is_fixed):
        """Vertex of an exported landmark. Scenes written before landmarks
        became vertices store a face and barycentric weights instead; those
        snap the same way as pickLandmark: the corner with the largest weight
        that is not on the outline (None if all three are)."""
        if "vertex" in l:
            return int(l["vertex"])
        corners = self.F[int(l["face"])]
        bary = np.asarray(l["bary"], dtype=np.float64)
        ok = ~is_fixed[corners]
        if not ok.any():
            return None
        return int(corners[ok][np.argmax(bary[ok])])

    def landmark_matrix(self):
        """Sparse k x n matrix B whose row k reads a field at landmark k's
        vertex, so (B u)[k] = u[vertex k]."""
        k = len(self.landmarks)
        cols = [lm.vertex for lm in self.landmarks]
        return sp.csr_matrix((np.ones(k), (np.arange(k), cols)), shape=(k, self.n))

    def to_original_heights(self, h):
        return np.asarray(h) * self.R

    def write_result(self, path, h):
        """Heights in original units, in the format "Import SDS result" reads."""
        result = {
            "version": 1,
            "n": self.n,
            "xy": self.V_input[:, :2].tolist(),
            "h": self.to_original_heights(h).tolist(),
        }
        with open(path, "w") as f:
            json.dump(result, f)

    def write_obj(self, path, h):
        V = self.V_input.copy()
        V[:, 2] = self.to_original_heights(h)
        with open(path, "w") as f:
            for v in V:
                f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
            for t in self.F + 1:
                f.write(f"f {t[0]} {t[1]} {t[2]}\n")
