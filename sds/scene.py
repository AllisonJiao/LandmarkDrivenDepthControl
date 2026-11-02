"""Scene exported from Monster Mash ("Export for SDS"): the default rest mesh,
the exact flat operators used by the C++ landmark solver, and the landmarks.

Coordinates are normalized so that the sketch spans about [-1, 1]. Under the
scaling x -> x / R the 2D cotangent Laplacian is unchanged, the mass matrix
scales by 1 / R^2 and so does u = sgn(h) h^2, which leaves the pressure
p = M^-1 L u unchanged.

Version 2 scenes (the M1 contract's export) also carry the run they were
exported after: method, radius, ok, the solved heights and the fields in the
solver's and the outward-positive signs. vertex_table() and
write_vertex_csv() give them per vertex, in the app's units.
"""

import csv
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
    def __init__(self, path, allow_failed=False):
        with open(path) as f:
            d = json.load(f)

        # a scene exported after a failed solve has zero fields that look
        # like a valid "nothing changed" run
        self.ok = bool(d.get("ok", True))
        if not self.ok and not allow_failed:
            raise ValueError(f"{path}: exported after a failed landmark solve "
                             "(pass allow_failed=True to load it anyway)")

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
        self.is_fixed = is_fixed
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

        # the exported run, in the app's units (absent before version 2)
        self.mesh_name = d.get("mesh_name", "")
        self.method = d.get("method")
        self.radius = d.get("radius")  # mean edge lengths in version 2
        self.radius_px = d.get("radius_px")
        self.weight = d.get("weight")  # None: hard constraints
        self.c = float(d.get("c", 2.0))  # default pressure
        self.mean_edge = d.get("mean_edge")
        self.component = (np.asarray(d["component"], dtype=np.int64)
                          if "component" in d else None)
        self.fields = {k: np.asarray(d[k], dtype=np.float64)
                       for k in ("h", "du", "dp", "du_out", "dp_out", "dh_out")
                       if k in d}

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

    def header(self):
        """Description of the exported run, as the M1 contract's CSV header."""
        return {
            "mesh_name": self.mesh_name,
            "method": self.method,
            "k": len(self.landmarks),
            "R_edges": self.radius,
            "R_px": self.radius_px,
            "weight": "hard" if self.weight is None else self.weight,
            "c": self.c,
            "max_abs_h0": float(np.abs(self.V_input[:, 2]).max()),
            # vertex:factor per landmark, ":anchor" on anchors
            "handles": ";".join(f"{l.vertex}:{l.factor:g}" + (":anchor" if l.anchor else "")
                                for l in self.landmarks),
        }

    def vertex_table(self):
        """Per-vertex fields of the exported run in the app's units, with the
        M1 contract's normalizations dh_rel = dh_out / max|h0| and
        dp_rel = dp_out / c. Columns are numpy arrays, so
        pandas.DataFrame(scene.vertex_table()) works."""
        if "dh_out" not in self.fields:
            raise ValueError("scene has no solve fields (exported before version 2)")
        f = self.fields
        h0 = self.V_input[:, 2]
        return {
            "vertex": np.arange(self.n),
            "x": self.V_input[:, 0],
            "y": self.V_input[:, 1],
            "h0": h0,
            "h": f["h"],
            "dh_out": f["dh_out"],
            "du_out": f["du_out"],
            "dp_out": f["dp_out"],
            "dh_rel": f["dh_out"] / np.abs(h0).max(),
            "dp_rel": f["dp_out"] / self.c,
            "component": self.component,
            "is_fixed": self.is_fixed,
        }

    def write_vertex_csv(self, path):
        """vertex_table() as CSV after "# key: value" header lines, which
        pandas.read_csv(path, comment="#") skips."""
        table = self.vertex_table()
        with open(path, "w", newline="") as f:
            for key, value in self.header().items():
                f.write(f"# {key}: {value}\n")
            writer = csv.writer(f)
            writer.writerow(table.keys())
            writer.writerows(zip(*table.values()))

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
