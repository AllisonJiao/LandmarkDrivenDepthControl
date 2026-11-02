"""Handle files (part_b_plan.md B5.1a): where the handles are, how much each
changes, and which are anchors. The method is not part of a handle file; it is
run configuration, so one handle set can be solved under every method.

Format (JSON):
{
  "version": 1,
  "mesh": {"n": 2424, "faces_sha1": "..."},     # required for vertex handles
  "depth_unit": "factor",                        # factor | dh_px | dh_rel
  "handles": [
    {"vertex": 812, "depth": 1.5},
    {"at": [0.3, 0.27], "side": -1, "depth": 1.5},
    {"vertex": 905, "anchor": true}
  ],
  "sweep": {"handle": 0, "depths": [0.7, 1.0, 1.5]}   # optional
}

- "vertex": a vertex id of the base mesh; the "mesh" fingerprint must match.
- "at" + "side": a position in the base scene's NORMALIZED coordinates
  (Scene.V, sketch spans about [-1, 1]) and the sign of h0 of the sheet,
  the same convention as Landmark.side in the app. Resolves to the nearest
  interior vertex on that side. Use this when the mesh changes (refinement).
  Note: in the app front parts have h0 > 0 (side +1); in the synthetic ellipse
  of sds/make_test_scene.py the first sheet has h0 < 0 (side -1).
- "depth": in depth_unit; anchors default to factor 1 (hold the default
  height) when they give no depth.

CLI:
  python handles.py split  --scene mm_sds_scene.json --base base.json --handles h.json
  python handles.py show   --scene base.json --handles h.json
"""

import argparse
import hashlib
from dataclasses import asdict, dataclass, field
from typing import Optional

import numpy as np

from common import SOLVE_KEYS, load_base, read_json, write_json

DEPTH_UNITS = ("factor", "dh_px", "dh_rel")


@dataclass
class Handle:
    vertex: Optional[int] = None
    at: Optional[tuple] = None  # (x, y), normalized scene coordinates
    side: Optional[int] = None  # sign of h0 of the sheet, for "at" handles
    depth: Optional[float] = None
    anchor: bool = False


@dataclass
class HandleSet:
    handles: list = field(default_factory=list)
    depth_unit: str = "factor"
    mesh: Optional[dict] = None
    sweep: Optional[dict] = None


# --- file I/O ---------------------------------------------------------------

def load_handles(path):
    d = read_json(path)
    unit = d.get("depth_unit", "factor")
    if unit not in DEPTH_UNITS:
        raise ValueError(f"{path}: depth_unit must be one of {DEPTH_UNITS}, got {unit!r}")
    handles = []
    for i, h in enumerate(d.get("handles", [])):
        has_vertex, has_at = "vertex" in h, "at" in h
        if has_vertex == has_at:
            raise ValueError(f"{path}: handle {i} needs exactly one of 'vertex' or 'at'")
        if has_at and h.get("side") not in (1, -1):
            raise ValueError(f"{path}: handle {i} with 'at' needs 'side': 1 or -1")
        handles.append(Handle(
            vertex=int(h["vertex"]) if has_vertex else None,
            at=tuple(float(x) for x in h["at"]) if has_at else None,
            side=int(h["side"]) if "side" in h else None,
            depth=float(h["depth"]) if h.get("depth") is not None else None,
            anchor=bool(h.get("anchor", False)),
        ))
    return HandleSet(handles=handles, depth_unit=unit, mesh=d.get("mesh"),
                     sweep=d.get("sweep"))


def save_handles(hs, path):
    out = {"version": 1, "depth_unit": hs.depth_unit, "handles": []}
    if hs.mesh is not None:
        out["mesh"] = hs.mesh
    for h in hs.handles:
        out["handles"].append({k: (list(v) if isinstance(v, tuple) else v)
                               for k, v in asdict(h).items()
                               if v is not None and not (k == "anchor" and not v)})
    if hs.sweep is not None:
        out["sweep"] = hs.sweep
    write_json(path, out)


def mesh_fingerprint(scene):
    """Identifies the base mesh, so vertex ids are only used on their mesh."""
    faces = np.ascontiguousarray(scene.F, dtype=np.int64)
    return {"n": int(scene.n), "faces_sha1": hashlib.sha1(faces.tobytes()).hexdigest()}


def from_scene_landmarks(scene):
    """The landmarks stored in an app export, as vertex handles (factors)."""
    return HandleSet(
        handles=[Handle(vertex=l.vertex, depth=l.factor, anchor=l.anchor)
                 for l in scene.landmarks],
        depth_unit="factor", mesh=mesh_fingerprint(scene))


def split_scene(app_export, base_out, handles_out):
    """One app export -> a base scene file (without landmarks and solve fields)
    and a handle file (its landmarks as vertex handles)."""
    d = read_json(app_export)
    base = {k: v for k, v in d.items() if k not in SOLVE_KEYS}
    base["landmarks"] = []
    write_json(base_out, base)
    save_handles(from_scene_landmarks(load_base(app_export)), handles_out)


# --- resolving handles on a base scene ---------------------------------------

def nearest_vertex(scene, at, side):
    """Interior vertex closest in xy to `at` among those with sgn(h0) = side."""
    candidates = np.nonzero(~scene.is_fixed & (np.sign(scene.h0) == side))[0]
    if len(candidates) == 0:
        raise ValueError(f"no interior vertex with side {side}")
    d2 = ((scene.V[candidates, :2] - np.asarray(at)) ** 2).sum(1)
    return int(candidates[np.argmin(d2)])


def resolve_vertices(scene, hs):
    """Vertex id of every handle, validated like the C++ (setConstraints and
    pickLandmark): in range, not on the outline, h0 != 0, no duplicates."""
    if any(h.vertex is not None for h in hs.handles):
        fp = mesh_fingerprint(scene)
        if hs.mesh is None or hs.mesh != fp:
            raise ValueError("vertex handles were made for a different mesh "
                             f"(file: {hs.mesh}, base: {fp}); use 'at' handles instead")
    vertices = []
    for i, h in enumerate(hs.handles):
        v = h.vertex if h.vertex is not None else nearest_vertex(scene, h.at, h.side)
        if not 0 <= v < scene.n:
            raise ValueError(f"handle {i}: vertex {v} out of range")
        if scene.is_fixed[v] or scene.h0[v] == 0:
            raise ValueError(f"handle {i}: vertex {v} is on the outline")
        if v in vertices:
            raise ValueError(f"handle {i}: vertex {v} already has a handle")
        vertices.append(v)
    return np.array(vertices, dtype=np.int64)


def handle_depths(hs, depth_override=None):
    """Depth of each handle in hs.depth_unit. Anchors without a depth hold the
    default height. depth_override = (handle index, depth) replaces one depth,
    for sweeps."""
    hold = {"factor": 1.0, "dh_px": 0.0, "dh_rel": 0.0}[hs.depth_unit]
    depths = []
    for i, h in enumerate(hs.handles):
        d = h.depth
        if d is None:
            if not h.anchor:
                raise ValueError(f"handle {i} has no depth")
            d = hold
        depths.append(d)
    depths = np.array(depths, dtype=np.float64)
    if depth_override is not None:
        depths[depth_override[0]] = depth_override[1]
    return depths


def handle_factors(scene, vertices, depths, depth_unit):
    """Depths as height factors f (target |h| = f |h0|), the input of the C++
    solver (landmark_cli).
      factor  f as given
      dh_px   a change Dh of |h| in canvas pixels
      dh_rel  a change of |h| as a fraction d of max|h0| over the mesh
    For the two dh units, with Dh_n the change in normalized units
    (dh_px: Dh / scene.R; dh_rel: d * max|h0|), f = (|h0[v]| + Dh_n) / |h0[v]|.
    A positive change thickens the shape on either sheet (outward positive).
    f <= 0 is refused: the handle would cross the drawing plane."""
    depths = np.asarray(depths, dtype=np.float64)
    if depth_unit == "factor":
        factors = depths
    elif depth_unit in ("dh_px", "dh_rel"):
        base = np.abs(scene.h0[vertices])  # normalized; > 0 (validated)
        if depth_unit == "dh_px":
            dh = depths / scene.R
        else:
            dh = depths * np.abs(scene.h0).max()
        factors = (base + dh) / base
    else:
        raise ValueError(f"unknown depth unit {depth_unit!r}")
    if np.any(factors <= 0):
        raise ValueError("a target height <= 0 would push a handle through the drawing plane")
    return factors


# --- CLI ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("split", help="app export -> base scene + handle file")
    s.add_argument("--scene", required=True)
    s.add_argument("--base", required=True)
    s.add_argument("--handles", required=True)
    s = sub.add_parser("show", help="print the resolved handles on a base scene")
    s.add_argument("--scene", required=True)
    s.add_argument("--handles", required=True)
    args = ap.parse_args()

    if args.cmd == "split":
        split_scene(args.scene, args.base, args.handles)
        print(f"wrote {args.base} and {args.handles}")
    else:
        scene, hs = load_base(args.scene), load_handles(args.handles)
        vertices = resolve_vertices(scene, hs)
        depths = handle_depths(hs)
        print(f"{len(vertices)} handles on {args.scene} (n = {scene.n}), unit {hs.depth_unit}")
        for h, v, d in zip(hs.handles, vertices, depths):
            x, y = scene.V[v, :2]
            print(f"  vertex {v:6d} at ({x:+.3f}, {y:+.3f}) side {int(np.sign(scene.h0[v])):+d}"
                  f"  depth {d:g}{'  anchor' if h.anchor else ''}")


if __name__ == "__main__":
    main()
