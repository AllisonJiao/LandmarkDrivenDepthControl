"""Shared helpers for the B5 scripts: paths, loading a base scene, the heat-map
colour table, and unit conversions between the scene's normalized units and
the app's canvas pixels.

The scene loader is Part A's sds/scene.py (numpy/scipy only). It normalizes
x, y and h by Scene.R (the half-extent of the sketch), M by R^2 and therefore
u = sgn(h) h^2 by R^2; pressure is unchanged.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(ROOT / "sds"))

from scene import Scene  # noqa: E402  (sds/scene.py)

# Fields of an app export that describe one solve, not the base shape.
SOLVE_KEYS = ("landmarks", "method", "radius", "radius_px", "weight", "ok",
              "h", "du", "dp", "du_out", "dp_out", "dh_out")

# Same 9 stops as kHeatmapTable in src/viz/heatmap.cpp (matplotlib RdBu_r),
# so figures and the in-app heat map use identical colours.
HEATMAP_TABLE = np.array([
    [0.020, 0.188, 0.380], [0.129, 0.400, 0.675], [0.263, 0.576, 0.765],
    [0.573, 0.773, 0.871], [0.969, 0.969, 0.969], [0.957, 0.647, 0.510],
    [0.839, 0.376, 0.302], [0.698, 0.094, 0.169], [0.404, 0.000, 0.122],
])


def heatmap_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("mm_heatmap", HEATMAP_TABLE)
    cmap.set_bad((0.6, 0.6, 0.6))  # NaN (outline in dp), as in the app
    return cmap


def load_base(path):
    """A scene file used as the base shape. Its landmarks and solve fields are
    ignored by the B5 scripts (handles come from a handle file), so a scene
    exported after a failed solve is still a valid base."""
    return Scene(str(path), allow_failed=True)


def read_json(path):
    with open(path) as f:
        return json.load(f)


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=1)


# --- units: normalized (Scene) <-> app canvas pixels -----------------------

def length_to_px(scene, x):
    """Lengths and heights (x, y, h, dh) scale by R."""
    return np.asarray(x) * scene.R


def u_to_px(scene, u):
    """u = sgn(h) h^2 and du scale by R^2."""
    return np.asarray(u) * scene.R**2
