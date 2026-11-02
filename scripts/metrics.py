"""Metrics of one solve, computed from the C++ solver's output.

Inputs (results.py), all in normalized units:
  scene  the C++ output scene (sds/scene.py): V, F, h0, u0, is_fixed, ...
  ops    Operators: interior, L, m (masses), sheet labels, mean_edge, c
  res    Result: method, du, vertices (handle ids), r (targets in u)
  f      Fields: u, h, du, dp, du_out, dp_out, dh_out, dh_rel, dp_rel
"""

import numpy as np

from common import Scene
from results import Fields, Operators, Result

# column order of the metrics table
METRIC_NAMES = (
    "handle_error", "area_moved", "dp_max", "dp_min", "dp_energy",
    "neg_dp_area", "dirichlet", "sheet_crossing", "other_sheet_du",
)

"""
Example of the pattern every metric follows:
    res.vertices  the handle vertex ids, in handle order (anchors included)
    res.r         the target change of u at each of them (same order)
    res.du / f.*  per-vertex arrays of length n, so index them with ids
"""

def handle_error(scene: Scene, ops: Operators, res: Result, f: Fields) -> float:
    """max_k |du(v_k) - r_k| / max|r|: how far the worst handle is from its
    target change of u (0 = every handle exactly on target).
    """
    du_at_handles = res.du[res.vertices]  # (k,) achieved change of u
    error = np.abs(du_at_handles - res.r)  # (k,) 1-D, so abs, not a norm
    scale = np.abs(res.r).max()
    # every factor 1 (or anchors only): r = 0, report the absolute error
    return float(error.max() / scale) if scale > 0 else float(error.max())


def area_moved(scene: Scene, ops: Operators, res: Result, f: Fields) -> float:
    """Area (sum of M) where |dh_out| > 1 % of the largest |dh_out| at a
    handle, as a fraction of the handle sheets' area; independent of the mesh
    resolution."""
    
    result = 0
    
    max_dh = np.max(np.abs(f.dh_out[res.vertices]))
    
    if max_dh > 0:
        # find the sheets that contain handles
        handle_sheets = np.unique(ops.sheet[res.vertices])
        
        # create a mask for vertices on those sheets
        mask = np.isin(ops.sheet, handle_sheets)
        
        # compute the area where |dh_out| > 1% of max_dh
        moved_area = np.sum(ops.m[mask & (np.abs(f.dh_out) > 0.01 * max_dh)])
        
        # compute the total area of the handle sheets
        total_area = np.sum(ops.m[mask])
        
        if total_area > 0:
            result = moved_area / total_area

    return result


def dp_max(scene: Scene, ops: Operators, res: Result, f: Fields) -> float:
    """max dp_rel over the interior."""
    # ops.interior is an array of vertex ids (not a mask): test its size
    if ops.interior.size == 0:
        return float("nan")
    return float(np.max(f.dp_rel[ops.interior]))


def dp_min(scene: Scene, ops: Operators, res: Result, f: Fields) -> float:
    """min dp_rel over the interior."""
    if ops.interior.size == 0:
        return float("nan")
    return float(np.min(f.dp_rel[ops.interior]))


def dp_energy(scene: Scene, ops: Operators, res: Result, f: Fields) -> float:
    """dp' M dp (the quantity P3 minimizes), normalized units."""
    I = ops.interior
    if I.size == 0:
        return float("nan")
    # dp with itself: dp_out = -sgn(u0) dp, so dp_out . dp is negative on the
    # u0 > 0 sheet. The sign drops out of dp'M dp anyway.
    return float(f.dp[I] @ (ops.m[I] * f.dp[I]))


def neg_dp_area(scene: Scene, ops: Operators, res: Result, f: Fields) -> float:
    """Fraction of the interior area with dp_out < 0 (less inflation)."""
    I = ops.interior
    if I.size == 0:
        return float("nan")
    m = ops.m[I]
    dp = f.dp_out[I]
    # ignore round-off: P1/P2 have dp ~ -1e-13 far from the handle, which
    # would otherwise count as negative pressure
    tol = 1e-9 * np.abs(dp).max()
    # an area (sum of masses), not a vertex count, so it doesn't depend on
    # the mesh resolution
    return float(m[dp < -tol].sum() / m.sum())


def dirichlet(scene: Scene, ops: Operators, res: Result, f: Fields) -> float:
    """du' (-L) du, the smoothness of the displacement (Dirichlet energy)."""
    if ops.interior.size == 0:
        return float("nan")
    # du with itself (du_out . du would flip sign on the u0 < 0 sheet)
    return float(f.du @ (-(ops.L @ f.du)))


def sheet_crossing(scene: Scene, ops: Operators, res: Result, f: Fields) -> float:
    """min over the interior of u / u0 (Part A's crossing measure): 1 for an
    untouched vertex, < 1 where the sheet moved toward the drawing plane,
    negative where a vertex crossed it."""
    I = ops.interior
    if I.size == 0:
        return float("nan")
    # per vertex, not / max|u0|: that global scale made every method read the
    # same small value set by the vertices next to the outline
    return float(np.min(f.u[I] / scene.u0[I]))  # u0 != 0 on the interior


def other_sheet_du(scene: Scene, ops: Operators, res: Result, f: Fields) -> float:
    """max |du| on sheets that carry no handle, relative to max|r|; 0 for
    one-sided methods."""
    handle_sheets = np.unique(ops.sheet[res.vertices])
    others = (ops.sheet >= 0) & ~np.isin(ops.sheet, handle_sheets)
    if not others.any():
        return 0.0  # every sheet carries a handle
    worst = np.abs(f.du[others]).max()
    scale = np.abs(res.r).max()
    return float(worst / scale) if scale > 0 else float(worst)


METRICS = {name: globals()[name] for name in METRIC_NAMES}


_reported = set()


def all_metrics(scene: Scene, ops: Operators, res: Result, f: Fields) -> dict:
    """{name: value}; NaN for metrics that are not implemented yet, and for
    metrics that raise (the error is printed once per metric, so a broken
    metric shows up without stopping the figures)."""
    out = {}
    for name, fn in METRICS.items():
        try:
            out[name] = float(fn(scene, ops, res, f))
        except NotImplementedError:
            out[name] = float("nan")
        except Exception as e:  # noqa: BLE001 - reported, not hidden
            out[name] = float("nan")
            if name not in _reported:
                _reported.add(name)
                print(f"metric {name} failed: {type(e).__name__}: {e}")
    return out
