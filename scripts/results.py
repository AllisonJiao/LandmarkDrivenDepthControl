"""Data types shared by the bridge, the metrics and the plots: one solve's
Result and Fields, and per-mesh Operators.

The solves themselves come from the app's C++ solver (cpp_bridge.py), so
Operators is read from the C++ output scene: its L and Voronoi M are the ones
the solve used.
"""

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.csgraph

METHODS = ("G1", "G2", "P1", "P2", "P3")


@dataclass
class Operators:
    """Per-mesh quantities, normalized units (sds/scene.py)."""
    n: int
    interior: np.ndarray  # interior vertex ids (not fixed)
    L: sp.csr_matrix  # flat cotangent Laplacian (negative semi-definite)
    m: np.ndarray  # vertex masses (diagonal of M), all vertices
    sheet: np.ndarray  # sheet label per vertex, -1 on fixed
    mean_edge: float  # mean xy edge length (the unit of R)
    c: float  # default pressure the dp normalization uses


@dataclass
class Result:
    """One method's solve. du and r in normalized units."""
    method: str
    du: np.ndarray  # change of u = sgn(h) h^2, per vertex
    vertices: np.ndarray  # handle vertex ids (anchors included), in handle order
    r: np.ndarray  # target change of u at each handle: (f^2 - 1) u0[v]
    settings: dict = field(default_factory=dict)


@dataclass
class Fields:
    """Per-vertex fields of one solve. Normalized units unless noted."""
    u: np.ndarray
    h: np.ndarray
    du: np.ndarray
    dp: np.ndarray  # solver signs, 0 on the outline
    du_out: np.ndarray
    dp_out: np.ndarray
    dh_out: np.ndarray
    dh_rel: np.ndarray  # dh_out / max|h0| (unitless, M1 contract)
    dp_rel: np.ndarray  # dp_out / c (unitless, M1 contract), NaN on the outline


def sheet_labels(scene):
    """The exported "component" field (C++ sheet labels), or for scenes
    without it the connected components of the interior over mesh edges
    (every outline is fixed, so each part's front and back are separate)."""
    if scene.component is not None:
        return np.asarray(scene.component, dtype=np.int64)
    F = scene.F
    i = np.concatenate([F[:, 0], F[:, 1], F[:, 2]])
    j = np.concatenate([F[:, 1], F[:, 2], F[:, 0]])
    keep = ~scene.is_fixed[i] & ~scene.is_fixed[j]
    A = sp.coo_matrix((np.ones(keep.sum()), (i[keep], j[keep])), shape=(scene.n, scene.n))
    _, labels = scipy.sparse.csgraph.connected_components(A, directed=False)
    sheet = np.full(scene.n, -1, dtype=np.int64)
    _, sheet[scene.interior] = np.unique(labels[scene.interior], return_inverse=True)
    return sheet


def mean_edge_length(scene):
    """Mean xy edge length in normalized units: the exported "mean_edge"
    (canvas pixels) / R, or computed from the mesh."""
    if scene.mean_edge:
        return float(scene.mean_edge) / scene.R
    E = np.concatenate([scene.F[:, [0, 1]], scene.F[:, [1, 2]], scene.F[:, [2, 0]]])
    E = np.unique(np.sort(E, axis=1), axis=0)
    xy = scene.V[:, :2]
    return float(np.linalg.norm(xy[E[:, 0]] - xy[E[:, 1]], axis=1).mean())


def operators(scene):
    return Operators(n=scene.n, interior=scene.interior, L=scene.L.tocsr(), m=scene.M,
                     sheet=sheet_labels(scene), mean_edge=mean_edge_length(scene),
                     c=scene.c)
