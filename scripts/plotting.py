"""Figure helpers for the B5 scripts: one sheet of a scene as a triangulation,
the method-comparison grid, the height profile and the metrics tables.

Figures use the app's canvas coordinates (Scene.V_input, pixels, y down), so
they look like the sketch on screen, and the app's heat-map colours.
"""

import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import Normalize, SymLogNorm  # noqa: E402

from common import heatmap_cmap  # noqa: E402

FIELD_LABELS = {
    "dh_rel": r"$dh_{out}\,/\,\max|h_0|$",
    "dp_rel": r"$dp_{out}\,/\,c$",
}


# --- one sheet as a triangulation --------------------------------------------------

def sheet_faces(scene, sheet, s):
    """Faces of sheet s: every corner is on s or on the outline, and at least
    one corner is on s. Front and back sheets overlap in xy, so figures show
    one sheet at a time."""
    labels = sheet[scene.F]
    on_s = labels == s
    ok = on_s | (labels == -1)
    return scene.F[ok.all(1) & on_s.any(1)]


def sheet_triangulation(scene, sheet, s):
    xy = scene.V_input[:, :2]
    return mtri.Triangulation(xy[:, 0], xy[:, 1], sheet_faces(scene, sheet, s))


def face_values(tri, values):
    """Per-face mean over corners with a finite value (NaN marks the outline
    in dp), so outline-only faces stay grey."""
    v = values[tri.triangles]
    finite = np.isfinite(v)
    with np.errstate(invalid="ignore"):
        out = np.where(finite, v, 0).sum(1) / finite.sum(1)
    return np.where(finite.any(1), out, np.nan)


def symmetric_limit(arrays):
    vals = [np.abs(a[np.isfinite(a)]) for a in arrays if a is not None]
    vals = [v for v in vals if v.size]
    vmax = max((v.max() for v in vals), default=0.0)
    return vmax if vmax > 0 else 1.0


def make_norm(vmax, symlog=False, linthresh_frac=1e-2):
    if symlog:
        return SymLogNorm(linthresh=vmax * linthresh_frac, vmin=-vmax, vmax=vmax, base=10)
    return Normalize(vmin=-vmax, vmax=vmax)


def draw_handles(ax, scene, vertices, anchors):
    xy = scene.V_input[vertices, :2]
    for (x, y), is_anchor in zip(xy, anchors):
        if is_anchor:
            ax.plot(x, y, "s", ms=6, mfc="none", mec="k", mew=1.2)
        else:
            # hollow, so a one-vertex field (G1, P1 dp) stays visible under it
            ax.plot(x, y, "*", ms=12, mfc="none", mec="k", mew=0.9)


def style_axes(ax):
    ax.set_aspect("equal")
    ax.invert_yaxis()  # canvas y points down
    ax.set_xticks([])
    ax.set_yticks([])


# --- the method grid ---------------------------------------------------------------------

def plot_method_grid(scene, sheet, s, columns, path, rows=("dh_rel", "dp_rel"),
                     symlog_rows=("dp_rel",), symlog=False, wireframe=True, title=None):
    """columns: list of (method, Fields or None, vertices, anchors). None draws
    a 'not implemented' panel. Limits are symmetric and shared within a row."""
    tri = sheet_triangulation(scene, sheet, s)
    cmap = heatmap_cmap()
    ncol = len(columns)
    fig, axes = plt.subplots(len(rows), ncol, figsize=(2.6 * ncol + 1.2, 2.4 * len(rows) + 0.6),
                             squeeze=False, constrained_layout=True)
    for i, row in enumerate(rows):
        values = [None if f is None else getattr(f, row) for _, f, _, _ in columns]
        vmax = symmetric_limit(values)
        norm = make_norm(vmax, symlog and row in symlog_rows)
        mappable = None
        for j, (method, f, vertices, anchors) in enumerate(columns):
            ax = axes[i, j]
            style_axes(ax)
            if i == 0:
                ax.set_title(method)
            if f is None:
                ax.text(0.5, 0.5, "not implemented", ha="center", va="center",
                        transform=ax.transAxes, color="0.4")
                continue
            mappable = ax.tripcolor(tri, facecolors=face_values(tri, getattr(f, row)),
                                    cmap=cmap, norm=norm, edgecolors="none")
            if wireframe:
                ax.triplot(tri, lw=0.15, color="k", alpha=0.25)
            draw_handles(ax, scene, vertices, anchors)
        axes[i, 0].set_ylabel(FIELD_LABELS.get(row, row))
        if mappable is not None:
            fig.colorbar(mappable, ax=axes[i, :], shrink=0.85, pad=0.01)
    if title:
        fig.suptitle(title, fontsize=10)
    fig.savefig(path, dpi=200)
    plt.close(fig)


# --- the height profile --------------------------------------------------------------------

def plot_profile(scene, sheet, s, h0, curves, through, path, title=None, samples=400):
    """|h| along the horizontal line through vertex `through` on sheet s.
    curves: list of (label, h per vertex in pixels). h0 in pixels."""
    tri = sheet_triangulation(scene, sheet, s)
    used = np.unique(tri.triangles)
    x0, y0 = scene.V_input[through, :2]
    xs = np.linspace(scene.V_input[used, 0].min(), scene.V_input[used, 0].max(), samples)
    ys = np.full_like(xs, y0)

    fig, ax = plt.subplots(figsize=(6.4, 3.2), constrained_layout=True)
    base = mtri.LinearTriInterpolator(tri, np.abs(h0))(xs, ys)
    ax.plot(xs, base, "k--", lw=1.2, label="default")
    for label, h in curves:
        ax.plot(xs, mtri.LinearTriInterpolator(tri, np.abs(h))(xs, ys), lw=1.4, label=label)
    ax.axvline(x0, color="0.6", lw=0.8)
    ax.set_xlabel("x (canvas px)")
    ax.set_ylabel("|h| (px)")
    ax.legend(fontsize=8, ncol=3)
    if title:
        ax.set_title(title, fontsize=10)
    fig.savefig(path, dpi=200)
    plt.close(fig)


# --- tables ----------------------------------------------------------------------------------

def write_metrics(rows, header, csv_path, md_path, names):
    """rows: list of (label, {metric: value}). header: {key: value} written as
    '# key: value' lines (pandas.read_csv(path, comment='#') skips them)."""
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        for k, v in header.items():
            f.write(f"# {k}: {v}\n")
        w = csv.writer(f)
        w.writerow(["run", *names])
        for label, m in rows:
            w.writerow([label, *(m.get(n, float("nan")) for n in names)])
    with open(md_path, "w", encoding="utf-8") as f:
        for k, v in header.items():
            f.write(f"- **{k}:** {v}\n")
        f.write("\n| run | " + " | ".join(names) + " |\n")
        f.write("|---" * (len(names) + 1) + "|\n")
        for label, m in rows:
            cells = ["—" if not np.isfinite(m.get(n, np.nan)) else f"{m[n]:.4g}" for n in names]
            f.write(f"| {label} | " + " | ".join(cells) + " |\n")
