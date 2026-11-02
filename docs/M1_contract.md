# M1 contract between Part A (solver) and Part B (UI, measurements)

Amended 2026-11-01. Changed rows are marked **[changed]**; crossed-out text is
the original proposal.

## B0. Agreed items

| Item | Agreement |
|---|---|
| Method enum | `enum class HandleMethod { G1, G2, P1, P2, P3 }` (values 0–4) in `landmarksolver.h` |
| Setters exported to JS | `setHandleMethod(int)`, `setFalloffRadius(double)` in edge lengths. **[changed]** An edge length is the mean x,y length over all edges of the default rest mesh (`LandmarkSolver::meanEdge`, recomputed after every reconstruction). Default R = 4 edges. |
| Fields B can read | `lastDu()`, `lastDp()` in the solver's signs, and **[changed]** `lastDuOutward()`, `lastDpOutward()`, `lastDhOutward()` in the contract's signs. All are n×1, 0 on the outline, and 0 until the first solve after a reconstruction. **[changed]** `lastSolveOk()` is false before the first solve and after a failed one (rejected landmark, singular setup); the fields are 0 then and must not be recorded. |
| Sign convention | "Outward positive": thickening is positive on both sheets. **[changed]** Part A computes it: du_out = sgn(u0)·du, dp_out = −sgn(u0)·dp, dh_out = sgn(h0)·(h − h0). A sheet crossing shows as a large negative dh_out. This is the convention of `docs/math_explanation.md` (K du = M dp). |
| Normalization | dp reported as dp_out / c (c = `defaultInflationAmount` = 2, exported as `c`); dh as dh_out / max\|h0\| over the mesh. dp / c is relative to the default pressure exactly only for one-part shapes with the default inflation amount. |
| ~~CSV schema~~ Export **[changed]** | ~~milestone plan §3.1, plus a header line~~ One JSON file per run from `exportScene` (fields below). `sds/scene.py` turns it into the per-vertex table and CSV: `Scene(path).vertex_table()`, `.write_vertex_csv(path)`, with the header (mesh name, method, k, R in edges and pixels, weight, c, max\|h0\|, handles) as `# key: value` lines. |
| Snapped handle | `Landmark::vertexId`, which every method reads. **[changed]** ~~−1 when not snapped~~ Always set in M1 (snapping is always on). −1 is invalid: `setConstraints` rejects it, and also out-of-range and outline vertices. After the merge the struct holds `vertexId`, `faceId` and `bary`; `faceId`/`bary` are for drawing and a later barycentric version. |

~~Until A's fields exist, use P3 only.~~ All five methods provide the fields.

## Export fields (JSON version 2)

| Field | Meaning | Units |
|---|---|---|
| `mesh_name` | name passed from JS, `""` if none | |
| `V`, `F`, `fixed` | default rest mesh, faces, outline vertices | canvas pixels |
| `L` (`i`, `j`, `v`), `M` | flat cotangent Laplacian, Voronoi mass diagonal | |
| `component` | sheet of each vertex (each part's front and back are separate), −1 on the outline | |
| `mean_edge`, `c` | edge length unit, default pressure | pixels, pressure |
| `landmarks[]` | `vertex`, `factor`, `side`, `anchor`, `baseHeight` | |
| `method`, `radius`, `radius_px` | the last solve's method and R | edges, pixels |
| `weight` | P3's soft weight, `null` = hard (the baselines ignore it) | |
| `ok` | whether the last solve succeeded | |
| `h` | solved rest heights | pixels |
| `du`, `dp` | solver signs | |
| `du_out`, `dp_out`, `dh_out` | outward positive | |

Mesh name from JS: `Module.ccall('exportSdsScene', null, ['string'], [name])`.

## Conventions to expect

- **No clamping (plan Δ11).** G2 and P2 with two handles much closer than R
  can push part of a sheet through the drawing plane. Report it as
  min sgn(u0)·u, or as negative dh_out.
- **The baselines hit their targets exactly.** dh_out at a handle is
  (f − 1)·|h0|.

## Left to Part B

- R slider range in edge lengths (`src/ui/myshell.html`, `[D3 radius]`).
- Hide or disable the weight slider for G1–P2 (`src/ui/main.js`, `[D6 weight]`).
- After a failed solve: keep showing the previous shape or fall back to the
  default inflation.
- Pass a mesh name to `exportSdsScene`.
- Merge: the `Landmark` struct above, Part A's snapping (outline check,
  duplicate rejection, `sketchPos` kept on rebind), and remove or force on
  `snapToVertex`.
