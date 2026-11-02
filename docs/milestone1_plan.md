# Milestone 1 plan: Landmark-driven depth control (v2.2)

> **How to read the change markers**
> - **[Δn]**: changed in v2. The reason is row *n* of the table below.
> - **[new]**: added in v2, with no counterpart in v1.
> - **[clarified]**: same content as v1, made more precise.
> - ~~struck text~~: removed from v1.
> - **[v2.1]**: added on 2026-11-01, after the method wiring and vertex
>   landmarks (D1) were implemented. **[done]** marks plan items that now
>   exist in the code.
> - **[v2.2]**: added later on 2026-11-01, after B1–B3, the four baselines and
>   the amended M1 contract (`docs/M1_contract.md`) were implemented.
>
> Unmarked text is unchanged from v1.

Goal: a 5-minute presentation that shows (1) a usable handle-editing tool,
(2) a fair comparison between our pressure-space solve and simple baselines,
and (3) clear figures that explain *why* each result looks the way it does.

**Rule for every method:** each method produces **both** a displacement map
(`dh`, and `du` in the solver's field) **and** a pressure-change map (`dp`),
whichever of the two it edits directly.

---

## What changed from v1 **[new section]**

Every claim below was checked against the current code or measured on the
two-sheet test ellipse (`sds/test_ellipse_plain.json`, handle ×1.5, normalized
units).

| # | Change | Why |
|---|---|---|
| 1 | **P3 energy uses interior rows only.** Today it sums over every row, including the outline rows (`buildEnergy` in `src/landmarksolver.cpp`). | With every row, (a) the solver does not minimize the `dp'M dp` we report: 7.65 vs the true optimum of 5.73. (b) The front and back sheets are coupled: a front handle moves the back sheet by **37 %** of the handle's displacement. (c) The slide claim "dp ∝ L⁻¹Bᵀ" is false for it. |
| 2 | **Refinement claim corrected.** | P3's *peak* dp grows like log(1/h): in the v1 table 9.6 → 12.1 → 14.4 is about +2.4 per halving. Only its *energy* converges. G2 and P2 at a fixed R also converge. The new claim and the new metric (area moved) are in §3.2. |
| 3 | **G2 falloff is multiplicative:** `du_j = (k²−1)·u0_j·w(d_ij/R)` | The additive form `r·w` moves outline vertices when the disk crosses the outline, and it can push u below zero for k < 1. Multiplying by `u0` is zero on the outline and keeps `u > 0` for any k > 0. |
| 4 | **P3-R (region-restricted optimal pressure) moves from future work to Should.** | It has a closed form (§0) and costs about 10 lines on top of the shared helpers. With R = the handle vertex it equals P1 exactly; with R = everything it equals P3. On the same disk it always beats P2 (15.4 vs 21.8 at R = 0.15, 7.9 vs 10.4 at R = 0.3), which is the "why optimization" figure. |
| 5 | **Falloffs (G2, P2, P3-R) stay on the handle's own sheet.** | The front and back sheets share xy, so a 2D-distance disk also covers the back sheet, and overlapping parts as well. |
| 6 | **Several handles and anchors in G2/P1/P2/P3-R from day 1** (one k×k solve). | The UI already supports several landmarks and anchors. Single-handle versions would silently ignore anchors. |
| 7 | **Snapping skips outline vertices and duplicates.** | A handle on a fixed vertex is an infeasible constraint (0 = r). Two handles on one vertex make the KKT system singular. |
| 8 | **Keys apply to the web build only.** | `keyPressEvent` is desktop-only (`#ifndef __EMSCRIPTEN__`, `src/mainwindow.cpp:653`), and the desktop build is not set up. `+`/`−` already shift animation frames in `src/ui/main.js:213`, so the factor steps use `[`/`]`. `W` is free on the web. |
| 9 | **Front/back decided: edits are one-sided** (the landmark's own sheet). | This matches the existing per-side factors. With #1 every method is exactly one-sided; the back sheet's du is 0, which is a checkable invariant. |
| 10 | **The figures come from the existing JSON scene export plus Python**, not a new CSV schema. | `exportScene` already writes V, F, fixed, L, M and the landmarks, and downloads through the Save menu; `sds/scene.py` already loads it. A Python reference implementation of every method lets the figures start before the C++ lands, and it cross-checks the C++. |
| 11 | **New metric: sheet crossing,** `min sgn(u0)·u`. | No method clamps. A crossing vertex flips to the other side of the drawing plane. It did not occur on the ellipse even at ×0.3. Keep it for multi-part shapes and handles near the outline. |
| 12 | **Already done, so removed from Must:** the CMake fix (committed); L, M and the fixed set stored in `prepare()`; the scene export and download path; several landmarks and anchors in P3. | — |

---

## Status, 2026-11-01 **[v2.2]**

Part A (the solver) is done apart from the open items in the last row. The
agreements with Part B are in **`docs/M1_contract.md`**. §7 keeps the math
and the test results.

| State | Item | Where |
|---|---|---|
| **Done** | Method wiring: `HandleMethod` enum (G1 G2 P1 P2 P3), a dispatch into `setConstraintsP3/Baseline` and `solveP3/Baseline`. **[v2.2]** JS setters `setHandleMethod(int)` and `setFalloffRadius(double)`, with R in mean edge lengths | `landmarksolver.*`, `mainwindow.*`, `main.cpp` |
| **Done** | **D1 vertex landmarks** with the outline guard (D7) and the duplicate guard (A4); the Python side reads vertex ids. **[v2.2]** `setConstraints` also rejects a `vertexId` of −1, out of range, or on the outline | §7.1 |
| **Done [v2.2]** | B1–B3: interior map, K_II = −L_II with `SimplicialLDLT`, sheet components, `forward`, `implied` | §7.2 |
| **Done [v2.2]** | G1, G2, P1, P2 through one k×k basis solve (`factorBasis`, `falloff`) | §7.3 |
| **Done [v2.2]** | Decisions D2, D3, D4, D5, D7 | §7.4 |
| **Done [v2.2]** | The amended M1 contract: outward-positive fields, `lastSolveOk()`, scene export version 2 (`component`, `h`, `ok`, …), `Scene.vertex_table()` / `write_vertex_csv()` | `docs/M1_contract.md`, §3.1 |
| **Part B** | D6 (weight slider), the R slider range, the mesh name for the export, the `Landmark` struct merge, heat map, figures | `docs/M1_contract.md`, §6 |
| **Open (A)** | P3 energy over interior rows only (Δ1), still being decided; P3-R; refinement data; the C++/Python cross-check | §2, §6 |

All five methods run in the app. On the heart (11,549 vertices) a drag frame
takes 0.1–0.4 ms for the baselines.

---

## 0. Story for the talk

> "Monster Mash inflates with a constant pressure. We let the user drag a point
> on the surface and solve backwards for a non-constant pressure field. We
> compare this with editing the geometry directly and with changing the pressure
> locally."

Every method answers the same question: **the handle must reach height ×k.
Which field do we change, and where?**

### Notation
- `u = sgn(h)·h²` is the field in which inflation is linear. On the interior
  vertices I, `L_II u = M_II p`, with `u = 0` on the fixed set (every part's
  outline, `defEng.bnds`) **[clarified]**. L is the flat cotangent Laplacian
  and M is the Voronoi mass matrix.
- The handle at vertex i reaches `k·hᵢ` exactly when `duᵢ = r = (k²−1)·u0ᵢ`.
  **[new]** With non-snapped handles, `B u0` interpolates u while the screen
  interpolates h, so the displayed ratio is only approximately k. Snapping
  removes this.
- **Pressure → displacement:** `forward(dp) = L_II⁻¹ M_II dp`, one
  back-substitution with a prefactored `L_II` **[clarified: named helper]**.
- **Displacement → pressure:** `implied(du) = M_II⁻¹ L_I· du`, one sparse
  matrix-vector product, on interior vertices only **[clarified: named
  helper]**.
- **[new]** All comparisons use P3 in **hard** mode (weight slider at the
  right end).

### Methods
**[Δ1 Δ3 Δ4 Δ5 Δ6]**: the G2, P2, P3-R and P3 rows changed or are new, and
the "Optimization?" column became "Solve".

| | Method | Edits | The other field | Solve |
|---|---|---|---|---|
| **G1** | Move the handle vertex only: `du = r·eᵢ` | du (1 vertex) | `dp = implied(du)`: a 1-ring stencil with zero net force | none |
| **G2** | **[Δ3 Δ5]** Multiplicative falloff: `du_j = (k²−1)·u0_j·w(d_ij/R)`, `w(t) = (1−t²)³`, on the handle's sheet. ~~`du_j = r·w(d_ij/R)`~~ | du (disk) | `dp = implied(du)` | **[Δ6]** k×k when there are several handles |
| **P1** | Pressure spike: `dp = α·eᵢ` | dp (1 vertex) | `du = α·forward(eᵢ)`, `α = r / forward(eᵢ)ᵢ` | **[Δ6]** k×k |
| **P2** | Pressure bump: `dp = α·φ`, `φ_j = w(d_ij/R)` **[Δ5]** on the handle's sheet | dp (disk) | `du = α·forward(φ)` | **[Δ6]** k×k |
| **P3-R** | **[Δ4 new row]** Optimal pressure inside the disk: `dp = α·mask_R(L_II⁻¹eᵢ)` | dp (disk) | `du = α·forward(dp)` | k×k plus one back-substitution per handle |
| **P3** | Minimize `dp'M dp` subject to `B du = r` **[Δ1]** (interior rows) | dp (everywhere) | existing solver plus `implied(du)` | KKT |

### Unified view (one slide)
Write the pressure change as `dp = Φ·α`, where the basis Φ has m columns.
- **P1:** Φ = eᵢ. **P2:** Φ = a bump. Both have m = k, so a k×k solve gives a
  unique answer.
- **P3:** Φ = the identity, so m = n > k, and a selection rule is needed:
  minimize `dp'M dp`.
  - **[clarified: derivation added]** Stationarity gives
    `M dp = Gᵀ Bᵀ λ = M L⁻¹ Bᵀ λ`, so **`dp = L⁻¹Bᵀλ`**. P3's pressure has
    the shape of the handle's harmonic Green's function, which is the shape of
    P1's *displacement*.
  - **[new]** P3's displacement `L⁻¹ M L⁻¹ Bᵀ λ` is the biharmonic Green's
    function, which is smooth in 2D (`r² log r`), whereas P1's is a log spike.
    P3 is "P1 made one Laplacian smoother".
- **[Δ4]** **P3-R:** restrict dp to a region R. The same derivation (M is
  diagonal and L is symmetric) gives `dp_R = (L⁻¹Bᵀ)_R λ`: the global
  optimum's pressure, **masked to R and rescaled**. R = {i} gives P1 and
  R = everything gives P3, so **one radius slider morphs P1 into P3**.
  ~~"Future work: Φ = the vertices inside a user region R … a small version of
  P3."~~
- **[new]** **P2 vs P3-R at the same R** separates "what locality buys" from
  "what optimization buys": same support, fixed bump vs optimal shape.
- Optimization is needed exactly when the pressure has more degrees of freedom
  than there are constraints.

### Known numbers
**Unit-square grid, handle at ×1.5, base inflation = 1 (v1, unchanged):**

| Grid N | P1 peak dp | G1 peak dp | P3 peak dp | P3 `dp'M dp` |
|---|---|---|---|---|
| 21 | 50 | 122 | 9.6 | 1.25 |
| 41 | 171 | 490 | 12.1 | 1.39 |
| 81 | 593 | 1961 | 14.4 | 1.46 |

**[Δ2]** Peak dp grows as 1/h² for G1, as 1/(h² log(1/h)) for P1, and as
**log(1/h)** for P3, where h is the edge length. P3's energy converges: the
increments 0.14 and 0.07 halve with each refinement, so the limit is about
1.53. Under refinement P1's displacement narrows into G1's spike
(logarithmically slowly), while P3's bump keeps its width. ~~"P3's pressure
energy levels off"~~ was correct but incomplete.

**[new] Two-sheet ellipse, ×1.5 handle off-centre (measured for v2):**

| | `dp'M dp` | peak dp | back-sheet du / r |
|---|---|---|---|
| P1 = P3-R (R = 0) | 129.5 | 274.8 | 0 |
| P2 / P3-R, R = 0.15 | 21.8 / 15.4 | 46.5 / 32.8 | 0 |
| P2 / P3-R, R = 0.30 | 10.4 / 7.9 | 16.1 / 16.8 | 0 |
| P3, interior rows | **5.73** | 12.2 | 0 |
| P3, every row (current C++) | 7.65 | 15.6 | **0.374** |

---

## 1. Handle control upgrades

1. **Wireframe overlay.** **[Δ8]** The `W` key (free in `main.js`) plus a
   checkbox.
   - **[new]** Compute `E = igl::edges(F)` once per reconstruction, not every
     frame.
   - **[new]** `GLMeshDraw(…, GL_LINES)` already draws a two-column index
     matrix (`third_party/miscutils/opengltools.cpp`). It needs its **own**
     `GLMeshData`, because refilling the mesh VAO would overwrite the triangle
     indices, and a flat-colour shader, because both existing shaders are
     matcap or textured.
   - Set `glPolygonOffset` on the fill pass. Where:
     `MainWindow::drawModelOpenGL` (`src/mainwindow.cpp:455`).
2. **Snap handles to vertices** in `pickLandmark` (`src/landmarks.cpp:58`):
   **[v2.1 done]**, except the highlight (see §7.1).
   - **[Δ7]** Choose the face corner with the largest barycentric weight that
     is **not fixed**, and reject the pick if all three corners are fixed.
     This needs the fixed set passed in (`isFixed` from `LandmarkSolver`).
   - ~~Set the weights one-hot.~~ **[v2.1]** The landmark stores a vertex id
     (`Landmark::vertexId`) instead of a face and weights. **[new]** Set
     `sketchPos` to the vertex's xy so `rebindLandmarks` re-snaps after a
     redraw.
   - **[Δ7]** If the vertex already has a handle, select that handle instead
     of adding a duplicate. **[v2.1]** `rebindLandmarks` also drops a second
     landmark that snaps to a taken vertex.
   - Highlight the snapped vertex.
3. **Numeric readout:** ×factor (exists), absolute `h`, vertex id,
   **[new]** method. **[Δ8]** **`[` / `]`** step the factor by ±0.1 in
   Landmarks mode, so demos are reproducible. ~~± keys~~
4. ~~Desktop key for Landmarks mode (e.g. `K`)~~ **[Δ8]**: dropped.
5. *(Could)* Save landmarks in the project file (`loadsave.cpp`), or
   hard-code the demo scenes.

---

## 2. Methods implementation

### 2.1 Structure
- **[Δ4]** Add `enum class HandleMethod { G1, G2, P1, P2, P3R, P3 }` with
  `setHandleMethod(int)` and **[clarified]** ~~`setHandleRadius(double)`~~
  **[v2.2]** `setFalloffRadius(double)`
  exported next to `setLandmarkWeight` (`src/main.cpp:137`), plus a dropdown
  and an R slider in `myshell.html` / `main.js`.
  **[v2.1 done]**, without P3R for now. ~~The method and the radius go
  through a single `setHandleMethod(int method, double radius)`.~~
  **[v2.2]** Two setters, `setHandleMethod(int)` and
  `setFalloffRadius(double)`, with R in mean edge lengths (M1 contract). The
  R slider shows only for G2 and P2.
- **[Δ12]** `prepare()` already stores L, M and `fixedIdx`. Add
  (**[v2.2 done]**; the components come from `igl::vertex_components` on
  K_II and are stored as `sheetLabels`):
  - the interior list and its inverse map, `L_II`, and a `SimplicialLDLT` of
    `−L_II` (SPD);
  - **[Δ5 new]** **sheet components**: connected components of the mesh with
    the fixed vertices removed, giving one id per sheet of each part.
    Falloffs, export and plots use them.
- **[Δ1 new]** **Energy:** `Q = Lᵀ S M⁻¹ S L`, where S zeroes the outline rows
  (a one-line change in `buildEnergy`). Keep the every-row version only as an
  optional "P3c" figure about boundary conditions.
- Helpers: `forward(dp) → du` and `implied(du) → dp`, both with zeros on the
  fixed set. **[v2.2 done]**
- **[Δ6 new]** **`basisSolve(Φ)`**, shared by P1, P2 and P3-R
  (**[v2.2 done]** as `factorBasis()` plus one per-frame case; G2 shares it
  too, with the column u₀∘w):
  - When the landmark set changes, compute `g = forward(Φ)` (one
    back-substitution per handle).
  - Per drag frame: solve `(B g) α = r` (k×k), then `du = g α`, `dp = Φ α`.
  - This is the same factor-once, back-substitute-per-drag pattern as P3, so
    every method updates live while dragging.
- Every method writes `VRestOut` (the rest of the pipeline is unchanged) and
  stores `lastDu` / `lastDp`. **[v2.2]** Also the outward-positive fields and
  `lastSolveOk()` (M1 contract).

### 2.2 Per-method notes
- **G1:** **[Δ6]** `du = Σ r_k e_{i_k}` over the snapped handles. No solve.
- **G2:** **[Δ3 Δ5 Δ6]** `du_j = u0_j · Σ_k c_k w(d_jk/R)`, with c from a k×k
  system (overlapping falloffs). For one handle, `c = k²−1`. Use 2D distance,
  restricted to the handle's component. Geodesic distance is optional.
  ~~"Apply the falloff in u, not in h"~~: still true; the multiplicative form
  is in u.
- **P1:** `Φ = eᵢ` per handle. For a non-snapped handle, use `Φ = Σ bⱼ eⱼ`
  over the face's interior corners. **[v2.2]** Not needed in M1: landmarks
  are always snapped (M1 contract).
- **P2:** `Φ_k = w(d/R)` **[Δ5]** on the handle's component. As R → 0 it
  becomes P1.
- **P3-R:** **[Δ4 new]** `Φ_k = mask_R(L_II⁻¹ B_kᵀ)`, where R is the union of
  the handles' disks on their own components. Then the same `basisSolve`.
- **P3:** the existing solver **[Δ1]** with the interior-row Q, plus
  `dp = implied(du)`.
- **[Δ11 new]** **No clamping in any method.** Sheet crossing is reported, not
  hidden.
- ~~"Several handles in P1/P2 (optional for M1)"~~ **[Δ6]**: now required.
- ~~"Two-sided mesh: decide whether an edit changes thickness or only the
  front"~~ **[Δ9]**: decided, one-sided.

### 2.3 Sanity checks **[clarified: moved here from Student A's list and extended]**
- With no handles, or every factor = 1: du = 0 for every method.
- Handle error ≈ 0 for every method.
- `forward(implied(du)) = du` for any du with du = 0 on the fixed set.
- G1's net force `Σ M dp = 0` **[clarified]** when the handle is not next to
  the outline.
- R → 0: P2 = **[Δ4]** P3-R = P1. **[Δ4]** R → ∞: P3-R = P3.
- **[Δ9 new]** Front handle: the back sheet's du = 0 exactly, for every
  method.
- **[Δ10 new]** C++ P3 = the Python projection
  (`PressureModel.project_landmarks()` from θ = 0, `sds/pressure.py`) on the
  same exported scene. That projection is an independent implementation of
  interior-row P3.

**[v2.2] Results** on the two-sheet test ellipse (2,424 vertices, mean edge
6.7 px):
- Factor 1: du = 0 exactly, every method.
- Handle error: 0 to about 1e-14 for G1, G2, P1 and P2.
- `implied(forward(e_v)) = e_v` to 1.3e-15.
- Net force ≈ 1e-12 for G1 and ≈ 1e-11 for G2 when the disk stays off the
  outline. When it crosses the outline, the net force equals −Σ (L du) over
  the outline rows exactly.
- R below one edge: G2 = G1 and P2 = P1 bit for bit.
- Front handle: the back sheet's du and dp are exactly 0 for G1, G2, P1 and
  P2. P3 waits for Δ1.
- **[new]** The outward fields satisfy K·du_out = M·dp_out to about 1e-16
  relative, for all five methods.
- **[new]** Sheet crossing does occur: G2 and P2 with a ×1.0 anchor on the
  handle's 1-ring reach min u/u₀ = −1.4 (G2) and −1.8 (P2) at R = 7.4 edges.
  P1 never crosses.
- Not run yet: C++ P3 vs the Python projection.

---

## 3. Result demonstration

### 3.1 Export **[Δ10]**
~~"Write one CSV per method run. Per vertex: x, y, h0, h, dh, u0, u, du, dp,
is_boundary, sheet."~~

Extend `LandmarkSolver::exportScene` (it already writes V, F, fixed, L, M
and the landmarks, downloaded through the Save menu) with `method`, `R`,
`h`, `du`, `dp` and `component`. Write one JSON file per run; `sds/scene.py`
loads it.
- Use dp, not absolute p. On the stitched mesh, `M⁻¹L u0` is not constant
  near part seams.
- Mask the fixed vertices in dp.
- **[Δ10 new]** `scripts/handle_methods.py` recomputes every method from a
  scene in Python (dense G, as in `sds/pressure.py`). The figures do not wait
  for the C++, and the two implementations check each other.

**[v2.2 done]** `exportScene` writes version 2: `method`, `radius` (mean
edge lengths), `radius_px`, `weight` (`null` = hard), `ok`, `h`, `du`, `dp`,
the outward fields `du_out`, `dp_out`, `dh_out`, plus `component`,
`mean_edge`, `c` and `mesh_name`. In `sds/scene.py`,
`Scene.vertex_table()` and `write_vertex_csv()` give the per-vertex table
with the contract's header lines, and a scene with `ok: false` is refused.
The CSV the M1 contract asks for is generated from this JSON.

### 3.2 Offline figures, `scripts/plot_handles.py`
- **Main grid, 2 × 6** **[Δ4]** (G1 G2 P1 P2 P3-R P3; 2 × 5 without P3-R).
  Row 1 is `dh` and row 2 is `dp`.
  - **[Δ5 new]** Plot **only the handle's component**. Front and back overlap
    in xy, so a `tripcolor` of the whole mesh draws both sheets on top of
    each other.
  - Use `RdBu_r` with symmetric limits shared within each row, a star at the
    handle, and a symmetric-log variant for dp.
- **Height profile:** a 1D cut through the handle **[new]**
  (`LinearTriInterpolator` on the handle's component), showing h0 and h for
  every method.
- **Metrics table:**
  - handle error;
  - **[Δ2]** **area moved** `Σ M_j·[|dh_j| > 1 % |dh_handle|]`, which does
    not depend on the resolution. ~~number of vertices moved~~
  - max/min dp, `dp'M dp`, and the area fraction with negative dp;
  - `du'(−L)du`;
  - **[Δ11 new]** sheet crossing `min sgn(u0)·u`;
  - **[Δ9 new]** back-sheet `max |du|`.
- **Refinement plot:** peak |dp|, `dp'M dp` and **[Δ2]** area moved against h.
  - **[new]** Data: the ellipse from `make_test_scene.py` at 4 resolutions
    (no app needed), plus one real shape with `triangleOpts = "pqa25QYY"`
    (`src/commonStructs.h:129`, where `a` is the maximum triangle area) set to
    a12 and a6.
  - Show G1, P1 and P3, plus **[Δ2]** G2/P2 at a fixed R as dashed lines.
  - **[Δ2 new]** Expected: G1's and P1's area moved shrinks toward a spike;
    G2 and P2 converge but need R; P3 converges with no R.
- **R sweep** (*Should*): **[Δ4]** `dp'M dp` against R for P2 and P3-R.
  P3-R lies below P2 at every R.

### 3.3 In-app heat map (live demo)
**[new]** `GLMeshFillBuffers` already uploads a `color` attribute when `C` is
non-empty, but neither shader reads it. Add `attribute vec3 color` and a
varying to `third_party/shaderMatcap/shaderMatcap.h`, and mix it into the
matcap colour. Fill `C` in `drawModelOpenGL` (`src/mainwindow.cpp:481`).
A toggle cycles **Off → dh → dp**. Use the same colormap as the plots, and
draw the range as on-screen text.

---

## 4. Five-minute presentation outline

| Time | Slide | Content |
|---|---|---|
| 0:00–0:30 | Problem | Monster Mash inflates with a constant pressure, which makes depth hard to control locally. |
| 0:30–1:15 | Idea | `L u = M p`. Drag a handle and invert for the pressure. The two maps: forward and implied. |
| 1:15–2:15 | **Live demo / video** | Wireframe on, snapped handle, drag it, switch methods, **[Δ4]** R slider (P1 → P3-R → P3), dh/dp heat map. |
| 2:15–3:30 | Baselines | **[Δ4]** The 2 × 6 grid and the profile. "Local in height means a sharp, zero-net-force pressure. Local in pressure means a kinked height change. P3 is smooth in both." |
| 3:30–4:15 | Why optimization | The unified Φ view. **[Δ4]** P2 vs P3-R at the same R. **[Δ2]** Refinement: point methods collapse to a spike; radius methods need R; P3 needs neither. ~~"only P3 stays bounded"~~. Caveats: **[Δ2]** P3's peak dp grows like log(1/h), negative pressure, **[new]** P3's global falloff (P3-R answers it), z-only edits. |
| 4:15–5:00 | Next steps | A user-drawn region for P3-R **[Δ4]**, saving landmarks, outline (xy) edits, **[new]** the pressure-space diffusion prior (`sds/`, prototype exists). |

Record a **backup demo video**. Use 2–3 shapes: a simple blob, a multi-part
character from `data/examples/` **[clarified: path]**, and optionally a long
thin part.

---

## 5. Priorities

**[Δ12 new]** **Done:** CMake fix; L, M and the fixed set in `prepare()`; the
scene JSON export and download path; P3 with several landmarks, anchors and a
soft weight; the Python pressure model and landmark projection.
**[v2.1]** Also done: the method dropdown and R slider with their setter;
vertex landmarks with the outline and duplicate guards (C++ and Python);
`du`, `dp`, `method` and `radius` in the export.
**[v2.2]** Also done: B1–B3 with the sheet components; G1, G2, P1 and P2
through one basis solve; the export fields; the M1 contract's setters,
outward fields and solve flag.

**Must:**
- **[Δ1]** interior-row P3 energy, checked against Python;
- snapping **[Δ7]** with the outline and duplicate guards **[v2.1 done]**;
  wireframe;
- G1, G2, P1, P2 and P3 **[Δ6]** through `basisSolve`, with du/dp
  **[v2.2 done]**;
- **[Δ10]** the export fields (~~CSV export~~) **[v2.2 done]**;
- **[Δ10]** the Python reference, 2 × 5 grid, profile and metrics;
- slides and demo video.

~~"commit the CMake fix"~~ **[Δ12]**: done.

**Should:** **[Δ4]** P3-R and the R sweep; in-app heat map; refinement plot;
**[Δ8]** `[`/`]` keys (~~± keys~~).

**Could:** save/load landmarks; ghost overlay of the default shape;
triharmonic P3; geodesic falloff; **[Δ1 new]** the P3c boundary-condition
figure; **[Δ9 new]** a mirror-thickness option applied to every method alike.
~~desktop key~~ **[Δ8]**. ~~several handles in P1/P2~~ **[Δ6]**: moved to
Must.

---

## 6. Two-student task division

The split follows the boundary between the solver and the figures. On day 1,
agree on **[Δ10]** **the export fields (§3.1)** (~~CSV schema~~), **[Δ5 new]**
**component ids**, and **the `HandleMethod` enum and its setter signatures**.

### Student A: solver and methods (C++)
1. **[Δ1 new]** Interior-row `buildEnergy`. Check it against
   `project_landmarks` on an exported scene.
2. Extend `prepare()` with the interior map, the `L_II` LDLT and **[Δ5]** the
   components. Add `forward` / `implied` / **[Δ6]** `basisSolve`.
   ~~"Commit the CMake fix"~~ **[Δ12]**. **[v2.2 done]**
3. Implement G1, G2, P1, P2 **[Δ4]** (and P3-R), with `lastDu` / `lastDp` for
   every method, P3 included. **[v2.2 done]**, except P3-R.
4. Add **[Δ10]** the export fields and the JS setters for the method and R.
   **[v2.2 done]**, as amended in `docs/M1_contract.md`.
5. Run the sanity checks **[clarified]** (§2.3). Produce the refinement data
   **[new]** (ellipse resolutions plus `triangleOpts`). **[v2.2]** Checks
   done except the C++/Python comparison; refinement data open.

### Student B: interaction, visualization, presentation
1. Add the wireframe, snapping **[Δ7]** with its guards, the readout and
   **[Δ8]** the `[`/`]` keys. ~~desktop mode key~~
2. Add the method dropdown and R slider, wired to A's setters (use stubs until
   they land).
3. Add the heat map (**[new]** shader colour attribute, C buffer, range text).
4. **[Δ10 new]** Write `scripts/handle_methods.py` (the Python reference,
   built on `sds/scene.py`) and `scripts/plot_handles.py`.
5. Choose the demo shapes, record the backup video and draft the slides.

### Shared
- **Midpoint sync:** A's first exported scene goes through B's scripts.
  **[Δ10 new]** Compare the C++ and Python du/dp; they must agree to about
  1e-6.
- **Final 1–2 days:** both review the figures, rehearse and write the speaker
  notes. A presents the method and unified-view slides; B presents the demo
  and results.
- **Ownership** (to avoid conflicts):
  - A owns `landmarksolver.*` and the export code.
  - B owns `drawModelOpenGL`, `landmarks.*` (picking and snapping),
    `src/ui/*` and **[new]** `third_party/shaderMatcap/`.
  - **[Δ7 new]** Snapping needs `isFixed` from A's solver; agree on that
    getter on day 1. **[v2.1 done]**: `fixedMask()`.
  - **[v2.2]** D1 put the snapping in `landmarks.*` on A's branch, and the M1
    contract keeps A's snapping logic in the merge (outline check, duplicate
    rejection, `sketchPos` kept on rebind). B owns the rest of `landmarks.*`
    and `src/ui/*`.
  - Coordinate the `main.cpp` exports directly.

---

## 7. Implementation guide: what is done and the math for what is next **[v2.1]**

**[v2.2]** Everything in §7.1–7.4 is now in the code, and its
`TODO(human)` tags are gone, so this section names functions instead of line
numbers. It keeps the math and adds what was measured.

### 7.0 Notation for this section
- n vertices. The **fixed set** is `fixedIdx`, or `isFixed` as one flag per
  vertex: every part's outline. The **interior** I is everything else.
- L is libigl's cotangent Laplacian. It is negative semi-definite, so write
  **K = −L**, which is positive. mᵢ = M(i,i) is the vertex area.
- u₀ = sgn(h₀)·h₀² is the signed default field. It is **negative on the back
  sheet**.
- Handle k sits at vertex vₖ with factor fₖ. Its target is
  **rₖ = (fₖ² − 1)·u₀(vₖ)**, and `solve()` already computes r.
- Every method returns a du with **B du = r** and du = 0 on the fixed set.
  `solve()` then writes h = sgn(u)·√|u| with u = u₀ + du, and stores
  dp = `implied(du)`.

### 7.1 Done: D1 vertex landmarks
Every method's constraint reads the field at the handles. With a landmark
on a vertex, each handle is exactly one unknown, du(vₖ), and row k of B is
e_{vₖ}ᵀ. G1 and P1 are defined by a single vertex (eᵥ in du or in dp), so a
point inside a face would smear them over three vertices.

| Change | Where | Role in the algorithm |
|---|---|---|
| `vertexId` replaces `faceId` + `bary` | `src/landmarks.h` | One unknown du(vₖ) per handle. G1's column is eᵥ, and P1's column is forward(eᵥ), the response to pressure at one vertex. |
| Rows of B are e_{vₖ}ᵀ | `setConstraints` in `src/landmarksolver.cpp` | B only picks out values: (B u₀)ₖ = u₀(vₖ). The baselines' k×k matrix becomes a lookup, A(m,k) = cₖ(vₘ), and G1's A is the identity. |
| `pickCorner`: the largest weight among corners not on the outline (D7) | `src/landmarks.cpp:47` | Snaps to the corner nearest the click on screen. An outline vertex has u₀ = 0, so r = 0 and forward(eᵥ) = 0: its column of A would be zero, and A singular. |
| `fixedMask()`, built in `prepare()` | `prepare` in `src/landmarksolver.cpp` | The fixed set as one flag per vertex. It is also the first piece of B1: the interior map for K_II. |
| `surfacePoint = V.row(v)` | `src/landmarks.cpp` | The marker is drawn exactly where the constraint applies. |
| Rebind re-snaps; keeps `sketchPos`; sets `baseHeight = \|h₀(v)\|` | `rebindLandmarks`, `src/landmarks.cpp:137` | Vertex ids change with every reconstruction. Keeping the click position avoids drift, and measuring the height at the vertex keeps rₖ consistent with the drawn target. |
| No two landmarks on one vertex | click: `src/mainwindow.cpp:1437`; rebind: `src/landmarks.cpp` | Two identical rows of B make P3's hard-mode KKT system and the baselines' A singular. |
| Export `"vertex"`; Python reads it | `exportScene`; `sds/scene.py` (`landmark_matrix`, `_landmark_vertex`); `sds/make_test_scene.py` | The SDS projection Δp ← Δp − M⁻¹GᵀBᵀ(BGM⁻¹GᵀBᵀ)⁻¹(BGΔp − r) uses the same B as the C++ code. Old face/bary scene files still load: they snap with the same rule. |

**Verified on 2026-11-01:**
- Old (face/bary) and new (vertex) scene files load. Each row of B has one
  nonzero, (B u₀)ₖ = u₀(vₖ), no landmark is on the outline, and each
  landmark's side matches sgn(h₀) at its vertex.
- `sds/check.py` on a two-landmark scene meets the targets exactly
  (×1.500 and ×0.700, constraint error 3e-8).
- `./build-web.sh` succeeds. In the app, clicks snap to vertices and P3
  solves in about 1 ms on 11,549 vertices.

**Not exercised:**
- The click-side duplicate branch never fired. A vertex's snap area reaches
  about edge/√3 ≈ 9 px, just inside the 10 px hit test, so it only fires in
  rare cases such as faces along the outline.
- The rebind-side duplicate guard is compiled only.
- The snapped-vertex highlight from §1 is not done.

### 7.2 Done: per-mesh operators (B1–B3) **[v2.2]**
**B1. `prepare()`.** After the `isFixed` lines and before `ready = true`. The
members are `interiorIdx`, `toInterior`, `massI` and `ldlt`, plus
`sheetLabels` (D2) and `meanEdge` (D3).

**Math:** K_II = −L restricted to interior rows and columns. It is symmetric
positive definite: it is the finite-element stiffness matrix with zero values
on the outline, and every sheet touches the outline.

```text
toInterior[0..n) = -1;  interior = []
for i in 0..n-1: if !isFixed[i]: toInterior[i] = interior.size; interior.append(i)
for each nonzero (i, j, value) of L:                # InnerIterator loop, as in exportScene
    if toInterior[i] >= 0 and toInterior[j] >= 0:
        triplets.add(toInterior[i], toInterior[j], -value)    # K_II = -L_II
K_II.setFromTriplets(...);  mI[a] = M.diagonal()(interior[a])
ldlt.compute(K_II);  if ldlt.info() != Eigen::Success: return   # ready stays false
```

- **Members:** `interiorIdx`, `toInterior`, `mI`, and
  `Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt`.
- **Includes:** `SimplicialLDLT` arrives through the `-mod` headers that are
  already included. Don't add `<Eigen/SparseCholesky>`, which brings back the
  include-order errors.

**B2. `forward()`**: from pressure to displacement.

  L_II du_I = M_II dp_I, so **du_I = L_II⁻¹ M_II dp_I = −K_II⁻¹ (m ∘ dp)_I**,
  with du = 0 on the fixed set.

```text
b[a] = mI[a] * dp(interior[a])
x = ldlt.solve(b)
du = zeros(n);  du(interior[a]) = -x[a]          # minus sign: we factored -L
```

One back-substitution: about 1 ms at 11.5k vertices.

**B3. `implied()`**: from displacement to pressure.

  **dp_I = M_II⁻¹ (L du)_I**, and dp = 0 on the fixed set. Since du = 0 on the
  fixed set, (L du)_I = L_II du_I, so `implied` is exactly the inverse of
  `forward`.

```text
Ldu = L * du                                      # sparse mat-vec over all rows
dp = zeros(n);  for interior i: dp(i) = Ldu(i) / M(i,i)
```

- **Self-check:** `implied(forward(e_v))` must equal `e_v` to about 1e-12.
  This catches sign and indexing mistakes before any shape is drawn.
- **Sign:** L is negative semi-definite, so a bump in u has a Laplacian of the
  opposite sign. A dp that pushes a vertex away from the drawing plane has sign
  −sgn(u₀). For display or export, −sgn(u₀)·dp makes "pushes outward" positive
  on both sheets. **[v2.2]** Decided in D4: `lastDpOutward()` returns it.
- **Verified:** `implied(forward(e_v)) = e_v` to 1.3e-15, `forward` is 0 on
  the fixed set, and the LDLT of K_II succeeds on the heart (11,549
  vertices).

### 7.3 Done: the four baselines as one basis solve **[v2.2]**
G1, G2, P1 and P2 differ only in a **column** cₖ ∈ ℝⁿ per handle: the du that
coefficient 1 on handle k produces.

  du = Σₖ αₖ cₖ = C α. Every target needs du(vₘ) = rₘ, so **A α = r** with
  **A(m,k) = cₖ(vₘ)**, a k×k matrix.

This is §2.1's `basisSolve`: `(B g) α = r` there is `A α = r` here.

| Method | Column cₖ | Pressure basis | A |
|---|---|---|---|
| **G1** | e_{vₖ} | — | identity |
| **G2** | u₀ ∘ wₖ | — | wₖ(vₘ)·u₀(vₘ) |
| **P1** | forward(e_{vₖ}) | e_{vₖ} | −m_{vₖ}(K_II⁻¹)_{vₘvₖ} |
| **P2** | forward(wₖ) | wₖ | read off the columns |

- **Setup, once per landmark set:** in `setConstraintsBaseline`:
  1. Build C (n×k), the member `basis`.
  2. Read off A(m,k) = C(vₘ, k), i.e. A = B·C.
  3. Factor A with `Eigen::FullPivLU<Eigen::MatrixXd>` (`basisLU`). Return
     false if `!isInvertible()`. Steps 2 and 3 are `factorBasis()`.
- **Per drag frame:** one case in `solveBaseline`, shared by G2, P1 and P2:
  `alpha = basisLU.solve(r); duOut = basis * alpha`. That takes
  microseconds.
- **Why LU:** A is not symmetric. For P1, A = −S·D, with S a block of K_II⁻¹
  and D the handle masses.
- **Placeholder:** **[v2.2]** removed; every case returns its own result.

**G1** **[v2.2 done]**
- **Math:** du = Σ rₖ e_{vₖ}. A = I, so α = r.
- **Setup:** store `handles`.
- **Solve:** `duOut(handles[k]) = r(k)`.
- **Expect:** only k vertices move. dp is a 1-ring stencil, and Σ mᵢ dpᵢ ≈ 0
  when the handle isn't next to the outline (L's rows sum to zero).
- **Measured:** du is nonzero only at the handles, dp(v) = r·L_vv/m_v, and
  the net force is about 1e-12.

**G2** **[v2.2 done]**
- **Falloff:** wₖ(j) = (1 − (dₖⱼ/R)²)³ when dₖⱼ < R, vertex j is on handle
  k's sheet (D2) and j is not fixed. Otherwise wₖ(j) = 0. dₖⱼ is the xy
  distance on `VDefault`, in canvas pixels; R is set in mean edge lengths and
  converted with `meanEdge` (D3). This is `falloff(v)`, which P2 shares.
- **Column:** cₖ = u₀ ∘ wₖ.

```text
for k: for j in 0..n-1:
    if isFixed[j] or !sameSheet(j, k): skip
    t = dist(j, v_k) / R
    if t < 1: C(j, k) = u0(j) * (1 - t*t)^3
A(m, k) = C(v_m, k);  factor
```

- **One-handle check:** α = r / u₀(v) = f² − 1, so
  hⱼ = h₀ⱼ·√(1 + (f² − 1)wⱼ). This scales the thickness from f at the handle
  down to 1 at the rim (Δ3).
- **Why this falloff:** (1 − t²)³ has w, w′ and w″ all zero at t = 1, so the
  implied dp has no ridge at the rim.
- **Measured:** α = f² − 1 and h = h₀√(1 + αw) to 1e-14. At R = 3.7, 7.4 and
  14.9 edges, 28, 97 and 344 of the sheet's 1,152 vertices move (P2 moves all
  of them). dp has both signs: a push inside and a pull ring.

**P1** **[v2.2 done]**
- **Math:** dp = Σ αₖ e_{vₖ}, pressure at the handle vertices only. The column
  cₖ = forward(e_{vₖ}) is the discrete **Green's function**: the shape's
  response to pressure at one vertex.
- **Why A is invertible:** for distinct interior handles, A is a principal
  block of K_II⁻¹ (positive definite) times a positive diagonal.
- **Setup:** one `forward` per handle (k back-substitutions), then A, then the
  factorization.
- **Solve:** the same as G1's.
- **Expect:** du moves the whole sheet, and `implied(du)` is nonzero **only**
  at the handle vertices. That makes a strong correctness check.
- **Measured:** `implied(du)` = α at the handles and about 1e-14 elsewhere.
  du is exactly 0 on the other sheet, and A's entries between the two sheets
  are exactly 0. A ×1.0 anchor next to a handle gets a negative α (suction).

**P2** **[v2.2 done]**
- **Math:** the pressure basis φₖ = wₖ: G2's falloff without the u₀ factor,
  on the interior and the handle's sheet. The column is cₖ = forward(φₖ).
- **Setup:** build φₖ, run `forward`, then A, then the factorization.
- **Solve:** the same as G1's.
- **Expect:** `implied(du)` ≈ Σ αₖ φₖ, a smooth bump. When R is below one edge
  length, P2 becomes P1.
- **Measured:** `implied(du)` = Σ αₖ φₖ to 1e-14, and R below one edge gives P1
  bit for bit. At R = 3.7, 7.4 and 14.9 edges the peak dp is 24, 12 and 6.3
  (P1: 100), while the net force grows from 3.9k (P1) to 28k: a wider bump
  lowers the peak but adds more air.

**Later, P3-R** (Δ4): the same machinery with cₖ = forward(φₖ), where
φₖ = mask_R(K_II⁻¹ e_{vₖ}).

### 7.4 Decisions **[v2.2]**
- **D2: which sheet a falloff may touch.** **[v2.2 done]** The robust option,
  computed with `igl::vertex_components` on K_II instead of a flood fill:
  K_II's off-diagonal entries are exactly the edges between interior
  vertices. Stored as `sheetLabels` (−1 on the outline) and exported as
  `component`. The options as first written:
  - **One-part shapes:** `sameSheet(j, k) = sgn(u0(j)) == lms[k].side`.
  - **Several parts (robust):** connected components with the fixed vertices
    removed, computed once in `prepare`:
    ```text
    comp[0..n) = -1;  id = 0
    for s in interior with comp[s] == -1:
        flood fill from s over mesh edges, never entering fixed vertices; comp = id
        id++
    sameSheet(j, k) = comp[j] == comp[v_k]
    ```
    `igl::adjacency_list(F, adj)` gives the neighbours. Every part outline is
    fixed, so each part's front and back become separate components. This
    also keeps a disk on a leg from spilling onto an overlapping body at the
    same xy. The components are also the `component` export field (§3.1).
- **D3: the units of R.** **[v2.2 done]** R is in mean edge lengths, as the
  M1 contract specifies: `setFalloffRadius(edges)`, converted with the mean
  x,y edge length `meanEdge` from `prepare()`. Default 4 edges. The slider
  range belongs to Part B (`[D3 radius]` in `src/ui/myshell.html`).
  - **Units:** `VDefault` x and y are canvas pixels.
  - **Edge length:** `triangleOpts = "pqa25QYY"` allows at most 25 px² on the
    half-resolution image, so 100 px² on the canvas, and edges are about
    10–15 px. R ≈ 10 px makes G2 ≈ G1 and P2 ≈ P1.
  - ~~**Option:** scale R by the size of the sheet, for example
    ρ·√(sheet area), so one slider value means the same locality on any
    drawing.~~ Edge lengths give that, and also stay comparable under mesh
    refinement.
  - **Keep the HTML default equal to `LandmarkSolver::radiusEdges` (4).**
- **D4: signs on the back sheet.** **[v2.2 done]** The solver keeps signed
  u. The reported fields are outward positive (`lastDuOutward()`,
  `lastDpOutward()`, `lastDhOutward()`; M1 contract decision 1). The test
  below passes for all five methods.
  - **What already works:** u₀ < 0 on the back sheet, so r < 0 for f > 1. That
    makes u more negative, and h = sgn(u)·√|u| moves *away* from the plane.
    G2's signed u₀ and the α of P1/P2 carry the sign automatically.
  - **What breaks it:** using |u₀| or h₀ where the signed u₀ belongs. Also, if
    u crosses 0, the square root flips that vertex to the other sheet (Δ11).
    At the handle, u = f²u₀ keeps its sign; elsewhere that isn't guaranteed
    (P methods with f < 1).
  - **Test:** put a back-sheet handle at f = 1.5. In every method it must move
    away from the plane.
- **D5: anchors** **[v2.2 done]**: nothing extra. An anchor is a handle whose
  r comes from its held factor; its own column lets the k×k solve keep it
  still while the other handles move.
- **D6: the weight slider** **[Part B]** (`src/ui/main.js`, `[D6 weight]`):
  the baselines ignore `weight`. In the UI, for example
  `$('#landmarkWeight').prop('disabled', method != 4)`, or hide it.
- **D7: outline** **[done]**: `pickLandmark` and `rebindLandmarks` never put a
  handle on a fixed vertex. **[v2.2]** `setConstraints` also rejects such
  landmarks, and a `vertexId` of −1, when they come from somewhere else.

### 7.5 Test checklist
These extend §2.3. **[v2.2]** The status of each item is at the end of it;
the numbers are in §2.3.
1. All factors 1 → r = 0 → du = 0 for every method. **Passes.**
2. At each handle, the new height is exactly fₖ·h₀(vₖ): the target marker sits
   on the surface. **Passes** (to about 1e-14).
3. `implied(forward(e_v)) = e_v`. **Passes** (1.3e-15).
4. P1: `implied(du)` is nonzero only at the handles. P2: `implied(du)` ≈ α·φ.
   **Passes.**
5. G1: Σ mᵢ dpᵢ ≈ 0 when the handle is away from the outline. **Passes.**
6. R below one edge: G2 = G1 and P2 = P1. **Passes** (bit for bit).
7. A back-sheet handle moves away from the plane (D4). A front handle leaves
   the back sheet untouched (D2, Δ9). **Passes** for the baselines; P3's
   back sheet waits for Δ1.
8. **Export for SDS…**, then compare `du` / `dp` with the Python reference on
   the same scene (§3.1). **Open:** needs `scripts/handle_methods.py`.

### 7.6 Order of work
1. **B1**, then **B2/B3** with test 3. **[v2.2 done]**
2. **G1**, which needs no `forward`. **[v2.2 done]**
3. **P1**, the first real test of `forward`. **[v2.2 done]**
4. **P2** with **D2**. **[v2.2 done]**
5. **G2**. **[v2.2 done]**
6. **D3** and **D6** in the UI. **[v2.2]** D3's units are done; the slider
   range and D6 are Part B's.
7. Then the open items in the status table: interior-row P3 (Δ1), P3-R, the
   export fields **[v2.2 done]** and the figures.
