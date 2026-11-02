# Offline method comparison

Figures and metrics for the five handle methods. Inputs are a **base scene**
(mesh + default heights) and a **handle file** (where the handles are and how
much they change); the method is chosen at run time.

**The methods are solved by the app's own C++ solver** (`LandmarkSolver`),
built as a headless command-line tool (`src/tools/landmark_cli.cpp`) and run
with Node. Python resolves handles, calls the tool (`cpp_bridge.py`), computes
metrics and plots, so the figures show exactly what the app computes.

```bash
../build-web.sh                            # also builds build/web/cli/landmark_cli.js
pip install -r requirements.txt            # numpy, scipy, matplotlib
python make_grid_scenes.py                 # data/ellipse_nb{60,120,240,480}.json
python handles.py show --scene data/ellipse_nb120.json --handles configs/ellipse_x1.5.json
python cpp_bridge.py   --scene data/ellipse_nb120.json --handles configs/ellipse_x1.5.json
python plot_handles.py --scene data/ellipse_nb120.json --handles configs/ellipse_x1.5.json --out out/ellipse
python plot_refinement.py --scenes "data/ellipse_nb*.json" --handles configs/ellipse_x1.5.json --out out/refinement
```

Node: `node` on PATH (or set `MM_NODE`); the tool path can be overridden with
`MM_CLI`. Direct use of the tool:
`node build/web/cli/landmark_cli.js --scene in.json --out-dir DIR [--methods G1,G2,P1,P2,P3] [--radius 4] [--weight W]`
(the input is a scene with `landmarks: [{vertex, factor, anchor}]`; outputs
are scene-version-2 files, the same format as the app's **Export for SDS…**).

From the app: **Export for SDS…**, then

```bash
python handles.py split --scene mm_sds_scene.json --base data/heart_base.json --handles configs/heart.json
```

| File | Role |
|---|---|
| `common.py` | paths, base-scene loading (`sds/scene.py`), heat-map colours, unit conversion |
| `handles.py` | handle files: load/save, `split`, mesh fingerprint, position → vertex, validation, depth → factor |
| `cpp_bridge.py` | runs `landmark_cli` and returns its results as `Result` / `Fields` |
| `results.py` | `Result`, `Fields`, `Operators` / `operators()` (sheet labels, masses, L, mean edge) |
| `metrics.py` | one function per metric, on the C++ results; `all_metrics()` fills NaN for unwritten or failing ones |
| `plotting.py` | sheet triangulation, 2 × N grid, height profile, metrics tables |
| `plot_handles.py` | one scene + one handle set (+ sweep) → grid, profile, metrics |
| `plot_refinement.py` | the same handles on several resolutions → log-log metrics |
| `make_grid_scenes.py` | synthetic ellipse bases at several resolutions |
| `configs/` | handle files used for figures |

**Left to do (`TODO(B5)` in the code):**
- P3 is the app's current P3 (every row of L in the energy). Update the
  `p3` header in `plot_handles.py` and `plot_refinement.py` when Part A lands
  the interior-row version (no Python change is needed for the solve itself).
- A native (non-Node) build of the tool needs `struct Landmark` in an
  SDL-free header (see `src/CMakeLists.txt`); not needed with Node.

Depth units in handle files: `factor` (target |h| = f |h0|), `dh_px` (change
of |h| in canvas pixels), `dh_rel` (change of |h| as a fraction of max|h0|);
positive thickens on either sheet.

The synthetic ellipse bases are re-solved with libigl's Voronoi mass
(`make_grid_scenes.py --mass voronoi`, the default), so the C++ solver sees a
uniform default pressure +-c.

Handle positions (`"at"`) are in the base scene's normalized coordinates
(`Scene.V`, about [-1, 1]); `"side"` is the sign of h0 of the sheet. In the
synthetic ellipse the first sheet has h0 < 0 (`"side": -1`); in the app,
front parts have h0 > 0.

Windows note: the bridge passes an empty pipe as the tool's stdin; Emscripten's
NODERAWFS cannot stat the `NUL` device and node would exit before `main()`.
