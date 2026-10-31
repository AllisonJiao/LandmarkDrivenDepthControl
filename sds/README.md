# Score distillation in pressure space

Optimizes the inflation pressure of a Monster Mash sketch so that renders of the
inflated shape match a text prompt (DreamFusion-style score distillation with
Stable Diffusion). The shape stays an inflation of the drawing: x, y are fixed,
the outline stays clamped, and only the pressure changes, so the result can
still be edited with landmarks in the app.

## Setup (Windows, RTX 50 series)

```bash
python -m venv ../../tools/sdsenv
../../tools/sdsenv/Scripts/python -m pip install torch --index-url https://download.pytorch.org/whl/cu130
../../tools/sdsenv/Scripts/python -m pip install -r requirements.txt
```

Stable Diffusion weights (about 4 GB) download to the Hugging Face cache on the first run.

## Use

1. In the app: draw, Inflate, optionally place landmarks, then **Export for SDS…** (`mm_sds_scene.json`).
2. Run:
   ```bash
   ../../tools/sdsenv/Scripts/python run_sds.py --scene mm_sds_scene.json --prompt "a puffy cloud" --out runs/cloud
   ```
   `--rep spectral --k 20` restricts the pressure change to the first 20 Laplacian eigenmodes.
3. In the app: **Import SDS result…** with `runs/cloud/result.json`, then keep editing with landmarks.

Without the app, `make_test_scene.py` writes a synthetic inflated ellipse in the same format, and
`check.py --scene <file>` verifies the forward model, the silhouette and the renderer's gradients.

## Files

| File | Purpose |
|---|---|
| `scene.py` | loads the exported scene, normalizes coordinates, writes results |
| `pressure.py` | pressure change → heights (full or spectral), regularizer, landmark projection |
| `raster.py` | pure-PyTorch rasterizer with differentiable shading |
| `sds.py` | Stable Diffusion score distillation loss |
| `run_sds.py` | optimization loop, cameras, logging |
