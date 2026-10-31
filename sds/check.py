"""Sanity checks for the pressure model and the renderer.

  python check.py --scene test_ellipse.json
"""

import argparse

import numpy as np
import torch
from PIL import Image

from pressure import PressureModel
from raster import orbit_camera, rasterize, render
from scene import Scene


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scene", required=True)
    ap.add_argument("--rep", choices=["full", "spectral"], default="full")
    args = ap.parse_args()
    torch.manual_seed(0)

    scene = Scene(args.scene)
    model = PressureModel(scene, rep=args.rep)
    xy = torch.tensor(scene.V[:, :2], dtype=torch.float32, device="cuda")
    Fc = torch.tensor(scene.F, dtype=torch.long, device="cuda")

    # 1. no pressure change reproduces the exported default exactly
    h = model.heights()
    err = (h - model.h0).abs().max().item()
    print(f"[1] default heights reproduced: max |dh| = {err:.2e} (h range {h.min():.3f}..{h.max():.3f})")

    # 2. silhouette of the inflated shape equals the flat outline (near-orthographic front view)
    mvp, eye = orbit_camera(0, 0, 40.0, fov_deg=3.5)
    inflated = rasterize(torch.cat([xy, h[:, None]], 1), Fc, mvp, 256, 256) >= 0
    flat = rasterize(torch.cat([xy, torch.zeros_like(h)[:, None]], 1), Fc, mvp, 256, 256) >= 0
    iou = (inflated & flat).sum().item() / max((inflated | flat).sum().item(), 1)
    print(f"[2] front silhouette IoU inflated vs flat outline: {iou:.4f}")

    # 3. finite-difference check of the shading gradient w.r.t. the parameters
    mvp, eye = orbit_camera(30, 20, 3.2)
    W = torch.rand(3, 256, 256, device="cuda")

    def objective():
        img = render(torch.cat([xy, model.heights()[:, None]], 1), Fc, mvp, eye, 256, 256)
        return (img.double() * W.double()).sum()

    model.theta.data.normal_(0, 0.05)
    model.theta.grad = None
    objective().backward()
    grad = model.theta.grad.clone()
    idx = torch.topk(grad.abs(), 5).indices
    # global spectral modes move the whole sheet, and the square root makes
    # that strongly nonlinear near the outline, so they need a smaller step
    step = 1e-2 if args.rep == "full" else 1e-4
    for i in idx.tolist():
        with torch.no_grad():
            model.theta[i] += step
            fp = objective().item()
            model.theta[i] -= 2 * step
            fm = objective().item()
            model.theta[i] += step
        fd = (fp - fm) / (2 * step)
        print(f"[3] theta[{i}]: autograd {grad[i].item():+.4f}  finite diff {fd:+.4f}")

    # 4. landmark projection meets the targets
    model.theta.data.zero_()
    model.project_landmarks()
    if scene.landmarks:
        B = scene.landmark_matrix()
        h = model.heights().detach().cpu().numpy()
        for k, lm in enumerate(scene.landmarks):
            row = B[k].toarray().ravel()
            print(f"[4] landmark {k}: factor {lm.factor:.2f}, height {row @ scene.h0:+.4f} -> "
                  f"{row @ h:+.4f} (ratio {(row @ h) / (row @ scene.h0):.3f}), "
                  f"constraint error {model.landmark_error():.2e}")
    else:
        print("[4] no landmarks in the scene")

    # 5. front and side render for a visual orientation check
    tiles = []
    for azimuth, elevation in ((0, 0), (90, 0), (30, 20)):
        mvp, eye = orbit_camera(azimuth, elevation, 3.2)
        img = render(torch.cat([xy, model.heights().detach()[:, None]], 1), Fc, mvp, eye, 256, 256)
        tiles.append(img[0].cpu().numpy())
    Image.fromarray((np.concatenate(tiles, 1) * 255).astype(np.uint8)).save("check_views.png")
    print("[5] wrote check_views.png (front, side, three-quarter)")


if __name__ == "__main__":
    main()
