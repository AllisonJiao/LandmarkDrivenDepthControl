"""Optimize the inflation pressure of a Monster Mash scene with score
distillation from Stable Diffusion.

  python run_sds.py --scene mm_sds_scene.json --prompt "a puffy cloud" --out runs/cloud

Writes result.json (load it with "Import SDS result" in the app), result.obj,
progress images and log.csv into --out.
"""

import argparse
import csv
import os
import random
import time

import numpy as np
import torch
from PIL import Image

from pressure import PressureModel
from raster import orbit_camera, render
from scene import Scene


def view_suffix(azimuth):
    a = (azimuth + 180) % 360 - 180
    if abs(a) <= 45:
        return "front view"
    if abs(a) >= 135:
        return "back view"
    return "side view"


def sample_view(rng, front_weight):
    elevation = rng.uniform(-10, 45)
    if rng.random() < front_weight:
        azimuth = rng.gauss(0, 50)
    else:
        azimuth = rng.uniform(-180, 180)
    return azimuth, elevation


def mesh_vertices(xy, h):
    return torch.cat([xy, h[:, None]], dim=1)


@torch.no_grad()
def save_views(path, xy, h, Fc, distance, res=256):
    tiles = []
    for azimuth in (0, 45, 90, 180):
        mvp, eye = orbit_camera(azimuth, 15, distance)
        img = render(mesh_vertices(xy, h), Fc, mvp, eye, res, res)
        tiles.append(img[0].cpu().numpy())
    grid = np.concatenate(tiles, axis=1)
    Image.fromarray((np.clip(grid, 0, 1) * 255).astype(np.uint8)).save(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scene", required=True)
    ap.add_argument("--prompt", default="")
    ap.add_argument("--negative", default="")
    ap.add_argument("--rep", choices=["full", "spectral"], default="full")
    ap.add_argument("--k", type=int, default=20, help="spectral modes")
    ap.add_argument("--iters", type=int, default=1500)
    ap.add_argument("--lr", type=float, default=1e-2)
    ap.add_argument("--reg", type=float, default=1e-3, help="weight of 0.5 dp'M dp")
    ap.add_argument("--cfg", type=float, default=100.0)
    ap.add_argument("--model", default="stable-diffusion-v1-5/stable-diffusion-v1-5")
    ap.add_argument("--res", type=int, default=512)
    ap.add_argument("--distance", type=float, default=3.2)
    ap.add_argument("--front-weight", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="runs/latest")
    ap.add_argument("--render-only", action="store_true",
                    help="only render the default inflation (no diffusion model)")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    torch.manual_seed(args.seed)
    rng = random.Random(args.seed)
    device = "cuda"

    scene = Scene(args.scene)
    model = PressureModel(scene, rep=args.rep, k=args.k, device=device)
    model.project_landmarks()
    xy = torch.tensor(scene.V[:, :2], dtype=torch.float32, device=device)
    Fc = torch.tensor(scene.F, dtype=torch.long, device=device)
    print(f"{scene.n} vertices, {len(scene.F)} faces, {len(scene.landmarks)} landmarks, "
          f"rep={args.rep}, params={model.theta.numel()}")

    save_views(os.path.join(args.out, "views_0000.png"), xy, model.heights().detach(),
               Fc, args.distance)
    if args.render_only:
        return

    from sds import StableDiffusionSDS
    guidance = StableDiffusionSDS(args.model, device)
    embeddings = {s: guidance.encode_text(f"{args.prompt}, {s}", args.negative)
                  for s in ("front view", "side view", "back view")}

    opt = torch.optim.Adam([model.theta], lr=args.lr)
    log = open(os.path.join(args.out, "log.csv"), "w", newline="")
    writer = csv.writer(log)
    writer.writerow(["iter", "t_max", "sds", "reg", "landmark_err", "ms"])

    for it in range(1, args.iters + 1):
        start = time.perf_counter()
        t_max = 0.98 - (0.98 - 0.5) * it / args.iters
        azimuth, elevation = sample_view(rng, args.front_weight)
        mvp, eye = orbit_camera(azimuth, elevation, args.distance)

        h = model.heights()
        img = render(mesh_vertices(xy, h), Fc, mvp, eye, args.res, args.res)
        sds = guidance.loss(img, embeddings[view_suffix(azimuth)], args.cfg, (0.02, t_max))
        reg = model.regularizer()
        loss = sds + args.reg * reg

        opt.zero_grad()
        loss.backward()
        opt.step()
        model.project_landmarks()

        ms = (time.perf_counter() - start) * 1000
        writer.writerow([it, f"{t_max:.3f}", f"{sds.item():.4f}", f"{reg.item():.6f}",
                         f"{model.landmark_error():.2e}", f"{ms:.1f}"])
        if it % 100 == 0 or it == args.iters:
            log.flush()
            save_views(os.path.join(args.out, f"views_{it:04d}.png"), xy,
                       model.heights().detach(), Fc, args.distance)
            print(f"iter {it}: sds {sds.item():.2f}, reg {reg.item():.4f}, {ms:.0f} ms/it")

    log.close()
    h = model.heights().detach().cpu().numpy()
    scene.write_result(os.path.join(args.out, "result.json"), h)
    scene.write_obj(os.path.join(args.out, "result.obj"), h)
    np.save(os.path.join(args.out, "theta.npy"), model.theta.detach().cpu().numpy())
    print(f"wrote {args.out}/result.json")


if __name__ == "__main__":
    main()
