"""Pure-PyTorch triangle renderer (no CUDA compilation needed).

Pass 1 (no gradient) finds, for every pixel, the nearest triangle with a
z-buffer built by scatter_reduce. Pass 2 (differentiable) recomputes the
barycentric coordinates of the pixel centres from the current vertices and
shades interpolated vertex normals. Gradients therefore reach the vertices
through shading and interior positions, but not through silhouette edges,
which is enough here because the outline of an inflated sketch is fixed.
"""

import math

import torch
import torch.nn.functional as F


def look_at(eye, target, up):
    f = F.normalize(target - eye, dim=0)
    r = F.normalize(torch.linalg.cross(f, up), dim=0)
    u = torch.linalg.cross(r, f)
    view = torch.eye(4, device=eye.device)
    view[0, :3], view[1, :3], view[2, :3] = r, u, -f
    view[:3, 3] = -view[:3, :3] @ eye
    return view


def perspective(fov_deg, near=0.1, far=100.0, device="cuda"):
    t = 1.0 / math.tan(math.radians(fov_deg) / 2)
    P = torch.zeros(4, 4, device=device)
    P[0, 0], P[1, 1] = t, t
    P[2, 2], P[2, 3] = (far + near) / (near - far), 2 * far * near / (near - far)
    P[3, 2] = -1.0
    return P


def orbit_camera(azimuth_deg, elevation_deg, distance, fov_deg=40.0, device="cuda"):
    """Camera around the origin. The sketch's y axis points down (image
    coordinates) and its front sheet has z < 0, so azimuth 0 looks at the
    front along +z with -y as up; positive elevation looks from above."""
    az, el = math.radians(azimuth_deg), math.radians(elevation_deg)
    eye = torch.tensor([
        distance * math.cos(el) * math.sin(az),
        -distance * math.sin(el),
        -distance * math.cos(el) * math.cos(az),
    ], device=device)
    view = look_at(eye, torch.zeros(3, device=device),
                   torch.tensor([0.0, -1.0, 0.0], device=device))
    return perspective(fov_deg, device=device) @ view, eye


def vertex_normals(V, Fc):
    v0, v1, v2 = V[Fc[:, 0]], V[Fc[:, 1]], V[Fc[:, 2]]
    fn = torch.linalg.cross(v1 - v0, v2 - v0)  # area weighted
    N = torch.zeros_like(V)
    for j in range(3):
        N = N.index_add(0, Fc[:, j], fn)
    return F.normalize(N, dim=1)


def _to_screen(V, mvp, H, W):
    Vh = torch.cat([V, torch.ones_like(V[:, :1])], dim=1) @ mvp.T
    ndc = Vh[:, :3] / Vh[:, 3:4]
    x = (ndc[:, 0] * 0.5 + 0.5) * W
    y = (1.0 - (ndc[:, 1] * 0.5 + 0.5)) * H
    return x, y, ndc[:, 2]


def _edge(ax, ay, bx, by, px, py):
    return (px - ax) * (by - ay) - (py - ay) * (bx - ax)


@torch.no_grad()
def rasterize(V, Fc, mvp, H, W):
    """Triangle id per pixel (-1 for background), as a flat (H*W,) tensor."""
    x, y, z = _to_screen(V, mvp, H, W)
    fx, fy, fz = x[Fc], y[Fc], z[Fc]
    xmin = torch.clamp(torch.floor(fx.min(1).values - 0.5), 0, W - 1).long()
    xmax = torch.clamp(torch.ceil(fx.max(1).values - 0.5), 0, W - 1).long()
    ymin = torch.clamp(torch.floor(fy.min(1).values - 0.5), 0, H - 1).long()
    ymax = torch.clamp(torch.ceil(fy.max(1).values - 0.5), 0, H - 1).long()
    nx, ny = xmax - xmin + 1, ymax - ymin + 1
    counts = nx * ny
    face = torch.repeat_interleave(torch.arange(len(Fc), device=V.device), counts)
    start = torch.cumsum(counts, 0) - counts
    local = torch.arange(len(face), device=V.device) - start[face]
    px = xmin[face] + local % nx[face]
    py = ymin[face] + local // nx[face]
    cx, cy = px.float() + 0.5, py.float() + 0.5

    ax, ay, bx, by, qx, qy = (fx[face, 0], fy[face, 0], fx[face, 1], fy[face, 1],
                              fx[face, 2], fy[face, 2])
    area = _edge(ax, ay, bx, by, qx, qy)
    w0 = _edge(bx, by, qx, qy, cx, cy) / area
    w1 = _edge(qx, qy, ax, ay, cx, cy) / area
    w2 = 1 - w0 - w1
    inside = (w0 >= 0) & (w1 >= 0) & (w2 >= 0) & (area.abs() > 1e-12)
    depth = w0 * fz[face, 0] + w1 * fz[face, 1] + w2 * fz[face, 2]

    pix = (py * W + px)[inside]
    depth, face = depth[inside], face[inside]
    zbuf = torch.full((H * W,), float("inf"), device=V.device)
    zbuf.scatter_reduce_(0, pix, depth, reduce="amin")
    front = depth <= zbuf[pix]
    tri = torch.full((H * W,), -1, dtype=torch.long, device=V.device)
    tri[pix[front]] = face[front]
    return tri


def key_light(eye, up_amount=0.6, side_amount=0.3):
    """Light from above and to the side of the camera, so that surfaces facing
    the viewer are not shaded as flat white."""
    d = F.normalize(eye, dim=0)
    up = torch.tensor([0.0, -1.0, 0.0], device=eye.device)
    side = F.normalize(torch.linalg.cross(up, d), dim=0)
    return F.normalize(d + up_amount * up + side_amount * side, dim=0)


def render(V, Fc, mvp, eye, H=512, W=512, light_dir=None, ambient=0.25,
           albedo=0.85, background=1.0):
    """Shaded grey image (3, H, W) in [0, 1], differentiable w.r.t. V. The
    surface is at most `albedo` bright, so it stays visible on a white
    background."""
    tri = rasterize(V, Fc, mvp, H, W)
    covered = torch.nonzero(tri >= 0).squeeze(1)
    t = tri[covered]
    px = (covered % W).float() + 0.5
    py = (covered // W).float() + 0.5

    x, y, _ = _to_screen(V, mvp, H, W)
    f = Fc[t]
    ax, ay, bx, by, qx, qy = x[f[:, 0]], y[f[:, 0]], x[f[:, 1]], y[f[:, 1]], x[f[:, 2]], y[f[:, 2]]
    area = _edge(ax, ay, bx, by, qx, qy)
    w0 = _edge(bx, by, qx, qy, px, py) / area
    w1 = _edge(qx, qy, ax, ay, px, py) / area
    w2 = 1 - w0 - w1

    N = vertex_normals(V, Fc)
    n = F.normalize(w0[:, None] * N[f[:, 0]] + w1[:, None] * N[f[:, 1]] + w2[:, None] * N[f[:, 2]], dim=1)
    P = w0[:, None] * V[f[:, 0]] + w1[:, None] * V[f[:, 1]] + w2[:, None] * V[f[:, 2]]
    to_eye = F.normalize(eye[None, :] - P, dim=1)
    # two-sided: always shade the side facing the camera
    n = n * torch.sign((n * to_eye).sum(1, keepdim=True) + 1e-12)

    if light_dir is None:
        light_dir = key_light(eye)
    diffuse = torch.clamp((n * light_dir[None, :]).sum(1), min=0)
    shade = albedo * (ambient + (1 - ambient) * diffuse)

    img = torch.full((H * W,), background, device=V.device)
    img = img.index_put((covered,), shade)
    return img.view(1, H, W).expand(3, H, W)
