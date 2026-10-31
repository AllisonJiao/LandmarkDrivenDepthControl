"""Differentiable inflation: pressure change -> heights.

On the flat mesh the inflation is linear, L u = M p with u = 0 on the outline,
so a pressure change dp moves the pre-sqrt field by du = G dp on the interior
vertices with G = L_II^-1 M_II. Heights are h = sgn(u) sqrt(|u|).

Representations of dp:
  full      one value per interior vertex
  spectral  dp = Phi c with the first K Dirichlet eigenmodes of
            -L_II phi = lambda M_II phi (M-orthonormal, so |Phi c|_M = |c|)
"""

import numpy as np
import scipy.linalg
import scipy.sparse.linalg
import torch


class PressureModel(torch.nn.Module):
    def __init__(self, scene, rep="full", k=20, eps=1e-6, device="cuda"):
        super().__init__()
        I = scene.interior
        L_II = scene.L[I][:, I]
        M_I = scene.M[I]

        # dense G = L_II^-1 M_II (fine for meshes up to several thousand vertices)
        G = scipy.linalg.solve(L_II.toarray(), np.diag(M_I), assume_a="sym")

        # pressure units: the median default pressure, so parameters are O(1)
        p0 = (scene.L @ scene.u0)[I] / M_I
        self.p_scale = float(np.median(np.abs(p0)))

        if rep == "full":
            G_eff = G * self.p_scale
            w_inv = 1.0 / M_I
            n_params = len(I)
            self.Phi = None
        elif rep == "spectral":
            k = min(k, len(I) - 2)
            evals, Phi = scipy.sparse.linalg.eigsh(
                -L_II, k=k, M=scipy.sparse.diags(M_I), sigma=0, which="LM")
            order = np.argsort(evals)
            Phi = Phi[:, order]
            G_eff = G @ Phi * self.p_scale
            w_inv = np.ones(k)
            n_params = k
            self.register_buffer("Phi", torch.tensor(Phi, dtype=torch.float32))
        else:
            raise ValueError(f"unknown representation {rep}")
        self.rep = rep

        self.register_buffer("G_eff", torch.tensor(G_eff, dtype=torch.float32))
        self.register_buffer("M_I", torch.tensor(M_I, dtype=torch.float32))
        self.register_buffer("u0_I", torch.tensor(scene.u0[I], dtype=torch.float32))
        self.register_buffer("h0", torch.tensor(scene.h0, dtype=torch.float32))
        self.register_buffer("I", torch.tensor(I, dtype=torch.long))
        self.theta = torch.nn.Parameter(torch.zeros(n_params))
        self.eps = eps

        # landmark constraints A theta = r on the parameters, with A = B_I G_eff
        # and r = (f^2 - 1) * (B u0): a factor scales the height, so u by f^2
        B = scene.landmark_matrix()
        factors = np.array([lm.factor for lm in scene.landmarks])
        r = (factors**2 - 1) * (B @ scene.u0)
        A = B[:, I] @ G_eff
        KK = (A * w_inv) @ A.T  # A W^-1 A' (k x k), constant during a run
        self.register_buffer("A", torch.tensor(A, dtype=torch.float32))
        self.register_buffer("r", torch.tensor(r, dtype=torch.float32))
        # inverse of the metric the projection is measured in (M for full
        # pressure, identity for M-orthonormal spectral coefficients)
        self.register_buffer("w_inv", torch.tensor(w_inv, dtype=torch.float32))
        self.register_buffer("KK", torch.tensor(KK, dtype=torch.float32))

        self.to(device)

    def dp(self):
        """Pressure change on the interior vertices, in units of p_scale."""
        return self.theta if self.Phi is None else self.Phi @ self.theta

    def heights(self):
        s = torch.sign(self.u0_I)
        u = self.u0_I + self.G_eff @ self.theta
        # keep each sheet on its own side of the drawing plane
        v = torch.clamp(s * u, min=self.eps)
        h = self.h0.clone()
        h[self.I] = s * torch.sqrt(v)
        return h

    def regularizer(self):
        """0.5 dp' M dp (formulation A), normalized by the total area."""
        if self.Phi is None:
            return 0.5 * (self.M_I * self.theta**2).sum() / self.M_I.sum()
        return 0.5 * (self.theta**2).sum() / self.M_I.sum()

    def landmark_error(self):
        if self.A.shape[0] == 0:
            return 0.0
        return float((self.A @ self.theta.detach() - self.r).abs().max())

    @torch.no_grad()
    def project_landmarks(self):
        """Project theta onto {A theta = r} in the metric W (W^-1 = w_inv):
        theta <- theta - W^-1 A' (A W^-1 A')^-1 (A theta - r)."""
        if self.A.shape[0] == 0:
            return
        residual = self.A @ self.theta - self.r  # (k,) constraint violation
        lam = torch.linalg.solve(self.KK, residual)  # (A W^-1 A')^-1 residual
        self.theta.sub_(self.w_inv * (self.A.T @ lam))  # in place: theta is the optimized Parameter
