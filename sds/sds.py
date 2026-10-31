"""Score distillation sampling (DreamFusion) with Stable Diffusion."""

import torch
import torch.nn.functional as F
from diffusers import DDPMScheduler, StableDiffusionPipeline


class StableDiffusionSDS:
    def __init__(self, model_id, device="cuda"):
        pipe = StableDiffusionPipeline.from_pretrained(
            model_id, torch_dtype=torch.float16, safety_checker=None,
            requires_safety_checker=False)
        pipe.to(device)
        self.device = device
        self.vae = pipe.vae
        self.unet = pipe.unet
        self.tokenizer = pipe.tokenizer
        self.text_encoder = pipe.text_encoder
        for m in (self.vae, self.unet, self.text_encoder):
            m.requires_grad_(False)
        self.scheduler = DDPMScheduler.from_pretrained(model_id, subfolder="scheduler")
        self.alphas_cumprod = self.scheduler.alphas_cumprod.to(device)
        self.num_train_timesteps = self.scheduler.config.num_train_timesteps

    @torch.no_grad()
    def encode_text(self, prompt, negative_prompt=""):
        """(2, 77, D) embeddings: unconditional, conditional."""
        tokens = self.tokenizer([negative_prompt, prompt], padding="max_length",
                                max_length=self.tokenizer.model_max_length,
                                truncation=True, return_tensors="pt")
        return self.text_encoder(tokens.input_ids.to(self.device))[0]

    def loss(self, image, text_emb, guidance_scale=100.0, t_range=(0.02, 0.98)):
        """SDS surrogate loss for one (3, H, W) image in [0, 1]; its gradient
        w.r.t. the image latents is w(t) (eps_hat - eps)."""
        x = F.interpolate(image[None], size=(512, 512), mode="bilinear", align_corners=False)
        x = (x * 2 - 1).to(torch.float16)
        latents = self.vae.encode(x).latent_dist.sample() * self.vae.config.scaling_factor

        lo = int(t_range[0] * self.num_train_timesteps)
        hi = int(t_range[1] * self.num_train_timesteps)
        t = torch.randint(lo, hi + 1, (1,), device=self.device)
        noise = torch.randn_like(latents)
        with torch.no_grad():
            noisy = self.scheduler.add_noise(latents, noise, t)
            eps = self.unet(torch.cat([noisy] * 2), torch.cat([t] * 2),
                            encoder_hidden_states=text_emb.to(torch.float16)).sample
            eps_uncond, eps_cond = eps.chunk(2)
            eps_hat = eps_uncond + guidance_scale * (eps_cond - eps_uncond)

        w = 1 - self.alphas_cumprod[t]
        grad = torch.nan_to_num(w * (eps_hat - noise)).float()
        target = (latents.float() - grad).detach()
        return 0.5 * F.mse_loss(latents.float(), target, reduction="sum")
