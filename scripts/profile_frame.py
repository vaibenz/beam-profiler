"""Measure a single beam image: D4 sigma second moments plus a Gaussian fit cross check.

    python scripts/profile_frame.py beam.png --dark dark.png --pixel-um 3.0
    python scripts/profile_frame.py --simulate --z 0
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from beamprofiler import gaussfit, moments, simulate  # noqa: E402


def load(path, channel):
    import cv2
    img = cv2.imread(str(path), cv2.IMREAD_UNCHANGED)
    if img is None:
        raise SystemExit(f"Could not read {path}")
    if img.ndim == 3:
        img = img[:, :, {"blue": 0, "green": 1, "red": 2}[channel]]
    return img.astype(np.float64)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("image", nargs="?")
    p.add_argument("--dark")
    p.add_argument("--channel", default="red")
    p.add_argument("--pixel-um", type=float, default=3.0)
    p.add_argument("--simulate", action="store_true")
    p.add_argument("--z", type=float, default=0.0, help="z position for the simulated frame, mm")
    p.add_argument("--plot", default=None, help="save a profile plot to this path")
    a = p.parse_args()

    if a.simulate:
        cam = simulate.CameraSpec(pixel_um=a.pixel_um)
        rng = np.random.default_rng(0)
        img = np.mean([simulate.frame(a.z, simulate.BeamSpec(), cam, rng=rng) for _ in range(8)], axis=0)
        dark = simulate.dark_frame(cam, rng=rng)
    else:
        if not a.image:
            raise SystemExit("Give an image path or --simulate")
        img = load(a.image, a.channel)
        dark = load(a.dark, a.channel) if a.dark else None

    corr = moments.subtract_background(img, dark)
    bw = moments.measure(corr)
    gf = gaussfit.fit(corr, guess=bw)
    s = a.pixel_um
    print(f"Centroid           ({bw.cx:.2f}, {bw.cy:.2f}) px")
    print(f"D4 sigma x, y      {bw.dx * s:.1f}, {bw.dy * s:.1f} um")
    print(f"Principal axes     {bw.d_major * s:.1f} x {bw.d_minor * s:.1f} um at {bw.phi_deg:.1f} deg")
    print(f"Ellipticity        {bw.ellipticity:.3f}")
    print(f"Gaussian 1/e^2     {gf.d_u * s:.1f} x {gf.d_v * s:.1f} um at {gf.theta_deg:.1f} deg")
    print(f"Fit RMS residual   {100 * gf.rms_residual:.2f} % of peak")
    agree = abs(max(gf.d_u, gf.d_v) / bw.d_major - 1)
    print(f"D4 sigma vs Gaussian disagreement {100 * agree:.1f} %"
          + ("  (beam is not a clean TEM00, trust D4 sigma)" if agree > 0.05 else ""))

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(9, 3.8))
        half = int(1.5 * max(bw.dx, bw.dy))
        cx, cy = int(round(bw.cx)), int(round(bw.cy))
        crop = corr[max(cy - half, 0):cy + half, max(cx - half, 0):cx + half]
        axs[0].imshow(crop, cmap="inferno")
        axs[0].set_title("Background corrected frame")
        axs[0].axis("off")
        xs = (np.arange(corr.shape[1]) - bw.cx) * s
        row = corr[cy]
        axs[1].plot(xs, row / row.max(), ".", ms=3, label="row through centroid")
        model = np.exp(-2 * xs ** 2 / (gf.d_u * s / 2) ** 2)
        axs[1].plot(xs, model, "-", label="Gaussian fit")
        axs[1].set_xlim(-1.5 * bw.dx * s, 1.5 * bw.dx * s)
        axs[1].set_xlabel("x (µm)")
        axs[1].set_ylabel("normalised intensity")
        axs[1].legend(frameon=False)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=150)
        print(f"Saved {a.plot}")


if __name__ == "__main__":
    main()
