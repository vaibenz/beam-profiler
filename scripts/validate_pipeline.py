"""Monte Carlo check of the whole measurement chain against known ground truth.

Runs many simulated through focus scans with independent camera noise, fits
each one, and reports the bias and spread of the recovered waist, waist
position and M^2. This is how I know how far to trust a single hardware scan.

    python scripts/validate_pipeline.py --runs 40
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from beamprofiler import caustic, moments, simulate  # noqa: E402


def one_scan(rng, beam, cam, z, frames, dark):
    dx, dy = [], []
    for zi in z:
        img = np.mean([simulate.frame(zi, beam, cam, rng=rng) for _ in range(frames)], axis=0)
        bw = moments.measure(moments.subtract_background(img, dark))
        dx.append(bw.dx * cam.pixel_um)
        dy.append(bw.dy * cam.pixel_um)
    return caustic.fit(z, dx, beam.wavelength_nm), caustic.fit(z, dy, beam.wavelength_nm)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--runs", type=int, default=40)
    p.add_argument("--frames", type=int, default=4)
    p.add_argument("--out", default="results/validation.json")
    a = p.parse_args()

    beam, cam = simulate.BeamSpec(), simulate.CameraSpec()
    rng = np.random.default_rng(42)
    dark = simulate.dark_frame(cam, rng=rng)
    zr_x = np.pi * (beam.w0x_um * 1e-3) ** 2 / (beam.m2x * beam.wavelength_nm * 1e-6)
    z = caustic.iso_scan_positions(0.3, zr_x, n=20)

    res = {k: [] for k in ("d0x", "d0y", "z0x", "z0y", "m2x", "m2y")}
    for i in range(a.runs):
        fx, fy = one_scan(rng, beam, cam, z, a.frames, dark)
        res["d0x"].append(fx.d0_um); res["d0y"].append(fy.d0_um)
        res["z0x"].append(fx.z0_mm); res["z0y"].append(fy.z0_mm)
        res["m2x"].append(fx.m2); res["m2y"].append(fy.m2)

    truth = {"d0x": 2 * beam.w0x_um, "d0y": 2 * beam.w0y_um, "z0x": beam.z0x_mm,
             "z0y": beam.z0y_mm, "m2x": beam.m2x, "m2y": beam.m2y}
    report = {}
    print(f"{a.runs} simulated ISO spaced scans, {len(z)} positions each, {a.frames} frames per position\n")
    print(f"{'quantity':10s}{'truth':>10s}{'mean':>10s}{'bias %':>9s}{'std':>9s}")
    for k, vals in res.items():
        v = np.array(vals)
        t = truth[k]
        bias = (v.mean() - t) / t * 100 if t else float("nan")
        report[k] = {"truth": t, "mean": float(v.mean()), "std": float(v.std(ddof=1)),
                     "bias_pct": None if t == 0 else float(bias)}
        bias_s = f"{bias:+8.2f}" if t else f"{'n/a':>8s}"
        print(f"{k:10s}{t:10.3f}{v.mean():10.3f}{bias_s} {v.std(ddof=1):8.3f}")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps({"runs": a.runs, "positions": len(z), "frames": a.frames,
                                       "results": report}, indent=2))
    print(f"\nSaved {a.out}")


if __name__ == "__main__":
    main()
