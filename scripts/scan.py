"""Through focus scan: step the camera along the beam, measure D4 sigma at each z,
fit the caustic in x and y, and report waist, Rayleigh length and M^2.

Hardware:
    python scripts/scan.py --port /dev/ttyACM0 --camera 0 --z-start -40 --z-stop 40 --points 21

Simulation (no hardware, known ground truth):
    python scripts/scan.py --simulate
"""

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from beamprofiler import camera, caustic, moments, simulate, stage  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--simulate", action="store_true")
    p.add_argument("--port", default="/dev/ttyACM0")
    p.add_argument("--camera", type=int, default=0)
    p.add_argument("--exposure", type=float, default=None)
    p.add_argument("--gamma", type=float, default=1.0, help="from scripts/linearity.py")
    p.add_argument("--pixel-um", type=float, default=3.0, help="sensor pixel pitch")
    p.add_argument("--mm-per-step", type=float, default=1.25 / 2048)
    p.add_argument("--wavelength-nm", type=float, default=650.0)
    p.add_argument("--z-start", type=float, default=-40.0)
    p.add_argument("--z-stop", type=float, default=40.0)
    p.add_argument("--points", type=int, default=21)
    p.add_argument("--frames", type=int, default=8, help="frames averaged per position")
    p.add_argument("--settle-s", type=float, default=0.5)
    p.add_argument("--iso-refine", action="store_true",
                   help="after the coarse scan, rescan at ISO 11146 spaced positions")
    p.add_argument("--out", default="results/scan")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args()


def open_hardware(args):
    if args.simulate:
        st = stage.SimulatedStage(args.mm_per_step)
        src = camera.SimulatedSource(st, seed=args.seed, cam=simulate.CameraSpec(pixel_um=args.pixel_um))
    else:
        st = stage.SerialStage(args.port, args.mm_per_step)
        src = camera.WebcamSource(args.camera, exposure=args.exposure, gamma=args.gamma)
    return st, src


def take_dark(src, n):
    if isinstance(src, camera.SimulatedSource):
        src.shutter_open = False
        dark, _ = camera.average(src, n)
        src.shutter_open = True
        return dark
    input("Block the laser, then press Enter to record a dark frame...")
    dark, _ = camera.average(src, n)
    input("Unblock the laser, then press Enter to start the scan...")
    return dark


def scan(st, src, positions, dark, args):
    rows = []
    for z in positions:
        st.move_to_mm(z)
        if not args.simulate:
            time.sleep(args.settle_s)
        img, sat = camera.average(src, args.frames)
        bw = moments.measure(moments.subtract_background(img, dark))
        um = bw.scaled(args.pixel_um)
        row = dict(z_mm=st.position_mm(), dx_um=um["dx"], dy_um=um["dy"],
                   d_major_um=um["d_major"], d_minor_um=um["d_minor"],
                   phi_deg=bw.phi_deg, cx_px=bw.cx, cy_px=bw.cy,
                   peak=float(img.max()), saturated_fraction=sat, converged=int(bw.converged))
        rows.append(row)
        flag = "  SATURATED, lower exposure" if sat > 0 else ""
        print(f"z={row['z_mm']:+8.3f} mm  dx={row['dx_um']:7.1f} um  dy={row['dy_um']:7.1f} um{flag}")
    return rows


def main():
    args = parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    st, src = open_hardware(args)
    try:
        st.home()
        dark = take_dark(src, 16)
        positions = np.linspace(args.z_start, args.z_stop, args.points)
        rows = scan(st, src, positions, dark, args)
        z = np.array([r["z_mm"] for r in rows])
        fx = caustic.fit(z, [r["dx_um"] for r in rows], args.wavelength_nm)
        fy = caustic.fit(z, [r["dy_um"] for r in rows], args.wavelength_nm)

        if args.iso_refine:
            zr = min(fx.zr_mm, fy.zr_mm)
            z0 = 0.5 * (fx.z0_mm + fy.z0_mm)
            print(f"\nRefining around z0={z0:.2f} mm with zR={zr:.2f} mm")
            positions = caustic.iso_scan_positions(z0, zr, n=max(args.points, 16))
            rows = scan(st, src, positions, dark, args)
            z = np.array([r["z_mm"] for r in rows])
            fx = caustic.fit(z, [r["dx_um"] for r in rows], args.wavelength_nm)
            fy = caustic.fit(z, [r["dy_um"] for r in rows], args.wavelength_nm)
        st.release()
    finally:
        src.close()
        st.close()

    with open(out.with_suffix(".csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    summary = {"x": fx.as_dict(), "y": fy.as_dict(),
               "astigmatic_offset_mm": fy.z0_mm - fx.z0_mm,
               "waist_ellipticity": min(fx.d0_um, fy.d0_um) / max(fx.d0_um, fy.d0_um),
               "simulated": args.simulate}
    if args.simulate:
        b = simulate.BeamSpec()
        summary["truth"] = {"d0x_um": 2 * b.w0x_um, "d0y_um": 2 * b.w0y_um,
                            "z0x_mm": b.z0x_mm, "z0y_mm": b.z0y_mm, "m2x": b.m2x, "m2y": b.m2y}
    out.with_suffix(".json").write_text(json.dumps(summary, indent=2))

    for axis, f in (("x", fx), ("y", fy)):
        print(f"\n{axis}: d0={f.d0_um:.1f} um  z0={f.z0_mm:+.2f} mm  zR={f.zr_mm:.2f} mm  "
              f"theta={f.theta_mrad:.2f} mrad  M2={f.m2:.3f}  R2={f.r2:.4f}  "
              f"ISO sampling {'ok' if f.iso_sampling_ok else 'NOT met, rerun with --iso-refine'}")

    plot(rows, fx, fy, out.with_suffix(".png"), args)


def plot(rows, fx, fy, path, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    z = np.array([r["z_mm"] for r in rows])
    zz = np.linspace(z.min(), z.max(), 400)
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(z, [r["dx_um"] for r in rows], "o", color="#1f5fa8", label="x measured")
    ax.plot(zz, caustic.model(zz, fx), "-", color="#1f5fa8", label=f"x fit, M²={fx.m2:.2f}")
    ax.plot(z, [r["dy_um"] for r in rows], "s", color="#c0392b", label="y measured")
    ax.plot(zz, caustic.model(zz, fy), "-", color="#c0392b", label=f"y fit, M²={fy.m2:.2f}")
    ax.set_xlabel("z (mm)")
    ax.set_ylabel("D4σ diameter (µm)")
    ax.set_title("Beam caustic" + (" (simulated data)" if args.simulate else ""))
    ax.grid(alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"\nSaved {path}")


if __name__ == "__main__":
    main()
