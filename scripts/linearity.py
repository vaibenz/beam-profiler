"""Measure the webcam's response curve so frames can be linearised.

Webcams apply a gamma curve, so a pixel value of 200 is not twice the light of
100. Second moments computed on gamma encoded data overestimate beam width.
This script holds the laser fixed, sweeps the exposure, and fits

    DN = k * exposure ^ (1 / gamma)

at the brightest unsaturated pixels. Pass the resulting gamma to scan.py.

    python scripts/linearity.py --camera 0 --exposures -10 -9 -8 -7 -6 -5
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from beamprofiler import camera  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--camera", type=int, default=0)
    p.add_argument("--exposures", type=float, nargs="+", required=True,
                   help="OpenCV exposure settings; on V4L2 these are often log2 seconds")
    p.add_argument("--log2", action="store_true", default=True,
                   help="treat exposure values as log2 seconds (default)")
    a = p.parse_args()

    times, levels = [], []
    for e in a.exposures:
        src = camera.WebcamSource(a.camera, exposure=e)
        frame, sat = camera.average(src, 8)
        src.close()
        mask = frame > 0.9 * frame.max()
        if sat > 0:
            print(f"exposure {e}: saturated, skipped")
            continue
        times.append(2.0 ** e if a.log2 else e)
        levels.append(float(frame[mask].mean()))
        print(f"exposure {e}: mean DN {levels[-1]:.1f}")

    slope, _ = np.polyfit(np.log(times), np.log(levels), 1)
    gamma = 1 / slope
    print(f"\nlog DN vs log exposure slope {slope:.3f}, so gamma = {gamma:.2f}")
    print("A linear sensor gives slope 1.0. Pass --gamma to scan.py.")


if __name__ == "__main__":
    main()
