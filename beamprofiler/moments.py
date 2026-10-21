"""Second moment beam width measurement following ISO 11146.

For an intensity image I(x, y) the first moments give the centroid and the
second moments give the beam variances. The D4 sigma diameter is four times
the standard deviation of the intensity distribution, which equals the 1/e^2
diameter for an ideal TEM00 Gaussian beam.

ISO 11146 asks for the integration to be restricted to a rectangular aperture
three times the beam diameter, centred on the centroid. Because the diameter
is not known in advance, the aperture is refined iteratively until the result
stops changing.
"""

from dataclasses import dataclass, asdict

import numpy as np


@dataclass
class BeamWidth:
    cx: float          # centroid x, pixels
    cy: float          # centroid y, pixels
    dx: float          # D4 sigma diameter along image x, pixels
    dy: float          # D4 sigma diameter along image y, pixels
    d_major: float     # D4 sigma diameter along the major principal axis, pixels
    d_minor: float     # D4 sigma diameter along the minor principal axis, pixels
    phi_deg: float     # angle of the major axis from image x, degrees
    ellipticity: float # d_minor / d_major, 1.0 is round
    iterations: int
    converged: bool

    def scaled(self, pixel_um: float) -> dict:
        """Return the widths converted to micrometres."""
        d = asdict(self)
        for k in ("cx", "cy", "dx", "dy", "d_major", "d_minor"):
            d[k] = d[k] * pixel_um
        return d


def subtract_background(frame, dark=None, border=10):
    """Remove the dark frame, then any residual offset estimated from the image border.

    Negative values are kept on purpose. Clipping them to zero biases the
    second moments upwards because noise in the wings only adds energy.
    """
    img = np.asarray(frame, dtype=np.float64)
    if dark is not None:
        img = img - np.asarray(dark, dtype=np.float64)
    if border > 0:
        edge = np.concatenate([
            img[:border].ravel(), img[-border:].ravel(),
            img[:, :border].ravel(), img[:, -border:].ravel(),
        ])
        img = img - np.median(edge)
    return img


def _moments(img, x0, x1, y0, y1):
    sub = img[y0:y1, x0:x1]
    ys, xs = np.mgrid[y0:y1, x0:x1]
    total = sub.sum()
    if total <= 0:
        raise ValueError("No signal inside the integration window")
    cx = (sub * xs).sum() / total
    cy = (sub * ys).sum() / total
    sxx = (sub * (xs - cx) ** 2).sum() / total
    syy = (sub * (ys - cy) ** 2).sum() / total
    sxy = (sub * (xs - cx) * (ys - cy)).sum() / total
    return cx, cy, sxx, syy, sxy


def _diameters(sxx, syy, sxy):
    root = np.sqrt((sxx - syy) ** 2 + 4 * sxy ** 2)
    d_major = 2 * np.sqrt(2) * np.sqrt(max(sxx + syy + root, 0.0))
    d_minor = 2 * np.sqrt(2) * np.sqrt(max(sxx + syy - root, 0.0))
    phi = 0.5 * np.degrees(np.arctan2(2 * sxy, sxx - syy))
    return d_major, d_minor, phi


def measure(img, aperture_factor=3.0, max_iter=30, tol=1e-3):
    """Iterative ISO 11146 second moment measurement on a background corrected image."""
    img = np.asarray(img, dtype=np.float64)
    h, w = img.shape
    x0, x1, y0, y1 = 0, w, 0, h
    prev = None
    converged = False
    for it in range(1, max_iter + 1):
        cx, cy, sxx, syy, sxy = _moments(img, x0, x1, y0, y1)
        dx, dy = 4 * np.sqrt(max(sxx, 0)), 4 * np.sqrt(max(syy, 0))
        if prev is not None and abs(dx - prev[0]) <= tol * prev[0] and abs(dy - prev[1]) <= tol * prev[1]:
            converged = True
            break
        prev = (dx, dy)
        hx, hy = aperture_factor * dx / 2, aperture_factor * dy / 2
        x0, x1 = max(int(np.floor(cx - hx)), 0), min(int(np.ceil(cx + hx)) + 1, w)
        y0, y1 = max(int(np.floor(cy - hy)), 0), min(int(np.ceil(cy + hy)) + 1, h)
    d_major, d_minor, phi = _diameters(sxx, syy, sxy)
    return BeamWidth(
        cx=cx, cy=cy, dx=dx, dy=dy, d_major=d_major, d_minor=d_minor,
        phi_deg=phi, ellipticity=d_minor / d_major if d_major > 0 else float("nan"),
        iterations=it, converged=converged,
    )
