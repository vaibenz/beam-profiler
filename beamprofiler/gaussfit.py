"""Rotated elliptical 2D Gaussian fit, used as a cross check on the second moments.

The model is I = A exp(-2 u^2 / wu^2 - 2 v^2 / wv^2) + B, where u and v are
the coordinates rotated by theta about the centre. wu and wv are 1/e^2 radii,
so the 1/e^2 diameters are 2 wu and 2 wv and should agree with D4 sigma for a
clean TEM00 beam. A large disagreement, or a large residual, means the beam
is not Gaussian (side lobes, clipping, higher order modes) and the D4 sigma
number is the one to trust.
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit


@dataclass
class GaussFit:
    amplitude: float
    cx: float
    cy: float
    d_u: float        # 1/e^2 diameter along the fitted u axis, pixels
    d_v: float        # 1/e^2 diameter along the fitted v axis, pixels
    theta_deg: float
    offset: float
    rms_residual: float   # RMS of fit residual as a fraction of amplitude


def _model(coords, a, cx, cy, wu, wv, theta, b):
    x, y = coords
    ct, st = np.cos(theta), np.sin(theta)
    u = (x - cx) * ct + (y - cy) * st
    v = -(x - cx) * st + (y - cy) * ct
    return a * np.exp(-2 * u ** 2 / wu ** 2 - 2 * v ** 2 / wv ** 2) + b


def fit(img, guess=None, roi_factor=3.0):
    """Fit the Gaussian model. `guess` is an optional BeamWidth from moments.measure."""
    img = np.asarray(img, dtype=np.float64)
    h, w = img.shape
    if guess is not None:
        cx0, cy0 = guess.cx, guess.cy
        wu0, wv0 = guess.dx / 2, guess.dy / 2
        half = roi_factor * max(guess.dx, guess.dy) / 2
        x0, x1 = max(int(cx0 - half), 0), min(int(cx0 + half) + 1, w)
        y0, y1 = max(int(cy0 - half), 0), min(int(cy0 + half) + 1, h)
    else:
        cy0, cx0 = np.unravel_index(np.argmax(img), img.shape)
        wu0 = wv0 = min(h, w) / 10
        x0, x1, y0, y1 = 0, w, 0, h
    sub = img[y0:y1, x0:x1]
    ys, xs = np.mgrid[y0:y1, x0:x1]
    p0 = [sub.max() - np.median(sub), cx0, cy0, wu0, wv0, 0.0, np.median(sub)]
    popt, _ = curve_fit(_model, (xs.ravel(), ys.ravel()), sub.ravel(), p0=p0, maxfev=20000)
    a, cx, cy, wu, wv, theta, b = popt
    resid = sub - _model((xs, ys), *popt)
    return GaussFit(
        amplitude=a, cx=cx, cy=cy, d_u=2 * abs(wu), d_v=2 * abs(wv),
        theta_deg=np.degrees(theta) % 180, offset=b,
        rms_residual=float(np.sqrt(np.mean(resid ** 2)) / a),
    )
