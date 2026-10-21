"""Beam caustic fit and M squared estimate following ISO 11146.

Along the propagation axis the second moment diameter of any beam obeys

    d(z)^2 = d0^2 + theta^2 (z - z0)^2

so d^2 is an exact quadratic a + b z + c z^2 in z. A linear least squares fit
of that quadratic gives

    z0    = -b / (2c)                  waist position
    d0    = sqrt(a - b^2 / (4c))       waist diameter
    theta = sqrt(c)                    full far field divergence angle
    zR    = d0 / theta                 Rayleigh length
    M^2   = pi d0 theta / (4 lambda)   beam quality factor (1 for ideal TEM00)
"""

from dataclasses import dataclass, asdict

import numpy as np


@dataclass
class CausticFit:
    d0_um: float
    z0_mm: float
    theta_mrad: float
    zr_mm: float
    m2: float
    r2: float
    n_points: int
    n_within_zr: int
    n_beyond_2zr: int
    iso_sampling_ok: bool

    def as_dict(self):
        return asdict(self)


def fit(z_mm, d_um, wavelength_nm):
    z = np.asarray(z_mm, dtype=np.float64)
    d = np.asarray(d_um, dtype=np.float64)
    if len(z) < 5:
        raise ValueError("Need at least 5 points for a caustic fit (ISO 11146 asks for 10)")

    # Fit d^2 in metres squared against z in metres to keep units honest.
    zm, dm = z * 1e-3, d * 1e-6
    c, b, a = np.polyfit(zm, dm ** 2, 2)
    if c <= 0:
        raise ValueError("Fitted curvature is not positive, the scan does not bracket a waist")
    z0 = -b / (2 * c)
    d0_sq = a - b ** 2 / (4 * c)
    if d0_sq <= 0:
        raise ValueError("Fitted waist diameter is imaginary, check the data")
    d0 = np.sqrt(d0_sq)
    theta = np.sqrt(c)
    zr = d0 / theta
    lam = wavelength_nm * 1e-9
    m2 = np.pi * d0 * theta / (4 * lam)

    pred = a + b * zm + c * zm ** 2
    ss_res = np.sum((dm ** 2 - pred) ** 2)
    ss_tot = np.sum((dm ** 2 - np.mean(dm ** 2)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    dist = np.abs(zm - z0)
    n_in = int(np.sum(dist <= zr))
    n_out = int(np.sum(dist >= 2 * zr))
    iso_ok = len(z) >= 10 and n_in >= 5 and n_out >= 5

    return CausticFit(
        d0_um=d0 * 1e6, z0_mm=z0 * 1e3, theta_mrad=theta * 1e3, zr_mm=zr * 1e3,
        m2=m2, r2=r2, n_points=len(z), n_within_zr=n_in, n_beyond_2zr=n_out,
        iso_sampling_ok=iso_ok,
    )


def model(z_mm, fit_result: CausticFit):
    """Diameter in micrometres predicted by a fit, for plotting."""
    z = np.asarray(z_mm, dtype=np.float64)
    dz_m = (z - fit_result.z0_mm) * 1e-3
    d0_m = fit_result.d0_um * 1e-6
    th = fit_result.theta_mrad * 1e-3
    return np.sqrt(d0_m ** 2 + th ** 2 * dz_m ** 2) * 1e6


def iso_scan_positions(z0_mm, zr_mm, n=16, span_zr=4.0):
    """Suggested z positions for a second pass: half within one zR, half beyond two zR."""
    half = n // 2
    inner = np.linspace(-zr_mm, zr_mm, half)
    outer_each = (n - half) // 2
    far = np.linspace(2 * zr_mm, span_zr * zr_mm, outer_each)
    outer = np.concatenate([-far[::-1], far])
    return np.sort(z0_mm + np.concatenate([inner, outer]))
