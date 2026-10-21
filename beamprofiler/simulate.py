"""Synthetic camera frames of a focused laser beam.

Used for the unit tests and for validating the whole pipeline end to end
before the hardware is built: frames are generated from a beam with known
M^2, waist and astigmatism, passed through a simple camera model (shot noise,
read noise, dark offset, 8 bit quantisation and saturation), and the
measured caustic is compared against the truth.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class BeamSpec:
    wavelength_nm: float = 650.0
    w0x_um: float = 60.0      # 1/e^2 waist radius, x
    w0y_um: float = 55.0      # 1/e^2 waist radius, y
    z0x_mm: float = 0.0       # waist position, x
    z0y_mm: float = 0.6       # waist position, y (offset gives astigmatism)
    m2x: float = 1.3
    m2y: float = 1.4
    tilt_deg: float = 0.0     # rotation of the beam axes on the sensor


@dataclass
class CameraSpec:
    width: int = 640
    height: int = 480
    pixel_um: float = 3.0
    full_well_e: float = 6000.0
    read_noise_e: float = 8.0
    dark_dn: float = 4.0
    bits: int = 8


def radius_um(z_mm, w0_um, z0_mm, m2, wavelength_nm):
    zr_mm = np.pi * (w0_um * 1e-3) ** 2 / (m2 * wavelength_nm * 1e-6)
    return w0_um * np.sqrt(1 + ((z_mm - z0_mm) / zr_mm) ** 2)


def ideal_image(z_mm, beam: BeamSpec, cam: CameraSpec, centre=None):
    """Noise free intensity on the sensor, normalised to a peak of 1."""
    wx = radius_um(z_mm, beam.w0x_um, beam.z0x_mm, beam.m2x, beam.wavelength_nm) / cam.pixel_um
    wy = radius_um(z_mm, beam.w0y_um, beam.z0y_mm, beam.m2y, beam.wavelength_nm) / cam.pixel_um
    cx, cy = centre if centre is not None else (cam.width / 2 + 3.3, cam.height / 2 - 2.1)
    ys, xs = np.mgrid[0:cam.height, 0:cam.width]
    t = np.radians(beam.tilt_deg)
    u = (xs - cx) * np.cos(t) + (ys - cy) * np.sin(t)
    v = -(xs - cx) * np.sin(t) + (ys - cy) * np.cos(t)
    return np.exp(-2 * u ** 2 / wx ** 2 - 2 * v ** 2 / wy ** 2), 2 * wx, 2 * wy


def frame(z_mm, beam: BeamSpec, cam: CameraSpec, peak_fraction=0.7, rng=None):
    """One noisy 8 bit frame. peak_fraction sets exposure relative to full scale."""
    rng = np.random.default_rng() if rng is None else rng
    ideal, _, _ = ideal_image(z_mm, beam, cam)
    max_dn = 2 ** cam.bits - 1
    gain = cam.full_well_e / max_dn
    electrons = ideal * peak_fraction * cam.full_well_e
    electrons = rng.poisson(electrons).astype(np.float64)
    electrons += rng.normal(0, cam.read_noise_e, electrons.shape)
    dn = electrons / gain + cam.dark_dn
    return np.clip(np.round(dn), 0, max_dn).astype(np.uint8)


def dark_frame(cam: CameraSpec, n=16, rng=None):
    rng = np.random.default_rng() if rng is None else rng
    max_dn = 2 ** cam.bits - 1
    gain = cam.full_well_e / max_dn
    frames = [
        np.clip(np.round(rng.normal(0, cam.read_noise_e, (cam.height, cam.width)) / gain + cam.dark_dn), 0, max_dn)
        for _ in range(n)
    ]
    return np.mean(frames, axis=0)
