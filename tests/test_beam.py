import numpy as np
import pytest

from beamprofiler import caustic, gaussfit, moments, simulate


def gaussian(shape, cx, cy, wx, wy, theta_deg=0.0):
    ys, xs = np.mgrid[0:shape[0], 0:shape[1]]
    t = np.radians(theta_deg)
    u = (xs - cx) * np.cos(t) + (ys - cy) * np.sin(t)
    v = -(xs - cx) * np.sin(t) + (ys - cy) * np.cos(t)
    return np.exp(-2 * u ** 2 / wx ** 2 - 2 * v ** 2 / wy ** 2)


def test_d4sigma_equals_1_over_e2_diameter_for_tem00():
    img = gaussian((400, 400), 200.3, 190.7, 30.0, 30.0)
    bw = moments.measure(img)
    assert bw.dx == pytest.approx(60.0, rel=1e-3)
    assert bw.dy == pytest.approx(60.0, rel=1e-3)
    assert bw.cx == pytest.approx(200.3, abs=1e-3)
    assert bw.converged


def test_principal_axes_of_rotated_ellipse():
    img = gaussian((500, 500), 250, 250, 40.0, 20.0, theta_deg=30)
    bw = moments.measure(img)
    assert bw.d_major == pytest.approx(80.0, rel=2e-3)
    assert bw.d_minor == pytest.approx(40.0, rel=2e-3)
    assert bw.phi_deg == pytest.approx(30.0, abs=0.2)
    assert bw.ellipticity == pytest.approx(0.5, rel=5e-3)


def test_background_offset_is_removed():
    img = gaussian((300, 300), 150, 150, 20, 20) * 200 + 12.0
    bw = moments.measure(moments.subtract_background(img))
    assert bw.dx == pytest.approx(40.0, rel=1e-2)


def test_gaussian_fit_agrees_with_moments():
    img = gaussian((300, 300), 140, 160, 25, 18, theta_deg=10) * 180 + 3
    corr = moments.subtract_background(img)
    gf = gaussfit.fit(corr, guess=moments.measure(corr))
    assert gf.d_u == pytest.approx(50.0, rel=1e-3)
    assert gf.d_v == pytest.approx(36.0, rel=1e-3)
    assert gf.rms_residual < 1e-4


def test_caustic_recovers_known_beam_exactly():
    lam, w0, m2, z0 = 650.0, 60.0, 1.3, 1.5
    z = np.linspace(-40, 40, 21)
    d = 2 * simulate.radius_um(z, w0, z0, m2, lam)
    f = caustic.fit(z, d, lam)
    assert f.d0_um == pytest.approx(2 * w0, rel=1e-6)
    assert f.z0_mm == pytest.approx(z0, abs=1e-6)
    assert f.m2 == pytest.approx(m2, rel=1e-6)
    assert f.r2 == pytest.approx(1.0, abs=1e-9)


def test_caustic_rejects_scan_without_waist():
    z = np.linspace(10, 40, 10)
    d = 50 + 2 * z  # monotonic, not a hyperbola
    with pytest.raises(ValueError):
        caustic.fit(z, d - 0.01 * z ** 2, 650.0)


def test_iso_positions_split_inside_and_outside_zr():
    pos = caustic.iso_scan_positions(z0_mm=2.0, zr_mm=5.0, n=16)
    dist = np.abs(pos - 2.0)
    assert np.sum(dist <= 5.0) >= 8
    assert np.sum(dist >= 10.0) >= 8


def test_end_to_end_simulated_scan_recovers_m2_within_5_percent():
    beam, cam = simulate.BeamSpec(), simulate.CameraSpec()
    rng = np.random.default_rng(7)
    dark = simulate.dark_frame(cam, rng=rng)
    zr = np.pi * (beam.w0x_um * 1e-3) ** 2 / (beam.m2x * beam.wavelength_nm * 1e-6)
    z = caustic.iso_scan_positions(0.0, zr, n=16)
    dx = []
    for zi in z:
        img = np.mean([simulate.frame(zi, beam, cam, rng=rng) for _ in range(4)], axis=0)
        dx.append(moments.measure(moments.subtract_background(img, dark)).dx * cam.pixel_um)
    f = caustic.fit(z, dx, beam.wavelength_nm)
    assert f.m2 == pytest.approx(beam.m2x, rel=0.05)
    assert f.d0_um == pytest.approx(2 * beam.w0x_um, rel=0.05)
