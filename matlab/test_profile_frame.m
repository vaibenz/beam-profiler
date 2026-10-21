% Sanity check for profile_frame.m on a synthetic elliptical Gaussian.
% Run from this folder:  test_profile_frame
[xs, ys] = meshgrid(0:299, 0:299);
img = 200 * exp(-2 * (xs - 150.4) .^ 2 / 25 ^ 2 - 2 * (ys - 140.8) .^ 2 / 18 ^ 2);
out = profile_frame(img, 3.0);
assert(abs(out.dx_um - 150) / 150 < 2e-3, 'D4 sigma x wrong');
assert(abs(out.dy_um - 108) / 108 < 2e-3, 'D4 sigma y wrong');
assert(abs(out.gauss_du_um - 150) / 150 < 5e-3, 'Gaussian u wrong');
disp('profile_frame.m checks passed');
