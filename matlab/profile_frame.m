function out = profile_frame(img, pixel_um)
%PROFILE_FRAME D4 sigma second moments and a 2D Gaussian fit on one beam image.
%   out = profile_frame(double(imread('beam.png')), 3.0)
%   Pass a background corrected, single channel image. Base MATLAB only
%   (fminsearch), so it also runs in GNU Octave.

    if nargin < 2, pixel_um = 1; end
    img = double(img);
    [h, w] = size(img);
    [xs, ys] = meshgrid(0:w-1, 0:h-1);

    % Iterative ISO 11146 aperture: 3 x D4 sigma around the centroid.
    x0 = 1; x1 = w; y0 = 1; y1 = h; prev = [inf inf];
    for it = 1:30
        sub = img(y0:y1, x0:x1); X = xs(y0:y1, x0:x1); Y = ys(y0:y1, x0:x1);
        s = sum(sub(:));
        cx = sum(sub(:) .* X(:)) / s;
        cy = sum(sub(:) .* Y(:)) / s;
        sxx = sum(sub(:) .* (X(:) - cx) .^ 2) / s;
        syy = sum(sub(:) .* (Y(:) - cy) .^ 2) / s;
        sxy = sum(sub(:) .* (X(:) - cx) .* (Y(:) - cy)) / s;
        d = 4 * sqrt([sxx syy]);
        if all(abs(d - prev) <= 1e-3 * prev), break; end
        prev = d;
        x0 = max(floor(cx - 1.5 * d(1)) + 1, 1); x1 = min(ceil(cx + 1.5 * d(1)) + 1, w);
        y0 = max(floor(cy - 1.5 * d(2)) + 1, 1); y1 = min(ceil(cy + 1.5 * d(2)) + 1, h);
    end
    root = sqrt((sxx - syy) ^ 2 + 4 * sxy ^ 2);
    out.cx = cx; out.cy = cy;
    out.dx_um = d(1) * pixel_um; out.dy_um = d(2) * pixel_um;
    out.d_major_um = 2 * sqrt(2) * sqrt(sxx + syy + root) * pixel_um;
    out.d_minor_um = 2 * sqrt(2) * sqrt(sxx + syy - root) * pixel_um;
    out.phi_deg = 0.5 * atan2d(2 * sxy, sxx - syy);

    % Gaussian cross check on the converged aperture.
    sub = img(y0:y1, x0:x1); X = xs(y0:y1, x0:x1); Y = ys(y0:y1, x0:x1);
    model = @(p) p(1) * exp(-2 * ((X - p(2)) * cos(p(6)) + (Y - p(3)) * sin(p(6))) .^ 2 / p(4) ^ 2 ...
                           - 2 * (-(X - p(2)) * sin(p(6)) + (Y - p(3)) * cos(p(6))) .^ 2 / p(5) ^ 2) + p(7);
    cost = @(p) sum(sum((sub - model(p)) .^ 2));
    p0 = [max(sub(:)), cx, cy, d(1) / 2, d(2) / 2, 0, 0];
    opts = optimset('MaxFunEvals', 2e4, 'MaxIter', 2e4, 'TolX', 1e-8, 'TolFun', 1e-10);
    p = fminsearch(cost, p0, opts);
    out.gauss_du_um = 2 * abs(p(4)) * pixel_um;
    out.gauss_dv_um = 2 * abs(p(5)) * pixel_um;
    out.gauss_rms_residual = sqrt(mean((sub(:) - reshape(model(p), [], 1)) .^ 2)) / p(1);

    fprintf('D4 sigma: %.1f x %.1f um   Gaussian 1/e^2: %.1f x %.1f um   residual %.2f %%\n', ...
        out.dx_um, out.dy_um, out.gauss_du_um, out.gauss_dv_um, 100 * out.gauss_rms_residual);
end
