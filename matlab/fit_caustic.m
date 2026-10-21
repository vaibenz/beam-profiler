function result = fit_caustic(csv_path, wavelength_nm)
%FIT_CAUSTIC ISO 11146 caustic fit of a scan CSV written by scripts/scan.py.
%   result = fit_caustic('results/scan.csv', 650)
%
%   d(z)^2 = d0^2 + theta^2 (z - z0)^2 is quadratic in z, so polyfit on d^2
%   gives the waist, its position, the divergence and M^2 directly.
%   Uses base MATLAB only (no toolboxes) and also runs in GNU Octave.

    if nargin < 2, wavelength_nm = 650; end
    data = csvread(csv_path, 1, 0);         % columns: z_mm, dx_um, dy_um, ...
    z = data(:, 1);
    lam = wavelength_nm * 1e-9;

    labels = {'x', 'y'};
    result = struct();
    figure('Visible', 'off'); hold on; grid on;
    colours = {[0.12 0.37 0.66], [0.75 0.22 0.17]};
    zz = linspace(min(z), max(z), 400)';
    for k = 1:2
        d = data(:, k + 1);
        p = polyfit(z * 1e-3, (d * 1e-6) .^ 2, 2);   % p = [c b a]
        c = p(1); b = p(2); a = p(3);
        z0 = -b / (2 * c);
        d0 = sqrt(a - b ^ 2 / (4 * c));
        theta = sqrt(c);
        m2 = pi * d0 * theta / (4 * lam);
        r.d0_um = d0 * 1e6;
        r.z0_mm = z0 * 1e3;
        r.theta_mrad = theta * 1e3;
        r.zr_mm = d0 / theta * 1e3;
        r.m2 = m2;
        result.(labels{k}) = r;
        fprintf('%s: d0 = %.1f um, z0 = %+.2f mm, zR = %.2f mm, theta = %.2f mrad, M2 = %.3f\n', ...
            labels{k}, r.d0_um, r.z0_mm, r.zr_mm, r.theta_mrad, r.m2);

        model = sqrt(d0 ^ 2 + theta ^ 2 * (zz * 1e-3 - z0) .^ 2) * 1e6;
        plot(z, d, 'o', 'Color', colours{k});
        plot(zz, model, '-', 'Color', colours{k});
    end
    xlabel('z (mm)'); ylabel('D4\sigma diameter (\mum)');
    title('Beam caustic, MATLAB fit');
    legend({'x measured', 'x fit', 'y measured', 'y fit'}, 'Location', 'north');
    [folder, name] = fileparts(csv_path);
    out = fullfile(folder, [name '_matlab.png']);
    print(out, '-dpng', '-r150');
    fprintf('Saved %s\n', out);
end
