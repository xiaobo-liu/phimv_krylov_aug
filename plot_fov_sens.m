function fig = plot_fov_sens(results, sens_name)
%PLOT_FOV_SENS Plot field-of-values sets for one right-hand-side sensitivity test.

lg_linewidth = 1.5;
lg_linewidth_bnd = 1.2;
lg_fontsize = 9;

axlabel_linewidth = 1.0;
axlabel_fontsize = 9;

color_K = [0.635 0.078 0.184];
color_K_bnd = [0.85 0.325 0.098];
color_MK = [0 0.4470 0.7410];
color_W = [0 0 0];
color_W_bnd = [0.5 0.5 0.5];
color_A = [0.13 0.55 0.13];

num_cases = numel(results);
num_cols = ceil(sqrt(num_cases));
num_rows = ceil(num_cases / num_cols);

fig = figure;
clf(fig)
set(fig, 'Color', 'w', 'Units', 'inches', ...
    'Position', [1, 1, 3.8*num_cols, 3.3*num_rows]);

tiledlayout(fig, num_rows, num_cols, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:num_cases
    result = results{k};

    nexttile
    hold on
    % Filled F(A) region, drawn first so the other field-of-values curves stay visible.
    fill(real(result.range_A), imag(result.range_A), color_A, ...
        'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(real(result.range_K), imag(result.range_K), '-', ...
        'Color', color_K, 'LineWidth', lg_linewidth);
    plot(real(result.range_MK), imag(result.range_MK), '-', ...
        'Color', color_MK, 'LineWidth', lg_linewidth);
    plot(real(result.range_W), imag(result.range_W), '-', ...
        'Color', color_W, 'LineWidth', lg_linewidth);
    plot(real(result.range_A), imag(result.range_A), '-', ...
        'Color', color_A, 'LineWidth', lg_linewidth);
    plot(real(result.bnd_K), imag(result.bnd_K), '--', ...
        'Color', color_K_bnd, 'LineWidth', lg_linewidth_bnd);
    plot(real(result.bnd_W), imag(result.bnd_W), '--', ...
        'Color', color_W_bnd, 'LineWidth', lg_linewidth_bnd);
    grid on
    xlabel('Real part', 'interpreter', 'latex');
    ylabel('Imaginary part', 'interpreter', 'latex');
    title(case_title(result, sens_name), 'interpreter', 'latex');
    if k == 1
        legend({'$\mathcal F(K)$', '$\mathcal F_{M_{\rm K}}(K)$', ...
            '$\mathcal F(W)$', '$\mathcal F(A)$', ...
            'bound on $\mathcal F(K)$', 'bound on $\mathcal F(W)$'}, ...
            'interpreter', 'latex', 'FontSize', lg_fontsize, ...
            'Location', 'northwest');
    end
    set(gca, 'linewidth', axlabel_linewidth)
    set(gca, 'fontsize', axlabel_fontsize)
end

end

function txt = case_title(result, sens_name)
%CASE_TITLE Return the panel title for one parameter value.

if strcmp(sens_name, 'delta')
    txt = sprintf('$\\delta=%g$', result.delta);
elseif strcmp(sens_name, 'beta')
    txt = sprintf('$\\beta=%g$', result.beta);
else
    txt = sprintf('$\\beta=%g$, $\\delta=%g$', ...
        result.beta, result.delta);
end

end
