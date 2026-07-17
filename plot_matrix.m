function figs = plot_matrix(result)
%PLOT_MATRIX Plot field-of-values regions and convergence curves.

lg_linewidth = 1.6;
lg_linewidth_bnd = 1.2;
lg_markersize = 5;
lg_fontsize = 11;

axlabel_linewidth = 1.0;
axlabel_fontsize = 10;

color_K = [0.635 0.078 0.184];
color_K_bnd = [0.85 0.325 0.098];
color_MK = [0 0.4470 0.7410];
color_W = [0 0 0];
color_W_bnd = [0.5 0.5 0.5];
color_A = [0.13 0.55 0.13];

%% Field-of-values comparison

figs.fov = figure;
clf(figs.fov)
set(figs.fov, 'Color', 'w', 'Units', 'inches', 'Position', [1,1,5.5,4.5]);
hold on

% Filled F(A) region, drawn first so the other field-of-values curves stay visible.
fill(real(result.range_A), imag(result.range_A), color_A, ...
    'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(real(result.range_K), imag(result.range_K), '-', 'Color', color_K, ...
    'LineWidth', lg_linewidth);
plot(real(result.range_MK), imag(result.range_MK), '-', 'Color', color_MK, ...
    'LineWidth', lg_linewidth);
plot(real(result.range_W), imag(result.range_W), '-', 'Color', color_W, ...
    'LineWidth', lg_linewidth);
plot(real(result.range_A), imag(result.range_A), '-', 'Color', color_A, ...
    'LineWidth', lg_linewidth);
plot(real(result.bnd_K), imag(result.bnd_K), '--', 'Color', color_K_bnd, ...
    'LineWidth', lg_linewidth_bnd);
plot(real(result.bnd_W), imag(result.bnd_W), '--', 'Color', color_W_bnd, ...
    'LineWidth', lg_linewidth_bnd);
plot(real(result.eig_K), imag(result.eig_K), 'x', 'Color', color_K, ...
    'LineWidth', lg_linewidth, 'MarkerSize', 4); % Eigenvalues of the KIOPS matrix K.

% axis equal
grid on
box on
xlabel('Real part', 'Interpreter', 'latex');
ylabel('Imaginary part', 'Interpreter', 'latex');
ax = gca;
ax.LineWidth = axlabel_linewidth;
ax.FontSize = axlabel_fontsize;
% title(sprintf('%s, $n=%d$, $s=%d$', result.mat_label, result.n, result.s), 'interpreter', 'latex');
lgd = legend({'$\mathcal F(K)$', '$\mathcal F_{M_{\rm K}}(K)$', '$\mathcal F(W)$', ...
    '$\mathcal F(A)$', 'bound on $\mathcal F(K)$', ...
    'bound on $\mathcal F(W)$', ...
    '$\lambda(K)$'}, ...
    'Interpreter', 'latex', 'FontSize', lg_fontsize, 'Location', 'northwest');
drawnow;
increase_legend_width(lgd, 1.1);

%% error and estimate comparison

figs.error = figure;
clf(figs.error)
set(figs.error, 'Color', 'w', 'Units', 'inches', 'Position', [1,1,5.5,4.5]);
hold on

ga(1) = semilogy(result.mvals, result.err_kiops, 'o-', 'Color', color_K, ...
    'LineWidth', lg_linewidth, 'MarkerSize', lg_markersize);
ga(2) = semilogy(result.mvals, result.err_W, 's-', 'Color', color_W, ...
    'LineWidth', lg_linewidth, 'MarkerSize', lg_markersize);
ga(3) = semilogy(result.mvals, result.err_X, '^-', 'Color', color_MK, ...
    'LineWidth', lg_linewidth, 'MarkerSize', lg_markersize);
ga(4) = semilogy(result.mvals, result.est_K, '--', 'Color', color_K, ...
    'LineWidth', lg_linewidth_bnd);
ga(5) = semilogy(result.mvals, result.est_MK, '-.', 'Color', color_MK, ...
    'LineWidth', lg_linewidth_bnd);
ga(6) = semilogy(result.mvals, result.est_W, '-.', 'Color', color_W, ...
    'LineWidth', lg_linewidth_bnd);

% Error panel in logarithmic scale.
set(gca, 'YScale', 'log')
grid on
box on
xlabel('Krylov dimension $m$', 'Interpreter', 'latex');
ylabel('Relative error and field-of-values estimate', 'Interpreter', 'latex');
ax = gca;
ax.LineWidth = axlabel_linewidth;
ax.FontSize = axlabel_fontsize;
% title(sprintf('%s, $n=%d$, $s=%d$', result.mat_label, result.n, result.s), 'interpreter', 'latex');
set_m_axis(result.mvals)

positive_vals = [result.err_kiops(:); result.err_W(:); result.err_X(:); ...
    result.est_K(:); result.est_MK(:); result.est_W(:)];
positive_vals = positive_vals(isfinite(positive_vals) & positive_vals > 0);
if ~isempty(positive_vals)
    ylim([max(1e-16, min(positive_vals)/5), max(positive_vals)*5])
end

lgd = legend(ga, {'KIOPS basis ($K$)', 'Block formulation ($W$)', ...
    'Orthonormal basis ($K_{\rm X}$)', 'bound based on $\mathcal F(K)$', ...
    'bound based on $\mathcal F_{M_{\rm K}}(K)$', ...
    'bound based on $\mathcal F(W)$'}, ...
    'Interpreter', 'latex', 'FontSize', lg_fontsize, 'Location', 'southwest');
drawnow;
increase_legend_width(lgd, 1.1);

%% Illustration of Proposition 5.5 for Poisson example
if isfield(result,'est_cheb')
    figs.error_sym = figure;
    clf(figs.error_sym)
    set(figs.error_sym, 'Color', 'w', 'Units', 'inches', 'Position', [1,1,5.5,4.5]);
    hold on
    
    ga(1) = semilogy(result.mvals, result.err_W, 's-', 'Color', color_W, ...
        'LineWidth', lg_linewidth, 'MarkerSize', lg_markersize);
    ga(2) = semilogy(result.mvals, result.est_K, '--', 'Color', color_K, ...
        'LineWidth', lg_linewidth_bnd);
    ga(3) = semilogy(result.mvals, result.est_W, '-.', 'Color', color_W, ...
        'LineWidth', lg_linewidth_bnd);
    ga(4) = semilogy(result.mvals, result.est_cheb, '-.', 'Color', color_MK, ...
        'LineWidth', lg_linewidth_bnd);
    
    % Error panel in logarithmic scale.
    set(gca, 'YScale', 'log')
    grid on
    box on
    xlabel('Krylov dimension $m$', 'Interpreter', 'latex');
    ylabel('Relative error and field-of-values estimate', 'Interpreter', 'latex');
    ax = gca;
    ax.LineWidth = axlabel_linewidth;
    ax.FontSize = axlabel_fontsize;
    % title(sprintf('%s, $n=%d$, $s=%d$', result.mat_label, result.n, result.s), 'interpreter', 'latex');
    set_m_axis(result.mvals)
    
    positive_vals = [result.err_W(:); result.est_K(:); result.est_cheb(:); result.est_W(:)];
    positive_vals = positive_vals(isfinite(positive_vals) & positive_vals > 0);
    if ~isempty(positive_vals)
        ylim([max(1e-16, min(positive_vals)/5), max(positive_vals)*5])
    end
    
    lgd = legend(ga, {'error', 'bound based on $\mathcal F(K)$', ...
        'bound based on $\mathcal F(W)$', 'bound based on Proposition 5.5'}, ...
        'Interpreter', 'latex', 'FontSize', lg_fontsize, 'Location', 'southwest');
    drawnow;
    increase_legend_width(lgd, 1.1);
end


end

function set_m_axis(mvals)
%SET_M_AXIS Set a stable x-axis range for one or more Krylov dimensions.

if min(mvals) < max(mvals)
    xlim([min(mvals), max(mvals)])
else
    xlim([mvals(1)-0.5, mvals(1)+0.5])
end

end
