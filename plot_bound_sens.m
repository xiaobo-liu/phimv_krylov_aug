function fig = plot_bound_sens(results, sens_name)
%PLOT_BOUND_SENS Plot convergence bounds for one right-hand-side sensitivity test.

lg_linewidth = 1.5;
lg_linewidth_bnd = 1.2;
lg_markersize = 4;
lg_fontsize = 9;

axlabel_linewidth = 1.0;
axlabel_fontsize = 9;

color_K = [0.635 0.078 0.184];
color_MK = [0 0.4470 0.7410];

num_cases = numel(results);
num_cols = ceil(sqrt(num_cases));
num_rows = ceil(num_cases / num_cols);

fig = figure;
clf(fig)
set(fig, 'Color', 'w', 'Units', 'inches', ...
    'Position', [1, 1, 3.8*num_cols, 3.3*num_rows]);

tlo = tiledlayout(fig, num_rows, num_cols, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:num_cases
    result = results{k};

    nexttile(tlo)
    hold on
    ga(1) = semilogy(result.mvals, result.err_kiops, 'o-', ...
        'Color', color_K, 'LineWidth', lg_linewidth, ...
        'MarkerSize', lg_markersize);
    ga(2) = semilogy(result.mvals, result.est_K, '--', ...
        'Color', color_K, 'LineWidth', lg_linewidth_bnd);
    ga(3) = semilogy(result.mvals, result.est_MK, '-.', ...
        'Color', color_MK, 'LineWidth', lg_linewidth_bnd);

    set(gca, 'YScale', 'log')
    grid on
    box on
    xlabel('Krylov dimension $m$', 'Interpreter', 'latex');
    ylabel('Relative error and field-of-values estimate', 'Interpreter', 'latex');
    title(case_title(result, sens_name), 'Interpreter', 'latex');
    ax = gca;
    ax.LineWidth = axlabel_linewidth;
    ax.FontSize = axlabel_fontsize;
    set_m_axis(result.mvals)

    positive_vals = [result.err_kiops(:); result.est_K(:); result.est_MK(:)];
    positive_vals = positive_vals(isfinite(positive_vals) & positive_vals > 0);
    if ~isempty(positive_vals)
        % Keep all positive curves visible on the logarithmic scale.
        ylim([max(1e-16, min(positive_vals)/5), max(positive_vals)*5])
    end
    if k == 1 % only one (same) legend needed in the plots
        lgd = legend(ga, {'error', 'bound based on $\mathcal F(K)$', ...
            'bound based on $\mathcal F_{M_{\rm K}}(K)$'}, ...
            'Interpreter', 'latex', ...
            'FontSize', lg_fontsize, 'Location', 'southwest');
        drawnow;
        increase_legend_width(lgd, 1.1);
    end
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
