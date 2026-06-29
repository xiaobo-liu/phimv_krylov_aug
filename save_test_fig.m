function save_test_fig(fig, basename)
%SAVE_TEST_fig  Save one test figure as PDF and EPS.

fig_dir = fileparts(basename);
if ~isempty(fig_dir) && ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

drawnow

exportgraphics(fig, [basename, '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [basename, '.eps'], 'ContentType', 'vector');

end
