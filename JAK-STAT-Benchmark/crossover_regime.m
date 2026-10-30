%% crossover_regime.m
% Load and visualize the regime diagnostic (Part 9) exported from Pathway.ipynb:
% Lipschitz vs Hoelder-1/2 stability of the weak ln p2 direction at the fitted eta*.

clear; clc; close all;

%% Load data
dataDir  = fileparts(mfilename('fullpath'));
dataFile = fullfile(dataDir, 'crossover_analysis.mat');
S = load(dataFile);

Q        = S.quotient;          % diagnostic after removing the exact structural symmetry
V        = S.verify;            % numerical verification (second-order Taylor model)
epsCross = Q.eps_cross;         % gamma*nu/c2 (noise-sigma units)
sqrtN    = S.sqrt_n;            % norm of one full noise realisation
etaNames = cellstr(S.eta_names);

fprintf('Loaded %s\n', dataFile);
fprintf('nu = %.4f | gamma = %.4e | c2 = %.4e | M_N = %.4e | eps_cross = %.4f | sqrt(n) = %.3f\n', ...
        Q.nu, Q.gamma, Q.c2, Q.M_N, epsCross, sqrtN);
fprintf('weak direction (eta):');
for j = 1:numel(etaNames), fprintf(' %s:%+.3f', strtrim(etaNames{j}), S.weak_direction_eta(j)); end
fprintf('\n\n');

regimeTbl = table(S.eps_grid(:), string(cellstr(S.eps_labels(:))), Q.Xi(:), string(strtrim(cellstr(Q.regime(:)))), ...
                  'VariableNames', {'eps', 'label', 'Xi', 'regime'});
disp(regimeTbl);

%% Figure 1: worst-case displacement along the weak direction
figure('Name', 'Regime diagnostic', 'Color', 'w');
set(gcf, 'position', [203, 221, 700, 498]);
hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
set(gca, 'XScale', 'log', 'YScale', 'log');
e  = V.eps(:); d = V.sup_displacement(:);
[~, i10] = min(abs(e - 10));
plot(e, d, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
plot(e, d(1) * e / e(1), '--', 'Color', [0.12 0.47 0.71], 'LineWidth', 1.5);
plot(e, d(i10) * sqrt(e / e(i10)), '--', 'Color', [0.84 0.15 0.16], 'LineWidth', 1.5);
xline(epsCross, '--k', 'LineWidth', 1.2);
xline(sqrtN, '-.', 'Color', [0.17 0.63 0.17], 'LineWidth', 1.2);
xlabel('\epsilon = ||\delta||_2'); ylabel('||\eta_\delta - \eta^*||_2');
lgd = legend({'quadratic model (worst case)', 'Lipschitz', 'Hölder 1/2', '\epsilon_{cross}', '\surd{n} (noise)'}, ...
             'Location', 'northwest');
lgd.Box = 'off'; lgd.ItemTokenSize = [10, 6];

%% Figure 2: perturbation vs residual along the weak direction ln(p2*x1_0), one axes
gam = Q.gamma; nu = Q.nu; c2 = Q.c2;
e   = V.eps(:); d = V.sup_displacement(:);
[~, i10] = min(abs(e - 10));
ef  = logspace(-3, 2, 200);

figure('Name', 'Perturbation vs residual', 'Color', 'w');
set(gcf, 'position', [203, 221, 900, 560]);
xl = [1e-3 1e2]; yl = [1e-4 1e3];
% shaded regimes on a background axes (opaque light colours: transparency forces OpenGL, which mis-renders yyaxis log axes in batch mode)
ax0 = axes; hold(ax0, 'on');
patch(ax0, [xl(1) epsCross/3 epsCross/3 xl(1)],         [yl(1) yl(1) yl(2) yl(2)], [0.938 0.963 0.980], 'EdgeColor', 'none');
patch(ax0, [epsCross/3 3*epsCross 3*epsCross epsCross/3], [yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5],    'FaceAlpha', 0.12, 'EdgeColor', 'none');
patch(ax0, [3*epsCross xl(2) xl(2) 3*epsCross],         [yl(1) yl(1) yl(2) yl(2)], [0.989 0.941 0.941], 'EdgeColor', 'none');
set(ax0, 'XScale', 'log', 'YScale', 'log', 'XLim', xl, 'YLim', yl, 'Visible', 'off');
% data axes on top
ax = axes('Position', ax0.Position, 'Color', 'none'); hold(ax, 'on'); box off;
set(ax, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
yyaxis(ax, 'left'); set(ax, 'XScale', 'log', 'YScale', 'log', 'YColor', 'k'); xlim(xl); ylim(yl);
h1 = plot(e, d, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
h2 = plot(e, d(1) * e / e(1), 'k--', 'LineWidth', 1.5);
h3 = plot(e, d(i10) * sqrt(e / e(i10)), 'k:', 'LineWidth', 1.5);
ylabel('sup_\delta |\xi^T\Deltau| along ln(p_2x_1(0))');
yyaxis(ax, 'right'); set(ax, 'YScale', 'log', 'YColor', 'k'); ylim([1e-7 1e0]);
h4 = plot(ef, gam * nu * ones(size(ef)), '-', 'Color', [0.12 0.47 0.71], 'LineWidth', 2);
h5 = plot(ef, c2 * ef, '-', 'Color', [0.84 0.15 0.16], 'LineWidth', 2);
ylabel('coefficient of t^2 in l(t)');
xline(epsCross, '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xline(sqrtN, '-.', 'Color', [0.17 0.63 0.17], 'LineWidth', 1.2, 'HandleVisibility', 'off');
xlim(xl); xlabel('\epsilon = ||\delta||_2'); title('(i)');
lgd = legend([h1 h2 h3 h4 h5], {'worst-case displacement', 'slope 1', 'slope 1/2', '\gamma\nu (residual curvature)', 'c_2\epsilon (perturbation)'}, 'Location', 'northwest');
lgd.Box = 'off';
linkaxes([ax0 ax], 'x');

%% Figure 3: why ln(p2*x1_0) is the weak direction
Vt   = S.svd_Vt_J;                 % rows = right singular vectors of J_eta
sv   = S.singular_values_J(:);
z    = Vt(end, :);   if z(end) < 0, z = -z; end
w    = Vt(end-1, :); if w(2) < 0,   w = -w; end
wPred = S.weak_direction_predicted_eta(:)';
P    = load(fullfile(dataDir, 'profile_likelihood.mat'));
etaLabels = {'ln p_1', 'ln p_2', 'ln p_3', 'ln p_4', 'ln s_1', 'ln s_2', 'ln x_1(0)'};

figure('Name', 'Weak direction', 'Color', 'w');
set(gcf, 'position', [203, 221, 1500, 480]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (i) singular values
nexttile; hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
cols = repmat([0.5 0.5 0.5], 7, 1); cols(6, :) = [1 0.5 0.05]; cols(7, :) = [0.84 0.15 0.16];
b = bar(1:7, max(sv, 1e-16), 'FaceColor', 'flat'); b.CData = cols;
set(gca, 'YScale', 'log'); ylim([1e-15 1e3]); xlim([0.3 7.7]);
xlabel('index'); ylabel('singular value of J_\eta'); title('(i)');

% (ii) components of z and w
nexttile; hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
bb = bar(1:7, [z(:) w(:)], 'grouped');
bb(1).FaceColor = [0.84 0.15 0.16]; bb(2).FaceColor = [1 0.5 0.05];
plot((1:7) + 0.14, wPred, 'k_', 'MarkerSize', 14, 'LineWidth', 2);
yline(0, 'k', 'LineWidth', 0.6, 'HandleVisibility', 'off');
set(gca, 'XTick', 1:7, 'XTickLabel', etaLabels); xtickangle(30); ylabel('component'); title('(ii)');
lgd = legend({'z (\sigma_7): exact symmetry', 'w (\sigma_6): weak direction', 'e_{p_2} - (e_{p_2}\cdotz) z'}, 'Location', 'northeast'); lgd.Box = 'off';

% (iii) profiles of the six combination coordinates
nexttile; hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
combos = {'p2x1_0', 's1x1_0', 's2x1_0'}; cc = [0.84 0.15 0.16; 0.12 0.47 0.71; 0.17 0.63 0.17];
for k = 1:3
    C = P.PL_combo.(combos{k});
    plot(C.log10_grid(:) - C.log10_fit, C.delta_chi2(:), '.-', 'Color', cc(k, :), 'LineWidth', 1.5, 'MarkerSize', 10);
end
for k = 1:3
    nm = {'p1', 'p3', 'p4'}; C = P.PL.(nm{k});
    plot(C.log10_grid(:) - C.log10_fit, C.delta_chi2(:), '.-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1, 'MarkerSize', 8);
end
yline(P.threshold_pointwise, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
xlim([-1.5 3.5]); ylim([-0.5 12]);
xlabel('log_{10} shift from the optimum'); ylabel('\Delta\chi^2'); title('(iii)');
lgd = legend({'ln(p_2x_1(0))', 'ln(s_1x_1(0))', 'ln(s_2x_1(0))', 'ln p_1', 'ln p_3', 'ln p_4'}, 'Location', 'northeast'); lgd.Box = 'off';

%% Figure 4: eigenvalues and eigenvectors of F_eta
Id = load(fullfile(dataDir, 'identifiability_analysis.mat'), 'F_eta');
[UF, LF] = eig((Id.F_eta + Id.F_eta') / 2);
[lamF, idx] = sort(diag(LF), 'descend');    % lambda_i = sigma_i^2 of J_eta
UF = UF(:, idx);
for k = 1:size(UF, 2)                       % sign convention: largest component positive
    [~, imax] = max(abs(UF(:, k)));
    UF(:, k) = UF(:, k) * sign(UF(imax, k));
end
if UF(end, end) < 0, UF(:, end) = -UF(:, end); end   % z: x1_0 component positive, as in Figure 3
nEig = numel(lamF);

figure('Name', 'F_eta eigen-decomposition', 'Color', 'w');
set(gcf, 'position', [203, 221, 1144, 480]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (i) eigenvalues of F_eta
nexttile; hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
set(gca, 'YScale', 'log');
yl = [1e-14 1e6];
hL = plot(1:nEig, lamF, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'k');
hT = yline(1e-3 * lamF(1), ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
xlim([0.5 nEig + 0.5]); ylim(yl);
xticks(1:nEig); yticks(10.^(-14:4:6));
xlabel('index'); ylabel('eigenvalue'); title('(i)');
lgd = legend([hL hT], {'\lambda(F_\eta)', '10^{-3}\lambda_1'}, 'Location', 'southwest');
lgd.Box = 'off'; lgd.ItemTokenSize = [10, 6];

% (ii) eigenvector heatmap: rows = parameters, columns = eigenvectors
nexttile; hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
imagesc(1:nEig, 1:nEig, UF);
% diverging blue - white - red (RdBu-like): dark blue at -1, white at 0, dark red at +1
nc = 128; t = linspace(0, 1, nc)';
cBlue = [0.02 0.19 0.38]; cRed = [0.40 0.00 0.05]; cWhite = [0.97 0.97 0.97];
cmap = [cBlue + t .* (cWhite - cBlue); cWhite + t .* (cRed - cWhite)];
colormap(gca, cmap); caxis([-1 1]);
cb = colorbar; cb.LineWidth = 1.5;
for i = 1:nEig
    for j = 1:nEig
        if abs(UF(i, j)) > 0.6, tc = 'w'; else, tc = 'k'; end
        text(j, i, sprintf('%.2f', UF(i, j)), 'Color', tc, 'FontSize', 11, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
set(gca, 'YDir', 'reverse');
axis([0.5 nEig+0.5 0.5 nEig+0.5]);
xticks(1:nEig); yticks(1:nEig); yticklabels(etaLabels);
xlabel('eigenvector index (desc. eigenvalue)'); title('(ii)');

%% Figure 5: v_7 is the exact symmetry direction, v_6 is the weak direction log10(p2*x1_0)
v7 = UF(:, nEig);                                  % sign: x1_0 component positive
v6 = UF(:, nEig-1); if v6(2) < 0, v6 = -v6; end     % sign: p2 component positive
% Map eta = ln(theta) to the identifiable combinations q = (p1, p2*x1_0, p3, p4, s1*x1_0, s2*x1_0): ln q = A*eta
A = [1 0 0 0 0 0 0;
     0 1 0 0 0 0 1;
     0 0 1 0 0 0 0;
     0 0 0 1 0 0 0;
     0 0 0 0 1 0 1;
     0 0 0 0 0 1 1];
dq7 = A * v7 / log(10);                            % change of log10 q per unit step along v_7 (= 0: invisible)
dq6 = A * v6 / log(10);                            % change of log10 q per unit step along v_6 (only p2*x1_0)
comboLabels = {'p_1', 'p_2x_1(0)', 'p_3', 'p_4', 's_1x_1(0)', 's_2x_1(0)'};
zCol = [0.84 0.15 0.16]; wCol = [1 0.5 0.05];

figure('Name', 'Symmetry and weak direction', 'Color', 'w');
set(gcf, 'position', [203, 221, 700, 480]);
hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
bq = bar(1:numel(dq6), [dq7 dq6], 'grouped');
bq(1).FaceColor = zCol; bq(2).FaceColor = wCol;
yline(0, 'k', 'LineWidth', 0.6, 'HandleVisibility', 'off');
xlim([0.5 numel(dq6) + 0.5]); ylim([-0.1 0.6]);
set(gca, 'XTick', 1:numel(dq6), 'XTickLabel', comboLabels); xtickangle(30);
ylabel('\Delta log_{10} q per unit step');
lgd = legend({'along v_7', 'along v_6'}, 'Location', 'northeast');
lgd.Box = 'off'; lgd.ItemTokenSize = [10, 6];

%% Figure 6: F_eta eigenvalue, c_2 and gamma*nu of the weak (v_6) and exact-symmetry (v_7) directions
% c_2 = min ||Pi_perp H[xi,xi]||: second-order output change that the other parameters cannot compensate
% gamma*nu = lambda_min(V' C V): curvature the residual adds along the direction
%   v_7: S.tau_1e6  (null space of J_eta at tau = 1e-6 = exact symmetry direction only)
%   v_6: S.quotient (null direction left after removing the exact symmetry = weak direction)
barVals = [lamF(nEig-1), S.quotient.c2, S.quotient.gamma * S.quotient.nu;   % v_6
           lamF(nEig),   S.tau_1e6.c2,  S.tau_1e6.gamma  * S.tau_1e6.nu];   % v_7

figure('Name', 'Eigenvalue, c2 and gamma*nu of v6, v7', 'Color', 'w');
set(gcf, 'position', [203, 221, 700, 480]);
hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
set(gca, 'YScale', 'log');
hb = bar(1:2, barVals, 'grouped', 'BaseValue', 1e-16);
hb(1).FaceColor = [0 0 0];
hb(2).FaceColor = [0.84 0.15 0.16];
hb(3).FaceColor = [0.12 0.47 0.71];
for k = 1:3                                        % value labels on top of each bar
    xt = hb(k).XEndPoints;
    text(xt, barVals(:, k)' * 3, arrayfun(@(v) sprintf('%.1e', v), barVals(:, k)', 'UniformOutput', false), 'FontSize', 14, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end
xlim([0.4 2.6]); ylim([1e-16 1e2]);
yticks(10.^(-16:4:0));
set(gca, 'XTick', 1:2, 'XTickLabel', {'v_6 (weak)', 'v_7 (symmetry)'});
ylabel('value');
lgd = legend(hb, {'\lambda(F_\eta)', 'c_2', '\gamma\nu'}, 'Location', 'northeast');
lgd.Box = 'off'; lgd.ItemTokenSize = [10, 6];
