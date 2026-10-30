% File: identifiability.m
% 读取 Abeta_CN.ipynb 保存的结果，比较 F_eta = J'J 与 K_eta = J'J + sum_n r_n H_n：
%   Figure 1: F_eta 与 K_eta 的特征值谱
%   Figure 2: 特征向量热图 (F_eta, K_eta 符号对齐) 及扰动 |V_K - V_F|
%   Figure 3: 累积 Abeta 预测区间 (K_eta vs F_eta)
%   Figure 4: 扰动耦合 log10 rho_ji = log10 |E_ji| / |lambda_i - lambda_j|
%   Figure 5: 最小不确定性时间 s_K, s_F 与区间宽度最小值
clc; clear; close all;

% === 1. 加载数据 ===
S = load('Abeta_CN_interior_recon_analysis_results', 'eigenvector_difference', ...
    'perturbation_check', 'prediction', 'prediction_comparison', 'minimum_uncertainty_time', ...
    'patient_age', 'patient_abeta');
ED = S.eigenvector_difference;
PT = S.perturbation_check;
PC = S.prediction_comparison;

lam_F = ED.eigenvalues_F(:);           % 降序
lam_K = ED.eigenvalues_K(:);           % 降序
V_F   = ED.eigenvectors_F;             % 列为特征向量
V_K   = ED.eigenvectors_K_aligned;     % 与 V_F 逐列符号对齐
dV    = ED.absolute_difference;        % |V_K - V_F|
rho   = PT.rho;                        % rho(j,i)，对角线为 NaN
n_par = numel(lam_F);
n_roi = n_par/2;                       % 前 68 个为 log(lambda)，后 68 个为 log(K)

set_K = PC.settings.K_eta;             % 预测设置: threshold, sigma, tau0, level
set_F = PC.settings.F_eta;
thr_K = set_K.threshold;               % 预测中 U/V 划分的特征值阈值
thr_F = set_F.threshold;
n_F = sum(lam_F > thr_F);
n_K = sum(lam_K > thr_K);
split_labels = { ...                   % 阈值划分的模态: F_eta 蓝色, K_eta 橙色
    sprintf('$F_\\eta$: $\\lambda_i > 10^{%d}$, $i \\le %d$', round(log10(thr_F)), n_F), ...
    sprintf('$K_\\eta$: $\\lambda_i > 10^{%d}$, $i \\le %d$', round(log10(thr_K)), n_K)};

% === 2. 作图设置 ===
fs      = 18;                          % Figure 1, 3, 4
fs_hm   = 14;                          % Figure 2 (热图)
lw      = 1.5;
lw_hm   = 2.5;                         % 热图上模态划分虚线
c_F     = hex2rgb('#2a78d6');          % F_eta: blue
c_K     = hex2rgb('#eb6834');          % K_eta: orange
c_ref   = [0.55 0.55 0.55];
cmap_div = make_cmap({'#104281','#5598e7','#f0efec','#e66767','#a11d1d'}, 256);
cmap_seq = make_cmap({'#f3f1fb','#9085e9','#4a3aa7','#231a5e'}, 256);
bf = @(s) ['\boldmath\textbf{' s '}'];  % LaTeX 粗体 (文字与公式)

%% Figure 1: 特征值谱
idx = (1:n_par)';
figure(1); clf;
set(gcf, 'Color', 'w', 'Position', [100 100 760 560]);
lam_F_plot = lam_F; lam_F_plot(lam_F <= 0) = NaN;     % F 的非正值为舍入误差
lam_K_plot = lam_K; lam_K_plot(lam_K <= 0) = NaN;
neg_K = lam_K < 0;

h1 = semilogy(idx, lam_F_plot, 'o-', 'Color', c_F, 'MarkerFaceColor', c_F, ...
    'MarkerSize', 5, 'LineWidth', lw); hold on
h2 = semilogy(idx, lam_K_plot, 's--', 'Color', c_K, 'MarkerFaceColor', c_K, ...
    'MarkerSize', 5, 'LineWidth', lw);
hs = [h1 h2];
labels = {'$\lambda_i(F_\eta)$', '$\lambda_i(K_\eta)$'};
if any(neg_K)
    h3 = semilogy(idx(neg_K), abs(lam_K(neg_K)), 's', 'Color', c_K, ...
        'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', lw);
    hs(end+1) = h3;
    labels{end+1} = '$|\lambda_i(K_\eta)|,\ \lambda_i<0$';
end
if thr_K == thr_F
    hs(end+1) = yline(thr_F, ':', 'Color', c_ref, 'LineWidth', lw);
    labels{end+1} = sprintf('$\\epsilon = 10^{%d}$', round(log10(thr_F)));
else
    hs(end+1) = yline(thr_F, ':', 'Color', c_F, 'LineWidth', lw);
    labels{end+1} = sprintf('$\\epsilon_F = 10^{%d}$', round(log10(thr_F)));
    hs(end+1) = yline(thr_K, ':', 'Color', c_K, 'LineWidth', lw);
    labels{end+1} = sprintf('$\\epsilon_K = 10^{%d}$', round(log10(thr_K)));
end

set(gca, 'FontName', 'Helvetica', 'FontSize', fs, 'FontWeight', 'bold', ...
    'LineWidth', lw, 'YMinorTick', 'on', 'TickDir', 'out')
xlim([0 n_par+1])
y_lo = min([lam_F(lam_F > 0); abs(lam_K(lam_K ~= 0))]);
ylim([10^floor(log10(y_lo)) 10])
xlabel('Index (descending)')
ylabel('Eigenvalue')
title(bf('Spectra of $F_\eta$ and $K_\eta$'), 'Interpreter', 'latex')
lgd = legend(hs, cellfun(bf, labels, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'Location', 'southwest');
lgd.Box = 'off';
box off

%% Figure 2: 特征向量热图及扰动
figure(2); clf;
set(gcf, 'Color', 'w', 'Position', [100 100 1500 480]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
mats   = {V_F, V_K, dV};
titles = {'Eigenvectors of $F_\eta$', ...
          'Eigenvectors of $K_\eta$ (sign-aligned)', ...
          'Perturbation $|V_{K}-V_{F}|$'};
for k = 1:3
    ax = nexttile(tl);
    imagesc(ax, mats{k});
    axis(ax, 'ij');
    if k < 3
        colormap(ax, cmap_div); ax.CLim = [-1 1];
    else
        colormap(ax, cmap_seq); ax.CLim = [0 max(dV(:))];
    end
    hold(ax, 'on');
    yline(ax, n_roi+0.5, '-', 'Color', 'k', 'LineWidth', 1.0);          % lambda / K 分块
    h_split = [mode_line(ax, @xline, n_F+0.5, c_F, lw_hm), ...         % lambda_i(F) > epsilon
               mode_line(ax, @xline, n_K+0.5, c_K, lw_hm)];            % lambda_i(K) > epsilon
    if k == 1, h_split_1 = h_split; end
    set(ax, 'YTick', [n_roi/2, 3*n_roi/2], ...
        'YTickLabel', {'$\log\lambda$', '$\log K$'}, ...
        'TickLabelInterpreter', 'latex', 'XTick', [1 20:20:n_par], ...
        'FontName', 'Helvetica', 'FontSize', fs_hm, 'LineWidth', 1.0, 'TickDir', 'out');
    if k > 1, ax.YTickLabel = {}; end
    xlabel(ax, 'Eigenvector index $i$', 'Interpreter', 'latex')
    title(ax, titles{k}, 'Interpreter', 'latex')
    if k >= 2
        cb = colorbar(ax);
        cb.FontSize = fs_hm - 2;
        cb.LineWidth = 1.0;
        if k == 2, cb.Ticks = -1:0.5:1; end
    end
end
ax1 = nexttile(tl, 1);
ylabel(ax1, 'Log parameter', 'Interpreter', 'latex')
lgd = legend(ax1, h_split_1, split_labels, 'Interpreter', 'latex', ...
    'Orientation', 'horizontal', 'FontSize', fs_hm);
lgd.Layout.Tile = 'south';
lgd.Box = 'off';

%% Figure 3: 预测区间对比
ages  = PC.ages(:);
bK    = PC.K;
bF    = PC.F;
center  = S.prediction.center(:);      % K_eta 与 F_eta 共用同一拟合中心
age_obs = S.patient_age(:);
tot_obs = sum(S.patient_abeta, 2);
level   = set_K.level;

figure(3); clf;
set(gcf, 'Color', 'w', 'Position', [100 100 1650 520]);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) 置信带
ax = nexttile(tl); hold(ax, 'on');
hF = fill_band(ax, ages, bF.lower(:), bF.upper(:), c_F, lw);
hK = fill_band(ax, ages, bK.lower(:), bK.upper(:), c_K, lw);
hm = plot(ax, ages, center, 'k-', 'LineWidth', lw);
hd = plot(ax, age_obs, tot_obs, 'o', 'MarkerSize', 9, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
xlim(ax, [ages(1) ages(end)])
ylabel(ax, 'Accumulated A\beta')
title(ax, sprintf('%.0f%% prediction bands', 100*level))
lgd = legend(ax, [hK hF hm hd], cellfun(bf, {'$K_\eta$', '$F_\eta$', 'Model', 'Data'}, ...
    'UniformOutput', false), 'Interpreter', 'latex', 'Location', 'northwest');
lgd.Box = 'off';

% % (b) 区间宽度
% ax = nexttile(tl); hold(ax, 'on');
% plot(ax, ages, bK.width(:), '-', 'Color', c_K, 'LineWidth', lw);
% plot(ax, ages, bF.width(:), '-', 'Color', c_F, 'LineWidth', lw);
% xlim(ax, [ages(1) ages(end)])
% ylabel(ax, 'Interval width')
% title(ax, 'Interval width (upper - lower)')
% lgd = legend(ax, cellfun(bf, {'$K_\eta$', '$F_\eta$'}, 'UniformOutput', false), ...
%     'Interpreter', 'latex', 'Location', 'northwest');
% lgd.Box = 'off';

% (c) 方差分解: 稳定子空间 sigma^2 U Lambda^{-1} U' 与正则化子空间 tau0^2 V V'
ax = nexttile(tl); hold(ax, 'on');
plot(ax, ages, bK.variance_stable(:),   '-',  'Color', c_K, 'LineWidth', lw);
plot(ax, ages, bK.variance_unstable(:), '--', 'Color', c_K, 'LineWidth', lw);
plot(ax, ages, bF.variance_stable(:),   '-',  'Color', c_F, 'LineWidth', lw);
plot(ax, ages, bF.variance_unstable(:), '--', 'Color', c_F, 'LineWidth', lw);
set(ax, 'YScale', 'log')
xlim(ax, [ages(1) ages(end)])
ylabel(ax, 'Variance of accumulated A\beta')
title(ax, 'Variance contributions')
lgd = legend(ax, cellfun(bf, { ...
    sprintf('$K_\\eta$ stable ($%d$)',      bK.stable_dimension), ...
    sprintf('$K_\\eta$ regularized ($%d$)', bK.regularized_dimension), ...
    sprintf('$F_\\eta$ stable ($%d$)',      bF.stable_dimension), ...
    sprintf('$F_\\eta$ regularized ($%d$)', bF.regularized_dimension)}, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'Location', 'southeast');
lgd.Box = 'off';

for ax = findobj(gcf, 'Type', 'axes')'
    set(ax, 'FontName', 'Helvetica', 'FontSize', fs, 'FontWeight', 'bold', ...
        'TitleFontWeight', 'bold', 'LineWidth', lw, 'TickDir', 'out')
    xlabel(ax, 'Age (years)')
    box(ax, 'off')
end

%% Figure 4: 扰动耦合 log10 rho_ji
log_rho = log10(rho);
log_rho(~isfinite(log_rho)) = NaN;
L_rho = 10;                            % 颜色范围 [-10, 10]，超出部分饱和

figure(4); clf;
set(gcf, 'Color', 'w', 'Position', [100 100 760 620]);
ax = axes;
imagesc(ax, log_rho, 'AlphaData', ~isnan(log_rho));
axis(ax, 'ij'); axis(ax, 'square');
colormap(ax, cmap_div); ax.CLim = [-L_rho L_rho];   % 以 rho = 1 为中心
ax.Color = 'w';
hold(ax, 'on');
h_split = [mode_line(ax, @xline, n_F+0.5, c_F, lw_hm), mode_line(ax, @xline, n_K+0.5, c_K, lw_hm)];
mode_line(ax, @yline, n_F+0.5, c_F, lw_hm);
mode_line(ax, @yline, n_K+0.5, c_K, lw_hm);
set(ax, 'FontName', 'Helvetica', 'FontSize', fs, 'FontWeight', 'bold', ...
    'LineWidth', lw, 'TickDir', 'out', 'XTick', [1 20:20:n_par], 'YTick', [1 20:20:n_par])
xlabel(ax, 'Perturbed eigenvector i')
ylabel(ax, 'Coupled eigenvector j')
title(ax, bf('$\log_{10}\rho_{ji}=\log_{10}\left(|E_{ji}|/|\lambda_i-\lambda_j|\right)$'), ...
    'Interpreter', 'latex')
cb = colorbar(ax);
cb.FontSize = fs - 2;
cb.FontWeight = 'bold';
cb.LineWidth = lw;
cb.Label.String = 'log_{10} \rho_{ji}';
cb.Label.FontSize = fs;
cb.Label.FontWeight = 'bold';
lgd = legend(ax, h_split, cellfun(bf, split_labels, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'Orientation', 'horizontal', 'Location', 'southoutside');
lgd.Box = 'off';

%% Figure 5: 最小不确定性时间 (近零子空间敏感性, 与 tau0 无关)
% s_M(t) = || V_M' * grad_eta y(t) ||_2,  y = sum_i x_i,  V_M 为 M = K_eta, F_eta 阈值以下的特征向量
MU = S.minimum_uncertainty_time;
mu_ages = MU.ages(:);
s_Kt    = MU.s_K(:);
s_Ft    = MU.s_F(:);
g_norm  = MU.total_gradient_norm(:);
g_norm(g_norm <= 0) = NaN;                     % t = 0 处梯度为零
w_Kt    = MU.width_K(:);
w_Ft    = MU.width_F(:);

figure(5); clf;
set(gcf, 'Color', 'w', 'Position', [100 100 1650 520]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) 近零敏感性 s_K, s_F
ax = nexttile(tl); hold(ax, 'on');
plot(ax, mu_ages, s_Kt, '-', 'Color', c_K, 'LineWidth', lw);
plot(ax, mu_ages, s_Ft, '-', 'Color', c_F, 'LineWidth', lw);
set(ax, 'YScale', 'log')
mark_minimum(ax, mu_ages, s_Kt, MU.argmin_age_s_K, c_K, lw, fs);
mark_minimum(ax, mu_ages, s_Ft, MU.argmin_age_s_F, c_F, lw, fs);
ylabel(ax, '||V^T\nabla_\eta y(t)||_2')
title(ax, 'Near-null sensitivity')
lgd = legend(ax, cellfun(bf, { ...
    sprintf('$s_K$ ($%d$ null modes)', MU.null_dimension_K), ...
    sprintf('$s_F$ ($%d$ null modes)', MU.null_dimension_F)}, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'Location', 'west');
lgd.Box = 'off';

% (b) 梯度落在近零子空间的比例
ax = nexttile(tl); hold(ax, 'on');
plot(ax, mu_ages, s_Kt./g_norm, '-', 'Color', c_K, 'LineWidth', lw);
plot(ax, mu_ages, s_Ft./g_norm, '-', 'Color', c_F, 'LineWidth', lw);
set(ax, 'YScale', 'log')
ylabel(ax, 'Near-null fraction')
title(ax, 'Relative near-null fraction')
lgd = legend(ax, cellfun(bf, {'$s_K/\|\nabla_\eta y\|_2$', '$s_F/\|\nabla_\eta y\|_2$'}, ...
    'UniformOutput', false), 'Interpreter', 'latex', 'Location', 'southeast');
lgd.Box = 'off';

% (c) 实际区间宽度及其最小值位置
ax = nexttile(tl); hold(ax, 'on');
plot(ax, mu_ages, w_Kt, '-', 'Color', c_K, 'LineWidth', lw);
plot(ax, mu_ages, w_Ft, '-', 'Color', c_F, 'LineWidth', lw);
mark_minimum(ax, mu_ages, w_Kt, MU.argmin_age_width_K, c_K, lw, fs);
mark_minimum(ax, mu_ages, w_Ft, MU.argmin_age_width_F, c_F, lw, fs);
ylabel(ax, 'Interval width')
title(ax, sprintf('%.0f%% interval width', 100*level))
lgd = legend(ax, cellfun(bf, {'$K_\eta$', '$F_\eta$'}, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'Location', 'northwest');
lgd.Box = 'off';

for ax = findobj(gcf, 'Type', 'axes')'
    set(ax, 'FontName', 'Helvetica', 'FontSize', fs, 'FontWeight', 'bold', ...
        'TitleFontWeight', 'bold', 'LineWidth', lw, 'TickDir', 'out')
    xlim(ax, [mu_ages(1) mu_ages(end)])
    xlabel(ax, 'Age (years)')
    box(ax, 'off')
end

%% 汇总
fprintf('Modes above threshold:  F_eta %d (> %.0e),  K_eta %d (> %.0e)  (of %d)\n', ...
    n_F, thr_F, n_K, thr_K, n_par);
fprintf('Negative K_eta eigenvalues (round-off level): %d, min %.3e\n', sum(neg_K), min(lam_K));
fprintf('max |V_K - V_F|: first %d modes %.3e,  remaining modes %.3e\n', ...
    n_F, max(dV(:, 1:n_F), [], 'all'), max(dV(:, n_F+1:end), [], 'all'));
fprintf('max rho_ji within the first %d modes: %.3e;  over all modes: %.3e\n', ...
    n_F, max(rho(1:n_F, 1:n_F), [], 'all', 'omitnan'), max(rho, [], 'all', 'omitnan'));
fprintf('Prediction settings K_eta: threshold = %g, sigma = %g, tau0 = %g, level = %g\n', ...
    set_K.threshold, set_K.sigma, set_K.tau0, set_K.level);
fprintf('Prediction settings F_eta: threshold = %g, sigma = %g, tau0 = %g, level = %g\n', ...
    set_F.threshold, set_F.sigma, set_F.tau0, set_F.level);
fprintf('Minimum-uncertainty age: s_K %.3f, s_F %.3f (shift %+.3f yr);  width_K %.3f, width_F %.3f (shift %+.3f yr)\n', ...
    MU.argmin_age_s_K, MU.argmin_age_s_F, MU.argmin_age_s_K - MU.argmin_age_s_F, ...
    MU.argmin_age_width_K, MU.argmin_age_width_F, MU.argmin_age_width_K - MU.argmin_age_width_F);
fprintf('Mean / max interval width:  K_eta %.4f / %.4f,  F_eta %.4f / %.4f\n', ...
    mean(bK.width), max(bK.width), mean(bF.width), max(bF.width));

%% 辅助函数
function rgb = hex2rgb(hex)
    hex = char(hex);
    if hex(1) == '#', hex = hex(2:end); end
    rgb = reshape(sscanf(hex, '%2x'), 1, 3)/255;
end

function cmap = make_cmap(hex_list, n)
    anchors = cell2mat(cellfun(@hex2rgb, hex_list(:), 'UniformOutput', false));
    x = linspace(0, 1, size(anchors, 1));
    cmap = interp1(x, anchors, linspace(0, 1, n));
end

function h = mode_line(ax, line_fun, pos, color, lw)
    % 白色衬底 + 彩色虚线，保证在深色热图上可见
    line_fun(ax, pos, '-', 'Color', 'w', 'LineWidth', lw + 1, 'HandleVisibility', 'off');
    h = line_fun(ax, pos, '--', 'Color', color, 'LineWidth', lw);
end

function mark_minimum(ax, x, y, x_min, color, lw, fs)
    % 最小值位置: 竖直点线 + 标记 + 年龄标注
    [~, i] = min(abs(x - x_min));
    xline(ax, x_min, ':', 'Color', color, 'LineWidth', lw, 'HandleVisibility', 'off');
    plot(ax, x(i), y(i), 'o', 'MarkerSize', 9, 'MarkerFaceColor', color, ...
        'MarkerEdgeColor', 'w', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    text(ax, x_min, y(i), sprintf('  %.2f', x_min), 'FontSize', fs - 4, 'FontWeight', 'bold', ...
        'Color', [0.2 0.2 0.2], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

function h = fill_band(ax, x, lo, hi, color, lw)
    h = fill(ax, [x; flipud(x)], [hi; flipud(lo)], color, ...
        'FaceAlpha', 0.22, 'EdgeColor', 'none');
    plot(ax, x, lo, '-', 'Color', color, 'LineWidth', lw);
    plot(ax, x, hi, '-', 'Color', color, 'LineWidth', lw);
end
