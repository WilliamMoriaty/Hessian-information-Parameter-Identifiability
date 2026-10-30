% identifiability.m
% Identifiability figures for patient 3 (results saved by HIV.ipynb):
%   Figure 1  eigenvalues and eigenvectors of F_eta and K_eta
%   Figure 2  eigenvector perturbation heat map, rho_ji = |E_ji|/|lambda_i-lambda_j|
%   Figure 3  residual-Hessian indices I_p^RH, I_{p,resp}^RH, I_{p,curv}^RH
%   Figure 4  CrI_{eta_p} (sigma^2/kappa_p^K) versus MCMC credible intervals

clear; clc;

S = load('HIV_patient3_results.mat');
names   = S.FREE_PARAM_TXT(:);
n       = numel(names);
idx     = 1:n;
log_tex = cellfun(@(s) ['log ' param_tex(s)], names, 'UniformOutput', false);
font    = {'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.2};

%% Figure 1: eigenvalues and eigenvectors of F_eta and K_eta
% F_eta eigenpairs come from the SVD of J_eta, K_eta from eigh (both descending).
lam_F = S.eigvals_logtheta(:);
U_F   = align_signs(S.eigvecs_logtheta, []);
lam_K = S.K_eta_results.eigenvalues_descending(:);
U_K   = align_signs(S.K_eta_results.eigenvectors_descending, U_F);

fig1 = figure(1);
clf();
set(gcf, "Position", [100, 100, 1350, 400])
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
h = gobjects(0); lab = {};
h(end+1) = semilogy(idx, abs(lam_F), '-o', 'Color', 'k', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 7, 'LineWidth', 1.5); lab{end+1} = '\lambda(F_\eta)';
hold on
h(end+1) = semilogy(idx, abs(lam_K), '--s', 'Color', 'r', 'MarkerSize', 10, ...
    'LineWidth', 1.5); lab{end+1} = '\lambda(K_\eta)';
if any(lam_K < 0)
    h(end+1) = semilogy(idx(lam_K < 0), abs(lam_K(lam_K < 0)), 's', 'Color', 'r', ...
        'MarkerFaceColor', 'w', 'MarkerSize', 13, 'LineWidth', 1.5);
    lab{end+1} = '\lambda(K_\eta) < 0 (abs. shown)';
end
hold off
xlim([0.5, n+0.5])
set(gca, 'XTick', idx)
xlabel('index')
ylabel('eigenvalue')
title('Spectra')
set(gca, font{:})
lgd = legend(h, lab, 'Location', 'southwest', 'Box', 'off', 'FontSize', 11);
box off

nexttile
eig_heatmap(U_F, log_tex, 'Eigenvectors of F_\eta', font)
nexttile
eig_heatmap(U_K, log_tex, 'Eigenvectors of K_\eta', font)
cb = colorbar;
cb.Layout.Tile = 'east';

%% Figure 2: eigenvector perturbation heat map
rho  = S.perturbation_check.rho_scaled_coupling;   % NaN on the diagonal
Lrho = log10(rho);

fig2 = figure(2);
clf();
set(gcf, "Position", [150, 150, 560, 470])
imagesc(Lrho, 'AlphaData', ~isnan(Lrho));
set(gca, 'Color', [0.85, 0.85, 0.85])              % diagonal shown in gray
colormap(gca, redblue(256));
cb = colorbar;
cb.Label.String = 'log_{10}\rho_{ji}';
% White text on the dark ends of the map, as in Figure 1.
Lrange = [min(Lrho(:)), max(Lrho(:))];
for j = 1:n
    for i = 1:n
        if i ~= j
            s = 2*(Lrho(j,i) - Lrange(1))/diff(Lrange) - 1;
            if abs(s) > 0.6, tc = 'w'; else, tc = 'k'; end
            text(i, j, sprintf('%.1e', rho(j,i)), 'HorizontalAlignment', 'center', ...
                'Color', tc, 'FontSize', 13, 'FontWeight', 'bold');
        end
    end
end
axis square
set(gca, 'XTick', idx, 'YTick', idx, 'TickLength', [0, 0])
xlabel('mode i')
ylabel('mode j')
title('\rho_{ji} = |E_{ji}| / |\lambda_i(F_\eta) - \lambda_j(F_\eta)|')
set(gca, font{:})
cb.FontSize = 14;
cb.Label.FontSize = 18;
cb.Label.FontWeight = 'bold';

%% Figure 3: residual-Hessian correction indices
RH = S.residual_hessian_coordinate_stability;
I  = [RH.I_RH(:), RH.I_RH_response(:), RH.I_RH_curvature(:)];

fig3 = figure(3);
clf();
set(gcf, "Position", [200, 200, 760, 430])
b = bar(idx, I, 'grouped', 'EdgeColor', 'none');
b(1).FaceColor = [0.49, 0.18, 0.56];
b(2).FaceColor = [0.00, 0.45, 0.74];
b(3).FaceColor = [0.93, 0.69, 0.13];
for k = 1:3
    y  = b(k).YEndPoints;
    va = repmat({'bottom'}, size(y));
    va(y < 0) = {'top'};
    for m = 1:numel(y)
        text(b(k).XEndPoints(m), y(m), sprintf('%.2f', y(m)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', va{m}, ...
            'FontSize', 10, 'FontWeight', 'bold');
    end
end
yline(0, 'k-', 'LineWidth', 1.0);
ylim([min(I(:)) - 1, max(I(:)) + 1])
set(gca, 'XTick', idx, 'XTickLabel', log_tex)
ylabel('log correction')
set(gca, font{:})
legend(b, {'I^{RH}_p', 'I^{RH}_{p,resp}', 'I^{RH}_{p,curv}'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
box off

%% Figure 4: CrI_{eta_p} versus MCMC
% Quantile rows: [alpha/2, 0.25, 0.5, 0.75, 1-alpha/2]; columns: parameters.
% Coordinate: Delta eta_p = log(theta_p / theta_p^*).
C     = S.credible_interval_comparison;
Qm    = C.mcmc_quantiles;
Qk    = C.paper_quantiles;
ok_k  = logical(C.paper_interval_computable(:));
c_cri = [0.49, 0.18, 0.56];
c_mc  = [0.20, 0.63, 0.17];
off   = 0.17;
w     = 0.26;

fig4 = figure(4);
clf();
set(gcf, "Position", [250, 250, 740, 470])
hold on
for p = 1:n
    if ok_k(p)
        h_cri = draw_box(Qk(:,p), p - off, w, c_cri);
    end
    h_mc = draw_box(Qm(:,p), p + off, w, c_mc);
end
h_star = plot(zeros(n,1), idx, 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 12);
xline(0, 'k--', 'LineWidth', 1.0, 'Alpha', 0.5);
hold off
set(gca, 'YDir', 'reverse', 'YTick', idx, 'YTickLabel', log_tex)
ylim([0.4, n + 0.6])
xlabel('\Delta\eta_p = log(\theta_p / \theta_p^*)')
title(sprintf('%.0f%% credible intervals', 100*C.level))
set(gca, font{:})
legend([h_cri, h_mc, h_star], {'CrI_{\eta_p}: \sigma^2/\kappa_p^K', 'MCMC', '\eta^*'}, ...
    'Location', 'northeast', 'Box', 'off');
box off

fprintf('%.0f%% interval width, CrI_eta_p / MCMC:\n', 100*C.level);
for p = 1:n
    fprintf('  %-7s %.4f\n', names{p}, C.width_ratio_paper_over_mcmc(p));
end

% exportgraphics(fig1, 'HIV_3_eigen_F_K.png', 'Resolution', 300);
% exportgraphics(fig2, 'HIV_3_perturbation.png', 'Resolution', 300);
% exportgraphics(fig3, 'HIV_3_I_RH.png', 'Resolution', 300);
% exportgraphics(fig4, 'HIV_3_CrI_MCMC.png', 'Resolution', 300);

%% Local functions
function s = param_tex(name)
    switch name
        case 'lambda', s = '\lambda';
        case 'delta',  s = '\delta';
        case 'p',      s = '\pi';
        otherwise,     s = name;
    end
end

function V = align_signs(V, ref)
    % Fix the arbitrary sign of each eigenvector: largest-magnitude entry
    % positive, or non-negative inner product with the reference column.
    for c = 1:size(V, 2)
        if isempty(ref)
            [~, k] = max(abs(V(:,c)));
            s = sign(V(k,c));
        else
            s = sign(ref(:,c)' * V(:,c));
        end
        if s < 0
            V(:,c) = -V(:,c);
        end
    end
end

function eig_heatmap(V, row_labels, ttl, font)
    imagesc(V);
    clim([-1, 1]);
    colormap(gca, redblue(256));
    [nr, nc] = size(V);
    for r = 1:nr
        for c = 1:nc
            if abs(V(r,c)) > 0.6, tc = 'w'; else, tc = 'k'; end
            text(c, r, sprintf('%.2f', V(r,c)), 'HorizontalAlignment', 'center', ...
                'Color', tc, 'FontSize', 11, 'FontWeight', 'bold');
        end
    end
    set(gca, 'XTick', 1:nc, 'YTick', 1:nr, 'YTickLabel', row_labels, 'TickLength', [0, 0])
    xlabel('eigenvector index')
    title(ttl)
    set(gca, font{:})
end

function cmap = redblue(m)
    % Diverging blue-white-red map (ColorBrewer RdBu, reversed).
    anchors = [  5  48  97;  33 102 172;  67 147 195; 146 197 222; 209 229 240; 247 247 247; ...
               253 219 199; 244 165 130; 214  96  77; 178  24  43; 103   0  31] / 255;
    cmap = interp1(linspace(0, 1, size(anchors, 1)), anchors, linspace(0, 1, m));
end

function hb = draw_box(q, y, w, c)
    % Horizontal box: whiskers q(1)-q(5), box q(2)-q(4), median q(3).
    plot([q(1), q(2)], [y, y], '-', 'Color', c, 'LineWidth', 1.5);
    plot([q(4), q(5)], [y, y], '-', 'Color', c, 'LineWidth', 1.5);
    plot([q(1), q(1)], y + [-w, w]/4, '-', 'Color', c, 'LineWidth', 1.5);
    plot([q(5), q(5)], y + [-w, w]/4, '-', 'Color', c, 'LineWidth', 1.5);
    hb = patch([q(2), q(4), q(4), q(2)], y + [-w, -w, w, w]/2, c, ...
        'FaceAlpha', 0.3, 'EdgeColor', c, 'LineWidth', 1.5);
    plot([q(3), q(3)], y + [-w, w]/2, '-', 'Color', c, 'LineWidth', 2.5);
end
