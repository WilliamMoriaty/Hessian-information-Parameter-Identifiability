%% identifiability.m
% Load the local identifiability analysis exported from Pathway.ipynb
% (JAK-STAT model, eta = ln(theta) coordinates unless noted otherwise).

clear; clc; close all;

%% Load data
dataDir = fileparts(mfilename('fullpath'));
S = load(fullfile(dataDir, 'identifiability_analysis.mat'));

paramNames = cellstr(S.param_names);        % {'p1','p2','p3','p4','s1','s2','x1_0'}
etaNames   = cellstr(S.eta_names);          % {'ln p1', ..., 'ln x1_0'}
thetaFit   = S.theta_fit;                   % best fit, linear scale
etaFit     = S.eta_fit;                     % best fit, eta = ln(theta)
nParams    = numel(paramNames);

% Residuals and sensitivities at the best fit (31 data rows)
r_star   = S.r_star;                        % weighted residuals, sum(r.^2) = chi2_min
rowTime  = S.row_time;                      % time of each data row
rowObs   = S.row_observable;                % 1 = y1 (pSTAT), 2 = y2 (tSTAT)
J_eta    = S.J_eta;                         % 31 x 7 Jacobian dr/d(eta)
H_eta    = S.H_eta;                         % 31 x 7 x 7 second derivatives d2r/d(eta)2
J_theta  = S.J_theta;                       % same, in theta coordinates
H_theta  = S.H_theta;

% 7 x 7 matrices: F = J'*J (Fisher / Gauss-Newton), C = sum_i r_i*H_i (curvature), K = F + C (full Hessian)
F_eta   = S.F_eta;    C_eta   = S.C_eta;    K_eta   = S.K_eta;
F_theta = S.F_theta;  C_theta = S.C_theta;  K_theta = S.K_theta;

% Analysis results
T2    = S.theorem2;                         % Theorem 2 quantities (rank, null space, perturbation bounds)
EPeta = S.eig_perturbation_eta;             % eigen-perturbation F -> K, eta coordinates
EPth  = S.eig_perturbation_theta;           % eigen-perturbation F -> K, theta coordinates
CM    = S.coordinate_metrics;               % per-parameter metrics: CM.<param>.{kappa_F, L_K, ...}
Yao   = S.yao;                              % Yao et al. orthogonalization ranking
Qua   = S.quaiser;                          % Quaiser & Moenkemoeller eigenvalue method
Cmp   = S.comparison;                       % rankings of the three methods

% Same analyses in the 6 identifiable combinations q = (p1, p2*x1_0, p3, p4, s1*x1_0, s2*x1_0)
comboNames = cellstr(S.combo_names);        % {'p1','p2x1_0','p3','p4','s1x1_0','s2x1_0'}
log10Combo = S.log10_combo_fit;             % best fit of the combinations, log10 scale
nCombo     = numel(comboNames);
YaoC  = S.yao_combo;                        % Yao ranking on the combinations
QuaC  = S.quaiser_combo;                    % Quaiser method on the combinations
CmpC  = S.comparison_combo;                 % rankings of the three methods on the combinations
CMc   = S.coordinate_metrics_combo;         % per-combination metrics: CMc.<combo>.{L_K, L_F, ...}

fprintf('Loaded %s\n', fullfile(dataDir, 'identifiability_analysis.mat'));
fprintf('chi2_min = sum(r_star.^2) = %.4f\n\n', sum(r_star.^2));

%% Eigenvalues of F and K
eigTbl = table(EPeta.lamF(:), EPeta.lamK(:), EPth.lamF(:), EPth.lamK(:), ...
               'VariableNames', {'lamF_eta', 'lamK_eta', 'lamF_theta', 'lamK_theta'});
disp('Eigenvalues of F and K (descending):');
disp(eigTbl);

%% Theorem 2: rank and null space
fprintf('Theorem 2: rank(F) = %d, dim ker(K) = %d, sin(Theta_U) = %.3e, tan(theta_max) = %.3e\n', ...
        T2.r, T2.ker_dim_K, T2.sin_Theta_U, T2.tan_theta_max);
nullTbl = table(string(etaNames(:)), T2.V(:), 'VariableNames', {'eta', 'V_null'});
disp('Null-space direction V of F:');
disp(nullTbl);

%% Per-parameter coordinate metrics
M = cellfun(@(n) CM.(n), paramNames, 'UniformOutput', false);
M = [M{:}];                                 % 1 x 7 struct array
metricsTbl = struct2table(M(:));
metricsTbl.name = string(metricsTbl.name);
disp('Coordinate metrics:');
disp(metricsTbl(:, {'name', 'kappa_F', 'kappa_K', 'L_F', 'L_K', 'R_RH', ...
                    'ident_F', 'ident_K', 'sd_ln', 'sd_log10'}));

%% Ranking methods
fprintf('Yao (eps = %.3g): estimable = {%s}, excluded = {%s}\n', Yao.eps, ...
        strjoin(cellstr(Yao.estimable), ', '), strjoin(cellstr(Yao.excluded), ', '));
fprintf('Quaiser (eps = %.3g): identifiable = {%s}, unidentifiable = {%s}\n\n', Qua.eps, ...
        strjoin(cellstr(Qua.identifiable), ', '), strjoin(cellstr(Qua.unidentifiable), ', '));

rankTbl = table((1:nParams)', string(cellstr(Cmp.rank_yao_least_first))', ...
                string(cellstr(Cmp.rank_quaiser_least_first))', ...
                string(cellstr(Cmp.rank_haowang_LK_least_first))', ...
                Cmp.lam_min_path_yao(:), Cmp.lam_min_path_quaiser(:), Cmp.lam_min_path_haowang(:), ...
                'VariableNames', {'Step', 'Yao', 'Quaiser', 'HaoWang_LK', ...
                                  'lam_min_Yao', 'lam_min_Quaiser', 'lam_min_HaoWang'});
disp('Rankings (least -> most identifiable) and lambda_min along each removal path:');
disp(rankTbl);

%% Ranking methods on the identifiable combinations
fprintf('\n[Combinations] Yao (eps = %.3g): estimable = {%s}, excluded = {%s}\n', YaoC.eps, ...
        joinNames(YaoC.estimable), joinNames(YaoC.excluded));
fprintf('[Combinations] Quaiser (eps = %.3g): identifiable = {%s}, unidentifiable = {%s}\n', QuaC.eps, ...
        joinNames(QuaC.identifiable), joinNames(QuaC.unidentifiable));

stepsC = QuaC.steps;
if iscell(stepsC), stepsC = [stepsC{:}]; end
for s = 1:numel(stepsC)
    % u1 field names such as 'ln(p2*x1_0)' are not valid MATLAB identifiers, so read them via fieldnames
    u1Names = fieldnames(stepsC(s).u1);
    u1Vals  = cell2mat(struct2cell(stepsC(s).u1));
    u1Str   = strjoin(compose('%s:%+.3f', string(u1Names), u1Vals), ', ');
    if isempty(stepsC(s).removed)
        action = 'stop';
    else
        action = ['remove ' stepsC(s).removed];
    end
    fprintf('  n_I = %d: lambda_min = %.4g, u1 = [%s] -> %s\n', stepsC(s).n_I, stepsC(s).lam_min, u1Str, action);
end
fprintf('[Combinations] parameters to fix until lambda_min >= eps: Yao %d, Quaiser %d, HaoWang %d\n\n', ...
        CmpC.n_fix_yao, CmpC.n_fix_quaiser, CmpC.n_fix_haowang);

rankTblC = table((1:nCombo)', string(cellstr(CmpC.rank_yao_least_first))', ...
                 string(cellstr(CmpC.rank_quaiser_least_first))', ...
                 string(cellstr(CmpC.rank_haowang_LK_least_first))', ...
                 CmpC.lam_min_path_yao(:), CmpC.lam_min_path_quaiser(:), CmpC.lam_min_path_haowang(:), ...
                 'VariableNames', {'Step', 'Yao', 'Quaiser', 'HaoWang_LK', ...
                                   'lam_min_Yao', 'lam_min_Quaiser', 'lam_min_HaoWang'});
disp('[Combinations] Rankings (least -> most identifiable) and lambda_min along each removal path:');
disp(rankTblC);

%% Plot settings
dispNames = paramNames;
dispNames(strcmp(paramNames, 'x1_0')) = {'x_1(0)'};
x        = 1:nParams;
blueCol  = [0 0.4470 0.7410];
redCol   = [0.85 0.1 0.1];
orangeCol = [0.8500 0.3250 0.0980];

%% Figure 1: coordinate Lipschitz coefficients L_p^K and L_p^F
LK = [M.L_K];                               % residual-corrected, ||g_p^K|| / kappa_p^K
LF = [M.L_F];                               % Gauss-Newton, 1 / sqrt(kappa_p^F)
LKplot = LK;  LKplot(~isfinite(LK)) = NaN;  % Inf = non-identifiable, marked with an infinity sign
LFplot = LF;  LFplot(~isfinite(LF)) = NaN;
yTop = 1.5*max([LKplot LFplot]);            % same range on both axes so the bars are comparable

w = 0.35;                                   % bar width; left/right bars are offset by -/+ w/2

figure('Name', 'Coordinate metrics', 'Color', 'w');
set(gcf, 'position', [203, 221, 800, 420]);
hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);

% Left axis: L_p^K
yyaxis left; hold on;
hLK = bar(x - w/2, LKplot, w, 'FaceColor', blueCol);
ylim([0 yTop]);
text(x(~isfinite(LK)) - w/2, 0.85*yTop*ones(1, nnz(~isfinite(LK))), '\infty', 'Color', blueCol, ...
     'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
ylabel('L_p^K');

% Right axis: L_p^F
yyaxis right; hold on;
hLF = bar(x + w/2, LFplot, w, 'FaceColor', orangeCol);
ylim([0 yTop]);
text(x(~isfinite(LF)) + w/2, 0.85*yTop*ones(1, nnz(~isfinite(LF))), '\infty', 'Color', orangeCol, ...
     'FontSize', 24, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
ylabel('L_p^F');

ax = gca;
ax.YAxis(1).Color = blueCol;
ax.YAxis(2).Color = orangeCol;
xlim([0.5 nParams + 0.5]);
xticks(x); xticklabels(dispNames);
xlabel('Parameter');
lgd = legend([hLK hLF], {'L_p^K', 'L_p^F'}, 'Location', 'northoutside', 'Orientation', 'horizontal');
lgd.Box = 'off';
lgd.ItemTokenSize = [10, 6];

%% Figure 2: Yao orthogonal method
order = double(Yao.order_1based);           % parameters in selection order
Rll   = Yao.Rll;                            % normalized orthogonal residual |R_LL| at each step
isEst = x <= double(Yao.n_estimable);

figure('Name', 'Yao method', 'Color', 'w');
set(gcf, 'position', [290, 357, 830, 420]);
hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
set(gca, 'YScale', 'log');
b1 = bar(x(isEst),  Rll(isEst),  0.6, 'FaceColor', blueCol, 'BaseValue', 1e-16);
b2 = bar(x(~isEst), Rll(~isEst), 0.6, 'FaceColor', redCol,  'BaseValue', 1e-16);
hEps = yline(Yao.eps, '--k', 'LineWidth', 1.5);
ylim([1e-16 10]);
yticks(10.^(-16:4:0));
xticks(x); xticklabels(dispNames(order));
xlabel('Parameter (Yao ranking)');
ylabel('|R_{LL}|');
lgd = legend([b1 b2 hEps], {'Estimable', 'Excluded', sprintf('\\epsilon = %.2g', Yao.eps)}, ...
             'Location', 'eastoutside');
lgd.Box = 'off';
lgd.ItemTokenSize = [10, 6];

%% Figure 3: Quaiser & Moenkemoeller eigenvalue method
steps = Qua.steps;
if iscell(steps), steps = [steps{:}]; end
nSteps = numel(steps);
lamMin = [steps.lam_min];                   % smallest eigenvalue of H at each iteration
U1     = nan(nParams, nSteps);              % |u_1| component per parameter (NaN = already removed)
stepLbl = strings(1, nSteps);
for s = 1:nSteps
    for k = 1:nParams
        fn = ['ln_' paramNames{k}];
        if isfield(steps(s).u1, fn)
            U1(k, s) = abs(steps(s).u1.(fn));
        end
    end
    if isempty(steps(s).removed)
        stepLbl(s) = "stop";
    else
        stepLbl(s) = dispNames{strcmp(etaNames, steps(s).removed)} + " removed";
    end
end

figure('Name', 'Quaiser method', 'Color', 'w');
set(gcf, 'position', [203, 221, 1000, 400]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (i) |lambda_min| per iteration (step 1 is -1.8e-12, i.e. numerically zero, so plot the magnitude)
nexttile; hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
set(gca, 'YScale', 'log');
hLam = plot(1:nSteps, abs(lamMin), 'o-', 'Color', blueCol, 'LineWidth', 1.5, ...
            'MarkerSize', 8, 'MarkerFaceColor', blueCol);
epsQ = 0.01;                                % plotted threshold (the saved analysis used 1e-3*lambda_max = 15.7; same removals)
hEps = yline(epsQ, '--r', 'LineWidth', 1.5);
text(1:nSteps, abs(lamMin)/8, stepLbl, 'FontSize', 16, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
xlim([0.4 nSteps + 0.6]);
ylim([1e-16 1e4]);
yticks(10.^(-16:4:4));
xticks(1:nSteps);
xlabel('Iteration');
ylabel('|\lambda_{min}|');
title('(i)');
lgd = legend([hLam hEps], {'|\lambda_{min}|', sprintf('\\epsilon = %g', epsQ)}, 'Location', 'southeast');
lgd.Box = 'off';
lgd.ItemTokenSize = [10, 6];

% (ii) components of the smallest eigenvector u_1 at each iteration
nexttile; hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
bar(x, U1, 'grouped');
ylim([0 1.2]);
xticks(x); xticklabels(dispNames);
xlabel('Parameter');
ylabel('|u_{1,k}|');
title('(ii)');
lgd = legend(arrayfun(@(s) sprintf('Iteration %d', s), 1:nSteps, 'UniformOutput', false), ...
             'Location', 'northeast');
lgd.Box = 'off';
lgd.ItemTokenSize = [10, 6];

%% Figure 4: eigenvalues of K_eta vs F_eta, and residual-Hessian ratio R_p^RH
lamK = EPeta.lamK;                          % eigenvalues of K_eta (descending)
lamF = EPeta.lamF;                          % eigenvalues of F_eta (descending)
RRH  = [M.R_RH];                            % R_p^RH = L_p^K / L_p^F; NaN = non-identifiable
nEig = numel(lamK);

figure('Name', 'Eigenvalues and R_RH', 'Color', 'w');
set(gcf, 'position', [203, 221, 1000, 400]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (i) eigenvalue spectra (the smallest eigenvalue is ~ -1e-12, i.e. numerically zero, so plot magnitudes)
nexttile; hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
set(gca, 'YScale', 'log');
hK = plot(1:nEig, abs(lamK), 'o-', 'Color', blueCol, 'LineWidth', 1.5, ...
          'MarkerSize', 8, 'MarkerFaceColor', blueCol);
hF = plot(1:nEig, abs(lamF), 's--', 'Color', orangeCol, 'LineWidth', 1.5, ...
          'MarkerSize', 8, 'MarkerFaceColor', orangeCol);
xlim([0.5 nEig + 0.5]);
ylim([1e-14 1e6]);
yticks(10.^(-14:4:6));
xticks(1:nEig);
xlabel('Eigenvalue index i');
ylabel('|\lambda_i|');
title('(i)');
lgd = legend([hK hF], {'K_\eta', 'F_\eta'}, 'Location', 'northeast');
lgd.Box = 'off';
lgd.ItemTokenSize = [10, 6];

% (ii) R_p^RH per parameter (< 1: the residual Hessian reduces the coordinate Lipschitz coefficient)
nexttile; hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
bar(x, RRH, 0.6, 'FaceColor', blueCol);
yline(1, '--k', 'LineWidth', 1.5);
ylim([0 1.2]);
text(x(isnan(RRH)), zeros(1, nnz(isnan(RRH))), 'n/a', 'Color', [0.5 0.5 0.5], ...
     'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
xlim([0.5 nParams + 0.5]);
xticks(x); xticklabels(dispNames);
xlabel('Parameter');
ylabel('R_p^{RH}');
title('(ii)');

%% Figure 5: Yao orthogonal method on the identifiable combinations
comboDisp = {'p_1', 'p_2x_1(0)', 'p_3', 'p_4', 's_1x_1(0)', 's_2x_1(0)'};   % same order as comboNames
xc     = 1:nCombo;
orderC = double(YaoC.order_1based);         % combinations in selection order
RllC   = YaoC.Rll;                          % normalized orthogonal residual |R_LL| at each step
isEstC = xc <= double(YaoC.n_estimable);

figure('Name', 'Yao method (combinations)', 'Color', 'w');
set(gcf, 'position', [290, 357, 830, 420]);
hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
set(gca, 'YScale', 'log');
hBars = bar(xc(isEstC), RllC(isEstC), 0.6, 'FaceColor', blueCol, 'BaseValue', 1e-2);
lgdTxt = {'Estimable'};
if any(~isEstC)                             % all 6 combinations are estimable in the current data
    hBars(end+1) = bar(xc(~isEstC), RllC(~isEstC), 0.6, 'FaceColor', redCol, 'BaseValue', 1e-2);
    lgdTxt{end+1} = 'Excluded';
end
hEps = yline(YaoC.eps, '--k', 'LineWidth', 1.5);
ylim([1e-2 10]);
yticks(10.^(-2:1:1));
xticks(xc); xticklabels(comboDisp(orderC));
xlabel('Parameter (Yao ranking)');
ylabel('|R_{LL}|');
lgd = legend([hBars hEps], [lgdTxt, {sprintf('\\epsilon = %.2g', YaoC.eps)}], 'Location', 'eastoutside');
lgd.Box = 'off';
lgd.ItemTokenSize = [10, 6];

%% Figure 6: Quaiser & Moenkemoeller eigenvalue method on the identifiable combinations
comboKeys = comboKey(comboNames);           % e.g. 'p2x1_0' -> 'p2x10', to match 'ln(p2*x1_0)'
stepsC    = QuaC.steps;
if iscell(stepsC), stepsC = [stepsC{:}]; end
nStepsC   = numel(stepsC);
lamMinC   = [stepsC.lam_min];
U1C       = nan(nCombo, nStepsC);           % |u_1| component per combination (NaN = already removed)
stepLblC  = strings(1, nStepsC);
for s = 1:nStepsC
    % u1 field names such as 'ln(p2*x1_0)' are not valid identifiers, so read them via fieldnames
    u1Names = fieldnames(stepsC(s).u1);
    u1Vals  = cell2mat(struct2cell(stepsC(s).u1));
    for j = 1:numel(u1Names)
        U1C(strcmp(comboKeys, comboKey(u1Names(j))), s) = abs(u1Vals(j));
    end
    if isempty(stepsC(s).removed)
        stepLblC(s) = "stop";
    else
        stepLblC(s) = comboDisp{strcmp(comboKeys, comboKey({stepsC(s).removed}))} + " removed";
    end
end

figure('Name', 'Quaiser method (combinations)', 'Color', 'w');
set(gcf, 'position', [203, 221, 1000, 400]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (i) lambda_min per iteration, same plotted threshold epsQ as Figure 3 (saved eps = 6.97 gives the same removals)
nexttile; hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
set(gca, 'YScale', 'log');
hLam = plot(1:nStepsC, abs(lamMinC), 'o-', 'Color', blueCol, 'LineWidth', 1.5, ...
            'MarkerSize', 8, 'MarkerFaceColor', blueCol);
hEps = yline(epsQ, '--r', 'LineWidth', 1.5);
text(1:nStepsC, abs(lamMinC)/8, stepLblC, 'FontSize', 16, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
xlim([0.4 nStepsC + 0.6]);
ylim([1e-10 1e4]);
yticks(10.^(-10:2:4));
xticks(1:nStepsC);
xlabel('Iteration');
ylabel('|\lambda_{min}|');
title('(i)');
lgd = legend([hLam hEps], {'|\lambda_{min}|', sprintf('\\epsilon = %g', epsQ)}, 'Location', 'southeast');
lgd.Box = 'off';
lgd.ItemTokenSize = [10, 6];

% (ii) components of the smallest eigenvector u_1 at each iteration
nexttile; hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);
bar(xc, U1C, 'grouped');
ylim([0 1.2]);
xticks(xc); xticklabels(comboDisp);
xlabel('Parameter');
ylabel('|u_{1,k}|');
title('(ii)');
lgd = legend(arrayfun(@(s) sprintf('Iteration %d', s), 1:nStepsC, 'UniformOutput', false), ...
             'Location', 'northeast');
lgd.Box = 'off';
lgd.ItemTokenSize = [10, 6];

%% Figure 7: coordinate Lipschitz coefficients L_p^K and L_p^F of the identifiable combinations
Mc = cellfun(@(n) CMc.(n), comboNames, 'UniformOutput', false);
Mc = [Mc{:}];                               % 1 x 6 struct array, same order as comboDisp
LKc = [Mc.L_K];
LFc = [Mc.L_F];
% All combinations are identifiable (finite L), but L_F spans 0.04 to ~840, so both axes use the same log range
yLimC = [1e-2 1e4];

figure('Name', 'Coordinate metrics (combinations)', 'Color', 'w');
set(gcf, 'position', [203, 221, 800, 420]);
hold on; box off;
set(gca, 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold', 'linewidth', 1.5);

% Left axis: L_p^K
yyaxis left; hold on;
set(gca, 'YScale', 'log');
hLK = bar(xc - w/2, LKc, w, 'FaceColor', blueCol, 'BaseValue', yLimC(1));
ylim(yLimC);
yticks(10.^(-2:2:4));
ylabel('L_p^K');

% Right axis: L_p^F
yyaxis right; hold on;
set(gca, 'YScale', 'log');
hLF = bar(xc + w/2, LFc, w, 'FaceColor', orangeCol, 'BaseValue', yLimC(1));
ylim(yLimC);
yticks(10.^(-2:2:4));
ylabel('L_p^F');

ax = gca;
ax.YAxis(1).Color = blueCol;
ax.YAxis(2).Color = orangeCol;
xlim([0.5 nCombo + 0.5]);
xticks(xc); xticklabels(comboDisp);
xlabel('Parameter');
lgd = legend([hLK hLF], {'L_p^K', 'L_p^F'}, 'Location', 'northoutside', 'Orientation', 'horizontal');
lgd.Box = 'off';
lgd.ItemTokenSize = [10, 6];

%% Local functions
function str = joinNames(c)
    % Join a cell array of names; empty (e.g. no excluded parameters) -> '(none)'
    if isempty(c)
        str = '(none)';
    else
        str = strjoin(cellstr(c), ', ');
    end
end

function key = comboKey(names)
    % Normalize a parameter/combination name for matching: 'ln(p2*x1_0)', 'ln_p1', 'p2x1_0' -> 'p2x10', 'p1'
    key = regexprep(lower(cellstr(names)), '[^a-z0-9]', '');
    key = regexprep(key, '^ln', '');
end
