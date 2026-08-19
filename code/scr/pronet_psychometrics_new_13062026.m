% =========================================================================
% PRONET MASTER PIPELINE: Fitting, Plotting & Statistics
% =========================================================================

clear; clc; close all;

fprintf('======================================================\n');
fprintf('  STARTING PRONET MASTER PIPELINE                     \n');
fprintf('======================================================\n\n');

addpath('/mnt/beegfs/users/antonia.ceric/git/pronet/code/matlab/func');

global HPC_PATH DATA_PATH
HPC_PATH  = '/mnt/beegfs/workspace/2024-0404-PRONET/';
DATA_PATH = fullfile(HPC_PATH, 'DATA', 'proc');

FIT_DIR = fullfile(DATA_PATH, 'psychometric_fits_062026');
if ~exist(FIT_DIR, 'dir'), mkdir(FIT_DIR); end

STATS_OUTDIR = fullfile(FIT_DIR, 'stats_plots');
if ~exist(STATS_OUTDIR, 'dir'), mkdir(STATS_OUTDIR); end

%% ========================================================================
%  STEP 1: RUN DATA FITTING
%  ========================================================================
fprintf('>>> STEP 1: Fitting Behavioral Data...\n\n');
try
    fprintf('--- Fitting Intonation Task ---\n'); 
    pronet_fitting_intonation_custom_ml();
catch ME
    fprintf('WARNING in Intonation fit: %s\n', ME.message); 
end

try
    fprintf('--- Fitting Phoneme Task ---\n'); 
    pronet_fitting_phoneme_custom_ml();
catch ME
    fprintf('WARNING in Phoneme fit: %s\n', ME.message); 
end

%% ========================================================================
%  STEP 2: GENERATE OVERLAY PLOTS
%  ========================================================================
fprintf('\n>>> STEP 2: Generating Overlay Plots...\n\n');
try
    plot_pronet_fitted_curves('intonation'); 
catch ME
    fprintf('WARNING plotting Intonation: %s\n', ME.message); 
end

try
    plot_pronet_fitted_curves('phoneme');    
catch ME
    fprintf('WARNING plotting Phoneme: %s\n', ME.message); 
end

%% ========================================================================
%  STEP 3: STATISTICAL ANALYSIS & PLOTS
%  ========================================================================
fprintf('\n>>> STEP 3: Running Statistical Analysis...\n\n');

col_pho = [0 0 0.8]; col_int = [0.8 0 0];
VP = {'02','03','04','05','06','07','08','09','10','11','12','13','15','16','17',...
      '18','19','20','21','22','23','24','25','26','27','28','29','31','32'};

% ---- Manual exclusion lists (23-Jun-2025) ----
manual.phoneme.all = {'05','11','24','27','32'};
manual.phoneme.rt  = {'26'};
manual.phoneme.vx  = {'29'};
manual.phoneme.lt  = {};

manual.intonation.all = {'13'};
manual.intonation.rt  = {'26'};
manual.intonation.vx  = {'29'};
manual.intonation.lt  = {'32'};

% IMPORTANT: Number of trials per level for exact LL calculation
N_TRIALS_SCANNER = 62; 
x_pos = [1, 2, 4, 5]; w = 0.6; jitter = 0.15;
sem = @(x) std(x)/sqrt(length(x));

fprintf('=======================================================\n');
fprintf('=== PART 3A: SCANNER DATA (VX CONDITION)            ===\n');
fprintf('=======================================================\n');

load(fullfile(FIT_DIR, 'intonation_resp_fit_all_conditions.mat'), 'AllFitInfo'); 
Fit_Intonation = [AllFitInfo.vx.average]; 

load(fullfile(FIT_DIR, 'phoneme_resp_fit_all_conditions.mat'), 'AllFitInfo');  
Fit_Phoneme = [AllFitInfo.vx.average]; 

incl_pho = setdiff(unique({Fit_Phoneme.subject}), unique([manual.phoneme.all, manual.phoneme.vx]));
incl_int = setdiff(unique({Fit_Intonation.subject}), unique([manual.intonation.all, manual.intonation.vx]));

% Calculate R2, Adj R2, and AIC
[~, r2_lin_pho_vx, r2_sig_pho_vx, adj_lin_pho_vx, adj_sig_pho_vx, aic_lin_pho_vx, aic_sig_pho_vx] = calculate_scanner_metrics(Fit_Phoneme, @(s) ismember(s, incl_pho), N_TRIALS_SCANNER);
[~, r2_lin_int_vx, r2_sig_int_vx, adj_lin_int_vx, adj_sig_int_vx, aic_lin_int_vx, aic_sig_int_vx] = calculate_scanner_metrics(Fit_Intonation, @(s) ismember(s, incl_int), N_TRIALS_SCANNER);

% T-Tests for R2 (Two-Tailed)
[~, p_pho_vx_r2]     = ttest(r2_lin_pho_vx, r2_sig_pho_vx);
[~, p_int_vx_r2]     = ttest(r2_lin_int_vx, r2_sig_int_vx);
[~, p_pho_vx_adj]    = ttest(adj_lin_pho_vx, adj_sig_pho_vx);
[~, p_int_vx_adj]    = ttest(adj_lin_int_vx, adj_sig_int_vx);

% Standard Two-Tailed T-Tests for AIC
[~, p_pho_vx_aic]    = ttest(aic_lin_pho_vx, aic_sig_pho_vx);
[~, p_int_vx_aic]    = ttest(aic_lin_int_vx, aic_sig_int_vx);

% Robustness Checks for AIC (Hypothesis: Linear > Sigmoid)
[~, p_pho_vx_aic_1tail] = ttest(aic_lin_pho_vx, aic_sig_pho_vx, 'Tail', 'right');
[~, p_int_vx_aic_1tail] = ttest(aic_lin_int_vx, aic_sig_int_vx, 'Tail', 'right');
p_pho_vx_wilcox = signrank(aic_lin_pho_vx, aic_sig_pho_vx, 'tail', 'right');
p_int_vx_wilcox = signrank(aic_lin_int_vx, aic_sig_int_vx, 'tail', 'right');

fprintf('--- STANDARD R-SQUARED (Higher is better) ---\n');
fprintf('VX Phoneme    | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n', mean(r2_lin_pho_vx), mean(r2_sig_pho_vx), p_pho_vx_r2);
fprintf('VX Intonation | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n\n', mean(r2_lin_int_vx), mean(r2_sig_int_vx), p_int_vx_r2);

fprintf('--- ADJUSTED R-SQUARED (Higher is better) ---\n');
fprintf('VX Phoneme    | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n', mean(adj_lin_pho_vx), mean(adj_sig_pho_vx), p_pho_vx_adj);
fprintf('VX Intonation | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n\n', mean(adj_lin_int_vx), mean(adj_sig_int_vx), p_int_vx_adj);

fprintf('--- AIC (Lower is better) | STANDARD TWO-TAILED TEST ---\n');
fprintf('VX Phoneme    | Linear M=%.2f vs Sigmoid M=%.2f | p = %.4f\n', mean(aic_lin_pho_vx), mean(aic_sig_pho_vx), p_pho_vx_aic);
fprintf('VX Intonation | Linear M=%.2f vs Sigmoid M=%.2f | p = %.4f\n\n', mean(aic_lin_int_vx), mean(aic_sig_int_vx), p_int_vx_aic);

fprintf('--- AIC ROBUSTNESS CHECKS (Hypothesis: Linear > Sigmoid) ---\n');
fprintf('VX Phoneme    | One-Tailed p = %.4f | Wilcoxon Sign-Rank p = %.4f\n', p_pho_vx_aic_1tail, p_pho_vx_wilcox);
fprintf('VX Intonation | One-Tailed p = %.4f | Wilcoxon Sign-Rank p = %.4f\n\n', p_int_vx_aic_1tail, p_int_vx_wilcox);

% PLOT 1: Scanner Standard R2
fig1 = figure('Name', 'Scanner VX: Standard R2', 'Color', 'w', 'Position', [100 200 700 600]); hold on;
bar(x_pos(1), mean(r2_lin_pho_vx), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(2), mean(r2_lin_int_vx), w, 'FaceColor', col_int, 'EdgeColor', 'none');
bar(x_pos(3), mean(r2_sig_pho_vx), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(4), mean(r2_sig_int_vx), w, 'FaceColor', col_int, 'EdgeColor', 'none');
errorbar(x_pos, [mean(r2_lin_pho_vx), mean(r2_lin_int_vx), mean(r2_sig_pho_vx), mean(r2_sig_int_vx)], ...
    [sem(r2_lin_pho_vx), sem(r2_lin_int_vx), sem(r2_sig_pho_vx), sem(r2_sig_int_vx)], 'k.', 'LineWidth', 1.5, 'CapSize', 0);
plot(x_pos(1) + randn(size(r2_lin_pho_vx))*jitter, r2_lin_pho_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(2) + randn(size(r2_lin_int_vx))*jitter, r2_lin_int_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
plot(x_pos(3) + randn(size(r2_sig_pho_vx))*jitter, r2_sig_pho_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(4) + randn(size(r2_sig_int_vx))*jitter, r2_sig_int_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
y_max = 1.05; y_max2 = 1.15;
plot([x_pos(1), x_pos(3)], [y_max y_max], 'k-', 'LineWidth', 1.5); text(2.5, y_max + 0.03, get_pval_stars(p_pho_vx_r2), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_pho);
plot([x_pos(2), x_pos(4)], [y_max2 y_max2], 'k-', 'LineWidth', 1.5); text(3.5, y_max2 + 0.03, get_pval_stars(p_int_vx_r2), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_int);
title('Scanner (VX) - Model Comparison via Standard R^2', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Model Fit (Standard R^2)', 'FontSize', 14); ylim([0 1.25]); xlim([0 6]); xticks([1.5, 4.5]); xticklabels({'Linear Fit', 'Sigmoid Fit'});
set(gca, 'FontSize', 12, 'TickDir', 'out'); grid on; box on;
legend([bar(NaN, NaN, 'FaceColor', col_pho, 'EdgeColor', 'none'), bar(NaN, NaN, 'FaceColor', col_int, 'EdgeColor', 'none')], ...
    {'Phoneme', 'Intonation'}, 'Location', 'southeast', 'FontSize', 12, 'Box', 'off');
saveas(fig1, fullfile(STATS_OUTDIR, 'vx_linear_vs_sigmoid_R2_bar.png')); close(fig1);

% PLOT 2: Scanner AIC (Uses Parametric p-values for plotting)
fig2 = figure('Name', 'Scanner VX: AIC', 'Color', 'w', 'Position', [150 200 700 600]); hold on;
bar(x_pos(1), mean(aic_lin_pho_vx), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(2), mean(aic_lin_int_vx), w, 'FaceColor', col_int, 'EdgeColor', 'none');
bar(x_pos(3), mean(aic_sig_pho_vx), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(4), mean(aic_sig_int_vx), w, 'FaceColor', col_int, 'EdgeColor', 'none');
errorbar(x_pos, [mean(aic_lin_pho_vx), mean(aic_lin_int_vx), mean(aic_sig_pho_vx), mean(aic_sig_int_vx)], ...
    [sem(aic_lin_pho_vx), sem(aic_lin_int_vx), sem(aic_sig_pho_vx), sem(aic_sig_int_vx)], 'k.', 'LineWidth', 1.5, 'CapSize', 0);
plot(x_pos(1) + randn(size(aic_lin_pho_vx))*jitter, aic_lin_pho_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(2) + randn(size(aic_lin_int_vx))*jitter, aic_lin_int_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
plot(x_pos(3) + randn(size(aic_sig_pho_vx))*jitter, aic_sig_pho_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(4) + randn(size(aic_sig_int_vx))*jitter, aic_sig_int_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
max_aic_vx = max([aic_lin_pho_vx(:); aic_sig_pho_vx(:); aic_lin_int_vx(:); aic_sig_int_vx(:)]); 
y_max_aic = max_aic_vx * 1.05; y_max2_aic = max_aic_vx * 1.15;

% USING PARAMETRIC T-TEST P-VALUES NOW
plot([x_pos(1), x_pos(3)], [y_max_aic y_max_aic], 'k-', 'LineWidth', 1.5); text(2.5, y_max_aic + (max_aic_vx*0.02), get_pval_stars(p_pho_vx_aic), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_pho);
plot([x_pos(2), x_pos(4)], [y_max2_aic y_max2_aic], 'k-', 'LineWidth', 1.5); text(3.5, y_max2_aic + (max_aic_vx*0.02), get_pval_stars(p_int_vx_aic), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_int);

title('Scanner (VX) - Model Comparison via AIC', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Information Criterion (AIC)', 'FontSize', 14); 
ylim([0 max_aic_vx * 1.5]); xlim([0 6]); xticks([1.5, 4.5]); xticklabels({'Linear Fit', 'Sigmoid Fit'});
set(gca, 'FontSize', 12, 'TickDir', 'out'); grid on; box on;
legend([bar(NaN, NaN, 'FaceColor', col_pho, 'EdgeColor', 'none'), bar(NaN, NaN, 'FaceColor', col_int, 'EdgeColor', 'none')], ...
    {'Phoneme', 'Intonation'}, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');
saveas(fig2, fullfile(STATS_OUTDIR, 'vx_linear_vs_sigmoid_AIC_bar.png')); close(fig2);

% % PLOT 2: Scanner AIC (Uses Wilcoxon p-values for plotting)
% fig2 = figure('Name', 'Scanner VX: AIC', 'Color', 'w', 'Position', [150 200 700 600]); hold on;
% bar(x_pos(1), mean(aic_lin_pho_vx), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
% bar(x_pos(2), mean(aic_lin_int_vx), w, 'FaceColor', col_int, 'EdgeColor', 'none');
% bar(x_pos(3), mean(aic_sig_pho_vx), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
% bar(x_pos(4), mean(aic_sig_int_vx), w, 'FaceColor', col_int, 'EdgeColor', 'none');
% errorbar(x_pos, [mean(aic_lin_pho_vx), mean(aic_lin_int_vx), mean(aic_sig_pho_vx), mean(aic_sig_int_vx)], ...
%     [sem(aic_lin_pho_vx), sem(aic_lin_int_vx), sem(aic_sig_pho_vx), sem(aic_sig_int_vx)], 'k.', 'LineWidth', 1.5, 'CapSize', 0);
% plot(x_pos(1) + randn(size(aic_lin_pho_vx))*jitter, aic_lin_pho_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
% plot(x_pos(2) + randn(size(aic_lin_int_vx))*jitter, aic_lin_int_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
% plot(x_pos(3) + randn(size(aic_sig_pho_vx))*jitter, aic_sig_pho_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
% plot(x_pos(4) + randn(size(aic_sig_int_vx))*jitter, aic_sig_int_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
% max_aic_vx = max([aic_lin_pho_vx(:); aic_sig_pho_vx(:); aic_lin_int_vx(:); aic_sig_int_vx(:)]); 
% y_max_aic = max_aic_vx * 1.05; y_max2_aic = max_aic_vx * 1.15;
% % USING WILCOXON P-VALUES HERE
% plot([x_pos(1), x_pos(3)], [y_max_aic y_max_aic], 'k-', 'LineWidth', 1.5); text(2.5, y_max_aic + (max_aic_vx*0.02), get_pval_stars(p_pho_vx_wilcox), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_pho);
% plot([x_pos(2), x_pos(4)], [y_max2_aic y_max2_aic], 'k-', 'LineWidth', 1.5); text(3.5, y_max2_aic + (max_aic_vx*0.02), get_pval_stars(p_int_vx_wilcox), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_int);
% title('Scanner (VX) - Model Comparison via AIC', 'FontSize', 16, 'FontWeight', 'bold');
% ylabel('Information Criterion (AIC)', 'FontSize', 14); 
% ylim([0 max_aic_vx * 1.5]); xlim([0 6]); xticks([1.5, 4.5]); xticklabels({'Linear Fit', 'Sigmoid Fit'});
% set(gca, 'FontSize', 12, 'TickDir', 'out'); grid on; box on;
% legend([bar(NaN, NaN, 'FaceColor', col_pho, 'EdgeColor', 'none'), bar(NaN, NaN, 'FaceColor', col_int, 'EdgeColor', 'none')], ...
%     {'Phoneme (Wilcoxon)', 'Intonation (Wilcoxon)'}, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');
% saveas(fig2, fullfile(STATS_OUTDIR, 'vx_linear_vs_sigmoid_AIC_bar.png')); close(fig2);


fprintf('\n=======================================================\n');
fprintf('=== PART 3B: PRETEST DATA (R2 & AIC)                ===\n');
fprintf('=======================================================\n');

MOCS_DIR = fullfile(HPC_PATH, 'TMS/Pre-Session/results/');
r2_lin_pho_pre = []; r2_sig_pho_pre = []; r2_lin_int_pre = []; r2_sig_int_pre = [];
adj_lin_pho_pre = []; adj_sig_pho_pre = []; adj_lin_int_pre = []; adj_sig_int_pre = [];
aic_lin_pho_pre = []; aic_sig_pho_pre = []; aic_lin_int_pre = []; aic_sig_int_pre = [];
x_levels = 1:5;

for i = 1:length(VP)
    sub_id = VP{i};
    fpath = fullfile(MOCS_DIR, sprintf('%s_PRONET_MOCS_results.mat', sub_id));
    if ~exist(fpath, 'file'), continue; end
    tmp = load(fpath);
    tasks_pre = {'phoneme', 'intonation'}; file_names = {'phoneme', 'intonat'};
    
    for t = 1:2
        task = tasks_pre{t}; task_in = file_names{t};
        if ismember(sub_id, manual.(task).all), continue; end
        if ~isfield(tmp, 'm1') || ~isfield(tmp.m1, task_in), continue; end
        
        yes_tot = tmp.m1.(task_in).yes + tmp.f2.(task_in).yes;
        n_tot   = tmp.m1.(task_in).n + tmp.f2.(task_in).n;
        y_true  = yes_tot ./ n_tot;
        
        % LINEAR FIT
        p_lin = polyfit(x_levels, y_true, 1); y_lin = polyval(p_lin, x_levels);
        r2_lin = calc_r2(y_true, y_lin);
        adj_lin = 1 - ((1 - r2_lin) * (5 - 1) / (5 - 2 - 1));
        
        y_lin_clamped = max(min(y_lin, 0.999), 0.001);
        ll_lin = sum(yes_tot .* log(y_lin_clamped) + (n_tot - yes_tot) .* log(1 - y_lin_clamped));
        aic_lin = 4 - 2 * ll_lin;
        
        % SIGMOID FIT
        p_m1 = tmp.m1.(task_in).PFfit.params; p_f2 = tmp.f2.(task_in).PFfit.params;
        y_sig_m1 = eval_logistic(p_m1, x_levels); y_sig_f2 = eval_logistic(p_f2, x_levels);
        y_sig = (y_sig_m1 + y_sig_f2)/2;
        r2_sig = calc_r2(y_true, y_sig);
        adj_sig = 1 - ((1 - r2_sig) * (5 - 1) / (5 - 2 - 1));
        
        y_sig_clamped = max(min(y_sig, 0.999), 0.001);
        ll_sig = sum(yes_tot .* log(y_sig_clamped) + (n_tot - yes_tot) .* log(1 - y_sig_clamped));
        aic_sig = 4 - 2 * ll_sig;
        
        if strcmp(task, 'phoneme')
            r2_lin_pho_pre(end+1) = r2_lin; r2_sig_pho_pre(end+1) = r2_sig;
            adj_lin_pho_pre(end+1) = adj_lin; adj_sig_pho_pre(end+1) = adj_sig;
            aic_lin_pho_pre(end+1) = aic_lin; aic_sig_pho_pre(end+1) = aic_sig;
        else
            r2_lin_int_pre(end+1) = r2_lin; r2_sig_int_pre(end+1) = r2_sig;
            adj_lin_int_pre(end+1) = adj_lin; adj_sig_int_pre(end+1) = adj_sig;
            aic_lin_int_pre(end+1) = aic_lin; aic_sig_int_pre(end+1) = aic_sig;
        end
    end
end

[~, p_pho_pre_r2]  = ttest(r2_lin_pho_pre, r2_sig_pho_pre);
[~, p_int_pre_r2]  = ttest(r2_lin_int_pre, r2_sig_int_pre);
[~, p_pho_pre_adj] = ttest(adj_lin_pho_pre, adj_sig_pho_pre);
[~, p_int_pre_adj] = ttest(adj_lin_int_pre, adj_sig_int_pre);
[~, p_pho_pre_aic] = ttest(aic_lin_pho_pre, aic_sig_pho_pre);
[~, p_int_pre_aic] = ttest(aic_lin_int_pre, aic_sig_int_pre);

p_pho_pre_wilcox = signrank(aic_lin_pho_pre, aic_sig_pho_pre, 'tail', 'right');
p_int_pre_wilcox = signrank(aic_lin_int_pre, aic_sig_int_pre, 'tail', 'right');

fprintf('--- STANDARD R-SQUARED (Higher is better) ---\n');
fprintf('Pretest Phoneme    | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n', mean(r2_lin_pho_pre), mean(r2_sig_pho_pre), p_pho_pre_r2);
fprintf('Pretest Intonation | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n\n', mean(r2_lin_int_pre), mean(r2_sig_int_pre), p_int_pre_r2);

fprintf('--- ADJUSTED R-SQUARED (Higher is better) ---\n');
fprintf('Pretest Phoneme    | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n', mean(adj_lin_pho_pre), mean(adj_sig_pho_pre), p_pho_pre_adj);
fprintf('Pretest Intonation | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n\n', mean(adj_lin_int_pre), mean(adj_sig_int_pre), p_int_pre_adj);

fprintf('--- AIC (Lower is better) ---\n');
fprintf('Pretest Phoneme    | Linear M=%.2f vs Sigmoid M=%.2f | p = %.4f\n', mean(aic_lin_pho_pre), mean(aic_sig_pho_pre), p_pho_pre_aic);
fprintf('Pretest Intonation | Linear M=%.2f vs Sigmoid M=%.2f | p = %.4f\n\n', mean(aic_lin_int_pre), mean(aic_sig_int_pre), p_int_pre_aic);

% PLOT 3: Pretest Standard R2
fig3 = figure('Name', 'Pretest: Standard R2', 'Color', 'w', 'Position', [200 250 700 600]); hold on;
bar(x_pos(1), mean(r2_lin_pho_pre), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(2), mean(r2_lin_int_pre), w, 'FaceColor', col_int, 'EdgeColor', 'none');
bar(x_pos(3), mean(r2_sig_pho_pre), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(4), mean(r2_sig_int_pre), w, 'FaceColor', col_int, 'EdgeColor', 'none');
errorbar(x_pos, [mean(r2_lin_pho_pre), mean(r2_lin_int_pre), mean(r2_sig_pho_pre), mean(r2_sig_int_pre)], ...
    [sem(r2_lin_pho_pre), sem(r2_lin_int_pre), sem(r2_sig_pho_pre), sem(r2_sig_int_pre)], 'k.', 'LineWidth', 1.5, 'CapSize', 0);
plot(x_pos(1) + randn(size(r2_lin_pho_pre))*jitter, r2_lin_pho_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(2) + randn(size(r2_lin_int_pre))*jitter, r2_lin_int_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
plot(x_pos(3) + randn(size(r2_sig_pho_pre))*jitter, r2_sig_pho_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(4) + randn(size(r2_sig_int_pre))*jitter, r2_sig_int_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
plot([x_pos(1), x_pos(3)], [y_max y_max], 'k-', 'LineWidth', 1.5); text(2.5, y_max + 0.03, get_pval_stars(p_pho_pre_r2), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_pho);
plot([x_pos(2), x_pos(4)], [y_max2 y_max2], 'k-', 'LineWidth', 1.5); text(3.5, y_max2 + 0.03, get_pval_stars(p_int_pre_r2), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_int);
title('Pretest - Model Comparison via Standard R^2', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Model Fit (Standard R^2)', 'FontSize', 14); ylim([0 1.25]); xlim([0 6]); xticks([1.5, 4.5]); xticklabels({'Linear Fit', 'Sigmoid Fit'});
set(gca, 'FontSize', 12, 'TickDir', 'out'); grid on; box on;
legend([bar(NaN, NaN, 'FaceColor', col_pho, 'EdgeColor', 'none'), bar(NaN, NaN, 'FaceColor', col_int, 'EdgeColor', 'none')], ...
    {'Phoneme', 'Intonation'}, 'Location', 'southeast', 'FontSize', 12, 'Box', 'off');
saveas(fig3, fullfile(STATS_OUTDIR, 'pretest_linear_vs_sigmoid_R2_bar.png')); close(fig3);

% PLOT 4: Pretest AIC (Uses Parametric p-values)
fig4 = figure('Name', 'Pretest: AIC', 'Color', 'w', 'Position', [250 250 700 600]); hold on;
bar(x_pos(1), mean(aic_lin_pho_pre), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(2), mean(aic_lin_int_pre), w, 'FaceColor', col_int, 'EdgeColor', 'none');
bar(x_pos(3), mean(aic_sig_pho_pre), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(4), mean(aic_sig_int_pre), w, 'FaceColor', col_int, 'EdgeColor', 'none');
errorbar(x_pos, [mean(aic_lin_pho_pre), mean(aic_lin_int_pre), mean(aic_sig_pho_pre), mean(aic_sig_int_pre)], ...
    [sem(aic_lin_pho_pre), sem(aic_lin_int_pre), sem(aic_sig_pho_pre), sem(aic_sig_int_pre)], 'k.', 'LineWidth', 1.5, 'CapSize', 0);
plot(x_pos(1) + randn(size(aic_lin_pho_pre))*jitter, aic_lin_pho_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(2) + randn(size(aic_lin_int_pre))*jitter, aic_lin_int_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
plot(x_pos(3) + randn(size(aic_sig_pho_pre))*jitter, aic_sig_pho_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(4) + randn(size(aic_sig_int_pre))*jitter, aic_sig_int_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
max_aic_pre = max([aic_lin_pho_pre(:); aic_sig_pho_pre(:); aic_lin_int_pre(:); aic_sig_int_pre(:)]); 
y_max_aic_pre = max_aic_pre * 1.05; y_max2_aic_pre = max_aic_pre * 1.15;

% USING PARAMETRIC T-TEST P-VALUES NOW
plot([x_pos(1), x_pos(3)], [y_max_aic_pre y_max_aic_pre], 'k-', 'LineWidth', 1.5); text(2.5, y_max_aic_pre + (max_aic_pre*0.02), get_pval_stars(p_pho_pre_aic), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_pho);
plot([x_pos(2), x_pos(4)], [y_max2_aic_pre y_max2_aic_pre], 'k-', 'LineWidth', 1.5); text(3.5, y_max2_aic_pre + (max_aic_pre*0.02), get_pval_stars(p_int_pre_aic), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_int);

title('Pretest - Model Comparison via AIC', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Information Criterion (AIC)', 'FontSize', 14); 
ylim([0 max_aic_pre * 1.5]); xlim([0 6]); xticks([1.5, 4.5]); xticklabels({'Linear Fit', 'Sigmoid Fit'});
set(gca, 'FontSize', 12, 'TickDir', 'out'); grid on; box on;
legend([bar(NaN, NaN, 'FaceColor', col_pho, 'EdgeColor', 'none'), bar(NaN, NaN, 'FaceColor', col_int, 'EdgeColor', 'none')], ...
    {'Phoneme', 'Intonation'}, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');
saveas(fig4, fullfile(STATS_OUTDIR, 'pretest_linear_vs_sigmoid_AIC_bar.png')); close(fig4);

% % PLOT 4: Pretest AIC (Uses Wilcoxon p-values)
% fig4 = figure('Name', 'Pretest: AIC', 'Color', 'w', 'Position', [250 250 700 600]); hold on;
% bar(x_pos(1), mean(aic_lin_pho_pre), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
% bar(x_pos(2), mean(aic_lin_int_pre), w, 'FaceColor', col_int, 'EdgeColor', 'none');
% bar(x_pos(3), mean(aic_sig_pho_pre), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
% bar(x_pos(4), mean(aic_sig_int_pre), w, 'FaceColor', col_int, 'EdgeColor', 'none');
% errorbar(x_pos, [mean(aic_lin_pho_pre), mean(aic_lin_int_pre), mean(aic_sig_pho_pre), mean(aic_sig_int_pre)], ...
%     [sem(aic_lin_pho_pre), sem(aic_lin_int_pre), sem(aic_sig_pho_pre), sem(aic_sig_int_pre)], 'k.', 'LineWidth', 1.5, 'CapSize', 0);
% plot(x_pos(1) + randn(size(aic_lin_pho_pre))*jitter, aic_lin_pho_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
% plot(x_pos(2) + randn(size(aic_lin_int_pre))*jitter, aic_lin_int_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
% plot(x_pos(3) + randn(size(aic_sig_pho_pre))*jitter, aic_sig_pho_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
% plot(x_pos(4) + randn(size(aic_sig_int_pre))*jitter, aic_sig_int_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
% max_aic_pre = max([aic_lin_pho_pre(:); aic_sig_pho_pre(:); aic_lin_int_pre(:); aic_sig_int_pre(:)]); 
% y_max_aic_pre = max_aic_pre * 1.05; y_max2_aic_pre = max_aic_pre * 1.15;
% % USING WILCOXON P-VALUES HERE
% plot([x_pos(1), x_pos(3)], [y_max_aic_pre y_max_aic_pre], 'k-', 'LineWidth', 1.5); text(2.5, y_max_aic_pre + (max_aic_pre*0.02), get_pval_stars(p_pho_pre_wilcox), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_pho);
% plot([x_pos(2), x_pos(4)], [y_max2_aic_pre y_max2_aic_pre], 'k-', 'LineWidth', 1.5); text(3.5, y_max2_aic_pre + (max_aic_pre*0.02), get_pval_stars(p_int_pre_wilcox), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_int);
% title('Pretest - Model Comparison via AIC', 'FontSize', 16, 'FontWeight', 'bold');
% ylabel('Information Criterion (AIC)', 'FontSize', 14); 
% ylim([0 max_aic_pre * 1.5]); xlim([0 6]); xticks([1.5, 4.5]); xticklabels({'Linear Fit', 'Sigmoid Fit'});
% set(gca, 'FontSize', 12, 'TickDir', 'out'); grid on; box on;
% legend([bar(NaN, NaN, 'FaceColor', col_pho, 'EdgeColor', 'none'), bar(NaN, NaN, 'FaceColor', col_int, 'EdgeColor', 'none')], ...
%     {'Phoneme (Wilcoxon)', 'Intonation (Wilcoxon)'}, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');
% saveas(fig4, fullfile(STATS_OUTDIR, 'pretest_linear_vs_sigmoid_AIC_bar.png')); close(fig4);


fprintf('\n=======================================================\n');
fprintf('=== PART 3C: 2x2 RM-ANOVA (TASK x SESSION)          ===\n');
fprintf('=======================================================\n');

load(fullfile(FIT_DIR, 'phoneme_resp_fit_all_conditions.mat'), 'AllFitInfo');  
AllFitInfo_phoneme = AllFitInfo;
load(fullfile(FIT_DIR, 'intonation_resp_fit_all_conditions.mat'), 'AllFitInfo'); 
AllFitInfo_intonation = AllFitInfo;

data_PSE = []; data_Slope = []; valid_subs = 0;

for i = 1:length(VP)
    sub_id = VP{i};
    if ismember(sub_id, manual.phoneme.all) || ismember(sub_id, manual.phoneme.vx) || ...
       ismember(sub_id, manual.intonation.all) || ismember(sub_id, manual.intonation.vx)
        continue; 
    end
    
    fpath = fullfile(MOCS_DIR, sprintf('%s_PRONET_MOCS_results.mat', sub_id));
    if ~exist(fpath, 'file'), continue; end
    tmp = load(fpath);
    if ~isfield(tmp, 'm1') || ~isfield(tmp.m1, 'phoneme') || ~isfield(tmp.m1, 'intonat'), continue; end
    
    pre_pho_pse   = (tmp.m1.phoneme.PFfit.params(1) + tmp.f2.phoneme.PFfit.params(1)) / 2;
    pre_pho_slope = (tmp.m1.phoneme.PFfit.params(2) + tmp.f2.phoneme.PFfit.params(2)) / 2;
    pre_int_pse   = (tmp.m1.intonat.PFfit.params(1) + tmp.f2.intonat.PFfit.params(1)) / 2;
    pre_int_slope = (tmp.m1.intonat.PFfit.params(2) + tmp.f2.intonat.PFfit.params(2)) / 2;
    
    idx_pho = find(strcmp({AllFitInfo_phoneme.vx.average.subject}, sub_id) | strcmp({AllFitInfo_phoneme.vx.average.subject}, ['sub-' sub_id]), 1);
    idx_int = find(strcmp({AllFitInfo_intonation.vx.average.subject}, sub_id) | strcmp({AllFitInfo_intonation.vx.average.subject}, ['sub-' sub_id]), 1);
    if isempty(idx_pho) || isempty(idx_int), continue; end
    
    vx_pho_pse   = AllFitInfo_phoneme.vx.average(idx_pho).mlfit.alpha;
    vx_pho_slope = AllFitInfo_phoneme.vx.average(idx_pho).mlfit.beta;
    vx_int_pse   = AllFitInfo_intonation.vx.average(idx_int).mlfit.alpha;
    vx_int_slope = AllFitInfo_intonation.vx.average(idx_int).mlfit.beta;
    
    data_PSE   = [data_PSE; pre_pho_pse, vx_pho_pse, pre_int_pse, vx_int_pse];
    data_Slope = [data_Slope; pre_pho_slope, vx_pho_slope, pre_int_slope, vx_int_slope];
    valid_subs = valid_subs + 1;
end

if valid_subs > 2
    VarNames = {'Pho_Pre', 'Pho_VX', 'Int_Pre', 'Int_VX'};
    T_PSE   = array2table(data_PSE, 'VariableNames', VarNames);
    T_Slope = array2table(data_Slope, 'VariableNames', VarNames);
    
    Task    = categorical([1 1 2 2]', [1 2], {'Phoneme', 'Intonation'});
    Session = categorical([1 2 1 2]', [1 2], {'Pretest', 'VX'});
    WithinDesign = table(Task, Session);
    
    rm_PSE   = fitrm(T_PSE, 'Pho_Pre-Int_VX ~ 1', 'WithinDesign', WithinDesign);
    rm_Slope = fitrm(T_Slope, 'Pho_Pre-Int_VX ~ 1', 'WithinDesign', WithinDesign);
    
    disp('--- RANOVA: PSE (Bias/Threshold) ---'); disp(ranova(rm_PSE, 'WithinModel', 'Task*Session'));
    disp('--- RANOVA: SLOPE ---'); disp(ranova(rm_Slope, 'WithinModel', 'Task*Session'));
end

fprintf('\n=======================================================\n');
fprintf('=== PART 3D: POST-HOC TESTS (PSE INTERACTION)       ===\n');
fprintf('=======================================================\n');

if valid_subs > 2
    shift_pho = data_PSE(:,2) - data_PSE(:,1);
    shift_int = data_PSE(:,4) - data_PSE(:,3);
    
    [~, p_post_pho, ~, stats_pho] = ttest(data_PSE(:,1), data_PSE(:,2));
    [~, p_post_int, ~, stats_int] = ttest(data_PSE(:,3), data_PSE(:,4));
    
    fprintf('--- PHONEME ---\n');
    fprintf('Pretest M = %.3f  |  Scanner (VX) M = %.3f\n', mean(data_PSE(:,1)), mean(data_PSE(:,2)));
    fprintf('Mean Shift = %+.3f  |  t(%d) = %.2f  |  p = %.4f\n\n', mean(shift_pho), stats_pho.df, stats_pho.tstat, p_post_pho);
    
    fprintf('--- INTONATION ---\n');
    fprintf('Pretest M = %.3f  |  Scanner (VX) M = %.3f\n', mean(data_PSE(:,3)), mean(data_PSE(:,4)));
    fprintf('Mean Shift = %+.3f  |  t(%d) = %.2f  |  p = %.4f\n\n', mean(shift_int), stats_int.df, stats_int.tstat, p_post_int);
end

fprintf('\n=======================================================\n');
fprintf('=== PART 3E: INTERACTION PLOT (PSE)                 ===\n');
fprintf('=======================================================\n');

if valid_subs > 2
    means = mean(data_PSE, 1);
    sems  = std(data_PSE, 0, 1) ./ sqrt(size(data_PSE, 1));

    pho_means = [means(1), means(2)]; pho_sems  = [sems(1),  sems(2)];
    int_means = [means(3), means(4)]; int_sems  = [sems(3),  sems(4)];
    x_vals  = [1, 2];

    fig5 = figure('Name', 'Interaction Plot: Task x Session', 'Color', 'w', 'Position', [300 300 600 500]); hold on;
    errorbar(x_vals, int_means, int_sems, '-o', 'Color', col_int, 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'CapSize', 8);
    errorbar(x_vals, pho_means, pho_sems, '-o', 'Color', col_pho, 'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'CapSize', 8);
    xlim([0.5, 2.5]); xticks(x_vals); xticklabels({'Pretest', 'Scanner (VX)'});
    ylabel('Threshold / PSE (Stimulus Level)', 'FontSize', 14);
    title('Interaction Effect on Perceptual Threshold', 'FontSize', 16, 'FontWeight', 'bold');
    y_min_int = min([pho_means - pho_sems, int_means - int_sems]) - 0.2;
    y_max_int = max([pho_means + pho_sems, int_means + int_sems]) + 0.2;
    ylim([y_min_int, y_max_int]);
    set(gca, 'FontSize', 12, 'TickDir', 'out'); grid on; box on;
    legend({'Intonation', 'Phoneme'}, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');
    saveas(fig5, fullfile(STATS_OUTDIR, 'interaction_plot_PSE_lines.png'));
    close(fig5);
end

fprintf('\n=======================================================\n');
fprintf('=== PART 3F: ASSUMPTION CHECKS (NORMALITY)          ===\n');
fprintf('=======================================================\n');

if valid_subs > 2
    % 1. Check Normality for ANOVA variables (PSE)
    fprintf('--- Normality Check: PSE Variables (Kolmogorov-Smirnov) ---\n');
    [h_pp, p_pp] = kstest(zscore(data_PSE(:,1))); % Pretest Phoneme
    [h_pv, p_pv] = kstest(zscore(data_PSE(:,2))); % Scanner Phoneme
    [h_ip, p_ip] = kstest(zscore(data_PSE(:,3))); % Pretest Intonation
    [h_iv, p_iv] = kstest(zscore(data_PSE(:,4))); % Scanner Intonation
    
    fprintf('Pretest Phoneme    | p = %.4f %s\n', p_pp, check_normality_string(h_pp));
    fprintf('Scanner Phoneme    | p = %.4f %s\n', p_pv, check_normality_string(h_pv));
    fprintf('Pretest Intonation | p = %.4f %s\n', p_ip, check_normality_string(h_ip));
    fprintf('Scanner Intonation | p = %.4f %s\n\n', p_iv, check_normality_string(h_iv));
    
    % 2. Check Normality for Paired t-test differences (Scanner AIC)
    fprintf('--- Normality Check: Paired Differences (Scanner AIC) ---\n');
    diff_pho = aic_lin_pho_vx - aic_sig_pho_vx;
    diff_int = aic_lin_int_vx - aic_sig_int_vx;
    
    [h_diff_pho, p_diff_pho] = kstest(zscore(diff_pho));
    [h_diff_int, p_diff_int] = kstest(zscore(diff_int));
    
    fprintf('Diff Phoneme (Lin-Sig)    | p = %.4f %s\n', p_diff_pho, check_normality_string(h_diff_pho));
    fprintf('Diff Intonation (Lin-Sig) | p = %.4f %s\n', p_diff_int, check_normality_string(h_diff_int));
end

fprintf('\n=======================================================\n');
fprintf('=== PART 4: TMS ANALYSIS (VX vs. LT vs. RT)         ===\n');
fprintf('=======================================================\n');

% 1. Setup Colors & Conditions
cond_colors = struct(...
    'pre', [0.2 0.2 0.2], ...   % Black/Dark Grey for Pretest
    'vx', [0.60 0.60 0.60], ... % Light Grey
    'lt', [0.60 0.30 0.75], ... % Purple
    'rt', [0.95 0.70 0.10]  ... % Orange/Yellow
);
x_range = linspace(1, 5, 100);

tasks_plot = {'phoneme', 'intonation'};
titles = {'Phoneme', 'Intonation'};

for t = 1:2
    task_name = tasks_plot{t};
    fprintf('\n--- Analyzing Task: %s ---\n', upper(task_name));
    
    % Load Data
    load(fullfile(FIT_DIR, [task_name '_resp_fit_all_conditions.mat']), 'AllFitInfo');
    
    % Setup exclusion lists for this task
    incl_vx = setdiff(VP, [manual.(task_name).all, manual.(task_name).vx]);
    incl_lt = setdiff(VP, [manual.(task_name).all, manual.(task_name).lt]);
    incl_rt = setdiff(VP, [manual.(task_name).all, manual.(task_name).rt]);
    
    % --- Step 4A: Check Model Fits (AIC) for LT and RT ---
    [~, ~, ~, ~, ~, aic_lin_lt, aic_sig_lt] = calculate_scanner_metrics([AllFitInfo.lt.average], @(s) ismember(s, incl_lt), N_TRIALS_SCANNER);
    [~, ~, ~, ~, ~, aic_lin_rt, aic_sig_rt] = calculate_scanner_metrics([AllFitInfo.rt.average], @(s) ismember(s, incl_rt), N_TRIALS_SCANNER);
    
    [~, p_aic_lt] = ttest(aic_lin_lt, aic_sig_lt);
    [~, p_aic_rt] = ttest(aic_lin_rt, aic_sig_rt);
    
    fprintf('AIC Comparison LT: Linear M=%.2f vs Sigmoid M=%.2f | p = %.4f\n', mean(aic_lin_lt), mean(aic_sig_lt), p_aic_lt);
    fprintf('AIC Comparison RT: Linear M=%.2f vs Sigmoid M=%.2f | p = %.4f\n', mean(aic_lin_rt), mean(aic_sig_rt), p_aic_rt);

    % --- Step 4B: Extract PSEs and Slopes for ANOVA ---
    % Find intersection of subjects valid in ALL conditions for paired stats
    valid_subs_tms = intersect(incl_vx, intersect(incl_lt, incl_rt));
    pse_vx = []; pse_lt = []; pse_rt = [];
    slope_vx = []; slope_lt = []; slope_rt = [];
    
    for s = 1:length(valid_subs_tms)
        sub_id = valid_subs_tms{s};
        
        idx_v = find(strcmp({AllFitInfo.vx.average.subject}, sub_id) | strcmp({AllFitInfo.vx.average.subject}, ['sub-' sub_id]), 1);
        idx_l = find(strcmp({AllFitInfo.lt.average.subject}, sub_id) | strcmp({AllFitInfo.lt.average.subject}, ['sub-' sub_id]), 1);
        idx_r = find(strcmp({AllFitInfo.rt.average.subject}, sub_id) | strcmp({AllFitInfo.rt.average.subject}, ['sub-' sub_id]), 1);
        
        % Extract PSE (alpha)
        pse_vx(end+1) = AllFitInfo.vx.average(idx_v).mlfit.alpha;
        pse_lt(end+1) = AllFitInfo.lt.average(idx_l).mlfit.alpha;
        pse_rt(end+1) = AllFitInfo.rt.average(idx_r).mlfit.alpha;
        
        % Extract Slope (beta)
        slope_vx(end+1) = AllFitInfo.vx.average(idx_v).mlfit.beta;
        slope_lt(end+1) = AllFitInfo.lt.average(idx_l).mlfit.beta;
        slope_rt(end+1) = AllFitInfo.rt.average(idx_r).mlfit.beta;
    end
    
    % Run 1-Way Repeated Measures ANOVA on Scanner Conditions (PSE & SLOPE)
    if length(valid_subs_tms) > 2
        % Setup ANOVA variables
        Condition = categorical([1 2 3]', [1 2 3], {'VX', 'LT', 'RT'});
        
        % --- ANOVA for PSE ---
        T_PSE_TMS = array2table([pse_vx', pse_lt', pse_rt'], 'VariableNames', {'VX', 'LT', 'RT'});
        rm_pse_tms = fitrm(T_PSE_TMS, 'VX-RT ~ 1', 'WithinDesign', table(Condition));
        anova_pse_tms = ranova(rm_pse_tms, 'WithinModel', 'Condition');
        
        fprintf('\n--- ANOVA Results for TMS effect on PSE (%s) ---\n', task_name);
        disp(anova_pse_tms);
        % Index 3 is the row for '(Intercept):Condition'
        if anova_pse_tms.pValue(3) < 0.05
            fprintf('>> SIGNIFICANT TMS EFFECT ON PSE! (p = %.4f)\n', anova_pse_tms.pValue(3));
        else
            fprintf('>> No significant TMS effect on PSE (p = %.4f)\n', anova_pse_tms.pValue(3));
        end
        
        % --- ANOVA for SLOPE ---
        T_SLOPE_TMS = array2table([slope_vx', slope_lt', slope_rt'], 'VariableNames', {'VX', 'LT', 'RT'});
        rm_slope_tms = fitrm(T_SLOPE_TMS, 'VX-RT ~ 1', 'WithinDesign', table(Condition));
        anova_slope_tms = ranova(rm_slope_tms, 'WithinModel', 'Condition');
        
        fprintf('\n--- ANOVA Results for TMS effect on SLOPE (%s) ---\n', task_name);
        disp(anova_slope_tms);
        if anova_slope_tms.pValue(3) < 0.05
            fprintf('>> SIGNIFICANT TMS EFFECT ON SLOPE! (p = %.4f)\n', anova_slope_tms.pValue(3));
        else
            fprintf('>> No significant TMS effect on SLOPE (p = %.4f)\n', anova_slope_tms.pValue(3));
        end
    end

    % --- Step 4C: Generate Combined Average Plot ---
    fig_comb = figure('Color', 'w', 'Position', [400, 300, 600, 450]); hold on;
    title(sprintf('%s - Average Curves (Pre vs. TMS)', titles{t}), 'FontSize', 16, 'FontWeight', 'bold');
    
    conditions_to_plot = {'pre', 'vx', 'lt', 'rt'};
    plot_handles = [];
    
    for c = 1:length(conditions_to_plot)
        cond = conditions_to_plot{c};
        subj_y = [];
        
        if strcmp(cond, 'pre')
            % Pretest curve extraction
            for s = 1:length(valid_subs_tms)
                fpath = fullfile(MOCS_DIR, sprintf('%s_PRONET_MOCS_results.mat', valid_subs_tms{s}));
                if ~exist(fpath, 'file'), continue; end
                tmp = load(fpath);
                task_pre_name = task_name;
                if strcmp(task_name, 'intonation'), task_pre_name = 'intonat'; end
                p1 = tmp.m1.(task_pre_name).PFfit.params; 
                p2 = tmp.f2.(task_pre_name).PFfit.params;
                y_avg = (eval_logistic(p1, x_range) + eval_logistic(p2, x_range))/2;
                subj_y(end+1, :) = y_avg;
            end
        else
            % Scanner curve extraction
            for s = 1:length(valid_subs_tms)
                sub_id = valid_subs_tms{s};
                idx = find(strcmp({AllFitInfo.(cond).average.subject}, sub_id) | strcmp({AllFitInfo.(cond).average.subject}, ['sub-' sub_id]), 1);
                alpha = AllFitInfo.(cond).average(idx).mlfit.alpha;
                beta  = AllFitInfo.(cond).average(idx).mlfit.beta;
                subj_y(end+1, :) = eval_logistic([alpha, beta, 0, 0], x_range);
            end
        end
        
        % Calculate and plot mean curve
        if ~isempty(subj_y)
            mean_y = mean(subj_y, 1);
            if strcmp(cond, 'pre')
                p = plot(x_range, mean_y, '--', 'Color', cond_colors.(cond), 'LineWidth', 3);
            else
                p = plot(x_range, mean_y, '-', 'Color', cond_colors.(cond), 'LineWidth', 3);
            end
            plot_handles(end+1) = p;
        end
    end
    
    % Formatting
    if strcmp(task_name, 'intonation')
        ylabel('P("Question")', 'FontSize', 14); xlabel('Prosody Morph Levels', 'FontSize', 14);
    else
        ylabel('P("Paar")', 'FontSize', 14); xlabel('Phoneme Morph Levels', 'FontSize', 14);
    end
    
    ylim([0 1]); xlim([1 5]); xticks(1:5); yticks(0:0.2:1);
    set(gca, 'FontSize', 12, 'TickDir', 'out'); grid on; box on;
    legend(plot_handles, {'Pretest', 'Scanner (VX)', 'Scanner (LT)', 'Scanner (RT)'}, 'Location', 'southeast', 'FontSize', 12, 'Box', 'off');
    
    saveas(fig_comb, fullfile(STATS_OUTDIR, sprintf('combined_TMS_avg_curves_%s.png', task_name)));
    close(fig_comb);
end

fprintf('Master Pipeline completely finished! Outputs saved in:\n %s\n', FIT_DIR);

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function [subjects, r2_lin, r2_sig, adj_lin, adj_sig, aic_lin, aic_sig] = calculate_scanner_metrics(FitStruct, include_func, n_trials)
    all_subjects = unique({FitStruct.subject});
    keep = cellfun(include_func, all_subjects);
    subjects = all_subjects(keep);
    
    r2_lin = nan(size(subjects)); r2_sig = nan(size(subjects));
    adj_lin = nan(size(subjects)); adj_sig = nan(size(subjects));
    aic_lin = nan(size(subjects)); aic_sig = nan(size(subjects));

    for i = 1:numel(subjects)
        idx = strcmp({FitStruct.subject}, subjects{i});
        entry = FitStruct(idx);
        
        x = entry(1).x;
        y_true = entry(1).y; % Raw accuracy (proportions)
        k_true = round(y_true .* n_trials);
        n_array = repmat(n_trials, 1, length(x));

        % Linear Fit 
        p_lin = polyfit(x, y_true, 1);
        y_lin = polyval(p_lin, x);
        r2_lin(i) = calc_r2(y_true, y_lin);
        adj_lin(i) = 1 - ((1 - r2_lin(i)) * (5 - 1) / (5 - 2 - 1)); % n=5, k=2
        
        y_lin_clamped = max(min(y_lin, 0.999), 0.001);
        ll_lin = sum(k_true .* log(y_lin_clamped) + (n_array - k_true) .* log(1 - y_lin_clamped));
        aic_lin(i) = 4 - 2 * ll_lin; 

        % Sigmoid Fit 
        alpha = entry(1).mlfit.alpha;
        beta  = entry(1).mlfit.beta;
        gamma = entry(1).mlfit.gamma;
        lambda = entry(1).mlfit.lambda;
        
        y_sig = eval_logistic([alpha, beta, gamma, lambda], x);
        r2_sig(i) = calc_r2(y_true, y_sig);
        adj_sig(i) = 1 - ((1 - r2_sig(i)) * (5 - 1) / (5 - 2 - 1)); % n=5, k=2
        
        y_sig_clamped = max(min(y_sig, 0.999), 0.001);
        ll_sig = sum(k_true .* log(y_sig_clamped) + (n_array - k_true) .* log(1 - y_sig_clamped));
        aic_sig(i) = 4 - 2 * ll_sig; 
    end
end

function stars = get_pval_stars(pval)
    if pval < 0.001, stars = '***';
    elseif pval < 0.01, stars = '**';
    elseif pval < 0.05, stars = '*';
    else, stars = 'n.s.'; end
end

function r2 = calc_r2(y_true, y_pred)
    ss_res = sum((y_true - y_pred).^2);
    ss_tot = sum((y_true - mean(y_true)).^2);
    if ss_tot == 0
        if ss_res == 0, r2 = 1; else, r2 = 0; end
    else
        r2 = 1 - (ss_res / ss_tot);
    end
end

function y = eval_logistic(params, x)
    alpha = params(1); beta = params(2); 
    gamma = 0; lambda = 0;
    if length(params) >= 3, gamma = params(3); end
    if length(params) >= 4, lambda = params(4); end
    y = gamma + (1 - gamma - lambda) .* (1 ./ (1 + exp(-beta .* (x - alpha))));
end

function str = check_normality_string(h_val)
    % h = 0 means we CANNOT reject the null hypothesis (It IS normal)
    % h = 1 means we REJECT the null hypothesis (It is NOT normal)
    if h_val == 0
        str = '(Normal distribution assumed)';
    else
        str = '(WARNING: Significant deviation from normality)';
    end
end