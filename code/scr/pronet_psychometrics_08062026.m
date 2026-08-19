
%% Statistical Analysis - Linear vs. 4-Param Sigmoid Fit & RM-ANOVA (Pretest vs VX)
% Features manual exclusions, dynamic Adjusted R2 calculation for Pretest, 
% Adjusted R2 model fit bar plots, RM-ANOVAs, and high-resolution Psychometric Curves.
% Updated: June 2026 (Full Pipeline with Adjusted R2 & Individual Curve Plots)

clear; clc;

%% === 1. SETUP & PATHS ===
HPC_PATH  = '/mnt/beegfs/workspace/2024-0404-PRONET/';
DATA_PATH = fullfile(HPC_PATH, 'DATA', 'proc', 'psychometrics');
COMP_PATH = '/home/antonia.ceric/git/pronet/data/proc/psychometrics/fitting_comparison/';
OUTDIR    = fullfile(DATA_PATH, 'plots_06_2026'); 
if ~exist(OUTDIR, 'dir'), mkdir(OUTDIR); end

% Colors
col_pho = [0 0 0.8]; % Blue
col_int = [0.8 0 0]; % Red

% Subject List
VP = {'02','03','04','05','06','07','08','09','10','11','12','13','15','16','17',...
      '18','19','20','21','22','23','24','25','26','27','28','29','31','32'};

fprintf('\n=======================================================\n');
fprintf('=== PART 1: SCANNER DATA (VX CONDITION - ADJUSTED R2) ===\n');
fprintf('=======================================================\n');

% Load precomputed 4-parameter fit results
load(fullfile(COMP_PATH, 'intonation_resp_fit_all_conditions_comparison.mat'), 'AllFitInfo'); 
Fit_Intonation_VX = AllFitInfo.vx.separate;

load(fullfile(COMP_PATH, 'phoneme_resp_fit_all_conditions_comparison.mat'), 'AllFitInfo');  
Fit_Phoneme_VX = AllFitInfo.vx.separate;

% Load manual exclusions
load(fullfile(DATA_PATH, 'manual_exclusions.mat'), 'manual'); 

% Define included subject lists
incl_pho = setdiff(unique({Fit_Phoneme_VX.subject}), unique([manual.phoneme.all, manual.phoneme.vx]));
incl_int = setdiff(unique({Fit_Intonation_VX.subject}), unique([manual.intonation.all, manual.intonation.vx]));

% Extract Adjusted R² values for Phoneme
lin_pho_vx = []; sig_pho_vx = [];
for i = 1:length(incl_pho)
    idx = find(strcmp({Fit_Phoneme_VX.subject}, incl_pho{i}), 1);
    if ~isempty(idx)
        lin_pho_vx(end+1) = mean(Fit_Phoneme_VX(idx).linfit.adjrsq, 'omitnan');
        sig_pho_vx(end+1) = mean(Fit_Phoneme_VX(idx).sigfit4.adjrsq, 'omitnan');
    end
end

% Extract Adjusted R² values for Intonation
lin_int_vx = []; sig_int_vx = [];
for i = 1:length(incl_int)
    idx = find(strcmp({Fit_Intonation_VX.subject}, incl_int{i}), 1);
    if ~isempty(idx)
        lin_int_vx(end+1) = mean(Fit_Intonation_VX(idx).linfit.adjrsq, 'omitnan');
        sig_int_vx(end+1) = mean(Fit_Intonation_VX(idx).sigfit4.adjrsq, 'omitnan');
    end
end

% Paired t-tests
[~, p_pho_vx] = ttest(lin_pho_vx, sig_pho_vx);
[~, p_int_vx] = ttest(lin_int_vx, sig_int_vx);

fprintf('VX Phoneme    | Linear M=%.3f vs Sigmoid4 (Adj R²)=%.3f | p = %.4f\n', mean(lin_pho_vx), mean(sig_pho_vx), p_pho_vx);
fprintf('VX Intonation | Linear M=%.3f vs Sigmoid4 (Adj R²)=%.3f | p = %.4f\n', mean(lin_int_vx), mean(sig_int_vx), p_int_vx);

% --- PLOT 1: SCANNER (VX) ---
fig1 = figure('Name', 'Scanner VX: Lin vs Sig4', 'Color', 'w', 'Position', [200 200 700 600]); 
hold on;
x_pos = [1, 2, 4, 5]; w = 0.6;
bar(x_pos(1), mean(lin_pho_vx), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(2), mean(lin_int_vx), w, 'FaceColor', col_int, 'EdgeColor', 'none');
bar(x_pos(3), mean(sig_pho_vx), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(4), mean(sig_int_vx), w, 'FaceColor', col_int, 'EdgeColor', 'none');
sem = @(x) std(x)/sqrt(length(x));
errorbar(x_pos, [mean(lin_pho_vx), mean(lin_int_vx), mean(sig_pho_vx), mean(sig_int_vx)], ...
    [sem(lin_pho_vx), sem(lin_int_vx), sem(sig_pho_vx), sem(sig_int_vx)], 'k.', 'LineWidth', 1.5, 'CapSize', 0);
jitter = 0.15;
plot(x_pos(1) + randn(size(lin_pho_vx))*jitter, lin_pho_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(2) + randn(size(lin_int_vx))*jitter, lin_int_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
plot(x_pos(3) + randn(size(sig_pho_vx))*jitter, sig_pho_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(4) + randn(size(sig_int_vx))*jitter, sig_int_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
y_max = 1.05; plot([x_pos(1), x_pos(3)], [y_max y_max], 'k-', 'LineWidth', 1.5);
text(2.5, y_max + 0.03, get_pval_stars(p_pho_vx), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_pho);
y_max2 = 1.15; plot([x_pos(2), x_pos(4)], [y_max2 y_max2], 'k-', 'LineWidth', 1.5);
text(3.5, y_max2 + 0.03, get_pval_stars(p_int_vx), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_int);
title('Scanner (VX Condition) - Linear vs. 4-Param Sigmoid', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Model Fit (Adjusted R^2)', 'FontSize', 14);
ylim([0 1.25]); xlim([0 6]); xticks([1.5, 4.5]); xticklabels({'Linear Fit', 'Sigmoid Fit (4 Param)'});
set(gca, 'FontSize', 12, 'TickDir', 'out'); grid on; box on;
h1 = bar(NaN, NaN, 'FaceColor', col_pho, 'EdgeColor', 'none'); h2 = bar(NaN, NaN, 'FaceColor', col_int, 'EdgeColor', 'none');
legend([h1, h2], {'Phoneme', 'Intonation'}, 'Location', 'southeast', 'FontSize', 12, 'Box', 'off');
saveas(fig1, fullfile(OUTDIR, 'vx_linear_vs_sigmoid4_bar.png')); close(fig1);


fprintf('\n=======================================================\n');
fprintf('=== PART 2: PRETEST DATA (ADJUSTED R2 CONSISTENCY)  ===\n');
fprintf('=======================================================\n');

MOCS_DIR = fullfile(HPC_PATH, 'TMS/Pre-Session/results/');
excl_pre.phoneme    = {'08', '05'};
excl_pre.intonation = {'13', '05'}; 

lin_pho_pre = []; sig_pho_pre = [];
lin_int_pre = []; sig_int_pre = [];
x_levels = 1:5;

for i = 1:length(VP)
    sub_id = VP{i};
    fpath = fullfile(MOCS_DIR, sprintf('%s_PRONET_MOCS_results.mat', sub_id));
    if ~exist(fpath, 'file'), continue; end
    
    tmp = load(fpath);
    tasks_pre = {'phoneme', 'intonation'};
    file_names = {'phoneme', 'intonat'};
    
    for t = 1:2
        task = tasks_pre{t}; task_in = file_names{t};
        if ismember(sub_id, excl_pre.(task)), continue; end
        if ~isfield(tmp, 'm1') || ~isfield(tmp.m1, task_in), continue; end
        
        yes_tot = tmp.m1.(task_in).yes + tmp.f2.(task_in).yes;
        n_tot   = tmp.m1.(task_in).n + tmp.f2.(task_in).n;
        y_true  = yes_tot ./ n_tot;
        
        % Linear Fit (k=2)
        p_lin = polyfit(x_levels, y_true, 1);
        y_lin = polyval(p_lin, x_levels);
        r2_lin = calc_adjr2(y_true, y_lin, 2); 
        
        % Sigmoid 4-Param Fit (k=4)
        p_m1 = tmp.m1.(task_in).PFfit.params; p_f2 = tmp.f2.(task_in).PFfit.params;
        y_sig_m1 = eval_logistic(p_m1, x_levels); y_sig_f2 = eval_logistic(p_f2, x_levels);
        y_sig = (y_sig_m1 + y_sig_f2) / 2;
        r2_sig = calc_adjr2(y_true, y_sig, 4); 
        
        if strcmp(task, 'phoneme')
            lin_pho_pre(end+1) = r2_lin; sig_pho_pre(end+1) = r2_sig;
        else
            lin_int_pre(end+1) = r2_lin; sig_int_pre(end+1) = r2_sig;
        end
    end
end

[~, p_pho_pre] = ttest(lin_pho_pre, sig_pho_pre);
[~, p_int_pre] = ttest(lin_int_pre, sig_int_pre);

fprintf('Pretest Phoneme    | Linear M=%.3f vs Sigmoid4 (Adj R²)=%.3f | p = %.4f\n', mean(lin_pho_pre), mean(sig_pho_pre), p_pho_pre);
fprintf('Pretest Intonation | Linear M=%.3f vs Sigmoid4 (Adj R²)=%.3f | p = %.4f\n', mean(lin_int_pre), mean(sig_int_pre), p_int_pre);

% --- PLOT 2: PRETEST ---
fig2 = figure('Name', 'Pretest: Lin vs Sig4', 'Color', 'w', 'Position', [250 250 700 600]); 
hold on;
bar(x_pos(1), mean(lin_pho_pre), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(2), mean(lin_int_pre), w, 'FaceColor', col_int, 'EdgeColor', 'none');
bar(x_pos(3), mean(sig_pho_pre), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(4), mean(sig_int_pre), w, 'FaceColor', col_int, 'EdgeColor', 'none');
errorbar(x_pos, [mean(lin_pho_pre), mean(lin_int_pre), mean(sig_pho_pre), mean(sig_int_pre)], ...
    [sem(lin_pho_pre), sem(lin_int_pre), sem(sig_pho_pre), sem(sig_int_pre)], 'k.', 'LineWidth', 1.5, 'CapSize', 0);
plot(x_pos(1) + randn(size(lin_pho_pre))*jitter, lin_pho_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(2) + randn(size(lin_int_pre))*jitter, lin_int_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
plot(x_pos(3) + randn(size(sig_pho_pre))*jitter, sig_pho_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(4) + randn(size(sig_int_pre))*jitter, sig_int_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
plot([x_pos(1), x_pos(3)], [y_max y_max], 'k-', 'LineWidth', 1.5);
text(2.5, y_max + 0.03, get_pval_stars(p_pho_pre), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_pho);
plot([x_pos(2), x_pos(4)], [y_max2 y_max2], 'k-', 'LineWidth', 1.5);
text(3.5, y_max2 + 0.03, get_pval_stars(p_int_pre), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_int);
title('Pretest - Linear vs. 4-Param Sigmoid Fit', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Model Fit (Adjusted R^2)', 'FontSize', 14); 
ylim([0 1.25]); xlim([0 6]); xticks([1.5, 4.5]); xticklabels({'Linear Fit', 'Sigmoid Fit (4 Param)'});
set(gca, 'FontSize', 12, 'TickDir', 'out'); grid on; box on;
h1 = bar(NaN, NaN, 'FaceColor', col_pho, 'EdgeColor', 'none'); h2 = bar(NaN, NaN, 'FaceColor', col_int, 'EdgeColor', 'none');
legend([h1, h2], {'Phoneme', 'Intonation'}, 'Location', 'southeast', 'FontSize', 12, 'Box', 'off');
saveas(fig2, fullfile(OUTDIR, 'pretest_linear_vs_sigmoid_bar.png')); close(fig2);


fprintf('\n=======================================================\n');
fprintf('=== PART 3: 2x2 RM-ANOVA (TASK x SESSION)           ===\n');
fprintf('=======================================================\n');

% Setup Arrays for RM-ANOVA & High-Resolution Curves
data_PSE = []; data_Slope = []; valid_subs = 0;
x_fine = linspace(1, 5, 100); % Feines X-Raster für glatte Kurvenplots
curves_pho_pre = []; curves_pho_vx = [];
curves_int_pre = []; curves_int_vx = [];

for i = 1:length(VP)
    sub_id = VP{i};
    if ismember(sub_id, excl_pre.phoneme) || ismember(sub_id, excl_pre.intonation) || ...
       ismember(sub_id, manual.phoneme.all) || ismember(sub_id, manual.phoneme.vx) || ...
       ismember(sub_id, manual.intonation.all) || ismember(sub_id, manual.intonation.vx)
        continue; 
    end
    
    fpath = fullfile(MOCS_DIR, sprintf('%s_PRONET_MOCS_results.mat', sub_id));
    if ~exist(fpath, 'file'), continue; end
    tmp = load(fpath);
    if ~isfield(tmp, 'm1') || ~isfield(tmp.m1, 'phoneme') || ~isfield(tmp.m1, 'intonat'), continue; end
    
    % PRETEST Extract & Evaluate Curves
    p_pho_m1 = tmp.m1.phoneme.PFfit.params; p_pho_f2 = tmp.f2.phoneme.PFfit.params;
    pre_pho_pse   = (p_pho_m1(1) + p_pho_f2(1)) / 2;
    pre_pho_slope = (p_pho_m1(2) + p_pho_f2(2)) / 2;
    y_pho_pre_fine = (eval_logistic(p_pho_m1, x_fine) + eval_logistic(p_pho_f2, x_fine)) / 2;
    
    p_int_m1 = tmp.m1.intonat.PFfit.params; p_int_f2 = tmp.f2.intonat.PFfit.params;
    pre_int_pse   = (p_int_m1(1) + p_int_f2(1)) / 2;
    pre_int_slope = (p_int_m1(2) + p_int_f2(2)) / 2;
    y_int_pre_fine = (eval_logistic(p_int_m1, x_fine) + eval_logistic(p_int_f2, x_fine)) / 2;
    
    % VX Extract & Evaluate Curves
    idx_pho = find(strcmp({Fit_Phoneme_VX.subject}, sub_id), 1);
    idx_int = find(strcmp({Fit_Intonation_VX.subject}, sub_id), 1);
    if isempty(idx_pho) || isempty(idx_int), continue; end
    
    vx_pho_pse   = mean(Fit_Phoneme_VX(idx_pho).sigfit4.bias, 'omitnan');
    vx_pho_slope = mean(Fit_Phoneme_VX(idx_pho).sigfit4.slope, 'omitnan');
    p_pho_vx_v1 = [Fit_Phoneme_VX(idx_pho).sigfit4.bias(1), Fit_Phoneme_VX(idx_pho).sigfit4.slope(1), Fit_Phoneme_VX(idx_pho).sigfit4.guess_rate(1), Fit_Phoneme_VX(idx_pho).sigfit4.lapse_rate(1)];
    p_pho_vx_v2 = [Fit_Phoneme_VX(idx_pho).sigfit4.bias(2), Fit_Phoneme_VX(idx_pho).sigfit4.slope(2), Fit_Phoneme_VX(idx_pho).sigfit4.guess_rate(2), Fit_Phoneme_VX(idx_pho).sigfit4.lapse_rate(2)];
    y_pho_vx_fine = (eval_logistic(p_pho_vx_v1, x_fine) + eval_logistic(p_pho_vx_v2, x_fine)) / 2;
    
    vx_int_pse   = mean(Fit_Intonation_VX(idx_int).sigfit4.bias, 'omitnan');
    vx_int_slope = mean(Fit_Intonation_VX(idx_int).sigfit4.slope, 'omitnan');
    p_int_vx_v1 = [Fit_Intonation_VX(idx_int).sigfit4.bias(1), Fit_Intonation_VX(idx_int).sigfit4.slope(1), Fit_Intonation_VX(idx_int).sigfit4.guess_rate(1), Fit_Intonation_VX(idx_int).sigfit4.lapse_rate(1)];
    p_int_vx_v2 = [Fit_Intonation_VX(idx_int).sigfit4.bias(2), Fit_Intonation_VX(idx_int).sigfit4.slope(2), Fit_Intonation_VX(idx_int).sigfit4.guess_rate(2), Fit_Intonation_VX(idx_int).sigfit4.lapse_rate(2)];
    y_int_vx_fine = (eval_logistic(p_int_vx_v1, x_fine) + eval_logistic(p_int_vx_v2, x_fine)) / 2;
    
    % Collect all Matrix Data
    data_PSE   = [data_PSE; pre_pho_pse, vx_pho_pse, pre_int_pse, vx_int_pse];
    data_Slope = [data_Slope; pre_pho_slope, vx_pho_slope, pre_int_slope, vx_int_slope];
    
    curves_pho_pre = [curves_pho_pre; y_pho_pre_fine];
    curves_pho_vx  = [curves_pho_vx;  y_pho_vx_fine];
    curves_int_pre = [curves_int_pre; y_int_pre_fine];
    curves_int_vx  = [curves_int_vx;  y_int_vx_fine];
    valid_subs = valid_subs + 1;
end

fprintf('Included N = %d (Valid in all 4 conditions)\n\n', valid_subs);

if valid_subs > 2
    VarNames = {'Pho_Pre', 'Pho_VX', 'Int_Pre', 'Int_VX'};
    T_PSE   = array2table(data_PSE, 'VariableNames', VarNames);
    T_Slope = array2table(data_Slope, 'VariableNames', VarNames);
    Task    = categorical([1 1 2 2]', [1 2], {'Phoneme', 'Intonation'});
    Session = categorical([1 2 1 2]', [1 2], {'Pretest', 'VX'});
    WithinDesign = table(Task, Session);
    
    rm_PSE   = fitrm(T_PSE, 'Pho_Pre-Int_VX ~ 1', 'WithinDesign', WithinDesign);
    rm_Slope = fitrm(T_Slope, 'Pho_Pre-Int_VX ~ 1', 'WithinDesign', WithinDesign);
    anova_PSE   = ranova(rm_PSE, 'WithinModel', 'Task*Session');
    anova_Slope = ranova(rm_Slope, 'WithinModel', 'Task*Session');
    
    disp('--- RANOVA: PSE (Bias/Threshold) ---'); disp(anova_PSE);
    disp('--- RANOVA: SLOPE ---'); disp(anova_Slope);
else
    fprintf('Not enough valid subjects to run RM-ANOVA.\n');
end


%% === PART 4: POST-HOC TESTS ===
fprintf('\n=======================================================\n');
fprintf('=== POST-HOC TESTS: PSE INTERACTION                 ===\n');
fprintf('=======================================================\n');
pse_pho_pre = data_PSE(:,1); pse_pho_vx  = data_PSE(:,2);
pse_int_pre = data_PSE(:,3); pse_int_vx  = data_PSE(:,4);
[~, p_post_pho, ~, stats_pho] = ttest(pse_pho_pre, pse_pho_vx);
[~, p_post_int, ~, stats_int] = ttest(pse_int_pre, pse_int_vx);

fprintf('--- PHONEME ---\nPretest M = %.3f  |  Scanner (VX) M = %.3f\nMittlerer Shift = %+.3f  |  t(%d) = %.2f  |  p = %.4f\n\n', mean(pse_pho_pre), mean(pse_pho_vx), mean(pse_pho_vx - pse_pho_pre), stats_pho.df, stats_pho.tstat, p_post_pho);
fprintf('--- INTONATION ---\nPretest M = %.3f  |  Scanner (VX) M = %.3f\nMittlerer Shift = %+.3f  |  t(%d) = %.2f  |  p = %.4f\n\n', mean(pse_int_pre), mean(pse_int_vx), mean(pse_int_vx - pse_int_pre), stats_int.df, stats_int.tstat, p_post_int);

fprintf('\n=======================================================\n');
fprintf('=== POST-HOC TESTS: SLOPE INTERACTION               ===\n');
fprintf('=======================================================\n');
slope_pho_pre = data_Slope(:,1); slope_pho_vx  = data_Slope(:,2);
slope_int_pre = data_Slope(:,3); slope_int_vx  = data_Slope(:,4);
[~, p_post_slope_pho, ~, stats_slope_pho] = ttest(slope_pho_pre, slope_pho_vx);
[~, p_post_slope_int, ~, stats_slope_int] = ttest(slope_int_pre, slope_int_vx);

fprintf('--- PHONEME (SLOPE) ---\nPretest M = %.3f  |  Scanner (VX) M = %.3f\nMittlerer Shift = %+.3f  |  t(%d) = %.2f  |  p = %.4f\n\n', mean(slope_pho_pre), mean(slope_pho_vx), mean(slope_pho_vx - slope_pho_pre), stats_slope_pho.df, stats_slope_pho.tstat, p_post_slope_pho);
fprintf('--- INTONATION (SLOPE) ---\nPretest M = %.3f  |  Scanner (VX) M = %.3f\nMittlerer Shift = %+.3f  |  t(%d) = %.2f  |  p = %.4f\n\n', mean(slope_int_pre), mean(slope_int_vx), mean(slope_int_vx - slope_int_pre), stats_slope_int.df, stats_slope_int.tstat, p_post_slope_int);


%% === PART 5: VISUALIZE PSYCHOMETRIC CURVES (INDIVIDUAL PLOTS) ===
fprintf('\n=======================================================\n');
fprintf('=== PART 5: GENERATING INDIVIDUAL CURVE PLOTS       ===\n');
fprintf('=======================================================\n');

% Gemeinsame Plot-Einstellungen
lw_ind  = 1.2; % Liniendicke Einzelprobanden
lw_mean = 4;   % Liniendicke Gruppenmittelwert
col_ind = [0.65 0.65 0.65]; % Grau für Einzelprobanden

% --- 1. Plot: Phoneme Pretest ---
fig3 = figure('Name', 'Phoneme - Pretest', 'Color', 'w', 'Position', [100 100 500 450]);
hold on;
plot(x_fine, curves_pho_pre', 'Color', col_ind, 'LineStyle', ':', 'LineWidth', lw_ind);
plot(x_fine, mean(curves_pho_pre, 1), 'Color', col_pho, 'LineWidth', lw_mean);
title('Phoneme - Pretest', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('P("Paar")', 'FontSize', 12); 
xlabel('Phoneme Levels', 'FontSize', 12);
ylim([0 1]); xlim([1 5]); 
xticks(1:5); yticks(0:0.1:1);
grid on; box on; set(gca, 'TickDir', 'in', 'FontSize', 11);
saveas(fig3, fullfile(OUTDIR, 'curve_phoneme_pretest.png'));
close(fig3);

% --- 2. Plot: Phoneme Scanner (VX) ---
fig4 = figure('Name', 'Phoneme - VX', 'Color', 'w', 'Position', [150 150 500 450]);
hold on;
plot(x_fine, curves_pho_vx', 'Color', col_ind, 'LineStyle', ':', 'LineWidth', lw_ind);
plot(x_fine, mean(curves_pho_vx, 1), 'Color', col_pho, 'LineWidth', lw_mean);
title('Phoneme - VX', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('P(''Target'')', 'FontSize', 12); 
xlabel('Phoneme Levels', 'FontSize', 12);
ylim([0 1]); xlim([1 5]); 
xticks(1:5); yticks(0:0.1:1);
grid on; box on; set(gca, 'TickDir', 'in', 'FontSize', 11);
saveas(fig4, fullfile(OUTDIR, 'curve_phoneme_vx.png'));
close(fig4);

% --- 3. Plot: Intonation Pretest ---
fig5 = figure('Name', 'Intonation - Pretest', 'Color', 'w', 'Position', [200 200 500 450]);
hold on;
plot(x_fine, curves_int_pre', 'Color', col_ind, 'LineStyle', ':', 'LineWidth', lw_ind);
plot(x_fine, mean(curves_int_pre, 1), 'Color', col_int, 'LineWidth', lw_mean);
title('Intonation - Pretest', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('P(''Question'')', 'FontSize', 12); 
xlabel('Prosody Levels', 'FontSize', 12);
ylim([0 1]); xlim([1 5]); 
xticks(1:5); yticks(0:0.1:1);
grid on; box on; set(gca, 'TickDir', 'in', 'FontSize', 11);
saveas(fig5, fullfile(OUTDIR, 'curve_intonation_pretest.png'));
close(fig5);

% --- 4. Plot: Intonation Scanner (VX) ---
fig6 = figure('Name', 'Intonation - VX', 'Color', 'w', 'Position', [250 250 500 450]);
hold on;
plot(x_fine, curves_int_vx', 'Color', col_ind, 'LineStyle', ':', 'LineWidth', lw_ind);
plot(x_fine, mean(curves_int_vx, 1), 'Color', col_int, 'LineWidth', lw_mean);
title('Intonation - VX', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('P(''Question?'')', 'FontSize', 12); 
xlabel('Prosody Levels', 'FontSize', 12);
ylim([0 1]); xlim([1 5]); 
xticks(1:5); yticks(0:0.1:1);
grid on; box on; set(gca, 'TickDir', 'in', 'FontSize', 11);
saveas(fig6, fullfile(OUTDIR, 'curve_intonation_vx.png'));
close(fig6);

%% === PART 6: ASSUMPTION CHECKS (NORMALITY & OUTLIERS) ===
fprintf('\n=======================================================\n');
fprintf('=== PART 6: ASSUMPTION CHECKS (NORMALITY & OUTLIERS)===\n');
fprintf('=======================================================\n');

% Differenzen (Shifts) berechnen (bereits oben getan, hier nochmal sauber referenziert)
shift_pse_pho = pse_pho_vx - pse_pho_pre;
shift_pse_int = pse_int_vx - pse_int_pre;
shift_slope_pho = slope_pho_vx - slope_pho_pre;
shift_slope_int = slope_int_vx - slope_int_pre;

% Matrix aller Shifts für Schleife
shifts = [shift_pse_pho, shift_pse_int, shift_slope_pho, shift_slope_int];
shift_names = {'PSE Phoneme', 'PSE Intonation', 'Slope Phoneme', 'Slope Intonation'};

fprintf('1. NORMALITY OF DIFFERENCES (Shapiro-Wilk Test)\n');
fprintf('   H0: Data is normally distributed (p > 0.05 is good).\n');

for s = 1:4
    data = shifts(:, s);
    % Eigener, simpler Shapiro-Wilk Test Ersatz (Lilliefors Test), da SW nicht in Base-MATLAB ist
    [h, p_val] = lillietest(data); 
    
    if p_val > 0.05
        res_str = 'PASS (Normal)';
    else
        res_str = 'WARN (Non-Normal)';
    end
    fprintf('   %-18s | p = %.4f | %s\n', shift_names{s}, p_val, res_str);
end

fprintf('\n2. OUTLIER DETECTION (Z-Scores > |2.5|)\n');
fprintf('   Checks if individual subjects drive the effect excessively.\n');

for s = 1:4
    data = shifts(:, s);
    z_scores = (data - mean(data)) / std(data);
    outlier_idx = find(abs(z_scores) > 2.5);
    
    if isempty(outlier_idx)
        fprintf('   %-18s | PASS (0 Outliers found)\n', shift_names{s});
    else
        fprintf('   %-18s | WARN (%d Outlier(s) found: Subj-Indices ', shift_names{s}, length(outlier_idx));
        fprintf('%d ', outlier_idx);
        fprintf(')\n');
    end
end

fprintf('\n3. SPHERICITY ASSUMPTION\n');
fprintf('   PASS (Automatically met, as all Within-Subject factors only have 2 levels).\n\n');

fprintf('Alle Plots wurden erfolgreich gespeichert in:\n%s\n', OUTDIR);
fprintf('Pipeline Finished.\n');


%% === HELPER FUNCTIONS ===
function stars = get_pval_stars(pval)
    if pval < 0.001, stars = '***';
    elseif pval < 0.01, stars = '**';
    elseif pval < 0.05, stars = '*';
    else, stars = 'n.s.'; end
end

% Adjusted R2 Funktion 
function adjr2 = calc_adjr2(y_true, y_pred, k)
    n = length(y_true);
    ss_res = sum((y_true - y_pred).^2);
    ss_tot = sum((y_true - mean(y_true)).^2);
    
    if ss_tot == 0
        if ss_res == 0, r2 = 1; else, r2 = 0; end
    else
        r2 = 1 - (ss_res / ss_tot);
    end
    
    % Berechnung des Adjusted R2 basierend auf Parametern (k) und Datenpunkten (n)
    if n > k
        adjr2 = 1 - ((1 - r2) * (n - 1) / (n - k));
    else
        adjr2 = NaN; % Fallback falls n <= k
    end
end

function y = eval_logistic(params, x)
    alpha = params(1); beta = params(2); 
    gamma = 0; lambda = 0;
    if length(params) >= 3, gamma = params(3); end
    if length(params) >= 4, lambda = params(4); end
    y = gamma + (1 - gamma - lambda) .* (1 ./ (1 + exp(-beta .* (x - alpha))));
end

%% Statistical Analysis - Linear vs. 4-Param Sigmoid Fit & RM-ANOVA (Pretest vs VX)
% % Features manual exclusions, dynamic R2 calculation for Pretest, 
% % clean publication-style bar plots, and 2x2 RM-ANOVAs for Slope & PSE.
% % Updated: June 2026 (Using sigfit4 for VX condition)
% 
% clear; clc;
% 
% %% === 1. SETUP & PATHS ===
% HPC_PATH = '/mnt/beegfs/workspace/2024-0404-PRONET/';
% DATA_PATH = fullfile(HPC_PATH, 'DATA', 'proc', 'psychometrics');
% COMP_PATH = '/home/antonia.ceric/git/pronet/data/proc/psychometrics/fitting_comparison/';
% OUTDIR    = fullfile(DATA_PATH, 'psychometric_plots_08062026');
% if ~exist(OUTDIR, 'dir'), mkdir(OUTDIR); end
% 
% % Colors
% col_pho = [0 0 0.8]; % Blue
% col_int = [0.8 0 0]; % Red
% 
% % Subject List
% VP = {'02','03','04','05','06','07','08','09','10','11','12','13','15','16','17',...
%       '18','19','20','21','22','23','24','25','26','27','28','29','31','32'};
% 
% fprintf('\n=======================================================\n');
% fprintf('=== PART 1: SCANNER DATA (VX CONDITION - 4 PARAM)   ===\n');
% fprintf('=======================================================\n');
% 
% % Load precomputed 4-parameter fit results (comparison files)
% load(fullfile(COMP_PATH, 'intonation_resp_fit_all_conditions_comparison.mat'), 'AllFitInfo'); 
% Fit_Intonation_VX = AllFitInfo.vx.separate;
% 
% % Annahme: Die Phoneme-Datei liegt im selben Ordner und heißt analog
% load(fullfile(COMP_PATH, 'phoneme_resp_fit_all_conditions_comparison.mat'), 'AllFitInfo');  
% Fit_Phoneme_VX = AllFitInfo.vx.separate;
% 
% % Load manual exclusions
% load(fullfile(DATA_PATH, 'manual_exclusions.mat'), 'manual'); 
% 
% % Define included subject lists for VX specifically
% incl_pho = setdiff(unique({Fit_Phoneme_VX.subject}), unique([manual.phoneme.all, manual.phoneme.vx]));
% incl_int = setdiff(unique({Fit_Intonation_VX.subject}), unique([manual.intonation.all, manual.intonation.vx]));
% 
% % Extract paired Adjusted R² values for Phoneme (Averaging over voices m1 and f2)
% lin_pho_vx = []; sig_pho_vx = [];
% for i = 1:length(incl_pho)
%     idx = find(strcmp({Fit_Phoneme_VX.subject}, incl_pho{i}), 1);
%     if ~isempty(idx)
%         lin_pho_vx(end+1) = mean(Fit_Phoneme_VX(idx).linfit.adjrsq, 'omitnan');
%         sig_pho_vx(end+1) = mean(Fit_Phoneme_VX(idx).sigfit4.adjrsq, 'omitnan');
%     end
% end
% 
% % Extract paired Adjusted R² values for Intonation (Averaging over voices m1 and f2)
% lin_int_vx = []; sig_int_vx = [];
% for i = 1:length(incl_int)
%     idx = find(strcmp({Fit_Intonation_VX.subject}, incl_int{i}), 1);
%     if ~isempty(idx)
%         lin_int_vx(end+1) = mean(Fit_Intonation_VX(idx).linfit.adjrsq, 'omitnan');
%         sig_int_vx(end+1) = mean(Fit_Intonation_VX(idx).sigfit4.adjrsq, 'omitnan');
%     end
% end
% 
% % Paired t-tests
% [~, p_pho_vx] = ttest(lin_pho_vx, sig_pho_vx);
% [~, p_int_vx] = ttest(lin_int_vx, sig_int_vx);
% 
% fprintf('VX Phoneme    | Linear M=%.3f vs Sigmoid4 M=%.3f | p = %.4f\n', mean(lin_pho_vx), mean(sig_pho_vx), p_pho_vx);
% fprintf('VX Intonation | Linear M=%.3f vs Sigmoid4 M=%.3f | p = %.4f\n', mean(lin_int_vx), mean(sig_int_vx), p_int_vx);
% 
% % --- PLOT 1: SCANNER (VX) ---
% fig1 = figure('Name', 'Scanner VX: Lin vs Sig4', 'Color', 'w', 'Position', [200 200 700 600]); 
% hold on;
% 
% % Positions: Linear = Left (1, 2), Sigmoid = Right (4, 5)
% x_pos = [1, 2, 4, 5]; 
% w = 0.6;
% 
% % Bars
% bar(x_pos(1), mean(lin_pho_vx), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
% bar(x_pos(2), mean(lin_int_vx), w, 'FaceColor', col_int, 'EdgeColor', 'none');
% bar(x_pos(3), mean(sig_pho_vx), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
% bar(x_pos(4), mean(sig_int_vx), w, 'FaceColor', col_int, 'EdgeColor', 'none');
% 
% % Error Bars (SEM)
% sem = @(x) std(x)/sqrt(length(x));
% errorbar(x_pos, [mean(lin_pho_vx), mean(lin_int_vx), mean(sig_pho_vx), mean(sig_int_vx)], ...
%     [sem(lin_pho_vx), sem(lin_int_vx), sem(sig_pho_vx), sem(sig_int_vx)], ...
%     'k.', 'LineWidth', 1.5, 'CapSize', 0);
% 
% % Jittered Individual Points
% jitter = 0.15;
% plot(x_pos(1) + randn(size(lin_pho_vx))*jitter, lin_pho_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
% plot(x_pos(2) + randn(size(lin_int_vx))*jitter, lin_int_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
% plot(x_pos(3) + randn(size(sig_pho_vx))*jitter, sig_pho_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
% plot(x_pos(4) + randn(size(sig_int_vx))*jitter, sig_int_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
% 
% % Significance Brackets
% y_max = 1.05;
% plot([x_pos(1), x_pos(3)], [y_max y_max], 'k-', 'LineWidth', 1.5);
% text(2.5, y_max + 0.03, get_pval_stars(p_pho_vx), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_pho);
% 
% y_max2 = 1.15;
% plot([x_pos(2), x_pos(4)], [y_max2 y_max2], 'k-', 'LineWidth', 1.5);
% text(3.5, y_max2 + 0.03, get_pval_stars(p_int_vx), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_int);
% 
% % Formatting
% title('Scanner (VX Condition) - Linear vs. 4-Param Sigmoid', 'FontSize', 16, 'FontWeight', 'bold');
% ylabel('Model Fit (Adjusted R^2)', 'FontSize', 14);
% ylim([0 1.25]); xlim([0 6]);
% xticks([1.5, 4.5]); xticklabels({'Linear Fit', 'Sigmoid Fit (4 Param)'});
% set(gca, 'FontSize', 12, 'TickDir', 'out');
% grid on; box on;
% 
% % Dummy handles for legend
% h1 = bar(NaN, NaN, 'FaceColor', col_pho, 'EdgeColor', 'none'); 
% h2 = bar(NaN, NaN, 'FaceColor', col_int, 'EdgeColor', 'none');
% legend([h1, h2], {'Phoneme', 'Intonation'}, 'Location', 'southeast', 'FontSize', 12, 'Box', 'off');
% 
% saveas(fig1, fullfile(OUTDIR, 'vx_linear_vs_sigmoid4_bar.png'));
% close(fig1);
% 
% 
% fprintf('\n=======================================================\n');
% fprintf('=== PART 2: PRETEST DATA                            ===\n');
% fprintf('=======================================================\n');
% 
% MOCS_DIR = fullfile(HPC_PATH, 'TMS/Pre-Session/results/');
% 
% excl_pre.phoneme    = {'08', '05'};
% excl_pre.intonation = {'13', '05'}; 
% 
% lin_pho_pre = []; sig_pho_pre = [];
% lin_int_pre = []; sig_int_pre = [];
% x_levels = 1:5;
% 
% for i = 1:length(VP)
%     sub_id = VP{i};
%     fpath = fullfile(MOCS_DIR, sprintf('%s_PRONET_MOCS_results.mat', sub_id));
%     if ~exist(fpath, 'file'), continue; end
% 
%     tmp = load(fpath);
%     tasks_pre = {'phoneme', 'intonation'};
%     file_names = {'phoneme', 'intonat'};
% 
%     for t = 1:2
%         task = tasks_pre{t}; task_in = file_names{t};
%         if ismember(sub_id, excl_pre.(task)), continue; end
%         if ~isfield(tmp, 'm1') || ~isfield(tmp.m1, task_in), continue; end
% 
%         % Raw data pooling (Ground Truth)
%         yes_tot = tmp.m1.(task_in).yes + tmp.f2.(task_in).yes;
%         n_tot   = tmp.m1.(task_in).n + tmp.f2.(task_in).n;
%         y_true  = yes_tot ./ n_tot;
% 
%         % 1. Linear Fit Evaluation
%         p_lin = polyfit(x_levels, y_true, 1);
%         y_lin = polyval(p_lin, x_levels);
%         r2_lin = calc_r2(y_true, y_lin);
% 
%         % 2. Sigmoid Fit Evaluation (Avg of m1 and f2 params)
%         p_m1 = tmp.m1.(task_in).PFfit.params; p_f2 = tmp.f2.(task_in).PFfit.params;
%         y_sig_m1 = eval_logistic(p_m1, x_levels); y_sig_f2 = eval_logistic(p_f2, x_levels);
%         y_sig = (y_sig_m1 + y_sig_f2) / 2;
%         r2_sig = calc_r2(y_true, y_sig);
% 
%         if strcmp(task, 'phoneme')
%             lin_pho_pre(end+1) = r2_lin; sig_pho_pre(end+1) = r2_sig;
%         else
%             lin_int_pre(end+1) = r2_lin; sig_int_pre(end+1) = r2_sig;
%         end
%     end
% end
% 
% % Paired t-tests for Pretest
% [~, p_pho_pre] = ttest(lin_pho_pre, sig_pho_pre);
% [~, p_int_pre] = ttest(lin_int_pre, sig_int_pre);
% 
% fprintf('Pretest Phoneme    | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n', mean(lin_pho_pre), mean(sig_pho_pre), p_pho_pre);
% fprintf('Pretest Intonation | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n', mean(lin_int_pre), mean(sig_int_pre), p_int_pre);
% 
% % --- PLOT 2: PRETEST ---
% fig2 = figure('Name', 'Pretest: Lin vs Sig', 'Color', 'w', 'Position', [250 250 700 600]); 
% hold on;
% 
% bar(x_pos(1), mean(lin_pho_pre), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
% bar(x_pos(2), mean(lin_int_pre), w, 'FaceColor', col_int, 'EdgeColor', 'none');
% bar(x_pos(3), mean(sig_pho_pre), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
% bar(x_pos(4), mean(sig_int_pre), w, 'FaceColor', col_int, 'EdgeColor', 'none');
% 
% errorbar(x_pos, [mean(lin_pho_pre), mean(lin_int_pre), mean(sig_pho_pre), mean(sig_int_pre)], ...
%     [sem(lin_pho_pre), sem(lin_int_pre), sem(sig_pho_pre), sem(sig_int_pre)], ...
%     'k.', 'LineWidth', 1.5, 'CapSize', 0);
% 
% plot(x_pos(1) + randn(size(lin_pho_pre))*jitter, lin_pho_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
% plot(x_pos(2) + randn(size(lin_int_pre))*jitter, lin_int_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
% plot(x_pos(3) + randn(size(sig_pho_pre))*jitter, sig_pho_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
% plot(x_pos(4) + randn(size(sig_int_pre))*jitter, sig_int_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
% 
% % Significance Brackets
% plot([x_pos(1), x_pos(3)], [y_max y_max], 'k-', 'LineWidth', 1.5);
% text(2.5, y_max + 0.03, get_pval_stars(p_pho_pre), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_pho);
% 
% plot([x_pos(2), x_pos(4)], [y_max2 y_max2], 'k-', 'LineWidth', 1.5);
% text(3.5, y_max2 + 0.03, get_pval_stars(p_int_pre), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_int);
% 
% % Formatting
% title('Pretest - Linear vs. 4-Param Sigmoid Fit', 'FontSize', 16, 'FontWeight', 'bold');
% ylabel('Model Fit (R^2)', 'FontSize', 14);
% ylim([0 1.25]); xlim([0 6]);
% xticks([1.5, 4.5]); xticklabels({'Linear Fit', 'Sigmoid Fit'});
% set(gca, 'FontSize', 12, 'TickDir', 'out');
% grid on; box on;
% 
% % Dummy handles for legend
% h1 = bar(NaN, NaN, 'FaceColor', col_pho, 'EdgeColor', 'none'); 
% h2 = bar(NaN, NaN, 'FaceColor', col_int, 'EdgeColor', 'none');
% legend([h1, h2], {'Phoneme', 'Intonation'}, 'Location', 'southeast', 'FontSize', 12, 'Box', 'off');
% 
% saveas(fig2, fullfile(OUTDIR, 'pretest_linear_vs_sigmoid_bar.png'));
% close(fig2);
% fprintf('Plots generated and saved in: %s\n', OUTDIR);
% 
% 
% fprintf('\n=======================================================\n');
% fprintf('=== PART 3: 2x2 RM-ANOVA (TASK x SESSION)           ===\n');
% fprintf('=======================================================\n');
% 
% % Arrays to hold valid data
% data_PSE = [];
% data_Slope = [];
% valid_subs = 0;
% 
% for i = 1:length(VP)
%     sub_id = VP{i};
% 
%     % 1. Check ALL exclusions (Pretest & VX for both tasks)
%     if ismember(sub_id, excl_pre.phoneme) || ismember(sub_id, excl_pre.intonation) || ...
%        ismember(sub_id, manual.phoneme.all) || ismember(sub_id, manual.phoneme.vx) || ...
%        ismember(sub_id, manual.intonation.all) || ismember(sub_id, manual.intonation.vx)
%         continue; % Skip if excluded in ANY condition
%     end
% 
%     % 2. Extract Pretest Data (Average of m1 and f2 parameters)
%     fpath = fullfile(MOCS_DIR, sprintf('%s_PRONET_MOCS_results.mat', sub_id));
%     if ~exist(fpath, 'file'), continue; end
%     tmp = load(fpath);
%     if ~isfield(tmp, 'm1') || ~isfield(tmp.m1, 'phoneme') || ~isfield(tmp.m1, 'intonat'), continue; end
% 
%     % PSE = Parameter(1) | Slope = Parameter(2)
%     pre_pho_pse   = (tmp.m1.phoneme.PFfit.params(1) + tmp.f2.phoneme.PFfit.params(1)) / 2;
%     pre_pho_slope = (tmp.m1.phoneme.PFfit.params(2) + tmp.f2.phoneme.PFfit.params(2)) / 2;
%     pre_int_pse   = (tmp.m1.intonat.PFfit.params(1) + tmp.f2.intonat.PFfit.params(1)) / 2;
%     pre_int_slope = (tmp.m1.intonat.PFfit.params(2) + tmp.f2.intonat.PFfit.params(2)) / 2;
% 
%     % 3. Extract VX Data (Using sigfit4 and averaging across voices)
%     idx_pho = find(strcmp({Fit_Phoneme_VX.subject}, sub_id), 1);
%     idx_int = find(strcmp({Fit_Intonation_VX.subject}, sub_id), 1);
% 
%     if isempty(idx_pho) || isempty(idx_int), continue; end
% 
%     % Extracting bias and slope from sigfit4 and averaging across voices
%     vx_pho_pse   = mean(Fit_Phoneme_VX(idx_pho).sigfit4.bias, 'omitnan');
%     vx_pho_slope = mean(Fit_Phoneme_VX(idx_pho).sigfit4.slope, 'omitnan');
%     vx_int_pse   = mean(Fit_Intonation_VX(idx_int).sigfit4.bias, 'omitnan');
%     vx_int_slope = mean(Fit_Intonation_VX(idx_int).sigfit4.slope, 'omitnan');
% 
%     % Store valid subject data
%     data_PSE   = [data_PSE; pre_pho_pse, vx_pho_pse, pre_int_pse, vx_int_pse];
%     data_Slope = [data_Slope; pre_pho_slope, vx_pho_slope, pre_int_slope, vx_int_slope];
%     valid_subs = valid_subs + 1;
% end
% 
% fprintf('Included N = %d (Valid in all 4 conditions)\n\n', valid_subs);
% 
% % 4. Create Tables and RM-ANOVA Design
% if valid_subs > 2
%     % Format: Pho_Pre | Pho_VX | Int_Pre | Int_VX
%     VarNames = {'Pho_Pre', 'Pho_VX', 'Int_Pre', 'Int_VX'};
%     T_PSE   = array2table(data_PSE, 'VariableNames', VarNames);
%     T_Slope = array2table(data_Slope, 'VariableNames', VarNames);
% 
%     % Define Within-Subject Design
%     Task    = categorical([1 1 2 2]', [1 2], {'Phoneme', 'Intonation'});
%     Session = categorical([1 2 1 2]', [1 2], {'Pretest', 'VX'});
%     WithinDesign = table(Task, Session);
% 
%     % Fit Models
%     rm_PSE   = fitrm(T_PSE, 'Pho_Pre-Int_VX ~ 1', 'WithinDesign', WithinDesign);
%     rm_Slope = fitrm(T_Slope, 'Pho_Pre-Int_VX ~ 1', 'WithinDesign', WithinDesign);
% 
%     % Run ANOVAs
%     anova_PSE   = ranova(rm_PSE, 'WithinModel', 'Task*Session');
%     anova_Slope = ranova(rm_Slope, 'WithinModel', 'Task*Session');
% 
%     disp('--- RANOVA: PSE (Bias/Threshold) ---');
%     disp(anova_PSE);
% 
%     disp('--- RANOVA: SLOPE ---');
%     disp(anova_Slope);
% else
%     fprintf('Not enough valid subjects to run RM-ANOVA.\n');
% end
% 
% % --- 6. POST-HOC TESTS (PSE INTERACTION) ---
% fprintf('\n=======================================================\n');
% fprintf('=== POST-HOC TESTS: PSE INTERACTION                 ===\n');
% fprintf('=======================================================\n');
% 
% % Daten extrahieren (Spalten aus data_PSE: Pho_Pre, Pho_VX, Int_Pre, Int_VX)
% pse_pho_pre = data_PSE(:,1);
% pse_pho_vx  = data_PSE(:,2);
% pse_int_pre = data_PSE(:,3);
% pse_int_vx  = data_PSE(:,4);
% 
% % Verschiebungen (Shift) berechnen: Positiv = Verschiebung nach rechts, Negativ = nach links
% shift_pho = pse_pho_vx - pse_pho_pre;
% shift_int = pse_int_vx - pse_int_pre;
% 
% % Gepaarte t-Tests
% [~, p_post_pho, ~, stats_pho] = ttest(pse_pho_pre, pse_pho_vx);
% [~, p_post_int, ~, stats_int] = ttest(pse_int_pre, pse_int_vx);
% 
% % Ausgabe Phoneme
% fprintf('--- PHONEME ---\n');
% fprintf('Pretest M = %.3f  |  Scanner (VX) M = %.3f\n', mean(pse_pho_pre), mean(pse_pho_vx));
% fprintf('Mittlerer Shift = %+.3f  |  t(%d) = %.2f  |  p = %.4f\n\n', mean(shift_pho), stats_pho.df, stats_pho.tstat, p_post_pho);
% 
% % Ausgabe Intonation
% fprintf('--- INTONATION ---\n');
% fprintf('Pretest M = %.3f  |  Scanner (VX) M = %.3f\n', mean(pse_int_pre), mean(pse_int_vx));
% fprintf('Mittlerer Shift = %+.3f  |  t(%d) = %.2f  |  p = %.4f\n\n', mean(shift_int), stats_int.df, stats_int.tstat, p_post_int);
% 
% % --- 7. POST-HOC TESTS (SLOPE INTERACTION) ---
% fprintf('\n=======================================================\n');
% fprintf('=== POST-HOC TESTS: SLOPE INTERACTION               ===\n');
% fprintf('=======================================================\n');
% 
% % Daten extrahieren (Spalten aus data_Slope: Pho_Pre, Pho_VX, Int_Pre, Int_VX)
% slope_pho_pre = data_Slope(:,1);
% slope_pho_vx  = data_Slope(:,2);
% slope_int_pre = data_Slope(:,3);
% slope_int_vx  = data_Slope(:,4);
% 
% % Verschiebungen (Shift) berechnen: Positiv = steilere Kurve, Negativ = flachere Kurve
% shift_slope_pho = slope_pho_vx - slope_pho_pre;
% shift_slope_int = slope_int_vx - slope_int_pre;
% 
% % Gepaarte t-Tests
% [~, p_post_slope_pho, ~, stats_slope_pho] = ttest(slope_pho_pre, slope_pho_vx);
% [~, p_post_slope_int, ~, stats_slope_int] = ttest(slope_int_pre, slope_int_vx);
% 
% % Ausgabe Phoneme
% fprintf('--- PHONEME (SLOPE) ---\n');
% fprintf('Pretest M = %.3f  |  Scanner (VX) M = %.3f\n', mean(slope_pho_pre), mean(slope_pho_vx));
% fprintf('Mittlerer Shift = %+.3f  |  t(%d) = %.2f  |  p = %.4f\n\n', mean(shift_slope_pho), stats_slope_pho.df, stats_slope_pho.tstat, p_post_slope_pho);
% 
% % Ausgabe Intonation
% fprintf('--- INTONATION (SLOPE) ---\n');
% fprintf('Pretest M = %.3f  |  Scanner (VX) M = %.3f\n', mean(slope_int_pre), mean(slope_int_vx));
% fprintf('Mittlerer Shift = %+.3f  |  t(%d) = %.2f  |  p = %.4f\n\n', mean(shift_slope_int), stats_slope_int.df, stats_slope_int.tstat, p_post_slope_int);
% 
% fprintf('Pipeline Finished.\n');
% 
% %% === HELPER FUNCTIONS ===
% 
% function stars = get_pval_stars(pval)
%     if pval < 0.001, stars = '***';
%     elseif pval < 0.01, stars = '**';
%     elseif pval < 0.05, stars = '*';
%     else, stars = 'n.s.'; end
% end
% 
% function r2 = calc_r2(y_true, y_pred)
%     ss_res = sum((y_true - y_pred).^2);
%     ss_tot = sum((y_true - mean(y_true)).^2);
%     if ss_tot == 0
%         if ss_res == 0, r2 = 1; else, r2 = 0; end
%     else
%         r2 = 1 - (ss_res / ss_tot);
%     end
% end
% 
% function y = eval_logistic(params, x)
%     alpha = params(1); beta = params(2); 
%     gamma = 0; lambda = 0;
%     if length(params) >= 3, gamma = params(3); end
%     if length(params) >= 4, lambda = params(4); end
%     y = gamma + (1 - gamma - lambda) .* (1 ./ (1 + exp(-beta .* (x - alpha))));
% end