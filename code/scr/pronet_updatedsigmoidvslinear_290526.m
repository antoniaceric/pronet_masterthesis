%% Statistical Analysis - Linear vs. Sigmoid Fit & RM-ANOVA (Pretest vs VX)
% Date: May 2026

clear; clc;

%% === 1. SETUP & PATHS ===
HPC_PATH = '/mnt/beegfs/workspace/2024-0404-PRONET/';
DATA_PATH = fullfile(HPC_PATH, 'DATA', 'proc', 'psychometrics');
OUTDIR    = fullfile(DATA_PATH, 'plots_may_2026');
if ~exist(OUTDIR, 'dir'), mkdir(OUTDIR); end

% Colors
col_pho = [0 0 0.8]; % Blue
col_int = [0.8 0 0]; % Red

% Subject List
VP = {'02','03','04','05','06','07','08','09','10','11','12','13','15','16','17',...
      '18','19','20','21','22','23','24','25','26','27','28','29','31','32'};

fprintf('\n=======================================================\n');
fprintf('=== PART 1: SCANNER DATA (VX CONDITION)             ===\n');
fprintf('=======================================================\n');

% Load precomputed fit results
load(fullfile(DATA_PATH, 'intonation_resp_fit_all_conditions.mat'), 'AllFitInfo'); 
AllFitInfo_intonation = AllFitInfo;
Fit_Intonation = [AllFitInfo.vx.average];

load(fullfile(DATA_PATH, 'phoneme_resp_fit_all_conditions.mat'), 'AllFitInfo');  
AllFitInfo_phoneme = AllFitInfo;
Fit_Phoneme = [AllFitInfo.vx.average];

% Load manual exclusions
load(fullfile(DATA_PATH, 'manual_exclusions.mat'), 'manual'); 

% Define included subject lists for VX specifically
incl_pho = setdiff(unique({Fit_Phoneme.subject}), unique([manual.phoneme.all, manual.phoneme.vx]));
incl_int = setdiff(unique({Fit_Intonation.subject}), unique([manual.intonation.all, manual.intonation.vx]));

% Extract paired R² values
[~, lin_pho_vx, sig_pho_vx] = extract_adjrsq_clean(Fit_Phoneme, @(s) ismember(s, incl_pho));
[~, lin_int_vx, sig_int_vx] = extract_adjrsq_clean(Fit_Intonation, @(s) ismember(s, incl_int));

% Paired t-tests
[~, p_pho_vx] = ttest(lin_pho_vx, sig_pho_vx);
[~, p_int_vx] = ttest(lin_int_vx, sig_int_vx);

fprintf('VX Phoneme    | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n', mean(lin_pho_vx), mean(sig_pho_vx), p_pho_vx);
fprintf('VX Intonation | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n', mean(lin_int_vx), mean(sig_int_vx), p_int_vx);

% --- PLOT 1: SCANNER (VX) ---
fig1 = figure('Name', 'Scanner VX: Lin vs Sig', 'Color', 'w', 'Position', [200 200 700 600]); 
hold on;

% Positions: Linear = Left (1, 2), Sigmoid = Right (4, 5)
x_pos = [1, 2, 4, 5]; 
w = 0.6;

% Bars
bar(x_pos(1), mean(lin_pho_vx), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(2), mean(lin_int_vx), w, 'FaceColor', col_int, 'EdgeColor', 'none');
bar(x_pos(3), mean(sig_pho_vx), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(4), mean(sig_int_vx), w, 'FaceColor', col_int, 'EdgeColor', 'none');

% Error Bars (SEM)
sem = @(x) std(x)/sqrt(length(x));
errorbar(x_pos, [mean(lin_pho_vx), mean(lin_int_vx), mean(sig_pho_vx), mean(sig_int_vx)], ...
    [sem(lin_pho_vx), sem(lin_int_vx), sem(sig_pho_vx), sem(sig_int_vx)], ...
    'k.', 'LineWidth', 1.5, 'CapSize', 0);

% Jittered Individual Points
jitter = 0.15;
plot(x_pos(1) + randn(size(lin_pho_vx))*jitter, lin_pho_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(2) + randn(size(lin_int_vx))*jitter, lin_int_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
plot(x_pos(3) + randn(size(sig_pho_vx))*jitter, sig_pho_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(4) + randn(size(sig_int_vx))*jitter, sig_int_vx, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);

% Significance Brackets
y_max = 1.05;
plot([x_pos(1), x_pos(3)], [y_max y_max], 'k-', 'LineWidth', 1.5);
text(2.5, y_max + 0.03, get_pval_stars(p_pho_vx), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_pho);

y_max2 = 1.15;
plot([x_pos(2), x_pos(4)], [y_max2 y_max2], 'k-', 'LineWidth', 1.5);
text(3.5, y_max2 + 0.03, get_pval_stars(p_int_vx), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_int);

% Formatting
title('Scanner (VX Condition) - Linear vs. Sigmoid Fit', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Model Fit (Adjusted R^2)', 'FontSize', 14);
ylim([0 1.25]); xlim([0 6]);
xticks([1.5, 4.5]); xticklabels({'Linear Fit', 'Sigmoid Fit'});
set(gca, 'FontSize', 12, 'TickDir', 'out');
grid on; box on;

% Dummy handles for legend
h1 = bar(NaN, NaN, 'FaceColor', col_pho, 'EdgeColor', 'none'); 
h2 = bar(NaN, NaN, 'FaceColor', col_int, 'EdgeColor', 'none');
legend([h1, h2], {'Phoneme', 'Intonation'}, 'Location', 'southeast', 'FontSize', 12, 'Box', 'off');

saveas(fig1, fullfile(OUTDIR, 'vx_linear_vs_sigmoid_bar.png'));
close(fig1);


fprintf('\n=======================================================\n');
fprintf('=== PART 2: PRETEST DATA                            ===\n');
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
        
        % Raw data pooling (Ground Truth)
        yes_tot = tmp.m1.(task_in).yes + tmp.f2.(task_in).yes;
        n_tot   = tmp.m1.(task_in).n + tmp.f2.(task_in).n;
        y_true  = yes_tot ./ n_tot;
        
        % 1. Linear Fit Evaluation
        p_lin = polyfit(x_levels, y_true, 1);
        y_lin = polyval(p_lin, x_levels);
        r2_lin = calc_r2(y_true, y_lin);
        
        % 2. Sigmoid Fit Evaluation (Avg of m1 and f2 params)
        p_m1 = tmp.m1.(task_in).PFfit.params; p_f2 = tmp.f2.(task_in).PFfit.params;
        y_sig_m1 = eval_logistic(p_m1, x_levels); y_sig_f2 = eval_logistic(p_f2, x_levels);
        y_sig = (y_sig_m1 + y_sig_f2) / 2;
        r2_sig = calc_r2(y_true, y_sig);
        
        if strcmp(task, 'phoneme')
            lin_pho_pre(end+1) = r2_lin; sig_pho_pre(end+1) = r2_sig;
        else
            lin_int_pre(end+1) = r2_lin; sig_int_pre(end+1) = r2_sig;
        end
    end
end

% Paired t-tests for Pretest
[~, p_pho_pre] = ttest(lin_pho_pre, sig_pho_pre);
[~, p_int_pre] = ttest(lin_int_pre, sig_int_pre);

fprintf('Pretest Phoneme    | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n', mean(lin_pho_pre), mean(sig_pho_pre), p_pho_pre);
fprintf('Pretest Intonation | Linear M=%.3f vs Sigmoid M=%.3f | p = %.4f\n', mean(lin_int_pre), mean(sig_int_pre), p_int_pre);

% --- PLOT 2: PRETEST ---
fig2 = figure('Name', 'Pretest: Lin vs Sig', 'Color', 'w', 'Position', [250 250 700 600]); 
hold on;

bar(x_pos(1), mean(lin_pho_pre), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(2), mean(lin_int_pre), w, 'FaceColor', col_int, 'EdgeColor', 'none');
bar(x_pos(3), mean(sig_pho_pre), w, 'FaceColor', col_pho, 'EdgeColor', 'none');
bar(x_pos(4), mean(sig_int_pre), w, 'FaceColor', col_int, 'EdgeColor', 'none');

errorbar(x_pos, [mean(lin_pho_pre), mean(lin_int_pre), mean(sig_pho_pre), mean(sig_int_pre)], ...
    [sem(lin_pho_pre), sem(lin_int_pre), sem(sig_pho_pre), sem(sig_int_pre)], ...
    'k.', 'LineWidth', 1.5, 'CapSize', 0);

plot(x_pos(1) + randn(size(lin_pho_pre))*jitter, lin_pho_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(2) + randn(size(lin_int_pre))*jitter, lin_int_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);
plot(x_pos(3) + randn(size(sig_pho_pre))*jitter, sig_pho_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.6 0.6 1], 'MarkerSize', 5);
plot(x_pos(4) + randn(size(sig_int_pre))*jitter, sig_int_pre, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.6 0.6], 'MarkerSize', 5);

% Significance Brackets
plot([x_pos(1), x_pos(3)], [y_max y_max], 'k-', 'LineWidth', 1.5);
text(2.5, y_max + 0.03, get_pval_stars(p_pho_pre), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_pho);

plot([x_pos(2), x_pos(4)], [y_max2 y_max2], 'k-', 'LineWidth', 1.5);
text(3.5, y_max2 + 0.03, get_pval_stars(p_int_pre), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', col_int);

% Formatting
title('Pretest - Linear vs. Sigmoid Fit', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Model Fit (R^2)', 'FontSize', 14);
ylim([0 1.25]); xlim([0 6]);
xticks([1.5, 4.5]); xticklabels({'Linear Fit', 'Sigmoid Fit'});
set(gca, 'FontSize', 12, 'TickDir', 'out');
grid on; box on;

% Dummy handles for legend (recreated)
h1 = bar(NaN, NaN, 'FaceColor', col_pho, 'EdgeColor', 'none'); 
h2 = bar(NaN, NaN, 'FaceColor', col_int, 'EdgeColor', 'none');
legend([h1, h2], {'Phoneme', 'Intonation'}, 'Location', 'southeast', 'FontSize', 12, 'Box', 'off');

saveas(fig2, fullfile(OUTDIR, 'pretest_linear_vs_sigmoid_bar.png'));
close(fig2);
fprintf('Plots generated and saved in: %s\n', OUTDIR);


fprintf('\n=======================================================\n');
fprintf('=== PART 3: 2x2 RM-ANOVA (TASK x SESSION)           ===\n');
fprintf('=======================================================\n');

% Arrays to hold valid data
data_PSE = [];
data_Slope = [];
valid_subs = 0;

for i = 1:length(VP)
    sub_id = VP{i};
    
    % 1. Check ALL exclusions (Pretest & VX for both tasks)
    if ismember(sub_id, excl_pre.phoneme) || ismember(sub_id, excl_pre.intonation) || ...
       ismember(sub_id, manual.phoneme.all) || ismember(sub_id, manual.phoneme.vx) || ...
       ismember(sub_id, manual.intonation.all) || ismember(sub_id, manual.intonation.vx)
        continue; % Skip if excluded in ANY condition
    end
    
    % 2. Extract Pretest Data (Average of m1 and f2 parameters)
    fpath = fullfile(MOCS_DIR, sprintf('%s_PRONET_MOCS_results.mat', sub_id));
    if ~exist(fpath, 'file'), continue; end
    tmp = load(fpath);
    if ~isfield(tmp, 'm1') || ~isfield(tmp.m1, 'phoneme') || ~isfield(tmp.m1, 'intonat'), continue; end
    
    % PSE = Parameter(1) | Slope = Parameter(2)
    pre_pho_pse   = (tmp.m1.phoneme.PFfit.params(1) + tmp.f2.phoneme.PFfit.params(1)) / 2;
    pre_pho_slope = (tmp.m1.phoneme.PFfit.params(2) + tmp.f2.phoneme.PFfit.params(2)) / 2;
    pre_int_pse   = (tmp.m1.intonat.PFfit.params(1) + tmp.f2.intonat.PFfit.params(1)) / 2;
    pre_int_slope = (tmp.m1.intonat.PFfit.params(2) + tmp.f2.intonat.PFfit.params(2)) / 2;
    
    % 3. Extract VX Data
    idx_pho = find(strcmp({AllFitInfo_phoneme.vx.average.subject}, sub_id) | strcmp({AllFitInfo_phoneme.vx.average.subject}, ['sub-' sub_id]), 1);
    idx_int = find(strcmp({AllFitInfo_intonation.vx.average.subject}, sub_id) | strcmp({AllFitInfo_intonation.vx.average.subject}, ['sub-' sub_id]), 1);
    
    if isempty(idx_pho) || isempty(idx_int), continue; end
    
    vx_pho_pse   = AllFitInfo_phoneme.vx.average(idx_pho).sigfit.bias;
    vx_pho_slope = AllFitInfo_phoneme.vx.average(idx_pho).sigfit.slope;
    vx_int_pse   = AllFitInfo_intonation.vx.average(idx_int).sigfit.bias;
    vx_int_slope = AllFitInfo_intonation.vx.average(idx_int).sigfit.slope;
    
    % Store valid subject data
    data_PSE   = [data_PSE; pre_pho_pse, vx_pho_pse, pre_int_pse, vx_int_pse];
    data_Slope = [data_Slope; pre_pho_slope, vx_pho_slope, pre_int_slope, vx_int_slope];
    valid_subs = valid_subs + 1;
end

fprintf('Included N = %d (Valid in all 4 conditions)\n\n', valid_subs);

% 4. Create Tables and RM-ANOVA Design
if valid_subs > 2
    % Format: Pho_Pre | Pho_VX | Int_Pre | Int_VX
    VarNames = {'Pho_Pre', 'Pho_VX', 'Int_Pre', 'Int_VX'};
    T_PSE   = array2table(data_PSE, 'VariableNames', VarNames);
    T_Slope = array2table(data_Slope, 'VariableNames', VarNames);
    
    % Define Within-Subject Design
    Task    = categorical([1 1 2 2]', [1 2], {'Phoneme', 'Intonation'});
    Session = categorical([1 2 1 2]', [1 2], {'Pretest', 'VX'});
    WithinDesign = table(Task, Session);
    
    % Fit Models
    rm_PSE   = fitrm(T_PSE, 'Pho_Pre-Int_VX ~ 1', 'WithinDesign', WithinDesign);
    rm_Slope = fitrm(T_Slope, 'Pho_Pre-Int_VX ~ 1', 'WithinDesign', WithinDesign);
    
    % Run ANOVAs
    anova_PSE   = ranova(rm_PSE, 'WithinModel', 'Task*Session');
    anova_Slope = ranova(rm_Slope, 'WithinModel', 'Task*Session');
    
    disp('--- RANOVA: PSE (Bias/Threshold) ---');
    disp(anova_PSE);
    
    disp('--- RANOVA: SLOPE ---');
    disp(anova_Slope);
else
    fprintf('Not enough valid subjects to run RM-ANOVA.\n');
end

% --- 6. POST-HOC TESTS (PSE INTERACTION) ---
fprintf('\n=======================================================\n');
fprintf('=== POST-HOC TESTS: PSE INTERACTION                 ===\n');
fprintf('=======================================================\n');

% Daten extrahieren (Spalten aus data_PSE: Pho_Pre, Pho_VX, Int_Pre, Int_VX)
pse_pho_pre = data_PSE(:,1);
pse_pho_vx  = data_PSE(:,2);
pse_int_pre = data_PSE(:,3);
pse_int_vx  = data_PSE(:,4);

% Verschiebungen (Shift) berechnen: Positiv = Verschiebung nach rechts, Negativ = nach links
shift_pho = pse_pho_vx - pse_pho_pre;
shift_int = pse_int_vx - pse_int_pre;

% Gepaarte t-Tests
[~, p_post_pho, ~, stats_pho] = ttest(pse_pho_pre, pse_pho_vx);
[~, p_post_int, ~, stats_int] = ttest(pse_int_pre, pse_int_vx);

% Ausgabe Phoneme
fprintf('--- PHONEME ---\n');
fprintf('Pretest M = %.3f  |  Scanner (VX) M = %.3f\n', mean(pse_pho_pre), mean(pse_pho_vx));
fprintf('Mittlerer Shift = %+.3f  |  t(%d) = %.2f  |  p = %.4f\n\n', mean(shift_pho), stats_pho.df, stats_pho.tstat, p_post_pho);

% Ausgabe Intonation
fprintf('--- INTONATION ---\n');
fprintf('Pretest M = %.3f  |  Scanner (VX) M = %.3f\n', mean(pse_int_pre), mean(pse_int_vx));
fprintf('Mittlerer Shift = %+.3f  |  t(%d) = %.2f  |  p = %.4f\n\n', mean(shift_int), stats_int.df, stats_int.tstat, p_post_int);

end

fprintf('Pipeline Finished.\n');

%% === HELPER FUNCTIONS ===
function [subjects, lin_mean, sig_mean] = extract_adjrsq_clean(FitStruct, include_func)
    all_subjects = unique({FitStruct.subject});
    keep = cellfun(include_func, all_subjects);
    subjects = all_subjects(keep);
    lin_mean = nan(size(subjects)); sig_mean = nan(size(subjects));

    for i = 1:numel(subjects)
        idx = strcmp({FitStruct.subject}, subjects{i});
        entries = FitStruct(idx);
        l_vals = []; s_vals = [];
        for e = 1:numel(entries)
            if isfield(entries(e), 'linfit') && isfield(entries(e), 'sigfit')
                if isnumeric(entries(e).linfit.adjrsq), l_vals = [l_vals; entries(e).linfit.adjrsq(:)]; end
                if isnumeric(entries(e).sigfit.adjrsq), s_vals = [s_vals; entries(e).sigfit.adjrsq(:)]; end
            end
        end
        l_vals = l_vals(isfinite(l_vals)); s_vals = s_vals(isfinite(s_vals));
        if ~isempty(l_vals), lin_mean(i) = mean(l_vals, 'omitnan'); end
        if ~isempty(s_vals), sig_mean(i) = mean(s_vals, 'omitnan'); end
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