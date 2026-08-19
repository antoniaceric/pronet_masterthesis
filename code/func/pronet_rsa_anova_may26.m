function stats_results = pronet_rsa_anova_may26(all_results, VP, out_dir, pipeline_name)
% PRONET_RSA_ANOVA_UNIVERSAL
% Top-Down, Spatial Post-Hocs abd Bottom-Up ANOVAs for Perceptual Similarity Modell.

    %% 1. CONFIGURATION
    tasks    = {'phoneme', 'intonat'};
    sessions = {'lt', 'rt', 'vx'};
    metrics  = {'R_general'}; % Nur noch das General Model
    
    hemis    = {'l', 'r'};                                        
    regions  = {'IFG', 'PMC', 'PAC', 'pSTS', 'aSTS'}; 
    
    rois_list = {'lIFG', 'lPMC', 'lPAC', 'lpSTS', 'laSTS', ...
                 'rIFG', 'rPMC', 'rPAC', 'rpSTS', 'raSTS'};
    
    stats_results = struct();
    
    fprintf('\n=======================================================\n');
    fprintf('   RSA STATISTICAL ANALYSIS: %s        \n', pipeline_name);
    fprintf('=======================================================\n');

    %% ====================================================================
    %  PART 1: TOP-DOWN | 3-Way RM-ANOVA (Hemisphere x Region x Session)
    %  ====================================================================
    fprintf('\n--- [PART 1] TOP-DOWN: 3-Way RM-ANOVA (Hemi x Region x Session) ---\n');
    
    var_names_3way = {}; W_H={}; W_R={}; W_S={};
    for h=1:2
        for r=1:5
            for s=1:3
                var_names_3way{end+1} = sprintf('V_%d_%d_%d', h,r,s);
                W_H{end+1}=hemis{h}; W_R{end+1}=regions{r}; W_S{end+1}=sessions{s};
            end
        end
    end
    wd_3way = table(W_H', W_R', W_S', 'VariableNames', {'Hemisphere', 'Region', 'Session'});

    for m_idx = 1:numel(metrics)
        metric = metrics{m_idx};
        fprintf('\n*******************************************************\n');
        fprintf('   ANALYZING METRIC: %s\n', upper(metric));
        fprintf('*******************************************************\n');
        
        for t = 1:numel(tasks)
            task = tasks{t};
            mat_3way = nan(numel(VP), 30);
            
            for i = 1:numel(VP)
                s_key = ['sub' VP{i}]; col_idx = 1;
                if isfield(all_results, s_key)
                    for h=1:2
                        for r=1:5
                            roi = [hemis{h} regions{r}]; 
                            for s=1:3
                                try
                                    mat_3way(i, col_idx) = all_results.(s_key).(task).(roi).(sessions{s}).(metric);
                                catch; end
                                col_idx = col_idx + 1;
                            end
                        end
                    end
                end
            end
            
            valid_rows = ~any(isnan(mat_3way), 2);
            n_valid = sum(valid_rows);
            
            if n_valid >= 3
                try
                    t_data = array2table(mat_3way(valid_rows,:), 'VariableNames', var_names_3way);
                    rm = fitrm(t_data, sprintf('%s-%s~1',var_names_3way{1},var_names_3way{end}), 'WithinDesign', wd_3way);
                    ra = ranova(rm, 'WithinModel', 'Hemisphere*Region*Session');
                    
                    fprintf('\n>>> TASK: %s | 3-Way RM-ANOVA (N=%d) <<<\n', upper(task), n_valid);
                    disp(ra(ra.pValue < 0.05, {'SumSq','DF','F','pValue','pValueGG','pValueHF'}));
                    
                    stats_results.TopDown.(metric).(task).anova_table = ra;
                catch ME
                     fprintf('3-Way ANOVA failed for %s (%s): %s\n', task, metric, ME.message);
                end
         
                %  PART 1B: EXPLORATIVE EFFECT BREAKDOWNS
                fprintf('\n   --- [PART 1B] EXPLORATIVE BREAKDOWNS (%s) ---\n', upper(task));
                
                i_L_lt = strcmp(W_H, 'l') & strcmp(W_S, 'lt'); i_L_rt = strcmp(W_H, 'l') & strcmp(W_S, 'rt'); i_L_vx = strcmp(W_H, 'l') & strcmp(W_S, 'vx');
                i_R_lt = strcmp(W_H, 'r') & strcmp(W_S, 'lt'); i_R_rt = strcmp(W_H, 'r') & strcmp(W_S, 'rt'); i_R_vx = strcmp(W_H, 'r') & strcmp(W_S, 'vx');
                
                L_lt = mean(mat_3way(valid_rows, i_L_lt), 2); L_rt = mean(mat_3way(valid_rows, i_L_rt), 2); L_vx = mean(mat_3way(valid_rows, i_L_vx), 2);
                R_lt = mean(mat_3way(valid_rows, i_R_lt), 2); R_rt = mean(mat_3way(valid_rows, i_R_rt), 2); R_vx = mean(mat_3way(valid_rows, i_R_vx), 2);
                
                [p_L_lt, t_L_lt] = calc_ttest(L_lt, L_vx); [p_L_rt, t_L_rt] = calc_ttest(L_rt, L_vx);
                [p_R_lt, t_R_lt] = calc_ttest(R_lt, R_vx); [p_R_rt, t_R_rt] = calc_ttest(R_rt, R_vx);
                
                if p_L_lt < 0.05, fprintf('      Links  [LT vs VX]: p = %.4f (t = %.2f)\n', p_L_lt, t_L_lt); end
                if p_L_rt < 0.05, fprintf('      Links  [RT vs VX]: p = %.4f (t = %.2f)\n', p_L_rt, t_L_rt); end
                if p_R_lt < 0.05, fprintf('      Rechts [LT vs VX]: p = %.4f (t = %.2f)\n', p_R_lt, t_R_lt); end
                if p_R_rt < 0.05, fprintf('      Rechts [RT vs VX]: p = %.4f (t = %.2f)\n', p_R_rt, t_R_rt); end

                % ====================================================================
                %  PART 1C: EXHAUSTIVE SPATIAL POST-HOCS (FDR CORRECTED)
                % ====================================================================
                fprintf('\n   --- [PART 1C] EXHAUSTIVE SPATIAL POST-HOCS (%s) ---\n', upper(task));
                
                % 1. OVERALL HEMISPHERE (Left vs Right)
                idx_L = strcmp(W_H, 'l');
                idx_R = strcmp(W_H, 'r');
                val_L = mean(mat_3way(valid_rows, idx_L), 2);
                val_R = mean(mat_3way(valid_rows, idx_R), 2);
                [p_hemi, t_hemi] = calc_ttest(val_L, val_R);
                fprintf('   > HEMISPHERE (Left vs Right): p = %.4f (t = %.2f)\n', p_hemi, t_hemi);
                
                % 2. PAIRWISE REGIONS (10 Comparisons)
                fprintf('\n   > REGION (Pairwise comparisons, FDR corrected):\n');
                reg_pvals = []; reg_tstats = []; reg_pairs = {};
                reg_vals = nan(n_valid, 5);
                for r1 = 1:5
                    reg_vals(:, r1) = mean(mat_3way(valid_rows, strcmp(W_R, regions{r1})), 2);
                end
                
                for r1 = 1:4
                    for r2 = (r1+1):5
                        [p_r, t_r] = calc_ttest(reg_vals(:, r1), reg_vals(:, r2));
                        reg_pvals(end+1) = p_r;
                        reg_tstats(end+1) = t_r;
                        reg_pairs{end+1} = sprintf('%s vs %s', regions{r1}, regions{r2});
                    end
                end
                
                reg_fdr = mafdr(reg_pvals, 'BHFDR', true);
                sig_reg = false;
                for i = 1:10
                    if reg_fdr(i) < 0.05
                        sig_reg = true;
                        fprintf('      * %s: p_FDR = %.4f (p_uncorr = %.4f, t = %.2f)\n', reg_pairs{i}, reg_fdr(i), reg_pvals(i), reg_tstats(i));
                    end
                end
                if ~sig_reg
                    fprintf('      No significant pairwise regional differences after FDR correction.\n');
                end
                
                % 3. HEMISPHERE x REGION (Left vs Right for each specific Region)
                fprintf('\n   > HEMISPHERE x REGION (Left vs Right per Region, FDR corrected):\n');
                hx_pvals = []; hx_tstats = []; hx_names = {};
                for r1 = 1:5
                    idx_Lr = strcmp(W_H, 'l') & strcmp(W_R, regions{r1});
                    idx_Rr = strcmp(W_H, 'r') & strcmp(W_R, regions{r1});
                    val_Lr = mean(mat_3way(valid_rows, idx_Lr), 2);
                    val_Rr = mean(mat_3way(valid_rows, idx_Rr), 2);
                    
                    [p_hx, t_hx] = calc_ttest(val_Lr, val_Rr);
                    hx_pvals(end+1) = p_hx;
                    hx_tstats(end+1) = t_hx;
                    hx_names{end+1} = regions{r1};
                end
                
                hx_fdr = mafdr(hx_pvals, 'BHFDR', true);
                sig_hx = false;
                for i = 1:5
                    if hx_fdr(i) < 0.05
                        sig_hx = true;
                        fprintf('      * %s (Left vs Right): p_FDR = %.4f (p_uncorr = %.4f, t = %.2f)\n', hx_names{i}, hx_fdr(i), hx_pvals(i), hx_tstats(i));
                    end
                end
                if ~sig_hx
                    fprintf('      No significant L/R differences within regions after FDR correction.\n');
                end

            end
        end
    end
    
    %% ====================================================================
    %  PART 2: BOTTOM-UP | 1-Way RM-ANOVA per ROI (FDR Corrected)
    %  ====================================================================
    fprintf('\n--- [PART 2] BOTTOM-UP: 1-Way RM-ANOVA per ROI (10 ANOVAs) ---\n');
    wd_1way = table(sessions', 'VariableNames', {'Session'});
    
    for t = 1:numel(tasks)
        task = tasks{t};
        fprintf('\n>>> TASK: %s | ROI-Specific Session Effects <<<\n', upper(task));
        
        anova_pvals_gg = nan(1, numel(rois_list));
        roi_stats = struct();
        
        for r = 1:numel(rois_list)
            roi = rois_list{r};
            mat_roi = nan(numel(VP), 3); 
            
            for i = 1:numel(VP)
                s_key = ['sub' VP{i}];
                try
                    mat_roi(i, 1) = all_results.(s_key).(task).(roi).lt.(metrics{1});
                    mat_roi(i, 2) = all_results.(s_key).(task).(roi).rt.(metrics{1});
                    mat_roi(i, 3) = all_results.(s_key).(task).(roi).vx.(metrics{1});
                catch; end
            end
            
            valid_rows = ~any(isnan(mat_roi), 2);
            if sum(valid_rows) >= 3
                try
                    t_data = array2table(mat_roi(valid_rows,:), 'VariableNames', {'lt', 'rt', 'vx'});
                    rm = fitrm(t_data, 'lt-vx ~ 1', 'WithinDesign', wd_1way);
                    ra = ranova(rm, 'WithinModel', 'Session');
                    
                    pval_gg = ra{'(Intercept):Session', 'pValueGG'};
                    anova_pvals_gg(r) = pval_gg;
                    roi_stats.(roi).ranova = ra;
                    roi_stats.(roi).mat_roi = mat_roi;
                catch
                end
            end
        end
        
        valid_anova_mask = ~isnan(anova_pvals_gg);
        p_fdr_anova = nan(size(anova_pvals_gg));
        
        if any(valid_anova_mask)
            p_fdr_anova(valid_anova_mask) = mafdr(anova_pvals_gg(valid_anova_mask), 'BHFDR', true);
        end
        
        stats_results.BottomUp.(metrics{1}).(task).pvals_uncorr = anova_pvals_gg;
        stats_results.BottomUp.(metrics{1}).(task).pvals_fdr = p_fdr_anova;
        
        sig_rois = {};
        for r = 1:numel(rois_list)
            if p_fdr_anova(r) < 0.05
                sig_rois{end+1} = rois_list{r};
                fprintf('   ! %s: Significant Session Effect (p_FDR = %.4f, p_uncorr = %.4f)\n', ...
                    rois_list{r}, p_fdr_anova(r), anova_pvals_gg(r));
            end
        end
        
        if isempty(sig_rois)
             fprintf('   No ROIs survived FDR correction for the overall Session effect.\n');
        end

        %  PART 3: POST-HOCS | Only for Significant ROIs (FDR Corrected)
        if ~isempty(sig_rois)
            fprintf('\n   --- [PART 3] Post-Hoc T-Tests for Significant ROIs in %s ---\n', upper(task));
            
            test_results = struct('roi', {}, 'contrast', {}, 'p_uncorr', {}, 'tstat', {}, 'dir_s', {});
            p_vals_ttest = [];
            
            for sr = 1:numel(sig_rois)
                roi = sig_rois{sr};
                mat_roi = roi_stats.(roi).mat_roi;
                
                [p, tstat, diff_mean, valid] = calc_ttest(mat_roi(:,1), mat_roi(:,3));
                if valid
                    dir_s = '>'; if diff_mean < 0, dir_s = '<'; end
                    test_results(end+1) = struct('roi', roi, 'contrast', 'lt_vs_vx', 'p_uncorr', p, 'tstat', tstat, 'dir_s', dir_s);
                    p_vals_ttest(end+1) = p;
                end
                
                [p, tstat, diff_mean, valid] = calc_ttest(mat_roi(:,2), mat_roi(:,3));
                if valid
                    dir_s = '>'; if diff_mean < 0, dir_s = '<'; end
                    test_results(end+1) = struct('roi', roi, 'contrast', 'rt_vs_vx', 'p_uncorr', p, 'tstat', tstat, 'dir_s', dir_s);
                    p_vals_ttest(end+1) = p;
                end
                
                [p, tstat, diff_mean, valid] = calc_ttest(mat_roi(:,1), mat_roi(:,2));
                if valid
                    dir_s = '>'; if diff_mean < 0, dir_s = '<'; end
                    test_results(end+1) = struct('roi', roi, 'contrast', 'lt_vs_rt', 'p_uncorr', p, 'tstat', tstat, 'dir_s', dir_s);
                    p_vals_ttest(end+1) = p;
                end
            end
            
            if ~isempty(p_vals_ttest)
                p_fdr_ttest = mafdr(p_vals_ttest, 'BHFDR', true);
                stats_results.BottomUp.(metrics{1}).(task).posthocs = test_results;
                
                any_sig_ttest = false;
                for k = 1:numel(test_results)
                    if p_fdr_ttest(k) < 0.05
                        any_sig_ttest = true;
                        fprintf('      * %s | %s: %s (p_FDR = %.4f, t = %.2f)  [p_uncorr = %.4f]\n', ...
                            test_results(k).roi, test_results(k).contrast, test_results(k).dir_s, ...
                            p_fdr_ttest(k), test_results(k).tstat, test_results(k).p_uncorr);
                    end
                end
                
                if ~any_sig_ttest
                    fprintf('      No Post-Hoc comparisons survived FDR correction.\n');
                end
            end
        end
    end

    %% ====================================================================
    %  PART 4: VISUALIZATION OF ANOVA EFFECTS (INTERACTION PLOTS)
    %  ====================================================================
    fprintf('\n*******************************************************\n');
    fprintf('   PART 4: GENERATING ANOVA VISUALIZATIONS\n');
    fprintf('*******************************************************\n');
    
    try
        plot_rsa_anova_interaction(all_results, VP, 'intonat', 'R_general', 'l');
        plot_rsa_anova_interaction(all_results, VP, 'intonat', 'R_general', 'r');
    catch ME
        fprintf('  -> Visualization skipped or failed: %s\n', ME.message);
    end

    stats_dir = fullfile(out_dir, 'stats_april26');
    if ~exist(stats_dir, 'dir'), mkdir(stats_dir); end
    
    outfile = fullfile(stats_dir, sprintf('ANOVA_Stats_%s.mat', pipeline_name));
    save(outfile, 'stats_results');
    fprintf('\nALL ANALYSES COMPLETE. Full results saved to:\n%s\n', outfile);

end 
%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function [p, tstat, diff_mean, valid] = calc_ttest(v1, v2)
    mask = ~isnan(v1) & ~isnan(v2);
    if sum(mask) < 3
        p = NaN; tstat = NaN; diff_mean = NaN; valid = false; return;
    end
    
    [~, p, ~, stats] = ttest(v1(mask), v2(mask));
    tstat = stats.tstat;
    diff_mean = mean(v1(mask) - v2(mask));
    valid = true;
end

function plot_rsa_anova_interaction(data_struct, VP, task, metric, target_hemi)
    regions = {'IFG', 'PMC', 'PAC', 'pSTS', 'aSTS'};
    sessions = {'vx', 'lt', 'rt'};
    
    session_colors = [0.5 0.5 0.5; 0.8 0.2 0.2; 0.2 0.4 0.8]; 
    
    n_reg = numel(regions);
    n_sess = numel(sessions);
    
    mean_data = nan(n_sess, n_reg);
    se_data = nan(n_sess, n_reg);
    
    for s = 1:n_sess
        sess = sessions{s};
        mat_sess = nan(numel(VP), n_reg);
        
        for r = 1:n_reg
            roi = [target_hemi regions{r}];
            for i = 1:numel(VP)
                s_key = ['sub' VP{i}];
                try mat_sess(i, r) = data_struct.(s_key).(task).(roi).(sess).(metric); catch; end
            end
        end
        
        valid_rows = ~any(isnan(mat_sess), 2);
        clean_data = mat_sess(valid_rows, :);
        
        mean_data(s, :) = mean(clean_data, 1);
        se_data(s, :) = std(clean_data, 0, 1) / sqrt(sum(valid_rows));
    end
    
    figure('Color', 'w', 'Position', [100, 100, 700, 450], 'Name', sprintf('ANOVA: %s | %s', task, metric), 'Visible', 'off');
    hold on;
    
    x_offsets = [-0.15, 0, 0.15]; 
    for s = 1:n_sess
        x_pos = (1:n_reg) + x_offsets(s);
        errorbar(x_pos, mean_data(s, :), se_data(s, :), ...
            '-o', 'Color', session_colors(s,:), 'LineWidth', 2.5, ...
            'MarkerSize', 7, 'MarkerFaceColor', session_colors(s,:), 'CapSize', 0);
    end
    
    set(gca, 'XTick', 1:n_reg, 'XTickLabel', regions, 'FontSize', 12, 'LineWidth', 1.2);
    xlim([0.5, n_reg + 0.5]);
    
    hemi_str = 'Left Hemisphere'; if strcmpi(target_hemi, 'r'), hemi_str = 'Right Hemisphere'; end
    title(sprintf('Interaction: Region x Session\n%s | %s | %s', upper(task), strrep(metric, '_', ' '), hemi_str), 'FontSize', 14);
    ylabel('Fisher Z (Mean ± SE)', 'FontSize', 12, 'FontWeight', 'bold');
    
    legend(upper(sessions), 'Location', 'best', 'FontSize', 11, 'Box', 'off');
    grid on; set(gca, 'GridColor', [0.9 0.9 0.9]); box off; hold off;
end

% function stats_results = pronet_rsa_anova_may26(all_results, VP, out_dir, pipeline_name)
% % PRONET_RSA_ANOVA_UNIVERSAL
% % Rechnet Top-Down und Bottom-Up ANOVAs für ALLE DREI Modelle 
% % (General, Unique Acoustic, Unique Categorical).
% 
%     %% 1. CONFIGURATION
%     tasks    = {'phoneme', 'intonat'};
%     sessions = {'lt', 'rt', 'vx'};
%     metrics  = {'R_general', 'R_unique_Acoustic', 'R_unique_Categorical'};
% 
%     hemis    = {'l', 'r'};                                    
%     regions  = {'IFG', 'PMC', 'PAC', 'pSTS', 'aSTS'}; 
% 
%     rois_list = {'lIFG', 'lPMC', 'lPAC', 'lpSTS', 'laSTS', ...
%                  'rIFG', 'rPMC', 'rPAC', 'rpSTS', 'raSTS'};
% 
%     stats_results = struct();
% 
%     fprintf('\n=======================================================\n');
%     fprintf('   RSA STATISTICAL ANALYSIS: %s        \n', pipeline_name);
%     fprintf('=======================================================\n');
% 
%     %% ====================================================================
%     %  PART 1: TOP-DOWN | 3-Way RM-ANOVA (Hemisphere x Region x Session)
%     %  ====================================================================
%     fprintf('\n--- [PART 1] TOP-DOWN: 3-Way RM-ANOVA (Hemi x Region x Session) ---\n');
% 
%     var_names_3way = {}; W_H={}; W_R={}; W_S={};
%     for h=1:2
%         for r=1:5
%             for s=1:3
%                 var_names_3way{end+1} = sprintf('V_%d_%d_%d', h,r,s);
%                 W_H{end+1}=hemis{h}; W_R{end+1}=regions{r}; W_S{end+1}=sessions{s};
%             end
%         end
%     end
%     wd_3way = table(W_H', W_R', W_S', 'VariableNames', {'Hemisphere', 'Region', 'Session'});
% 
%     for m_idx = 1:numel(metrics)
%         metric = metrics{m_idx};
%         fprintf('\n*******************************************************\n');
%         fprintf('   ANALYZING METRIC: %s\n', upper(metric));
%         fprintf('*******************************************************\n');
% 
%         for t = 1:numel(tasks)
%             task = tasks{t};
%             mat_3way = nan(numel(VP), 30);
% 
%             for i = 1:numel(VP)
%                 s_key = ['sub' VP{i}]; col_idx = 1;
%                 if isfield(all_results, s_key)
%                     for h=1:2
%                         for r=1:5
%                             roi = [hemis{h} regions{r}]; 
%                             for s=1:3
%                                 try
%                                     mat_3way(i, col_idx) = all_results.(s_key).(task).(roi).(sessions{s}).(metric);
%                                 catch; end
%                                 col_idx = col_idx + 1;
%                             end
%                         end
%                     end
%                 end
%             end
% 
%             valid_rows = ~any(isnan(mat_3way), 2);
%             n_valid = sum(valid_rows);
% 
%             if n_valid >= 3
%                 try
%                     t_data = array2table(mat_3way(valid_rows,:), 'VariableNames', var_names_3way);
%                     rm = fitrm(t_data, sprintf('%s-%s~1',var_names_3way{1},var_names_3way{end}), 'WithinDesign', wd_3way);
%                     ra = ranova(rm, 'WithinModel', 'Hemisphere*Region*Session');
% 
%                     fprintf('\n>>> TASK: %s | 3-Way RM-ANOVA (N=%d) <<<\n', upper(task), n_valid);
%                     disp(ra(ra.pValue < 0.05, {'SumSq','DF','F','pValue','pValueGG','pValueHF'}));
% 
%                     % Save
%                     stats_results.TopDown.(metric).(task).anova_table = ra;
%                 catch ME
%                      fprintf('3-Way ANOVA failed for %s (%s): %s\n', task, metric, ME.message);
%                 end
% 
% %% ====================================================================
%                 %  PART 1B: EXPLORATIVE EFFECT BREAKDOWNS (Unkorrigierte Post-Hocs)
%                 %  ====================================================================
%                 fprintf('\n   --- [PART 1B] EXPLORATIVE BREAKDOWNS (%s) ---\n', upper(task));
% 
%                 % 1. HEMISPHERE x SESSION (averaged over all regions)
%                 fprintf('   > 1. HEMISPHERE x SESSION (mean over regions):\n');
%                 i_L_lt = strcmp(W_H, 'l') & strcmp(W_S, 'lt'); i_L_rt = strcmp(W_H, 'l') & strcmp(W_S, 'rt'); i_L_vx = strcmp(W_H, 'l') & strcmp(W_S, 'vx');
%                 i_R_lt = strcmp(W_H, 'r') & strcmp(W_S, 'lt'); i_R_rt = strcmp(W_H, 'r') & strcmp(W_S, 'rt'); i_R_vx = strcmp(W_H, 'r') & strcmp(W_S, 'vx');
% 
%                 L_lt = mean(mat_3way(valid_rows, i_L_lt), 2); L_rt = mean(mat_3way(valid_rows, i_L_rt), 2); L_vx = mean(mat_3way(valid_rows, i_L_vx), 2);
%                 R_lt = mean(mat_3way(valid_rows, i_R_lt), 2); R_rt = mean(mat_3way(valid_rows, i_R_rt), 2); R_vx = mean(mat_3way(valid_rows, i_R_vx), 2);
% 
%                 [p_L_lt, t_L_lt] = calc_ttest(L_lt, L_vx); [p_L_rt, t_L_rt] = calc_ttest(L_rt, L_vx);
%                 [p_R_lt, t_R_lt] = calc_ttest(R_lt, R_vx); [p_R_rt, t_R_rt] = calc_ttest(R_rt, R_vx);
% 
%                 if p_L_lt < 0.05, fprintf('      Links  [LT vs VX]: p = %.4f (t = %.2f)\n', p_L_lt, t_L_lt); end
%                 if p_L_rt < 0.05, fprintf('      Links  [RT vs VX]: p = %.4f (t = %.2f)\n', p_L_rt, t_L_rt); end
%                 if p_R_lt < 0.05, fprintf('      Rechts [LT vs VX]: p = %.4f (t = %.2f)\n', p_R_lt, t_R_lt); end
%                 if p_R_rt < 0.05, fprintf('      Rechts [RT vs VX]: p = %.4f (t = %.2f)\n', p_R_rt, t_R_rt); end
% 
%                 % 2. REGION x SESSION (averaged over hemispheres)
%                 fprintf('\n   > 2. REGION x SESSION (mean over hemispheres):\n');
%                 for r_idx = 1:5
%                     reg = regions{r_idx};
%                     i_lt = strcmp(W_R, reg) & strcmp(W_S, 'lt');
%                     i_rt = strcmp(W_R, reg) & strcmp(W_S, 'rt');
%                     i_vx = strcmp(W_R, reg) & strcmp(W_S, 'vx');
% 
%                     val_lt = mean(mat_3way(valid_rows, i_lt), 2);
%                     val_rt = mean(mat_3way(valid_rows, i_rt), 2);
%                     val_vx = mean(mat_3way(valid_rows, i_vx), 2);
% 
%                     [p_lt, t_lt] = calc_ttest(val_lt, val_vx);
%                     [p_rt, t_rt] = calc_ttest(val_rt, val_vx);
% 
%                     if p_lt < 0.05, fprintf('      Region %s [LT vs VX]: p = %.4f (t = %.2f)\n', reg, p_lt, t_lt); end
%                     if p_rt < 0.05, fprintf('      Region %s [RT vs VX]: p = %.4f (t = %.2f)\n', reg, p_rt, t_rt); end
%                 end
% 
%                 % 3. HEMISPHERE x REGION x SESSION (per roi)
%                 fprintf('\n   > 3. 3-WAY ISOLATED (per roi):\n');
%                 for h_idx = 1:2
%                     for r_idx = 1:5
%                         roi_name = [hemis{h_idx}, regions{r_idx}];
%                         i_lt = strcmp(W_H, hemis{h_idx}) & strcmp(W_R, regions{r_idx}) & strcmp(W_S, 'lt');
%                         i_rt = strcmp(W_H, hemis{h_idx}) & strcmp(W_R, regions{r_idx}) & strcmp(W_S, 'rt');
%                         i_vx = strcmp(W_H, hemis{h_idx}) & strcmp(W_R, regions{r_idx}) & strcmp(W_S, 'vx');
% 
%                         val_lt = mat_3way(valid_rows, i_lt);
%                         val_rt = mat_3way(valid_rows, i_rt);
%                         val_vx = mat_3way(valid_rows, i_vx);
% 
%                         [p_lt, t_lt] = calc_ttest(val_lt, val_vx);
%                         [p_rt, t_rt] = calc_ttest(val_rt, val_vx);
% 
%                         if p_lt < 0.05, fprintf('      ROI %s [LT vs VX]: p = %.4f (t = %.2f)\n', roi_name, p_lt, t_lt); end
%                         if p_rt < 0.05, fprintf('      ROI %s [RT vs VX]: p = %.4f (t = %.2f)\n', roi_name, p_rt, t_rt); end
%                     end
%                 end
%         end
%     end
%         %% ====================================================================
%         %  PART 2: BOTTOM-UP | 1-Way RM-ANOVA per ROI (FDR Corrected)
%         %  ====================================================================
%         fprintf('\n--- [PART 2] BOTTOM-UP: 1-Way RM-ANOVA per ROI (10 ANOVAs) ---\n');
%         wd_1way = table(sessions', 'VariableNames', {'Session'});
% 
%         for t = 1:numel(tasks)
%             task = tasks{t};
%             fprintf('\n>>> TASK: %s | ROI-Specific Session Effects <<<\n', upper(task));
% 
%             anova_pvals_gg = nan(1, numel(rois_list));
%             roi_stats = struct();
% 
%             for r = 1:numel(rois_list)
%                 roi = rois_list{r};
%                 mat_roi = nan(numel(VP), 3); 
% 
%                 for i = 1:numel(VP)
%                     s_key = ['sub' VP{i}];
%                     try
%                         mat_roi(i, 1) = all_results.(s_key).(task).(roi).lt.(metric);
%                         mat_roi(i, 2) = all_results.(s_key).(task).(roi).rt.(metric);
%                         mat_roi(i, 3) = all_results.(s_key).(task).(roi).vx.(metric);
%                     catch; end
%                 end
% 
%                 valid_rows = ~any(isnan(mat_roi), 2);
%                 if sum(valid_rows) >= 3
%                     try
%                         t_data = array2table(mat_roi(valid_rows,:), 'VariableNames', {'lt', 'rt', 'vx'});
%                         rm = fitrm(t_data, 'lt-vx ~ 1', 'WithinDesign', wd_1way);
%                         ra = ranova(rm, 'WithinModel', 'Session');
% 
%                         pval_gg = ra{'(Intercept):Session', 'pValueGG'};
%                         anova_pvals_gg(r) = pval_gg;
%                         roi_stats.(roi).ranova = ra;
%                         roi_stats.(roi).mat_roi = mat_roi;
%                     catch
%                     end
%                 end
%             end
% 
%             valid_anova_mask = ~isnan(anova_pvals_gg);
%             p_fdr_anova = nan(size(anova_pvals_gg));
% 
%             if any(valid_anova_mask)
%                 p_fdr_anova(valid_anova_mask) = mafdr(anova_pvals_gg(valid_anova_mask), 'BHFDR', true);
%             end
% 
%             stats_results.BottomUp.(metric).(task).pvals_uncorr = anova_pvals_gg;
%             stats_results.BottomUp.(metric).(task).pvals_fdr = p_fdr_anova;
% 
%             sig_rois = {};
%             for r = 1:numel(rois_list)
%                 if p_fdr_anova(r) < 0.05
%                     sig_rois{end+1} = rois_list{r};
%                     fprintf('   ? %s: Significant Session Effect (p_FDR = %.4f, p_uncorr = %.4f)\n', ...
%                         rois_list{r}, p_fdr_anova(r), anova_pvals_gg(r));
%                 end
%             end
% 
%             if isempty(sig_rois)
%                  fprintf('   No ROIs survived FDR correction for the overall Session effect.\n');
%             end
% 
%             %% ====================================================================
%             %  PART 3: POST-HOCS | Only for Significant ROIs (FDR Corrected)
%             %  ====================================================================
%             if ~isempty(sig_rois)
%                 fprintf('\n   --- [PART 3] Post-Hoc T-Tests for Significant ROIs in %s ---\n', upper(task));
% 
%                 test_results = struct('roi', {}, 'contrast', {}, 'p_uncorr', {}, 'tstat', {}, 'dir_s', {});
%                 p_vals_ttest = [];
% 
%                 for sr = 1:numel(sig_rois)
%                     roi = sig_rois{sr};
%                     mat_roi = roi_stats.(roi).mat_roi;
% 
%                     [p, tstat, diff_mean, valid] = calc_ttest(mat_roi(:,1), mat_roi(:,3));
%                     if valid
%                         dir_s = '>'; if diff_mean < 0, dir_s = '<'; end
%                         test_results(end+1) = struct('roi', roi, 'contrast', 'lt_vs_vx', 'p_uncorr', p, 'tstat', tstat, 'dir_s', dir_s);
%                         p_vals_ttest(end+1) = p;
%                     end
% 
%                     [p, tstat, diff_mean, valid] = calc_ttest(mat_roi(:,2), mat_roi(:,3));
%                     if valid
%                         dir_s = '>'; if diff_mean < 0, dir_s = '<'; end
%                         test_results(end+1) = struct('roi', roi, 'contrast', 'rt_vs_vx', 'p_uncorr', p, 'tstat', tstat, 'dir_s', dir_s);
%                         p_vals_ttest(end+1) = p;
%                     end
% 
%                     [p, tstat, diff_mean, valid] = calc_ttest(mat_roi(:,1), mat_roi(:,2));
%                     if valid
%                         dir_s = '>'; if diff_mean < 0, dir_s = '<'; end
%                         test_results(end+1) = struct('roi', roi, 'contrast', 'lt_vs_rt', 'p_uncorr', p, 'tstat', tstat, 'dir_s', dir_s);
%                         p_vals_ttest(end+1) = p;
%                     end
%                 end
% 
%                 if ~isempty(p_vals_ttest)
%                     p_fdr_ttest = mafdr(p_vals_ttest, 'BHFDR', true);
%                     stats_results.BottomUp.(metric).(task).posthocs = test_results;
% 
%                     any_sig_ttest = false;
%                     for k = 1:numel(test_results)
%                         if p_fdr_ttest(k) < 0.05
%                             any_sig_ttest = true;
%                             fprintf('      * %s | %s: %s (p_FDR = %.4f, t = %.2f)  [p_uncorr = %.4f]\n', ...
%                                 test_results(k).roi, test_results(k).contrast, test_results(k).dir_s, ...
%                                 p_fdr_ttest(k), test_results(k).tstat, test_results(k).p_uncorr);
%                         end
%                     end
% 
%                     if ~any_sig_ttest
%                         fprintf('      No Post-Hoc comparisons survived FDR correction.\n');
%                     end
%                 end
%             end
%         end
%     end 
% 
% %% ====================================================================
%     %  PART 3: MODEL COMPARISON | 3-Way RM-ANOVA (Acoustic vs. Categorical)
%     %  ====================================================================
%     fprintf('\n*******************************************************\n');
%     fprintf('   PART 3: MODEL COMPARISON ANOVA (ACOUSTIC - CATEGORICAL)\n');
%     fprintf('*******************************************************\n');
% 
%     for t = 1:numel(tasks)
%         task = tasks{t};
%         % Datenmatrix für Differenzwerte: (Subjekte x 30 Bedingungen)
%         mat_diff = nan(numel(VP), 30);
% 
%         for i = 1:numel(VP)
%             s_key = ['sub' VP{i}]; col_idx = 1;
%             if isfield(all_results, s_key)
%                 for h=1:2
%                     for r=1:5
%                         roi = [hemis{h} regions{r}]; 
%                         for s=1:3
%                             try
%                                 % Differenz berechnen: Acoustic - Categorical
%                                 val_aco = all_results.(s_key).(task).(roi).(sessions{s}).R_unique_Acoustic;
%                                 val_cat = all_results.(s_key).(task).(roi).(sessions{s}).R_unique_Categorical;
%                                 mat_diff(i, col_idx) = val_aco - val_cat;
%                             catch; end
%                             col_idx = col_idx + 1;
%                         end
%                     end
%                 end
%             end
%         end
% 
%         valid_rows = ~any(isnan(mat_diff), 2);
%         n_valid = sum(valid_rows);
% 
%         if n_valid >= 3
%             try
%                 t_data = array2table(mat_diff(valid_rows,:), 'VariableNames', var_names_3way);
%                 rm = fitrm(t_data, sprintf('%s-%s~1',var_names_3way{1},var_names_3way{end}), 'WithinDesign', wd_3way);
%                 ra = ranova(rm, 'WithinModel', 'Hemisphere*Region*Session');
% 
%                 fprintf('\n>>> TASK: %s | Model Comparison ANOVA (Aco-Cat) (N=%d) <<<\n', upper(task), n_valid);
%                 disp(ra(ra.pValue < 0.05, {'SumSq','DF','F','pValue','pValueGG'}));
% 
%                 stats_results.ModelComparison.(task).anova_table = ra;
% 
%                 % --- [PART 3B] DETAILED BREAKDOWN OF MODEL PREFERENCE (%s) ---
%                 fprintf('\n   --- [PART 3B] BREAKDOWN OF MODEL PREFERENCE (%s) ---\n', upper(task));
% 
%                 % 1. HEMISPHERE-LEVEL (mean over regions)
%                 fprintf('   > 1. HEMISPHERE x SESSION (Preference Shifts):\n');
%                 for h_idx = 1:2
%                     h_name = hemis{h_idx};
%                     idx_vx = strcmp(W_H, h_name) & strcmp(W_S, 'vx');
%                     idx_lt = strcmp(W_H, h_name) & strcmp(W_S, 'lt');
%                     idx_rt = strcmp(W_H, h_name) & strcmp(W_S, 'rt');
% 
%                     val_vx = mean(mat_diff(valid_rows, idx_vx), 2);
%                     val_lt = mean(mat_diff(valid_rows, idx_lt), 2);
%                     val_rt = mean(mat_diff(valid_rows, idx_rt), 2);
% 
%                     % Test gegen Null (Baseline Präferenz)
%                     [p_base] = calc_ttest(val_vx, zeros(size(val_vx)));
%                     if p_base < 0.10
%                         type = 'SIG'; if p_base >= 0.05, type = 'TREND'; end
%                         pref = 'ACOUSTIC'; if mean(val_vx) < 0, pref = 'CATEGORICAL'; end
%                         fprintf('      [%s] Hemi %s [Baseline VX]: %s Preference for %s (p=%.4f)\n', type, h_name, type, pref, p_base);
%                     end
% 
%                     % Test Shifts
%                     [p_lt] = calc_ttest(val_lt, val_vx);
%                     [p_rt] = calc_ttest(val_rt, val_vx);
%                     if p_lt < 0.10
%                         type = 'SIG'; if p_lt >= 0.05, type = 'TREND'; end
%                         fprintf('      [%s] Hemi %s [LT vs VX]: TMS shifted Aco/Cat balance (p=%.4f)\n', type, h_name, p_lt);
%                     end
%                     if p_rt < 0.10
%                         type = 'SIG'; if p_rt >= 0.05, type = 'TREND'; end
%                         fprintf('      [%s] Hemi %s [RT vs VX]: TMS shifted Aco/Cat balance (p=%.4f)\n', type, h_name, p_rt);
%                     end
%                 end
% 
%                 % 2. REGION-LEVEL (mean over hemispheres)
%                 fprintf('\n   > 2. REGION x SESSION (Preference Shifts):\n');
%                 for r_idx = 1:5
%                     reg = regions{r_idx};
%                     idx_vx = strcmp(W_R, reg) & strcmp(W_S, 'vx');
%                     idx_lt = strcmp(W_R, reg) & strcmp(W_S, 'lt');
%                     idx_rt = strcmp(W_R, reg) & strcmp(W_S, 'rt');
% 
%                     val_vx = mean(mat_diff(valid_rows, idx_vx), 2);
%                     val_lt = mean(mat_diff(valid_rows, idx_lt), 2);
%                     val_rt = mean(mat_diff(valid_rows, idx_rt), 2);
% 
%                     [p_lt] = calc_ttest(val_lt, val_vx);
%                     [p_rt] = calc_ttest(val_rt, val_vx);
% 
%                     if p_lt < 0.10
%                         type = 'SIG'; if p_lt >= 0.05, type = 'TREND'; end
%                         fprintf('      [%s] Region %s [LT vs VX]: Shift in Aco/Cat balance (p=%.4f)\n', type, reg, p_lt);
%                     end
%                     if p_rt < 0.10
%                         type = 'SIG'; if p_rt >= 0.05, type = 'TREND'; end
%                         fprintf('      [%s] Region %s [RT vs VX]: Shift in Aco/Cat balance (p=%.4f)\n', type, reg, p_rt);
%                     end
%                 end
% 
%                 % 3. ISOLATED ROI-LEVEL (10 specific ROIs)
%                 fprintf('\n   > 3. ISOLATED ROI COMPARISONS:\n');
%                 for h_idx = 1:2
%                     for r_idx = 1:5
%                         roi_name = [hemis{h_idx}, regions{r_idx}];
%                         idx_vx = strcmp(W_H, hemis{h_idx}) & strcmp(W_R, regions{r_idx}) & strcmp(W_S, 'vx');
%                         idx_lt = strcmp(W_H, hemis{h_idx}) & strcmp(W_R, regions{r_idx}) & strcmp(W_S, 'lt');
%                         idx_rt = strcmp(W_H, hemis{h_idx}) & strcmp(W_R, regions{r_idx}) & strcmp(W_S, 'rt');
% 
%                         v_vx = mat_diff(valid_rows, idx_vx);
%                         v_lt = mat_diff(valid_rows, idx_lt);
%                         v_rt = mat_diff(valid_rows, idx_rt);
% 
%                         % T-Tests berechnen (mit T-Statistik und Richtung)
%                         [p_lt, t_lt] = calc_ttest(v_lt, v_vx);
%                         [p_rt, t_rt] = calc_ttest(v_rt, v_vx);
% 
%                         % LT vs VX
%                         if p_lt < 0.10
%                             type = 'SIG'; if p_lt >= 0.05, type = 'TREND'; end
%                             dir_str = 'toward CATEGORICAL'; if t_lt > 0, dir_str = 'toward ACOUSTIC'; end
%                             fprintf('      [%s] ROI %s [LT vs VX]: Shift %s (p=%.4f, t=%.2f)\n', type, roi_name, dir_str, p_lt, t_lt);
% 
%                             % Speichern im Struct
%                             stats_results.ModelComparison.(task).PostHoc.ROI.(roi_name).lt_vs_vx.p = p_lt;
%                             stats_results.ModelComparison.(task).PostHoc.ROI.(roi_name).lt_vs_vx.t = t_lt;
%                             stats_results.ModelComparison.(task).PostHoc.ROI.(roi_name).lt_vs_vx.dir = dir_str;
%                         end
% 
%                         % RT vs VX
%                         if p_rt < 0.10
%                             type = 'SIG'; if p_rt >= 0.05, type = 'TREND'; end
%                             dir_str = 'toward CATEGORICAL'; if t_rt > 0, dir_str = 'toward ACOUSTIC'; end
%                             fprintf('      [%s] ROI %s [RT vs VX]: Shift %s (p=%.4f, t=%.2f)\n', type, roi_name, dir_str, p_rt, t_rt);
% 
%                             % Speichern im Struct
%                             stats_results.ModelComparison.(task).PostHoc.ROI.(roi_name).rt_vs_vx.p = p_rt;
%                             stats_results.ModelComparison.(task).PostHoc.ROI.(roi_name).rt_vs_vx.t = t_rt;
%                             stats_results.ModelComparison.(task).PostHoc.ROI.(roi_name).rt_vs_vx.dir = dir_str;
%                         end
%                     end
%                 end
%             catch ME
%                  fprintf('Model Comparison ANOVA failed for %s: %s\n', task, ME.message);
%             end
% 
%         end
%     end
% 
%     %% ====================================================================
%     %  PART 4: VISUALIZATION OF ANOVA EFFECTS (INTERACTION PLOTS)
%     %  ====================================================================
%     fprintf('\n*******************************************************\n');
%     fprintf('   PART 4: GENERATING ANOVA VISUALIZATIONS\n');
%     fprintf('*******************************************************\n');
% 
%     % Hier rufen wir die Visualisierung für die spannendsten Interaktionen auf.
%     % (Passe die Strings an, falls du andere Tasks/Metriken primär sehen willst).
%     try
%         % 1. Visualisierung: Single Model (z.B. Acoustic Model für Intonation)
%         % Spiegelt deine 3-Way ANOVA Hemi:Region:Session Interaktion wider
%         plot_rsa_anova_interaction(all_results, VP, 'intonat', 'R_unique_Acoustic', 'l');
%         plot_rsa_anova_interaction(all_results, VP, 'intonat', 'R_unique_Acoustic', 'r');
% 
%         % 2. Visualisierung: Model Comparison (Aco-Cat Differenz für Phoneme)
%         % Spiegelt den PAC-Ausreißer unter Left-TMS wider
%         plot_rsa_anova_interaction_diff(all_results, VP, 'phoneme', 'l');
% 
%     catch ME
%         fprintf('  -> Visualization skipped or failed: %s\n', ME.message);
%     end
% 
% 
%     stats_dir = fullfile(out_dir, 'stats_april26');
%     if ~exist(stats_dir, 'dir'), mkdir(stats_dir); end
% 
%     outfile = fullfile(stats_dir, sprintf('ANOVA_Stats_%s.mat', pipeline_name));
%     save(outfile, 'stats_results');
%     fprintf('\nALL ANALYSES COMPLETE. Full results saved to:\n%s\n', outfile);
% 
% 
%     end % new
% %% ========================================================================
% %  LOCAL HELPER FUNCTIONS
% %  ========================================================================
% 
% function [p, tstat, diff_mean, valid] = calc_ttest(v1, v2)
%     mask = ~isnan(v1) & ~isnan(v2);
%     if sum(mask) < 3
%         p = NaN; tstat = NaN; diff_mean = NaN; valid = false; return;
%     end
% 
%     [~, p, ~, stats] = ttest(v1(mask), v2(mask));
%     tstat = stats.tstat;
%     diff_mean = mean(v1(mask) - v2(mask));
%     valid = true;
% end
% 
% % =========================================================================
% %  VISUALIZATION HELPER FUNCTIONS (ANOVA Interactions)
% % =========================================================================
% 
% function plot_rsa_anova_interaction(data_struct, VP, task, metric, target_hemi)
% % Creates an interaction plot (Region x Session) for a specific hemisphere.
% % Error bars represent Standard Error (Within-Group variance).
% 
%     regions = {'IFG', 'PMC', 'PAC', 'pSTS', 'aSTS'};
%     sessions = {'vx', 'lt', 'rt'};
% 
%     % Farben: VX (Grau), LT (Rot/Triangle-like theme), RT (Blau/Square-like theme)
%     session_colors = [0.5 0.5 0.5; 
%                       0.8 0.2 0.2; 
%                       0.2 0.4 0.8]; 
% 
%     n_reg = numel(regions);
%     n_sess = numel(sessions);
% 
%     mean_data = nan(n_sess, n_reg);
%     se_data = nan(n_sess, n_reg);
% 
%     for s = 1:n_sess
%         sess = sessions{s};
%         mat_sess = nan(numel(VP), n_reg);
% 
%         for r = 1:n_reg
%             roi = [target_hemi regions{r}];
%             for i = 1:numel(VP)
%                 s_key = ['sub' VP{i}];
%                 try mat_sess(i, r) = data_struct.(s_key).(task).(roi).(sess).(metric); catch; end
%             end
%         end
% 
%         valid_rows = ~any(isnan(mat_sess), 2);
%         clean_data = mat_sess(valid_rows, :);
% 
%         mean_data(s, :) = mean(clean_data, 1);
%         se_data(s, :) = std(clean_data, 0, 1) / sqrt(sum(valid_rows));
%     end
% 
%     figure('Color', 'w', 'Position', [100, 100, 700, 450], 'Name', sprintf('ANOVA: %s | %s', task, metric));
%     hold on;
% 
%     % X-Offsets to prevent error bar overlap
%     x_offsets = [-0.15, 0, 0.15]; 
% 
%     for s = 1:n_sess
%         x_pos = (1:n_reg) + x_offsets(s);
%         errorbar(x_pos, mean_data(s, :), se_data(s, :), ...
%             '-o', 'Color', session_colors(s,:), 'LineWidth', 2.5, ...
%             'MarkerSize', 7, 'MarkerFaceColor', session_colors(s,:), 'CapSize', 0);
%     end
% 
%     set(gca, 'XTick', 1:n_reg, 'XTickLabel', regions, 'FontSize', 12, 'LineWidth', 1.2);
%     xlim([0.5, n_reg + 0.5]);
% 
%     hemi_str = 'Left Hemisphere'; if strcmpi(target_hemi, 'r'), hemi_str = 'Right Hemisphere'; end
%     title(sprintf('Interaction: Region x Session\n%s | %s | %s', upper(task), strrep(metric, '_', ' '), hemi_str), 'FontSize', 14);
%     ylabel('Fisher Z (Mean ± SE)', 'FontSize', 12, 'FontWeight', 'bold');
% 
%     legend(upper(sessions), 'Location', 'best', 'FontSize', 11, 'Box', 'off');
%     grid on; set(gca, 'GridColor', [0.9 0.9 0.9]); box off; hold off;
% end
% 
% % -------------------------------------------------------------------------
% 
% function plot_rsa_anova_interaction_diff(data_struct, VP, task, target_hemi)
% % Creates an interaction plot for the Model Comparison (Acoustic - Categorical).
% 
%     regions = {'IFG', 'PMC', 'PAC', 'pSTS', 'aSTS'};
%     sessions = {'vx', 'lt', 'rt'};
%     session_colors = [0.5 0.5 0.5; 0.8 0.2 0.2; 0.2 0.4 0.8]; 
% 
%     n_reg = numel(regions);
%     n_sess = numel(sessions);
% 
%     mean_data = nan(n_sess, n_reg);
%     se_data = nan(n_sess, n_reg);
% 
%     for s = 1:n_sess
%         sess = sessions{s};
%         mat_sess = nan(numel(VP), n_reg);
% 
%         for r = 1:n_reg
%             roi = [target_hemi regions{r}];
%             for i = 1:numel(VP)
%                 s_key = ['sub' VP{i}];
%                 try 
%                     v_aco = data_struct.(s_key).(task).(roi).(sess).R_unique_Acoustic;
%                     v_cat = data_struct.(s_key).(task).(roi).(sess).R_unique_Categorical;
%                     mat_sess(i, r) = v_aco - v_cat; 
%                 catch; end
%             end
%         end
% 
%         valid_rows = ~any(isnan(mat_sess), 2);
%         clean_data = mat_sess(valid_rows, :);
% 
%         mean_data(s, :) = mean(clean_data, 1);
%         se_data(s, :) = std(clean_data, 0, 1) / sqrt(sum(valid_rows));
%     end
% 
%     figure('Color', 'w', 'Position', [150, 150, 700, 450], 'Name', sprintf('ANOVA: Aco-Cat | %s', task));
%     hold on;
% 
%     % Null-Linie (0 = keine Präferenz)
%     yline(0, '--k', 'LineWidth', 1.5, 'Alpha', 0.5);
% 
%     x_offsets = [-0.15, 0, 0.15]; 
%     for s = 1:n_sess
%         x_pos = (1:n_reg) + x_offsets(s);
%         errorbar(x_pos, mean_data(s, :), se_data(s, :), ...
%             '-o', 'Color', session_colors(s,:), 'LineWidth', 2.5, ...
%             'MarkerSize', 7, 'MarkerFaceColor', session_colors(s,:), 'CapSize', 0);
%     end
% 
%     set(gca, 'XTick', 1:n_reg, 'XTickLabel', regions, 'FontSize', 12, 'LineWidth', 1.2);
%     xlim([0.5, n_reg + 0.5]);
% 
%     hemi_str = 'Left Hemisphere'; if strcmpi(target_hemi, 'r'), hemi_str = 'Right Hemisphere'; end
%     title(sprintf('Model Comparison (Acoustic - Categorical)\n%s | %s', upper(task), hemi_str), 'FontSize', 14);
%     ylabel('\Delta Fisher Z (Mean ± SE)', 'FontSize', 12, 'FontWeight', 'bold');
% 
%     % Text für die y-Achse zur Orientierung
%     text(0.6, max(ylim)*0.9, 'Acoustic > Categorical', 'FontSize', 10, 'Color', [0.3 0.3 0.3]);
%     text(0.6, min(ylim)*0.9, 'Categorical > Acoustic', 'FontSize', 10, 'Color', [0.3 0.3 0.3]);
% 
%     legend(['Baseline (0)', upper(sessions)], 'Location', 'best', 'FontSize', 11, 'Box', 'off');
%     grid on; set(gca, 'GridColor', [0.9 0.9 0.9]); box off; hold off;
% end
% 

