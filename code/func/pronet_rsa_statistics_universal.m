
function [stats_tables] = pronet_rsa_statistics_universal(all_subject_results, VP, rois, tasks, out_dir, pipeline_name)
% Performs session-wise group-level statistics (non-parametric one-sample t-test)
% Generates STANDARD plots and STACKED BASELINE OVERLAY plots.

fprintf('\n--- Starting Session-wise Statistical Analysis: %s ---\n', pipeline_name);

n_subjects = numel(VP);
n_rois = numel(rois);
n_tasks = numel(tasks);
sessions = {'lt', 'rt', 'vx'};
n_sessions = numel(sessions);
n_perms = 100000; 

plot_y_limit = 0.25; 
plot_y_limit_stacked = 0.25;

% --- 1. Reshape Data ---
flat_data = struct(); 

for t = 1:n_tasks
    task = tasks{t};
    for s = 1:n_sessions
        sess = sessions{s};
        R_matrix = zeros(n_subjects, n_rois, 1); % Signal (Corrected)
        B_matrix = zeros(n_subjects, n_rois, 1); % Baseline (Noise Floor)
        
        for i = 1:n_subjects
            subj_id = ['sub' VP{i}];
            for r = 1:n_rois
                roi_name = rois{r};
                if isfield(all_subject_results, subj_id) && ...
                   isfield(all_subject_results.(subj_id), task) && ...
                   isfield(all_subject_results.(subj_id).(task), roi_name) && ...
                   isfield(all_subject_results.(subj_id).(task).(roi_name), sess)
               
                    dat = all_subject_results.(subj_id).(task).(roi_name).(sess);
                    
                    % Corrected Signal
                    R_matrix(i, r, 1) = dat.R_general;
                    
                    % Null Baseline
                    if isfield(dat, 'Baseline_Z_general')
                        B_matrix(i, r, 1) = dat.Baseline_Z_general;
                    end
                end
            end
        end
        flat_data.(task).(sess).R_matrix = R_matrix;
        flat_data.(task).(sess).B_matrix = B_matrix;
    end
end
fprintf('Data reshaping complete.\n');

% --- 2. Perform Statistics (On Corrected Signal) ---
stats_tables = struct();
measure_names = {'Perceptual_Similarity'};

for t = 1:n_tasks
    task = tasks{t};
    for s = 1:n_sessions
        sess = sessions{s};
        fprintf('Calculating stats: Task %s | Session %s\n', upper(task), upper(sess));
        
        stats_table = table('RowNames', rois);
        for m = 1:1  
            measure = measure_names{m};
            data_to_test = flat_data.(task).(sess).R_matrix(:, :, m);
            
            % Run Permutation Test
            [t_stats, p_adj, ~, ~] = permuttest(data_to_test, 0, ...
                                               'nperm', n_perms, ...
                                               'tail', 'right', ...
                                               'correct', true, ...
                                               'verbose', 0);
            
            stats_table.([measure '_tstat']) = t_stats';
            stats_table.([measure '_p_fwer']) = p_adj';
        end
        stats_tables.(task).(sess) = stats_table;
    end
end
fprintf('Statistical testing complete.\n');

%% --- 3. Visualize Results ---
fig_outdir = fullfile(out_dir, 'figs_april26', 'group_bars');
if ~exist(fig_outdir, 'dir'), mkdir(fig_outdir); end

for t = 1:n_tasks
    task = tasks{t};
    for s = 1:n_sessions
        sess = sessions{s};
        stats_table = stats_tables.(task).(sess);
        
        % =========================================================
        % FIGURE 1: STANDARD PLOT (Corrected Signal Only)
        % =========================================================
        hf1 = figure('Color', 'w', 'Position', [100 100 600 600], ...
            'Name', sprintf('Standard | %s - %s', upper(task), upper(sess)), 'Visible', 'off'); 
        
        for m = 1:1
            measure = measure_names{m};
            subplot(1, 1, m);
            
            R_data = flat_data.(task).(sess).R_matrix(:, :, m);
            mean_R = mean(R_data, 1);
            sem_R = std(R_data, 0, 1) / sqrt(n_subjects);
            p_adj = stats_table.([measure '_p_fwer']);
            
            bar(mean_R, 'FaceColor', [0.5 0.7 0.9]); hold on;
            errorbar(1:n_rois, mean_R, sem_R, sem_R, 'k', 'LineStyle', 'none');
            
            for r = 1:n_rois
                if p_adj(r) < 0.05
                    star = '*'; if p_adj(r) < 0.01, star = '**'; end; if p_adj(r) < 0.001, star = '***'; end
                    text(r, mean_R(r) + sem_R(r) + (plot_y_limit * 0.05), star, ...
                        'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
                end
            end
            
            hold off;
            set(gca, 'XTick', 1:n_rois, 'XTickLabel', rois, 'XTickLabelRotation', 45);
            ylabel('Fisher Z (Corrected Signal)'); title(strrep(measure, '_', ' ')); ylim([0 plot_y_limit]); 
        end
        sgtitle(sprintf('Standard RSA: %s (%s)\n%s', upper(task), upper(sess), strrep(pipeline_name, '_', ' ')), 'FontSize', 14, 'FontWeight', 'bold');

        exportgraphics(hf1, fullfile(fig_outdir, sprintf('group_rsa_standard_%s_%s_%s.jpg', pipeline_name, task, sess)), 'Resolution', 300);
        saveas(hf1, fullfile(fig_outdir, sprintf('group_rsa_standard_%s_%s_%s.svg', pipeline_name, task, sess)), 'svg');
        close(hf1);
        
        % =========================================================
        % FIGURE 2: OBSERVED VALUES BARS & BASELINE LINE
        % =========================================================
        hf2 = figure('Color', 'w', 'Position', [100 100 600 600], ...
            'Name', sprintf('Observed & Baseline | %s - %s', upper(task), upper(sess)), 'Visible', 'off'); 
        
        for m = 1:1
            measure = measure_names{m};
            subplot(1, 1, m);
            
            R_data = flat_data.(task).(sess).R_matrix(:, :, m);
            B_data = flat_data.(task).(sess).B_matrix(:, :, m);
            
            mean_R = mean(R_data, 1);
            mean_B = mean(B_data, 1);
            p_adj = stats_table.([measure '_p_fwer']);
            
            uncorr_data = R_data + B_data; 
            mean_uncorr = mean_R + mean_B;
            sem_uncorr = std(uncorr_data, 0, 1) / sqrt(n_subjects);
            
            hold on;
            b1 = bar(1:n_rois, mean_uncorr, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
            errorbar(1:n_rois, mean_uncorr, sem_uncorr, sem_uncorr, 'k', 'LineStyle', 'none');
            l1 = plot(1:n_rois, mean_B, 'k--', 'LineWidth', 2);
            
            for r = 1:n_rois
                if p_adj(r) < 0.05
                    star = '*'; if p_adj(r) < 0.01, star = '**'; end; if p_adj(r) < 0.001, star = '***'; end
                    text(r, mean_uncorr(r) + sem_uncorr(r) + (plot_y_limit_stacked * 0.05), star, ...
                        'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
                end
            end
            
            hold off;
            set(gca, 'XTick', 1:n_rois, 'XTickLabel', rois, 'XTickLabelRotation', 45);
            ylabel('Fisher Z'); title(strrep(measure, '_', ' ')); 
            ylim([0 plot_y_limit_stacked]); 
            
            if m == 1
                lgd = legend([b1, l1], {'Observed Values', 'Null Baseline'});
                set(lgd, 'Location', 'best', 'FontSize', 10, 'EdgeColor', 'none', 'Color', 'none');
            end
        end
        sgtitle(sprintf('Observed Values & Baseline: %s (%s)\n%s', upper(task), upper(sess), strrep(pipeline_name, '_', ' ')), 'FontSize', 14, 'FontWeight', 'bold');

        exportgraphics(hf2, fullfile(fig_outdir, sprintf('group_rsa_overlay_%s_%s_%s.jpg', pipeline_name, task, sess)), 'Resolution', 300);
        saveas(hf2, fullfile(fig_outdir, sprintf('group_rsa_overlay_%s_%s_%s.svg', pipeline_name, task, sess)), 'svg');
        close(hf2);
    end
end

%% --- 4. Pairwise ROI Visualization (VX vs. RT & VX vs. LT) ---
target_measure_idx = 1; 
measure_name = 'Perceptual_Similarity';

idx_lt = 1; 
idx_rt = 2; 
idx_vx = 3;

comparison_pairs = { [idx_vx, idx_rt], [idx_vx, idx_lt] };
pair_names       = { 'VX_vs_RT',       'VX_vs_LT' };

for t = 1:n_tasks
    cur_task = tasks{t};
    roi_outdir = fullfile(out_dir, 'figs_april26', 'pairwise_comparisons', cur_task);
    if ~exist(roi_outdir, 'dir'), mkdir(roi_outdir); end
    
    for r = 1:n_rois
        roi_name = rois{r};
        for p = 1:length(comparison_pairs)
            curr_pair_indices = comparison_pairs{p};
            curr_pair_label   = pair_names{p};
            
            hf_roi = figure('Color', 'w', ... 
                'Position', [200 200 350 600], ... 
                'Name', sprintf('%s | %s: %s (%s)', pipeline_name, cur_task, roi_name, curr_pair_label), ...
                'Visible', 'off', 'InvertHardcopy', 'off'); 
            
            hold on;
            pair_mean = zeros(1, 2);
            pair_sem  = zeros(1, 2);
            pair_p    = zeros(1, 2);
            pair_tags = {sessions{curr_pair_indices(1)}, sessions{curr_pair_indices(2)}};
                
            for i = 1:2
                s_idx = curr_pair_indices(i);
                sess_name = sessions{s_idx};
                
                roi_data = flat_data.(cur_task).(sess_name).R_matrix(:, r, target_measure_idx);
                
                pair_mean(i) = mean(roi_data, 1);
                pair_sem(i)  = std(roi_data, 0, 1) / sqrt(n_subjects);
                
                % Fallback falls die Variable fehlt
                if isfield(stats_tables.(cur_task).(sess_name), [measure_name '_p_fwer'])
                    pair_p(i)    = stats_tables.(cur_task).(sess_name).([measure_name '_p_fwer'])(r);
                else
                    pair_p(i) = 1; 
                end
            end
                
            b = bar(1:2, pair_mean, 0.6, 'FaceColor', 'flat', 'EdgeColor', [0 0 0], 'LineWidth', 1.5);
            colors = zeros(3, 3);
            colors(idx_lt, :) = [0.60 0.30 0.75]; 
            colors(idx_rt, :) = [0.95 0.70 0.10]; 
            colors(idx_vx, :) = [0.60 0.60 0.60]; 
            
            b.CData(1,:) = colors(curr_pair_indices(1), :);
            b.CData(2,:) = colors(curr_pair_indices(2), :);
                
            errorbar(1:2, pair_mean, pair_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
                
            for i = 1:2
                if pair_p(i) < 0.05
                    star = '*';
                    if pair_p(i) < 0.01, star = '**'; end
                    if pair_p(i) < 0.001, star = '***'; end
                    star_y = pair_mean(i) + pair_sem(i) + (plot_y_limit * 0.05);
                    text(i, star_y, star, ...
                        'HorizontalAlignment', 'center', 'Color', 'k', ...
                        'FontSize', 16, 'FontWeight', 'bold');
                end
            end
                
            set(gca, 'XTick', 1:2, 'XTickLabel', pair_tags, ...
                'Box', 'off', 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1.2, ...
                'Color', 'w', 'XColor', 'k', 'YColor', 'k'); 
                
            ylabel('Fisher Z (Corrected Signal)'); 
            title(sprintf('%s: %s\n%s', strrep(pipeline_name,'_',' '), roi_name, curr_pair_label), 'Interpreter', 'none');
            
            ylim([0 plot_y_limit]); 

            % --- Save ---
            out_filename_png = sprintf('%s_%s_%s_%s_transp.png', pipeline_name, cur_task, roi_name, curr_pair_label);
            outpath_png = fullfile(roi_outdir, out_filename_png);
            exportgraphics(hf_roi, outpath_png, 'Resolution', 300, 'BackgroundColor', 'w');
            exportgraphics(hf_roi, fullfile(roi_outdir, sprintf('%s_%s_%s_%s.eps', pipeline_name, cur_task, roi_name, curr_pair_label)), 'ContentType', 'vector', 'BackgroundColor', 'none');
            saveas(hf_roi, fullfile(roi_outdir, sprintf('%s_%s_%s_%s.svg', pipeline_name, cur_task, roi_name, curr_pair_label)), 'svg');
            
            % Transparency Hack
            try
                [img, ~, ~] = imread(outpath_png);
                isWhite = (img(:,:,1) == 255) & (img(:,:,2) == 255) & (img(:,:,3) == 255);
                alphaMap = uint8(~isWhite * 255);
                imwrite(img, outpath_png, 'Alpha', alphaMap);
            catch
            end
            
            % SVG Bars only
            title(''); ylabel(''); axis off; 
            saveas(hf_roi, fullfile(roi_outdir, sprintf('%s_%s_%s_%s_barsonly.svg', pipeline_name, cur_task, roi_name, curr_pair_label)), 'svg');
            close(hf_roi);
        end
    end
end
fprintf('Figures saved successfully to: %s\n', fullfile(out_dir, 'figs_april26'));
end

% function [stats_tables] = pronet_rsa_statistics_universal(all_subject_results, VP, rois, tasks, out_dir, pipeline_name)
% % Performs session-wise group-level statistics (non-parametric one-sample t-test)
% % Generates STANDARD plots and STACKED BASELINE OVERLAY plots.
% 
% fprintf('\n--- Starting Session-wise Statistical Analysis: %s ---\n', pipeline_name);
% 
% n_subjects = numel(VP);
% n_rois = numel(rois);
% n_tasks = numel(tasks);
% sessions = {'lt', 'rt', 'vx'};
% n_sessions = numel(sessions);
% n_perms = 100000; 
% 
% plot_y_limit = 0.25; 
% plot_y_limit_stacked = 0.25;
% 
% % --- 1. Reshape Data ---
% flat_data = struct(); 
% 
% for t = 1:n_tasks
%     task = tasks{t};
%     for s = 1:n_sessions
%         sess = sessions{s};
%         R_matrix = zeros(n_subjects, n_rois, 3); % Signal (Corrected)
%         B_matrix = zeros(n_subjects, n_rois, 3); % Baseline (Noise Floor)
% 
%         for i = 1:n_subjects
%             subj_id = ['sub' VP{i}];
%             for r = 1:n_rois
%                 roi_name = rois{r};
%                 if isfield(all_subject_results, subj_id) && ...
%                    isfield(all_subject_results.(subj_id), task) && ...
%                    isfield(all_subject_results.(subj_id).(task), roi_name) && ...
%                    isfield(all_subject_results.(subj_id).(task).(roi_name), sess)
% 
%                     dat = all_subject_results.(subj_id).(task).(roi_name).(sess);
% 
%                     % Corrected Signal
%                     R_matrix(i, r, 1) = dat.R_general;
%                     R_matrix(i, r, 2) = dat.R_unique_Acoustic;
%                     R_matrix(i, r, 3) = dat.R_unique_Categorical;
% 
%                     % Null Baseline (Falls vorhanden)
%                     if isfield(dat, 'Baseline_Z_general')
%                         B_matrix(i, r, 1) = dat.Baseline_Z_general;
%                         B_matrix(i, r, 2) = dat.Baseline_Z_unique_Acoustic;
%                         B_matrix(i, r, 3) = dat.Baseline_Z_unique_Categorical;
%                     end
%                 end
%             end
%         end
%         flat_data.(task).(sess).R_matrix = R_matrix;
%         flat_data.(task).(sess).B_matrix = B_matrix;
%     end
% end
% fprintf('Data reshaping complete.\n');
% 
% % --- 2. Perform Statistics (On Corrected Signal) ---
% stats_tables = struct();
% measure_names = {'General', 'Unique_Acoustic', 'Unique_Categorical'};
% 
% for t = 1:n_tasks
%     task = tasks{t};
%     for s = 1:n_sessions
%         sess = sessions{s};
%         fprintf('Calculating stats: Task %s | Session %s\n', upper(task), upper(sess));
% 
%         stats_table = table('RowNames', rois);
%         for m = 1:3 
%             measure = measure_names{m};
%             data_to_test = flat_data.(task).(sess).R_matrix(:, :, m);
% 
%             % Run Permutation Test
%             [t_stats, p_adj, ~, ~] = permuttest(data_to_test, 0, ...
%                                                'nperm', n_perms, ...
%                                                'tail', 'right', ...
%                                                'correct', true, ...
%                                                'verbose', 0);
% 
%             stats_table.([measure '_tstat']) = t_stats';
%             stats_table.([measure '_p_fwer']) = p_adj';
%         end
%         stats_tables.(task).(sess) = stats_table;
%     end
% end
% fprintf('Statistical testing complete.\n');
% 
% %% --- 3. Visualize Results ---
% fig_outdir = fullfile(out_dir, 'figs_april26', 'group_bars');
% if ~exist(fig_outdir, 'dir'), mkdir(fig_outdir); end
% 
% for t = 1:n_tasks
%     task = tasks{t};
%     for s = 1:n_sessions
%         sess = sessions{s};
%         stats_table = stats_tables.(task).(sess);
% 
%         % =========================================================
%         % FIGURE 1: STANDARD PLOT (Corrected Signal Only)
%         % =========================================================
%         hf1 = figure('Color', 'w', 'Position', [100 100 1400 600], ...
%             'Name', sprintf('Standard | %s - %s', upper(task), upper(sess)), 'Visible', 'off'); 
% 
%         for m = 1:3
%             measure = measure_names{m};
%             subplot(1, 3, m);
% 
%             R_data = flat_data.(task).(sess).R_matrix(:, :, m);
%             mean_R = mean(R_data, 1);
%             sem_R = std(R_data, 0, 1) / sqrt(n_subjects);
%             p_adj = stats_table.([measure '_p_fwer']);
% 
%             bar(mean_R, 'FaceColor', [0.5 0.7 0.9]); hold on;
%             errorbar(1:n_rois, mean_R, sem_R, sem_R, 'k', 'LineStyle', 'none');
% 
%             for r = 1:n_rois
%                 if p_adj(r) < 0.05
%                     star = '*'; if p_adj(r) < 0.01, star = '**'; end; if p_adj(r) < 0.001, star = '***'; end
%                     text(r, mean_R(r) + sem_R(r) + (plot_y_limit * 0.05), star, ...
%                         'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
%                 end
%             end
% 
%             hold off;
%             set(gca, 'XTick', 1:n_rois, 'XTickLabel', rois, 'XTickLabelRotation', 45);
%             ylabel('Fisher Z (Corrected Signal)'); title(strrep(measure, '_', ' ')); ylim([0 plot_y_limit]); 
%         end
%         sgtitle(sprintf('Standard RSA: %s (%s) | %s', upper(task), upper(sess), strrep(pipeline_name, '_', ' ')), 'FontSize', 16, 'FontWeight', 'bold');
% 
%         % Save Fig 1
%         exportgraphics(hf1, fullfile(fig_outdir, sprintf('group_rsa_standard_%s_%s_%s.jpg', pipeline_name, task, sess)), 'Resolution', 300);
%         saveas(hf1, fullfile(fig_outdir, sprintf('group_rsa_standard_%s_%s_%s.svg', pipeline_name, task, sess)), 'svg');
%         close(hf1);
% 
%         % =========================================================
%         % FIGURE 2: OBSERVED VALUES BARS & BASELINE LINE
%         % =========================================================
%         hf2 = figure('Color', 'w', 'Position', [100 100 1400 600], ...
%             'Name', sprintf('Observed & Baseline | %s - %s', upper(task), upper(sess)), 'Visible', 'off'); 
% 
%         for m = 1:3
%             measure = measure_names{m};
%             subplot(1, 3, m);
% 
%             R_data = flat_data.(task).(sess).R_matrix(:, :, m);
%             B_data = flat_data.(task).(sess).B_matrix(:, :, m);
% 
%             mean_R = mean(R_data, 1);
%             mean_B = mean(B_data, 1);
%             p_adj = stats_table.([measure '_p_fwer']);
% 
%             % Uncorrected data = Corrected Signal + Baseline
%             uncorr_data = R_data + B_data; 
%             mean_uncorr = mean_R + mean_B;
%             sem_uncorr = std(uncorr_data, 0, 1) / sqrt(n_subjects);
% 
%             hold on;
% 
%             % 1. Die beobachteten Werte (inkl. Rauschen) als Balken zeichnen
%             b1 = bar(1:n_rois, mean_uncorr, 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
% 
%             % 2. Fehlerbalken für die Gesamtdaten hinzufügen
%             errorbar(1:n_rois, mean_uncorr, sem_uncorr, sem_uncorr, 'k', 'LineStyle', 'none');
% 
%             % 3. Das Rauschen (Baseline) als gestrichelte Linie darüberlegen
%             l1 = plot(1:n_rois, mean_B, 'k--', 'LineWidth', 2);
% 
%             % Signifikanz-Sterne
%             for r = 1:n_rois
%                 if p_adj(r) < 0.05
%                     star = '*'; if p_adj(r) < 0.01, star = '**'; end; if p_adj(r) < 0.001, star = '***'; end
%                     text(r, mean_uncorr(r) + sem_uncorr(r) + (plot_y_limit_stacked * 0.05), star, ...
%                         'HorizontalAlignment', 'center', 'Color', 'r', 'FontSize', 14);
%                 end
%             end
% 
%             hold off;
%             set(gca, 'XTick', 1:n_rois, 'XTickLabel', rois, 'XTickLabelRotation', 45);
%             ylabel('Fisher Z'); title(strrep(measure, '_', ' ')); 
%             ylim([0 plot_y_limit_stacked]); 
% 
%             % Legende (Frei schwebend im dritten Subplot)
%             if m == 3
%                 lgd = legend([b1, l1], {'Observed Values (incl. Noise)', 'Null Baseline (Noise)'});
%                 set(lgd, 'Position', [0.82 0.82 0.15 0.08], ...
%                          'Units', 'normalized', ...
%                          'FontSize', 10, ...
%                          'EdgeColor', 'none', ...
%                          'Color', 'none');
%             end
%         end
%         sgtitle(sprintf('Observed Values & Baseline: %s (%s) | %s', upper(task), upper(sess), strrep(pipeline_name, '_', ' ')), 'FontSize', 16, 'FontWeight', 'bold');
% 
%         % Speichern
%         exportgraphics(hf2, fullfile(fig_outdir, sprintf('group_rsa_overlay_%s_%s_%s.jpg', pipeline_name, task, sess)), 'Resolution', 300);
%         saveas(hf2, fullfile(fig_outdir, sprintf('group_rsa_overlay_%s_%s_%s.svg', pipeline_name, task, sess)), 'svg');
%         close(hf2);
%     end
% end
% 
% %% --- 4. Pairwise ROI Visualization (VX vs. RT & VX vs. LT) ---
% target_measure_idx = 1; % 1 = General Model
% measure_name = 'General';
% 
% idx_lt = 1; 
% idx_rt = 2; 
% idx_vx = 3;
% 
% comparison_pairs = { [idx_vx, idx_rt], [idx_vx, idx_lt] };
% pair_names       = { 'VX_vs_RT',       'VX_vs_LT' };
% 
% for t = 1:n_tasks
%     cur_task = tasks{t};
%     roi_outdir = fullfile(out_dir, 'figs_april26', 'pairwise_comparisons', cur_task);
%     if ~exist(roi_outdir, 'dir'), mkdir(roi_outdir); end
% 
%     for r = 1:n_rois
%         roi_name = rois{r};
%         for p = 1:length(comparison_pairs)
%             curr_pair_indices = comparison_pairs{p};
%             curr_pair_label   = pair_names{p};
% 
%             hf_roi = figure('Color', 'w', ... 
%                 'Position', [200 200 350 600], ... 
%                 'Name', sprintf('%s | %s: %s (%s)', pipeline_name, cur_task, roi_name, curr_pair_label), ...
%                 'Visible', 'off', 'InvertHardcopy', 'off'); 
% 
%             hold on;
%             pair_mean = zeros(1, 2);
%             pair_sem  = zeros(1, 2);
%             pair_p    = zeros(1, 2);
%             pair_tags = {sessions{curr_pair_indices(1)}, sessions{curr_pair_indices(2)}};
% 
%             for i = 1:2
%                 s_idx = curr_pair_indices(i);
%                 sess_name = sessions{s_idx};
% 
% 
%                 roi_data = flat_data.(cur_task).(sess_name).R_matrix(:, r, target_measure_idx);
% 
%                 pair_mean(i) = mean(roi_data, 1);
%                 pair_sem(i)  = std(roi_data, 0, 1) / sqrt(n_subjects);
%                 pair_p(i)    = stats_tables.(cur_task).(sess_name).([measure_name '_p_fwer'])(r);
%             end
% 
%             b = bar(1:2, pair_mean, 0.6, 'FaceColor', 'flat', 'EdgeColor', [0 0 0], 'LineWidth', 1.5);
%             colors = zeros(3, 3);
%             colors(idx_lt, :) = [0.60 0.30 0.75]; 
%             colors(idx_rt, :) = [0.95 0.70 0.10]; 
%             colors(idx_vx, :) = [0.60 0.60 0.60]; 
% 
%             b.CData(1,:) = colors(curr_pair_indices(1), :);
%             b.CData(2,:) = colors(curr_pair_indices(2), :);
% 
%             errorbar(1:2, pair_mean, pair_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 10);
% 
%             for i = 1:2
%                 if pair_p(i) < 0.05
%                     star = '*';
%                     if pair_p(i) < 0.01, star = '**'; end
%                     if pair_p(i) < 0.001, star = '***'; end
%                     star_y = pair_mean(i) + pair_sem(i) + (plot_y_limit * 0.05);
%                     text(i, star_y, star, ...
%                         'HorizontalAlignment', 'center', 'Color', 'k', ...
%                         'FontSize', 16, 'FontWeight', 'bold');
%                 end
%             end
% 
%             set(gca, 'XTick', 1:2, 'XTickLabel', pair_tags, ...
%                 'Box', 'off', 'TickDir', 'out', 'FontSize', 12, 'LineWidth', 1.2, ...
%                 'Color', 'w', 'XColor', 'k', 'YColor', 'k'); 
% 
%             ylabel('Fisher Z (Corrected Signal)'); 
%             title(sprintf('%s: %s\n%s', strrep(pipeline_name,'_',' '), roi_name, curr_pair_label), 'Interpreter', 'none');
% 
%             % Apply the fixed Y-Limit
%             ylim([0 plot_y_limit]); 
% 
%             % --- Save ---
%             out_filename_png = sprintf('%s_%s_%s_%s_transp.png', pipeline_name, cur_task, roi_name, curr_pair_label);
%             outpath_png = fullfile(roi_outdir, out_filename_png);
%             exportgraphics(hf_roi, outpath_png, 'Resolution', 300, 'BackgroundColor', 'w');
%             exportgraphics(hf_roi, fullfile(roi_outdir, sprintf('%s_%s_%s_%s.eps', pipeline_name, cur_task, roi_name, curr_pair_label)), 'ContentType', 'vector', 'BackgroundColor', 'none');
%             saveas(hf_roi, fullfile(roi_outdir, sprintf('%s_%s_%s_%s.svg', pipeline_name, cur_task, roi_name, curr_pair_label)), 'svg');
% 
%             % Transparency Hack
%             [img, ~, ~] = imread(outpath_png);
%             isWhite = (img(:,:,1) == 255) & (img(:,:,2) == 255) & (img(:,:,3) == 255);
%             alphaMap = uint8(~isWhite * 255);
%             imwrite(img, outpath_png, 'Alpha', alphaMap);
% 
%             % SVG Bars only
%             title(''); ylabel(''); axis off; 
%             saveas(hf_roi, fullfile(roi_outdir, sprintf('%s_%s_%s_%s_barsonly.svg', pipeline_name, cur_task, roi_name, curr_pair_label)), 'svg');
%             close(hf_roi);
%         end
%     end
% end
% fprintf('Figures saved successfully (incl. EPS & SVG for Adobe Illustrator) to: %s\n', fullfile(out_dir, 'figs_april26'));
% end