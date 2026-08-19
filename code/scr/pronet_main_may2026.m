%% pronet_main_may2026.m
% Master Script to run RSA Pipeline (4 selected versions)
% Includes modeling, statistics, ANOVA, Noise Ceilings.
% written by Antonia Ceric, 2026

clear; clc;
addpath /mnt/beegfs/users/antonia.ceric/git/pronet_masterthesis/func;
addpath /mnt/beegfs/users/antonia.ceric/git/pronet/pronet_masterthesis/func/pretest_modelling/;
run('pronet_startup.m'); 

global DATA_PATH;

% --- 1. GLOBAL CONFIGURATION ---
VP = {'02','03','04','05','06','07','08','09','10',...
      '11','12','13','15','16','17','18','19','20',...
      '21','22','23','24','25','26','27','28','29','31','32'};

TASKS = {'phoneme', 'intonat'};
ROIS  = {'lIFG', 'lPMC', 'lPAC', 'lpSTS', 'laSTS', ...
         'rIFG', 'rPMC', 'rPAC', 'rpSTS', 'raSTS'};
SESSIONS = {'lt', 'rt', 'vx'}; 

% Hardcoded Exclusions
excl.phoneme = {'08', '05'};
excl.intonat = {'13', '05'};

% Noise Ceiling Config
METRIC = 'R_general'; 
plot_y_limit = 1.3; % Universal Y-Axis for plots

% =========================================================================
% --- 2. PIPELINE CONFIGURATIONS ---
% =========================================================================
clear pipelines;

% ---------------------------------------------------------
% 4. avg Rsq (183 reg) [Late Fusion]
% ---------------------------------------------------------
p = 1;
pipelines(p).name       = 'LateFusion_183reg_CatOnly';
pipelines(p).logic_type = 'late_fusion'; 
pipelines(p).ldc_dir    = fullfile(DATA_PATH,'rsa_crossnobis_Speakerwise_Feb2026_110326');
pipelines(p).ldc_fname  = 'sub%s_task-%s_roi-%s_LDC_Speakerwise_Feb2026.mat';
pipelines(p).ldc_type   = 'avg_m1_f2'; 
pipelines(p).model_dir  = fullfile(DATA_PATH, 'model_rdms_pretest_Feb2026_speakerwisefit');
pipelines(p).model_fname= 'sub-%s_modelRDM_pretest_Speakerwise.mat';
pipelines(p).model_type = 'both'; 
pipelines(p).res_file   = fullfile(DATA_PATH,'rsa_results_april26_LateFusion_183reg_062026/rsa_whitened_perm_results_183reg_LateFusion_183reg.mat');
pipelines(p).out_dir    = fullfile(DATA_PATH,'rsa_results_april26_LateFusion_183reg_062026');

% ---------------------------------------------------------
% 5. avg spkwise (183 reg) [Early Fusion]
% ---------------------------------------------------------
p = 2;
pipelines(p).name       = 'EarlyFusion_183reg_CatOnly';
pipelines(p).logic_type = 'early_fusion'; 
pipelines(p).ldc_dir    = fullfile(DATA_PATH,'rsa_crossnobis_Speakerwise_Feb2026_110326');
pipelines(p).ldc_fname  = 'sub%s_task-%s_roi-%s_LDC_Speakerwise_Feb2026.mat';
pipelines(p).ldc_type   = 'avg_m1_f2'; 
pipelines(p).model_dir  = fullfile(DATA_PATH, 'model_rdms_pretest_Feb2026_speakerwisefit');
pipelines(p).model_fname= 'sub-%s_modelRDM_pretest_Speakerwise.mat';
pipelines(p).model_type = 'compute_avg'; 
pipelines(p).res_file   = fullfile(DATA_PATH,'rsa_results_april26_EarlyFusion_183reg_062026/rsa_whitened_perm_results_Averaged_183reg.mat');
pipelines(p).out_dir    = fullfile(DATA_PATH,'rsa_results_april26_EarlyFusion_183reg_062026');

% ---------------------------------------------------------
% 8. averaged Rsq (93 reg, bigblk) [Late Fusion]
% ---------------------------------------------------------
p = 3;
pipelines(p).name       = 'LateFusion_93reg_bigblk_CatOnly';
pipelines(p).logic_type = 'late_fusion';
pipelines(p).ldc_dir    = fullfile(DATA_PATH,'rsa_crossnobis_bigblocks_March2026');
pipelines(p).ldc_fname  = 'sub%s_task-%s_roi-%s_LDC_Speakerwise_bigblocks_March2026.mat';
pipelines(p).ldc_type   = 'avg_m1_f2'; 
pipelines(p).model_dir  = fullfile(DATA_PATH, 'model_rdms_pretest_Feb2026_speakerwisefit');
pipelines(p).model_fname= 'sub-%s_modelRDM_pretest_Speakerwise.mat';
pipelines(p).model_type = 'both';
pipelines(p).res_file   = fullfile(DATA_PATH,'rsa_results_april26_LateFusion_93reg_062026/rsa_whitened_perm_results_LateFusion.mat');
pipelines(p).out_dir    = fullfile(DATA_PATH,'rsa_results_april26_LateFusion_93reg_062026');

% ---------------------------------------------------------
% 9. avg spkwise (93 reg, bigblk) [Early Fusion]
% ---------------------------------------------------------
p = 4;
pipelines(p).name       = 'EarlyFusion_93reg_bigblk_CatOnly';
pipelines(p).logic_type = 'early_fusion';
pipelines(p).ldc_dir    = fullfile(DATA_PATH,'rsa_crossnobis_bigblocks_March2026');
pipelines(p).ldc_fname  = 'sub%s_task-%s_roi-%s_LDC_Speakerwise_bigblocks_March2026.mat';
pipelines(p).ldc_type   = 'avg_m1_f2'; 
pipelines(p).model_dir  = fullfile(DATA_PATH, 'model_rdms_pretest_Feb2026_speakerwisefit');
pipelines(p).model_fname= 'sub-%s_modelRDM_pretest_Speakerwise.mat';
pipelines(p).model_type = 'compute_avg'; 
pipelines(p).res_file   = fullfile(DATA_PATH,'rsa_results_april26_EarlyFusion_93reg_062026/rsa_whitened_perm_results_EarlyFusion.mat');
pipelines(p).out_dir    = fullfile(DATA_PATH,'rsa_results_april26_EarlyFusion_93reg_062026');


% =========================================================================
% --- 3. MASTER EXECUTION LOOP (MODELLING & STATS) ---
% =========================================================================

for p_idx = 1:numel(pipelines)
    cfg = pipelines(p_idx);
    
    fprintf('\n=================================================================\n');
    fprintf('STARTING PIPELINE %d/%d: %s\n', p_idx, numel(pipelines), cfg.name);
    fprintf('=================================================================\n');
    
    % Directory Check and Creation
    if ~exist(cfg.out_dir, 'dir')
        mkdir(cfg.out_dir);
    end
    
    [res_dir, ~, ~] = fileparts(cfg.res_file);
    if ~exist(res_dir, 'dir')
        mkdir(res_dir); 
    end
    
    all_subject_results = struct();
    
    try
        for i = 1:numel(VP)
            subj = VP{i};
            subj_key = ['sub' subj];
            
            % --- A. LOAD MODEL ---
            model_file = fullfile(cfg.model_dir, sprintf(cfg.model_fname, subj));
            if ~exist(model_file, 'file')
                fprintf(2, 'Missing Model for %s\n', subj);
                continue; 
            end
            tmp_mod = load(model_file);
            
            for t = 1:numel(TASKS)
                task = TASKS{t};
                if ismember(subj, excl.(task)), continue; end
                
                model_name = [task '_pretest'];
                if ~isfield(tmp_mod.modelRDMs, model_name), continue; end
                
                % --- EXTRACT OR COMPUTE MODEL MATRICES ---
                if strcmp(cfg.model_type, 'both')
                    m_cat_m1 = tmp_mod.modelRDMs.(model_name).m1.category;
                    m_aco_m1 = tmp_mod.modelRDMs.(model_name).m1.acoustic;
                    m_cat_f2 = tmp_mod.modelRDMs.(model_name).f2.category;
                    m_aco_f2 = tmp_mod.modelRDMs.(model_name).f2.acoustic;
                elseif strcmp(cfg.model_type, 'compute_avg')
                    m_cat_m1 = tmp_mod.modelRDMs.(model_name).m1.category;
                    m_cat_f2 = tmp_mod.modelRDMs.(model_name).f2.category;
                    m_aco_m1 = tmp_mod.modelRDMs.(model_name).m1.acoustic;
                    m_aco_f2 = tmp_mod.modelRDMs.(model_name).f2.acoustic;
                    m_cat = (m_cat_m1 + m_cat_f2) / 2;
                    m_aco = (m_aco_m1 + m_aco_f2) / 2;
                else
                    m_cat = tmp_mod.modelRDMs.(model_name).(cfg.model_type).category;
                    m_aco = tmp_mod.modelRDMs.(model_name).(cfg.model_type).acoustic;
                end
                
                for r = 1:numel(ROIS)
                    roi = ROIS{r};
                    
                    % --- B. LOAD LDC ---
                    ldc_file = fullfile(cfg.ldc_dir, sprintf(cfg.ldc_fname, subj, task, roi));
                    if ~exist(ldc_file, 'file'), continue; end
                    tmp_ldc = load(ldc_file, 'Sub');
                    
                    if ~isfield(tmp_ldc.Sub, task) || ~isfield(tmp_ldc.Sub.(task), roi), continue; end
                    
                    for s = 1:numel(SESSIONS)
                        sess = SESSIONS{s};
                        if ~isfield(tmp_ldc.Sub.(task).(roi), sess), continue; end
                        
                        sess_data = tmp_ldc.Sub.(task).(roi).(sess);
                        
                        % --- C. CORE MATHEMATICS (FUSION LOGIC) ---
                        if strcmp(cfg.logic_type, 'single')
                            if strcmp(cfg.ldc_type, 'agg')
                                if ~isfield(sess_data, 'RDM'), continue; end
                                n_rdm = sess_data.RDM;
                            else
                                if ~isfield(sess_data, cfg.ldc_type), continue; end
                                n_rdm = sess_data.(cfg.ldc_type).RDM;
                            end
                            
                            if isempty(n_rdm), continue; end
                            if isvector(n_rdm), n_rdm = squareform(n_rdm); end
                            
                            % m_aco entfernt
                            [R_cor, ~, ~] = pronet_single_perm_test_pretest(n_rdm, m_cat);
                            all_subject_results.(subj_key).(task).(roi).(sess) = R_cor;
                            
                        elseif strcmp(cfg.logic_type, 'early_fusion')
                            if ~isfield(sess_data, 'm1') || ~isfield(sess_data, 'f2'), continue; end
                            
                            rdm_m1 = sess_data.m1.RDM;
                            rdm_f2 = sess_data.f2.RDM;
                            if isempty(rdm_m1) || isempty(rdm_f2), continue; end
                            
                            if isvector(rdm_m1), rdm_m1 = squareform(rdm_m1); end
                            if isvector(rdm_f2), rdm_f2 = squareform(rdm_f2); end
                            
                            % Average Neural RDMs BEFORE fitting
                            neural_avg = (rdm_m1 + rdm_f2) / 2;
                            
                            % m_aco entfernt
                            [R_cor, ~, ~] = pronet_single_perm_test_pretest(neural_avg, m_cat);
                            
                            all_subject_results.(subj_key).(task).(roi).(sess) = R_cor;
                            
                        elseif strcmp(cfg.logic_type, 'late_fusion')
                            if ~isfield(sess_data, 'm1') || ~isfield(sess_data, 'f2'), continue; end
                            
                            rdm_m1 = sess_data.m1.RDM;
                            rdm_f2 = sess_data.f2.RDM;
                            if isempty(rdm_m1) || isempty(rdm_f2), continue; end
                            
                            if isvector(rdm_m1), rdm_m1 = squareform(rdm_m1); end
                            if isvector(rdm_f2), rdm_f2 = squareform(rdm_f2); end
                            
                            % Fit M1 and F2 SEPARATELY (m_aco_m1 und m_aco_f2 entfernt)
                            [R_cor_m1, ~, ~] = pronet_single_perm_test_pretest(rdm_m1, m_cat_m1);
                            [R_cor_f2, ~, ~] = pronet_single_perm_test_pretest(rdm_f2, m_cat_f2);
                            
                            % Average the resulting Fisher-Z transformed R values
                            avg_res.R_general = (R_cor_m1.R_general + R_cor_f2.R_general) / 2;
                            
                            % Average baselines
                            avg_res.Baseline_Z_general = (R_cor_m1.Baseline_Z_general + R_cor_f2.Baseline_Z_general) / 2;
                            
                            all_subject_results.(subj_key).(task).(roi).(sess) = avg_res;
                        end
                    end % Session
                end % ROI
            end % Task
        end % Subject
        
        % --- D. SAVE RAW RESULTS ---
        save(cfg.res_file, 'all_subject_results');
        fprintf(' Finished and saved raw results: %s\n', cfg.res_file);
        
        % --- E. RUN STATISTICS & ANOVA ---
        fprintf('Running Statistics and Generating Plots for %s...\n', cfg.name);
        
        try
            pronet_rsa_statistics_universal(all_subject_results, VP, ROIS, TASKS, cfg.out_dir, cfg.name);
        catch ME
            fprintf(2, 'Stats/Plotting failed for %s: %s\n', cfg.name, ME.message);
        end
        
        try
            pronet_rsa_anova_may26(all_subject_results, VP, cfg.out_dir, cfg.name);
        catch ME
            fprintf(2, 'ANOVA failed for %s: %s\n', cfg.name, ME.message);
        end
        
    catch ME
        fprintf(2, 'CRITICAL ERROR IN PIPELINE %s: %s\n', cfg.name, ME.message);
    end
end


% =========================================================================
% --- 4. NOISE CEILING COMPUTATION LOOP ---
% =========================================================================
fprintf('\n======================================================\n');
fprintf('NOISE CEILING COMPUTATION STARTED \n');
fprintf('======================================================\n');

for p_idx = 1:numel(pipelines)
    cfg = pipelines(p_idx);
    fprintf('\n===================================================================\n');
    fprintf('=== COMPUTING CEILINGS [%d/%d]: %s ===\n', p_idx, numel(pipelines), upper(cfg.name));
    fprintf('===================================================================\n');

    % Setup Output Directory for Noise Ceilings
    out_nc_dir = fullfile(cfg.out_dir, 'noise_ceilings');
    if ~exist(out_nc_dir, 'dir'), mkdir(out_nc_dir); end

    % Load Result File
    if ~exist(cfg.res_file, 'file')
        warning('Result file missing, skipping pipeline: %s', cfg.res_file);
        continue;
    end
    tmp_res = load(cfg.res_file);
    res_fields = fieldnames(tmp_res);
    model_res = tmp_res.(res_fields{1}); 

    results_table = {}; 

    for t = 1:numel(TASKS)
        task = TASKS{t};
        curr_excl = excl.(task);

        for r = 1:numel(ROIS)
            roi = ROIS{r};
            stacks = struct('vx', [], 'lt', [], 'rt', []);
            stacks_alt = struct('vx', [], 'lt', [], 'rt', []); % Used only for Late Fusion
            n_subs = struct('vx', 0, 'lt', 0, 'rt', 0);

            % 4A. LOAD AND PREPARE BRAIN DATA
            for i = 1:numel(VP)
                subj = VP{i};
                if ismember(subj, curr_excl), continue; end

                fname = sprintf(cfg.ldc_fname, subj, task, roi);
                fpath = fullfile(cfg.ldc_dir, fname);
                if ~exist(fpath, 'file'), continue; end

                try
                    tmp = load(fpath, 'Sub'); 
                    if isfield(tmp.Sub, task) && isfield(tmp.Sub.(task), roi)
                        for s = 1:numel(SESSIONS)
                            sess = SESSIONS{s};
                            rdm = [];

                            if strcmp(cfg.logic_type, 'late_fusion')
                                % Late Fusion: Keep M1 and F2 completely separate
                                if isfield(tmp.Sub.(task).(roi).(sess), 'm1') && isfield(tmp.Sub.(task).(roi).(sess), 'f2')
                                    rdm_m1 = tmp.Sub.(task).(roi).(sess).m1.RDM;
                                    rdm_f2 = tmp.Sub.(task).(roi).(sess).f2.RDM;
                                    if isvector(rdm_m1), rdm_m1 = squareform(rdm_m1); end
                                    if isvector(rdm_f2), rdm_f2 = squareform(rdm_f2); end

                                    if ~any(isnan(rdm_m1(:))) && ~any(isnan(rdm_f2(:)))
                                        if isempty(stacks.(sess))
                                            stacks.(sess) = rdm_m1;
                                            stacks_alt.(sess) = rdm_f2;
                                        else
                                            stacks.(sess) = cat(3, stacks.(sess), rdm_m1);
                                            stacks_alt.(sess) = cat(3, stacks_alt.(sess), rdm_f2);
                                        end
                                        n_subs.(sess) = n_subs.(sess) + 1;
                                    end
                                end
                            else
                                % Standard extraction (Early Fusion or Single)
                                switch cfg.ldc_type
                                    case 'm1'
                                        if isfield(tmp.Sub.(task).(roi).(sess), 'm1'), rdm = tmp.Sub.(task).(roi).(sess).m1.RDM; end
                                    case 'f2'
                                        if isfield(tmp.Sub.(task).(roi).(sess), 'f2'), rdm = tmp.Sub.(task).(roi).(sess).f2.RDM; end
                                    case 'agg'
                                        if isfield(tmp.Sub.(task).(roi).(sess), 'RDM'), rdm = tmp.Sub.(task).(roi).(sess).RDM; end
                                    case 'avg_m1_f2' % Early Fusion averages the brains immediately
                                        if isfield(tmp.Sub.(task).(roi).(sess), 'm1') && isfield(tmp.Sub.(task).(roi).(sess), 'f2')
                                            rdm_m1 = tmp.Sub.(task).(roi).(sess).m1.RDM;
                                            rdm_f2 = tmp.Sub.(task).(roi).(sess).f2.RDM;
                                            rdm = (rdm_m1 + rdm_f2) / 2;
                                        end
                                end

                                if ~isempty(rdm)
                                    if isvector(rdm), rdm = squareform(rdm); end
                                    if ~any(isnan(rdm(:)))
                                        if isempty(stacks.(sess)), stacks.(sess) = rdm;
                                        else, stacks.(sess) = cat(3, stacks.(sess), rdm); end
                                        n_subs.(sess) = n_subs.(sess) + 1;
                                    end
                                end
                            end
                        end
                    end
                catch
                end
            end

            % 4B. CALCULATE BOUNDS AND RETRIEVE MODEL R2
            for s = 1:numel(SESSIONS)
                sess = SESSIONS{s};
                valid_n = n_subs.(sess);
                if valid_n >= 3

                    % Calculate Fisher Z Bounds
                    if strcmp(cfg.logic_type, 'late_fusion')
                        % Late Fusion: Calculate bounds for M1 and F2 separately, then average them
                        [lb_m1, ub_m1] = calculate_bounds_fisher_Z(stacks.(sess));
                        [lb_f2, ub_f2] = calculate_bounds_fisher_Z(stacks_alt.(sess));
                        lb_Z = (lb_m1 + lb_f2) / 2;
                        ub_Z = (ub_m1 + ub_f2) / 2;
                    else
                        % Standard calculation (Early Fusion uses the pre-averaged stack)
                        [lb_Z, ub_Z] = calculate_bounds_fisher_Z(stacks.(sess));
                    end

                    model_Zs = [];
                    sub_keys = fieldnames(model_res);
                    for i = 1:numel(sub_keys)
                        s_idx = strrep(sub_keys{i},'sub','');
                        if ismember(s_idx, curr_excl), continue; end 
                        try
                            % Direct extraction (res_path removed for compatibility with new Modeling loop!)
                            val = model_res.(sub_keys{i}).(task).(roi).(sess).(METRIC);
                            if ~isnan(val), model_Zs(end+1) = val; end
                        catch
                        end
                    end

                    real_mean = mean(model_Zs); 
                    real_sem = std(model_Zs) / sqrt(numel(model_Zs));
                    fprintf('[%s] %-7s | %-5s | %s | N=%d | Mod: %.4f | Ceil: [%.4f - %.4f]\n', ...
                        cfg.name, task, roi, upper(sess), valid_n, real_mean, lb_Z, ub_Z);
                    results_table(end+1,:) = {task, roi, sess, valid_n, real_mean, real_sem, lb_Z, ub_Z};
                end
            end
        end 
    end 

    % 4C. SAVE AND PLOT
    if ~isempty(results_table)
        T_out = cell2table(results_table, 'VariableNames', {'Task', 'ROI', 'Session', 'N', 'Model_FisherZ_Mean', 'Model_FisherZ_SEM', 'Lower_FisherZ', 'Upper_FisherZ'});
        save(fullfile(out_nc_dir, sprintf('noise_ceiling_results_%s.mat', cfg.name)), 'T_out');
        writetable(T_out, fullfile(out_nc_dir, sprintf('noise_ceiling_results_%s.csv', cfg.name)));

        generate_plots(T_out, TASKS, ROIS, SESSIONS, METRIC, out_nc_dir, cfg.name, plot_y_limit);
    else
        fprintf('No valid data found for %s\n', cfg.name);
    end
end
fprintf('\n===================================================================\n');
fprintf('ALL PIPELINES COMPLETE! \n');
fprintf('===================================================================\n');


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function generate_plots(T, tasks, rois, sessions, metric_name, out_dir, prefix, y_limit)
    roi_colors = lines(numel(rois)); 
    markers = struct('vx', 'o', 'lt', '^', 'rt', 's');
    x_offset_map = struct('vx', -0.2, 'lt', 0, 'rt', 0.2); 
    box_width = 0.15; 
    
    for t = 1:numel(tasks)
        task = tasks{t};
        fig_name = sprintf('Ceiling Check: %s', task);
        hf = figure('Name', fig_name, 'Color', 'w', 'Position', [50 50 1400 600], 'Visible', 'off');
        hold on; h_leg = [];
        
        for r = 1:numel(rois)
            roi = rois{r};
            this_col = roi_colors(r, :);
            
            for s = 1:numel(sessions)
                sess = sessions{s};
                x_pos = r + x_offset_map.(sess);
                
                row = T(strcmp(T.Task, task) & strcmp(T.ROI, roi) & strcmp(T.Session, sess), :);
                if isempty(row), continue; end
                
                lb = row.Lower_FisherZ; ub = row.Upper_FisherZ;
                mod_mean = row.Model_FisherZ_Mean; mod_sem = row.Model_FisherZ_SEM;
                
                % Noise Ceiling Box
                patch([x_pos-box_width/2, x_pos+box_width/2, x_pos+box_width/2, x_pos-box_width/2], ...
                      [lb, lb, ub, ub], [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
                  
                % Real Model Fit
                h = errorbar(x_pos, mod_mean, mod_sem, 'Color', this_col, 'LineStyle', 'none', ...
                    'LineWidth', 1.5, 'Marker', markers.(sess), 'MarkerSize', 8, ...
                    'MarkerFaceColor', this_col, 'CapSize', 0);
                
                if r == 1, h_leg(end+1) = h; end
            end
        end
        
        title_str = sprintf('%s Fit vs. Ceiling: %s Task', strrep(prefix,'_',' '), upper(task));
        title(title_str, 'FontSize', 14);
        ylabel('Fisher Z (corrected)', 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Region of Interest (ROI)', 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'XTick', 1:numel(rois), 'XTickLabel', rois); grid on;
        
        % Automatically map the correct label to the order the sessions are processed in
        leg_labels = cell(1, numel(sessions));
        for s = 1:numel(sessions)
            if strcmp(sessions{s}, 'vx'), leg_labels{s} = 'Vertex (Circle)';
            elseif strcmp(sessions{s}, 'lt'), leg_labels{s} = 'Left TMS (Triangle)';
            elseif strcmp(sessions{s}, 'rt'), leg_labels{s} = 'Right TMS (Square)';
            end
        end
        legend(h_leg, leg_labels, 'Location', 'bestoutside', 'FontSize', 11);
                
        xlim([0.5, numel(rois) + 0.5]);
        
        % Fixed Y-Axis for all plots
        ylim([0 y_limit]); 
        
        % Save standard image
        fname = sprintf('Fig_NoiseCeiling_%s_%s.png', prefix, task);
        exportgraphics(hf, fullfile(out_dir, fname), 'Resolution', 300);
        
        % Save Vector Formats for Adobe Illustrator
        exportgraphics(hf, fullfile(out_dir, sprintf('Fig_NoiseCeiling_%s_%s.eps', prefix, task)), 'ContentType', 'vector');
        saveas(hf, fullfile(out_dir, sprintf('Fig_NoiseCeiling_%s_%s.svg', prefix, task)), 'svg');
        
        close(hf);
    end
end

function [lower_bound_Z, upper_bound_Z] = calculate_bounds_fisher_Z(neural_RDMs)
    if ndims(neural_RDMs) ~= 3 || size(neural_RDMs, 3) < 3
        lower_bound_Z = NaN; upper_bound_Z = NaN; return;
    end

    [n_rows, ~, n_subs] = size(neural_RDMs);
    n_conds = n_rows; 

    % Vectorize
    tmp_vec = get_vec(neural_RDMs(:,:,1));
    n_pairs = length(tmp_vec);
    
    rdm_vecs = zeros(n_pairs, n_subs);
    for s = 1:n_subs
        rdm_vecs(:,s) = get_vec(neural_RDMs(:,:,s));
    end
    rdm_vecs = double(rdm_vecs);

    % Whitening Matrix 'W'
    c_mat = zeros(n_pairs, n_conds);
    row_idx = 1;
    for i = 1:n_conds-1
        for j = i+1:n_conds
            c_mat(row_idx, i) = 1; c_mat(row_idx, j) = -1;
            row_idx = row_idx + 1;
        end
    end
    
    sigma_k = eye(n_conds); 
    v = c_mat * sigma_k * c_mat';
    V = v .* v;

    [K, L_matrix] = eig(V);      
    L = diag(L_matrix);         
    inv_l = 1 ./ sqrt(abs(L));      
    inv_l(~isfinite(inv_l)) = 0; 
    W = K * diag(inv_l) * K';    

    % Apply whitening
    rdm_vecs_w = W * rdm_vecs;
    M_intercept_w = W * ones(n_pairs, 1);

    z_lower = zeros(n_subs, 1);
    z_upper = zeros(n_subs, 1);
    
    mean_all_w = mean(rdm_vecs_w, 2);
    TSS_upper = sum(mean_all_w.^2); 
    
    for s = 1:n_subs
        subj_w = rdm_vecs_w(:, s);
        M_subj = [subj_w, M_intercept_w];
        
        % UPPER BOUND
        beta_up = lsqnonneg(M_subj, mean_all_w); 
        RSS_up = sum((mean_all_w - (M_subj * beta_up)).^2);
        r2_up_val = max(1 - (RSS_up / TSS_upper), 0);
        z_upper(s) = atanh(sqrt(r2_up_val)); 
        
        % LOWER BOUND (Leave-One-Out)
        others_idx = [1:s-1, s+1:n_subs];
        d_lower = mean(rdm_vecs_w(:, others_idx), 2);
        TSS_lower = sum(d_lower.^2); 
        
        beta_low = lsqnonneg(M_subj, d_lower);
        RSS_low = sum((d_lower - (M_subj * beta_low)).^2);
        r2_low_val = max(1 - (RSS_low / TSS_lower), 0);
        z_lower(s) = atanh(sqrt(r2_low_val)); 
    end
    
    % Average the Fisher Z scores
    lower_bound_Z = mean(z_lower, 'omitnan');
    upper_bound_Z = mean(z_upper, 'omitnan');
end

function v = get_vec(rdm)
    v = rdm(tril(true(size(rdm)), -1));
end
