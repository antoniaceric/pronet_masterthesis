function pronet_build_modelRDMs_pretest_Feb2026_aggfit(subject_id)
    % PRONET_BUILD_MODELRDMS_PRETEST_FEB2026_AGGFIT
    %
    % Logic:
    % 1. Pools raw trial counts from pre-test for both voices.
    % 2. Fits a SINGLE psychometric curve (The Aggregate Ground Truth).
    % 3. Evaluates this curve at standard ordinal levels 1-5.
    % 4. Generates a unified 5x5 RDM.

    global HPC_PATH
    
    %% 1. SETUP & PATHS
    MOCS_DIR  = '/mnt/beegfs/workspace/2024-0404-PRONET/TMS/Pre-Session/results';
    outputDir = fullfile(HPC_PATH, 'DATA', 'proc', 'model_rdms_pretest_Feb2026_aggregatefit');
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end

    sub_num = regexprep(subject_id, 'sub-', '');

    %% 2. LOAD PRE-TEST DATA
    master_file = fullfile(MOCS_DIR, sprintf('%s_PRONET_MOCS_results.mat', sub_num));
    if ~exist(master_file, 'file')
        warning('Pre-test master file missing: %s', master_file); return;
    end
    dat = load(master_file);

    tasks = {'phoneme','intonat'};
    modelRDMs = struct();
    
    %% 3. FIXED ACOUSTIC RDM (Ordinal steps 1 to 5)
    x_ord = (1:5)';
    rdm_acoustic = normalizeRDM(squareform(pdist(x_ord, 'euclidean').^2));

    %% 4. MAIN LOOP
    for t = 1:numel(tasks)
        task = tasks{t};
        key  = sprintf('%s_pretest', task);

        if isfield(dat, 'm1') && isfield(dat.m1, task) && isfield(dat, 'f2') && isfield(dat.f2, task)
            
            % --- A) Pool raw Pre-test data ---
            k = dat.m1.(task).yes(:) + dat.f2.(task).yes(:);
            n = dat.m1.(task).n(:)   + dat.f2.(task).n(:);

            % --- B) Fit Global Curve (Maximum Likelihood) ---
            initial = [3 1];
            fit_params = fminsearch(@(p) psych_likelihood_scalar(p, x_ord, k, n), initial);
            global_params = [fit_params(1), abs(fit_params(2)), 0.01, 0.01];

            % --- C) Evaluate Curve and Build RDM ---
            y = eval_psychometric(global_params, x_ord);
            rdm_cat = squareform(pdist(y, 'euclidean').^2);

            modelRDMs.(key).aggregate.category = normalizeRDM(rdm_cat);
            modelRDMs.(key).aggregate.acoustic = rdm_acoustic;
            modelRDMs.(key).aggregate.y = y;
            modelRDMs.(key).aggregate.global_params = global_params;
        else
            warning('Task %s data missing for m1 or f2 in pre-test for sub-%s', task, sub_num);
        end
    end

    %% 5. SAVE
    save(fullfile(outputDir, sprintf('sub-%s_modelRDM_pretest_Aggregate.mat', sub_num)), 'modelRDMs');
    fprintf('Saved pure PRETEST AGGREGATE model RDMs for sub-%s\n', sub_num);
end

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function negLogLik = psych_likelihood_scalar(params, x, k, n)
    mu = params(1); sigma = abs(params(2)) + eps;
    p = 0.01 + (1 - 0.02) .* normcdf(x, mu, sigma);
    p = min(max(p, 1e-4), 1 - 1e-4);
    negLogLik = -sum(k .* log(p) + (n - k).*log(1-p));
end

function y = eval_psychometric(p, x)
    mu = p(1); sigma = p(2); gamma = p(3); lambda = p(4);
    y = gamma + (1-gamma-lambda).*normcdf(x(:), mu, sigma);
end

function R = normalizeRDM(R)
    r = max(R(:)) - min(R(:));
    if r > 0, R = (R - min(R(:))) / r; else, R = zeros(size(R)); end
end
%% OLD    %% 1. SETUP & PATHS
% %     MOCS_DIR  = '/mnt/beegfs/workspace/2024-0404-PRONET/TMS/Pre-Session/results';
% %     BEHAV_DIR = '/home/antonia.ceric/git/pronet/data/proc';
% % 
% %     outputDir = fullfile(HPC_PATH, 'DATA', 'proc', 'model_rdms_pretest_Feb2026_aggregatefit');
% %     if ~exist(outputDir, 'dir'), mkdir(outputDir); end
% % 
% %     % Standardize Subject ID
% %     sub_num = regexprep(subject_id, 'sub-', ''); 
% %     sub_full = ['sub-', sub_num];
% % 
% %     %% 2. LOAD MASTER PRE-TEST FILE
% %     master_file = fullfile(MOCS_DIR, sprintf('%s_PRONET_MOCS_results.mat', sub_num));
% %     if ~exist(master_file, 'file')
% %         warning('Pre-test master file not found: %s', master_file); return;
% %     end
% %     dat = load(master_file);
% % 
% %     %% 3. READ VX SCANNER LOGFILE (Extracting true morph values)
% %     log_pattern = fullfile(BEHAV_DIR, sub_full, 'behav', '*_vx.txt');
% %     log_files = dir(log_pattern);
% % 
% %     if isempty(log_files)
% %         warning('No VX-logfile found for %s', sub_full); return; 
% %     end
% % 
% %     file_path = fullfile(log_files(1).folder, log_files(1).name);
% %     T = readtable(file_path, 'FileType', 'text', 'ReadVariableNames', false, 'HeaderLines', 1);
% % 
% %     T.Properties.VariableNames = {'SubjID', 'Session', 'Block', 'Trial', ...
% %                                   'Task', 'Button', 'RT', 'SoundTime', ...
% %                                   'Response', 'SOA', 'Speaker', 'P_Morph', ...
% %                                   'P_Level', 'I_Morph', 'I_Level', 'Wavefile'};
% % 
% %     scanner_morphs = struct();
% %     voices = {'m1', 'f2'};
% % 
% %     for v = 1:length(voices)
% %         s = voices{v};
% %         idx_spk = strcmpi(T.Speaker, s);
% % 
% %         % Extract Phoneme Morphs directly from strings (e.g., 'p09' -> 9)
% %         p_strs = T.P_Morph(idx_spk & strcmpi(T.Task, 'phoneme'));
% %         if iscell(p_strs)
% %             p_nums = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), p_strs);
% %         else, p_nums = p_strs; end
% %         scanner_morphs.(s).phoneme = unique(p_nums(~isnan(p_nums)));
% % 
% %         % Extract Intonation Morphs (e.g., 'i35' -> 35)
% %         i_strs = T.I_Morph(idx_spk & strcmpi(T.Task, 'intonat'));
% %         if iscell(i_strs)
% %             i_nums = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), i_strs);
% %         else, i_nums = i_strs; end
% %         scanner_morphs.(s).intonat = unique(i_nums(~isnan(i_nums)));
% %     end
% % 
% %     %% 4. POOL DATA, FIT GLOBAL CURVE, & BUILD SINGLE RDM
% %     tasks = {'intonat', 'phoneme'};
% %     modelRDMs = struct();
% % 
% %     % The physical Acoustic RDM. Since morph steps were consistently 5 units apart,
% %     % an ordinal 1:5 scale provides a perfectly proportional distance model.
% %     x_ord = [1; 2; 3; 4; 5];
% %     rdm_acoustic = normalizeRDM(squareform(pdist(x_ord, 'euclidean').^2));
% % 
% %     for t = 1:numel(tasks)
% %         task = tasks{t};
% % 
% %         if isfield(dat, 'm1') && isfield(dat.m1, task) && isfield(dat, 'f2') && isfield(dat.f2, task)
% % 
% %             % --- A. Pool the Pre-test Data across both speakers ---
% %             k_agg = dat.m1.(task).yes + dat.f2.(task).yes;
% %             n_agg = dat.m1.(task).n   + dat.f2.(task).n;
% % 
% %             % --- B. Fit the Global Psychometric Curve ---
% %             % We use a Maximum Likelihood fit to define the ground truth parameters
% %             initial_guess = [3, 1]; 
% %             fit_params = fminsearch(@(p) psych_likelihood(p, x_ord, k_agg, n_agg), initial_guess);
% %             % Format: [Threshold, Slope, Guess(0.01), Lapse(0.01)]
% %             global_params = [fit_params(1), abs(fit_params(2)), 0.01, 0.01];
% % 
% %             % --- C. Extract Speaker-Specific Scanner Morphs ---
% %             m1_raw = scanner_morphs.m1.(task);
% %             f2_raw = scanner_morphs.f2.(task);
% % 
% %             if length(m1_raw) ~= 5 || length(f2_raw) ~= 5
% %                 warning('Task %s: Missing morphs in scanner log for %s', task, sub_num); continue;
% %             end
% % 
% %             % --- D. Evaluate Speaker-Specific Points on the Global Curve ---
% %             % We divide the raw morphs by 5 to map them into the fit's domain (1 to 5)
% %             y_m1 = eval_psychometric(global_params, m1_raw / 5);
% %             y_f2 = eval_psychometric(global_params, f2_raw / 5);
% % 
% %             % Combine probabilities to create ONE unified categorical vector
% %             y_combined = (y_m1(:) + y_f2(:)) / 2;
% % 
% %             % Build the single Categorical RDM (Squared Euclidean)
% %             rdm_cat_raw = squareform(pdist(y_combined, 'euclidean').^2);
% % 
% %             % --- E. Store under .aggregate ---
% %             key = sprintf('%s_pretest', task);
% %             modelRDMs.(key).aggregate.category = normalizeRDM(rdm_cat_raw);
% %             modelRDMs.(key).aggregate.acoustic = rdm_acoustic; 
% % 
% %             % Metadata for validation
% %             modelRDMs.(key).aggregate.global_params = global_params;
% %             modelRDMs.(key).aggregate.m1_morphs_used = m1_raw;
% %             modelRDMs.(key).aggregate.f2_morphs_used = f2_raw;
% %             modelRDMs.(key).aggregate.y_probs_combined = y_combined;
% %         else
% %             warning('Task %s data missing for m1 or f2 in pre-test for %s', task, sub_num);
% %         end
% %     end
% % 
% %     %% 5. SAVE
% %     if ~isempty(fieldnames(modelRDMs))
% %         save_path = fullfile(outputDir, sprintf('sub-%s_modelRDM_pretest_Aggregate.mat', sub_num));
% %         save(save_path, 'modelRDMs');
% %         fprintf('Successfully saved Aggregate RDMs for sub-%s\n', sub_num);
% %     end
% % end
% % 
% % % =========================================================================
% % % HELPER FUNCTIONS
% % % =========================================================================
% % 
% % function negLogLik = psych_likelihood(params, x, k, n)
% %     % Maximum Likelihood estimation for psychometric curve parameters
% %     x = x(:); k = k(:); n = n(:);
% %     mu = params(1); sigma = abs(params(2)) + eps;
% %     gamma = 0.01; lambda = 0.01; 
% % 
% %     p_model = gamma + (1 - gamma - lambda) .* normcdf(x, mu, sigma);
% %     p_model = min(max(p_model, 0.0001), 0.9999); % Avoid log(0)
% % 
% %     negLogLik = -sum(k .* log(p_model) + (n - k) .* log(1 - p_model));
% % end
% % 
% % function y = eval_psychometric(p, x)
% %     % Standard Cumulative Gaussian evaluation
% %     mu = p(1); sigma = p(2);
% %     gamma = p(3); lambda = p(4);
% %     y = gamma + (1 - gamma - lambda) .* normcdf(x, mu, sigma);
% % end
% % 
% % function R = normalizeRDM(R)
% %     % Scales RDM values to exactly [0, 1] range
% %     range_val = max(R(:)) - min(R(:));
% %     if range_val > 0
% %         R = (R - min(R(:))) / range_val; 
% %     else
% %         R = zeros(size(R)); 
% %     end
% % end