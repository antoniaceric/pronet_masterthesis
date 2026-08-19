function pronet_build_modelRDMs_pretest_Feb2026_speakerwisefit(subject_id)
    % PRONET_BUILD_MODELRDMS_PRETEST_FEB2026_SPEAKERWISE
    %
    % Logic:
    % 1. Uses Pre-test psychometric parameters as the pure "Ground Truth".
    % 2. Evaluates the curves for m1 and f2 at standard ordinal levels (1-5).
    % 3. Builds separate RDMs for m1 and f2.
    % 4. Averages the resulting distance matrices.

    global HPC_PATH

    %% 1. SETUP & PATHS
    MOCS_DIR  = '/mnt/beegfs/workspace/2024-0404-PRONET/TMS/Pre-Session/results';
    outputDir = fullfile(HPC_PATH, 'DATA', 'proc', 'model_rdms_pretest_Feb2026_speakerwisefit');
    if ~exist(outputDir, 'dir'), mkdir(outputDir); end

    sub_num = regexprep(subject_id, 'sub-', '');

    %% 2. LOAD PRE-TEST DATA
    master_file = fullfile(MOCS_DIR, sprintf('%s_PRONET_MOCS_results.mat', sub_num));
    if ~exist(master_file, 'file')
        warning('Pre-test master file missing: %s', master_file); return;
    end
    dat = load(master_file);

    voices = {'m1','f2'};
    tasks  = {'phoneme','intonat'};
    modelRDMs = struct();

    %% 3. FIXED ACOUSTIC RDM (Ordinal steps 1 to 5)
    x_ord = (1:5)';
    rdm_acoustic = normalizeRDM(squareform(pdist(x_ord, 'euclidean').^2));

    %% 4. MAIN LOOP
    for t = 1:numel(tasks)
        task = tasks{t};
        key  = sprintf('%s_pretest', task);
        
        for v = 1:numel(voices)
            spk = voices{v};
            
            if isfield(dat, spk) && isfield(dat.(spk), task) && isfield(dat.(spk).(task), 'PFfit')
                % Extract pre-test parameters [mu, sigma, gamma, lambda]
                params = dat.(spk).(task).PFfit.params; 
                
                % Evaluate curve at the 5 standard steps
                y = eval_psychometric(params, x_ord);

                % Build Speaker-specific categorical RDM (Squared Euclidean)
                rdm_cat = squareform(pdist(y, 'euclidean').^2);
                
                modelRDMs.(key).(spk).category = normalizeRDM(rdm_cat);
                modelRDMs.(key).(spk).acoustic = rdm_acoustic;
                modelRDMs.(key).(spk).y = y; % Store for debugging
            end
        end

        %% 5. AVERAGE MATRICES
        if isfield(modelRDMs, key) && isfield(modelRDMs.(key), 'm1') && isfield(modelRDMs.(key), 'f2')
            R_avg = (modelRDMs.(key).m1.category + modelRDMs.(key).f2.category) / 2;
            modelRDMs.(key).avg.category = normalizeRDM(R_avg);
            modelRDMs.(key).avg.acoustic = rdm_acoustic;
        end
    end

    %% 6. SAVE
    save(fullfile(outputDir, sprintf('sub-%s_modelRDM_pretest_Speakerwise.mat', sub_num)), 'modelRDMs');
    fprintf('Saved pure PRETEST SPEAKERWISE model RDMs for sub-%s\n', sub_num);
end

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function y = eval_psychometric(p, x)
    mu = p(1); sigma = p(2); gamma = p(3); lambda = p(4);
    y = gamma + (1 - gamma - lambda) .* normcdf(x(:), mu, sigma);
end

function R = normalizeRDM(R)
    range_val = max(R(:)) - min(R(:));
    if range_val > 0, R = (R - min(R(:))) / range_val; else, R = zeros(size(R)); end
end

%% OLD     %% 1. SETUP & PATHS
% %     MOCS_DIR  = '/mnt/beegfs/workspace/2024-0404-PRONET/TMS/Pre-Session/results';
% %     BEHAV_DIR = '/home/antonia.ceric/git/pronet/data/proc';
% % 
% %     outputDir = fullfile(HPC_PATH, 'DATA', 'proc', 'model_rdms_pretest_Feb2026_speakerwisefit');
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
% %     % Read table (ignoring broken headers to avoid variable name issues)
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
% %         % Extract Phoneme Morphs (e.g., 'p09' -> 9) directly from strings to bypass the Level Bug
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
% %     %% 4. EVALUATE FIT & BUILD RDMs
% %     tasks = {'intonat', 'phoneme'};
% %     modelRDMs = struct();
% % 
% %     % Fixed Acoustic RDM: Since morph steps in the scanner were always 5 units apart,
% %     % the physical distance is perfectly represented by ordinal steps 1:5.
% %     rdm_acoustic = normalizeRDM(squareform(pdist([1; 2; 3; 4; 5], 'euclidean').^2));
% % 
% %     for t = 1:numel(tasks)
% %         task = tasks{t};
% % 
% %         for v = 1:length(voices)
% %             s = voices{v};
% % 
% %             % Check if the pre-test data contains the required fields
% %             if isfield(dat, s) && isfield(dat.(s), task) && isfield(dat.(s).(task), 'PFfit')
% %                 params = dat.(s).(task).PFfit.params;
% % 
% %                 % Retrieve the actual scanner morphs (e.g., 10, 15, 20, 25, 30)
% %                 morph_raw = scanner_morphs.(s).(task);
% %                 if length(morph_raw) ~= 5
% %                     warning('Subject %s, Task %s: Found %d morphs instead of 5.', sub_num, task, length(morph_raw));
% %                     continue; 
% %                 end
% % 
% %                 % Map to Pre-test Domain: Level 1 in pre-test = Morph 5, Level 2 = Morph 10, etc.
% %                 x_fit_scale = morph_raw / 5;
% % 
% %                 % [THE GROUND TRUTH STEP]: Query the pre-test curve for the scanner morphs
% %                 y_probs = eval_psychometric(params, x_fit_scale);
% % 
% %                 % Compute Categorical RDM (Squared Euclidean)
% %                 rdm_cat_raw = squareform(pdist(y_probs(:), 'euclidean').^2);
% % 
% %                 % Store individual speaker RDMs
% %                 key = sprintf('%s_pretest', task);
% %                 modelRDMs.(key).(s).category = normalizeRDM(rdm_cat_raw);
% %                 modelRDMs.(key).(s).acoustic = rdm_acoustic; 
% % 
% %                 % Save metadata for transparency
% %                 modelRDMs.(key).(s).raw_morphs = morph_raw; 
% %                 modelRDMs.(key).(s).x_fit_query = x_fit_scale;
% %                 modelRDMs.(key).(s).y_response  = y_probs;
% %             end
% %         end
% % 
% %         % --- SPEAKER-WISE MATRIX AVERAGING ---
% %         key = sprintf('%s_pretest', task);
% %         if isfield(modelRDMs, key) && isfield(modelRDMs.(key), 'm1') && isfield(modelRDMs.(key), 'f2')
% % 
% %             % Combine m1 and f2 distance matrices
% %             avg_cat_matrix = (modelRDMs.(key).m1.category + modelRDMs.(key).f2.category) / 2;
% % 
% %             % Store in the .avg field (normalize again to ensure 0-1 range)
% %             modelRDMs.(key).avg.category = normalizeRDM(avg_cat_matrix);
% %             modelRDMs.(key).avg.acoustic = rdm_acoustic;
% %         end
% %     end
% % 
% %     %% 5. SAVE
% %     if ~isempty(fieldnames(modelRDMs))
% %         save_path = fullfile(outputDir, sprintf('sub-%s_modelRDM_pretest_Speakerwise.mat', sub_num));
% %         save(save_path, 'modelRDMs');
% %         fprintf('Successfully saved model RDMs for sub-%s\n', sub_num);
% %     else
% %         fprintf('No valid RDM data generated for sub-%s. Check logfiles/master file.\n', sub_num);
% %     end
% % end
% % 
% % % =========================================================================
% % % HELPER FUNCTIONS
% % % =========================================================================
% % 
% % function y = eval_psychometric(p, x)
% %     % Evaluates a cumulative Gaussian distribution using pre-test parameters
% %     % p: [Threshold (mu), Slope (sigma), Guess (gamma), Lapse (lambda)]
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