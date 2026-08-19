%% PRONET RSA main pipeline WITH pretest behavioral data
clear; clc;

%% SETUP (paths, globals, SPM, RSA toolbox)
addpath /mnt/beegfs/users/antonia.ceric/git/pronet/code/matlab/func;
addpath /mnt/beegfs/users/antonia.ceric/git/pronet/code/matlab/func/pretest_modelling;
run('pronet_startup.m');

%% === Subject list (VP) ===
VP = {'02','03','04','05','06','07','08','09','10',...
      '11','12','13','15','16','17','18','19','20',...
      '21','22','23','24','25','26','27','28','29','31','32'};

% VP = {'02'}

%% BEHAVIORAL
% === 0. Model RDMs (Feb 2026 Versions) ===
fprintf('\n[0] Building Behavioral Model RDMs...\n');

for i = 1:length(VP)
    subj = VP{i};
    fprintf(' - Processing sub-%s...\n', subj);

    % Version A: Speaker-wise Fits (m1 vs f2 separate, log-aligned)
    % Saves to: model_rdms_pretest_Feb2026_speakerwisefit/

        pronet_build_modelRDMs_pretest_Feb2026_speakerwisefit(subj);


    % Version B: Aggregate Fit (All sessions/speakers pooled, log-aligned)
    % Saves to: model_rdms_pretest_Feb2026_aggregatefit/

        pronet_build_modelRDMs_pretest_Feb2026_aggfit(subj);

end

fprintf('\nBehavioral Modelling Complete.\n');

%% === 0.1 Plot both version model RDMs (Feb 2026 Versions) ===
fprintf('\n[0.3] Plotting Model RDMs...\n');

% Plot the Speaker-wise version (Generates plots for m1, f2, and avg)
pronet_plot_modelRDMs_Feb2026(VP, HPC_PATH, 'speakerwise');

% Plot the Aggregate version (Generates one plot for the pooled data)
pronet_plot_modelRDMs_Feb2026(VP, HPC_PATH, 'aggregate');


%% CHECK which RDM do we use now?

% best = pronet_select_best_pretest_RDM(VP, HPC_PATH, false);


%% NOTE
% % Do not need all of these when MOCS load('/mnt/beegfs/workspace/2024-0404-PRONET/TMS/Pre-Session/results/02_PRONET_MOCS_results.mat')
% % already has concatenated fits of the 3 pre-test sessions
% 
% %% 0.4 model RDMs
% % fprintf(' - Building model RDMs...\n');
% % 
% %      for i = 1:length(VP)
% %         subject_id = VP{i};
% % 
% %         fprintf('Processing subject: %s\n', subject_id);
% %         pronet_build_modelRDMs_pretest(subject_id); 
% %     end
% 
% %% Plot Model RDMs
% % fprintf('\n--- Plotting Group RDMs ---\n');
% % 
% % plot_model_rdms_pretest(VP, HPC_PATH);
% % 
% % fprintf('Done plotting.\n');

%% NEURAL
 === 1. GLM fitting (Speaker-wise Feb2026) 183 regressors ===
fprintf('\n[1] Preparing Speaker-wise GLM jobs (Feb2026)...\n');

Jobs_GLM = {};
for i = 1:numel(VP)
    Jobs_GLM{end+1} = VP{i}; %#ok<AGROW>
end

InitSh  = 'module load matlab'; 
InitCmd = 'addpath /mnt/beegfs/users/antonia.ceric/git/pronet/code/matlab/func; pronet_startup;';

if ~isempty(Jobs_GLM)
    fprintf('Submitting %d Speaker-wise GLM jobs to SLURM...\n', numel(Jobs_GLM));
    slurmitout(@pronet_GLM_Speakerwise_Feb2026, Jobs_GLM, InitCmd, InitSh);
else
    fprintf('No GLM jobs to run.\n');
end

fprintf('\nGLM Fitting submitted to cluster.\n');
%%
% % === 1. GLM fitting (Speaker-wise Feb2026) 108 ===
% fprintf('\n[1] Preparing Speaker-wise GLM jobs (March2026)...\n');
% 
% Jobs_GLM = {};
% for i = 1:numel(VP)
%     Jobs_GLM{end+1} = VP{i}; %#ok<AGROW>
% end
% 
% InitSh  = 'module load matlab'; 
% InitCmd = 'addpath /mnt/beegfs/users/antonia.ceric/git/pronet/code/matlab/func; pronet_startup;';
% 
% if ~isempty(Jobs_GLM)
%     fprintf('Submitting %d Speaker-wise GLM jobs to SLURM...\n', numel(Jobs_GLM));
%     slurmitout(@pronet_GLM_Speakerwise_bigblocks_March2026, Jobs_GLM, InitCmd, InitSh);
% else
%     fprintf('No GLM jobs to run.\n');
% end
% 
% fprintf('\nGLM Fitting submitted to cluster.\n');

%% === 2. ROI voxel patterns (Feb2026) ===
% fprintf('\n[2] Preparing Speaker-wise ROI pattern jobs (Feb2026)...\n');
% 
% Jobs_ROI = {};
% for i = 1:numel(VP)
%     J = struct();
%     J.VP        = VP{i};        
%     J.DATA_PATH = DATA_PATH;    
%     J.hcp_dir   = hcp_dir;      
% 
%     Jobs_ROI{end+1} = J; %#ok<AGROW>
% end
% 
% InitSh  = 'module load matlab'; 
% InitCmd = 'addpath /mnt/beegfs/users/antonia.ceric/git/pronet/code/matlab/func; pronet_startup;';
% 
% if ~isempty(Jobs_ROI)
%     fprintf('Submitting %d ROI extraction jobs to SLURM...\n', numel(Jobs_ROI));
%     slurmitout(@pronet_roi_patterns_Feb2026, Jobs_ROI, InitCmd, InitSh);
% else
%     fprintf('No ROI jobs to run.\n');
% end

%% === 4. Crossnobis (SESSION-WISE & SPEAKER-WISE) ===
% fprintf('\n[4] Preparing Speaker-wise Crossnobis jobs ...\n');
% 
% ldc_spk_dir = fullfile(DATA_PATH, 'rsa_crossnobis_bigblocks_March2026');
% if ~exist(ldc_spk_dir, 'dir'), mkdir(ldc_spk_dir); end
% 
% Jobs_LDC = {};
% for i = 1:numel(VP)
%     for task_index = 1:2
%         for roi_index = 1:10
%             J = struct();
%             J.VP         = VP{i};
%             J.DATA_PATH  = DATA_PATH;
%             J.task_index = task_index;
%             J.roi_index  = roi_index;
% 
%             Jobs_LDC{end+1} = J; %#ok<AGROW>
%         end
%     end
% end
% 
% InitSh  = 'module load matlab'; 
% InitCmd = 'addpath /mnt/beegfs/users/antonia.ceric/git/pronet/code/matlab/func; pronet_startup;';
% 
% 
% if ~isempty(Jobs_LDC)
%     fprintf('Submitting %d Speaker-wise LDC jobs to SLURM...\n', numel(Jobs_LDC));
%     slurmitout(@pronet_crossnobis_Feb2026, Jobs_LDC, InitCmd, InitSh);
% end
% 
% fprintf('\nSpeaker-wise Session-wise analysis submitted to cluster.\n');


%% Plot Crossnobis

% for i = 1:numel(VP)
%     fprintf('Plotting Speakerwise RDMs for sub-%s...\n', VP{i});
%     pronet_crossnobis_plot_Speakerwise_Feb2026(VP{i}, DATA_PATH, 'save', true);
% end

%% === 5. RSA Model Fitting (SESSION-WISE) ===
% fprintf('\n[5] Running RSA Model Fitting (SESSION-WISE)...\n');
% 
% %1. Define directories
% ldc_dir   = fullfile(DATA_PATH, 'rsa_crossnobis_sessionwise'); 
% model_dir = fullfile(DATA_PATH, 'model_rdms_pretest_Feb2026_aggregatefit');
% res_dir   = fullfile(DATA_PATH, 'rsa_results_new93reg'); 
% if ~exist(res_dir, 'dir'), mkdir(res_dir); end
% 
% %2. Check if Model RDMs exist
% have_models = ~isempty(dir(fullfile(model_dir, 'sub-*_modelRDM_pretest_Aggregate.mat')));
% 
% if ~have_models
%     warning('Skipping RSA fitting: Model RDM files are missing.');
% else
%     %Initialize storage structure
%     all_subject_results = struct();
% 
%     tasks_to_analyze = {'phoneme', 'intonat'}; 
%     rois_to_analyze  = {'lIFG', 'lPMC', 'lPAC', 'lpSTS', 'laSTS', ...
%                         'rIFG', 'rPMC', 'rPAC', 'rpSTS', 'raSTS'};
%     sessions = {'lt', 'rt', 'vx'}; % Define sessions
% 
%     %=========================================================
%     %3. HARDCODED EXCLUSIONS
%     %=========================================================
%     %Subjects to exclude entirely from specific tasks
% 
%     excl_phoneme = {'08', '05'};
%     excl_intonat = {'13', '05'};
% 
%     %=========================================================
% 
%     for i = 1:numel(VP)
%         current_subj = VP{i};
%         subj_key = ['sub' current_subj];
%         fprintf(' - Fitting models for sub-%s...\n', current_subj);
% 
%         try
%             %--- A. Load Model RDMs ---
%             model_file = fullfile(model_dir, sprintf('sub-%s_modelRDM_pretest_Aggregate.mat', current_subj));
%             if ~exist(model_file, 'file'), continue; end
%             model_data = load(model_file); 
% 
%             %--- B. Loop over tasks and ROIs ---
%             for t = 1:numel(tasks_to_analyze)
%                 task = tasks_to_analyze{t}; 
% 
%                 %[CHECK 1] Task Exclusion
%                 %Select the correct list
%                 if strcmp(task, 'phoneme')
%                     current_excl = excl_phoneme;
%                 else
%                     current_excl = excl_intonat;
%                 end
% 
%                 %Skip subject if in exclusion list
%                 if ismember(current_subj, current_excl)
%                     continue; 
%                 end
% 
%                 for r = 1:numel(rois_to_analyze)
%                     roi = rois_to_analyze{r}; 
% 
%                     %--- Load Specific SESSIONWISE File ---
%                     fname = sprintf('sub%s_task-%s_roi-%s_LDC_SESSIONWISE.mat', current_subj, task, roi);
%                     neural_file = fullfile(ldc_dir, fname);
% 
%                     if ~exist(neural_file, 'file'), continue; end
% 
%                     tmp = load(neural_file, 'Sub');
% 
%                     %Check structure path
%                     if ~isfield(tmp.Sub, task) || ~isfield(tmp.Sub.(task), roi)
%                         continue;
%                     end
% 
%                     ROI_Data = tmp.Sub.(task).(roi);
% 
%                     %--- Loop over Sessions (lt, rt, vx) ---
%                     for s = 1:3
%                         sess_name = sessions{s};
% 
%                         % %[CHANGE] 1. Construct Session-Specific Model Name
%                         % %Force use of '...' model for ALL sessions
%                         model_str = [task '_pretest'];
% 
%                         %2. Check/Load Specific Model
%                         if ~isfield(model_data.modelRDMs, model_str)
%                             continue; 
%                         end
% 
%                         model_RDM_Acoustic    = model_data.modelRDMs.(model_str).aggregate.acoustic;
%                         model_RDM_Categorical = model_data.modelRDMs.(model_str).aggregate.category;
% 
%                         %3. Check if Neural Data exists
%                         if isfield(ROI_Data, sess_name) && isfield(ROI_Data.(sess_name), 'RDM')
%                             neural_RDM = ROI_Data.(sess_name).RDM;
% 
%                             %--- Perform Fitting ---
%                             [R_corrected, R_real, R_perm] = pronet_single_perm_test_pretest(...
%                                 neural_RDM, model_RDM_Acoustic, model_RDM_Categorical);
% 
%                             %--- Store Results ---
%                             %Corrected
%                             all_subject_results.(subj_key).(task).(roi).(sess_name).R_general = R_corrected.R_general;
%                             all_subject_results.(subj_key).(task).(roi).(sess_name).R_unique_Acoustic = R_corrected.R_unique_Acoustic;
%                             all_subject_results.(subj_key).(task).(roi).(sess_name).R_unique_Categorical = R_corrected.R_unique_Categorical;
% 
%                             %Uncorrected
%                             all_subject_results.(subj_key).(task).(roi).(sess_name).R_general_uncorrected = R_real.R_general;
%                             all_subject_results.(subj_key).(task).(roi).(sess_name).R_unique_Acoustic_uncorrected = R_real.R_unique_Acoustic;
%                             all_subject_results.(subj_key).(task).(roi).(sess_name).R_unique_Categorical_uncorrected = R_real.R_unique_Categorical;
% 
%                             %Permuted
%                             all_subject_results.(subj_key).(task).(roi).(sess_name).R_general_permuted = R_perm.R_general;
%                             all_subject_results.(subj_key).(task).(roi).(sess_name).R_unique_Acoustic_permuted = R_perm.R_unique_Acoustic;
%                             all_subject_results.(subj_key).(task).(roi).(sess_name).R_unique_Categorical_permuted = R_perm.R_unique_Categorical;
%                         end
%                     end % End Session Loop
% 
%                 end % End ROI Loop
%             end % End Task Loop
% 
%         catch ME
%             warning('RSA fitting FAILED for sub-%s: %s', current_subj, ME.message);
%         end
%     end % End Subject Loop
% 
%     %--- 3. Save Final Group Results ---
%     results_file = fullfile(res_dir, 'rsa_whitened_perm_results_new93reg.mat');
%     save(results_file, 'all_subject_results');
%     fprintf('[5] Session-wise RSA fitting complete. Saved to %s\n', results_file);
% end
% 
%%  === 6. Statistical Analysis (SESSION-WISE) ===

fprintf('\n[6] Starting Session-wise Group Statistics...\n');
savehere_dir = fullfile(DATA_PATH, 'rsa_results_new93reg');
if ~exist(savehere_dir, 'dir'), mkdir(savehere_dir); end
results_dir = fullfile(DATA_PATH, 'rsa_results_new93reg');
results_file = fullfile(results_dir, 'rsa_whitened_perm_results_new93reg.mat');

if exist(results_file, 'file')
    load(results_file, 'all_subject_results'); 

    rois_to_analyze  = {'lIFG', 'lPMC', 'lPAC', 'lpSTS', 'laSTS', ...
                        'rIFG', 'rPMC', 'rPAC', 'rpSTS', 'raSTS'};
    tasks_to_analyze = {'phoneme', 'intonat'}; 

    stats_results = pronet_rsa_statistics_sessionwise_pretest(...
        all_subject_results, VP, rois_to_analyze, tasks_to_analyze, DATA_PATH);

    save(fullfile(savehere_dir, 'rsa_stats_tables_new93reg.mat'), 'stats_results');

else
    warning('[6] Session-wise results file not found: %s', results_file);
end

%% === 7. ANOVA Statistics (Pretest Models) ===
% fprintf('\n[7] Running ANOVA on Session-wise RSA Results (Pretest Models)...\n');
% 
% stats_input_dir = fullfile(DATA_PATH, 'rsa_results_new93reg');
% stats_input_file = fullfile(stats_input_dir, 'rsa_whitened_perm_results_new93reg.mat');
% 
% if exist(stats_input_file, 'file')
%     fprintf('Loading results from: %s\n', stats_input_file);
%     load(stats_input_file, 'all_subject_results'); % Loads the variable into the workspace
% 
%     stats_results = pronet_rsa_anova_pretest(...
%         all_subject_results, ...   % The struct containing the data
%         VP, ...                    % The subject list
%         DATA_PATH);                % Path for saving output
% 
% else
%     warning('Results file not found: %s\nPlease run Step 5 first.', stats_input_file);
% end