function pronet_GLM_Speakerwise_bigblocks_March2026(VP)
    % PRONET_GLM_SPEAKERWISE_BIGBLOCKS_MARCH2026 (MEGA-SESSION CONCAT)
    % 
    % 1. Creates Speaker-separated design matrices (m1 vs f2).
    % 2. 90-Regressor Model: Groups 6 original blocks into 3 Super-Blocks.
    % 3. CONCAT FIX: Automatically shifts onsets by (nV * TR).
    % 4. Fits a concatenated mega-session GLM (18 motion regressors).

    global HPC_PATH DATA_PATH
    
    % === Configuration ===
    subj_id = VP(1:2);
    stim_cond = {'lt','rt','vx'};
    task_list = {'phoneme','intonat'};
    n_levels = 5;
    n_super_blocks = 3; % Neu: 3 Super-Blöcke statt 6 Einzelblöcken
    speakers_to_split = {'m1', 'f2'};
    
    % Output Paths (Updated for March2026)
    DESIGNDIR = fullfile(DATA_PATH, 'design_matrices_bigblocks_March2026');
    GLM_BASE  = fullfile(DATA_PATH, ['sub-' subj_id], 'rsa_firstlevel_bigblocks_March2026');
    if ~exist(DESIGNDIR, 'dir'), mkdir(DESIGNDIR); end
    
    spm('defaults','FMRI'); spm_jobman('initcfg');
    getVols = @(dir,pat) spm_select('ExtFPList', dir, pat, Inf);

    subj_folder = fullfile(DATA_PATH, ['sub-', subj_id]);
    
    %% --- STEP 1: CONCATENATED DESIGN MATRIX GENERATION ---
    for t = 1:length(task_list)
        taskName = task_list{t};
        onsets = {}; durations = {}; names = {}; reg_count = 0;
        
        cum_time = 0; 
        TR_sec = [];  
        
        for s = 1:length(stim_cond)
            cond = stim_cond{s};
            
            sesdir = fullfile(subj_folder, ['ses-' cond]);
            vols = getVols(sesdir, '^wuabold-task_unwarped\.nii$');
            if isempty(vols)
                warning('No volumes found for %s %s', subj_id, cond);
                continue; 
            end
            nV = size(vols, 1);
            
            if isempty(TR_sec)
                TR_sec = spm_vol(vols(1,:)).private.timing.tspace; 
            end

            logdir = fullfile(subj_folder, 'behav');
            logfile = dir(fullfile(logdir, sprintf('sub-%s_*_%s.txt', subj_id, cond)));
            
            if isempty(logfile)
                cum_time = cum_time + (nV * TR_sec);
                continue; 
            end

            T = readtable(fullfile(logdir, logfile(1).name), 'delimiter', '\t', 'ReadVariableNames', false);
            T.Properties.VariableNames = {'SubjID','Session','Block','Trial','Task','Button','RT','SoundTime', ...
                                          'Response','SOA','Speaker','P_Morph','P_Level','I_Morph','I_Level','Wavefile'};

            if iscell(T.SoundTime)
                onset_times = cellfun(@str2double, T.SoundTime);
                bsl = str2double(T.SoundTime{4});
            else
                onset_times = T.SoundTime;
                bsl = T.SoundTime(4);
            end 
    
            % Mapping Logic (correct Morph-Level)
            ilevel_all = zeros(height(T),1);
            for spkIdx = 1:2
                spk = speakers_to_split{spkIdx};
                sel_spk = strcmp(T.Speaker, spk);
                morphVal = cellfun(@(x) str2double(x(2:end)), cellstr(T.I_Morph(sel_spk)));
                uniqVals = sort(unique(morphVal));
                for L = 1:length(uniqVals)
                    is_level = sel_spk & (cellfun(@(x) str2double(x(2:end)), cellstr(T.I_Morph)) == uniqVals(L));
                    ilevel_all(is_level) = L;
                end
            end

            % --- Regressors (with time shift) - SUPER-BLOCKS ---
            for spkIdx = 1:2
                spk = speakers_to_split{spkIdx};
                for super_blk = 1:n_super_blocks
                    
                    % Map super-block to original blocks
                    if super_blk == 1
                        orig_blocks = [1, 2];
                    elseif super_blk == 2
                        orig_blocks = [3, 4];
                    else
                        orig_blocks = [5, 6];
                    end
                    
                    for lvl = 1:n_levels
                        idx = find(strcmp(T.Task, taskName) & strcmp(T.Speaker, spk) & ...
                                   ilevel_all == lvl & ismember(T.Block, orig_blocks));
                                   
                        reg_count = reg_count + 1;
                        names{reg_count} = sprintf('%s_%s_%s_lvl%d_blk%d', taskName, cond, spk, lvl, super_blk);
                        
                        if isempty(idx)
                            onsets{reg_count} = []; durations{reg_count} = [];
                        else
                            o = ((onset_times(idx) - bsl) / 10000) + cum_time;
                            onsets{reg_count} = o(:);
                            durations{reg_count} = zeros(numel(o), 1);
                        end
                    end
                end
            end

            % Cue Regressor (Original 6 Blocks)
            reg_count = reg_count + 1;
            names{reg_count} = sprintf('cue_%s', cond);
            cue_idx = false(height(T),1);
            for blk = 1:6 % Hart auf 6 setzen für alle 6 originalen Cues
                task_blk_idx = find(strcmp(T.Task, taskName) & T.Block == blk);
                if ~isempty(task_blk_idx), [~, min_i] = min(onset_times(task_blk_idx)); cue_idx(task_blk_idx(min_i)) = true; end
            end
            
            c_o = ((onset_times(cue_idx) - bsl) / 10000) + cum_time;
            onsets{reg_count} = c_o(:); 
            durations{reg_count} = zeros(numel(c_o),1);

            cum_time = cum_time + (nV * TR_sec);
        end

        % SAVE
        task_design_dir = fullfile(DESIGNDIR, ['sub-', subj_id], sprintf('task-%s', taskName));
        if ~exist(task_design_dir, 'dir'), mkdir(task_design_dir); end
        save(fullfile(task_design_dir, sprintf('design_%s_bigblocks_concat.mat', taskName)), 'names', 'onsets', 'durations');
    end

    %% --- STEP 2: CONCATENATED MEGA-SESSION GLM FITTING ---
    scans_cat = {}; mot18 = []; TR_sec = [];
    col_start = [1 7 13];

    for c = 1:3
        cond = stim_cond{c};
        sesdir = fullfile(subj_folder, ['ses-' cond]);
        vols = getVols(sesdir, '^wuabold-task_unwarped\.nii$');
        nV = size(vols, 1);

        if nV == 0
            continue; 
        end

        scans_cat = [scans_cat; cellstr(vols)];

        if isempty(TR_sec), TR_sec = spm_vol(vols(1,:)).private.timing.tspace; end

        % Build 18-column motion matrix (6 params * 3 runs)
        rpfile = dir(fullfile(sesdir, 'rp_abold-task_unwarped.txt'));
        if isempty(rpfile), R = zeros(nV, 6); 
        else, R = dlmread(fullfile(rpfile.folder, rpfile.name)); end

        tmp = zeros(nV, 18);
        tmp(:, col_start(c):col_start(c)+5) = R(1:min(size(R,1), nV), :);
        mot18 = [mot18; tmp];
    end

    % Save concatenated motion file
    reg_file = fullfile(subj_folder, 'rp_18_bigblocks_March2026.txt');
    dlmwrite(reg_file, mot18, 'delimiter', '\t', 'precision', '%.6f');

    % Loop for Task Estimation
    for t = 1:numel(task_list)
        task = task_list{t};
        dsgn = fullfile(DESIGNDIR, ['sub-' subj_id], sprintf('task-%s', task), sprintf('design_%s_bigblocks_concat.mat', task));

        if ~exist(dsgn, 'file')
            warning('Design file not found: %s', dsgn);
            continue;
        end

        outdir = fullfile(GLM_BASE, ['task-' task]);
        if ~exist(outdir, 'dir'), mkdir(outdir); end

        % --- SPM BATCH ---
        matlabbatch = {};
        matlabbatch{1}.spm.stats.fmri_spec.dir = {outdir};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR_sec;
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans_cat;
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {dsgn};
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {reg_file};
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model spec', substruct('.','val', '{}',{1},'.','val', '{}',{1},'.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1; % Required for RSA whitening

        if ~isempty(scans_cat)
            spm_jobman('run', matlabbatch);
        end
    end
    fprintf('GLM MARCH2026 (BIG BLOCKS CONCAT) Complete: sub-%s\n', subj_id);
end