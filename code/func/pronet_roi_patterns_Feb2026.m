function pronet_roi_patterns_Feb2026(Job)
% PRONET_ROI_PATTERNS_FEB2026
%
% 1. Loads SPM.mat from Speakerwise_Feb2026 GLM.
% 2. Filters for Speaker-specific Stimuli (m1, f2).
% 3. Extracts 180 clean betas (90 per speaker).
% 4. Saves Y data, plus vectors for Condition, Partition, Session, AND Speaker.

    VP = Job.VP;
    DATA_PATH = Job.DATA_PATH;
    hcp_dir   = Job.hcp_dir;
    subj_id   = VP;
    
    global HPC_PATH
    
    %% === Paths ===
    spm('defaults', 'FMRI'); spm_jobman('initcfg');
    
    hcp_atlas = fullfile(hcp_dir, 'HCPex_2mm.nii');
    label_mat = fullfile(hcp_dir, 'HCPex_LabelID.mat');
    
    assert(isfile(hcp_atlas), 'Missing HCPex atlas: %s', hcp_atlas);
    assert(isfile(label_mat), 'Missing label file: %s', label_mat);
    load(label_mat); 

    %% === ROI Groups ===
    roi_groups = struct();
    roi_groups.lIFG  = [159,160,165,166];
    roi_groups.lPMC  = [44,45,47,41];
    roi_groups.lPAC  = [54,55,56,57];
    roi_groups.lpSTS = [64,66];
    roi_groups.laSTS = [63,65];
    roi_groups.rIFG  = [339,340,345,346];
    roi_groups.rPMC  = [224,225,227,221];
    roi_groups.rPAC  = [234,235,236,237];
    roi_groups.rpSTS = [244,246];
    roi_groups.raSTS = [243,245];
    group_names = fieldnames(roi_groups);

    task_list = {'phoneme','intonat'}; 

    %% === Reference Beta Search (Updated Path) ===
    Vref = []; Vref_dim = []; Vref_aff = [];
    for t = 1:numel(task_list)
        task = task_list{t};
        % Path points to the new Feb2026 Speakerwise GLM
        glm_dir  = fullfile(DATA_PATH, ['sub-' subj_id], 'rsa_firstlevel_bigblocks_March2026', ['task-' task]);
        spm_file = fullfile(glm_dir, 'SPM.mat');
        if exist(spm_file, 'file')
            tmp_SPM = load(spm_file); tmp_SPM = tmp_SPM.SPM;
            if ~isempty(tmp_SPM.Vbeta)
                beta1 = resolve_nii_path(tmp_SPM.Vbeta(1).fname, glm_dir);
                if exist(beta1,'file')
                    Vref = spm_vol(beta1); Vref_dim = Vref.dim; Vref_aff = Vref.mat; break;
                end
            end
        end
    end
    assert(~isempty(Vref), 'No valid reference beta found.');

    %% === Atlas Reslicing ===
    Vatlas0 = spm_vol(hcp_atlas); Y0 = spm_read_vols(Vatlas0);
    if ~isequal(size(Y0), Vref_dim) || max(abs(Vatlas0.mat(:) - Vref_aff(:))) > 1e-3
        fprintf('[Atlas] Reslicing atlas to match beta space...\n');
        tmp = fullfile(DATA_PATH, sprintf('tmp_ref_%s.nii', subj_id)); 
        Vtmp = Vref; Vtmp.fname = tmp; spm_write_vol(Vtmp, zeros(Vref.dim));
        spm_reslice({tmp, hcp_atlas}, struct('interp',0,'which',1,'mean',0));
        delete(tmp);
        [p0,nm0,ex0] = fileparts(hcp_atlas);
        Vatlas = spm_vol(fullfile(p0, ['r' nm0 ex0]));
    else
        Vatlas = Vatlas0;
    end
    Y_atlas = spm_read_vols(Vatlas);

    %% === Main Extraction Loop ===
    roi_patterns = struct();
    meta = struct(); 
    fprintf('\n=== ROI Extraction: sub-%s (Speakerwise) ===\n', subj_id);
    subj_key = ['sub' subj_id];

    for t = 1:numel(task_list)
        task    = task_list{t};
        glm_dir = fullfile(DATA_PATH, ['sub-' subj_id], 'rsa_firstlevel_bigblocks_March2026', ['task-' task]);
        spm_file= fullfile(glm_dir, 'SPM.mat');
        
        if ~exist(spm_file, 'file'), continue; end
        load(spm_file); % Loads SPM variable
        
        all_desc = reshape({SPM.Vbeta.descrip}, [], 1);
        all_idx  = (1:numel(all_desc))';
        
        % Robust Filtering
        mask_task   = contains(lower(all_desc), lower(task));
        mask_no_cue = ~contains(lower(all_desc), 'cue');
        mask_stim   = contains(lower(all_desc), 'lvl') & contains(lower(all_desc), 'blk');
        
        keep_idx = all_idx(mask_task & mask_no_cue & mask_stim);
        nKept    = numel(keep_idx);

        fprintf('  Task %-8s | Found %d valid speaker-wise betas\n', task, nKept);
        
        % Now expecting 180 (2 speakers * 5 levels * 6 blocks * 3 runs)
        if nKept ~= 180
            warning('Expected 180 betas for speaker-wise analysis. Found %d.', nKept);
            if nKept == 0, continue; end
        end

        % % Now expecting 90 (2 speakers * 5 levels * 3 blocks * 3 runs)
        % if nKept ~= 90
        %     warning('Expected 90 betas for big-block analysis. Found %d.', nKept);
        %     if nKept == 0, continue; end
        % end

        % Parse Metadata (Updated Regex for Speaker)
        condVec = nan(nKept,1); % Level
        sessVec = nan(nKept,1); % Session (lt, rt, vx)
        partVec = nan(nKept,1); % Block
        spkVec  = nan(nKept,1); % Speaker (1=m1, 2=f2)
        
        % Updated Regex: matches task_session_speaker_lvlX_blkX
        % Example: intonat_lt_m1_lvl1_blk1
        parse_re = '([a-z]+)_(lt|rt|vx)_(m1|f2)_lvl(\d+)_blk(\d+)';
        
        for i = 1:nKept
            d = all_desc{keep_idx(i)};
            tok = regexp(d, parse_re, 'tokens', 'once'); 
            
            if isempty(tok)
                warning('Regex failed to find speaker/level info in: %s', d);
                continue;
            end
            
            ses_str = tok{2};
            spk_str = tok{3};
            lvl     = str2double(tok{4});
            blk     = str2double(tok{5});
            
            % Encode Session
            if strcmp(ses_str,'lt'), sID=1; elseif strcmp(ses_str,'rt'), sID=2; else, sID=3; end
            % Encode Speaker
            if strcmp(spk_str,'m1'), spkID=1; else, spkID=2; end
            
            condVec(i) = lvl;
            partVec(i) = blk;
            sessVec(i) = sID;
            spkVec(i)  = spkID;
        end

        % Extract Voxels
        for g = 1:numel(group_names)
            gname  = group_names{g};
            labels = roi_groups.(gname);
            mask   = ismember(Y_atlas, labels);
            nvox   = nnz(mask);
            
            if nvox == 0, continue; end

            Y = nan(nKept, nvox, 'single');
            ok_row = false(nKept,1);
            
            for ii = 1:nKept
                bfile = resolve_nii_path(SPM.Vbeta(keep_idx(ii)).fname, glm_dir);
                Vbi = spm_vol(bfile);
                Yi  = spm_read_vols(Vbi);
                
                if nnz(isfinite(Yi) & mask) > 0
                    Y(ii,:) = single(Yi(mask));
                    ok_row(ii) = true;
                end
            end
            
            % Save Structure
            roi_patterns.(subj_key).(task).(gname).Y        = Y(ok_row,:);
            roi_patterns.(subj_key).(task).(gname).condVec  = condVec(ok_row);
            roi_patterns.(subj_key).(task).(gname).partVec  = partVec(ok_row);
            roi_patterns.(subj_key).(task).(gname).sessVec  = sessVec(ok_row);
            roi_patterns.(subj_key).(task).(gname).spkVec   = spkVec(ok_row); % Added Speaker Vector
            roi_patterns.(subj_key).(task).(gname).idx      = keep_idx(ok_row);
            roi_patterns.(subj_key).(task).(gname).descrip  = all_desc(keep_idx(ok_row));
        end
        meta.(subj_key).(task).nBetasKept = nKept;
    end

    %% === Save (Updated Filename) ===
    outdir = fullfile(DATA_PATH, 'rsa_roi_patterns_bigblocks_March2026');
    if ~exist(outdir, 'dir'), mkdir(outdir); end
    outfile = fullfile(outdir, sprintf('sub-%s_rsa_roi_pattern_bigblocks.mat', subj_id));
    save(outfile, 'roi_patterns', 'meta', 'roi_groups', '-v7.3');
    fprintf('Saved Speaker-wise patterns: %s\n', outfile);
end