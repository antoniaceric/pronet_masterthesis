function pronet_crossnobis_Feb2026(Job)
    % Computes Session-Specific & Speaker-Specific Crossnobis (lt, rt, vx | m1, f2)
    % AND automatically creates the Aggregated version (m1+f2 averaged).
    % Saves to: <DATA_PATH>/rsa_crossnobis_Speakerwise_Feb2026/

    VP = Job.VP;
    DATA_PATH = Job.DATA_PATH;
    task_index = Job.task_index;
    roi_index = Job.roi_index;

    %% SETUP
    subj_id = VP;
    subj_key = ['sub' subj_id];

    % Input File (Speakerwise Patterns)
    infile = fullfile(DATA_PATH, 'rsa_roi_patterns_bigblocks_March2026', sprintf('sub-%s_rsa_roi_pattern_bigblocks.mat', subj_id));
    
    % Output Directory (Separate from global)
    outdir = fullfile(DATA_PATH, 'rsa_crossnobis_bigblocks_March2026');
    if ~exist(outdir,'dir'), mkdir(outdir); end

    % --- Load Patterns ---
    if ~exist(infile, 'file')
        warning('Pattern file not found: %s', infile);
        return;
    end
    
    S = load(infile, 'roi_patterns', 'roi_groups');
    if ~isfield(S.roi_patterns, subj_key)
        warning('Subject %s not found in file', subj_key);
        return;
    end
    R          = S.roi_patterns;
    roi_names  = fieldnames(S.roi_groups);
    task_list  = {'phoneme','intonat'};
    sess_names = {'lt', 'rt', 'vx'}; 
    spk_names  = {'m1', 'f2'};

    % --- SPM Init ---
    spm('defaults','FMRI'); spm_jobman('initcfg');

    Sub = struct();
    Sub_Agg = struct(); % structure for averaged RDM

    % Select specific Task/ROI
    task = task_list{task_index};
    roi  = roi_names{roi_index};

    if isfield(R.(subj_key), task) && isfield(R.(subj_key).(task), roi)
        
        % Locate GLM
        glm_dir = fullfile(DATA_PATH, ['sub-' subj_id], ...
            'rsa_firstlevel_Speakerwise_Feb2026', ['task-' task]);

        if ~exist(glm_dir,'dir'), return; end

        res_list = spm_select('FPList', glm_dir, '^Res_.*\.nii$');
        if isempty(res_list), return; end

        % Get Data
        P = R.(subj_key).(task).(roi);

        Y = double(P.Y);              
        c = P.condVec(:);             
        p = P.partVec(:);             
        s_vec = P.sessVec(:); 
        spk_vec = P.spkVec(:);
        nvox = size(Y,2);
        
        if nvox >= 2
            % ====== Mask Logic ======
            mask_dir = '/home/antonia.ceric/git/pronet/data/proc/roi_masks';
            mask_file = fullfile(mask_dir, sprintf('HCPex_mask_%s.nii', roi));

            if ~exist(mask_file, 'file'), return; end

            SPM   = load(fullfile(glm_dir,'SPM.mat')); SPM = SPM.SPM;
            bfile = resolve_nii_path(SPM.Vbeta(P.idx(1)).fname, glm_dir);
            Vb    = spm_vol(bfile);
            Vm    = spm_vol(mask_file);

            % Reslice mask if geometry mismatches
            dim_mismatch = ~isequal(Vm.dim, Vb.dim);
            aff_mismatch = max(abs(Vm.mat(:) - Vb.mat(:))) > 1e-3;
            
            if dim_mismatch || aff_mismatch
                Vref = Vb; Vref.fname = fullfile(DATA_PATH, sprintf('tmp_ref_ses_%s_%s_%d.nii', subj_id, task, roi_index));
                spm_write_vol(Vref, zeros(Vref.dim));
                spm_reslice({Vref.fname, Vm.fname}, struct('interp',0,'which',1,'mean',0));
                [p_m,n_m,e_m] = fileparts(Vm.fname);
                Vm = spm_vol(fullfile(p_m, ['r' n_m e_m]));
                delete(Vref.fname);
            end

            mask3d   = spm_read_vols(Vm) > 0.5;
            mask_lin = find(mask3d);

            % ====== Whitening (Calculated on Global Data) ======
            Rres = read_residual_matrix(res_list, mask_lin);

            bad_nonfinite = any(~isfinite(Rres), 1);
            Rres = Rres - mean(Rres, 1, 'omitnan');
            zv   = std(Rres, 0, 1) < eps;
            drop = bad_nonfinite | zv; %drops NaNs because they do not go with crossnobis DISTANCE

            if any(drop)
                keep = ~drop;
                Rres = Rres(:, keep);
                Y    = Y(:, keep);
                nvox = size(Y,2);
            end
            
            if nvox >= 2
                % Covariance Shrinkage (1:1 ORIGINAL)
                Sig  = cov(Rres, 1);                  
                dSig = diag(diag(Sig));

                if any(~isfinite(Sig(:)))
                    sd = std(Rres, 0, 1); sd(sd < eps) = 1;
                    Y_white = Y ./ sd;
                else
                    varS    = mean((Sig(:) - mean(Sig(:))).^2);
                    varDiag = mean((diag(Sig) - diag(dSig)).^2);
                    lambda  = min(1, max(0, varS / (varS + varDiag + eps)));
                    Sig_sh  = (1 - lambda)*Sig + lambda*dSig;
                    Sig_sh  = (Sig_sh + Sig_sh')/2;
                    [V,D]   = eig(Sig_sh);
                    d       = diag(D);  d(d < 0) = 0;
                    Winvhalf = V * diag(1 ./ sqrt(d + eps)) * V';
                    Y_white  = Y * Winvhalf';
                end

                % % Demean per partition (Global)
                % up = unique(p(:))';
                % for pp = up
                %     idx = (p==pp);
                %     Y_white(idx,:) = Y_white(idx,:) - mean(Y_white(idx,:), 1, 'omitnan');
                % end

                % ====== Demean per partition AND Session ======
                up_sess = unique(s_vec(:))';
                for ss = up_sess
                    up_part = unique(p(:))';
                    for pp = up_part
                        % Jetzt wird strikt nach Block UND Session getrennt!
                        idx = (p == pp) & (s_vec == ss);
                        if sum(idx) > 0
                            Y_white(idx,:) = Y_white(idx,:) - mean(Y_white(idx,:), 1, 'omitnan');
                        end
                    end
                end

                % ====== SESSION-WISE & SPEAKER-WISE SPLIT ======
                for s = 1:3 % loop speakers
                    for sp = 1:2
                        s_name = sess_names{s};
                        spk_name = spk_names{sp};
                        
                        % Filter for this session AND speaker
                        idx_sess_spk = (s_vec == s) & (spk_vec == sp);
                        
                        if sum(idx_sess_spk) > 0
                            Ys = Y_white(idx_sess_spk, :);
                            cs = c(idx_sess_spk);
                            ps = p(idx_sess_spk);
                            
                            % Need at least 2 blocks for LDC
                            if numel(unique(ps)) >= 2
                                try
                                    dvec = rsa.distanceLDC(Ys, ps, cs, []);
                                    try, D = rsa.vector2RDM(dvec); catch, D = squareform(dvec); end
                                    
                                    % Save to nested structure
                                    Sub.(task).(roi).(s_name).(spk_name).RDM  = D;
                                    Sub.(task).(roi).(s_name).(spk_name).nVox = nvox;
                                catch
                                    % LDC failed
                                end
                            end
                        end
                    end
                    
                    %  ====== AGGREGATION ======
                    s_name = sess_names{s};
                    if isfield(Sub.(task).(roi), s_name)
                        sess_data = Sub.(task).(roi).(s_name);
                        % Average m1 and f2 if both exist
                        if isfield(sess_data, 'm1') && isfield(sess_data, 'f2')
                            if isfield(sess_data.m1, 'RDM') && isfield(sess_data.f2, 'RDM')
                                Sub_Agg.(task).(roi).(s_name).RDM = (sess_data.m1.RDM + sess_data.f2.RDM) / 2;
                                Sub_Agg.(task).(roi).(s_name).nVox = nvox;
                            end
                        end
                    end
                end
            end
        end
        
        %  ====== SAVING ======
        % 1. Speaker-wise
        fnMat_SPK = fullfile(outdir, sprintf('%s_task-%s_roi-%s_LDC_Speakerwise_bigblocks_March2026.mat', subj_key, task, roi));
        save(fnMat_SPK, 'Sub', 'VP', 'task', 'roi');
        fprintf('Saved Speaker-wise: %s\n', fnMat_SPK);

        % 2. Aggregated (Rename to 'Sub' so downstream scripts read it natively)
        if isfield(Sub_Agg, task)
            Sub_Backup = Sub;
            Sub = Sub_Agg; 
            
            fnMat_AGG = fullfile(outdir, sprintf('%s_task-%s_roi-%s_LDC_Aggregated_bigblocks_March2026.mat', subj_key, task, roi));
            save(fnMat_AGG, 'Sub', 'VP', 'task', 'roi');
            fprintf('Saved Aggregated:   %s\n', fnMat_AGG);
            
            Sub = Sub_Backup;
        end
    end
end

% ===== Helpers =====
function R = read_residual_matrix(res_list, mask_lin)
    nT  = size(res_list,1);
    V0  = spm_vol(deblank(res_list(1,:)));
    dim = V0.dim; nv = numel(mask_lin);
    R = nan(nT, nv, 'single');
    for t=1:nT
        V = spm_vol(deblank(res_list(t,:)));
        Y = spm_read_vols(V);
        R(t,:) = single(Y(mask_lin));
    end
end

function p = resolve_nii_path(fname, base_dir)
    f = char(fname); f = regexprep(f, ',\d+$', '');
    if isfile(f), p = f; return; end
    cand = fullfile(base_dir, f);
    if isfile(cand), p = cand; return; end
    [~,nm,ext] = fileparts(f);
    if isempty(ext)
        cand = fullfile(base_dir,[nm '.nii']); if isfile(cand), p = cand; return; end
        cand = fullfile(base_dir,[nm '.img']); if isfile(cand), p = cand; return; end
    elseif any(strcmpi(ext,{'.img','.hdr'}))
        cand = fullfile(base_dir,[nm '.img']); if isfile(cand), p = cand; return; end
        cand = fullfile(base_dir,[nm '.nii']); if isfile(cand), p = cand; return; end
    elseif strcmpi(ext,'.nii')
        cand = fullfile(base_dir,[nm '.nii']); if isfile(cand), p = cand; return; end
        cand = fullfile(base_dir,[nm '.img']); if isfile(cand), p = cand; return; end
    end
    error('Missing beta file: "%s"', fname);
end