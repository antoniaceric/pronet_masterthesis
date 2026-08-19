function pronet_crossnobis_plot_Speakerwise_Feb2026(subj_id, data_path, varargin)
% Plots Speaker-wise & Aggregated Session-wise RDMs side-by-side for each ROI.
% Style: Matches original pronet_crossnobis_plot_sessionwise
% Usage: pronet_crossnobis_plot_bigblocks_March2026('02', DATA_PATH, 'save', true);

    p = inputParser;
    addRequired(p,'subj_id',@ischar);
    addRequired(p,'data_path',@ischar);
    addParameter(p,'save',false,@islogical);
    parse(p,subj_id,data_path,varargin{:});
    opts = p.Results;

    % Dirs - UPDATE AUF NEUEN ORDNER
    ldc_dir = fullfile(opts.data_path, 'rsa_crossnobis_bigblocks_March2026');
    fig_dir = fullfile(ldc_dir, 'figs');
    if opts.save && ~exist(fig_dir,'dir'), mkdir(fig_dir); end

    % Settings
    tasks    = {'phoneme', 'intonat'};
    rois     = {'lIFG', 'lPMC', 'lPAC', 'lpSTS', 'laSTS', ...
                'rIFG', 'rPMC', 'rPAC', 'rpSTS', 'raSTS'};
    sessions = {'lt', 'rt', 'vx'};
  
    for t = 1:numel(tasks)
        task = tasks{t};
        
        hf = figure('Name', sprintf('Sub-%s %s', subj_id, task), ...
            'Color','w', 'Position', [100, 100, 1800, 1500]); % Plot etwas breiter gemacht
        
        % 9 Spalten statt 6 (m1, f2, Agg für jede der 3 Sessions)
        tlo = tiledlayout(numel(rois), 9, 'TileSpacing','compact','Padding','compact');
        
        for r = 1:numel(rois)
            roi = rois{r};
            
            % Dateinamen definieren - UPDATE AUF NEUE NAMEN
            fname_spk = sprintf('sub%s_task-%s_roi-%s_LDC_Speakerwise_bigblocks_March2026.mat', subj_id, task, roi);
            fname_agg = sprintf('sub%s_task-%s_roi-%s_LDC_Aggregated_bigblocks_March2026.mat', subj_id, task, roi);
            fpath_spk = fullfile(ldc_dir, fname_spk);
            fpath_agg = fullfile(ldc_dir, fname_agg);
            
            % Default to empty 5x5 matrices (9 Stück)
            RDM_data = cell(1,9);
            for i=1:9, RDM_data{i} = nan(5,5); end
            
            % 1. Speakerwise Daten laden
            dat_spk = struct();
            if exist(fpath_spk, 'file')
                tmp_spk = load(fpath_spk, 'Sub');
                if isfield(tmp_spk.Sub, task) && isfield(tmp_spk.Sub.(task), roi)
                    dat_spk = tmp_spk.Sub.(task).(roi);
                end
            end
            
            % 2. Aggregated Daten laden
            dat_agg = struct();
            if exist(fpath_agg, 'file')
                tmp_agg = load(fpath_agg, 'Sub');
                if isfield(tmp_agg.Sub, task) && isfield(tmp_agg.Sub.(task), roi)
                    dat_agg = tmp_agg.Sub.(task).(roi);
                end
            end
            
            % Daten in das 1x9 Array mappen
            % Order: LT_m1, LT_f2, LT_agg, RT_m1, RT_f2, RT_agg, VX_m1, VX_f2, VX_agg
            count = 1;
            for s_idx = 1:3
                s_nm = sessions{s_idx};
                
                % m1
                if isfield(dat_spk, s_nm) && isfield(dat_spk.(s_nm), 'm1') && isfield(dat_spk.(s_nm).m1, 'RDM')
                    RDM_data{count} = dat_spk.(s_nm).m1.RDM;
                end
                count = count + 1;
                
                % f2
                if isfield(dat_spk, s_nm) && isfield(dat_spk.(s_nm), 'f2') && isfield(dat_spk.(s_nm).f2, 'RDM')
                    RDM_data{count} = dat_spk.(s_nm).f2.RDM;
                end
                count = count + 1;
                
                % Aggregated
                if isfield(dat_agg, s_nm) && isfield(dat_agg.(s_nm), 'RDM')
                    RDM_data{count} = dat_agg.(s_nm).RDM;
                end
                count = count + 1;
            end
            
            % Labels für die 9 Spalten
            labels = {'LT (m1)', 'LT (f2)', 'LT (Agg)', ...
                      'RT (m1)', 'RT (f2)', 'RT (Agg)', ...
                      'VX (m1)', 'VX (f2)', 'VX (Agg)'};
                      
            for i = 1:9
                plot_rdm_tile(RDM_data{i}, roi, labels{i});
            end
        end
        
        sgtitle(sprintf('Subject %s - Task: %s', subj_id, task), 'FontSize', 16, 'Interpreter', 'none');
        
        if opts.save
            % Dateiname für das gespeicherte Bild angepasst
            outname = fullfile(fig_dir, sprintf('sub-%s_task-%s_Combined_RDMs_bigblocks_March2026.png', subj_id, task));
            exportgraphics(hf, outname, 'Resolution', 300);
            fprintf('Saved: %s\n', outname);
            close(hf);
        end
    end
end

function plot_rdm_tile(D, roi_name, sess_label)
    nexttile;
    imagesc(D);
    axis square off;
    colormap(gca, parula);
    
    % Original scale color logic
    vals = D(~isnan(D));
    if isempty(vals)
        caxis([0 1]); 
    else
        % Check for completely flat matrices to avoid caxis error
        if min(vals) == max(vals)
            caxis([min(vals)-0.1, max(vals)+0.1]);
        else
            caxis([min(vals) max(vals)]); 
        end
    end
    
    % Title style
    title(sprintf('%s - %s', roi_name, sess_label), 'FontSize', 8, 'Interpreter','none');
end


%% Feb version % function pronet_crossnobis_plot_Speakerwise_Feb2026(subj_id, data_path, varargin)
% % Plots Speaker-wise & Aggregated Session-wise RDMs side-by-side for each ROI.
% % Style: Matches original pronet_crossnobis_plot_sessionwise
% % Usage: pronet_crossnobis_plot_Speakerwise_Feb2026('02', DATA_PATH, 'save', true);
% 
%     p = inputParser;
%     addRequired(p,'subj_id',@ischar);
%     addRequired(p,'data_path',@ischar);
%     addParameter(p,'save',false,@islogical);
%     parse(p,subj_id,data_path,varargin{:});
%     opts = p.Results;
% 
%     % Dirs
%     ldc_dir = fullfile(opts.data_path, 'rsa_crossnobis_Speakerwise_Feb2026');
%     fig_dir = fullfile(ldc_dir, 'figs');
%     if opts.save && ~exist(fig_dir,'dir'), mkdir(fig_dir); end
% 
%     % Settings
%     tasks    = {'phoneme', 'intonat'};
%     rois     = {'lIFG', 'lPMC', 'lPAC', 'lpSTS', 'laSTS', ...
%                 'rIFG', 'rPMC', 'rPAC', 'rpSTS', 'raSTS'};
%     sessions = {'lt', 'rt', 'vx'};
% 
%     for t = 1:numel(tasks)
%         task = tasks{t};
% 
%         hf = figure('Name', sprintf('Sub-%s %s', subj_id, task), ...
%             'Color','w', 'Position', [100, 100, 1800, 1500]); % Plot etwas breiter gemacht
% 
%         % [UPDATE]: 9 Spalten statt 6 (m1, f2, Agg für jede der 3 Sessions)
%         tlo = tiledlayout(numel(rois), 9, 'TileSpacing','compact','Padding','compact');
% 
%         for r = 1:numel(rois)
%             roi = rois{r};
% 
%             % Dateinamen definieren
%             fname_spk = sprintf('sub%s_task-%s_roi-%s_LDC_Speakerwise_Feb2026.mat', subj_id, task, roi);
%             fname_agg = sprintf('sub%s_task-%s_roi-%s_LDC_Aggregated_Feb2026.mat', subj_id, task, roi);
%             fpath_spk = fullfile(ldc_dir, fname_spk);
%             fpath_agg = fullfile(ldc_dir, fname_agg);
% 
%             % Default to empty 5x5 matrices (9 Stück)
%             RDM_data = cell(1,9);
%             for i=1:9, RDM_data{i} = nan(5,5); end
% 
%             % 1. Speakerwise Daten laden
%             dat_spk = struct();
%             if exist(fpath_spk, 'file')
%                 tmp_spk = load(fpath_spk, 'Sub');
%                 if isfield(tmp_spk.Sub, task) && isfield(tmp_spk.Sub.(task), roi)
%                     dat_spk = tmp_spk.Sub.(task).(roi);
%                 end
%             end
% 
%             % 2. Aggregated Daten laden
%             dat_agg = struct();
%             if exist(fpath_agg, 'file')
%                 tmp_agg = load(fpath_agg, 'Sub');
%                 if isfield(tmp_agg.Sub, task) && isfield(tmp_agg.Sub.(task), roi)
%                     dat_agg = tmp_agg.Sub.(task).(roi);
%                 end
%             end
% 
%             % Daten in das 1x9 Array mappen
%             % Order: LT_m1, LT_f2, LT_agg, RT_m1, RT_f2, RT_agg, VX_m1, VX_f2, VX_agg
%             count = 1;
%             for s_idx = 1:3
%                 s_nm = sessions{s_idx};
% 
%                 % m1
%                 if isfield(dat_spk, s_nm) && isfield(dat_spk.(s_nm), 'm1') && isfield(dat_spk.(s_nm).m1, 'RDM')
%                     RDM_data{count} = dat_spk.(s_nm).m1.RDM;
%                 end
%                 count = count + 1;
% 
%                 % f2
%                 if isfield(dat_spk, s_nm) && isfield(dat_spk.(s_nm), 'f2') && isfield(dat_spk.(s_nm).f2, 'RDM')
%                     RDM_data{count} = dat_spk.(s_nm).f2.RDM;
%                 end
%                 count = count + 1;
% 
%                 % Aggregated
%                 if isfield(dat_agg, s_nm) && isfield(dat_agg.(s_nm), 'RDM')
%                     RDM_data{count} = dat_agg.(s_nm).RDM;
%                 end
%                 count = count + 1;
%             end
% 
%             % [UPDATE]: Labels für die 9 Spalten
%             labels = {'LT (m1)', 'LT (f2)', 'LT (Agg)', ...
%                       'RT (m1)', 'RT (f2)', 'RT (Agg)', ...
%                       'VX (m1)', 'VX (f2)', 'VX (Agg)'};
% 
%             for i = 1:9
%                 plot_rdm_tile(RDM_data{i}, roi, labels{i});
%             end
%         end
% 
%         sgtitle(sprintf('Subject %s - Task: %s', subj_id, task), 'FontSize', 16, 'Interpreter', 'none');
% 
%         if opts.save
%             % Name angepasst, damit er alte Bilder nicht direkt überschreibt
%             outname = fullfile(fig_dir, sprintf('sub-%s_task-%s_Combined_RDMs.png', subj_id, task));
%             exportgraphics(hf, outname, 'Resolution', 300);
%             fprintf('Saved: %s\n', outname);
%             close(hf);
%         end
%     end
% end
% 
% function plot_rdm_tile(D, roi_name, sess_label)
%     nexttile;
%     imagesc(D);
%     axis square off;
%     colormap(gca, parula);
% 
%     % Original scale color logic
%     vals = D(~isnan(D));
%     if isempty(vals)
%         caxis([0 1]); 
%     else
%         % Check for completely flat matrices to avoid caxis error
%         if min(vals) == max(vals)
%             caxis([min(vals)-0.1, max(vals)+0.1]);
%         else
%             caxis([min(vals) max(vals)]); 
%         end
%     end
% 
%     % Title style
%     title(sprintf('%s - %s', roi_name, sess_label), 'FontSize', 8, 'Interpreter','none');
% end