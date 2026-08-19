function pronet_plot_modelRDMs_Feb2026(VP, HPC_PATH, version_type)
% Plots Group Model RDMs for Feb 2026 versions.
% version_type: 'speakerwise' or 'aggregate'

    % Set Directory based on version
    if strcmpi(version_type, 'speakerwise')
        model_dir = fullfile(HPC_PATH, 'DATA', 'proc', 'model_rdms_pretest_Feb2026_speakerwisefit');
        suffix    = 'Speakerwise'; % UPDATE: Exakt so, wie es das neue Skript speichert
        sub_keys  = {'m1', 'f2', 'avg'}; % We will loop through these
    else
        model_dir = fullfile(HPC_PATH, 'DATA', 'proc', 'model_rdms_pretest_Feb2026_aggregatefit');
        suffix    = 'Aggregate'; % UPDATE: Großes 'A' wegen Linux Case-Sensitivity
        sub_keys  = {'aggregate'}; % Only one field here
    end

    out_dir = fullfile(model_dir, 'figs_group');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    tasks_keys   = {'phoneme', 'intonat'};
    tasks_titles = {'Phoneme Task', 'Intonation Task'};

    n_subs = length(VP);
    n_cols = ceil(sqrt(n_subs));
    n_rows = ceil(n_subs / n_cols);

    % Loop through Tasks
    for t = 1:numel(tasks_keys)
        curr_task  = tasks_keys{t};
        field_name = [curr_task '_pretest'];

        % Loop through the sub-fields (m1, f2, avg OR aggregate)
        for k = 1:numel(sub_keys)
            curr_sub_key = sub_keys{k};
            
            f = figure('Name', sprintf('Group Plot: %s (%s)', tasks_titles{t}, curr_sub_key), ...
                       'Color', 'w', 'Position', [50 50 1400 900], 'Visible', 'off');
            
            sgtitle(sprintf('%s Category RDMs: %s (%s)', version_type, tasks_titles{t}, curr_sub_key), ...
                    'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');

            for i = 1:n_subs
                subject_id = VP{i};
                sub_num = regexprep(subject_id, 'sub-', ''); 

                % Filename matches the new saving conventions
                if strcmpi(version_type,'speakerwise')
                    fname = sprintf('sub-%s_modelRDM_pretest_Speakerwise.mat', sub_num);
                else
                    fname = sprintf('sub-%s_modelRDM_pretest_Aggregate.mat', sub_num);
                end

                file_path = fullfile(model_dir, fname);
                
                if exist(file_path, 'file')
                    tmp = load(file_path);
                    % Check if task and speaker sub-key exists
                    if isfield(tmp, 'modelRDMs') && isfield(tmp.modelRDMs, field_name) && ...
                       isfield(tmp.modelRDMs.(field_name), curr_sub_key)

                        R = tmp.modelRDMs.(field_name).(curr_sub_key).category;

                        subplot(n_rows, n_cols, i);
                        imagesc(R);
                        axis square off; % Clean look
                        colormap(parula);
                        clim([0 1]); % Updated from caxis

                        title(subject_id, 'Interpreter', 'none', 'FontSize', 8);
                    end
                else
                    warning('File not found: %s', file_path); % Kleines Extra-Feedback
                end
            end
            
            % Save
            save_name = fullfile(out_dir, sprintf('GroupPlot_%s_%s_Category.png', curr_task, curr_sub_key));
            exportgraphics(f, save_name, 'Resolution', 300);
            close(f);
            
            fprintf(' -> Saved %s version (%s) plot: %s\n', version_type, curr_sub_key, save_name);
        end
    end
    fprintf('Done plotting %s RDMs.\n', version_type);
end