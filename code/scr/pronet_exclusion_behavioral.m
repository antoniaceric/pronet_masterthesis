%% Exclude participants based on behavioral performance
%%
% written by A. Ceric
% 03.06.2025
%
% This script identifies participants to exclude due to insufficient
% classification performance in at least one task. Insufficient
% classification is defined as an average speaker-wise difference
% between minimal and maximal classification proportions < 0.6 in any task.

%% Setup

HPC_PATH = '/mnt/beegfs/workspace/2024-0404-PRONET/';
DATA_PATH = fullfile(HPC_PATH, 'DATA/proc/psychometrics/');
OUTDIR = fullfile(HPC_PATH, 'DATA/proc/psychometrics/plots_fitted');

% Load precomputed fit data
AllFitInfo_into = load(fullfile(DATA_PATH, 'intonation_resp_fit_all_conditions_FITTED.mat'));    
AllFitInfo_phon = load(fullfile(DATA_PATH, 'phoneme_resp_fit_all_conditions_FITTED.mat'));         

% Threshold for exclusion
min_diff_required = 0.6;

% Helper: avg. speaker-wise (max - min) difference
get_avg_diff = @(Y) mean(max(Y, [], 2) - min(Y, [], 2));

%% Prepare conditions list
conditions = {'vx', 'rt', 'lt'};

%% Initialize container for subjects
subjects = containers.Map();  % key: subject ID string

%% Process Intonation task: average over conditions
for c = 1:length(conditions)
    cond = conditions{c};
    info = AllFitInfo_into.AllFitInfo.(cond).separate;  % access loaded struct
    
    for i = 1:length(info)
        subj = info(i).subject;
        Y = info(i).y;
        avg_diff = get_avg_diff(Y);

        if ~isKey(subjects, subj)
            % Initialize structure with fields for intonation and phoneme diffs
            subjects(subj) = struct( ...
                'intonation_per_cond', containers.Map(), ...
                'phoneme_per_cond', containers.Map());
        end

        s = subjects(subj);
        % Store difference for this condition (intonation)
        s.intonation_per_cond(cond) = avg_diff;
        subjects(subj) = s;
    end
end

%% Process Phoneme task: average over conditions
for c = 1:length(conditions)
    cond = conditions{c};
    info = AllFitInfo_phon.AllFitInfo.(cond).separate;  % access loaded struct
    
    for i = 1:length(info)
        subj = info(i).subject;
        Y = info(i).y;
        avg_diff = get_avg_diff(Y);

        if ~isKey(subjects, subj)
            subjects(subj) = struct( ...
                'intonation_per_cond', containers.Map(), ...
                'phoneme_per_cond', containers.Map());
        end

        s = subjects(subj);
        % Store difference for this condition (phoneme)
        s.phoneme_per_cond(cond) = avg_diff;
        subjects(subj) = s;
    end
end

%% Compute average across conditions per subject for each task
subject_keys = keys(subjects);
excluded_subjects = {};

intonation_avg = zeros(length(subject_keys),1);
phoneme_avg = zeros(length(subject_keys),1);

for i = 1:length(subject_keys)
    subj = subject_keys{i};
    s = subjects(subj);
    
    % Convert maps to arrays
    intonation_vals = cell2mat(values(s.intonation_per_cond));
    phoneme_vals = cell2mat(values(s.phoneme_per_cond));
    
    % Average across conditions
    intonation_avg(i) = mean(intonation_vals);
    phoneme_avg(i) = mean(phoneme_vals);
    
    % Check exclusion
    if intonation_avg(i) < min_diff_required || phoneme_avg(i) < min_diff_required
        excluded_subjects{end+1} = subj;
    end
end

%% Report excluded subjects
fprintf('\nExcluded %d subject(s):\n', length(excluded_subjects));
disp(excluded_subjects');

% Optionally save excluded subjects
save(fullfile(DATA_PATH, 'excluded_subjects_avg.mat'), 'excluded_subjects');

%% Plot proportion differences

figure('Color', 'w', 'Position', [100 100 1000 500]);

n = length(subject_keys);
is_excluded = false(n,1);
for i = 1:n
    is_excluded(i) = ismember(subject_keys{i}, excluded_subjects);
end

hold on;

% Dummy bars for legend
h_included_into = bar(n+1, NaN, 0.3, 'FaceColor', [0 0.447 0.741], 'EdgeColor', 'none');  % Blue
h_included_phon = bar(n+1, NaN, 0.3, 'FaceColor', [0 0.7 0.2], 'EdgeColor', 'none');      % Green
h_excl_into =     bar(n+1, NaN, 0.3, 'FaceColor', [1 0 0], 'EdgeColor', 'none');          % Bright Red
h_excl_phon =     bar(n+1, NaN, 0.3, 'FaceColor', [0.6 0 0], 'EdgeColor', 'none');        % Dark Red

for i = 1:n
    if is_excluded(i)
        % Excluded: red tones
        bar(i - 0.15, intonation_avg(i), 0.3, 'FaceColor', [1 0 0], 'EdgeColor', 'none');      % Bright red
        bar(i + 0.15, phoneme_avg(i), 0.3, 'FaceColor', [0.6 0 0], 'EdgeColor', 'none');       % Dark red
    else
        % Included: blue and green
        bar(i - 0.15, intonation_avg(i), 0.3, 'FaceColor', [0 0.447 0.741], 'EdgeColor', 'none');  % Blue
        bar(i + 0.15, phoneme_avg(i), 0.3, 'FaceColor', [0 0.7 0.2], 'EdgeColor', 'none');         % Green
    end
end

% Threshold line
yline(min_diff_required, 'r--', 'LineWidth', 1.5, ...
    'Label', 'Threshold = 0.6', 'LabelHorizontalAlignment', 'left', ...
    'LabelVerticalAlignment', 'bottom');

xticks(1:n);
xticklabels(subject_keys);
xtickangle(45);
ylabel('Avg. Speaker-wise Proportion Difference');
title('Classification Proportion Differences per Subject');

% Add updated legend
legend([h_included_into, h_included_phon, h_excl_into, h_excl_phon], ...
       {'Included - Intonation', 'Included - Phoneme', ...
        'Excluded - Intonation', 'Excluded - Phoneme'}, ...
        'Location', 'northeast');

ylim([0 1]);
grid on;
box on;

% Save plot
saveas(gcf, fullfile(OUTDIR, 'subject_proportion_differences.png'));

%% -------- PER-CONDITION ANALYSIS -------- %%

for c = 1:length(conditions)
    cond = conditions{c};

    % Initialize storage
    subject_cond_info = containers.Map();  % key: subject ID

    % Intonation Task - this condition
    into_info = AllFitInfo_into.AllFitInfo.(cond).separate;
    for i = 1:length(into_info)
        subj = into_info(i).subject;
        Y = into_info(i).y;
        diff = get_avg_diff(Y);

        if ~isKey(subject_cond_info, subj)
            subject_cond_info(subj) = struct('intonation', NaN, 'phoneme', NaN);
        end
        s = subject_cond_info(subj);
        s.intonation = diff;
        subject_cond_info(subj) = s;
    end

    % Phoneme Task - this condition
    phon_info = AllFitInfo_phon.AllFitInfo.(cond).separate;
    for i = 1:length(phon_info)
        subj = phon_info(i).subject;
        Y = phon_info(i).y;
        diff = get_avg_diff(Y);

        if ~isKey(subject_cond_info, subj)
            subject_cond_info(subj) = struct('intonation', NaN, 'phoneme', NaN);
        end
        s = subject_cond_info(subj);
        s.phoneme = diff;
        subject_cond_info(subj) = s;
    end

    % Process data
    cond_subjects = keys(subject_cond_info);
    n_cond = length(cond_subjects);
    intonation_diff = zeros(n_cond, 1);
    phoneme_diff = zeros(n_cond, 1);
    cond_excluded = false(n_cond, 1);

    for i = 1:n_cond
        subj = cond_subjects{i};
        s = subject_cond_info(subj);
        intonation_diff(i) = s.intonation;
        phoneme_diff(i) = s.phoneme;
        cond_excluded(i) = (s.intonation < min_diff_required || s.phoneme < min_diff_required);
    end

    %% Plot for this condition
    figure('Color', 'w', 'Position', [100 100 1000 500]);
    hold on;

    % Legend bars (dummy)
    h_included_into = bar(n_cond+1, NaN, 0.3, 'FaceColor', [0 0.447 0.741], 'EdgeColor', 'none');  % Blue
    h_included_phon = bar(n_cond+1, NaN, 0.3, 'FaceColor', [0 0.7 0.2], 'EdgeColor', 'none');      % Green
    h_excl_into =     bar(n_cond+1, NaN, 0.3, 'FaceColor', [1 0 0], 'EdgeColor', 'none');          % Bright Red
    h_excl_phon =     bar(n_cond+1, NaN, 0.3, 'FaceColor', [0.6 0 0], 'EdgeColor', 'none');        % Dark Red

    for i = 1:n_cond
        if cond_excluded(i)
            bar(i - 0.15, intonation_diff(i), 0.3, 'FaceColor', [1 0 0], 'EdgeColor', 'none');   % Bright red
            bar(i + 0.15, phoneme_diff(i), 0.3, 'FaceColor', [0.6 0 0], 'EdgeColor', 'none');    % Dark red
        else
            bar(i - 0.15, intonation_diff(i), 0.3, 'FaceColor', [0 0.447 0.741], 'EdgeColor', 'none');  % Blue
            bar(i + 0.15, phoneme_diff(i), 0.3, 'FaceColor', [0 0.7 0.2], 'EdgeColor', 'none');         % Green
        end
    end

    yline(min_diff_required, 'r--', 'LineWidth', 1.5, ...
        'Label', 'Threshold = 0.6', 'LabelHorizontalAlignment', 'left', ...
        'LabelVerticalAlignment', 'bottom');

    xticks(1:n_cond);
    xticklabels(cond_subjects);
    xtickangle(45);
    ylabel('Speaker-wise Proportion Difference');
    title(sprintf('Proportion differences for condition: %s', upper(cond)));

    legend([h_included_into, h_included_phon, h_excl_into, h_excl_phon], ...
        {'Included - Intonation', 'Included - Phoneme', ...
         'Excluded - Intonation', 'Excluded - Phoneme'}, ...
         'Location', 'northeastoutside', 'Box', 'off');

    ylim([0 1]);
    grid off;
    box off;

    % Save plot
    fname = fullfile(OUTDIR, ['subject_diff_condition_' cond '.png']);
    saveas(gcf, fname);
end

%% ---- Manual exclusion lists (23-Jun-2025) ----
% manual.phoneme.all = {'05','11','24','27','32'};
% manual.phoneme.rt  = {'26'};
% manual.phoneme.vx  = {'29'};
% manual.phoneme.lt  = {};
% 
% manual.intonation.all = {'13'};
% manual.intonation.rt  = {'26'};
% manual.intonation.vx  = {'29'};
% manual.intonation.lt  = {'32'};
% 
% save(fullfile(DATA_PATH,'manual_exclusions.mat'),'manual');

