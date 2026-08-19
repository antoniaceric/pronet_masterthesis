%% Plot Psychometrics - Linear and sigmoid fits per condition and task 
%  PRONET
%
% Written by A. Ceric (Modified by Coding-Assistent)
% Date: 09.04.2026

clear; clc;

%% === Setup ===
HPC_PATH = '/mnt/beegfs/workspace/2024-0404-PRONET/';
DATA_PATH = fullfile(HPC_PATH, 'DATA/proc/psychometrics/');
OUTDIR = fullfile(DATA_PATH, 'plots','APRIL26');

if ~exist(OUTDIR, 'dir'), mkdir(OUTDIR); end

% Alle drei Bedingungen werden analysiert
conditions = {'vx','lt','rt'};
tasks = {'intonation', 'phoneme'};

% Ursprüngliche Task-Farben (werden für die Einzelplots und den neuen reinen VX-Plot genutzt)
colors = struct('intonation', [0.8 0 0], 'phoneme', [0 0 0.8]);

% Neue Farb-Definitionen für den kombinierten Plot am Ende
cond_colors = struct(...
    'lt', [0.60 0.30 0.75], ... % Purple
    'rt', [0.95 0.70 0.10], ... % Orange/Yellow
    'vx', [0.60 0.60 0.60]  ... % Grey
);

n_voices = 2;

load(fullfile(DATA_PATH, 'manual_exclusions.mat'));

x_range = linspace(1, 5, 100);
lw_individual = 1; lw_average = 3;

%% === Loop over Tasks ===
for t = 1:length(tasks)
    task = tasks{t};
    load(fullfile(DATA_PATH, [task '_resp_fit_all_conditions.mat']));
    AllFitInfo_task = AllFitInfo;
    
    % Struktur zum Speichern der Durchschnittskurven für den kombinierten Plot
    mean_curves = struct();

    % === Plot per condition ===
    for ci = 1:length(conditions)
        cond = conditions{ci};
        entries = AllFitInfo_task.(cond).separate;
        
        % Lade die spezifischen Ausschlüsse für die aktuelle Aufgabe und Bedingung
        exclude_ids = unique([manual.(task).all, manual.(task).(cond)]);

        % --- Init plot ---
        figure('Color', 'w', 'Position', [300, 300, 600, 400]); hold on;
        title([capitalize(task) ' - ' upper(cond)], 'FontSize', 14);
        set(gca, 'FontSize', 11);
        grid on;

        % --- Per-subject voicewise average curves ---
        subj_y = [];  % for average
        for i = 1:numel(entries)
            sid = entries(i).subject;
            if ismember(sid, exclude_ids), continue; end

            a = entries(i).sigfit.bias(:);
            b = entries(i).sigfit.slope(:);
            y1 = 1 ./ (1 + exp((a(1) - x_range) * b(1)));
            y2 = 1 ./ (1 + exp((a(2) - x_range) * b(2)));
            y_avg = (y1 + y2) / 2;

            subj_y(end+1, :) = y_avg; %#ok<AGROW>
            plot(x_range, y_avg, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', lw_individual);
        end

        % --- Mean curve across subjects ---
        if ~isempty(subj_y)
            % Berechne den Durchschnitt und plotte ihn
            mean_y = mean(subj_y, 1);
            plot(x_range, mean_y, '-', 'Color', colors.(task), 'LineWidth', lw_average);
            
            % *** Speichere den Durchschnitt für spätere Plots ***
            mean_curves.(cond) = mean_y;
        end

        % --- Labels ---
        if strcmp(task, 'intonation')
            ylabel('P("Question")');
            xlabel('5 Prosody Levels');
        else
            ylabel('P("Paar")');
            xlabel('5 Phoneme Levels');
        end

        % --- Save ---
        fname = fullfile(OUTDIR, sprintf('sigmoid_fit_%s_%s.png', task, cond));
        saveas(gcf, fname); close(gcf);
    end
    
    % === Plot average across ALL conditions in ONE figure ===
    
    figure('Color', 'w', 'Position', [350, 350, 600, 400]); hold on;
    title([capitalize(task) ' - Averages (LT, RT, VX)'], 'FontSize', 14);
    set(gca, 'FontSize', 11);
    grid on;
    
    plot_handles = [];
    legend_labels = {};
    
    % Gehe durch die gespeicherten Durchschnittskurven und plotte sie
    for ci = 1:length(conditions)
        cond = conditions{ci};
        if isfield(mean_curves, cond)
            % Nutze die oben definierten speziellen Farben für lt, rt, vx
            p = plot(x_range, mean_curves.(cond), '-', 'Color', cond_colors.(cond), 'LineWidth', lw_average);
            
            % Speichere Handle und Label für die Legende
            plot_handles(end+1) = p; %#ok<AGROW>
            legend_labels{end+1} = upper(cond); %#ok<AGROW>
        end
    end
    
    % Füge eine Legende hinzu, damit man die Farben zuordnen kann
    if ~isempty(plot_handles)
        legend(plot_handles, legend_labels, 'Location', 'best');
    end
    
    % --- Labels für kombinierten Plot ---
    if strcmp(task, 'intonation')
        ylabel('P("Question")');
        xlabel('5 Prosody Levels');
    else
        ylabel('P("Paar")');
        xlabel('5 Phoneme Levels');
    end

    xlim([1 5]);
    xticks(1:5);

    % --- Speichern des kombinierten Plots ---
    fname_combined = fullfile(OUTDIR, sprintf('sigmoid_fit_%s_ALL_CONDITIONS_AVG.png', task));
    saveas(gcf, fname_combined); close(gcf);
    
    
    % === NEUER BLOCK: Plot ONLY average for VX (ohne Einzellinien) als EPS ===
    if isfield(mean_curves, 'vx')
        figure('Color', 'w', 'Position', [400, 400, 600, 400]); hold on;
        title([capitalize(task) ' - VX Average Only'], 'FontSize', 14);
        set(gca, 'FontSize', 11);
        grid on;
        
        % Zeichne nur die gespeicherte Durchschnittskurve für vx (in der Task-Farbe)
        plot(x_range, mean_curves.vx, '-', 'Color', colors.(task), 'LineWidth', lw_average);
        
        % --- Labels ---
        if strcmp(task, 'intonation')
            ylabel('P("Question")');
            xlabel('5 Prosody Levels');
        else
            ylabel('P("Paar")');
            xlabel('5 Phoneme Levels');
        end

        % --- Speichern als EPS ---
        fname_vx_only = fullfile(OUTDIR, sprintf('sigmoid_fit_%s_VX_AVG_ONLY.svg', task));
        saveas(gcf, fname_vx_only, 'svg'); close(gcf);
    end
    
    % === NEUER BLOCK: Pretest-Style Figure für VX ===
    % Erstellt für jeden Task einen Plot im exakten Stil der Pretest-Bilder
    % (Graue gestrichelte Einzellinien, dicker farbiger Durchschnitt, feste Achsen)
    
    cond_to_plot = 'vx';
    
    % Prüfen, ob die Bedingung existiert
    if isfield(AllFitInfo_task, cond_to_plot)
        entries_vx = AllFitInfo_task.(cond_to_plot).separate;
        exclude_ids_vx = unique([manual.(task).all, manual.(task).(cond_to_plot)]);

        fig_style = figure('Name', sprintf('%s - %s (Reference Style)', capitalize(task), upper(cond_to_plot)), ...
                     'Color', 'w', 'Position', [200 200 800 600]); 
        hold on;

        valid_subj_y_vx = [];
        
        % 1. Einzellinien (Hintergrund)
        for i = 1:numel(entries_vx)
            sid = entries_vx(i).subject;
            if ismember(sid, exclude_ids_vx), continue; end

            a = entries_vx(i).sigfit.bias(:);
            b = entries_vx(i).sigfit.slope(:);
            y1 = 1 ./ (1 + exp((a(1) - x_range) * b(1)));
            y2 = 1 ./ (1 + exp((a(2) - x_range) * b(2)));
            y_avg = (y1 + y2) / 2;

            valid_subj_y_vx(end+1, :) = y_avg; %#ok<AGROW>

            % Plot im Pretest-Stil: gestrichelt, grau, Linienbreite 1.2
            plot(x_range, y_avg, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
        end

        % 2. Durchschnittslinie (Vordergrund)
        if ~isempty(valid_subj_y_vx)
            mean_y_vx = mean(valid_subj_y_vx, 1);
            % Plot im Pretest-Stil: durchgezogen, Task-Farbe (Rot/Blau), Linienbreite 4
            plot(x_range, mean_y_vx, '-', 'Color', colors.(task), 'LineWidth', 4);
        end

        % 3. Formatierung exakt wie in den angehängten Bildern
        title(sprintf('%s - %s', capitalize(task), upper(cond_to_plot)), 'FontSize', 16, 'FontWeight', 'bold');
        
        if strcmp(task, 'intonation')
            ylabel('P("Question")', 'FontSize', 14);
            xlabel('5 Prosody Levels', 'FontSize', 14);
        else
            ylabel('P("Paar")', 'FontSize', 14);
            xlabel('5 Phoneme Levels', 'FontSize', 14);
        end

        ylim([0 1]); 
        xlim([1 5]);
        grid on; 
        box on;
        
        % Ticks manuell setzen (X-Achse nur ganze Zahlen)
        yticks(0:0.1:1);
        xticks(1:5);

        % Speichern
        fname_pretest_style = fullfile(OUTDIR, sprintf('sigmoid_fit_%s_%s_REFERENCE_STYLE.png', task, cond_to_plot));
        saveas(fig_style, fname_pretest_style);
        close(fig_style);
    end

end

%% Helper function
function out = capitalize(str)
    out = [upper(str(1)), str(2:end)];
end