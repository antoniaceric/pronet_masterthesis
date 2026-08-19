function pronet_psychometric_setup(VP)
global HPC_PATH DATA_PATH

% still needs:    subi = 1:length(VP)
%           subj = VP{subi};
%           subjID = subj(1:2); check if this works properly

% This script summarizes behavioral performance for the PRONET **intonation**
% and **phoneme** tasks based on the raw logfile data from each participant.
% It extracts accuracy, reaction times, trial counts, and missing responses
% per morph level and per speaker (m1, f2), then saves the results in
% MATLAB structs for later analysis (e.g., psychometric curve fitting or
% RSA behavioral model construction).
%
% === INPUTS ===
% 1. Subject lists:
%       /SCRIPTS/NEW_2024/ses-lt
%       /SCRIPTS/NEW_2024/ses-rt
%       /SCRIPTS/NEW_2024/ses-vx
%    Each file contains one subject ID per line (e.g. "sub-02").
%
% 2. Behavioral logfiles:
%       /DATA/proc/sub-XX/behav/*.txt
%    Each tab-delimited logfile contains columns (no header):
%       SubjID, Session, Block, Trial, Task, Button, RT, SoundTime,
%       Response, SOA, Speaker, P_Morph, P_Level, I_Morph, I_Level, Wavefile
%
% === PROCESSING STEPS ===
% Reads subject lists and loops over all participants.
% For each participant and logfile:
%     - Reads trialwise behavioral data.
%     - Corrects intonation morph levels (I_Level) per speaker.
%     - Separates data by task ('intonat' vs 'phoneme').
%     - Computes:
%         Accuracy per speaker × morph level
%         Mean reaction times per speaker × morph level
%         Number of correct and total trials
%         Number of missed trials
%     - Assigns session condition ('lt', 'rt', or 'vx') from filename.
%
% === OUTPUTS ===
% Saves two summary files to:
%       /DATA/proc/psychometrics/
%     summary_intonation_task.mat  (struct array INTONT)
%     summary_phoneme_task.mat     (struct array PHONEME)
%
% Each struct element corresponds to one subject × logfile and contains:
%     .subj       Subject ID
%     .cond       Condition ('lt','rt','vx')
%     .acc        Accuracy per speaker × level
%     .acc_all    Combined accuracy across speakers
%     .RT         Reaction times per speaker × level (cell arrays)
%     .corr       Number of correct responses
%     .num        Number of valid trials
%     .miss       Number of missed trials
%     .absLevel   Absolute morph level per speaker × level

stim_cond = {'lt', 'rt', 'vx'};
subjID = VP

OUTDIR = fullfile(HPC_PATH, 'DATA','proc','psychometrics');
INTONT = struct([]);
PHONEME = struct([]);
BEHAV_DIR = fullfile(DATA_PATH, ['sub-', subjID], 'behav');
    
    % List all txt logfiles for this subject
    behav_list = dir(fullfile(BEHAV_DIR, '*.txt'));
    logfiles = {behav_list.name};
    
    for si = 1:length(logfiles)
        logfile = fullfile(behav_list(si).folder, logfiles{si});
        
        % Read data table once per logfile
        T = readtable(logfile, 'delimiter', '\t', 'ReadVariableNames', false);

        new_header = {'SubjID','Session', 'Block','Trial', 'Task', 'Button', ...
          'RT', 'SoundTime', 'Response', 'SOA', 'Speaker', ...
          'P_Morph', 'P_Level', 'I_Morph', 'I_Level', 'Wavefile'};

        T.Properties.VariableNames = new_header;
        
        % Extract variables once
        id_all        = T.SubjID;
        ses_all       = T.Session;
        blc_all       = T.Block;
        trl_all       = T.Trial;
        task_all      = T.Task;
        button_all    = T.Button;
        rt_all        = T.RT;
        resp_all      = T.Response;
        speaker_all   = T.Speaker;
        pmorph_all    = T.P_Morph;
        plevel_all    = T.P_Level;
        imorph_all    = T.I_Morph;
        ilevel_all    = T.I_Level;
        wavname_all   = T.Wavefile;

        % Correct ilevel_all as you do in your original code
        ilevel_all_corrected = zeros(size(imorph_all));
        speakers = {'m1', 'f2'};
        for spkIdx = 1:length(speakers)
            spk = speakers{spkIdx};
            spk_idx = strcmp(speaker_all, spk);
            morph_strs = imorph_all(spk_idx);
            morph_vals = cellfun(@(x) str2double(x(2:end)), morph_strs);
            unique_morphs = sort(unique(morph_vals));
            for level = 1:length(unique_morphs)
                morph_val = unique_morphs(level);
                match_idx = spk_idx & cellfun(@(x) str2double(x(2:end)) == morph_val, imorph_all);
                ilevel_all_corrected(match_idx) = level;
            end
        end
        ilevel_all = ilevel_all_corrected;

        % Use these for filtering below:
        task = task_all;
        ilevel = ilevel_all;
        resp = resp_all;
        speaker = speaker_all;
        plevel = plevel_all;
        rt = rt_all;
        imorph = imorph_all;
        pmorph = pmorph_all;

        % Extract condition from filename:
        chunks = split(logfiles{si}, '_');
        cond = chunks{3}(1:2);  % adjust as needed

        %% INTONATION TASK PROCESSING
        Tasks = {'intonat'};
        Speakers = {'m1','f2'};
        nLevel = 5;
        
        p = zeros(length(Speakers), nLevel);
        RT = cell(length(Speakers), nLevel);
        yes = zeros(length(Speakers), nLevel);
        n = zeros(length(Speakers), nLevel);
        mo = zeros(length(Speakers), nLevel);
        miss = zeros(length(Speakers), nLevel);

        for spki = 1:length(Speakers)
            spk = Speakers{spki};

            ii = strcmp(task, 'intonat') & ilevel > 0 & strcmp(speaker, spk) & ~strcmp(resp, 'miss');
            rrtt = rt(ii);
            lvl = ilevel(ii);
            val = strcmp('Frage', resp(ii));
            morph = imorph(ii);

            for l = 1:nLevel
                jj = find(lvl == l);
                if ~isempty(jj)
                    p(spki, l) = sum(val(jj)) / numel(jj);
                    RT{spki, l} = rrtt(jj);
                    yes(spki, l) = sum(val(jj));
                    n(spki, l) = numel(jj);
                    mo(spki, l) = sscanf(char(unique(morph(jj))), '%*c%d');
                end
            end

            % MISS trials
            ii = strcmp(task, 'intonat') & ilevel > 0 & strcmp(speaker, spk) & strcmp(resp, 'miss');
            lvl = ilevel(ii);
            for l = 1:nLevel
                jj = find(lvl == l);
                miss(spki, l) = numel(jj);
            end
        end
        
        INTONT(subi, si).subj = subjID;
        INTONT(subi, si).acc = p;
        INTONT(subi, si).acc_all = sum(yes, 1) ./ sum(n, 1);
        INTONT(subi, si).RT = RT;
        INTONT(subi, si).corr = yes;
        INTONT(subi, si).num = n;
        INTONT(subi, si).absLevel = mo;
        INTONT(subi, si).miss = miss;
        INTONT(subi, si).cond = cond;

        %% PHONEME TASK PROCESSING
        Tasks = {'phoneme'};
        p = zeros(length(Speakers), nLevel);
        RT = cell(length(Speakers), nLevel);
        yes = zeros(length(Speakers), nLevel);
        n = zeros(length(Speakers), nLevel);
        mo = zeros(length(Speakers), nLevel);
        miss = zeros(length(Speakers), nLevel);

        for spki = 1:length(Speakers)
             spk = Speakers{spki};

            ii = strcmp(task, 'phoneme') & plevel > 0 & strcmp(speaker, spk) & ~strcmp(resp, 'miss');
            rrtt = rt(ii);
            lvl = plevel(ii);
            val = strcmp('Paar', resp(ii));
            morph = pmorph(ii);

            for l = 1:nLevel
                jj = find(lvl == l);
                if ~isempty(jj)
                    p(spki, l) = sum(val(jj)) / numel(jj);
                    RT{spki, l} = rrtt(jj);
                    yes(spki, l) = sum(val(jj));
                    n(spki, l) = numel(jj);
                    mo(spki, l) = sscanf(char(unique(morph(jj))), '%*c%d');
                end
            end

            % Miss trials
            ii = strcmp(task, 'phoneme') & plevel > 0 & strcmp(speaker, spk) & strcmp(resp, 'miss');
            lvl = plevel(ii);
            for l = 1:nLevel
                jj = find(lvl == l);
                miss(spki, l) = numel(jj);
            end
        end

        PHONEME(subi, si).subj = subjID;
        PHONEME(subi, si).acc = p;
        PHONEME(subi, si).acc_all = sum(yes, 1) ./ sum(n, 1);
        PHONEME(subi, si).RT = RT;
        PHONEME(subi, si).corr = yes;
        PHONEME(subi, si).num = n;
        PHONEME(subi, si).absLevel = mo;
        PHONEME(subi, si).miss = miss;
        PHONEME(subi, si).cond = cond;

    end % logfile loop
    
    disp([subjID ' finished! (' num2str((subi/length(VP))*100, '%.1f') '%)']);
end % subject loop

% Final save
if ~exist(OUTDIR, 'dir')
    mkdir(OUTDIR);
end

save(fullfile(OUTDIR, 'summary_intonation_task.mat'), 'INTONT');
save(fullfile(OUTDIR, 'summary_phoneme_task.mat'), 'PHONEME');

