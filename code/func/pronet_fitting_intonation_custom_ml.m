function pronet_fitting_intonation_custom_ml()
global HPC_PATH DATA_PATH

% === OVERVIEW ===
% Fits psychometric data using a custom Maximum Likelihood (ML) estimation 
% for a Logistic function, bypassing the Palamedes toolbox.
% Writes condition-wise fit summaries plus an all-conditions bundle to a 
% new specific folder: psychometric_fits_062026

%% SETUP
HPC_PATH  = '/mnt/beegfs/workspace/2024-0404-PRONET/';
DATA_PATH = fullfile(HPC_PATH,'DATA','proc');
OUTDIR    = fullfile(HPC_PATH,'DATA','proc','psychometric_fits_062026');

if ~exist(OUTDIR, 'dir')
    mkdir(OUTDIR);
end

conditions = {'vx', 'rt', 'lt'};
n_conditions = length(conditions);
voices = {'m1', 'f2'};
n_voices = length(voices);

% CRITICAL: ML requires integer trial counts, not just accuracy proportions. 
% Set the total number of trials per level here.
N_TRIALS_PER_LEVEL = 10; % <--- UPDATE THIS TO YOUR ACTUAL EXPERIMENT N

%% load behavioral data
load(fullfile(DATA_PATH, 'psychometrics', 'summary_intonation_task.mat'))

AllFitInfo = struct();  % init all conditions

for ci = 1:n_conditions
    cond = conditions{ci};

    % select current conditions entries
    cond_idx = arrayfun(@(s) strcmp(s.cond, cond), INTONT);
    INTONT_cond = INTONT(cond_idx);

    % Initialization
    FitInfo_separate = struct();  % Per Voice
    FitInfo_avg = struct();       % Averaged across voices

    for i = 1:length(INTONT_cond)
        subject = INTONT_cond(i).subj;
        ys = INTONT_cond(i).acc; 
        x = 1:size(ys, 2);
        
        n_array = repmat(N_TRIALS_PER_LEVEL, 1, length(x));

        %% ---- Fits per voice---- %%
        mlfit.alpha  = zeros(n_voices, 1); % Threshold
        mlfit.beta   = zeros(n_voices, 1); % Slope
        mlfit.gamma  = zeros(n_voices, 1); % Guess rate
        mlfit.lambda = zeros(n_voices, 1); % Lapse rate
        mlfit.LL     = zeros(n_voices, 1); % Log-Likelihood

        for vi = 1:n_voices
            y_acc = ys(vi,:);
            
            % Convert accuracy proportion back to "number of YES responses"
            yes_array = round(y_acc .* N_TRIALS_PER_LEVEL); 
            
            % Custom ML Fit
            [a, b, g, l, LL] = fit_logistic_ML(x, yes_array, n_array);
            
            mlfit.alpha(vi)  = a;
            mlfit.beta(vi)   = b;
            mlfit.gamma(vi)  = g;
            mlfit.lambda(vi) = l;
            mlfit.LL(vi)     = LL;
        end

        FitInfo_separate(i).subject = subject;
        FitInfo_separate(i).cond = cond;
        FitInfo_separate(i).x = x;
        FitInfo_separate(i).y = ys;
        FitInfo_separate(i).mlfit = mlfit;

        %% ---- Fit averaged across voices ---- %%
        y_avg_acc = mean(ys, 1);
        yes_avg_array = round(y_avg_acc .* N_TRIALS_PER_LEVEL);
        
        % Custom ML fit for average
        [a_avg, b_avg, g_avg, l_avg, LL_avg] = fit_logistic_ML(x, yes_avg_array, n_array);
        
        mlfit_avg.alpha  = a_avg;
        mlfit_avg.beta   = b_avg;
        mlfit_avg.gamma  = g_avg;
        mlfit_avg.lambda = l_avg;
        mlfit_avg.LL     = LL_avg;

        FitInfo_avg(i).subject = subject;
        FitInfo_avg(i).cond = cond;
        FitInfo_avg(i).x = x;
        FitInfo_avg(i).y = y_avg_acc;
        FitInfo_avg(i).mlfit = mlfit_avg;

        disp([subject ' - ' cond ' - ' num2str((i/length(INTONT_cond))*100, '%.1f') ' %'])
    end

    % Save separate voices
    fname_sep = fullfile(OUTDIR, ['intonation_resp_fit_' cond '_separate.mat']);
    save(fname_sep, 'FitInfo_separate');

    % Save averaged voices
    fname_avg = fullfile(OUTDIR, ['intonation_resp_fit_' cond '_average.mat']);
    save(fname_avg, 'FitInfo_avg');

    % Put all data
    AllFitInfo.(cond).separate = FitInfo_separate;
    AllFitInfo.(cond).average = FitInfo_avg;
end

% Save all
fname_all = fullfile(OUTDIR, 'intonation_resp_fit_all_conditions.mat');
save(fname_all, 'AllFitInfo');

%% === NEW SECTION: Generate and Save Fitted Curve Data Points ===
fprintf('\n--- Generating Fitted Curve Data Points (Custom ML) ---\n');

AllFitInfo_Fitted = AllFitInfo;
cond_fields = fieldnames(AllFitInfo_Fitted);

for c = 1:length(cond_fields)
    curr_cond = cond_fields{c};
    
    % Process 'separate' (per voice) data
    n_subs = length(AllFitInfo_Fitted.(curr_cond).separate);
    for s = 1:n_subs
        entry = AllFitInfo_Fitted.(curr_cond).separate(s);
        x_vals = entry.x; 
        
        y_fitted_matrix = zeros(size(entry.y));
        n_voices_local = size(entry.y, 1);
        
        for v = 1:n_voices_local
            alpha  = entry.mlfit.alpha(v);
            beta   = entry.mlfit.beta(v);
            gamma  = entry.mlfit.gamma(v);
            lambda = entry.mlfit.lambda(v);
            
            % Explicit logistic function equation
            y_fitted_matrix(v, :) = gamma + (1 - gamma - lambda) .* (1 ./ (1 + exp(-beta .* (x_vals - alpha))));
        end
        
        AllFitInfo_Fitted.(curr_cond).separate(s).y = y_fitted_matrix;
    end
    
    % Process 'average' data
    n_subs_avg = length(AllFitInfo_Fitted.(curr_cond).average);
    for s = 1:n_subs_avg
        entry = AllFitInfo_Fitted.(curr_cond).average(s);
        x_vals = entry.x;
        
        alpha  = entry.mlfit.alpha;
        beta   = entry.mlfit.beta;
        gamma  = entry.mlfit.gamma;
        lambda = entry.mlfit.lambda;
        
        y_fitted_vec = gamma + (1 - gamma - lambda) .* (1 ./ (1 + exp(-beta .* (x_vals - alpha))));
        AllFitInfo_Fitted.(curr_cond).average(s).y = y_fitted_vec;
    end
end

fname_fitted = fullfile(OUTDIR, 'intonation_resp_fit_all_conditions_FITTED.mat');
AllFitInfo = AllFitInfo_Fitted; 
save(fname_fitted, 'AllFitInfo'); 
fprintf('Saved custom ML fitted curve data to: %s\n', fname_fitted);

end

%% =======================================================================
%  LOCAL FUNCTION: Maximum Likelihood Estimation for Logistic Function
%  =======================================================================
function [alpha, beta, gamma, lambda, LL] = fit_logistic_ML(x, k, n)
    % Fixed parameters matching your old Palamedes setup
    gamma = 0.01;  % Guess rate
    lambda = 0.01; % Lapse rate
    
    % Initial guesses for optimization
    alpha_init = mean(x); % Threshold guess is the middle of the stimulus range
    beta_init = 1.0;      % Slope guess
    
    % Optimization options (suppress text output, ensure precision)
    options = optimset('Display', 'off', 'MaxFunEvals', 1000, 'MaxIter', 1000, 'TolX', 1e-9, 'TolFun', 1e-9);
    
    % Minimize the Negative Log-Likelihood
    [best_params, negLL] = fminsearch(@(params) neg_log_likelihood(params, x, k, n, gamma, lambda), [alpha_init, beta_init], options);
    
    alpha = best_params(1);
    beta = best_params(2);
    LL = -negLL; % Convert back to positive log-likelihood
end

function nll = neg_log_likelihood(params, x, k, n, gamma, lambda)
    alpha = params(1);
    beta = params(2);
    
    % Core Logistic Probability Function
    p = gamma + (1 - gamma - lambda) .* (1 ./ (1 + exp(-beta .* (x - alpha))));
    
    % Prevent log(0) errors by bounding probabilities slightly off 0 and 1
    p = max(min(p, 1 - 1e-10), 1e-10);
    
    % Binomial Log-Likelihood formula (ignoring the combinatorial constant)
    % LL = sum( k*log(p) + (n-k)*log(1-p) )
    nll = -sum(k .* log(p) + (n - k) .* log(1 - p)); 
end