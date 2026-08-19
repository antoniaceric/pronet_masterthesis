function pronet_fitting_intonation_comparison()
global HPC_PATH DATA_PATH

%% SETUP PATHS
INDIR = fullfile(DATA_PATH, 'psychometrics');
OUTDIR = fullfile(INDIR, 'fitting_comparison');

if ~exist(OUTDIR, 'dir')
    mkdir(OUTDIR);
    fprintf('Erstelle Output-Ordner: %s\n', OUTDIR);
end

%% CONFIG
conditions = {'vx', 'rt', 'lt'};
n_conditions = length(conditions);
voices = {'m1', 'f2'};
n_voices = length(voices);

%% LOAD DATA
load(fullfile(INDIR, 'summary_intonation_task.mat'), 'INTONT');

AllFitInfo = struct(); 

for ci = 1:n_conditions
    cond = conditions{ci};
    cond_idx = arrayfun(@(s) strcmp(s.cond, cond), INTONT);
    INTONT_cond = INTONT(cond_idx);

    FitInfo_separate = struct();  
    FitInfo_avg = struct();       

    for i = 1:length(INTONT_cond)
        subject = INTONT_cond(i).subj;
        ys = INTONT_cond(i).acc;
        
        % HIER: Deine 5 Steps (1 bis 5)
        x = 1:size(ys, 2);

        %% ---- Fits per voice ---- %%
        linfit = init_fit_struct(n_voices);
        sigfit2 = init_fit_struct(n_voices); 
        sigfit4 = init_fit_struct(n_voices); 
        sigfit4.guess_rate = zeros(n_voices, 1);
        sigfit4.lapse_rate = zeros(n_voices, 1);

        for vi = 1:n_voices
            y = ys(vi,:);
            [xData, yData] = prepareCurveData(x, y);

            % Linear
            [linfit.intercept(vi), linfit.slope(vi), linfit.rsq(vi), linfit.adjrsq(vi)] = run_linear_fit(xData, yData);
            
            % Sigmoid 2-Param (Baek Logic)
            [params, gof] = run_sigmoid_2param(xData, yData);
            sigfit2.bias(vi) = params(1); 
            sigfit2.slope(vi) = params(2); 
            sigfit2.rsq(vi) = gof.rsquare; 
            sigfit2.adjrsq(vi) = gof.adjrsquare;
            
            % Sigmoid 4-Param (Baek Logic)
            [params, gof] = run_sigmoid_4param(xData, yData);
            sigfit4.guess_rate(vi) = params(1);       
            sigfit4.lapse_rate(vi) = 1 - params(2);   
            sigfit4.bias(vi) = params(3);             
            sigfit4.slope(vi) = params(4);            
            sigfit4.rsq(vi) = gof.rsquare; 
            sigfit4.adjrsq(vi) = gof.adjrsquare;
        end

        FitInfo_separate(i).subject = subject;
        FitInfo_separate(i).cond = cond;
        FitInfo_separate(i).x = x;
        FitInfo_separate(i).y = ys;
        FitInfo_separate(i).linfit = linfit;
        FitInfo_separate(i).sigfit2 = sigfit2;
        FitInfo_separate(i).sigfit4 = sigfit4;

        %% ---- Fit averaged across voices ---- %%
        y_avg = mean(ys, 1);
        [xData, yData] = prepareCurveData(x, y_avg);

        % Linear
        [linfit_avg.intercept, linfit_avg.slope, linfit_avg.rsq, linfit_avg.adjrsq] = run_linear_fit(xData, yData);
        
        % Sigmoid 2-Param
        [params, gof] = run_sigmoid_2param(xData, yData);
        sigfit2_avg.bias = params(1); sigfit2_avg.slope = params(2);
        sigfit2_avg.rsq = gof.rsquare; sigfit2_avg.adjrsq = gof.adjrsquare;
        
        % Sigmoid 4-Param
        [params, gof] = run_sigmoid_4param(xData, yData);
        sigfit4_avg.guess_rate = params(1); sigfit4_avg.lapse_rate = 1 - params(2);
        sigfit4_avg.bias = params(3); sigfit4_avg.slope = params(4);
        sigfit4_avg.rsq = gof.rsquare; sigfit4_avg.adjrsq = gof.adjrsquare;

        FitInfo_avg(i).subject = subject;
        FitInfo_avg(i).cond = cond;
        FitInfo_avg(i).x = x;
        FitInfo_avg(i).y = y_avg;
        FitInfo_avg(i).linfit = linfit_avg;
        FitInfo_avg(i).sigfit2 = sigfit2_avg;
        FitInfo_avg(i).sigfit4 = sigfit4_avg;

        disp([subject ' - ' cond ' - ' num2str((i/length(INTONT_cond))*100, '%.0f') '%'])
    end

    % SAVE
    save(fullfile(OUTDIR, ['intonation_resp_fit_' cond '_separate_comparison.mat']), 'FitInfo_separate');
    save(fullfile(OUTDIR, ['intonation_resp_fit_' cond '_average_comparison.mat']), 'FitInfo_avg');

    AllFitInfo.(cond).separate = FitInfo_separate;
    AllFitInfo.(cond).average = FitInfo_avg;
end

save(fullfile(OUTDIR, 'intonation_resp_fit_all_conditions_comparison.mat'), 'AllFitInfo');
end

%% === HELPER FUNCTIONS (Baek Logic applied to 5 Steps) ===

function s = init_fit_struct(n)
    s.intercept = zeros(n,1); s.slope = zeros(n,1); s.rsq = zeros(n,1); s.adjrsq = zeros(n,1);
end

function [intercept, slope, rsq, adjrsq] = run_linear_fit(x, y)
    [fitresult, gof] = fit(x, y, fittype('poly1'), fitoptions('poly1'));
    p = coeffvalues(fitresult); slope = p(1); intercept = p(2); rsq = gof.rsquare; adjrsq = gof.adjrsquare;
end

function [params, gof] = run_sigmoid_2param(x, y)
    ft = fittype('0.01 + 0.98 ./ (1+exp((a-x)*b))', 'independent', 'x', 'dependent', 'y');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    
    slope_est = 1 / (max(x) - min(x)); 
    opts.StartPoint = [quantile(x, 0.5), slope_est];
    
    [fitresult, gof] = fit(x, y, ft, opts);
    params = coeffvalues(fitresult); 
end

function [params, gof] = run_sigmoid_4param(x, y)
    ft = fittype('a + (b-a) ./ (1+exp((c-x)*d))', 'independent', 'x', 'dependent', 'y');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    
    slope_est = 1 / (max(x) - min(x));
    opts.StartPoint = [0, 1, quantile(x, 0.5), slope_est];
    
    % Dynamic Bounds to match Baek's [-0.5, 6.5] for 1-5 data
    opts.Lower = [0, 0, min(x)-1.5, 0]; 
    opts.Upper = [1, 1, max(x)+1.5, inf];
    
    [fitresult, gof] = fit(x, y, ft, opts);
    params = coeffvalues(fitresult); 
end