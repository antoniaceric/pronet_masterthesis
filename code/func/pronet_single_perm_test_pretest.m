function [R_corrected, R_real, R_perm] = pronet_single_perm_test_pretest(neural_RDM, model_RDM_Categorical)
% Performs a subject-level permutation test for a SINGLE model
%
% OUTPUTS:
%   R_corrected - Struct with median-corrected R-squared values
%   R_real      - Struct with original, uncorrected R-squared values
%   R_perm      - Struct with null distributions (raw and zeroed)

% --- 1. Get the "real" R-squared values ---
R_gen_real = pronet_rsa_modeling_pretest(neural_RDM, model_RDM_Categorical);

R_real = struct('R_general_raw', R_gen_real);

% --- 2. Create the Null Distribution ---
n_conds = 5;
all_perms = perms(1:n_conds);

% Remove the identity permutation [1 2 3 4 5]
identity_perm = 1:n_conds;
[~, loc] = ismember(identity_perm, all_perms, 'rows');
all_perms(loc, :) = [];

null_perms = all_perms; 
n_perms = size(null_perms, 1); % = 119

null_R_gen = zeros(n_perms, 1);

for p = 1:n_perms
    perm_order = null_perms(p, :);
    shuffled_RDM = neural_RDM(perm_order, perm_order);
    
    % Just Categorical, output still "gen"
    null_R_gen(p) = pronet_rsa_modeling_pretest(shuffled_RDM, model_RDM_Categorical);
end

% --- 3. Calculate Median-Corrected R-squareds ---
median_null_gen = median(null_R_gen);
raw_null_R_gen = null_R_gen;

% zeroing out negative rsq values
zero_median_null_gen = max(0, median_null_gen);
null_R_gen = max(0, null_R_gen);
R_gen_real = max(0, R_gen_real);

% r-values via squareroot
null_r_gen = sqrt(null_R_gen);
median_null_r_gen = sqrt(zero_median_null_gen);
r_gen_real = sqrt(R_gen_real);

% --- Fisher z-transformation ---
z_median_null_gen = atanh(median_null_r_gen);
z_gen_real = atanh(r_gen_real);

% Subtract the median of the null from the real Z-value
Z_gen_corr = z_gen_real - z_median_null_gen;

R_corrected = struct(...
    'R_general', Z_gen_corr, ...
    'Baseline_Z_general', z_median_null_gen);

 R_perm = struct(...
    'raw_R_general', raw_null_R_gen, ...
    'zeroed_R_general', null_R_gen);

R_real.R_general_zeroed = R_gen_real;

end