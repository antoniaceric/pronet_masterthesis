function R_Categorical = pronet_rsa_modeling_pretest(neural_RDM, model_RDM_Categorical)
% Performs hierarchical linear modeling for RSA including pre-whitening
%
% This function calculates the whitening matrix 'W' internally.
% It relies on the helper functions:
%   - matlab_pairwise_contrast.m
%   - vectorize_rdm.m
%
% INPUTS:
%   neural_RDM            - 5x5 neural RDM
%   model_RDM_Acoustic    - 5x5 acoustic model RDM
%   model_RDM_Categorical - 5x5 categorical model RDM
%
% OUTPUTS:
%   R_general             - General explained variance (Full vs. Null)
%   R_unique_Acoustic     - Unique acoustic variance (Full vs. Reduced Cat)
%   R_unique_Categorical  - Unique categorical variance (Full vs. Reduced Ac)

% % --- 1. Handle Negative Neural Distances ---
% % Set all negative (Crossnobis) distances to 0,
% % as NNLS cannot fit negative target data.
% neural_RDM(neural_RDM < 0) = 0;

% --- 2. Vectorize RDMs (lower triangular) ---
% Call the helper function from your separate .m file
vec_d = vectorize_rdm(neural_RDM);
vec_M_C = vectorize_rdm(model_RDM_Categorical);

% --- 3. Ensure 'double' Data Type ---
% lsqnonneg and chol require double precision inputs.
vec_d = double(vec_d);
vec_M_C = double(vec_M_C);

% --- 4. Calculate Model-Based Whitening Matrix 'W' ---
n_conds = 5;
index_vector = (1:n_conds)'; 

c_mat = matlab_pairwise_contrast(index_vector);
sigma_k = eye(n_conds); 
v = c_mat * sigma_k * c_mat';
V = v .* v;

% Calculate the whitening matrix W = V^(-1/2)
[K, L_matrix] = eig(V);      % Get eigenvectors (K) and eigenvalues (L_matrix)
L = diag(L_matrix);         % Extract eigenvalues as a vector
L_sqrt = sqrt(abs(L));      % Take sqrt of absolute values
inv_l = 1 ./ L_sqrt;        % Calculate inverse
inv_l(~isfinite(inv_l)) = 0;% Set any non-finite (e.g., 1/0) to 0
W = K * diag(inv_l) * K';   % Reconstruct the whitening matrix W = V^(-1/2)

% --- 5. Apply whitening to data and models ---
d_w = W * vec_d;
M_C_w = W * vec_M_C;
M_intercept_w = W * ones(size(d_w)); % Intercept must also be whitened

% --- 6. Create design matrices for the models ---
M_full = [M_C_w, M_intercept_w];
M_null = M_intercept_w;

% --- 7. Fit models using NNLS ---
beta_full = lsqnonneg(M_full, d_w);
beta_null = lsqnonneg(M_null, d_w);

% --- Get d_hat (X(neuralRDM) x beta ---
d_hat_full = M_full * beta_full;
d_hat_null = M_null * beta_null;

% --- Get errors (d - d_hat)---
e_full = d_w - d_hat_full;
e_null = d_w - d_hat_null;

% --- Get RSS ---
RSS_full = sum(e_full.^2);
RSS_null = sum(e_null.^2);

% --- Get TSS ---
TSS = sum(d_w.^2);

% --- 8. Calculate R-squared (Explained Variance) ---
R_sq_full = 1 - (RSS_full / TSS);
R_sq_null = 1 - (RSS_null / TSS);

% % --- 7. Fit Models using fitlm ---
% %'Intercept', false because included M_intercept_w manually
% % Full Model
% mdl_full = fitlm(M_full, d_w, 'Intercept', false);
% R_sq_full = mdl_full.Rsquared.Ordinary;
% 
% % Reduced Model A (Acoustic Only)
% mdl_red_A = fitlm(M_red_A, d_w, 'Intercept', false);
% R_sq_red_A = mdl_red_A.Rsquared.Ordinary;
% 
% % Reduced Model C (Categorical Only)
% mdl_red_C = fitlm(M_red_C, d_w, 'Intercept', false);
% R_sq_red_C = mdl_red_C.Rsquared.Ordinary;
% 
% % Null Model (Intercept Only)
% R_sq_null = 0;

% % --- 9. Perform Variance Partitioning ---
% R_general = max(R_sq_full - max(R_sq_null, 0), 0);
% R_unique_Acoustic = max(R_sq_full - max (R_sq_red_C, 0), 0);
% R_unique_Categorical = max(R_sq_full - max(R_sq_red_A, 0),0);

% --- 9-2. Perform Variance Partitioning ---
R_Categorical = R_sq_full - R_sq_null;

end