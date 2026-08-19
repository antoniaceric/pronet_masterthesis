function indicator_matrix = matlab_pairwise_contrast(index_vector)
% MATLAB translation of the Python function pairwise_contrast from rsa
% toolbox
% This creates a contrast matrix for all unique pairs of conditions
% found in the index_vector.

    % Find unique conditions
    c_unique = unique(index_vector);
    n_unique = numel(c_unique);
    
    % Get dimensions
    rows = numel(index_vector);
    cols = (n_unique * (n_unique - 1)) / 2; % nchoosek(n_unique, 2)
    
    indicator_matrix = zeros(cols, rows);
    n_row = 1; % MATLAB is 1-indexed
    
    % Build the contrast matrix, one row for each pair
    for i = 1:n_unique
        for j = (i + 1):n_unique
            
            % Get logical indices for condition i
            select_i = (index_vector == c_unique(i));
            % Set positive weights, normalized by number of repeats
            indicator_matrix(n_row, select_i) = 1 / sum(select_i);
            
            % Get logical indices for condition j
            select_j = (index_vector == c_unique(j));
            % Set negative weights, normalized by number of repeats
            indicator_matrix(n_row, select_j) = -1 / sum(select_j);
            
            n_row = n_row + 1;
        end
    end
end % End of matlab_pairwise_contrast
