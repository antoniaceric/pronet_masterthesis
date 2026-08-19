function v = vectorize_rdm(rdm)
% Converts a square RDM into a column vector of its
% lower triangular elements (excluding the diagonal).
    mask = tril(true(size(rdm)), -1);
    v = rdm(mask);
end % End of vectorize_rdm