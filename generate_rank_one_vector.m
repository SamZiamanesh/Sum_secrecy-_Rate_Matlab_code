function calculate_rank_one = generate_rank_one_vector(V)

    [eigvecs, eigvals] = eig(V);
    
    [~, max_idx] = max(diag(eigvals));

    max_eigvec = eigvecs(:, max_idx);
    
    max_eigvec_row = max_eigvec.';
    
    weight = unifrnd(-1,1);
    
    calculate_rank_one = weight * max_eigvec_row;


end