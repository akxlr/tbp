function decomps = get_decomps(infile, K, n_decomps, outfile)

    T = tensor(dlmread(infile, ','));
    decomps = cell(n_decomps, 1);
    for i = 1:n_decomps
        decomps{i} = cp_nmu(T, K);
        % Make weights sum to 1
        weight_sum = sum(decomps{i}.lambda);
        decomps{i}.lambda = decomps{i}.lambda / weight_sum;
        decomps{i}.U{1} = decomps{i}.U{1} * weight_sum;
    end
    
end

