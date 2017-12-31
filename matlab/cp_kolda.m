function [w, T] = cp_kolda(X, dims, K)

    % X is a 1D array in fortran (matlab, not numpy) order
    % We pass dims to Tensor as well as reshape, since reshape ignores
    % trailing dimensions of size 1.
    if (length(dims) >= 2)
        T = tensor(reshape(X, dims), dims);
    else
        T = tensor(X, dims);
    end
    decomp = cp_nmu(T, K);
    % Make weights sum to 1
    weight_sum = sum(decomp.lambda);
    decomp.lambda = decomp.lambda / weight_sum;
    decomp.U{1} = decomp.U{1} * weight_sum;
    w = decomp.lambda;
    T = decomp.U;
end

