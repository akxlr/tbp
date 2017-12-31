function total_var = decomp_var(decomp)
    mat_sq = decomp.lambda(1) * (decomp.U{1}(:,1) .^ 2) * (decomp.U{2}(:,1)' .^ 2) + ...
             decomp.lambda(2) * (decomp.U{1}(:,2) .^ 2) * (decomp.U{2}(:,2)' .^ 2);
    
    vars = mat_sq - decomp2mat(decomp) .^ 2;
    total_var = sum(vars(:));
end