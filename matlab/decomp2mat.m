function mat = decomp2mat(decomp)
    mat = decomp.lambda(1) * decomp.U{1}(:,1) * decomp.U{2}(:,1)' + ...
          decomp.lambda(2) * decomp.U{1}(:,2) * decomp.U{2}(:,2)';
end