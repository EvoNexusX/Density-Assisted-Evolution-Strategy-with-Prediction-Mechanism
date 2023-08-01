function [old_pop,old_fit,bestmem_set] = add_old_pop(bestmem_set,pro,algRand,lambda)
    ss = get_ss(pro);
    old_pop = [];
    old_fit = [];
    
    bestmem_set = bestmem_set(2:end,:);
    if ~isempty(bestmem_set)
        if pro.D == 5
            for i = 1:size(bestmem_set,1)
                temp = bestmem_set(i,:)+ss*randn(algRand,lambda,pro.D);
                old_pop = [old_pop;temp;bestmem_set(i,:)];
            end
            old_fit = pro.GetFits(old_pop);
        end
    end

end