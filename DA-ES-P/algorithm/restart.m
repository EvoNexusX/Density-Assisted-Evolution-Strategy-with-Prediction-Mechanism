function groups = restart(track,pro, algRand)
    rest = pro.freq - rem(pro.evaluated, pro.freq);
    useless_pop = rand(10*rest,pro.D) .* (pro.upper - pro.lower) + pro.lower;
    rho = kernel(useless_pop,track,pro.D,0);
    [~,index] = sort(rho);
    get_pop = useless_pop(index(1),:);
    D= pro.D;
    min_popsize = 7 + floor(3*log(D));
    groups = struct();
    i = 1;
    groups(i).idx = i;
    groups(i).OPTS.first = 1;
%     groups(i).OPTS.pop = useless_pop(species(j).idx(1 : min_popsize), :);
%     groups(i).OPTS.val = pro.GetFits(useless_pop(species(j).idx(1 : min_popsize),:));
    groups(i).xmean = get_pop';
    groups(i).OPTS.pop = groups(i).xmean' + 0.5*randn(algRand, min_popsize, D);
    groups(i).OPTS.val = pro.GetFits(groups(i).OPTS.pop);
    groups(i).OPTS.count = 0;
    groups(i).OPTS.sigma = 0.5;
    groups(i).cc = std(groups(i).OPTS.val);
    [~,index] = sort(groups(i).OPTS.pop,"descend");
    groups(i).bestval = groups(i).OPTS.val(index(1));
    groups(i).bestmem = groups(i).OPTS.pop(index(1),:);
    groups(i).delta = 100;
    groups(i).iters = 0;
    num_groups = length(groups);
    groups(i).mean_distance = min(pdist2(groups(i).xmean',track));
% 
    %      
%     get_pop = useless_pop(index(1:rest),:);
%     i = 1;
%     D= pro.D;
%     min_popsize = 7 + floor(3*log(D));
%      
%     species = NBC_lp(get_pop);
%     num_species = length(species);
%     groups = struct();
%     
%     p = 0;
%     for j = 1 : num_species
%         if species(j).len < min_popsize
%             continue;
%         end
%         p = 1;
%         groups(i).idx = i;
%         groups(i).OPTS.first = 1;
%         groups(i).OPTS.pop = useless_pop(species(j).idx(1 : min_popsize), :);
%         groups(i).OPTS.val = pro.GetFits(useless_pop(species(j).idx(1 : min_popsize),:));
%         groups(i).xmean = mean(groups(i).OPTS.pop)';
%         x = groups(i).OPTS.pop - groups(i).xmean';
%         groups(i).OPTS.sigma = sqrt((1/(min_popsize*D))*sum(x(:).^2));
%         groups(i).cc = std(groups(i).OPTS.val);
%         groups(i).bestval = pro.GetFits(useless_pop(species(j).seed,:));
%         groups(i).bestmem = useless_pop(species(j).seed, :);
%         groups(i).delta = 0;
%         groups(i).iters = 0;
% 
%         i = i + 1;
%     end
%  
%     if p == 0
%         groups = [];
%         fits = pro.GetFits(useless_pop);
%     end
%  
 
end