function [peak,allpeak] = DA_ES_P(Fn,Run)
    %% Choose the Funciton and Run
    if nargin == 0
        Fn = 2;
        Run = 5;
    end
    %% Showing the problem's information
    fprintf("----Function %d Run %d is running----\n",Fn,Run);
    
    %% Initialize the problem's parameters
    pro = DMMOP(Fn);

    %% Setup random seed parameters
    algRand = RandStream.create('mt19937ar','seed', Run);
    RandStream.setGlobalStream(algRand);

    %% Initialize the related parameters
    D = pro.D;  % The dim of the problem
    min_popsize = 7 + floor(3*log(D)); % The minimum size of the population
    lambda = min_popsize; % The minimum size of the population
    data = []; % Store the fitness
    predict_fits = zeros(D,10);

    %% The set of the global optimal solutions
    bestmem_set = [];       
    bestval_set = [];
    temp_best_pop  = [];

    %% The set of the global optimal solutions 
    old_bestmem_set = [];

    %% Setup the network
    input_size = 40;
    output_size = 10;
    hidden_size = 300; % LSTM隐藏层大小
    layers = [ ...
    sequenceInputLayer(input_size)
    lstmLayer(hidden_size)
    fullyConnectedLayer(output_size)
    regressionLayer];
   
    %% Test Area ---- Skip the previous environments and test the specified environment
%     while ~pro.Terminate()
%         if pro.env==5
%             break
%         end
%         fprintf("%d\n",pro.env+1);
%         while ~pro.CheckChange(bestmem_set,bestval_set)
%             rest = pro.freq - rem(pro.evaluated, pro.freq);
%             useless_pop = rand(rest, D) .* (pro.upper - pro.lower) + pro.lower;
%             useless_fit = pro.GetFits(useless_pop);         
%         end
%     end 

%% Determine if the problem is terminated 
    while ~pro.Terminate()
        %% Initialization of the populations - Use IDBPI to generate populations based on density (IDBPI)
        path = sprintf('./IDBPI_POP/Init_Pop_dim%02d_Run%02d.mat',D,Run);
        load(path,"-mat",'pop_I');
        init_pop = pop_I; 

        %% Get the fitness of the populations
        [data,fits,layers,predict_fits] = GetFitness(data,pro,pop_I,layers,predict_fits);
    
        %% Add the position of the past environment's global optimal solutions to this environment
        [old_pop,old_fit,old_bestmem_set] = add_old_pop(old_bestmem_set,pro,algRand,lambda);
        init_pop = [init_pop;old_pop];
        fits = [fits;old_fit];

        %% FDBPI and evaluate the fitness
        pop_F = FDBPI(pro.lower, pro.upper, init_pop, fits, 0.1*pro.freq, D,algRand);
        fits_F= pro.GetFits(pop_F); 
        init_pop = [init_pop;pop_F];
        fits = [fits;fits_F];

        %% NBC_LP
        [fits, sort_index] = sort(fits, 'descend');
        init_pop = init_pop(sort_index, :);
        species = NBC_lp(init_pop);
        num_species = length(species);
        species_arr = [species.len];
        imp_spec_count = length(species_arr(species_arr>=min_popsize));% The number of the important species
        [~,sort_index] = sort(species_arr,'descend');

        %% Generate groups based on clusters and record the groups's postions
        groups = init_groups(pro,lambda,init_pop,species,fits,sort_index,algRand,bestmem_set);
        track  = [track_record(groups)];

        %% Initialize the archive solutions
        bestmem_set = zeros(1,D);
        bestval_set = -999;

        %% Sort the groups by fitness
        [~,index] = sort([groups.bestval],'descend');
        groups = groups(index);
        num_groups = length(groups);
        temp_pop = zeros(num_groups*lambda,pro.D);
        % Store all the population
        for j =1:num_groups
            temp_pop((j-1)*lambda+1:j*lambda,:) = groups(j).OPTS.pop;
        end
        num_groups = length(groups);
        
        %% Initialize the contribution of each groups
        itermax = min(20,ceil((0.25*pro.freq)/(num_groups*min_popsize)));
        ind = 1;
        while ind<=num_groups
            [new_groups,track] = KE_CMA_ES(pro, groups(ind), pro.lower, pro.upper, itermax, algRand,temp_pop,ind,temp_best_pop,imp_spec_count,track);
            if ~isempty(new_groups)
                groups(ind) = new_groups;
                temp_pop((ind-1)*lambda+1:ind*lambda,:) = groups(ind).OPTS.pop;
            else
                groups(ind)= [];
                ind = ind - 1;
                num_groups = num_groups-1;
            end
            ind = ind + 1;
        end
        num_groups = length(groups);
        val = [groups.bestval];  pop = cat(1, groups.bestmem);
        [bestval, ibest] = max(val);  bestmem = pop(ibest, :);
        bestval = round(bestval);
        [~,index] = sort([groups.bestval],'descend');
        groups = groups(index);% Sort the groups by the bestmem
        for j =1:num_groups
            temp_pop((j-1)*lambda+1:j*lambda,:) = groups(j).OPTS.pop;
        end
        track = [track;track_record(groups)];

        %% Iteration
        while ~pro.CheckChange(bestmem_set,bestval_set)
            %% Choose the best and second-best groups
            num_groups = length(groups);
            i = 1;
            while i<=num_groups
                if (groups(i).delta<=0.01 || groups(i).OPTS.sigma<0.05)&& abs(groups(i).bestval-bestval)>5
                    groups(i) = [];
                    i = i - 1;
                    num_groups = num_groups - 1;
                end
                i = i + 1;
            end
            % Restart if the group is empty
            if isempty(groups)
                rest = pro.freq - rem(pro.evaluated, pro.freq);
                if rest == 0 || rest == pro.freq
                    continue;
                else
                    if rest>=min_popsize
                        disp("----Restart!----");
                        groups = restart(track,pro,algRand);
                        if isempty(groups)
                            break;
                        end
    
                    else
                        useless_pop = rand(10*rest,pro.D) .* (pro.upper - pro.lower) + pro.lower;
                        rho = kernel(useless_pop,track,pro.D,0);
                        [~,index] = sort(rho);
                        get_pop = useless_pop(index(1:rest),:);
                        get_fit = pro.GetFits(get_pop);
                        for i = 1:size(get_pop,1)
                            if abs(get_fit(i)-bestval)<=1e-3
                                bestmem_set = [bestmem_set;get_pop(i,:)];
                                bestval_set = [bestval_set;get_fit(i,:)];
                            end
                        end
                    end
                    continue;
                end
            end
            num_groups = length(groups);
            fprintf("Evaluated:%d\n",rem(pro.evaluated,pro.freq));
            val = [groups.bestval]; pop = cat(1, groups.bestmem);
            [~, first_idx] = max(val);

            delta = [groups.delta];
            expected_gen = ceil((bestval - val)./(delta./itermax));
            expected_gen(first_idx) = Inf;
            if ~isempty(bestmem_set)
                gdis = pdist2(pop, bestmem_set);
                gdis = min(gdis, [], 2);
                temp_arr = [expected_gen', -gdis];
                [~, idx] = sortrows(temp_arr);
            else
                randnum = randperm(length(groups));
                temp_arr = [expected_gen', randnum'];
                [~, idx] = sortrows(temp_arr);
            end
            second_idx = idx(1);
         
            %% The best evolved populations
            ind = first_idx;
            itermax = 20;
            rest = pro.freq - rem(pro.evaluated, pro.freq);
            if pro.change == 1
                continue;
            end
            if min_popsize*itermax <= rest
                [new_groups,track] = KE_CMA_ES(pro, groups(ind), pro.lower, pro.upper, itermax, algRand,temp_pop,ind,temp_best_pop,imp_spec_count,track);
                track = [track;track_record(new_groups)];
            else
                itermax = min(itermax, floor(rest/min_popsize));
                if itermax>0
                    [new_groups,track] = KE_CMA_ES(pro, groups(ind), pro.lower, pro.upper, itermax, algRand,temp_pop,ind,temp_best_pop,imp_spec_count,track);
                    track = [track;track_record(new_groups)];
                end            
                rest = pro.freq - rem(pro.evaluated, pro.freq);
                if rest == 0 || rest == pro.freq
                    continue;
                else
                    if rest>=min_popsize
                        disp("----Restart!----");
                        groups = restart(track,pro,algRand);
                        if isempty(groups)
                            break;
                        end
                        track = [track;groups.OPTS.pop];
                    else
                        useless_pop = rand(10*rest,pro.D) .* (pro.upper - pro.lower) + pro.lower;
                        rho = kernel(useless_pop,track,pro.D,0);
                        [~,index] = sort(rho);
                        get_pop = useless_pop(index(1:rest),:);
                        get_fit = pro.GetFits(get_pop);
                        for i = 1:size(get_pop,1)
                            if abs(get_fit(i)-bestval)<=1e-3
                                bestmem_set = [bestmem_set;get_pop(i,:)];
                                bestval_set = [bestval_set;get_fit(i,:)];
                            end
                        end
                    end
                    continue;
                end
            end
            
            %% Update the groups
            idx = first_idx;
            if ~isempty(new_groups)
                groups(ind) = new_groups;
                if groups(ind).bestval > bestval
                    bestmem = groups(ind).bestmem;
                    bestval = groups(ind).bestval;
                end
            else
                groups(ind)= [];
%                 if first_idx<=second_idx
%                     second_idx = second_idx-1;
%                 end
%                 idx = second_idx;
                continue;
            end
            num_groups = length(groups);
            [~,index] = sort([groups.bestval],'descend');
            groups = groups(index);
            for j =1:num_groups
                temp_pop((j-1)*lambda+1:j*lambda,:) = groups(j).OPTS.pop;
            end
            % Restart if the group is empty
            if isempty(groups)
                disp("----Restart!----");
                groups = restart(track,pro,algRand);
                if isempty(groups)
                    break;
                end
                track = [track;groups.OPTS.pop];
                continue;
            end

            %% The second-best evolved populations
            ind = second_idx;
            itermax = 20;
            rest = pro.freq - rem(pro.evaluated, pro.freq);
            if pro.change == 1
                continue;
            end
            if min_popsize*itermax <= rest
                [new_groups,track] = KE_CMA_ES(pro, groups(ind), pro.lower, pro.upper, itermax, algRand,temp_pop,ind,temp_best_pop,imp_spec_count,track);
                track = [track;track_record(new_groups)];
            else
                itermax = min(itermax, floor(rest/min_popsize));
                if itermax>0
                    [new_groups,track] = KE_CMA_ES(pro, groups(ind), pro.lower, pro.upper, itermax, algRand,temp_pop,ind,temp_best_pop,imp_spec_count,track);
                    track = [track;track_record(new_groups)];
                end            
                rest = pro.freq - rem(pro.evaluated, pro.freq);
                if rest == 0 || rest == pro.freq
                    continue;
                else
                    if rest>=min_popsize
                        disp("----Restart!----");
                        groups = restart(track,pro,algRand);
                        if isempty(groups)
                            break;
                        end
                        track = [track;groups.OPTS.pop];
                        
                    else
                        useless_pop = rand(10*rest,pro.D) .* (pro.upper - pro.lower) + pro.lower;
                        track = [track;useless_pop];
                        rho = kernel(useless_pop,track,pro.D,0);
                        [~,index] = sort(rho);
                        get_pop = useless_pop(index(1:rest),:);
                        get_fit = pro.GetFits(get_pop);
                        for i = 1:size(get_pop,1)
                            if abs(get_fit(i)-bestval)<=1e-3
                                bestmem_set = [bestmem_set;get_pop(i,:)];
                                bestval_set = [bestval_set;get_fit(i,:)];
                            end
                        end
                    end
                    continue;
                end
            end
            
            %% Update the groups
            if ~isempty(new_groups)
                groups(ind) = new_groups;
                if groups(ind).bestval > bestval
                    bestmem = groups(ind).bestmem;
                    bestval = groups(ind).bestval;
                end
            else
                groups(ind)= [];
%                 if idx == second_idx && first_idx~=second_idx
%                     continue;
%                 end
            end
            num_groups = length(groups);
            [~,index] = sort([groups.bestval],'descend');
            groups = groups(index);
            for j =1:num_groups
                temp_pop((j-1)*lambda+1:j*lambda,:) = groups(j).OPTS.pop;
            end
 
            % Restart if the group is empty
            if isempty(groups)
               if rest <= min_popsize
                   useless_pop = rand(rest, D) .* (pro.upper - pro.lower) + pro.lower;
                   useless_fit = pro.GetFits(useless_pop);
                   continue;
               end
                disp("----Restart!----");
                groups = restart(track,pro,algRand);
                if isempty(groups)
                    break;
                end
                track = [track;groups.OPTS.pop];
                continue;
            end

            %% Save the archive solution if it converges to the global optimum solution
            if groups(idx).cc<=1e-7
                p = 1;  
                if ~isempty(bestmem_set)
                    dis_arr = pdist2(groups(idx).bestmem, bestmem_set);
                    if min(dis_arr) >= 1e-3
                        bestmem_set = [bestmem_set; groups(idx).bestmem];
                        bestval_set = [bestval_set; groups(idx).bestval];
                    end
                end
                temp_best_pop = [temp_best_pop;groups(idx).OPTS.pop];
                groups(idx) = [];
                % Restart if the group is empty
                if isempty(groups)
                    disp("----Restart!----");
                    groups = restart(track,pro,algRand);
                    if isempty(groups)
                        break;
                    end
                    track = [track;groups.OPTS.pop];
                    continue;
                end
                continue
            end

            % Restart if the group is empty
            if isempty(groups)
               if rest <= min_popsize
                   useless_pop = rand(rest, D) .* (pro.upper - pro.lower) + pro.lower;
                   useless_fit = pro.GetFits(useless_pop);
                   continue;
               end
                disp("----Restart!----");
                groups = restart(track,pro,algRand);
                if isempty(groups)
                    break;
                end
                track = [track;groups.OPTS.pop];
                continue;
            end
        end
        peak = length(bestval_set(abs(bestval_set-max(bestval_set))<1e-5));
        old_bestmem_set = [old_bestmem_set;bestmem_set];
        fprintf("---Function%d Run%d Env%d/60 Find %d个峰---\n",Fn,Run,pro.env,peak);
        temp_best_pop = [];
    end
    [peak, allpeak] = pro.GetPeak();
    PR = sum(peak, 2) ./ sum(allpeak, 2);
    fprintf("1e-3:")
    disp(peak(1,:));
    fprintf("1e-4:")
    disp(peak(2,:));
    fprintf("1e-5:")
    disp(peak(3,:));
    disp(allpeak);
    disp(PR);
    filename = sprintf('./PEAKS/F%d_Run%d.mat', Fn, Run);
    save(filename,'peak','allpeak',"PR");
end