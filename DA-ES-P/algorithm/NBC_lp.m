function [species, meandis]= NBC_lp(pop)
    [n, dim] = size(pop);       % edge
    min_num_edge = 10;
    
    factor = 4 - log(dim);
    
    nbc=zeros(n,3);             % Info Matrix,Start,End,Distance
    nbc(1:n, 1) = 1:n;          % Start
    nbc(1, 2) = -1;             % End
    nbc(1, 3) = 0;              % Distance
    
    for i = 2 : n
        arrdis = pdist2(pop(i, :), pop(1:i-1, :));
        [u, v] = min(arrdis);
        nbc(i,2) = v;
        nbc(i,3) = u;
    end

    % filename = sprintf('./nbc/nbc_F%d_runs%d.mat', func, runs);
    % load(filename);
    
    % nbc = nbc(1:n, :);
    
    % meandis=factor*mean(nbc(2:n,3));
    % nbc(nbc(:,3)>meandis,2)=-1;
    % nbc(nbc(:,3)>meandis,3)=0;

    meandis=factor*mean(nbc(2:n,3));
    if  length(nbc(nbc(:,3)>meandis)) >= min_num_edge
        nbc(nbc(:,3)>meandis,2)=-1;
        nbc(nbc(:,3)>meandis,3)=0;
    else
        pred = sortrows(nbc,3,'descend');
        meandis = pred(min_num_edge, 3);
        nbc(nbc(:,3)>=meandis,2)=-1;
        nbc(nbc(:,3)>=meandis,3)=0;
    end
    
    
    seeds=nbc(nbc(:,2)==-1,1);
    m=zeros(n,2); % Save the index of seed and speices
    m(1:n,1)=1:n;
    for i=1:n
        j=nbc(i,2);
        k=j;
        while j~=-1
            k=j;
            j=nbc(j,2);
        end
        if k==-1
            m(i,2)=i;
        else
            m(i,2)=k;
        end
    end  
    
    % Construct the result
    species = struct();
    for i=1:length(seeds)
        species(i).seed = seeds(i);
        species(i).idx = m(m(:, 2) == seeds(i), 1);
        species(i).len = length(species(i).idx);
    end

end

