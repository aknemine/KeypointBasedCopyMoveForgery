function output = partition_matches( matches )
%PARTITION_MATCHES Summary of this function goes here
%   Detailed explanation goes here
    A = cat(2, matches.source, matches.target);
    M = [];

    % Initialize first cluster center and initial grouping
    idx = randperm(size(A, 1), 1);

    M = [M; A(idx,:);];

    centroid = A(idx, 1:2);
    A(idx, :) = [];
    
    % While a line hasn't been 
    while size(A, 1) > 0
        % Prepare centroid
        c1 = centroid(1) * ones(size(A, 1), 1); c2 = centroid(2) * ones(size(A, 1), 1);
        C = cat(2, c1, c2, c1, c2);
        
        % Compute distance to centroid
        I = (A - C) .^ 2;
        I_dists = cat(2, I(:, 1) + I(:, 2), I(:, 3) + I(:, 4));
        min_dist = min(I_dists, [], 2);
        [idx, col] = find(min_dist==min(min_dist));
        idx = idx(1);
        [i, j] = find(I_dists(idx, :) == min(I_dists(idx, :)));
        j = j(1);
        if j == 1
            result = A(idx,:);
        else
            % invert it
            result = cat(2, A(idx,3:4), A(idx,1:2));
        end 
        
        % Update the output
        M = [M; result;];
        
        % Update the cluster center
        centroid = mean(M(:,1:2));
        
        % Remove the row that we already added
        A(idx, :) = [];
    end
    
    output = struct('source', M(:,1:2), 'target', M(:,3:4));
end
