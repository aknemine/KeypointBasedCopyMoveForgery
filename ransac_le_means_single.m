function output = ransac_le_means_single(matches, iters, error_threshold )
%RANSAC_LE_MEANS Summary of this function goes here
%   Detailed explanation goes here
    error_threshold = error_threshold ^ 2;
    O = ones(size(matches.source, 1), 1);
    best_inlier_indices = [];
    for i=1:iters
        indices = randperm(size(matches.source, 1), 6);
        
        % Get the subset of matches
        source = matches.source(indices, :);
        target = matches.target(indices, :);
        tmp_matches = partition_matches(struct('source', source, 'target', target));
        H = compute_h(tmp_matches.source, tmp_matches.target);
        
        % TODO: Should we re-aggregate w.r.t. to this centroid?
        X = matches.source(:, 1);  Y = matches.source(:, 2);
        query = cat(1, X', Y', O');
        
        XP = matches.target(:, 1);  YP = matches.target(:, 2);
        expected = cat(1, XP', YP', O');
        
        results = H * query;
        for j=1:size(matches.source, 1)
            f = results(3, j);
            if f ~= 1.0
                results(:, j) = results(:, j) / f;
            end
        end
        
        distances = sum((results - expected) .^ 2, 1);
        inlier_indices = find(distances < error_threshold);
         
        if numel(inlier_indices) > numel(best_inlier_indices) && numel(inlier_indices) >= 10
            best_inlier_indices = inlier_indices;
        end
         
    end
    
    source = matches.source(best_inlier_indices, :);
    target = matches.target(best_inlier_indices, :);
    
    output = struct('source', source, 'target', target, 'indices', best_inlier_indices);
%     tmp_matches = partition_matches(struct('source', source, 'target', target));
%     H = compute_h(tmp_matches.source, tmp_matches.target);
% 
%     X = matches.source(:, 1);  Y = matches.source(:, 2);
%     query = cat(1, X', Y', O');
% 
%     XP = matches.target(:, 1);  YP = matches.target(:, 2);
%     expected = cat(1, XP', YP', O');
%     
%     results = H * query;
%     for j=1:size(matches.source, 1)
%         f = results(3, j);
%         if f ~= 1.0
%             results(:, j) = results(:, j) / f;
%         end
%     end
% 
%     distances = sum((results - expected) .^ 2, 1);
%     inlier_indices = find(distances < error_threshold);
% 
%     output = struct('source', matches.source(inlier_indices,:), 'target', matches.target(inlier_indices,:));
end