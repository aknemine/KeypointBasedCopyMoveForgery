function output = g2nn(descriptors1, descriptors2, nn_threshold)
%G2NN Summary of this function goes here
%   Detailed explanation goes here
    % For matching descriptors of 1 to descriptors of 2
    distances = dist2(descriptors1.descriptors, descriptors2.descriptors);
   
    % THINK ABOUT -- DO WE NEED A MAX G?
    
    
    % don't allow matches on same pixels
    for i=1:size(distances, 1)
        distances(i, i) = Inf;
    end
    
    % Proceed with matching
    corresponding1 = []; corresponding2 = [];
    c = 1;

    for i=1:size(distances, 1)
        ds = distances(i, :);
        
        % Initialize first NN
        nn1 = min(ds);
        matching_indices = find(ds == nn1);
        j = matching_indices(1);

        % Delete NN1 from neighbors
        ds(matching_indices) = [];

        % Find new NN2
        nn2 = min(ds);
        
        % TODO: WHAT ABOUT EXACT MATCHES?
        while ((nn1 / nn2) < nn_threshold) && numel(matching_indices) == 1 && i < j
            % Add correspondence points
            corresponding1(c, :) = descriptors1.points(i, :);
            corresponding2(c, :) = descriptors2.points(j, :);
            c = c + 1;
            
            % Set NN1 = NN2
            nn1 = nn2;
            matching_indices = find(ds == nn1);
            j = matching_indices(1);
            
            % Delete NN1 from neighbors
            ds(matching_indices) = [];

            % Find new NN2
            nn2 = min(ds);
        end
    end
    
    % we need to return correspondences
    output = struct('source', corresponding1, 'target', corresponding2);
end
