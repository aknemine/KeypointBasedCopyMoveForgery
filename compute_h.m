function [ H ] = compute_h(im1_pts, im2_pts)
%COMPUTE_H Return the homography computed between the corresponding points.
%   The homography will transform the im1_pts to the im2_pts
%   i.e. im2_pts(1,:) ~= H * im1_pts(1,:); 
    
    n = size(im1_pts, 1);  % number of points that exist
    % Setup the matrix for solving least squares
    % Remove +3 for homog
    A = spalloc(n * 2 + 3, 8, 6 * 2 * n +1 + 3);
    b = zeros(n * 2 + 3, 1);
    eq = 1;
    for i=1:n
        % TODO: check if x, y need to be swapped
        x = im1_pts(i, 1);  y = im1_pts(i, 2);
        u = im2_pts(i, 1); v = im2_pts(i, 2);
        
        % new first equation
        A(eq, 1) = x; A(eq, 2) = y; A(eq, 3) = 1; 
        A(eq, 7) = -u * x; A(eq, 8) = -u * y;
        b(eq) = u;
        eq = eq + 1;
        
        % new second equation
        A(eq, 4) = x; A(eq, 5) = y; A(eq, 6) = 1;
        A(eq, 7) = -v * x; A(eq, 8) = -v * y;
        b(eq) = v;
        eq = eq + 1;
    end
    % Affine only
    A(eq, 7) = 1;
    b(eq) = 0;
    A(eq + 1, 8) = 1;
    b(eq + 1) = 0;
    A(eq + 2, 9) = 1;
    b(eq + 2) = 1;
    
    v = A \ b;
    H = [v(1) v(2) v(3); v(4) v(5) v(6); v(7) v(8) 1;];
end