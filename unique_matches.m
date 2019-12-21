function output = unique_matches( matches )
%UNIQUE_MATCHES Summary of this function goes here
%   Detailed explanation goes here
    A = cat(2, matches.source, matches.target);
    B = unique(A,'rows');
    output = struct('source', B(:,1:2), 'target', B(:,3:4));
end