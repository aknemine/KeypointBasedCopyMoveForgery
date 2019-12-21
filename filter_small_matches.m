function [ output ] = filter_small_matches(matches, threshold)
%FILTER_SMALL_POINTS Summary of this function goes here
%   Filters the matches by removing matches that have too small of a
%   distance.
    D = sqrt(sum((matches.source - matches.target) .^ 2, 2));
    idxs = find(D > threshold);
    output = struct('source', matches.source(idxs,:), 'target', matches.target(idxs,:));
end