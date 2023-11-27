function [y,ind] = maxN(in)
%maximum of image or matrix
N = length(size(in));

[y ind] = max(flatrow(in));
