function out = flatrow(in)
%flatten into a row vector

out = reshape(in,prod(size(in)),1);
