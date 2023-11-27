function out = minN(in)
%minimum nonzero values of image or matrix
N = length(size(in));

if (N == 2),
  out = min(min(in));
elseif (N == 3),
  out = min(min(min(in)));
elseif (N == 4),
  out = min(min(min(min(in))));
end;
