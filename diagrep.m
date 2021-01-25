function out = diagrep(A,n)
% 'diagrep' is user-made function, reducing memory from creating additional
% cells or matrix from repmat and blkdiag
%
% out = diagrep(A,n) where out is the output of the matrix with input of
% matrix A, repeated diagonally n times
ACell = repmat({A}, 1, n);
out = blkdiag(ACell{:});
end