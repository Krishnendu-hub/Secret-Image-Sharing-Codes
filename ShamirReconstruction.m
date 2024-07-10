function s = ShamirReconstruction(d,k)

if size(d,1)<k
    error('insufficient pieces of information parts for reconstruction')
end
% obtain k pieces of info
x = d(1:k,1)';
y = d(1:k,2)';
% generate Largrange Polynomial
% see detailed algorithm on
% http://en.wikipedia.org/wiki/Shamir's_Secret_Sharing
% and
% http://en.wikipedia.org/wiki/Lagrange_polynomial
delElement = @(x,i) [x(1:i-1),x(i+1:end)];
lj = @(a,j) prod(delElement(a-x,j)'./delElement(x(j)-x,j)');
Lj = @(y,j) y(j)*lj(0,j);
% reconstruct secret info
s = 0;
for i = 1:k
    s = s+Lj(y,i);
end
