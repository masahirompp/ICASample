function V = eigenvalue_decomposition(A)

% A=E*D*E' 
% eigenvalue decomposition
% 
% V=E*D^(-1/2)*E'
% decorrelation or orthogonalization matrix

[E,D]=eig(A);		% eigenvalue

% A=E*D*E' eigenvalue decomposition

% decorrelation matrix
V=(D^(-1/2))*E';

