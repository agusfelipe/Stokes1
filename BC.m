function [A, b, nDir, confined] = BC(X,dom,n)
% [A, b, nDir, confined] = BC(X,dom,n)
% Matrices to impose Dirichlet boundary conditions using Lagrange
% multipliers on a rectangular domain
% Input: 
%    X: nodal coordinates
%    dom: domain description [x1,x2,y1,y2]
%    n: number of velocity degrees of freedom
% Output:
%    A,b: matrix and r.h.s. vector to impose the boundary conditions using
%         Lagrange multipliers
%    nDir: number of prescribed degrees of freedom
%    confined: 


tol = 1e-6; 
y1 = dom(3); nodesY1 = find(abs(X(:,2)-y1) < tol); 
y2 = dom(4); nodesY2 = find(abs(X(:,2)-y2) < tol); nNodesY2 = length(nodesY2); 
x1 = dom(1); nodesX1 = find(abs(X(:,1)-x1) < tol & abs(X(:,2)-y1)>tol & abs(X(:,2)-y2) > tol);
x2 = dom(2); nodesX2 = find(abs(X(:,1)-x2) < tol & abs(X(:,2)-y1)>tol & abs(X(:,2)-y2) > tol);

nodesDirBC = [nodesY2; nodesY1; nodesX1; nodesX2]; 
% confined flow (velocity is imposed for all the nodes in the boundary)
confined = 1; 
% number of prescribed degrees of freedom
nDir = 2*length(nodesDirBC); 

% Degrees of freedom where the velocity is prescribed
% (horizontal and vertical component)
C = [2*nodesDirBC - 1; 2*nodesDirBC];
A = zeros(nDir,n); 
A(:,C) = eye(nDir); 
% Imposed value
b = [ones(nNodesY2,1); zeros(nDir-nNodesY2,1)]; 



