function [dofDir,valDir,dofUnk,confined] = BC_red(X,dom,ndofV)
% [dofDir,valDir,dofUnk,confined] = BC_red(X,dom,ndofV)

tol = 1e-6; 
y1 = dom(3); nodesY1 = find(abs(X(:,2)-y1) < tol); 
y2 = dom(4); nodesY2 = find(abs(X(:,2)-y2) < tol); 
x1 = dom(1); nodesX1 = find(abs(X(:,1)-x1) < tol & abs(X(:,2)-y1)>tol & abs(X(:,2)-y2) > tol);
x2 = dom(2); nodesX2 = find(abs(X(:,1)-x2) < tol & abs(X(:,2)-y1)>tol & abs(X(:,2)-y2) > tol);

confined = 1; 
dofDir = [
    2*nodesX1-1; 2*nodesX1
    2*nodesX2-1; 2*nodesX2
    2*nodesY1-1; 2*nodesY1
    2*nodesY2-1; 2*nodesY2
     ]; 
valDir = [
    zeros(size(nodesX1)); zeros(size(nodesX1))
    zeros(size(nodesX2)); zeros(size(nodesX2))
    zeros(size(nodesY1)); zeros(size(nodesY1))
    ones(size(nodesY2));   zeros(size(nodesY2))
    ];
dofUnk = setdiff(1:ndofV,dofDir); 



