% This program solves the cavity flow problem
clear; close all; clc

addpath('Func_ReferenceElement')

dom = [0,1,0,1]; 
mu = 1; 

% Element type and interpolation degree
% (0: quadrilaterals, 1: triangles, 11: triangles with bubble function)
%elemV = 0; degreeV = 2; degreeP = 1;
 elemV = 1; degreeV = 1; degreeP = 1;
% elemV = 11; degreeV = 1;  degreeP = 1; 
if elemV == 11
    elemP = 1; 
else
    elemP = elemV; 
end
referenceElement = SetReferenceElementStokes(elemV,degreeV,elemP,degreeP); 

nx = cinput('Number of elements in each direction',10);
ny = nx; 
[X,T,XP,TP] = CreateMeshes(dom,nx,ny,referenceElement);

figure; PlotMesh(T,X,elemV,'b-');
figure; PlotMesh(TP,XP,elemP,'r-');

% Matrices arising from the discretization
[K,G,f] = StokesSystem(X,T,XP,TP,referenceElement);
K = mu*K; 
[ndofP,ndofV] = size(G); 

% Prescribed velocity degrees of freedom
[dofDir,valDir,dofUnk,confined] = BC_red(X,dom,ndofV);
nunkV = length(dofUnk); 

% Total system of equations
if confined
   nunkP = ndofP-1;
   disp(' ')
   disp('Confined flow. Pressure on lower left corner is set to zero');
   G(1,:) = [];
else
   nunkP = ndofP;
end

f = f - K(:,dofDir)*valDir; 
Kred = K(dofUnk,dofUnk); 
Gred = G(:,dofUnk); 
fred = f(dofUnk); 

A = [Kred   Gred'; 
     Gred   zeros(nunkP)]; 
b = [fred; zeros(nunkP,1)];

sol = A\b; 

velo = zeros(ndofV,1);
velo(dofDir) = valDir; 
velo(dofUnk) = sol(1:nunkV); 
velo = reshape(velo,2,[])'; 
pres = sol(nunkV+1:end); 
if confined
    pres = [0; pres]; 
end

nPt = size(X,1); 
figure; 
quiver(X(1:nPt,1),X(1:nPt,2),velo(1:nPt,1),velo(1:nPt,2));
hold on 
plot(dom([1,2,2,1,1]),dom([3,3,4,4,3]),'k')
axis equal; axis tight

PlotStreamlines(X,velo,dom); 

if degreeP == 0
    PlotResults(X,T,pres,referenceElement.elemP,referenceElement.degreeP)
else
    PlotResults(XP,TP,pres,referenceElement.elemP,referenceElement.degreeP)
end
