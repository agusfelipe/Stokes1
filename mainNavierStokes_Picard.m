% This program solves the Navier-Stokes cavity problem
clear; close all; clc

addpath('Func_ReferenceElement')

dom = [0,1,0,1]; 
Re = 100; 
nu = 1/Re; 

% Element type and interpolation degree
% (0: quadrilaterals, 1: triangles, 11: triangles with bubble function)
elemV = 0; degreeV = 2; degreeP = 1;
% elemV = 1; degreeV = 2; degreeP = 1;
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
K = nu*K; 
[ndofP,ndofV] = size(G); 

% Prescribed velocity degrees of freedom
[dofDir,valDir,dofUnk,confined] = BC_red(X,dom,ndofV);
nunkV = length(dofUnk); 
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
A = [Kred   Gred'
     Gred   zeros(nunkP)]; 

% Initial guess
disp(' ')
IniVelo_file = input('.mat file with the initial velocity = ','s'); 
if isempty(IniVelo_file)
    velo = zeros(ndofV/2,2); 
    y2 = dom(4); 
    nodesY2 = find(abs(X(:,2)-y2) < 1e-6); 
    velo(nodesY2,1) = 1; 
else
    load(IniVelo_file); 
end
pres = zeros(nunkP,1);
veloVect = reshape(velo',ndofV,1);
sol0  = [veloVect(dofUnk);pres(1:nunkP)];

iter = 0; tol = 0.5e-08; 
while iter < 100
    fprintf('Iteration = %d\n',iter);
    
    C = ConvectionMatrix(X,T,referenceElement,velo);
    Cred = C(dofUnk,dofUnk); 
    
    Atot = A;
    Atot(1:nunkV,1:nunkV) = A(1:nunkV,1:nunkV) + Cred; 
    btot = [fred - C(dofUnk,dofDir)*valDir; zeros(nunkP,1)]; 
    
    % Computation of residual
    res = btot - Atot*sol0;
    % Computation of velocity and pressure increment
    solInc = Atot\res;
    
    % Update the solution
    veloInc = zeros(ndofV,1); 
    veloInc(dofUnk) = solInc(1:nunkV); 
    presInc = solInc(nunkV+1:end); 
    velo = velo + reshape(veloInc,2,[])'; 
    pres = pres + presInc; 
    
    % Check convergence
    delta1 = max(abs(veloInc)); 
    delta2 = max(abs(res)); 
    fprintf('Velocity increment=%8.6e, Residue max=%8.6e\n',delta1,delta2); 
    if delta1 < tol*max(max(abs(velo))) && delta2 < tol
        fprintf('\nConvergence achieved in iteration number %g\n',iter); 
        break
    end
    
    % Update variables for next iteration
    veloVect = reshape(velo',ndofV,1);
    sol0 = [veloVect(dofUnk); pres]; 
    iter = iter + 1; 
end

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
