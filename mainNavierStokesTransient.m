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

% figure; PlotMesh(T,X,elemV,'b-');
% figure; PlotMesh(TP,XP,elemP,'r-');

L = zeros(size(XP,1),size(XP,1));
f_q=zeros(size(XP,1),1);

% Matrices arising from the discretization
if degreeV == 2 % add stability only when needed (degreeV=1)
    [M,K,G,f] = StokesSystem(X,T,XP,TP,referenceElement);
else
    [M,K,G,f,L,f_q] = StokesSystemStable(X,T,XP,TP,referenceElement,nu);
end
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
   L(1,:) = [];
   L(:,1) = [];
   f_q(1) = [];
else
   nunkP = ndofP;
end

%f = f - K(:,dofDir)*valDir;
Mred = M(dofUnk,dofUnk); 
Kred = K(dofUnk,dofUnk); 
Gred = G(:,dofUnk); 
fred = f(dofUnk); 


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

iter = 0; 
tol = 1e-12;
dt = 0.01;
teta = 1;
step = 0;
nstep = 100;
velo_aux=velo;
veloVect_aux=veloVect;
pres_aux=pres;
tic
method = cinput('Select iterative method ([0] - Picard, [1] - Newton-Raphson):',0);
if method ==0
    while step < nstep
        step = step +1;
        while iter < 10
            fprintf('Iteration = %d\n',iter);
            iter = iter + 1;
            C = ConvectionMatrix(X,T,referenceElement,velo_aux);
            Cred = C(dofUnk,dofUnk); 
            fredn = fred - (K(dofUnk,dofDir)+C(dofUnk,dofDir))*valDir;

            Atot = [Mred+teta*dt*(Kred+Cred)   dt*teta*Gred'
             Gred   L]; 
            btot = [dt*(fredn-(Kred+Cred)*veloVect_aux(dofUnk))-Gred'*pres_aux; f_q]; 

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
        %     delta1 = max(abs(veloInc)); 
        %     delta2 = max(abs(res));
        %     Conv(iter+1,1) = iter+1;
        %     Conv(iter+1,2) = max(abs(res));
        %     fprintf('Velocity increment=%8.6e, Residue max=%8.6e\n',delta1,delta2); 
        %     if delta1 < tol*max(max(abs(velo))) && delta2 < tol
        %         fprintf('\nConvergence achieved in iteration number %g\n',iter); 
        %         break
        %     end

            % Update variables for next iteration
            veloVect = reshape(velo',ndofV,1);
            sol0 = [veloVect(dofUnk); pres]; 

        end
        veloInc = zeros(ndofV,1); 
        veloInc(dofUnk) = solInc(1:nunkV); 
        presInc = solInc(nunkV+1:end); 
        velo = velo + reshape(veloInc,2,[])'; 
        pres = pres + presInc; 
    end
else %Newton-Raphson
while iter < 100
    fprintf('Iteration = %d\n',iter);
    
    [C1,C2] = ConvectionMatrixNR(X,T,referenceElement,velo);
    Cred = C1(dofUnk,dofUnk); 
    Cred2= C2(dofUnk,dofUnk);
    
    Atot = A;
    Atot(1:nunkV,1:nunkV) = A(1:nunkV,1:nunkV) + Cred; 
    btot = [fred - C1(dofUnk,dofDir)*valDir; f_q]; 
    
    J = A;
    J(1:nunkV,1:nunkV) = A(1:nunkV,1:nunkV) + Cred + Cred2;
%     Z = ConvectionMatrix2(X,T,referenceElement,velo);
%     Zred = Z(dofUnk,dofUnk); 
%     J(1:nunkV,1:nunkV) = A(1:nunkV,1:nunkV) + Cred + Zred;
     
    res = Atot*sol0 - btot;
   
    % Computation of velocity and pressure increment
    solInc = -J\res;
    
    % Update the solution
    veloInc = zeros(ndofV,1); 
    veloInc(dofUnk) = solInc(1:nunkV); 
    presInc = solInc(nunkV+1:end); 
    velo = velo + reshape(veloInc,2,[])'; 
    pres = pres + presInc; 
    
    % Check convergence
    delta1 = max(abs(veloInc)); 
    delta2 = max(abs(res));
    Conv(iter+1,1) = iter+1;
    Conv(iter+1,2) = max(abs(res));
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

end
toc
hold on
plot(Conv(:,1)-1,log10(Conv(:,2)),'k:o','LineWidth',2,'MarkerSize',6); %b-s k:o
l = legend('Newton-Raphson','Picard''s');
ylabel('log_{10}(Maximum Residual)')
xlabel('No. of Iterations')

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
