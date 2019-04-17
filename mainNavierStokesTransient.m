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

nstep = cinput('Number of Time Steps',100);

pres = zeros(nunkP,1);
veloVect = reshape(velo',ndofV,1);
sol0  = [veloVect(dofUnk);pres(1:nunkP)];
solf = zeros(size(sol0,1),nstep);

iter = 0; 
tol = 1e-12;
dt = 0.01;

step = 0;

tic
method = cinput('Select method ([0] - Semi-Implicit, [1] - Chorin-Temam):',0);
if method == 0
    teta = cinput('Select Theta for method [1] - implicit, [0.5] - Crank Nicolson):',0.5);
end
while step < nstep
    step = step +1;
    fprintf('Time Step = %d\n',step);
    C = ConvectionMatrix(X,T,referenceElement,velo);
    Cred = C(dofUnk,dofUnk); 
    fredn = fred - (K(dofUnk,dofDir)+C(dofUnk,dofDir))*valDir;
    
    if method == 0    
        Atot = [Mred+teta*dt*(Kred+Cred)   dt*teta*Gred'
         Gred   L]; 
        btot = [dt*(fredn-(Kred+Cred)*veloVect(dofUnk)-Gred'*pres); f_q]; 
        % Computation of velocity and pressure increment
        solInc = Atot\btot;
        
        % Update the solution
        veloInc = zeros(ndofV,1); 
        veloInc(dofUnk) = solInc(1:nunkV); 
        presInc = solInc(nunkV+1:end); 
        velo = velo + reshape(veloInc,2,[])'; 
        pres = pres + presInc; 
    else
        % FIRST STEP             
        btot = dt*fredn+Mred*veloVect(dofUnk);
        Atot = Mred+dt*(Cred+Kred);
        zumba = Atot\btot;
        
        % SECOND STEP
        %if Stokes == 1
            btot = [Mred*zumba; f_q];
            Atot = [Mred Gred'*dt; Gred L];
            aux = Atot\btot;
%             clear btot;
            
        veloInc = zeros(ndofV,1); 
        veloInc(dofUnk) = aux(1:nunkV); 
        presInc = aux(nunkV+1:end); 
        velo = reshape(veloInc,2,[])'; 
        pres = presInc; 
        %else
%             btot = [G*aux1(1:nunk); 0]/dt;
%             aux = U2nd\(L2nd\btot);
%             pres = aux(1:nunkP);
%             btot = [M*aux1(1:nunk)-dt*G'*pres; btot2nd];
%             aux1 = UM\(LM\btot);
%             velo = reshape(aux1(1:nunk),2,numnp)';
        %end
    end
        

    % Update variables for next iteration
    veloVect = reshape(velo',ndofV,1);
    sol0 = [veloVect(dofUnk); pres];
    
    solf(:,step) = sol0;
    figure(1); clf;
    PlotStreamlines(X,velo,dom);
    hold on
    %pause(0.5)
end
toc

if confined
    pres = [0; pres];
    
end

nPt = size(X,1); 
figure(2); 
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
