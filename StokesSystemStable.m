function [K,G,f,L,f_q] = StokesSystemStable(X,T,XP,TP,referenceElement,nu)
% [K,G,f] = StokesSystem(X,T,XP,TP,referenceElement)
% Matrices K, G and r.h.s vector f obtained after discretizing a Stokes problem
%
% X,T: nodal coordinates and connectivities for velocity
% XP,TP: nodal coordinates and connectivities for pressure
% referenceElement: reference element properties (quadrature, shape functions...)


elem = referenceElement.elemV;
ngaus = referenceElement.ngaus;
wgp = referenceElement.GaussWeights;
N = referenceElement.N;
Nxi = referenceElement.Nxi;
Neta = referenceElement.Neta;
NP = referenceElement.NP; 
ngeom = referenceElement.ngeom;
h=XP(2)-XP(1);
tau1 = 1/3*h^2/(4*nu);


% Number of elements and number of nodes in each element
[nElem,nenV] = size(T);
nenP = size(TP,2); 

% Number of nodes
nPt_V = size(X,1);
if elem == 11
    nPt_V = nPt_V + nElem; 
end
nPt_P = size(XP,1);

% Number of degrees of freedom 
nedofV = 2*nenV; 
nedofP = nenP;
ndofV = 2*nPt_V; 
ndofP = nPt_P; 

K = zeros(ndofV,ndofV);
G = zeros(ndofP,ndofV); 
f = zeros(ndofV,1);
L = zeros(ndofP,ndofP);
f_q=zeros(ndofP,1);

% Loop on elements
for ielem = 1:nElem
    % Global number of the nodes in element ielem
    Te = T(ielem,:);
    TPe = TP(ielem,:); 
    % Coordinates of the nodes in element ielem
    Xe = X(Te(1:ngeom),:);
    % Degrees of freedom in element ielem
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofV);
    TPe_dof = TPe; 
    
    % Element matrices
    [Ke,Ge,fe,Le,f_qe] = EleMatStokes(Xe,ngeom,nedofV,nedofP,ngaus,wgp,N,Nxi,Neta,NP,tau1);
    
    % Assemble the element matrices
    K(Te_dof, Te_dof) = K(Te_dof, Te_dof) + Ke;
    G(TPe_dof,Te_dof) = G(TPe_dof,Te_dof) + Ge; 
    f(Te_dof) = f(Te_dof) + fe;
    L(TPe_dof,TPe_dof)= L(TPe_dof,TPe_dof) + Le;
    f_q(TPe_dof) = f_q(TPe_dof) + f_qe;
end






function [Ke,Ge,fe,Le,f_qe] = EleMatStokes(Xe,ngeom,nedofV,nedofP,ngaus,wgp,N,Nxi,Neta,NP,tau1)
%

Ke = zeros(nedofV,nedofV);
Ge = zeros(nedofP,nedofV);
fe = zeros(nedofV,1);
Le = zeros(nedofP,nedofP);
f_qe=zeros(nedofP,1);


% Loop on Gauss points 
for ig = 1:ngaus
    N_ig    = N(ig,:);
    Nxi_ig  = Nxi(ig,:);
    Neta_ig = Neta(ig,:);
    NP_ig = NP(ig,:); 
    Jacob = [
        Nxi_ig(1:ngeom)*(Xe(:,1))	Nxi_ig(1:ngeom)*(Xe(:,2))
        Neta_ig(1:ngeom)*(Xe(:,1))	Neta_ig(1:ngeom)*(Xe(:,2))
        ];
    dvolu = wgp(ig)*det(Jacob);
    res = Jacob\[Nxi_ig;Neta_ig];
    nx = res(1,:);
    ny = res(2,:);
    
	Ngp = [reshape([1;0]*N_ig,1,nedofV); reshape([0;1]*N_ig,1,nedofV)];
    % Gradient
    Nx = [reshape([1;0]*nx,1,nedofV); reshape([0;1]*nx,1,nedofV)];
    Ny = [reshape([1;0]*ny,1,nedofV); reshape([0;1]*ny,1,nedofV)];
    %NPx = [reshape([1;0]*nx,1,nedofP); reshape([0;1]*nx,1,nedofP)];
    %NPy = [reshape([1;0]*ny,1,nedofP); reshape([0;1]*ny,1,nedofP)];
    % Divergence
    dN = reshape(res,1,nedofV);

    Ke = Ke + (Nx'*Nx+Ny'*Ny)*dvolu; 
    Ge = Ge - NP_ig'*dN*dvolu; 
    %Le = Le + tau1*([nx;nx;ny;ny]'*[nx;ny;nx;ny])*dvolu;
    Le = Le - tau1*(nx'*nx+ny'*ny)*dvolu; 
    x_ig = N_ig(1:ngeom)*Xe; 
    f_igaus = SourceTerm(x_ig); 
    fe = fe + Ngp'*f_igaus*dvolu;
    f_qe = f_qe - tau1*([nx; ny]'*f_igaus)*dvolu;
end

