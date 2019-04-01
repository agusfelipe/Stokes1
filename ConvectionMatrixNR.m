function [C1,C2] = ConvectionMatrixNR(X,T,referenceElement,velo)

% X,T: nodal coordinates and connectivities for velocity
% XP,TP: nodal coordinates and connectivities for pressure
% referenceElement: reference element properties (quadrature, shape functions...)

elem = referenceElement.elemV;
ngaus = referenceElement.ngaus;
wgp = referenceElement.GaussWeights;
N = referenceElement.N;
Nxi = referenceElement.Nxi;
Neta = referenceElement.Neta;
ngeom = referenceElement.ngeom; 

% Number of elements and number of nodes in each element
[nElem,nenV] = size(T);

% Number of nodes
nPt_V = size(X,1);
if elem == 11
    nPt_V = nPt_V + nElem; 
end

% Number of degrees of freedom 
nedofV = 2*nenV; 
ndofV = 2*nPt_V; 

C1 = zeros(ndofV,ndofV);
C2 = zeros(ndofV,ndofV);

for ielem = 1:nElem
    % Global number of the nodes in element ielem
    Te = T(ielem,:);
    u_e = velo(Te(1:ngeom),:);
    Ce1 = zeros(nedofV,nedofV);
    Ce2 = zeros(nedofV,nedofV);
    Xe = X(Te(1:ngeom),:);
    Te_dof = reshape([2*Te-1; 2*Te],1,nedofV);
       
    for ig=1:ngaus
    N_ig    = N(ig,:);
    Nxi_ig  = Nxi(ig,:);
    Neta_ig = Neta(ig,:);
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
    
    v_igaus = N_ig*u_e;
    
    
    Jacob_v = [
        Nxi_ig(1:ngeom)*(u_e(:,1))	Nxi_ig(1:ngeom)*(u_e(:,2))
        Neta_ig(1:ngeom)*(u_e(:,1))	Neta_ig(1:ngeom)*(u_e(:,2))
        ];
    
    Ce1 = Ce1 + Ngp'*(v_igaus(1)*Nx+v_igaus(2)*Ny)*dvolu; 
    Ce2 = Ce2 + Ngp'*Jacob_v'*Ngp*dvolu; 
    end
    C1(Te_dof,Te_dof) = C1(Te_dof,Te_dof) + Ce1;
    C2(Te_dof,Te_dof) = C2(Te_dof,Te_dof) + Ce2;
    
end