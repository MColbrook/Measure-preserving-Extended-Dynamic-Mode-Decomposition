function [mpK,mpV,mpD] = mpEDMD(G,A)
% mpEDMD algorithm
% input: G (Gram) and A (<Kpsi_j.psi_k>) matrices
% output:Koopman matrix mpK and its eigendecomposition [mpV,mpD]
% reference: M. Colbrook, "The mpEDMD Algorithm for Data-Driven Computations of Measure-Preserving Dynamical Systems"

G=(G+G')/2; % make sure G is Hermitian
[VG,DG]=eig(G);
Gsqrt=VG*diag(sqrt(diag(DG)))*VG'; % G^{1/2}
GsqrtI=VG*diag(sqrt(1./diag(DG)))*VG'; % G^{-1/2}

[U,~,V] = svd(GsqrtI*A'*GsqrtI);
[mpV,mpD]= schur(V*U','complex'); % schur usedto ensure orthonormal basis

mpV=GsqrtI*mpV;
mpK=GsqrtI*V*U'*Gsqrt;
mpD=diag(diag(mpD));

end

