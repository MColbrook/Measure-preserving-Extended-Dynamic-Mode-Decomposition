%% Code for the turbulent flow example from 
%% ref: M. Colbrook, "The mpEDMD Algorithm for Data-Driven Computations of Measure-Preserving Dynamical Systems"
clear
close all

%% Data set (see link in README file)
load('mpEDMD_turbulent_data')

%% Algorithmic parameters - NB: parameters of paper were larger
N=100;%999;
LL=200; % number of timesteps for prediction
sam_num=10; % number of random initial conditions for TKE plot

%%
M=N;
ind1=1:M;
Fsamp=floor(12000/(ind1(end)+1));
DATA=DATA2(:,1:Fsamp:end);

%% Compute measurement matrices
[~,~,V]=svd(transpose(DATA(:,ind1))/sqrt(M),'econ');
PSI_X=transpose(DATA(:,ind1))*V(:,1:N);
PSI_Y=transpose(DATA(:,ind1+1))*V(:,1:N);

clear S V

G=(PSI_X(1:M,1:N)'*PSI_X(1:M,1:N))/M;    
A=(PSI_X(1:M,1:N)'*PSI_Y(1:M,1:N))/M;

% mpEDMD
[~,V_mpEDMD,D_mpEDMD] = mpEDMD(G,A);
E_mpEDMD=diag(D_mpEDMD);

% EDMD
K_EDMD=G\A;
[V_EDMD,D_EDMD]=eig(K_EDMD);
E_EDMD=diag(D_EDMD);

% piDMD
[U,~,V] = svd(transpose(PSI_Y)*transpose(PSI_X)');
K_piDMD=transpose(U*V');
clear U V
[V_piDMD,D_piDMD]=eig(K_piDMD);
E_piDMD=diag(D_piDMD);

%% Koopman mode decomposition
sigu=DATA(1:51150,ind1);
sigv=DATA(51151:end,ind1);

uMode_EDMD=pinv(PSI_X*V_EDMD)*transpose(sigu);
vMode_EDMD=pinv(PSI_X*V_EDMD)*transpose(sigv);

uMode_mpEDMD=pinv(PSI_X*V_mpEDMD)*transpose(sigu);
vMode_mpEDMD=pinv(PSI_X*V_mpEDMD)*transpose(sigv);

uMode_piDMD=pinv(PSI_X*V_piDMD)*transpose(sigu);
vMode_piDMD=pinv(PSI_X*V_piDMD)*transpose(sigv);

%%
DATA_EDMD=zeros(size(DATA,1),LL);   DATA_mpEDMD=DATA_EDMD;  DATA_piDMD=DATA_EDMD;
YYY_EDMD=[uMode_EDMD,vMode_EDMD];
YYY_mpEDMD=[uMode_mpEDMD,vMode_mpEDMD];
YYY_piDMD=[uMode_piDMD,vMode_piDMD];

uf=zeros(310,165,size(DATA,2));    vf=uf;
for j=1:size(DATA,2)
    uf(:,:,j)=real(reshape(transpose(DATA(1:51150,j)),[310,165]));
    vf(:,:,j)=real(reshape(transpose(DATA(51151:end,j)),[310,165]));
end
[TKE,Ruu,Psi11_y,k_x,NNx,lambda_x] = flow_en_stats(uf,vf,x);

SS=randsample(length(ind1),sam_num); % random initial conditions to sample over

for jjj=1:length(SS)
    jjj
    XXX_EDMD=PSI_X(SS(jjj),:)*V_EDMD;
    XXX_mpEDMD=PSI_X(SS(jjj),:)*V_mpEDMD;
    XXX_piDMD=PSI_X(SS(jjj),:)*V_piDMD;

    pf = parfor_progress(LL);
    pfcleanup = onCleanup(@() delete(pf));
    parfor j=1:LL
        DATA_EDMD(:,j)=transpose(XXX_EDMD*(sparse(diag(E_EDMD.^(j-1))))*YYY_EDMD);
        DATA_mpEDMD(:,j)=transpose(XXX_mpEDMD*(sparse(diag(E_mpEDMD.^(j-1))))*YYY_mpEDMD);
        DATA_piDMD(:,j)=transpose(XXX_piDMD*(sparse(diag(E_piDMD.^(j-1))))*YYY_piDMD);
        parfor_progress(pf);
    end
    
    % compute flow fields
    uf_EDMD=0*uf;    vf_EDMD=0*uf;
    uf_mpEDMD=0*uf;  vf_mpEDMD=0*uf;
    uf_piDMD=0*uf;   vf_piDMD=0*uf;

    pf = parfor_progress(LL);
    pfcleanup = onCleanup(@() delete(pf));
    
    for j=1:LL
        uf_EDMD(:,:,j)=real(reshape(transpose(DATA_EDMD(1:51150,j)),[310,165]));
        vf_EDMD(:,:,j)=real(reshape(transpose(DATA_EDMD(51151:end,j)),[310,165]));
        uf_mpEDMD(:,:,j)=real(reshape(transpose(DATA_mpEDMD(1:51150,j)),[310,165]));
        vf_mpEDMD(:,:,j)=real(reshape(transpose(DATA_mpEDMD(51151:end,j)),[310,165]));
        uf_piDMD(:,:,j)=real(reshape(transpose(DATA_piDMD(1:51150,j)),[310,165]));
        vf_piDMD(:,:,j)=real(reshape(transpose(DATA_piDMD(51151:end,j)),[310,165]));
        parfor_progress(pf);
    end
    
    [TKE_EDMD0,Ruu_EDMD0,Psi11_y_EDMD0,~,~,~] = flow_en_stats(uf_EDMD,vf_EDMD,x);
    [TKE_mpEDMD0,Ruu_mpEDMD0,Psi11_y_mpEDMD0,~,~,~] = flow_en_stats(uf_mpEDMD,vf_mpEDMD,x);
    [TKE_piDMD0,Ruu_piDMD0,Psi11_y_piDMD0,~,~,~] = flow_en_stats(uf_piDMD,vf_piDMD,x);

    if jjj==1
        TKE_EDMD=TKE_EDMD0; Ruu_EDMD=Ruu_EDMD0; Psi11_y_EDMD=Psi11_y_EDMD0;
        TKE_mpEDMD=TKE_mpEDMD0; Ruu_mpEDMD=Ruu_mpEDMD0; Psi11_y_mpEDMD=Psi11_y_mpEDMD0;
        TKE_piDMD=TKE_piDMD0; Ruu_piDMD=Ruu_piDMD0; Psi11_y_piDMD=Psi11_y_piDMD0;              
    else
        TKE_EDMD=(TKE_EDMD*(jjj-1)+TKE_EDMD0)/jjj; Ruu_EDMD=(Ruu_EDMD*(jjj-1)+Ruu_EDMD0)/jjj; Psi11_y_EDMD=(Psi11_y_EDMD*(jjj-1)+Psi11_y_EDMD0)/jjj;
        TKE_mpEDMD=(TKE_mpEDMD*(jjj-1)+TKE_mpEDMD0)/jjj; Ruu_mpEDMD=(Ruu_mpEDMD*(jjj-1)+Ruu_mpEDMD0)/jjj; Psi11_y_mpEDMD=(Psi11_y_mpEDMD*(jjj-1)+Psi11_y_mpEDMD0)/jjj;
        TKE_piDMD=(TKE_piDMD*(jjj-1)+TKE_piDMD0)/jjj; Ruu_piDMD=(Ruu_piDMD*(jjj-1)+Ruu_piDMD0)/jjj; Psi11_y_piDMD=(Psi11_y_piDMD*(jjj-1)+Psi11_y_piDMD0)/jjj;
    end

end

%% TKE plots
III=find(abs(y-5)==min(abs(y-5)));

figure
semilogy((1:LL)/1000,squeeze(mean(mean(TKE_mpEDMD(:,III,:),2),1))/2500,'g','linewidth',2); hold on
plot((1:LL)/1000,squeeze(mean(mean(TKE_piDMD(:,III,:),2),1))/2500,'b','linewidth',2)
plot((1:LL)/1000,squeeze(mean(mean(TKE_EDMD(:,III,:),2),1))/2500,'r','linewidth',2)
semilogy((1:LL)/1000,(1:LL)*0+mean(squeeze(mean(mean(TKE(:,III,:),2),1)))/2500,'k','linewidth',2)
semilogy((3500:LL)/1000,max(abs(E_EDMD)).^(2*(3500:LL))/1000000,'--k','linewidth',2)
ax = gca; ax.FontSize = 16;
legend({'\texttt{mpEDMD}','piDMD','EDMD','Mean TKE of flow'},'interpreter','latex','fontsize',16,'location','northwest')

III=find(abs(y-35)==min(abs(y-35)));

figure
semilogy((1:LL)/1000,squeeze(mean(mean(TKE_mpEDMD(:,III,:),2),1))/2500,'g','linewidth',2); hold on
plot((1:LL)/1000,squeeze(mean(mean(TKE_piDMD(:,III,:),2),1))/2500,'b','linewidth',2)
plot((1:LL)/1000,squeeze(mean(mean(TKE_EDMD(:,III,:),2),1))/2500,'r','linewidth',2)
semilogy((1:LL)/1000,(1:LL)*0+mean(squeeze(mean(mean(TKE(:,III,:),2),1)))/2500,'k','linewidth',2)
semilogy((3500:LL)/1000,max(abs(E_EDMD)).^(2*(3500:LL))/800000,'--k','linewidth',2)
ax = gca; ax.FontSize = 16;
legend({'\texttt{mpEDMD}','piDMD','EDMD','Mean TKE of flow'},'interpreter','latex','fontsize',16,'location','northwest')

figure
plot(y,mean(mean(TKE_mpEDMD,3),1)/2500,'g','linewidth',4); hold on
plot(y,mean(mean(TKE_piDMD,3),1)/2500,'b','linewidth',4)
plot(y,mean(mean(TKE,3),1)/2500,'--k','linewidth',2)
xlim([0,max(y)]);   ax = gca; ax.FontSize = 16;
legend({'\texttt{mpEDMD}','piDMD','Mean TKE of flow'},'interpreter','latex','fontsize',16,'location','best')

%% Wavenumber spectra plots
figure
toPlot = (k_x(NNx).*abs(Psi11_y(:,NNx))./(mean(mean(Ruu,2)) )  );
LLL=max(toPlot(:));
contourf(k_x(NNx), y , toPlot ,40,'edgecolor','none')
set(gca,'xscale','log');
colormap(turbo);    caxis([0 LLL]);    colorbar;
ax = gca; ax.FontSize = 16;

figure
toPlot = (k_x(NNx).*abs(Psi11_y_mpEDMD(:,NNx))./(mean(mean(Ruu_mpEDMD,2)) )  );
contourf(k_x(NNx), y , toPlot ,40,'edgecolor','none')
set(gca,'xscale','log');
colormap(turbo);    caxis([0 LLL]);    colorbar;
ax = gca; ax.FontSize = 16;

figure
toPlot = (k_x(NNx).*abs(Psi11_y_piDMD(:,NNx))./(mean(mean(Ruu_mpEDMD,2)) )  );
contourf(k_x(NNx), y , toPlot ,40,'edgecolor','none')
set(gca,'xscale','log');
colormap(turbo);    caxis([0 LLL]);    colorbar;
ax = gca; ax.FontSize = 16;

figure
toPlot = (k_x(NNx).*abs(Psi11_y_EDMD(:,NNx))./(mean(mean(Ruu_EDMD,2)) )  );
contourf(k_x(NNx), y , toPlot ,40,'edgecolor','none')
set(gca,'xscale','log');
colormap(turbo);    caxis([0 LLL]);    
colorbar;
ax = gca; ax.FontSize = 16;

%% Flow plots
figure
XI1=uf(:,:,min(500,size(uf,3)));
myContours = linspace(min(real(XI1(:))),max(real(XI1(:))), 21);
contourf(Xgrid+1282.7,Ygrid,(real(XI1)),myContours,'edgecolor','none')
set(gca,'ydir','normal')
colorbar
colormap(brighten(redblueTecplot(21),-0.6))
ax = gca; ax.FontSize = 16;

figure
XI1=uf_mpEDMD(:,:,min(4000,LL));
myContours = linspace(min(real(XI1(:))),max(real(XI1(:))), 21);
contourf(Xgrid+1282.7,Ygrid,(real(XI1)),myContours,'edgecolor','none')
set(gca,'ydir','normal')
colorbar
colormap(brighten(redblueTecplot(21),-0.6))
ax = gca; ax.FontSize = 16;

figure
XI1=uf_piDMD(:,:,min(4000,LL));
myContours = linspace(min(real(XI1(:))),max(real(XI1(:))), 21);
contourf(Xgrid+1282.7,Ygrid,(real(XI1)),myContours,'edgecolor','none')
set(gca,'ydir','normal')
colorbar
colormap(brighten(redblueTecplot(21),-0.6))
ax = gca; ax.FontSize = 16;

figure
XI1=uf_EDMD(:,:,min(4000,LL));
myContours = linspace(min(real(XI1(:))),max(real(XI1(:))), 21);
contourf(Xgrid+1282.7,Ygrid,(real(XI1)),myContours,'edgecolor','none')
set(gca,'ydir','normal')
colorbar
colormap(brighten(redblueTecplot(21),-0.6))
ax = gca; ax.FontSize = 16;

    
function [TKE,Ruu,Psi11_y,k_x,NNx,lambda_x] = flow_en_stats(uf,vf,x)
%% Define grid quantities
dX=abs(x(1)-x(2));
TKE = (uf.*uf+vf.*vf)/2; % turbulent kinectic energy

%% Two-point autocorrelation R_uu(deltaX,y) and R_vv(deltaX,y)
nY=165;   nX=310; nt=size(uf,3);
R11_y = zeros(nY, 2*nX-1);
parfor jj = 1:nY
    Rtmp1 = 0;
    for ii = 1:nt
        [r1,~] = xcorr(squeeze(uf(:,jj,ii)),'biased');
        Rtmp1 = Rtmp1 + r1;
    end
    R11_y(jj,:) = Rtmp1/nt;
end

%% FFT in the dx direction
[nR1, nR2] = size(R11_y);
Psi11_y=zeros(nR1,nR2);
parfor i = 1:nR1
    Psi11_y(i,:) = fftshift(fftn(((R11_y(i,:)))))./(nR2);
end
k_x = 2*pi*(-nR2/2:nR2/2-1)*((1/dX)/nR2);
NNx = (round(nR2/2,0)+1):nR2;
lambda_x = 2*pi./k_x(NNx);

Ruu = nanmean(uf.*uf,3);

end







