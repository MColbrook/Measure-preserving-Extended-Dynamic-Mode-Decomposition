%% Code for the Lorenz example from 
%% ref: M. Colbrook, "The mpEDMD Algorithm for Data-Driven Computations of Measure-Preserving Dynamical Systems"
clear
close all

%% Set parameters
M=10^4; % number of data points
delta_t=0.1; % time step
SIGMA=10;   BETA=8/3;   RHO=28;
ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
N=1000; % number of delay embeddings for reference

%% Produce the data
Y0=(rand(3,1)-0.5)*4;
[~,Y]=ode45(ODEFUN,[0.000001 (1:(10000+2*M))*delta_t],Y0,options);
Y=Y(10001:end,:); % sample after when on the attractor

%% Scalar-valued spectral measure
PSI_X1=zeros(M,N+1);   PSI_X1(:,1)=Y(1:M,1);
PSI_X2=zeros(M,N+1);   PSI_X2(:,1)=Y(1:M,2);
PSI_X3=zeros(M,N+1);   PSI_X3(:,1)=Y(1:M,3);

pf = parfor_progress(N);
pfcleanup = onCleanup(@() delete(pf));
for j=2:N+1
    PSI_X1(:,j)=Y((1:M)+j,1);    PSI_X2(:,j)=Y((1:M)+j,2);    PSI_X3(:,j)=Y((1:M)+j,3);
    parfor_progress(pf);
end
%% G and A matrices
G1=(PSI_X1(:,1:N)'*PSI_X1(:,1:N))/M;    A1=(PSI_X1(:,1:N)'*PSI_X1(:,2:(N+1)))/M;
G2=(PSI_X2(:,1:N)'*PSI_X2(:,1:N))/M;    A2=(PSI_X2(:,1:N)'*PSI_X2(:,2:(N+1)))/M;
G3=(PSI_X3(:,1:N)'*PSI_X3(:,1:N))/M;    A3=(PSI_X3(:,1:N)'*PSI_X3(:,2:(N+1)))/M;

%% Compute mpEDMD matrices and spec measures
c=zeros(N,1); c(1)=1;

[~,mpV1,mpD1] = mpEDMD(G1,A1);
piE1=diag(mpD1); TH1=angle(piE1);
MU1=abs(mpV1'*G1*c).^2;

[~,mpV2,mpD2] = mpEDMD(G2,A2);
piE2=diag(mpD2); TH2=angle(piE2);
MU2=abs(mpV2'*G2*c).^2;

[~,mpV3,mpD3] = mpEDMD(G3,A3);
piE3=diag(mpD3);  TH3=angle(piE3);
MU3=abs(mpV3'*G3*c).^2;

%% Cdf plots
[~,Ib] = sort(TH1(:),'ascend');
THp1=TH1(Ib); THp1=[THp1(:)-10^(-14),THp1(:)]'; THp1=THp1(:);
cdf1=0*THp1;
cc=0;
for j=1:length(TH1)
    cdf1(2*j-1)=cc;    cc=cc+MU1(Ib(j));    cdf1(2*j)=cc;
end

[~,Ib] = sort(TH2(:),'ascend');
THp2=TH2(Ib); THp2=[THp2(:)-10^(-14),THp2(:)]'; THp2=THp2(:);
cdf2=0*THp2;
cc=0;
for j=1:length(TH2)
    cdf2(2*j-1)=cc;    cc=cc+MU2(Ib(j));    cdf2(2*j)=cc;
end

[~,Ib] = sort(TH3(:),'ascend');
THp3=TH3(Ib); THp3=[THp3(:)-10^(-14),THp3(:)]'; THp3=THp3(:);
cdf3=0*THp3;
cc=0;
for j=1:length(TH3)
    cdf3(2*j-1)=cc;    cc=cc+MU3(Ib(j));    cdf3(2*j)=cc;
end

figure
plot(THp1,cdf1/sum(MU1),'linewidth',2)
hold on
plot(THp2,cdf2/sum(MU2),'linewidth',2)
plot(THp3,cdf3/sum(MU3),'linewidth',2)
ax = gca; ax.FontSize = 16; ylim([0,1]); xlim([-pi,pi]);
legend({'$g_1$','$g_2$','$g_3$'},'interpreter','latex','fontsize',16,'location','northwest')

%% Vector-valued spectral measure
GG=([PSI_X1(:,1:333),PSI_X2(:,1:333),PSI_X3(:,1:333)]'*[PSI_X1(:,1:333),PSI_X2(:,1:333),PSI_X3(:,1:333)])/M;
AA=([PSI_X1(:,1:333),PSI_X2(:,1:333),PSI_X3(:,1:333)]'*[PSI_X1(:,2:(333+1)),PSI_X2(:,2:(333+1)),PSI_X3(:,2:(333+1))])/M;

[~,mpVV,mpDD] = mpEDMD(GG,AA);
mpEE=diag(mpDD);
%%
phi=chebfun(@(t) exp(sin(t)),[-pi pi],'periodic'); % test function
NN=size(GG,1);
cc1=zeros(NN,1); cc1(1)=1;
dE1=mpVV*((mpVV'*GG*cc1).*phi(mpEE(:)));

cc2=zeros(999,1); cc2(NN/3+1)=1;
dE2=mpVV*((mpVV'*GG*cc2).*phi(mpEE(:)));

cc3=zeros(999,1); cc3(2*NN/3+1)=1;
dE3=mpVV*((mpVV'*GG*cc3).*phi(mpEE(:)));

%% Plot results (function on attractor)
C=[PSI_X1(1:min(2*10^4,M),1:NN/3),PSI_X2(1:min(2*10^4,M),1:NN/3),PSI_X3(1:min(2*10^4,M),1:NN/3)]*dE1;

figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter3(PSI_X1(1:min(2*10^4,M),1),PSI_X2(1:min(2*10^4,M),1),PSI_X3(1:min(2*10^4,M),1),3,real(C),'filled');
colormap turbo; colorbar; axis equal; view(axes1,[44.1307171302158 20.3999998682605]);
grid(axes1,'on'); axis(axes1,'tight'); hold(axes1,'off');
set(axes1,'DataAspectRatio',[1 1 1]);

C=[PSI_X1(1:min(2*10^4,M),1:NN/3),PSI_X2(1:min(2*10^4,M),1:NN/3),PSI_X3(1:min(2*10^4,M),1:NN/3)]*dE2;

figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter3(PSI_X1(1:min(2*10^4,M),1),PSI_X2(1:min(2*10^4,M),1),PSI_X3(1:min(2*10^4,M),1),3,real(C),'filled');
colormap turbo; colorbar; axis equal; view(axes1,[44.1307171302158 20.3999998682605]);
grid(axes1,'on'); axis(axes1,'tight'); hold(axes1,'off');
set(axes1,'DataAspectRatio',[1 1 1]);

C=[PSI_X1(1:min(2*10^4,M),1:NN/3),PSI_X2(1:min(2*10^4,M),1:NN/3),PSI_X3(1:min(2*10^4,M),1:NN/3)]*dE3;

figure1=figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter3(PSI_X1(1:min(2*10^4,M),1),PSI_X2(1:min(2*10^4,M),1),PSI_X3(1:min(2*10^4,M),1),3,real(C),'filled');
colormap turbo; colorbar; axis equal; view(axes1,[44.1307171302158 20.3999998682605]);
grid(axes1,'on'); axis(axes1,'tight'); hold(axes1,'off');
set(axes1,'DataAspectRatio',[1 1 1]);


function [Z] = W1(TH1,MU1,TH2,MU2)
% Computes W_1 distance of two discrete measures

break_pts=[TH1(:);TH2(:)];
ind=[TH1(:)*0;TH2(:)*0+1];
MU=[MU1(:);MU2(:)];

[break_pts,I] = sort(break_pts,'ascend');
ind=ind(I); MU=MU(I);

L=length(I);    F1=0;   F2=0;

if ind(1)==0
    F1=F1+MU(1);
else
    F2=F2+MU(1);
end

Z=0;
for j=2:L
    Z=Z+(break_pts(j)-break_pts(j-1))*abs(F1-F2);
    if ind(j)==0
        F1=F1+MU(j);
    else
        F2=F2+MU(j);
    end
end
end


