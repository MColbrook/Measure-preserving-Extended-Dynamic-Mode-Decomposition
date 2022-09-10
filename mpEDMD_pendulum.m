%% Code for the pendulum example from 
%% ref: M. Colbrook, "The mpEDMD Algorithm for Data-Driven Computations of Measure-Preserving Dynamical Systems"
clear
close all

%% parameters
N=100;  M1=200; M2=M1;

%% setup noise-free matrices
[PSI_X_nf,PSI_Y_nf,G_nf,A_nf,W,x1,x2] = pendulum_matrix(M1,M2,N,5);

%% compute spectra
sigma=0; % percentage of Gaussian noise
N=size(G_nf,1);

% create noise perturbations
R1=randn(size(PSI_X_nf,1),size(PSI_X_nf,2)); R1=R1/std(R1(:))*std(PSI_X_nf(:));
R2=randn(size(PSI_X_nf,1),size(PSI_X_nf,2)); R2=R2/std(R2(:))*std(PSI_Y_nf(:));

R1X=(R1'*PSI_X_nf)*W(1);    R1R1=(R1'*R1)*W(1);
R1Y=(R1'*PSI_Y_nf)*W(1);    R1R2=(R1'*R2)*W(1);
XR2=(PSI_X_nf'*R2)*W(1);

% G and A matrices - see paper for notation
G=G_nf+sigma*(R1X+R1X')+sigma^2*R1R1;   A=A_nf+sigma*(R1Y+XR2)+sigma^2*R1R2;

% mpEDMD
[mpK,mpV,mpD] = mpEDMD(G,A);    mpE=diag(mpD);

% EDMD
K=G\A;  [V,D]=eig(K);           E=diag(D);  

% residuals
Res=zeros(N,1); mpRes=Res;
for jj=1:N
    Res(jj)=real((V(:,jj)'*(G_nf-E(jj)*A_nf'-conj(E(jj))*A_nf+abs(E(jj))^2*G_nf)*V(:,jj))/(V(:,jj)'*G_nf*V(:,jj)));
    mpRes(jj)=real((mpV(:,jj)'*(G_nf-mpE(jj)*A_nf'-conj(mpE(jj))*A_nf+abs(mpE(jj))^2*G_nf)*mpV(:,jj))/(mpV(:,jj)'*G_nf*mpV(:,jj)));
end

%% plot eigenvalues
figure
plot(real(mpE),imag(mpE),'g*','markersize',10)
hold on
plot(real(E),imag(E),'rx','markersize',10)
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k')
ax=gca; ax.FontSize=14;
legend({'\texttt{mpEDMD}','EDMD'},'interpreter','latex','fontsize',14,'location','best')
axis equal
axis([-1.2,1.2,-1.2,1.2])

%% plot approximate eigenfunctions
lam=exp(1i*pi/4); % lambda value for approximate eigenfunctions

I1=find(abs(E-lam)==min(abs(E-lam)));       C=(PSI_X_nf*V(:,I1))/sqrt(V(:,I1)'*G_nf*V(:,I1));
I2=find(abs(mpE-lam)==min(abs(mpE-lam)));   mpC=(PSI_X_nf*mpV(:,I2))/sqrt(mpV(:,I2)'*G_nf*mpV(:,I2));

figure
toPlot = reshape((mpC),length(x1),length(x2));
contourf(x1, x2 , log10(abs(toPlot)),60,'edgecolor','none');
colormap((turbo));
caxis([-2,0]);  colorbar
axis equal on;   view(0,90);    ylim([-4,4]);   xlim([-pi,pi])
ax=gca; ax.FontSize=14;
set(gca,'YDir','normal')
    
figure
toPlot = reshape((C),length(x1),length(x2));
contourf(x1, x2 , log10(abs(toPlot)),60,'edgecolor','none');
colormap((turbo));
caxis([-2,0]);  colorbar
axis equal on;   view(0,90);    ylim([-4,4]);   xlim([-pi,pi])
ax=gca; ax.FontSize=14;
set(gca,'YDir','normal')


function [PSI_X,PSI_Y,G,A,W,x1,x2] = pendulum_matrix(M1,M2,N,L)
% computes data matrices for the pendulum
delta_t=0.5;
ODEFUN=@(t,y) [y(2);-sin(y(1))];
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

x1=linspace(-pi,pi,M1+1);   x1=x1+(x1(2)-x1(1))/2;
x1=x1(1:end-1); x2=linspace(-L,L,M2);
[X1,X2] = meshgrid(x1,x2);  X1=X1(:); X2=X2(:);
M=length(X1); % number of data points

PSI_X=zeros(M,N+1);
PSI_X(:,1)=exp(1i*X1).*X2.*exp(-X2.^2/2);
XX0=PSI_X(:,1);

pf = parfor_progress(M);
pfcleanup = onCleanup(@() delete(pf));
parfor j=1:M
    Y0=[X1(j);X2(j)];
    [~,Y]=ode45(ODEFUN,[0.000001 (1:N)*delta_t],Y0,options);
    PSI_X(j,:)=[XX0(j,1),transpose(exp(1i*Y(2:end,1)).*Y(2:end,2).*exp(-Y(2:end,2).^2/2))];
    parfor_progress(pf);
end

PSI_Y=PSI_X(:,2:end);     PSI_X=PSI_X(:,1:end-1);
W=zeros(M,1)+(x1(2)-x1(1))*(x2(2)-x2(1));
G=(PSI_X'*PSI_X)*W(1);     A=(PSI_X'*PSI_Y)*W(1);
end
