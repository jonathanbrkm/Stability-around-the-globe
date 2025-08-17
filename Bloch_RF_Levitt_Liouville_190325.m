clear all; close all; clc;

om=1;
T=2*pi/om;

N=1E2;
N1=N/4; N2=N/2;
tau1=T/4; tau2=T/2;
beta1=tau1; beta2=tau2;
phi1=0; phi2=pi/2;

a0_vec=[0:0.002:0.5].'; phi0_vec=[0:0.002:0.5].*2*pi.';
[P0,A0]=meshgrid(phi0_vec,a0_vec);

E0=1; delta=0;

% figure;
% sphere(25); view(150,20);
% set(gcf,'color','w');
% axis equal;
% shading interp; colormap gray; alpha 0.25; 
% line([-1 1],[0 0],[0 0],'LineWidth',1,'Color',[0 0 0]);
% line([0 0],[-1 1],[0 0],'LineWidth',1,'Color',[0 0 0]);
% line([0 0],[0 0],[-1 1],'LineWidth',1,'Color',[0 0 0]);
% xlabel('$x$','interpreter','latex','fontsize',20,'position',[1.1 0 0]);
% ylabel('$y$','interpreter','latex','fontsize',20,'position',[0 1.1 0]);
% zlabel('$z$','interpreter','latex','fontsize',20,'position',[0 0 1.1]); h=gca; h.ZAxis.Label.Rotation=0;
% hold on;
% 
% theta_aux1=linspace(0,2*pi,100);
% x1=cos(theta_aux1);
% z1=sin(theta_aux1);
% y1=zeros(size(theta_aux1));
% 
% theta_aux2=pi/2;
% x2=cos(theta_aux1)*sin(theta_aux2);
% y2=sin(theta_aux1)*sin(theta_aux2);
% z2=cos(theta_aux2)*ones(size(theta_aux1));
% 
% theta_aux3=pi/2;
% x3=cos(theta_aux3)*ones(size(theta_aux1));
% y3=cos(theta_aux1)*sin(theta_aux3);
% z3=sin(theta_aux1)*sin(theta_aux3);
% 
% plot3(x1,y1,z1,'--k','Color',[0.5 0.5 0.5]); hold on;
% plot3(x2,y2,z2,'--k','Color',[0.5 0.5 0.5]); hold on;
% plot3(x3,y3,z3,'--k','Color',[0.5 0.5 0.5]); hold on;

theta=pi/2;
t=[0:1/N:1-1/N].';

for i=1:length(a0_vec)
i
for j=1:length(phi0_vec)

nx1=sin(theta).*cos(phi1);
ny1=sin(theta).*sin(phi1);
nz1=cos(theta);
K1=[0,-nz1,ny1;nz1,0,-nx1;-ny1,nx1,0];

nx2=sin(theta).*cos(phi2);
ny2=sin(theta).*sin(phi2);
nz2=cos(theta);
K2=[0,-nz2,ny2;nz2,0,-nx2;-ny2,nx2,0];

Rx_90=expm(-(beta1./N1).*K1);
Ry_180=expm(-(beta2./N2).*K2);

r=zeros(3,N);

a0=A0(i,j).*exp(1i.*P0(i,j)); b0=sqrt(1-abs(a0).^2);
r(:,1)=[a0.*conj(b0)+b0.*conj(a0),1i.*(conj(a0).*b0-conj(b0).*a0),abs(b0).^2-abs(a0).^2].';

for f=1:N1
r(:,f+1)=Rx_90*r(:,f);
end

for f=N1:(N1+N2)
r(:,f+1)=Ry_180*r(:,f);
end

for f=(N1+N2):(N1+N2+N1)-1
r(:,f+1)=Rx_90*r(:,f);
end

xt=r(1,:); yt=r(2,:); zt=r(3,:);
% plot3(xt,yt,zt);
% drawnow;

X(i,j,:)=xt;
Y(i,j,:)=yt;
Z(i,j,:)=zt;

clear r;

end
end

Phi=atan2(Y,X);
Eta=Z;

Area=zeros(length(t),1);
%%
vidfile=VideoWriter('Bloch_RF_Levitt_Liouville_190325.mp4','MPEG-4');
vidfile.FrameRate=10;
open(vidfile);

figure;
set(gcf,'color','w');
for f=1:length(t)
plot(Phi(:,:,f),Eta(:,:,f),'.k');
xlabel('$\phi$','interpreter','latex','fontsize',20);
ylabel('$\eta$','interpreter','latex','fontsize',20);
title(['$t=$',num2str(round(t(f),3)),'$T$'],'interpreter','latex','fontsize',20);
axis([-4 4 -1 1]);
% axis off;
drawnow;
BW=getframe(gcf);
BW=BW.cdata(:,:,1);
Area(f)=bwarea(BW);
writeVideo(vidfile,getframe(gcf));
end
close(vidfile);
%%
figure;
set(gcf,'color','w'); grid on;
plot(t,Area./max(Area),'k','LineWidth',2); ylim([0.99 1.01]);
xlabel('$t/T$','interpreter','latex','fontsize',20);
ylabel('$\rho(t)$','interpreter','latex','fontsize',20);