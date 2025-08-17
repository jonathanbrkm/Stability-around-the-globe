clear all; close all; clc;

om=1;
T=2*pi/om;

N0=1E3;
N10=N0/4; N20=N0/2;
tau1=T/4; tau2=T/2;
beta10=tau1; beta20=tau2;
phi1=0; phi2=pi/2;

E00=1;
E0_vec=[0.8:0.001:0.9].*E00.'; 
delta=0;

figure;
sphere(25); view(150,20);
set(gcf,'color','w');
axis equal;
shading interp; colormap gray; alpha 0.25; 
line([-1 1],[0 0],[0 0],'LineWidth',1,'Color',[0 0 0]);
line([0 0],[-1 1],[0 0],'LineWidth',1,'Color',[0 0 0]);
line([0 0],[0 0],[-1 1],'LineWidth',1,'Color',[0 0 0]);
xlabel('$x$','interpreter','latex','fontsize',20,'position',[1.1 0 0]);
ylabel('$y$','interpreter','latex','fontsize',20,'position',[0 1.1 0]);
zlabel('$z$','interpreter','latex','fontsize',20,'position',[0 0 1.1]); h=gca; h.ZAxis.Label.Rotation=0;
hold on;

theta_aux1=linspace(0,2*pi,100);
x1=cos(theta_aux1);
z1=sin(theta_aux1);
y1=zeros(size(theta_aux1));

theta_aux2=pi/2;
x2=cos(theta_aux1)*sin(theta_aux2);
y2=sin(theta_aux1)*sin(theta_aux2);
z2=cos(theta_aux2)*ones(size(theta_aux1));

theta_aux3=pi/2;
x3=cos(theta_aux3)*ones(size(theta_aux1));
y3=cos(theta_aux1)*sin(theta_aux3);
z3=sin(theta_aux1)*sin(theta_aux3);

plot3(x1,y1,z1,'--k','Color',[0.5 0.5 0.5]); hold on;
plot3(x2,y2,z2,'--k','Color',[0.5 0.5 0.5]); hold on;
plot3(x3,y3,z3,'--k','Color',[0.5 0.5 0.5]); hold on;

N_vec=zeros(length(E0_vec),1); 
Nmax_factor=sqrt((delta.^2+(max(E0_vec)).^2)./(E00.^2));
N_max=round(N0.*Nmax_factor);
Nmin_factor=sqrt((delta.^2+(min(E0_vec)).^2)./(E00.^2));
N_min=round(N0.*Nmin_factor);

ti=0;
dt=1/N_max;
tf=1-dt;
t=[0:dt:tf];

N1_min=N_min/4;
N2_min=N_min/2;
N1_max=N_max/4;
N2_max=N_max/2;

X=NaN.*zeros(length(E0_vec),N_max);
Y=NaN.*zeros(length(E0_vec),N_max);
Z=NaN.*zeros(length(E0_vec),N_max);
Energy=NaN.*zeros(length(E0_vec),N_max);

theta_vec=zeros(length(E0_vec),1);
beta1_vec=zeros(length(E0_vec),1);
beta2_vec=zeros(length(E0_vec),1);

N1_vec=zeros(length(E0_vec),1);
N2_vec=zeros(length(E0_vec),1);


for d=1:length(E0_vec)
    d
E0=E0_vec(d);

N1=round(N10.*sqrt((E0.^2+delta.^2)./(E00.^2)));
N2=round(N20.*sqrt((E0.^2+delta.^2)./(E00.^2)));

N=N1+N2+N1;
N_vec(d)=N;
N1_vec(d)=N1;
N2_vec(d)=N2;

td=linspace(0,2*pi,N);
V=-(E0/2).*(exp(-1i.*phi1).*(td>=ti).*(td<=tau1)+exp(-1i.*phi2).*(td>tau1).*(td<=(tau1+tau2))+exp(-1i.*phi1).*(td>(tau1+tau2)).*(td<=(tau1+tau2+tau1)));
Omega_x=V+conj(V);
Omega_y=1i.*(V-conj(V));
Omega_z=delta;

beta1=beta10.*sqrt((E0.^2+delta.^2)./(E00.^2));
beta2=beta20.*sqrt((E0.^2+delta.^2)./(E00.^2));
beta1_vec(d)=beta1;
beta2_vec(d)=beta2;

theta=atan2(E0,delta);
theta_vec(d)=theta;

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

a0=0; b0=sqrt(1-abs(a0).^2);
% r(:,1)=[2.*a0.*b0,0,abs(b0).^2-abs(a0).^2].';
r(:,1)=[conj(a0).*b0+conj(b0).*a0,1i.*(conj(a0).*b0-a0.*conj(b0)),abs(b0).^2-abs(a0).^2].';

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
Energy(d,1:N)=Omega_x.*xt+Omega_y.*yt+Omega_z.*zt;

X(d,1:N)=xt;
Y(d,1:N)=yt;
Z(d,1:N)=zt;

clear r;

end

Eta=Z;
Phi=atan2(Y,X);
%%
x_initial=X(sub2ind(size(X),(1:length(E0_vec)).',N1_vec));
y_initial=Y(sub2ind(size(Y),(1:length(E0_vec)).',N1_vec));
z_initial=Z(sub2ind(size(Z),(1:length(E0_vec)).',N1_vec));

eta_initial=Eta(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec));
phi_initial=Phi(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec));

x_final=X(sub2ind(size(X),(1:length(E0_vec)).',N_vec));
y_final=Y(sub2ind(size(Y),(1:length(E0_vec)).',N_vec));
z_final=Z(sub2ind(size(Z),(1:length(E0_vec)).',N_vec));

eta_final=Eta(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec));
phi_final=Phi(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec));

% plot3(X,Y,Z,'.k','MarkerSize',0.5);
% hold on;
% plot3(x_initial,y_initial,z_initial,'.c','MarkerSize',10);
% hold on;
% plot3(x_final,y_final,z_final,'.m','MarkerSize',10);
% drawnow;
%%
M11=zeros(length(E0_vec),N_max-N1_max);
M12=zeros(length(E0_vec),N_max-N1_max);
M21=zeros(length(E0_vec),N_max-N1_max);
M22=zeros(length(E0_vec),N_max-N1_max);
lambda1=zeros(length(E0_vec),N_max-N1_max);
lambda2=zeros(length(E0_vec),N_max-N1_max);
DetM=zeros(length(E0_vec),N_max-N1_max);

for s=1:(N_max-N1_max)
s

x_initial=X(sub2ind(size(X),(1:length(E0_vec)).',N1_vec+s));
y_initial=Y(sub2ind(size(Y),(1:length(E0_vec)).',N1_vec+s));
z_initial=Z(sub2ind(size(Z),(1:length(E0_vec)).',N1_vec+s));

% width(s)=sum((1/length(E0_vec)).*((x_initial-x_initial(51)).^2+(y_initial-y_initial(51)).^2+(z_initial-z_initial(51)).^2));
% % width(s)=(max(x_initial)-min(x_initial)).^2+(max(y_initial)-min(y_initial)).^2+(max(z_initial)-min(z_initial)).^2
width(s)=sqrt((1./length(E0_vec)).*sum((x_initial-mean(x_initial)).^2+(y_initial-mean(y_initial)).^2+(z_initial-mean(z_initial)).^2));

% eta_initial=z_initial;
% phi_initial=atan2(y_initial,x_initial);

% plot(phi_initial,eta_initial,'.k','MarkerSize',20);
% drawnow;
% BW=getframe(gcf);
% BW=BW.cdata(:,:,1);
% Area(s)=bwarea(BW);

% M11(:,s)=gradient(eta_final)./gradient(eta_initial);
% M12(:,s)=gradient(eta_final)./gradient(phi_initial);
% M21(:,s)=gradient(phi_final)./gradient(eta_initial);
% M22(:,s)=gradient(phi_final)./gradient(phi_initial);
 
% for d=1:length(E0_vec)
% J=[M11(d,s),M12(d,s);M21(d,s),M22(d,s)];
% lambda=eig(J);
% lambda1(d,s)=lambda(1);
% lambda2(d,s)=lambda(2);
% end
% 
% DetM(:,s)=M11(:,s).*M22(:,s)-M12(:,s).*M21(:,s);

% set(gcf,'color','w');
% plot(E0_vec,M11(:,s),'r',E0_vec,M22(:,s),'b',E0_vec,M11(:,s).*M22(:,s),'k','LineWidth',2);
% xlabel('$\mathcal{E}_{0}$','interpreter','latex','fontsize',20);
% legend('$\frac{\partial\eta_{f}}{\partial\eta_{i}}$','$\frac{\partial\phi_{f}}{\partial\phi_{i}}$','$\frac{\partial\eta_{f}}{\partial\eta_{i}}\frac{\partial\phi_{f}}{\partial\phi_{i}}$','interpreter','latex','fontsize',14);
% drawnow;
end


%%
vidfile=VideoWriter('Bloch_RF_Levitt_field_ih_stability_matrix_190325.mp4','MPEG-4');
vidfile.FrameRate=10;
open(vidfile);
figure;
set(gcf,'color','w');

figure;
set(gcf,'color','w','units','normalized','outerposition',[0 0 1 1]);

for s=1:100:(N_max-N1_max)
sgtitle(['$t=$',num2str(round(0.25+s/N,2)),'$T$'],'interpreter','latex','fontsize',20);
subplot(2,2,1);
plot(E0_vec,M11(:,s),'r');
xlabel('$\mathcal{E}_{0}$','interpreter','latex','fontsize',14);
ylabel('$\partial\eta_{f}/\partial\eta_{i}$','interpreter','latex','fontsize',14);
drawnow;
subplot(2,2,2);
plot(E0_vec,M12(:,s),'g');
xlabel('$\mathcal{E}_{0}$','interpreter','latex','fontsize',14);
ylabel('$\partial\eta_{f}/\partial\phi_{i}$','interpreter','latex','fontsize',14);
drawnow;
subplot(2,2,3);
plot(E0_vec,M21(:,s),'b');
xlabel('$\mathcal{E}_{0}$','interpreter','latex','fontsize',14);
ylabel('$\partial\phi_{f}/\partial\eta_{i}$','interpreter','latex','fontsize',14);
drawnow;
subplot(2,2,4);
plot(E0_vec,M22(:,s),'k');
xlabel('$\mathcal{E}_{0}$','interpreter','latex','fontsize',14);
ylabel('$\partial\phi_{f}/\partial\phi_{i}$','interpreter','latex','fontsize',14);
drawnow;
writeVideo(vidfile,getframe(gcf));

end
close(vidfile);



%%
figure; set(gcf,'color','w');
plot_ind=200;

for s=1:plot_ind:(N_max-N1_max)-plot_ind

eta_LM=Eta(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec+s-1));
phi_LM=Phi(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec+s-1));

eta_LM_next=Eta(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec+s-1+plot_ind));
phi_LM_next=Phi(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec+s-1+plot_ind));

plot(phi_LM,eta_LM,'color',[0.7 0.7 0.7]);
hold on;
plot(phi_LM_next,eta_LM_next,'k');
xlabel('$\phi$','interpreter','latex','fontsize',20);
ylabel('$\eta$','interpreter','latex','fontsize',20)
title(['$t=$',num2str(round(0.25+s/N,2)),'$T$'],'interpreter','latex','fontsize',20);
axis([1.5 3 -1.1 0.4]);
drawnow; hold on;

end


%%
load Bloch_RF_Levitt_field_ih_histogram_eta_A_190325.mat;
load Bloch_RF_Levitt_field_ih_histogram_eta_B_190325.mat;

time_unit=round((N_vec-N1_vec)/100);
etaAi=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec));
etaBi=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec));

% figure; set(gcf,'color','w');
% vidfile=VideoWriter('Bloch_RF_Levitt_field_ih_histogram_eta_eta_190325.mp4','MPEG-4');
% vidfile.FrameRate=10;
% open(vidfile);

for s=1:101
etaAf=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec+(s-1).*time_unit));
etaBf=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec+(s-1).*time_unit));
Q=(etaBf-etaAf)./(etaBi-etaAi);
range(s)=max(Q)-min(Q);
% histogram(abs(Q),20,'facecolor','b');
% xlabel('$|\langle\partial\eta_{f}/\partial\eta_{i}\rangle|$','interpreter','latex','fontsize',20);
% ylabel('$N$','interpreter','latex','fontsize',20);
% title(['$t=$',num2str(0.25+(s-1).*0.0075),'$T$'],'interpreter','latex','fontsize',20);
% drawnow;
% writeVideo(vidfile,getframe(gcf));
end
% close(vidfile);


figure; set(gcf,'color','w');
etaAi=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaBi=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaAf=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaBf=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
Q1a=(etaBf-etaAf)./(etaBi-etaAi);

subplot(2,4,1);
histogram(abs(Q1a),20,'facecolor','b');
xlabel('$|\langle\partial\eta_{f}/\partial\eta_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-200 6400]);
title('$t_{f}=T/4$','interpreter','latex','fontsize',20);
drawnow;

etaAi=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaBi=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaAf=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',round(3*N_vec/8)));
etaBf=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',round(3*N_vec/8)));
Q2a=(etaBf-etaAf)./(etaBi-etaAi);

subplot(2,4,2);
histogram(abs(Q2a),20,'facecolor','b');
xlabel('$|\langle\partial\eta_{f}/\partial\eta_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-200 6400]);
title('$t_{f}=3T/8$','interpreter','latex','fontsize',20);
drawnow;

etaAi=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaBi=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaAf=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/2));
etaBf=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/2));
Q3=(etaBf-etaAf)./(etaBi-etaAi);

subplot(2,4,3);
histogram(abs(Q3),20,'facecolor','b');
xlabel('$|\langle\partial\eta_{f}/\partial\eta_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-200 6400]);
title('$t_{f}=T/2$','interpreter','latex','fontsize',20);
drawnow;

etaAi=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaBi=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaAf=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',round(5*N_vec/8)));
etaBf=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',round(5*N_vec/8)));
Q4=(etaBf-etaAf)./(etaBi-etaAi);

subplot(2,4,4);
histogram(abs(Q4),20,'facecolor','b');
xlabel('$|\langle\partial\eta_{f}/\partial\eta_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-200 6400]);
title('$t_{f}=5T/8$','interpreter','latex','fontsize',20);
drawnow;

etaAi=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaBi=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaAf=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',3*N_vec/4));
etaBf=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',3*N_vec/4));
Q5=(etaBf-etaAf)./(etaBi-etaAi);

subplot(2,4,5);
histogram(abs(Q5),20,'facecolor','b');
xlabel('$|\langle\partial\eta_{f}/\partial\eta_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-200 6400]);
title('$t_{f}=3T/4$','interpreter','latex','fontsize',20);
drawnow;

etaAi=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaBi=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaAf=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',round(7*N_vec/8)));
etaBf=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',round(7*N_vec/8)));
Q6=(etaBf-etaAf)./(etaBi-etaAi);

subplot(2,4,6);
histogram(abs(Q6),20,'facecolor','b');
xlabel('$|\langle\partial\eta_{f}/\partial\eta_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-200 6400]);
title('$t_{f}=7T/8$','interpreter','latex','fontsize',20);
drawnow;

etaAi=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaBi=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec/4));
etaAf=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec));
etaBf=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N_vec));
Q7=(etaBf-etaAf)./(etaBi-etaAi);

subplot(2,4,7);
histogram(abs(Q7),20,'facecolor','b');
xlabel('$|\langle\partial\eta_{f}/\partial\eta_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-200 6400]);
title('$t_{f}=T$','interpreter','latex','fontsize',20);
drawnow;

subplot(2,4,8)
plot(linspace(0.25,1,length(range)),range,'.-k');
xlabel('$t_{f}/T$','interpreter','latex','fontsize',20);
ylabel('$h_{\eta}(t_{f})$','interpreter','latex','fontsize',20);

%%
load Bloch_RF_Levitt_field_ih_histogram_phi_A_190325.mat;
load Bloch_RF_Levitt_field_ih_histogram_phi_B_190325.mat;

time_unit=round((N_vec-N1_vec)/100);
phiAi=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec));
phiBi=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec));

% figure; set(gcf,'color','w');
% vidfile=VideoWriter('Bloch_RF_Levitt_field_ih_histogram_phi_phi_190325.mp4','MPEG-4');
% vidfile.FrameRate=10;
% open(vidfile);

for s=1:101
phiAf=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec+(s-1).*time_unit));
phiBf=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec+(s-1).*time_unit));
Q=(phiBf-phiAf)./(phiBi-phiAi);
range(s)=max(Q)-min(Q);
% histogram(abs(Q),20,'facecolor','b');
% xlabel('$|\langle\partial\phi_{f}/\partial\phi_{i}\rangle|$','interpreter','latex','fontsize',20);
% ylabel('$N$','interpreter','latex','fontsize',20);
% title(['$t=$',num2str(0.25+(s-1).*0.0075),'$T$'],'interpreter','latex','fontsize',20);
% drawnow;
% writeVideo(vidfile,getframe(gcf));
end
% close(vidfile);

figure; set(gcf,'color','w');
phiAi=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiBi=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiAf=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiBf=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
Q1a=(phiBf-phiAf)./(phiBi-phiAi);

subplot(2,4,1);
histogram(abs(Q1a),20,'facecolor','b');
xlabel('$|\langle\partial\phi_{f}/\partial\phi_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-1E4 4E5]);
title('$t_{f}=T/4$','interpreter','latex','fontsize',20);
drawnow;

phiAi=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiBi=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiAf=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',round(3*N_vec/8)));
phiBf=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',round(3*N_vec/8)));
Q2a=(phiBf-phiAf)./(phiBi-phiAi);

subplot(2,4,2);
histogram(abs(Q2a),20,'facecolor','b');
xlabel('$|\langle\partial\phi_{f}/\partial\phi_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-1E4 4E5]);
title('$t_{f}=3T/8$','interpreter','latex','fontsize',20);
drawnow;

phiAi=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiBi=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiAf=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/2));
phiBf=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/2));
Q3=(phiBf-phiAf)./(phiBi-phiAi);

subplot(2,4,3);
histogram(abs(Q3),20,'facecolor','b');
xlabel('$|\langle\partial\phi_{f}/\partial\phi_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-1E4 4E5]);
title('$t_{f}=T/2$','interpreter','latex','fontsize',20);
drawnow;

phiAi=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiBi=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiAf=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',round(5*N_vec/8)));
phiBf=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',round(5*N_vec/8)));
Q4=(phiBf-phiAf)./(phiBi-phiAi);

subplot(2,4,4);
histogram(abs(Q4),20,'facecolor','b');
xlabel('$|\langle\partial\phi_{f}/\partial\phi_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-1E4 4E5]);
title('$t_{f}=5T/8$','interpreter','latex','fontsize',20);
drawnow;

phiAi=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiBi=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiAf=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',3*N_vec/4));
phiBf=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',3*N_vec/4));
Q5=(phiBf-phiAf)./(phiBi-phiAi);

subplot(2,4,5);
histogram(abs(Q5),20,'facecolor','b');
xlabel('$|\langle\partial\phi_{f}/\partial\phi_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-1E4 4E5]);
title('$t_{f}=3T/4$','interpreter','latex','fontsize',20);
drawnow;

phiAi=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiBi=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiAf=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',round(7*N_vec/8)));
phiBf=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',round(7*N_vec/8)));
Q6=(phiBf-phiAf)./(phiBi-phiAi);

subplot(2,4,6);
histogram(abs(Q6),20,'facecolor','b');
xlabel('$|\langle\partial\phi_{f}/\partial\phi_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-1E4 4E5]);
title('$t_{f}=7T/8$','interpreter','latex','fontsize',20);
drawnow;

phiAi=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiBi=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec/4));
phiAf=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec));
phiBf=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N_vec));
Q7=(phiBf-phiAf)./(phiBi-phiAi);

subplot(2,4,7);
histogram(abs(Q7),20,'facecolor','b');
xlabel('$|\langle\partial\phi_{f}/\partial\phi_{i}\rangle|$','interpreter','latex','fontsize',20);
ylabel('$N$','interpreter','latex','fontsize',20);
xlim([-1E4 4E5]);
title('$t_{f}=T$','interpreter','latex','fontsize',20);
drawnow;

subplot(2,4,8);
plot(linspace(0.25,1,length(range)),range,'.-k');
xlabel('$t_{f}/T$','interpreter','latex','fontsize',20);
ylabel('$h_{\phi}(t_{f})$','interpreter','latex','fontsize',20);


%%
load Bloch_RF_Levitt_field_ih_histogram_eta_A_190325.mat;
load Bloch_RF_Levitt_field_ih_histogram_eta_B_190325.mat;
load Bloch_RF_Levitt_field_ih_histogram_phi_A_190325.mat;
load Bloch_RF_Levitt_field_ih_histogram_phi_B_190325.mat;

time_unit=round((N_vec-N1_vec)/100);

etaAi=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec));
etaBi=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec));
etacAi=EtacA(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec));
etacBi=EtacB(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec));

phiAi=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec));
phiBi=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec));
phicAi=PhicA(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec));
phicBi=PhicB(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec));
%%
for s=1:101
etaAf=EtaA(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec+(s-1).*time_unit));
etaBf=EtaB(sub2ind(size(Eta),(1:length(E0_vec)).',N1_vec+(s-1).*time_unit));

phiAf=PhiA(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec+(s-1).*time_unit));
phiBf=PhiB(sub2ind(size(Phi),(1:length(E0_vec)).',N1_vec+(s-1).*time_unit));

Q11=(((etaBf-etaAf)./(etaBi-etaAi)));
Q12=(((etaBf-etaAf)./(phicBi-phicAi)));
Q21=(((phiBf-phiAf)./(etacBi-etacAi)));
Q22=(((phiBf-phiAf)./(phiBi-phiAi)));

DetQ(:,s)=Q11.*Q22-Q12.*Q21;
drawnow;
end

%%
S=300;
time_unit=round((N_vec-N1_vec)./S);
width=zeros(S,1);

for s=1:S
x_initial=X(sub2ind(size(X),(1:length(E0_vec)).',N1_vec+s.*time_unit));
y_initial=Y(sub2ind(size(Y),(1:length(E0_vec)).',N1_vec+s.*time_unit));
z_initial=Z(sub2ind(size(Z),(1:length(E0_vec)).',N1_vec+s.*time_unit));

width(s)=sqrt((1./length(E0_vec)).*sum((x_initial-mean(x_initial)).^2+(y_initial-mean(y_initial)).^2+(z_initial-mean(z_initial)).^2));
end

% figure;
% set(gcf,'color','w');
% plot(linspace(0.25,1,length(width)),width,'.-k'); xlim([0.25 1]); ylim([0 0.1]);
% xlabel('$t/T$','interpreter','latex','fontsize',20);
% ylabel('$\sigma(t)$','interpreter','latex','fontsize',20);

width_after=width;
clear width;

S=S/3;
width=zeros(S,1);
time_unit=round(N1_vec/S);

for s=1:S
x_initial=X(sub2ind(size(X),(1:length(E0_vec)).',N1_vec-(s-1).*time_unit));
y_initial=Y(sub2ind(size(Y),(1:length(E0_vec)).',N1_vec-(s-1).*time_unit));
z_initial=Z(sub2ind(size(Z),(1:length(E0_vec)).',N1_vec-(s-1).*time_unit));

width(s)=sqrt((1./length(E0_vec)).*sum((x_initial-mean(x_initial)).^2+(y_initial-mean(y_initial)).^2+(z_initial-mean(z_initial)).^2));
end

width_before=flip(width);

width_tot=[width_before;width_after];

figure;
set(gcf,'color','w');
plot(linspace(0,1,length(width_tot)),width_tot,'.-k'); xlim([0 1]); ylim([0 0.1]);
xlabel('$t/T$','interpreter','latex','fontsize',20);
ylabel('$\sigma(t)$','interpreter','latex','fontsize',20);

%%
Q_end=N0;
y_initial=Y(sub2ind(size(Y),(1:length(E0_vec)).',N1_vec));
z_initial=Z(sub2ind(size(Z),(1:length(E0_vec)).',N1_vec));

Q1a=zeros(length(E0_vec),length(E0_vec),Q_end/4);
Q1b=zeros(length(E0_vec),length(E0_vec),Q_end/2);
Q1c=zeros(length(E0_vec),length(E0_vec),Q_end/4);

Q2a=zeros(1,Q_end/4);
Q2b=zeros(1,Q_end/2);
Q2c=zeros(1,Q_end/4);

for i=1:length(E0_vec)
for k=1:length(E0_vec)
    [i,k]

for f=1:(Q_end/4)
Q1a(i,k,f)=E0_vec(k).*cos((E0_vec(i)+E0_vec(k)).*(f-1)./(2*pi));
end

for f=1:(Q_end/2)
Q1b(i,k,f)=E0_vec(k).*z_initial(i).*z_initial(k).*sin((E0_vec(i)+E0_vec(k)).*(f-251)./(2*pi));
end

for f=1:(Q_end/4)
Q1c(i,k,f)=E0_vec(k).*y_initial(i).*y_initial(k).*sin((E0_vec(i)-E0_vec(k)).*(f-751)./(2*pi));
end

end
end

for f=1:(Q_end/4)
Q2a(f)=-2./((length(E0_vec)).^3).*sum(sum(Q1a(:,:,f)));
Q2c(f)=-2./((length(E0_vec)).^3).*sum(sum(Q1c(:,:,f)));
end

for f=1:(Q_end/2)
Q2b(f)=-2./((length(E0_vec)).^3).*sum(sum(Q1b(:,:,f)));
end

Q2_tot=[Q2a,Q2b,Q2c];
% plot(linspace(0,1,Q_end),Q2_tot,'.b');

Q3a=sqrt(dt.*cumsum(abs(Q2a)));
Q3b=sqrt(Q3a(end)+dt.*cumsum(abs(Q2b))); qab=Q3b(1)-Q3a(end); Q3b_new=Q3b-qab;
Q3c=sqrt(Q3b(end)+dt.*cumsum(abs(Q2c))); qbc=Q3c(1)-Q3b_new(end); Q3c_new=Q3c-qbc;

Q3_tot=[Q3a,Q3b_new,Q3c_new];
plot(linspace(0,1,Q_end),Q3_tot,'.b');
