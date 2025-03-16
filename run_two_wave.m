clear all
close all
clc

% RUNS jihoon_no_extra_term_in_acou.m

% Kr: gb=0.1124cm/GW, seed duration=1ns, pump duration=250ns, 
% Brillouin frequency=1.5057GHz, \tilde(rho_0)=10^-4, temperature=273K,

spow=1;
kinv=1.6573;
g=0.39605;a0=1;eps=1/400;sigma=0.10871;tmax=80;dt=0.005;ximax=2;dxi=0.001; 

[a,b,f,xi,dxi,t,dt]=jihoon_no_extra_term_in_acou(spow,g,a0,eps,sigma,tmax,dt,ximax,dxi);

ttrans=-1/(2*g*a0^2)*log(g*eps^2*sigma)

figure('Position', [100, 100, 410, 410]);  
box on; hold on; pbaspect([1 1 1])
plot(xi,b(:,end).^2,'LineWidth',3,'Color',[0.2,0.2,0.6])
plot(xi,b(:,ceil(0.8*length(t))).^2,'LineWidth',3,'LineStyle','-','Color',[0.7,0.2,0.2])
plot(xi,b(:,ceil(0.6*length(t))).^2,'LineWidth',3,'LineStyle','-','Color',[0.2,0.6,0.2])
plot(xi,b(:,ceil(0.4*length(t))).^2,'LineWidth',3,'LineStyle','-','Color',[0.5,0.3,0.2])


set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\xi=(\omega t-kz)/2$','Interpreter','latex')
ylabel('$E^2/E_{L,0}^2$','Interpreter','latex')
legend(strcat('$b^2(\xi,\tau=$',string(round(tmax)),'$)$'), ...
 strcat('$b^2(\xi,\tau=$',string(round(0.8*tmax)),'$)$'), ...
 strcat('$b^2(\xi,\tau=$',string(round(0.6*tmax)),'$)$'), ...
 strcat('$b^2(\xi,\tau=$',string(round(0.4*tmax)),'$)$'), ...
 'Interpreter','latex')

figure('Position', [100, 100, 410, 410]);  
box on; hold on; pbaspect([1 1 1])
plot(xi,a(:,end).^2,'LineWidth',3,'Color',[0.2,0.2,0.6])
plot(xi,a(:,ceil(0.8*length(t))).^2,'LineWidth',3,'LineStyle','-','Color',[0.7,0.2,0.2])
plot(xi,a(:,ceil(0.6*length(t))).^2,'LineWidth',3,'LineStyle','-','Color',[0.2,0.6,0.2])
plot(xi,a(:,ceil(0.4*length(t))).^2,'LineWidth',3,'LineStyle','-','Color',[0.5,0.3,0.2])
set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\xi=(\omega t-kz)/2$','Interpreter','latex')
ylabel('$E^2/E_{L,0}^2$','Interpreter','latex')
legend(strcat('$a^2(\xi,\tau=$',string(round(tmax)),'$)$'), ...
 strcat('$a^2(\xi,\tau=$',string(round(0.8*tmax)),'$)$'), ...
 strcat('$a^2(\xi,\tau=$',string(round(0.6*tmax)),'$)$'), ...
 strcat('$a^2(\xi,\tau=$',string(round(0.4*tmax)),'$)$'), ...
 'Interpreter','latex')

% fluence 
fluence=zeros(1,length(t));
for i=1:length(fluence)
 fluence(i)=sum(b(:,i).^2)*dxi;
end 

fl=a0^2*(t-ttrans);
fl(1:ceil(ttrans/dt))=0;

figure('Position', [100, 100, 410, 410]);  
box on; hold on; grid on; grid minor; pbaspect([1 1 1])
plot(t,fluence,'LineWidth',4,'LineStyle','-', 'Color','blue')
plot(t,fl,'LineWidth',6,'LineStyle',':','Color',"#D95319");

set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\tau=kz$','Interpreter','latex')
ylabel('$U_b$','Interpreter','latex')
legend(strcat('$U_{b}, \mathrm{simulation}$'), ...
     strcat('$a_0^2\tau$'), ...
 'Interpreter','latex')

% peak intensity
[peaks, ~] = max(b, [], 1);
inten=2*spow*a0^(2-1/spow)*g^(-1/(2*spow)).*(t-ttrans).^(1-1/(2*spow))/(sigma);
inten(1:ceil(ttrans/dt))=0;

figure('Position', [100, 100, 410, 410]); 
box on; hold on; grid on; grid minor; pbaspect([1 1 1])
plot(t,peaks.^2,'LineWidth',4,'LineStyle','-', 'Color','blue')
plot(t,inten,'LineWidth',5,'LineStyle',':','Color',"#D95319");

set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\tau=kz$','Interpreter','latex')
ylabel('$E^2/E_{L,0}^2$','Interpreter','latex')
legend(strcat('$b^2_{\mathrm{max}}, \mathrm{simulation}$'), ...
    strcat('$(2a_0/\sigma)\sqrt{\tau/g}$'), ...
 'Interpreter','latex')


% peak coordinate
crd=zeros(1,length(t));
for i=1:length(crd)
[m,pk_ind]=max(b(:,i));crd(i)=pk_ind;
end

pk=sigma*(a0*sqrt(g*(t-ttrans))).^(1/spow);
pk(1:ceil(ttrans/dt))=0;

figure('Position', [100, 100, 410, 410]);  
box on; hold on; grid on; grid minor; pbaspect([1 1 1])
plot(t,(length(xi)-crd)*dxi-(length(xi)-crd(1))*dxi, ...
    'LineWidth',4,'LineStyle','-', 'Color','blue')
plot(t,pk,'LineWidth',5,'LineStyle',':','Color',"#D95319");

set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\tau=kz$','Interpreter','latex')
ylabel('$\xi=(\omega t-kz)/2$','Interpreter','latex')
legend( strcat('$\xi_{S}, \mathrm{simulation}$'), ...
    strcat('$\sigma (a_0\sqrt{g\tau})^{1/n}$'), ...
 'Interpreter','latex')


% extraction efficiency
eff=zeros(1,length(t));
for i=1:length(t)
 eff(i)=(sum(b(:,i).^2))*dxi;
end 

eta=(t-ttrans)./t;
eta(1:ceil(ttrans/dt))=0;

figure('Position', [100, 100, 410, 410]);  
box on; hold on; grid on; grid minor; pbaspect([1 1 1])
plot(t,eff./(a0^2*t),'LineWidth',4,'LineStyle','-', 'Color','blue')
plot(t,eta,'LineWidth',6,'LineStyle',':','Color',"#D95319");

set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\tau=kz$','Interpreter','latex')
ylabel('$U_{b,\mathrm{out}}/U_{a,\mathrm{in}}$','Interpreter','latex')
legend(strcat('$\eta, \mathrm{simulation}$'), ...
    strcat('$\frac{\tau-\tau_{\mathrm{tr}}}{\tau}$'), ...
 'Interpreter','latex')

front=zeros(1,length(t));
back=zeros(1,length(t));

% width
for i=1:length(front)
 ind = find(b(:,i).^2 >= max(b(:,i).^2)/2);
 front(i)=ind(1);
 back(i)=ind(end);
end

figure('Position', [100, 100, 410, 410]);  
box on; hold on; grid on; grid minor; pbaspect([1 1 1])

width_sim=(back-front)*dxi;
%wdth=(sigma*a0/2)*sqrt(g*(t-ttrans));
wdth=( 1 - (1/2)^( 1/(2*spow-1) ) )*sigma*(a0*sqrt(g*(t-ttrans))).^(1/spow);
wdth(1:ceil(ttrans/dt))=0;

plot(t,width_sim/ximax,'LineWidth',4,'LineStyle','-', 'Color','blue')
plot(t,wdth/ximax,'LineWidth',5,'LineStyle',':','Color',"#D95319");

set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\tau=kz$','Interpreter','latex')
ylabel('$\Delta \tau_{S,\mathrm{out}}/\Delta \tau_{L,\mathrm{in}}$','Interpreter','latex')
%legend(strcat('$R, \mathrm{simulation}$'), ...
 %   strcat('$\frac{\sigma a_0}{2}\sqrt{g(\tau-\tau_{\mathrm{tr}})}$'), ...
 %'Interpreter','latex')
legend(strcat('$R, \mathrm{simulation}$'), ...
    strcat('$R, \mathrm{theory}$'), ...
 'Interpreter','latex')

