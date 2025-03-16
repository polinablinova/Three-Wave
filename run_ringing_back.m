clear 
close all
clc

% RUNS ringing_back.m

g=1;a0=1.5;eps=1/400;sigma=1;
tmax=40;dt=0.005;ximax=20;dxi=0.002;
spow=1;
G=200;

[a,b,f,xi,dxi,t,dt]=ringing_back(G,spow,g,a0,eps,sigma,tmax,dt,ximax,dxi);


figure('Position', [100, 100, 500, 500]); 
box on; hold on; pbaspect([1 1 1])

plot(xi,b(:,end),'LineWidth',3)
plot(xi,a(:,end),'LineWidth',3)
plot(xi,f(:,end),'LineWidth',3)



set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\xi=(\omega t-kz)/2$','Interpreter','latex')
ylabel('$E/E_{L,0}$','Interpreter','latex')
legend(strcat('$b(\xi,\tau=$',string(round(tmax)),'$)$'), ...
 strcat('$a(\xi,\tau=$',string(round(tmax)),'$)$'), ...
 strcat('$f(\xi,\tau=$',string(round(tmax)),'$)$'),...
 'Interpreter','latex')

ttrans=-1/(2*g*a0^2)*log(g*eps^2*sigma*sqrt(pi/2)/(2*log(2)))
%% peak coordinate
crd=zeros(1,length(t));
for i=1:length(crd)
[m,pk_ind]=max(b(:,i));crd(i)=pk_ind;
end

% modified two-wave
pk=sigma*(a0*sqrt(g*(t-ttrans))).^(1/spow)-g*a0^2.*(t-ttrans)/G;
pk(1:ceil(ttrans/dt))=0;

% original two-wave solution
pk2=sigma*(a0*sqrt(g*(t-ttrans))).^(1/spow);
pk2(1:ceil(ttrans/dt))=0;

figure('Position', [100, 100, 350, 350]); 
box on; hold on; grid on; grid minor; pbaspect([1 1 1])

plot(t,pk,'LineWidth',3,'Color','black','LineStyle',":")
plot(t,pk2,'LineWidth',3,'Color','red','LineStyle',"-")
plot(t,(length(xi)-crd)*dxi-(length(xi)-crd(1))*dxi,'LineWidth',4,'LineStyle','--', 'Color','blue')


set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\tau=kz$','Interpreter','latex')
ylabel('$\xi=(\omega t-kz)/2$','Interpreter','latex')
legend(strcat('Intermediate Regime'), ...
 strcat('Two-Wave'), ...
 strcat('Simulation'), ...
 'Interpreter','latex')

%% peak intensity
[peaks, ~] = max(b, [], 1);

% modified two-wave
inten=a0^2./( sigma*(a0*sqrt(g))^(1/spow).*(t-ttrans).^(1/(2*spow)-1)./(2*spow)-g*a0^2/G);
inten(1:ceil(ttrans/dt))=0;

% original two-wave
inten2=2*spow*a0^(2-1/spow)*g^(-1/(2*spow)).*(t-ttrans).^(1-1/(2*spow))/(sigma);
inten2(1:ceil(ttrans/dt))=0;

figure('Position', [100, 100, 350, 350]);  
box on; hold on; grid on; grid minor; pbaspect([1 1 1])
plot(t,inten,'LineWidth',3,'Color','black','LineStyle',":")
plot(t,inten2,'LineWidth',3,'Color','red','LineStyle',"-")
plot(t,peaks(1:end-1).^2,'LineWidth',4,'LineStyle','--', 'Color','blue')

set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\tau=kz$','Interpreter','latex')
ylabel('$E^2/E_{L,0}^2$','Interpreter','latex')
legend(strcat('Intermediate Regime'), ...
 strcat('Two-Wave'), ...
 strcat('Simulation'), ...
 'Interpreter','latex')

%% plots
figure('Position', [100, 100, 410, 410]);  
box on; hold on; pbaspect([1 1 1])
plot(xi,b(:,end),'LineWidth',3,'Color',[0.2,0.2,0.6])
plot(xi,b(:,ceil(0.8*length(t))),'LineWidth',3,'LineStyle','-','Color',[0.7,0.2,0.2])
plot(xi,b(:,ceil(0.6*length(t))),'LineWidth',3,'LineStyle','-','Color',[0.2,0.6,0.2])
plot(xi,b(:,ceil(0.4*length(t))),'LineWidth',3,'LineStyle','-','Color',[0.5,0.3,0.2])


set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\xi=(\omega t-kz)/2$','Interpreter','latex')
ylabel('$b$','Interpreter','latex')
legend(strcat('$b(\xi,\tau=$',string(round(tmax)),'$)$'), ...
 strcat('$b(\xi,\tau=$',string(round(0.8*tmax)),'$)$'), ...
 strcat('$b(\xi,\tau=$',string(round(0.6*tmax)),'$)$'), ...
 strcat('$b(\xi,\tau=$',string(round(0.4*tmax)),'$)$'), ...
 'Interpreter','latex')

figure('Position', [100, 100, 410, 410]);  % [left, bottom, width, height]
box on; hold on; pbaspect([1 1 1])
plot(xi,a(:,end),'LineWidth',3,'Color',[0.2,0.2,0.6])
plot(xi,a(:,ceil(0.8*length(t))),'LineWidth',3,'LineStyle','-','Color',[0.7,0.2,0.2])
plot(xi,a(:,ceil(0.6*length(t))),'LineWidth',3,'LineStyle','-','Color',[0.2,0.6,0.2])
plot(xi,a(:,ceil(0.4*length(t))),'LineWidth',3,'LineStyle','-','Color',[0.5,0.3,0.2])
set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\xi=(\omega t-kz)/2$','Interpreter','latex')
ylabel('$a$','Interpreter','latex')
legend(strcat('$a(\xi,\tau=$',string(round(tmax)),'$)$'), ...
 strcat('$a(\xi,\tau=$',string(round(0.8*tmax)),'$)$'), ...
 strcat('$a(\xi,\tau=$',string(round(0.6*tmax)),'$)$'), ...
 strcat('$a(\xi,\tau=$',string(round(0.4*tmax)),'$)$'), ...
 'Interpreter','latex')


figure('Position', [100, 100, 410, 410]);  % [left, bottom, width, height]
box on; hold on; pbaspect([1 1 1])
plot(xi,f(:,end),'LineWidth',3,'Color',[0.2,0.2,0.6])
plot(xi,f(:,ceil(0.8*length(t))),'LineWidth',3,'LineStyle','-','Color',[0.7,0.2,0.2])
plot(xi,f(:,ceil(0.6*length(t))),'LineWidth',3,'LineStyle','-','Color',[0.2,0.6,0.2])
plot(xi,f(:,ceil(0.4*length(t))),'LineWidth',3,'LineStyle','-','Color',[0.5,0.3,0.2])
set(gca,'fontsize', 18) 
set(gca,'linewidth',1)
xlabel('$\xi=(\omega t-kz)/2$','Interpreter','latex')
ylabel('$f$','Interpreter','latex')
legend(strcat('$f(\xi,\tau=$',string(round(tmax)),'$)$'), ...
 strcat('$f(\xi,\tau=$',string(round(0.8*tmax)),'$)$'), ...
 strcat('$f(\xi,\tau=$',string(round(0.6*tmax)),'$)$'), ...
 strcat('$f(\xi,\tau=$',string(round(0.4*tmax)),'$)$'), ...
 'Interpreter','latex')

%% check convergence to original two-wave
close all
g=2;a0=1;eps=1/400;sigma=0.5;
tmax=50;dt=0.01;ximax=15;dxi=0.002;
spow=1;
G=5000;
[a,b,f,xi,dxi,t,dt]=ringing_back(G,spow,g,a0,eps,sigma,tmax,dt,ximax,dxi);
[a1,b1,f1,xi1,dxi1,t1,dt1]=jihoon_no_extra_term_in_acou(spow,g,a0,eps,sigma,tmax,dt,ximax,dxi);

figure('Position', [100, 100, 400, 400]);  % [left, bottom, width, height]
box on; hold on; pbaspect([1 1 1])

plot(xi,b(:,end),'LineWidth',3)
plot(xi,a(:,end),'LineWidth',3)

plot(xi1,b1(:,end),'LineWidth',3,'LineStyle','--')
plot(xi1,a1(:,end),'LineWidth',3)

set(gca,'fontsize', 18) 
set(gca,'linewidth',2)
xlabel('$\xi=(\omega t-kz)/2$','Interpreter','latex')
ylabel('$E/E_{L,0}$','Interpreter','latex')
legend(strcat('$b(\xi,\tau=$',string(round(tmax)),'$)$'), ...
 strcat('$a(\xi,\tau=$',string(round(tmax)),'$)$'), ...
 strcat('$b_{\mathrm{two-wave}}(\xi,\tau=$',string(round(tmax)),'$)$'),...
  strcat('$a_{\mathrm{two-wave}}(\xi,\tau=$',string(round(tmax)),'$)$'),...
 'Interpreter','latex')
