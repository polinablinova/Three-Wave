function [a,b,f,xi,dxi,t,dt] = ringing_back(G,spow,g,a0,eps,sigma,tmax,dt,ximax,dxi)

% BACKWARD EULER IMPLEMENTATION OF THE THREE-WAVE EQUATION

%%initializing block
xi=0:dxi:ximax;nxi=length(xi);


t=0:dt:tmax;nt=length(t);
tbeg=2;
tmid=floor(nt/2);
tend=nt;

a=zeros(nxi,nt);
b=zeros(nxi,nt);
f=zeros(nxi,nt);

a(1,:)=a0;
f(1,:)=0;


b(1:end,1)=eps*exp(-abs(xi-2*ximax/3).^(2*spow)/sigma^(2*spow));

%%loop
for j=1:nxi-1
      bb=b(j,1);db=b(j+1,1)-bb;
      ff=f(j,1);aa=a(j,1);
      ka1=-bb*ff;kf1=G*(aa*bb-g*ff);
        

      a(j+1,1)=a(j,1)+(ka1)*dxi;
      %f(j+1,1)=f(j,1)+(kf1)*dxi;
      f(j+1,1)=(f(j,1)+G*dxi*g*a(j+1,1)*b(j+1,1))/(1+G*dxi);

end

    ff=f(:,1);aa=a(:,1);
    kb1=dt*aa.*ff;
    b(:,2)=b(:,1)+(kb1);


for i=1:nt-1
    for j=1:nxi-1
    %push a and f
        bb=b(j,i+1);db=b(j+1,i+1)-bb;aa=a(j,i+1);ff=f(j,i+1);
    if G>1       
        ka1 = -dxi * bb * ff;
        %kf1 = G * dxi * (g * aa * bb - ff);
        a(j+1,i+1) = aa + ka1; 
        %f(j+1,i+1) = f(j,i) + kf1;
        f(j+1,i+1) = (ff + G*dxi*g*a(j+1,i+1)*b(j+1,i+1))/(1+G*dxi);
    else
        ka1=-dxi*bb*ff;
        kf1=G*dxi*(g*aa*bb-ff);

        ka2=-dxi*(bb+db/2)*(ff+kf1/2);
        kf2=G*dxi*(g*(aa+ka1/2)*(bb+db/2)-(ff+kf1/2));

        ka3=-dxi*(bb+db/2)*(ff+kf2/2);
        kf3=G*dxi*(g*(aa+ka2/2)*(bb+db/2)-(ff+kf2/2));

        ka4=-dxi*(bb+db)*(ff+kf3);
        kf4=G*dxi*(g*(aa+ka3)*(bb+db)-(ff+kf3));

        a(j+1,i+1)=a(j,i+1)+(ka1+2*ka2+2*ka3+ka4)/6;
        f(j+1,i+1)=f(j,i+1)+(kf1+2*kf2+2*kf3+kf4)/6;
    end 

    end
    %push b
    ff=f(:,i+1);aa=a(:,i+1);
    %df=f(:,i+1)-ff;da=a(:,i+1)-aa;
    %aa=a(:,i);ff=f(:,i);da=a(:,i+1)-aa;df=f(:,i+1)-ff;
    %kb1=dt*aa.*ff;kb2=dt*(aa+da/2).*(ff+df/2);
    %kb3=dt*(aa+da/2).*(ff+df/2);kb4=dt*(aa+da).*(ff+df);
    %b(:,i+2)=b(:,i+1)+(kb1+2*kb2+2*kb3+kb4)/6;
    kb1=dt*aa.*ff;
    b(:,i+2)=b(:,i+1)+(kb1);

end

