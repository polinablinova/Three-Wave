function [a,b,f,xi,dxi,t,dt] = jihoon_no_extra_term_in_acou(spow,g,a0,eps,sigma,tmax,dt,ximax,dxi)

%%initializing block
xi=0:dxi:ximax;nxi=length(xi);


t=0:dt:tmax;nt=length(t);
tbeg=2;
tmid=floor(nt/2);
tend=nt;

a=zeros(nxi,nt);
b=zeros(nxi,nt);
b2=zeros(nxi,1);
f=zeros(nxi,nt);

a(1,:)=a0;


b(1:end,1)=eps*exp(-abs(xi-2*ximax/3).^(2*spow)/sigma^(2*spow));

%%loop
for j=1:nxi-1
        bb=b(j,1);db=b(j+1,1)-bb;%db=0;
        aa=a(j,1);
        ka1=-dxi*bb*g*aa*bb;
        ka2=-dxi*(bb+db/2)*(g*(bb+db/2)*(aa+ka1/2));
        ka3=-dxi*(bb+db/2)*(g*(bb+db/2)*(aa+ka2/2));
        ka4=-dxi*(bb+db)*(g*(bb+db)*(aa+ka3));

        a(j+1,1)=a(j,1)+(ka1+2*ka2+2*ka3+ka4)/6;
        %a(j+1,1)=a(j,1)+(ka1);

end



for i=1:nt-1
    for j=1:nxi-1
    %push a and f
        bb=b(j,i);db=b(j+1,i)-bb;aa=a(j,i);ff=g*aa*bb;
        ka1=-dxi*bb*g*aa*bb;
        ka2=-dxi*(bb+db/2)*(g*(bb+db/2)*(aa+ka1/2));
        ka3=-dxi*(bb+db/2)*(g*(bb+db/2)*(aa+ka2/2));
        ka4=-dxi*(bb+db)*(g*(bb+db)*(aa+ka3));

         a(j+1,i+1)=a(j,i+1)+(ka1+2*ka2+2*ka3+ka4)/6;
         if(a(j+1,i+1)<0)
            a(j+1,i+1)=0;
         end
         % a(j+1,i+1)=a(j,i+1)+ka1;
    end
    %push b
    aa=a(:,i);da=a(:,i+1)-aa;da=0;bb=b(:,i);
    kb1=dt*aa.*(g*aa.*bb);
    b(:,i+1)=b(:,i)+(kb1);
end

f=g*a.*b;
end
