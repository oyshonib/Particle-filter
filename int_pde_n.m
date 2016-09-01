%PF_param;
function [u,v,theta] = int_pde_n(M,x,nsteps,dt,UgM,VgM,N,Hs,param)

linf = param(1);
alpha = param(2);
theta_inf = param(3);

k=0.4;
g=9.8;
beta=15;
f=10^-4.0;

xs = 0.5*(x(1:length(x)-1) + x(2:length(x)));

%set initial conditions
u=M(1,:)';
v=M(2,:)';
theta=M(3,:)';
uNew=u;
vNew=v;
thetaNew=theta;


fM=zeros(length(xs),1);
fH=zeros(length(xs),1);

for n = 1 : nsteps 
    %get Ug,Vg
     nwind=floor((n-1)*dt/3600)+1;
     prop = (n-1-3600/dt*(nwind-1))*dt/3600;  %proportion of interval to next geowind data
     
     Ug=prop*UgM(nwind+1)+(1-prop)*UgM(nwind);
     Vg=prop*VgM(nwind+1)+(1-prop)*VgM(nwind);
    
% Calculate l, Ri, fM,H, kM,H - tested
        v2u2 = ((v(2:N+2)-v(1:N+1)).^2+(u(2:N+2)-u(1:N+1)).^2)./(x(2:N+2) - x(1:N+1)).^2; %(u_x^2 + v_x^2)/(x_{j+1}-x_{j})^2
        l = 1./(1./(k*xs(1:N+1)) + 1/linf);
        Ri = g./(0.5*(theta(2:N+2)+theta(1:N+1))).*((theta(2:N+2) - theta(1:N+1))./(x(2:N+2) - x(1:N+1)))./v2u2;

    for j=1:N+1
        if (Ri(j) < 0)
            fM(j) = (1-beta*Ri(j))^0.5;
            fH(j) = (1-beta*Ri(j))^0.75;
        elseif (alpha*Ri(j) <= 1) && (alpha*Ri(j) >= 0)
            fM(j) = (1- alpha*Ri(j))^2;
            fH(j) = fM(j);
        else
            fM(j) = 0;
            fH(j) = fM(j);
        end
    end
    midk = l.^2.*sqrt(v2u2);
    kM = midk.*fM(1:N+1);
    kH = midk.*fH(1:N+1);
        
    
    %calculate new u,v and theta
    %NB: main nodes are from 2 to N+2, excludes ghost nodes
    for j=3:N+1   
        %update u
        uNew(j)=u(j)+dt*( f*(v(j)-Vg)+ (kM(j)*(((u(j+1)-u(j))/(x(j+1)-x(j)))) - kM(j-1)*(((u(j)-u(j-1))/(x(j)-x(j-1)))))/...
       (xs(j)-xs(j-1)));
        %update v
        vNew(j)=v(j)+dt*( -f*(u(j)-Ug)+ (kM(j)*(((v(j+1)-v(j))/(x(j+1)-x(j)))) - kM(j-1)*(((v(j)-v(j-1))/(x(j)-x(j-1)))))/...
       (xs(j)-xs(j-1)));
        %update theta
        thetaNew(j)=theta(j)+dt*( kH(j)*(((theta(j+1)-theta(j))/(x(j+1)-x(j)))) - kH(j-1)*(((theta(j)-theta(j-1))/(x(j)-x(j-1)))))/...
       (xs(j)-xs(j-1));
    end
    
    %2. Correct boundary conditions
    %u and v
    uNew(2)=0.0;   uNew(N+2)=Ug;
    vNew(2)=0.0;   vNew(N+2)=Vg;

    thetaNew(N+2)=theta_inf; 

     nwind2 = 1+floor((n-1)*dt/(600));
     prop2 = (n-1-600/dt*(nwind2-1))*dt/600;
     
     theta(1) = theta(3)+2*(prop2*Hs(1+nwind2) + (1-prop2)*Hs(nwind2))*(x(3)-x(1))/(kM(2)+kM(1));
    thetaNew(2)=theta(2)+dt*( kH(2)*(((theta(3)-theta(2))/(x(3)-x(2)))) - kH(1)*(((theta(2)-theta(1))/(x(2)-x(1)))))/...
    (xs(2)-xs(1));

    %set u,v and theta for next time loop
    u=uNew;    
    v=vNew;
    theta=thetaNew;
 
%      plot(x(2:N+2),u(2:N+2));
%      plot(x(2:N+2),v(2:N+2)); 
%      plot(theta(2:N+2), x(2:N+2),'bo-','markerfacecolor','b');
%      shg
%      pause(0.1);
%      theta(2)
%      pause(0.01)
end