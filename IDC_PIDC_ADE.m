clc;
clear;
%----this code is for advection-diffusion equation u'-nu*u_{xx}+u_x=g----
%----the main funciton PIDC is used for both IDC and PIDC depending the input parameter flag_PIDC=0 or 1----
global dt M Nt Kmax nu c x sp  % M is the quadrature nodes
sp=1000; 
sp=5; 
% sd: the parameter sigma used for the Gaussian parameter 
c=2/2;   % c=0.5: TR; c=1: Backward Euler            
nu=1e-3;
T=3;
Dt=1/10;
M=5;
Kmax=3;

Nt=T/Dt;
dx=1/64;
x=(0:dx:1)';
  
Nx=length(x);
Ix=speye(Nx);
e=ones(Nx,1);
A1 = spdiags([e -e], [-1,1], Nx, Nx)/(2*dx);
A2= spdiags([e -2*e e], -1:1, Nx, Nx)/dx^2;
A2(1,Nx)=A2(2,1);
A2(Nx,1)=A2(1,2);
A1(1,Nx)=A1(2,1);
A1(Nx,1)=A1(1,2);
A=nu*A2+A1;
nPoints=Nt*(M-1)+1;
t=linspace(0,T,nPoints);
dt=t(2)-t(1);
invA=Ix/(Ix-c*dt*A);
u0=0*sin(8*pi*(1-x).^2).^2; 
opts = odeset('RelTol',(1e-13),'AbsTol',(1e-13)*ones(Nx,1));
[tt,y] = ode45(@(t,u) odefun(t,u,A), t, u0, opts);
uExact=transpose(y);
%------------ PIDC-----------
flag_PIDC=1; 
Err_PIDC=zeros(Nt,Kmax);
[t,uIDC_ref]=PIDC(A,[0,T],u0,6,flag_PIDC);
for K=1:Kmax
    [t,uIDC]=PIDC(A,[0,T],u0,K,flag_PIDC);
    for n=1:Nt
        Err_PIDC(n,K)=max(max(abs((uExact(:,n*(M-1)+1-(M-1):n*(M-1)+1)-uIDC(:,n*(M-1)+1-(M-1):n*(M-1)+1))/...
            max(max(abs(uExact))))));
    end
end
%------------ IDC-----------
flag_PIDC=0;
Err_IDC=zeros(Nt,Kmax);
for K=1:Kmax
    [t,uIDC]=PIDC(A,[0,T],u0,K,flag_PIDC);
    for n=1:Nt
         Err_IDC(n,K)=max(max(abs((uExact(:,n*(M-1)+1-(M-1):n*(M-1)+1)-uIDC(:,n*(M-1)+1-(M-1):n*(M-1)+1))/...
         max(max(abs(uExact))))));
    end
end
% semilogy(1:Nt,Err_PIDC(:,1),'r-.o', 1:Nt,Err_IDC(:,1),'b-.+',...
%          1:Nt,Err_PIDC(:,3),'r-o',1:Nt,Err_IDC(:,3),'b-+',...
%          1:Nt,Err_PIDC(:,5),'r:o',1:Nt,Err_IDC(:,5),'b:+',...
%          'linewidth',1,'markersize',10);shg
figure(2);
d=1;
semilogy(1:d:Nt,Err_IDC(1:d:Nt,1),'b-.+',1:d:Nt,Err_IDC(1:d:Nt,2),'b-+',1:d:Nt,Err_IDC(1:d:Nt,3),'b:+',...
    1:d:Nt,Err_PIDC(1:d:Nt,1),'r-.o', 1:d:Nt,Err_PIDC(1:d:Nt,2),'r-o', 1:d:Nt,Err_PIDC(1:d:Nt,3),'r:o',...
         'linewidth',1,'markersize',10);shg
set(gca,'fontsize',14);
xlabel('Window Index','interpreter','latex','fontsize',19);
ylabel('Max Erorr on Each Window','interpreter','latex','fontsize',19);
% leg=legend( 'IDC $k_{\max}=0$ (initial)','IDC $k_{\max}=1$','IDC $k_{\max}=2$',...
%                   'PIDC $k_{\max}=0$ (initial)','PIDC $k_{\max}=1$','PIDC $k_{\max}=2$');
% set(leg,'interpreter','latex','fontsize',15);
title(['$\nu=',num2str(nu),'$'],'interpreter','latex','fontsize',20);
% if c==1
%     title(['backward Euler, $\nu=',num2str(nu),'$'],'interpreter','latex','fontsize',20);
% else
%     title(['trapezoidal rule, $\nu=',num2str(nu),'$'],'interpreter','latex','fontsize',20);
% end
xlim([1,Nt]);
set(gca,'xtick',[1:2:Nt]);
ylim([1e-7,1e-0])
%ylim([8e-6,25]);
set(gca,'ytick',10.^(-10:1:1));shg



function [t,u]=PIDC(A,tspan,u0,K,flag_PIDC)
global c Nt M
weights=LagrangeWeights(M); % Integral term
nPoints=Nt*(M-1)+1;
t=linspace(tspan(1),tspan(2),nPoints); % Time grid for all intervals
dt=t(2)-t(1); 
if flag_PIDC==1
    uStep=kron(ones(1,K),u0); % Initial guess
else
    uStep=u0; % Initial guess    
end
Nx=length(u0);
Ix=eye(Nx);
invA=Ix/(Ix-dt*c*A);
u=zeros(Nx,Nt,K,M);
for n=1:Nt
    u(:,n,1,1)=uStep(:,1);
    for m=1:M-1
        u(:,n,1,m+1)=u(:,n,1,1);
        %u(:,n,1,m+1)=invA*((Ix+dt*(1-c)*A)*u(:,n,1,m)+dt*c*f(t((n-1)*(M-1)+m+1))+dt*(1-c)*f(t((n-1)*(M-1)+m)));
    end
    for k=2:K % Sweep iterations
        if flag_PIDC==1
            u(:,n,k,1)=uStep(:,k);
        else
            u(:,n,k,1)=uStep(:,1);
        end
        for m=1:M-1
            integral=zeros(Nx,1);
            for mm=1:M
                %integral=integral+(M-1)*dt*weights(m, mm)*A*u(:,n,k-1,mm);
                integral=integral+(M-1)*dt*weights(m, mm)*odefun(t((n-1)*(M-1)+mm),u(:,n,k-1,mm),A);
            end
            u(:,n,k,m+1)=invA*(u(:,n,k,m)-c*dt*A*u(:,n,k-1,m+1)+(1-c)*dt*A*(u(:,n,k,m)-u(:,n,k-1,m))+integral);
        end
    end
    if flag_PIDC==1
        for k=1:K
            uStep(:,k)=u(:,n,k,M); % Initial values for next step
        end
    else
        uStep(:,1)=u(:,n,K,M);
    end
end
uAll=u0; % Combine each step
for n=1:Nt
    for m=2:M
        uAll=[uAll,u(:,n,K,m)];
    end
end
u=uAll;
end

function w=LagrangeWeights(M)
% LAGRANGEWEIGHTS Compute Lagrange weights
% w=LagrangeWeights(M) Returns analytical Lagrange weights for M uniform
% points.
weights{2}=[1/2, 1/2];
weights{3}=[
[ 5/24, 1/3, -1/24];
[-1/24, 1/3, 5/24]];
weights{4}=[
[ 9/72, 19/72, -5/72, 1/72];
[-1/72, 13/72, 13/72, -1/72];
[ 1/72, -5/72, 19/72, 9/72]];
weights{5}=[
[ 251/2880, 323/1440, -11/120, 53/1440, -19/2880];
[ -19/2880, 173/1440, 19/120, -37/1440, 11/2880];
[ 11/2880, -37/1440, 19/120, 173/1440, -19/2880];
[ -19/2880, 53/1440, -11/120, 323/1440, 251/2880]];
weights{6}=[
[ 19/288, 1427/7200, -399/3600, 241/3600, -173/7200, 3/800];
[ -3/800, 637/7200, 511/3600, -129/3600, 77/7200, -11/7200];
[ 11/7200, -31/2400, 401/3600, 401/3600, -31/2400, 11/7200];
[-11/7200, 77/7200, -129/3600, 511/3600, 637/7200, -3/800];
[ 3/800, -173/7200, 241/3600, -399/3600, 1427/7200, 19/288]];
w=weights{M};
end

function val=g(tn) % 'val' means 'value'
global   x sp
y=x;
t=tn;
val=exp(-sp*(y-0.5).^2).*10*(exp(-sp*(t-0.1).^2)+exp(-sp*(t-0.6).^2)+exp(-sp*(t-1.35).^2)+exp(-sp*(t-1.85).^2));
end


function  val = odefun(t,u,A)
val=A*u+g(t);
end