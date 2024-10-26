clc;
clear;
%----Code of ParaExp for Burgers equation u'-nu*u_{xx}+0.5*(u^2)_x=g----
global dx Dt J T nu NT_Tol Tol NTmax c Kmax
flag_BC=2;
nu=2e-2;
%nu=1;
T=2;
Dt=0.01;
J=20;
dt=Dt/J;

t=0:Dt:T;
Nt=length(t);
dx=1/100;
Tol0=max(dt,dx^2);
Tol=1e-10;
NT_Tol=1e-13;
NTmax=100;
Kmax=12;
x=(0:dx:1)';
c=1;
NT_It=zeros(1,Nt);
 
Nx=length(x);
Ix=speye(Nx);
e=ones(Nx,1);
A1 = 0.5*spdiags([-e e], [-1,1], Nx, Nx)/(2*dx);
A2= nu*spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx^2;
A2(1,Nx)=A2(2,1);
A2(Nx,1)=A2(1,2);
A1(1,Nx)=A1(2,1);
A1(Nx,1)=A1(1,2);
u_ref=zeros(Nx,Nt);
for j=1:Nx
    if abs(x(j)-0.5)<=0.2
        u_ref(j,1)=1;
    end
end
for n=1:Nt-1
    [val,N_It]=ProF_IE(u_ref(:,n),A1,A2);
    NT_It(n+1)=N_It;
    u_ref(:,n+1)=val;
end
Uk_IE=zeros(Nx,Nt);
Uk_EXP=zeros(Nx,Nt);
Uk1_IE=zeros(Nx,Nt);
Uk1_EXP=zeros(Nx,Nt);
Gk_IE=zeros(Nx,Nt);
Gk_IE(:,1)=u_ref(:,1);
for n=1:Nt-1
    Gk_IE(:,n+1)=ProG_EXP(Gk_IE(:,n),A2);
end
Gk_EXP=Gk_IE;
Uk_IE(:,:)=Gk_IE;
Uk_EXP(:,:)=Gk_IE;
Uk1_IE(:,1)=u_ref(:,1);
Uk1_EXP(:,1)=u_ref(:,1);
U_IE=zeros(Nx,Nt,Kmax);
U_EXP=zeros(Nx,Nt,Kmax);
k=1;
U_IE(:,:,k)=Gk_IE;
U_EXP(:,:,k)=Gk_EXP;
Err_IE=zeros(1,Kmax);
Err_EXP=zeros(1,Kmax);
k=1;
Err_IE(k)=max(max(abs(u_ref-Uk_IE(:,:,k))));
Err_EXP(k)=Err_IE(k);
fprintf('The error at %d-th iteration is (Err_IE, Err_EXP)=(%2.15f, %2.15f)\n',k,Err_IE(k),Err_EXP(k));

for k=1:Kmax-1
    for n=1:Nt-1
        [temp,~]=ProG_IE(Uk1_IE(:,n),A1,A2);
        Uk1_IE(:,n+1)=temp+ProF_IE(Uk_IE(:,n),A1,A2)-Gk_IE(:,n+1);
        Gk_IE(:,n+1)=temp;
    end
    Uk_IE=Uk1_IE;
    Err_IE(k+1)=max(max(abs(u_ref-Uk_IE)));
    for n=1:Nt-1
        temp=ProG_EXP(Uk1_EXP(:,n),A2);
        Uk1_EXP(:,n+1)=temp+ProF_IE(Uk_EXP(:,n),A1,A2)-Gk_EXP(:,n+1);
        Gk_EXP(:,n+1)=temp;
    end
    Uk_EXP=Uk1_EXP;
    U_IE(:,:,k+1)=Uk_IE;
    U_EXP(:,:,k+1)=Uk_EXP;
    Err_EXP(k+1)=max(max(abs(u_ref-Uk_EXP)));
    fprintf('The error at %d-th iteration is (Err_IE, Err_EXP)=(%2.15f, %2.15f)\n',k+1,Err_IE(1+k),Err_EXP(k+1));
    if max(Err_EXP(k+1),Err_IE(k+1))<=Tol
        break;
    end
end
figure(1);
semilogy(0:Kmax-1,Err_EXP,'b-.o',0:Kmax-1,Err_IE,'r-.+',...
    0:Kmax-1,Tol0*ones(Kmax),'k--','linewidth',1,'markersize',11);shg
 set(gca,'fontsize',16);
xlabel('Iteration Number','fontsize',17);
ylabel('Error','fontsize',17);
title(['$\nu=$',num2str(nu)],'interpreter','latex','fontsize',20);
leg=legend('Iterative ParaEXP','Standard Parareal');
set(leg,'fontsize',15);
ylim([Tol,2*max(max(Err_IE),max(Err_EXP))]);
set(gca,'ytick',10.^(-10:0));
set(gca,'xtick',0:Kmax-1);
xlim([0,10])
 ylim([1e-8,0.5])
% figure(2);
% plot(t,U_IE(round(0.5/dx),:,1),'r-.',t,U_IE(round(0.5/dx),:,3),'b-.',t,U_IE(round(0.5/dx),:,5),'m-.',t,u_ref(round(0.5/dx),:),'k','linewidth',1);shg
% plot(x,U_IE(:,Nt,1),'r-.',x,U_IE(:,Nt,3),'b-.',x,U_IE(:,Nt,5),'m-.',x,u_ref(:,Nt),'k','linewidth',1);shg
% figure(3);
% plot(x,U_EXP(:,Nt,1),'r-.',x,U_EXP(:,Nt,3),'b-.',x,U_EXP(:,Nt,5),'m-.',x,u_ref(:,Nt),'k','linewidth',1);shg
% plot(t,U_EXP(round(0.5/dx),:,1),'r-.',t,U_EXP(round(0.5/dx),:,3),'b-.',t,U_EXP(round(0.5/dx),:,5),'m-.',t,u_ref(round(0.5/dx),:),'k','linewidth',1);shg


function val=ProG_EXP(U,A2)
global Dt  
val=expm(-Dt*A2)*U;
end


function [val,N_It]=ProG_IE(U,A1,A2)
global Dt c NT_Tol NTmax
Ix=eye(length(U));
z=U;
for k=1:NTmax
    res=z-U+ Dt*A2*(c*z+(1-c)*U)+Dt*A1*(c*z.^2+(1-c)*U.^2);
    if norm(res,inf)<=NT_Tol
        break;
    else
        Jac=Ix+Dt*c*A2+2*Dt*c*A1*diag(z);
        z=z-Jac\res;
    end
end
val=z;
N_It=k;
end
 
function [val,N_It]=ProF_IE(U,A1,A2)
global J Dt c NT_Tol NTmax
dt=Dt/J;
z0=U;
z=U;
N_It=0;
Ix=eye(length(U));
for j=1:J
    for k=1:NTmax
        res=z-z0+ dt*A2*(c*z+(1-c)*z0)+dt*A1*(c*z.^2+(1-c)*z0.^2);
        if norm(res,inf)<=NT_Tol
            break;
        else
            Jac=Ix+dt*c*A2+2*dt*c*A1*diag(z);
            z=z-Jac\res;
        end
    end
    z0=z;
    N_It=N_It+k;
end
val=z;
end