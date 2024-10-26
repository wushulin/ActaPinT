clc;
clear; 
%----- SWR method for advection-diffusion equation u'-nu*u_{xx}+u_x=g ----
global L dx dt nu N Tol Kmax
Tol=1e-8;
Kmax=1000;
nnu=2.^(-7:0);
nu=nnu(1);
nu=0.004;
dx=1/50;
x=(0:dx:8.2)';
Nx=length(x)-1;% number of space inverval
N=4;
M=round(1+(Nx-2)/N);
T=5;
dt=0.01;
t=0:dt:T;
Nt=T/dt; % number of time interval
u0=exp(-10*(x-L/2).^2);
Ix=speye(Nx-1);
r1=nu*dt/dx^2;
r2=1*dt/dx;
e=ones(Nx-1,1);
Ac=spdiags([-e e], [-1,1], Nx-1, Nx-1);
Ad=spdiags([-e 2*e -e], -1:1, Nx-1, Nx-1);
% Ac(1,Nx-1)=Ac(2,1);
% Ac(Nx-1,1)=Ac(1,2);
% Ad(1,Nx-1)=Ad(2,1);
% Ad(Nx-1,1)=Ad(1,2);
A=r1*Ad+r2*Ac;
invA=Ix/(Ix+A);
U_ref=zeros(Nx+1,Nt+1);
U_ref(:,1)=u0;
for n=1:Nt
    U_ref(2:Nx,n+1)=invA*U_ref(2:Nx,n);
end
subdomain_x=zeros(M+2,N);
for s=1:N
    x_start=1+(s-1)*M-(s-1);
    x_end=x_start+M+1;
    subdomain_x(:,s)=x(x_start:x_end);
end
Ldx=2*dx;
y0=Ldx/nu;
C2=dt/dx;
yc=1.618386576;
if y0>=yc
    tp0=fsolve(@(p) y0-p*sqrt(p/(4+p)),2);
else
    tp0=fsolve(@(p) fun_R0(y0,p,y0)-fun_R0(ty(p,y0),p,y0),1.2);
end
tp=abs(tp0);
p=tp*nu/Ldx;
q=1/(1+dx*p);
q=0;
e=ones(M,1);
IM=speye(M);
A1_s1=spdiags([-e 2*e -e], -1:1, M, M);
A2_s1=spdiags([-e e], [-1,1], M, M);
A1_s1(M,M)=2-q;
A2_s1(M,M)=q;
A_s1=IM+r1*A1_s1+r2*A2_s1;
B_r=zeros(M,M);
B_r(M,1)=-(r1-r2)*q;
B_r(M,2)=r1-r2;

A1_sn=spdiags([-e 2*e -e], -1:1, M, M);
A2_sn=spdiags([-e e], [-1,1], M, M);
A1_sn(1,1)=2-q;
A1_sn(M,M)=2-q;
A2_sn(1,1)=q;
A2_sn(M,M)=q;
A_sn=IM+r1*A1_sn+r2*A2_sn;
B_l=zeros(M,M);
B_l(1,M)=-(r1+r2)*q;
B_l(1,M-1)=r1+r2;

A1_sN=spdiags([-e 2*e -e], -1:1, M, M);
A2_sN=spdiags([-e e], [-1,1], M, M);
A1_sN(1,1)=2-q;
A2_sN(1,1)=q;
A_sN=IM+r1*A1_sN+r2*A2_sN;
Q=zeros(N*M,N*M);
W=zeros(N*M,N*M);
for s=1:N
    if s==1
        Q((s-1)*M+1:s*M,(s-1)*M+1:s*M)=A_s1;
        W((s-1)*M+1:s*M,s*M+1:(s+1)*M)=B_r;
    elseif s==N
        Q((s-1)*M+1:s*M,(s-1)*M+1:s*M)=A_sN;
        W((s-1)*M+1:s*M,(s-2)*M+1:(s-1)*M)=B_l;
    else
        Q((s-1)*M+1:s*M,(s-1)*M+1:s*M)=A_sn;
        W((s-1)*M+1:s*M,(s-2)*M+1:(s-1)*M)=B_l;
        W((s-1)*M+1:s*M,s*M+1:(s+1)*M)=B_r;
    end
end
U0=zeros(N*M,1);
for s=1:N
    U0((s-1)*M+1:s*M)=get_u0(subdomain_x(2:M+1,s));
end

Uk=zeros(N*M,Nt+1);
Uk(:,1)=U0;
Uk1=Uk;
IU=speye(N*M);
invQ=IU/Q;
invQW=invQ*W;
Err=zeros(1,Kmax);
for k=1:Kmax
    for n=1:Nt
        Uk1(:,n+1)=invQ*Uk1(:,n)+invQW*Uk(:,n+1);
    end
    Err(k)=max(max(abs(Uk1-Uk)));
    Uk=Uk1;
    fprintf('%d-th SWR iteration: error=%2.15f\n',k,Err(k));
    if Err(k)<=Tol
        break;
    end
end
%mesh(t,x,U_ref);shg

function val=get_u0(x)
global L
val=exp(-10*(x-L/2).^2);
end

function val=fun_R0(y,tp,y0)
val=-(((y-tp)^2+y^2-y0^2)/((y+tp)^2+y^2-y0^2))*exp(-y);
end

function val=ty(tp,y0)
d=tp*(-tp^3-4*tp^2+(4+2*y0^2)*tp+8*y0^2);
val=sqrt((y0^2+2*tp+sqrt(d))/2);
end