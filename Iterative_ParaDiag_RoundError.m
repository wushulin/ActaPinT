clc;
clear;
%----- This code is used to plot the roundoff error of ParaDiag II ----
global dt T dx alp Nx Nt gamma Kmax
gamma=1/100;
aalp=[10^(-3), 10^(-6),10^(-11)];
AK=length(aalp);
T=4;
dt=1/16;
Nt=T/dt-1;  % for wave equation we use two-step time-integrator, so Nt=T/dt-1
t=0:dt:T;
dx=1/128;
Nx=1/dx;
e=ones(Nx,1);
A=1*spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
% A(1,Nx)=A(2,1);
% A(Nx,1)=A(1,2);
Ix=eye(Nx);
Kmax=10;
 
Err1=zeros(AK,Kmax); % direct implementation
Err2=zeros(AK,Kmax); % IR implementation

%------Stationary Iteration for K*U=b------
r1=Ix+(dt^2*A)/12+10*gamma*(dt^2*A)^2/12;
r2=-2*Ix+10*(dt^2*A)/12-20*gamma*(dt^2*A)^2/12;

b=zeros(Nx*Nt,1);
x=linspace(0,1,Nx)';
u0=sin(2*pi*x);
u1=(Ix+0.5*dt^2*A)\u0;
b(1:Nx)=-r2*u1-r1*u0;
b(Nx+1:2*Nx)=-r1*u1;
Uk_Ini=random('unif',-1,1,Nx*Nt,1);
for ja=1:AK
    Uk=Uk_Ini;
    alp=aalp(ja);
    for k=1:Kmax
        res=b-K_times_U(Uk,r1,r2);
        dUk=invP_times_U(res,r1,r2);
        Uk=dUk+Uk;
        Err2(ja,k)=norm(dUk,inf);
        fprintf('IR_Error at %d-th iteration is %2.15f\n',k,Err2(ja,k));
    end
end
for ja=1:AK
    Uk=Uk_Ini;
    alp=aalp(ja);
    bk=zeros(Nx*Nt,1);
    for k=1:Kmax
        bk(1:Nx)=alp*(r1*Uk((Nt-2)*Nx+1:(Nt-1)*Nx)+r2*Uk((Nt-1)*Nx+1:Nt*Nx));
        bk(Nx+1:2*Nx)=alp*r1*Uk((Nt-1)*Nx+1:Nt*Nx);
        Uk1=invP_times_U(b+bk,r1,r2);
        Err1(ja,k)=norm(Uk1-Uk,inf);
        Uk=Uk1;
        fprintf('Direct_Error at %d-th iteration is %2.15f\n',k,Err1(ja,k));
    end
end
figure(1);
semilogy(0:Kmax-1,Err1(1,:),'r-.+',0:Kmax-1,Err1(2,:),'k-.*',0:Kmax-1,Err1(3,:),'b-.o','markersize',10);shg;
xlim([0,Kmax-1]);
ylim([1e-16,max(max(Err1))]);
set(gca,'xtick',[0:Kmax-1]);
set(gca,'fontsize',14);
set(gca,'ytick',10.^(-16:2:0));
xlabel('iteration index','interpreter','latex','fontsize',20);
ylabel('Measured Error','interpreter','latex','fontsize',20);
title('Direct','fontsize',20);
leg=legend('$\alpha=10^{-3}$','$\alpha=10^{-6}$','$\alpha=10^{-11}$');
set(leg,'interpreter','latex','fontsize',15);
figure(2);
semilogy(0:Kmax-1,Err2(1,:),'r-.+',0:Kmax-1,Err2(2,:),'k-.*',0:Kmax-1,Err2(3,:),'b-.o','markersize',10);shg;
xlim([0,Kmax-1]);
ylim([1e-16,max(max(Err1))]);
set(gca,'xtick',[0:Kmax-1]);
set(gca,'fontsize',14);
set(gca,'ytick',10.^(-16:2:0));
xlabel('iteration index','interpreter','latex','fontsize',20);
ylabel('Measured Error','interpreter','latex','fontsize',20);
title('IR framework','fontsize',20);
 leg=legend('$\alpha=10^{-3}$','$\alpha=10^{-6}$','$\alpha=10^{-11}$');
set(leg,'interpreter','latex','fontsize',15);
 
function val=K_times_U(u,r1,r2)
global Nt
Nx=length(r1(:,1));
U=reshape(u,Nx,Nt);
Z=zeros(Nx,Nt);
for n=1:Nt
    if n==1
        Z(:,n)=r1*U(:,n);
    elseif n==2
        Z(:,n)=r1*U(:,n)+r2*U(:,n-1);
    else
        Z(:,n)=r1*U(:,n)+r2*U(:,n-1)+r1*U(:,n-2);
    end
end
val=reshape(Z,Nx*Nt,1);
end
function val=invP_times_U(u,r1,r2)
global Nt alp
Nx=length(r1(:,1));
U=reshape(u,Nx,Nt);
c=zeros(Nt,1);
c(2,1)=1;
tc=zeros(Nt,1);
tc(1,1)=1;
tc(3,1)=1;
Da=alp.^((0:Nt-1)'/Nt);
invDa=alp.^((0:-1:1-Nt)'/(Nt));
D=fft(Da.*c);   
tD=fft(Da.*tc);  

S1=fft(Da.*(U.')).';
S2=zeros(Nx,Nt);
for n=1:Nt 
    S2(:,n)=(tD(n)*r1+r2*D(n))\S1(:,n);
end 
val=reshape((invDa.*ifft(S2.')).',Nx*Nt,1); 
end