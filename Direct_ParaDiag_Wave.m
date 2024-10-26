clc;
clear;
%----Direct ParaDiag for 2nd-order wave equation u_{tt}=u_{xx}+g----
global  dx T
NN=[5,10,15,20,30];
%ee=0.01:0.022:1;
ee=1.11.^(-47:0);
It_e=length(ee);
It_N=length(NN);
dx=0.1/2;
T=0.2;
x=(0:dx:1)';
Nx=length(x);
e=ones(Nx,1);
A2= spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
% A2(1,Nx)=A2(2,1);
% A2(Nx,1)=A2(1,2);
Ix=speye(Nx);
A=[0*Ix,-Ix;A2,0*Ix];
Nx=2*Nx;
Ix=speye(Nx);

y0=[sin(2*pi*x);x*0];
error_dt=zeros(It_N);
error_dtn=zeros(It_N,It_e);
% varrho_opt=[0.010573790659212   0.0253392422702238   0.0682874151453803   0.095481730604595   0.119223400860857];
varrho_opt=[0.014775019077301   0.050815171349603   0.080940125907062   0.102404302317876   0.132337516569827];
error_dtn_opt=zeros(1,It_N);
for jNt=1:It_N
    Nt=NN(jNt);
    val=get_varrho_opt(Nt);
    if Nt==NN(end)
        varrho_opt(jNt)=0.4;
    else
        varrho_opt(jNt)=val;
    end
    It=eye(Nt);
    dt=T/Nt;
    invA=(Ix+0.5*dt*A)\(Ix-0.5*dt*A);
    y_ref=zeros(Nx,Nt);
    y_ref(:,1)=y0;
    y_TR=zeros(Nx,Nt);
    u_ref=zeros(Nx/2,Nt);
    u_TR=u_ref;
    for n=1:Nt
        y_ref(:,n)=expm(-dt*n*A)*y0;
        u_ref(:,n)=y_ref(1:Nx/2,n);
        if n==1
            y_TR(:,n)=invA*y0;
        else
            y_TR(:,n)=invA*y_TR(:,n-1);
        end
        u_TR(:,n)=y_TR(1:Nx/2,n);
    end
    error_dt(jNt)=max(max(abs(u_ref-u_TR)));   
    
    B=0.5*toeplitz([1;1;zeros(Nt-2,1)],[1,zeros(1,Nt-1)]);
    Bt=zeros(Nt,Nt);
    for je=1:It_e
        mu=1+ee(je);
        dtn=(T/sum(mu.^(1:Nt)))*mu.^(1:Nt);
        for n=1:Nt
            Bt(n,n)=1/dtn(n);
        end
        for n=2:Nt
            Bt(n,n-1)=-1/dtn(n);
        end
        b=zeros(Nt*Nx,1);
        b(1:Nx)=(1/dtn(1))*y0;
        tb=kron(B\It,Ix)*b;
        [V,D]=eig(B\Bt);
        Ya=kron(It/V,Ix)*tb;
        Yb=zeros(Nx*Nt,1);
        for n=1:Nt
            Yb((n-1)*Nx+1:n*Nx)=(D(n,n)*Ix+A)\Ya((n-1)*Nx+1:n*Nx);
        end
        Yc=kron(V,Ix)*Yb;    
        yc=reshape(Yc,Nx,Nt);
        for n=1:Nt
            u_TR(:,n)=yc(1:Nx/2,n);
        end
        y_ref=zeros(Nx,Nt);
        y_ref(:,1)=y0;
        for n=1:Nt
            y_ref(:,n)=expm(-sum(dtn(1:n))*A)*y0;
            u_ref(:,n)=y_ref(1:Nx/2,n);
        end
        error_dtn(jNt,je)=max(max(abs(u_ref-u_TR)));
    end
    
    Bt=zeros(Nt,Nt);
    mu=1+varrho_opt(jNt);
    dtn=(T/sum(mu.^(1:Nt)))*mu.^(1:Nt);
    for n=1:Nt
        Bt(n,n)=1/dtn(n);
    end
    for n=2:Nt
        Bt(n,n-1)=-1/dtn(n);
    end
    b=zeros(Nt*Nx,1);
    b(1:Nx)=(1/dtn(1))*y0;
    tb=kron(B\It,Ix)*b;
    [V,D]=eig(B\Bt);
    Ya=kron(It/V,Ix)*tb;
    Yb=zeros(Nx*Nt,1);
    for n=1:Nt
        Yb((n-1)*Nx+1:n*Nx)=(D(n,n)*Ix+A)\Ya((n-1)*Nx+1:n*Nx);
    end
    Yc=kron(V,Ix)*Yb;    

    yc=reshape(Yc,Nx,Nt);
    for n=1:Nt
        u_TR(:,n)=yc(1:Nx/2,n);
    end
    y_ref=zeros(Nx,Nt);
    y_ref(:,1)=y0;
    for n=1:Nt
        y_ref(:,n)=expm(-sum(dtn(1:n))*A)*y0;
        u_ref(:,n)=y_ref(1:Nx/2,n);
    end
    error_dtn_opt(jNt)=max(max(abs(u_ref-u_TR)));

    fprintf('It_Nt=%d: jNt=%d\n',It_N,jNt);
end
for jn=1:It_N
    loglog(ee,error_dtn(jn,:),'linewidth',1);shg
    hold on; 
end
loglog(varrho_opt,error_dtn_opt(1:It_N),'*','markersize',11);
hold off; 
set(gca,'fontsize',15);
xlabel('$\varrho$','interpreter','latex','fontsize',20);
ylabel('measured error','interpreter','latex','fontsize',20);
title('wave equation','interpreter','latex','fontsize',20);
ylim([min(min(error_dtn))/2,1]);shg
xlim([min(ee),max(ee)]);
function val=get_varrho_opt(Nt)
val=((1e-7)*15*2^(2*Nt-0.5)/((Nt^2-1)*factorial(Nt-1))).^(1/(1+Nt));
end