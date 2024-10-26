clc;
clear;
%----Direct ParaDiag for advection-diffusion equation u'-nu*u_{xx}+u_x=g----
global  dx nu T eigA
NN=[5,10,15,20,30];
ee=0.01:0.01:1;
It_e=length(ee);
It_N=length(NN);
dx=0.1/5;
adc=0;
if adc==0
    nu=1;
else
    nu=0.01;
end
T=0.2;
x=(0:dx:1)';
Nx=length(x);
e=ones(Nx,1);
A1 = spdiags([-e e], [-1,1], Nx, Nx)/(2*dx);
A2= nu*spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
A1(1,Nx)=A1(2,1);
A1(Nx,1)=A1(1,2);
A2(1,Nx)=A2(2,1);
A2(Nx,1)=A2(1,2);
A=adc*A1+A2;
Ix=speye(Nx);

eigA=eig(full(A));
y0=sin(2*pi*x);
error_dt=zeros(It_N);
error_dtn=zeros(It_N,It_e);
% varrho_opt=[0.010573790659212   0.0253392422702238   0.0682874151453803   0.095481730604595   0.119223400860857];
varrho_opt=[0.014775019077301   0.050815171349603   0.080940125907062   0.102404302317876   0.132337516569827];
%zeros(1,It_N);
error_dtn_opt=zeros(1,It_N);
for jNt=1:It_N
    Nt=NN(jNt);
    val=get_varrho_opt(Nt);
    if Nt==NN(end)
        varrho_opt(jNt)=0.16;
    else
        varrho_opt(jNt)=val;
    end
    It=eye(Nt);
    dt=T/Nt;
    invA=Ix/(Ix+dt*A);
    y_ref=zeros(Nx,Nt);
    y_ref(:,1)=y0;
    y_IE=zeros(Nx,Nt);
    for n=1:Nt
        y_ref(:,n)=expm(-dt*n*A)*y0;
        if n==1
            y_IE(:,n)=invA*y0;
        else
            y_IE(:,n)=invA*y_IE(:,n-1);
        end
    end
    error_dt(jNt)=max(max(abs(y_ref-y_IE)));   
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
        [V,D]=eig(Bt);
        Ya=kron(It/V,Ix)*b;
        Yb=zeros(Nx*Nt,1);
        for n=1:Nt
            Yb((n-1)*Nx+1:n*Nx)=(D(n,n)*Ix+A)\Ya((n-1)*Nx+1:n*Nx);
        end
        Yc=kron(V,Ix)*Yb;    
        
        y_ref=zeros(Nx,Nt);
        y_ref(:,1)=y0;
        for n=1:Nt
            y_ref(:,n)=expm(-sum(dtn(1:n))*A)*y0;
        end
        error_dtn(jNt,je)=max(abs(reshape(y_ref,Nx*Nt,1)-Yc)); 
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
    [V,D]=eig(Bt);
    Ya=kron(It/V,Ix)*b;
    Yb=zeros(Nx*Nt,1);
    for n=1:Nt
        Yb((n-1)*Nx+1:n*Nx)=(D(n,n)*Ix+A)\Ya((n-1)*Nx+1:n*Nx);
    end
    Yc=kron(V,Ix)*Yb;    

    y_ref=zeros(Nx,Nt);
    y_ref(:,1)=y0;
    for n=1:Nt
        y_ref(:,n)=expm(-sum(dtn(1:n))*A)*y0;
    end
    error_dtn_opt(jNt)=max(abs(reshape(y_ref,Nx*Nt,1)-Yc)); 

    fprintf('It_Nt=%d: jNt=%d\n',It_N,jNt);
end
for jn=1:It_N
    loglog(ee,error_dtn(jn,:),'linewidth',1);shg
    hold on; 
end
loglog(varrho_opt,error_dtn_opt(1:It_N),'*','markersize',11);
ylim([5e-2,2]);shg
hold off; 
set(gca,'fontsize',15);
xlabel('$\varrho$','interpreter','latex','fontsize',20);
ylabel('measured error','interpreter','latex','fontsize',20);
if adc==0
    title('heat equation','interpreter','latex','fontsize',20);
else
    title(['advection diffusion equation with $\nu=',num2str(nu),'$'],'interpreter','latex','fontsize',20);
end
ylim([min(min(min(error_dtn)),min(error_dtn_opt))/1.96,1]);shg
xlim([min(ee),max(ee)]);
function val=get_varrho_opt(Nt)
global eigA T
val=((eps*Nt^2*(1+2*Nt)*(Nt+max(abs(eigA))*T))/(get_phi(Nt)*get_C(Nt))).^(1/(1+Nt));
end
function val=get_phi(Nt)
if mod(Nt,2)==0
    val=factorial(Nt/2)*factorial(Nt/2-1);
else
    val=(factorial((Nt-1)/2))^2;
end
end
function val=get_C(Nt)
global T eigA
x=(abs(eigA))*T/Nt;
r=((x./(1+x)).^2).*((1+x).^(-Nt));
val=(Nt*(Nt^2-1)/24)*max(r);
end