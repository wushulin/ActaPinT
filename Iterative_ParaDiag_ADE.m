clc;
clear;
%----Iterative ParaDiag for advection-diffusion equation u'-nu*u_{xx}+u_x=g----
global dt T dx theta alp nu flag_eigPK Nx Nt flag_GMRES Kmax Tol
Kmax=20;
Tol=1e-12;
flag_GMRES=1;
flag_eigPK=0;
nu=1e-3;
nu=1e-6;
alp=1;
T=2;
dt=1/50;
Nt=T/dt; 
t=0:dt:T;
dx=1/100;
Nx=1/dx;
theta=1;
e=ones(Nx,1);
A1=spdiags([-e e], [-1,1], Nx, Nx)/(2*dx);
A2= spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
% A1(1,Nx)=A1(2,1);
% A1(Nx,1)=A1(1,2);
% A2(1,Nx)=A2(2,1);
% A2(Nx,1)=A2(1,2);
Ix=eye(Nx);
if nu<0
    A=A2;
else
    A=A1+nu*A2;
end
%------- plot the eigenvalues of inv(P)*K -------
if flag_eigPK==1
    eigA=eig(full(A));
    eigPK=zeros(Nt*Nx,1);
    B=toeplitz([0;1;zeros(Nt-2,1)],zeros(1,Nt));
    C=B;
    C(1,Nt)=alp;
    It=eye(Nt);
    for j=1:Nx
        lam=eigA(j);
        r1=1+dt*lam*theta;
        r2=-1+dt*(1-theta)*lam;
        K=r1*It+r2*B;
        P=r1*It+r2*C;
        eigPK((j-1)*Nt+1:j*Nt,1)=eig(P\K);
        if mod(j,10)==0
            fprintf('Nx=%d: j=%d\n',Nx,j);
        end
    end
    % r2=1/(1-alp)-1;
    % cc=0:0.01:2*pi;
    % c2=r2*exp(1i*cc);
%     plot(real(eigPK),imag(eigPK),'+',real(c2),imag(c2),'r-');shg
    figure(1);
    plot(real(eigPK),imag(eigPK),'k+','markersize',10);shg
    %hold on;
    
    set(gca,'fontsize',14);
    xlabel('Re$(\lambda(\mathcal{P}^{-1}\mathcal{K})$','interpreter','latex','fontsize',20);
    ylabel('Im$(\lambda(\mathcal{P}^{-1}\mathcal{K})$','interpreter','latex','fontsize',20);
    %hold on;
    if nu<0
        title('Heat equation','fontsize',20);
    else
        title('Advection diffusion equation','fontsize',20);
    end
 
end

%------GMRES solver for K*U=b------
r1=Ix+dt*theta*A;
r2=-Ix+dt*(1-theta)*A;
b=zeros(Nx*Nt,1);
x=linspace(0,1,Nx)';
u0=sin(2*pi*x);
b(1:Nx)=-r2*u0;
if flag_GMRES==1
    [X,FLAG,RELRES,ITER,RESVEC]=gmres(@(u)K_times_U(u,r1,r2),b,[],Tol,1000,@(u)invP_times_U(u,r1,r2),[],kron(ones(Nt,1),0.125*u0));
    u_ref=[u0,reshape(X,Nx,Nt)];
else
    Uk=kron(ones(Nt,1),0.125*u0);
    RESVEC=zeros(1,Kmax);
    for k=1:Kmax
        res=b-K_times_U(Uk,r1,r2);
        dUk=invP_times_U(res,r1,r2);
        Uk=dUk+Uk;
        RESVEC(k)=norm(dUk,inf);
        fprintf('Error at %d-th iteration is %2.15f\n',k,RESVEC(k));
        if RESVEC(k)<=Tol
            break;
        end
    end
end
figure(2);
if nu==1e-3
    semilogy(0:length(RESVEC)-1,RESVEC,'k-.s','markersize',10);shg;
else
    semilogy(0:length(RESVEC)-1,RESVEC,'r-.o','markersize',10);shg;
end

xlim([0,length(RESVEC)-1]);
ylim([1e-14,max(RESVEC)]);
set(gca,'xtick',[0:length(RESVEC)-1]);
set(gca,'fontsize',14);
set(gca,'ytick',10.^(-14:2:0));
xlabel('iteration index','interpreter','latex','fontsize',20);
if flag_GMRES==1
    ylabel('GMRES residual','interpreter','latex','fontsize',20);
    if nu<0
        title('Heat equation','fontsize',20);
    else
        title('Advection diffusion equation','fontsize',20);
    end
else
    ylabel('Error $\max_n\|{\bf u}_n^k-{\bf u}_n\|$','interpreter','latex','fontsize',20);
    leg=legend('Advection diffusion equation with $\nu=10^{-6}$');
    set(leg,'interpreter','latex','fontsize',20);
end
    
function val=K_times_U(u,r1,r2)
global Nt
Nx=length(r1(:,1));
U=reshape(u,Nx,Nt);
Z=zeros(Nx,Nt);
for n=1:Nt
    if n==1
        Z(:,n)=r1*U(:,n);
    else
        Z(:,n)=r1*U(:,n)+r2*U(:,n-1);
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
Da=alp.^((0:Nt-1)'/Nt);
invDa=alp.^((0:-1:1-Nt)'/(Nt));
D=fft(Da.*c);   
 
S1=fft(Da.*(U.')).';
S2=zeros(Nx,Nt);
for n=1:Nt 
        S2(:,n)=(r1+r2*D(n))\S1(:,n);
end 
val=reshape((invDa.*ifft(S2.')).',Nx*Nt,1); 

end