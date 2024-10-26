clc;
clear;
%----Iterative ParaDiag for 2nd-order wave equation u_{tt}=u_{xx}+g----
global dt T dx alp flag_eigPK Nx Nt gamma flag_GMRES Kmax Tol
flag_GMRES=0;
flag_eigPK=1;
gamma=1/120.01;
alp=0.02;
T=20;
dt=1/16;
Nt=T/dt-1;  % for wave equation we use two-step time-integrator, so Nt=T/dt-1
t=0:dt:T;
dx=1/128;
Nx=1/dx;
e=ones(Nx,1);
A=spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
% A(1,Nx)=A(2,1);
% A(Nx,1)=A(1,2);
Ix=eye(Nx);

%------- plot the eigenvalues of inv(P)*K -------
if flag_eigPK==1
    eigA=eig(full(A));
    eigPK=zeros(Nt*Nx,1);
    B=toeplitz([0;1;zeros(Nt-2,1)],zeros(1,Nt));
    tB=toeplitz([1;0;1;zeros(Nt-3,1)],[1,zeros(1,Nt-1)]);
    C=B; C(1,Nt)=alp;
    tC=tB;tC(1,Nt-1)=alp;tC(2,Nt)=alp;
    It=eye(Nt);
    for j=1:Nx
        lam=eigA(j);
        r1=1+(dt^2*lam)/12+10*gamma*(dt^2*lam)^2/12;
        r2=-2+10*(dt^2*lam)/12-20*gamma*(dt^2*lam)^2/12;
        K=r1*tB+r2*B;
        P=r1*tC+r2*C;
        eigPK((j-1)*Nt+1:j*Nt,1)=eig(P\K);
        if mod(j,10)==0
            fprintf('Nx=%d: j=%d\n',Nx,j);
        end
    end
    figure(1);
    r2=alp/(1-alp);
    cc=0:0.01:2*pi;
    c2=r2*exp(1i*cc);
    plot(real(eigPK-1),imag(eigPK),'k+',real(c2),imag(c2),'r--','markersize',10,'linewidth',1.5);shg;axis equal;
%     plot(real(eigPK),imag(eigPK),'+');shg
    set(gca,'fontsize',14);
    xlabel('Re$(\lambda(\mathcal{M})$','interpreter','latex','fontsize',20);
    ylabel('Im$(\lambda(\mathcal{M})$','interpreter','latex','fontsize',20);
    %title('Wave equation','fontsize',20);
    title(['$T=',num2str(T),'$'],'interpreter','latex','fontsize',20);
    xlim([-3.85*alp/(1-alp),1.2*alp/(1-alp)]);
%     leg=legend('$\lambda(\mathcal{M})$','circle with radius $\frac{\alpha}{1-\alpha}$');
%     set(leg,'interpreter','latex','fontsize',15);
end

% %------GMRES solver for K*U=b------
% r1=Ix+(dt^2*A)/12+10*gamma*(dt^2*A)^2/12;
% r2=-2*Ix+10*(dt^2*A)/12-20*gamma*(dt^2*A)^2/12;
% 
% b=zeros(Nx*Nt,1);
% x=linspace(0,1,Nx)';
% u0=sin(2*pi*x);
% u1=(Ix+0.5*dt^2*A)\u0;
%  
% b(1:Nx)=-r2*u1-r1*u0;
% b(Nx+1:2*Nx)=-r1*u1;
% if flag_GMRES==1
%     [X,FLAG,RELRES,ITER,RESVEC]=gmres(@(u)K_times_U(u,r1,r2),b,[],1e-14,1000,@(u)invP_times_U(u,r1,r2),[],kron(ones(Nt,1),0.125*u0));
%     u_ref=[u0,u1,reshape(X,Nx,Nt)];
% else
%     Uk=kron(ones(Nt,1),0.125*u0);
%     for k=1:Kmax
%         res=b-K_times_U(Uk,r1,r2);
%         dUk=invP_times_U(res,r1,r2);
%         Uk=dUk+Uk;
%         RESVEC(k)=norm(dUk,inf);
%         fprintf('Error at %d-th iteration is %2.15f\n',k,RESVEC(k));
%         if RESVEC(k)<=Tol
%             break;
%         end
%     end
% end
% figure(2);
% semilogy(0:length(RESVEC)-1,RESVEC,'k-.+');shg;
% xlim([0,length(RESVEC)-1]);
% ylim([1e-14,max(RESVEC)]);
% set(gca,'xtick',[0:5:length(RESVEC)-1]);
% set(gca,'fontsize',14);
% set(gca,'ytick',10.^(-14:2:0));
% xlabel('iteration index','interpreter','latex','fontsize',20);
% ylabel('GMRES residual','interpreter','latex','fontsize',20);
% title('Wave equation','fontsize',20);
%  
 
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