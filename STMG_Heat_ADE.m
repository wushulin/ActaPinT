clc;
clear;
%----- This is STMG code for heat equation and advection-diffusion equation-----

global eta l soomth_iter_num Nt Nx dt dx theta nu  Kmax
Kmax=30;
nu=0.001;
%nu=-1;
theta=1/2; % theta=0.5: trapezoidal rule; theta=1: implicit Euler BDF2
eeta=0.1:0.016:1.1;
%eeta=[0.5];
errcxt=zeros(Kmax,length(eeta));
soomth_iter_num=3;
l=8; % with the choice l=3 and
Nx=2^l-1; % N=(J+1)*2^(3*l)-1 the
Nt=Nx; % space coarsening condition
e=ones(Nx,1); % fails, try it
dx=1/(Nx+1); x=(0:dx:1)';
%A=1/dx^2*spdiags([e -2*e e],[-1 0 1],Nx,Nx);
A1=spdiags([e -e], [-1,1], Nx, Nx)/(2*dx);
A2= spdiags([e -2*e e], -1:1, Nx, Nx)/(dx^2);
Ix=eye(Nx);
if nu<0
    A=abs(nu)*A2;
else
    A=A1+nu*A2;
end

T=5; dt=T/Nt; t=(0:dt:T);
Q=speye(size(A))-dt*theta*A; % diagonal of all-at-once matrix
tQ=speye(size(A))+dt*(1-theta)*A;

u0=sin(2*pi*x);
u=zeros(Nx+2,Nt+1);
u(:,1)=u0;
for n=1:Nt % compute exact solution
    u(2:end-1,n+1)=Q\(tQ*u(2:end-1,n)); % exact linear-theta
end
u_theta=u;
Nx_c=(Nx+1)/2-1; % coarse grid size
Px=sparse(Nx,Nx_c); % prolongation by interpolation
for j=1:Nx_c
    Px(2*j,j)=1; Px(2*j-1,j)=0.5; Px(2*j+1,j)=0.5;
end
Rx=0.5*Px'; % restriction
Ac=Rx*A*Px; % coarse matrix by Getaerkin
Nt_c=(Nt+1)/2-1;
Pt=sparse(Nt,Nt_c); % prolongation by interpolation
for j=1:Nt_c
    Pt(2*j,j)=1; Pt(2*j-1,j)=0.5; Pt(2*j+1,j)=0.5;
end
Rt=0.5*Pt'; % restriction
%Pt(end,end)=1; % no final zero bc in time
Qc=speye(size(Ac))-2*dt*theta*Ac; % coarsening in time also
tQc=speye(size(Ac))+2*dt*(1-theta)*Ac; % coarsening in time also

rng('default') % random initieta guess
for je=1:length(eeta)
    eta=eeta(je);
    u(2:end-1,2:end)=rand(Nx,Nt); % with correct ic and bc
    errcxt(1,je)=max(max(abs(u_theta-u)));
    for k=1:Kmax
        for j=1:soomth_iter_num % use soomth_iter_num block Jacobi steps
            uold=u; % block Jacobi is paretalel
            res=zeros(Nx,Nt);
            for n=1:Nt
                res(:,n)=tQ*uold(2:end-1,n)-Q*uold(2:end-1,n+1);
            end
            for n=1:Nt
                u(2:end-1,n+1)=uold(2:end-1,n+1)+eta*(Q\res(:,n));
            end
        end
        %mesh(x,t,u_theta'-u'); xlabel('x'); ylabel('t');
        %title('Error after presmoothing');
        %pause
        res=zeros(Nx,Nt);
        for n=1:Nt % compute residual
            res(:,n)=tQ*uold(2:end-1,n)-Q*uold(2:end-1,n+1);
        end
        res_c=Rx*res; % restrict residueta in space
        res_c=res_c*Rt'; % restrict residueta in time
        err_c=zeros(Nx_c+2,Nt_c+1); % zero ic for correction
        for n=1:Nt_c % coarse correction
            err_c(2:end-1,n+1)=Qc\(tQc*err_c(2:end-1,n)+res_c(:,n)); % exact BE using 2*dt
        end % extend in space and time
        u(2:end-1,2:end)=u(2:end-1,2:end)+Px*err_c(2:end-1,2:end)*Pt';
        %mesh(x,t,u_theta'-u'); xlabel('x'); ylabel('t');
        %title('Error after space-time coarse correction');
        %pause
        for j=1:soomth_iter_num % use soomth_iter_num block Jacobi steps
            uold=u; % block Jacobi is paretalel
            res=zeros(Nx,Nt);
            for n=1:Nt
                res(:,n)=tQ*uold(2:end-1,n)-Q*uold(2:end-1,n+1);
            end
            for n=1:Nt
                u(2:end-1,n+1)=uold(2:end-1,n+1)+eta*(Q\res(:,n));
            end
        end
        errcxt(k+1,je)=max(max(abs(u_theta-u)));
        fprintf('EK=%d: je=%d, err at %d-th STMG iteration is=%2.15f\n',length(eeta), je,k+1,errcxt(k+1,je));
%         semilogy(0:k,errcxt(1:k+1,je),'-.o','markersize',11);shg
%         set(gca,'fontsize',14);
%         xlabel('STMG Iteration Index','fontsize',20);
%         ylabel('Global Error','fontsize',20);
%         ylim([eps,2]);
%         title(['je=',num2str(je),', Error after two-grid iteration k=',num2str(k)],'fontsize',20);
%         pause(0.5)
    end
end
d=1;
semilogy(eeta(1:d:length(eeta)),errcxt(5,1:d:length(eeta)),'r',...
    eeta(1:d:length(eeta)),errcxt(10,1:d:length(eeta)),'b:',...
    eeta(1:d:length(eeta)),errcxt(15,1:d:length(eeta)),'k--','linewidth',1);shg
xlim([0.15,1.1])
ylim([10e-5,1e+2])
set(gca,'fontname','Times New Roman','fontsize',14);
xlabel('Damping parameter $\eta$','interpreter','latex','fontsize',20);
ylabel('Error','fontsize',20);
if nu==-1
    title('Two-grid STMG for heat equation','interpreter','latex','fontsize',20);
else
    title(['Two-grid STMG for ADE with $\nu=',num2str(nu),'$'],'interpreter','latex','fontsize',20);
end
leg=legend('$k=5$', '$k=10$', '$k=15$');
set(leg,'interpreter','latex','fontsize',20);

% if nu==0.001
%     d=3;
%     semilogy(0:d:k,errcxt(1:d:k+1,je),'r--o','markersize',8);shg
%     hold on;
% elseif nu==0.01
%     d=3;
%     semilogy(0:d:k,errcxt(1:d:k+1,je),'m--s','markersize',8);shg
% elseif nu==0.1
%     d=3;
%     semilogy(0:d:k,errcxt(1:d:k+1,je),'b--d','markersize',8);shg
% else
%     d=3;
%     semilogy(0:d:k,errcxt(1:d:k+1,je),'k-.*','markersize',8);shg
% end 
% set(gca,'fontname','Times New Roman','fontsize',14);
% xlabel('Iteration Index','fontsize',20);
% ylabel('Error','fontsize',20);
% ylim([1e-10,2]);
% xlim([0,k]);
% title('1 Jacobi smoothing step','fontsize',20);
% leg=legend('ADE with $\nu=0.001$','ADE with $\nu=0.01$','ADE with $\nu=0.1$','Heat equation');
% set(leg,'interpreter','latex','fontsize',16);
