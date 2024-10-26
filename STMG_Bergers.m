clc;
clear;

%----- This is STMG code for Burgers equation u'-nu*u_{xx}+0.5*(u^2)_x=g----

global eta soomth_iter_num Nt Nx dt dx theta A1 A2 u0 A1c A2c NT_Kmax NT_Tol STMG_Kmax
nu=1;
%eeta=0.1:0.016:1.1;
eeta=[0.25];
theta=1; % theta=0.5: trapezoidal rule; theta=1: implicit Euler BDF2
soomth_iter_num=4;
NT_Kmax=100;
STMG_Kmax=16;
NT_Tol=1e-10;
l=8; % with the choice l=3 and
Nx=2^l-1; % N=(J+1)*2^(3*l)-1 the
Nt=Nx; % space coarsening condition
Nx_c=(Nx+1)/2-1; % coarse grid size
Nt_c=(Nt+1)/2-1;
e=ones(Nx,1); % fails, try it
dx=1/(Nx+1); x=(0:dx:1)';
A1=spdiags([e -e], [-1,1], Nx, Nx)/(4*dx);
A2=nu*spdiags([e -2*e e], -1:1, Nx, Nx)/(dx^2);
T=5; dt=T/Nt; t=(0:dt:T);
u0=sin(2*pi*x(2:Nx+1));
u=zeros(Nx,Nt+1);
u(:,1)=u0;
for n=1:Nt % compute exact solution
    z0=u(:,n);
    for k=1:NT_Kmax
        res=dt*(theta*f(z0,A1,A2)+(1-theta)*f(u(:,n),A1,A2))-(z0-u(:,n));
        if norm(res,inf)<=NT_Tol
            break;
        else
            Jac=eye(Nx)-dt*theta*df(z0,A1,A2);
            z0=z0+Jac\res;
        end
    end
    u(:,n+1)=z0; 
end
u_theta=reshape(u(:,2:Nt+1),Nt*Nx,1);


Px=sparse(Nx,Nx_c); % prolongation by interpolation
for j=1:Nx_c
    Px(2*j,j)=1; Px(2*j-1,j)=0.5; Px(2*j+1,j)=0.5;
end
Rx=0.5*Px'; % restriction
A1c=Rx*A1*Px;A2c=Rx*A2*Px; % coarse matrix by Getaerkin
Pt=sparse(Nt,Nt_c); % prolongation by interpolation
for j=1:Nt_c
    Pt(2*j,j)=1; Pt(2*j-1,j)=0.5; Pt(2*j+1,j)=0.5;
end
Rt=0.5*Pt'; % restriction
%Pt(end,end)=1; % no final zero bc in time

rng('default') % random initieta guess
U_old=random('unif',-1,1,Nx*Nt,1); % with correct ic and bc
err_STMG=zeros(STMG_Kmax,length(eeta));

for je=1:length(eeta)
    eta=eeta(je);
    k=1;
    err_STMG(k,je)=norm(U_old-u_theta,inf);
    fprintf('EK=%d: je=%d, err at %d-th STMG iteration is=%2.15f\n',length(eeta), je,k,err_STMG(k,je));
    for k=1:STMG_Kmax
        eta=0.28;
        U_1over3=S(U_old,soomth_iter_num);
        res=reshape(get_res(U_1over3),Nx,Nt);
        res_c=reshape(Rx*reshape(res,Nx,Nt)*Rt',Nt_c*Nx_c,1); % restrict residueta in space-time
        Uc_1over3=reshape(Rx*reshape(U_1over3,Nx,Nt)*Rt',Nx_c*Nt_c,1);
        Uc_2over3=solve_Kc(res_c,Uc_1over3);
        err_c=Uc_2over3-Uc_1over3;
        err=reshape(Px*reshape(err_c,Nx_c,Nt_c)*Pt',Nx*Nt,1);
        U_2over3=U_1over3+err;
        U_old=S(U_2over3,soomth_iter_num);
        %U_old=U_1over3+err;
        err_STMG(k+1,je)=norm(U_old-u_theta,inf);
        fprintf('EK=%d: je=%d, err at %d-th STMG iteration is=%2.15f\n',length(eeta), je,k+1,err_STMG(k+1,je));
    end
end
% d=1;
% semilogy(eeta(1:d:length(eeta)),err_STMG(5,1:d:length(eeta)),'r',...
%     eeta(1:d:length(eeta)),err_STMG(10,1:d:length(eeta)),'b:',...
%     eeta(1:d:length(eeta)),err_STMG(15,1:d:length(eeta)),'k--','linewidth',1);shg
% xlim([0.15,1.1])
% ylim([10e-5,1e+2])
% set(gca,'fontname','Times New Roman','fontsize',14);
% xlabel('Damping parameter $\eta$','interpreter','latex','fontsize',20);
% ylabel('Error','fontsize',20);
% title(['Two-grid STMG for Burgers equation with $\nu=',num2str(nu),'$'],'interpreter','latex','fontsize',20);

d=2;
semilogy(0:d:k,err_STMG(1:d:k+1),'b--*','markersize',10);shg
ylim([1e-6,3])
set(gca,'fontname','Times New Roman','fontsize',14);
xlabel('Iteration Index','fontsize',20);
ylabel('Error','fontsize',20);
xlim([0,k]);
title(['Two-grid STMG for Burgers equation with $\nu=',num2str(nu),'$'],'interpreter','latex','fontsize',20);
leg=legend('$\nu=0.1$','$\nu=1$');
set(leg,'interpreter','latex','fontsize',16);


function val=f(un,A1,A2)
val=A2*un+A1*un.^2;
end
function val=df(un,A1,A2)
val=A2+2*A1*diag(un);
end

function val=get_res(U)
global theta dt Nt Nx u0 A1 A2
res=zeros(Nx*Nt,1);
res(1:Nx)=dt*(theta*f(U(1:Nx),A1,A2)+(1-theta)*f(u0,A1,A2))-(U(1:Nx)-u0);
for n=2:Nt
    res((n-1)*Nx+1:n*Nx)=dt*(theta*f(U((n-1)*Nx+1:n*Nx),A1,A2)+(1-theta)*f(U((n-2)*Nx+1:(n-1)*Nx),A1,A2))...
        -(U((n-1)*Nx+1:n*Nx)-U((n-2)*Nx+1:(n-1)*Nx));
end
val=res;
end



function val=S(U_old,soomth_iter_num)
global theta eta  dt Nt Nx NT_Tol NT_Kmax A1 A2 
U0=U_old; 
for j=1:soomth_iter_num
    b0=get_res(U0);
    b=eta*(b0);
    for n=1:Nt
        C0=zeros(Nx,1);
        for k=1:NT_Kmax
            res=b((n-1)*Nx+1:n*Nx)-(C0-theta*dt*f(C0,A1,A2));
            if norm(res,inf)<=NT_Tol
                break;
            else
                Jac=eye(Nx)-theta*dt*df(C0,A1,A2);
                C0=C0+Jac\res;
            end
        end
        U0((n-1)*Nx+1:n*Nx)=U0((n-1)*Nx+1:n*Nx)+C0;
    end
end
val=U0;
end

function val=solve_Kc(res_c,U_c)
global Nx Nt dt theta A1c A2c NT_Tol NT_Kmax
Nx_c=(Nx+1)/2-1; % coarse grid size
Nt_c=(Nt+1)/2-1;
dt_c=2*dt;
b_c=zeros(Nx_c*Nt_c,1);
n=1;
b_c((n-1)*Nx_c+1:n*Nx_c)=res_c((n-1)*Nx_c+1:n*Nx_c)+U_c((n-1)*Nx_c+1:n*Nx_c)-theta*dt_c*f(U_c((n-1)*Nx_c+1:n*Nx_c),A1c,A2c);
val=zeros(Nt_c*Nx_c,1);
Un=U_c((n-1)*Nx_c+1:n*Nx_c);
for k=1:NT_Kmax
    res=b_c((n-1)*Nx_c+1:n*Nx_c)-(Un-theta*dt_c*f(Un,A1c,A2c));
    if norm(res,inf)<=NT_Tol
        break;
    else
        Jac=eye(Nx_c)-theta*dt_c*df(Un,A1c,A2c);
        Un=Un+Jac\res;
    end
    val((n-1)*Nx_c+1:n*Nx_c)=Un;
end
for n=2:Nt_c
    b_c((n-1)*Nx_c+1:n*Nx_c)=res_c((n-1)*Nx_c+1:n*Nx_c)+U_c((n-1)*Nx_c+1:n*Nx_c)-U_c((n-2)*Nx_c+1:(n-1)*Nx_c)-...
        (theta*dt_c*f(U_c((n-1)*Nx_c+1:n*Nx_c),A1c,A2c)+(1-theta)*dt_c*f(U_c((n-2)*Nx_c+1:(n-1)*Nx_c),A1c,A2c));
    Un=U_c((n-1)*Nx_c+1:n*Nx_c);
    for k=1:NT_Kmax
        res=b_c((n-1)*Nx_c+1:n*Nx_c)-(Un-U_c((n-2)*Nx_c+1:(n-1)*Nx_c)-...
            (theta*dt_c*f(Un,A1c,A2c)+(1-theta)*dt_c*f(U_c((n-2)*Nx_c+1:(n-1)*Nx_c),A1c,A2c)));
        if norm(res,inf)<=NT_Tol
            break;
        else
            Jac=eye(Nx_c)-theta*dt_c*df(Un,A1c,A2c);
            Un=Un+Jac\res;
        end
        val((n-1)*Nx_c+1:n*Nx_c)=Un;
    end
end
end



 