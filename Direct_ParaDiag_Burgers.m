clc;
clear;
%---- Direct ParaDiag method for Burgers equation u'-nu*u_{xx}+0.5*(u^2)_x=g ----
global nu T dt dx Kmax NT_Tol   Nx Nt
Kmax=30;
Nt=200;
nu=0.002;
dx=0.01;
x=(0:dx:1)';
Nx=length(x);
Ix=eye(Nx);
It=eye(Nt);
e=ones(Nx,1);
A1 = 0.5*spdiags([-e e], [-1,1], Nx, Nx)/(2*dx);
A2= nu*spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
A1(1,Nx)=A1(2,1);
A1(Nx,1)=A1(1,2);
A2(1,Nx)=A2(2,1);
A2(Nx,1)=A2(1,2);
u0=sin(2*pi*x);
U_TR=zeros(Nx,Nt+1);
NewtonIt_TR=zeros(1,Nt);
U_TR(:,1)=u0;
NT_Tol=max(dx^2,dt^2)/100;

TT=[0.1,0.2,0.4,0.8,1.6];
T_Num=length(TT);
Jac_Num_TR=zeros(1,T_Num);
Fevu_Num_TR=zeros(1,T_Num);
Jac_Num_ParaDiagI=zeros(1,T_Num);
Fevu_Num_ParaDiagI=zeros(1,T_Num);
It_ParaDiagI=zeros(1,T_Num);
B0=toeplitz([0;-0.5;zeros(Nt-2,1)],[0,0.5,zeros(1,Nt-2)]);
B0(Nt,Nt-1)=-1;B0(Nt,Nt)=1;
Error=zeros(T_Num,Kmax);
for jT=1:T_Num
    T=TT(jT);
    dt=T/Nt;
    for n=1:Nt
        z=U_TR(:,n);
        for k=1:Kmax
            res=U_TR(:,n)-0.5*dt*(A2*U_TR(:,n)+A1*U_TR(:,n).^2)-(z+0.5*dt*(A2*z+A1*z.^2));
            if norm(res,inf)<=NT_Tol
                break;
            else
                Jac=Ix+0.5*dt*(A2+2*A1*diag(z));
                z=z+Jac\res;
            end
        end
        U_TR(:,n+1)=z;
        NewtonIt_TR(n)=k;
    end
    Jac_Num_TR(jT)=sum(NewtonIt_TR);
    Fevu_Num_TR(jT)=Jac_Num_TR(jT)+Nt;


    B=B0/dt;
    [V,D]=eig(B);
    invV=It/V;
    b=zeros(Nx*Nt,1);
    b(1:Nx)=u0/(2*dt);
    Uk=random('unif',0,1,Nx*Nt,1);
%     for n=1:Nt
%         Uk((n-1)*Nx+1:n*Nx,1)=u0;
%     end

    fk=zeros(Nx*Nt,1);
    for n=1:Nt
        fk((n-1)*Nx+1:n*Nx)=f(Uk((n-1)*Nx+1:n*Nx),A1,A2);
    end
    % k=1;
    % Error(k)=norm(dt*(b-stepAC(B,Uk)-fk),inf);
    % fprintf('Error at %d-th iteration is %2.15f\n',k,Error(k));

    for k=1:Kmax
    %     fk=zeros(Nx*Nt,1);
    %     for n=1:Nt
    %         fk((n-1)*Nx+1:n*Nx)=f(Uk((n-1)*Nx+1:n*Nx),A1,A2);
    %     end
    %     Error(k)=norm(dt*(b-stepAC(B,Uk)-fk),inf);
            %Ak=zeros(Nx,Nx);
        avg_Uk=zeros(Nx,1);
        for n=1:Nt
            %Ak=Ak+df(Uk((n-1)*Nx+1:n*Nx),A1,A2)/Nt;
            avg_Uk=avg_Uk+Uk((n-1)*Nx+1:n*Nx);
        end
        Ak=df(avg_Uk/Nt,A1,A2);
        bk=zeros(Nx*Nt,1);
        for n=1:Nt
            bk((n-1)*Nx+1:n*Nx)=Ak*Uk((n-1)*Nx+1:n*Nx)-f(Uk((n-1)*Nx+1:n*Nx),A1,A2);
        end
        res=b+bk;
        g=stepAC(invV,res);
        w=zeros(Nx*Nt,1);
        for n=1:Nt
            w((n-1)*Nx+1:n*Nx)=(D(n,n)*Ix+Ak)\g((n-1)*Nx+1:n*Nx);
        end
        Uk1=stepAC(V,w);
        Error(jT,k)=norm(Uk-Uk1,inf);
        fprintf('jT=%d: Error at %d-th iteration is %2.15f\n',jT,k,Error(jT,k));
        Uk=Uk1;
        if Error(jT,k)<=NT_Tol
            Jac_Num_ParaDiagI(jT)=k;
            break;
        elseif Error(jT,k)>=1e+8
            Jac_Num_ParaDiagI(jT)=NaN;
            break;
        end
    end
    It_ParaDiagI(jT)=k;
    Fevu_Num_ParaDiagI(jT)=k;
end
figure(1);
for jT=1:T_Num
    k=It_ParaDiagI(jT);
    dt=TT(jT)/Nt;
    if jT==1
        semilogy(0:k-1,Error(jT,1:k),'r:+','markersize',11,'linewidth',1);
    elseif jT==2
        semilogy(0:k-1,Error(jT,1:k),'c:o','markersize',11,'linewidth',1);
    elseif jT==3
        semilogy(0:k-1,Error(jT,1:k),'b:s','markersize',11,'linewidth',1);
    elseif jT==4
        semilogy(0:k-1,Error(jT,1:k),'k:*','markersize',11,'linewidth',1);
    else
        semilogy(0:k-1,Error(jT,1:k),'m:d','markersize',11,'linewidth',1);
        %semilogy(1:k,Error(jT,1:k),'m-.d',1:k,ones(1,k)*max(dx^2,dt^2),'m-','markersize',11);
    end
    hold on;
end
hold off;
set(gca,'fontsize',15);
xlabel('Iteration Index','fontsize',20);
ylabel('Error','fontsize',20);
% title('advection diffusion equation','fontsize',20);
title(['$\nu=',num2str(nu),'$'],'fontsize',20,'interpreter','latex');
xlim([0,max(It_ParaDiagI)-1]);
if nu==0.1
    ylim([NT_Tol,2*max(max(Error))]);
else
    ylim([NT_Tol,100]);
end
% leg=legend('$(T,\Delta t)=(0.1, 0.0005)$','$(T,\Delta t)=(0.2, 0.001)$',...
%     '$(T,\Delta t)=(0.4, 0.002)$',...
%     '$(T,\Delta t)=(0.8, 0.004)$','$(T,\Delta t)=(1.6, 0.008)$');
leg=legend('$T=0.1$','$T=0.2$',...
    '$T=0.4$',...
    '$T=0.8$','$T=1.6$');

set(leg,'fontsize',15,'interpreter','latex');
figure(2);
loglog(TT,Jac_Num_TR,'r-.o',TT,Jac_Num_ParaDiagI,'b-.+','markersize',11,'linewidth',1);
set(gca,'fontsize',15);
xlabel('$T$','fontsize',20,'interpreter','latex');
ylabel('Total Jacobian Solves','fontsize',20);
title(['$\nu=',num2str(nu),'$'],'fontsize',20,'interpreter','latex');
xlim([min(TT)/1.5,max(TT)*1.5]);
set(gca,'xtick',TT);
leg=legend('trapezoidal rule','ParaDiag I');
set(leg,'fontsize',15);
% Nx=4;Nt=5;
% V=rand(Nt,Nt);
% U=random('unif',-1,1,Nx*Nt,1);
% Ix=eye(Nx);
% s1=kron(V,Ix)*U-stepAC(V,U);
function val=f(u,A1,A2)
val=A1*u.^2+A2*u;
end
function val=df(u,A1,A2)
val=2*A1*diag(u)+A2;
end

function val=stepAC(V,U)
global Nx Nt
val=zeros(Nx*Nt,1);
for n=1:Nt
    c=0;
    for j=1:Nt
        c=c+V(n,j)*U((j-1)*Nx+1:j*Nx);
    end
    val((n-1)*Nx+1:n*Nx)=c;
end
end

 