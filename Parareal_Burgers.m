clc;
clear;
%----Code of Parareal for Burgers equation u'-nu*u_{xx}+0.5*(u^2)_x=g----

global dt T dx J nu flag_eigPK Nx Nt Kmax Tol time_integrator_F
time_integrator_F=22;
Kmax=44;
Tol=1e-13;
flag_eigPK=0;
J=32;
nu=1;
T=4;
Nt=40;
dx=1/128;

Dt=T/Nt;

dt=Dt/J;
t=0:Dt:T;

Nx=1/dx;
e=ones(Nx,1);
A1=spdiags([-e e], [-1,1], Nx, Nx)/(4*dx);
A2=nu*spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
A1(1,Nx)=A1(2,1);
A1(Nx,1)=A1(1,2);
A2(1,Nx)=A2(2,1);
A2(Nx,1)=A2(1,2);
Ix=eye(Nx);

x=(linspace(0,1,Nx))';
%u0=exp(-30*(x-0.5).^2);
u0=sin(2*pi*x).^2;
u_ref=zeros(Nx,Nt+1);
u_ref(:,1)=u0;
for n=1:Nt
    u_ref(:,n+1)=Pro_F(A1,A2,u_ref(:,n),time_integrator_F);
end
U_Ini=random('unif',-1,1,Nx,Nt+1);
U_Ini(:,1)=u0;
for n=1:Nt
    U_Ini(:,n+1)=u0;
    %Pro_IE(INV_AC,U_Ini(:,n));
end
Uk=U_Ini;
Gk=Uk;
Uk1=Uk;
Err=zeros(1,Kmax);
k=1;
Err(k)=max(max(abs(Uk-u_ref)));
if time_integrator_F==1
    fprintf('time_integrator_F=implicit Euler: error at %d-th iteration is %2.15f\n',k,Err(k));
elseif time_integrator_F==21
    fprintf('time_integrator_F=TR: error at %d-th iteration is %2.15f\n',k,Err(k));
else % time_integrator_F==22
    fprintf('time_integrator_F=SDIRK22: error at %d-th iteration is %2.15f\n',k,Err(k));
end

for k=1:Kmax
    Fk=zeros(Nx,Nt);
    for n=1:Nt % parallel for Nt steps
        Fk(:,n)=Pro_F(A1,A2,Uk(:,n),time_integrator_F);
    end
    
    for n=1:Nt % sequential correcctions
        temp=Pro_IE(A1,A2,Dt,Uk1(:,n));
        Uk1(:,n+1)=Pro_IE(A1,A2,Dt,Uk1(:,n))+Fk(:,n)-Gk(:,n+1);
        Gk(:,n+1)=temp;
    end
    %Uk=Uk1;
    Err(k+1)=max(max(abs(Uk1-u_ref)));
    Uk=Uk1;
    if time_integrator_F==1
        fprintf('time_integrator_F=implicit Euler: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
    elseif time_integrator_F==21
        fprintf('time_integrator_F=TR: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
    else % time_integrator_F==22
        fprintf('time_integrator_F=SDIRK22: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
    end
    if Err(k+1)<=Tol
        break;
    end
end

% semilogy(0:k,rho_s(1:k+1)*Err(1),'b:',0:k,rho_l(1:k+1)*Err(1),'r-.',0:k,Err(1:k+1),'k--o','markersize',10,'linewidth',1);shg
% set(gca,'fontsize',14);
% xlabel('Iteration Index','fontsize',20);
% title(['$(T, N_t)=(',num2str(T),', ',num2str(Nt),'$)'],'interpreter','latex','fontsize',20);
% %ylabel('Error','fontsize',20);
% % if nu<0
% %     title('Heat equation', 'fontsize',20);
% % else
% %     title(['Advection diffusion equation with $\nu=',num2str(nu),'$'],'interpreter','latex','fontsize',20);
% % end
% leg=legend('bound $\varrho_s$','bound $\varrho_l$','measured error');
% set(leg,'interpreter','latex','fontsize',15);

 function val=Pro_F(A1,A2,Yn,time_integrator_F)
 global J dt
 if time_integrator_F==1 % F=implicit Euler
     z=Yn;
     for j=1:J
    	z=Pro_IE(A1,A2,dt,z);
     end
     val=z;
 elseif time_integrator_F==21 % F=trapezoidal rule
     z=Yn;
     for j=1:J
    	z=Pro_IE(A1,A2,0.5*dt,z-0.5*dt*(A2*z+A1*z.^2));
     end
     val=z;
 else  % F=SDIRK2
     gam=(2-sqrt(2))/2;
     z=Yn;
     for j=1:J
    	g1=Pro_IE(A1,A2,gam*dt,z);
        z=Pro_IE(A1,A2,gam*dt,z-(1-gam)*dt*(A2*g1+A1*g1.^2));
     end
     val=z;
 end
 end
function val=Pro_IE(A1,A2,tau,Yn)
global Tol
Ix=eye(length(Yn));
Kmax=1000;
z=Yn;
for k=1:Kmax
    res=z+tau*(A2*z+A1*z.^2)-Yn;
    if norm(res,inf)<=Tol/3
        break;
    else
        Jac=Ix+tau*(A2+2*A1*diag(z));
        z=z-Jac\res;
    end
end
val=z;
end

 