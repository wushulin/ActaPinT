clc;
clear;
%----Code of Parareal and MGRiT for Burgers equation u'-nu*u_{xx}+0.5*(u^2)_x=g----

global dt T dx J nu Nx Nt Kmax Tol time_integrator_F
time_integrator_F=22;
Kmax=60;
Tol=1e-14;
J=10;
nu=1/2;
T=5;
Dt=1/16;
Nt=T/Dt;
dx=1/160;

dt=Dt/J;
t=0:Dt:T;

Nx=1/dx;
e=ones(Nx,1);
A1=spdiags([-e e], [-1,1], Nx, Nx)/(4*dx);
A2=nu*spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
% A1(1,Nx)=A1(2,1);
% A1(Nx,1)=A1(1,2);
% A2(1,Nx)=A2(2,1);
% A2(Nx,1)=A2(1,2);
Ix=eye(Nx);

x=(linspace(0,1,Nx))';
u0=sin(8*pi*(1-x).^2).^2;
u_ref=zeros(Nx,Nt+1);
u_ref(:,1)=u0;
for n=1:Nt
    u_ref(:,n+1)=Pro_F(A1,A2,u_ref(:,n),time_integrator_F);
end
rng('default')
U_Ini=random('unif',0,1,Nx,Nt+1);
U_Ini(:,1)=u0;
Uk=U_Ini;
Gk=Uk;
Uk1=Uk;
Err_parareal=zeros(1,Kmax);
k=1;
Err_parareal(k)=max(max(abs(Uk-u_ref)));
if time_integrator_F==1
    fprintf('Parareal: time_integrator_F=implicit Euler, error at %d-th iteration is %2.15f\n',k,Err_parareal(k));
elseif time_integrator_F==21
    fprintf('Parareal: time_integrator_F=TR, error at %d-th iteration is %2.15f\n',k,Err_parareal(k));
else % time_integrator_F==22
    fprintf('Parareal: time_integrator_F=SDIRK22, error at %d-th iteration is %2.15f\n',k,Err_parareal(k));
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
    Err_parareal(k+1)=max(max(abs(Uk1-u_ref)));
    Uk=Uk1;
    if time_integrator_F==1
        fprintf('Parareal: time_integrator_F=implicit Euler, error at %d-th iteration is %2.15f\n',k+1,Err_parareal(k+1));
    elseif time_integrator_F==21
        fprintf('Parareal: time_integrator_F=TR, error at %d-th iteration is %2.15f\n',k+1,Err_parareal(k+1));
    else % time_integrator_F==22
        fprintf('Parareal: time_integrator_F=SDIRK22, error at %d-th iteration is %2.15f\n',k+1,Err_parareal(k+1));
    end
    if Err_parareal(k+1)<=Tol
        break;
    end
end
k_parareal=k;

Uk=U_Ini;
Uk(:,2)=u_ref(:,2);
Uk1=Uk;
Err_mgrit=zeros(1,Kmax);
k=1;
Err_mgrit(k)=max(max(abs(Uk-u_ref)));
if time_integrator_F==1
    fprintf('MGRiT: time_integrator_F=implicit Euler, error at %d-th iteration is %2.15f\n',k,Err_mgrit(k));
elseif time_integrator_F==21
    fprintf('MGRiT: time_integrator_F=TR: error at %d-th iteration is %2.15f\n',k,Err_mgrit(k));
else % time_integrator_F==22
    fprintf('MGRiT: time_integrator_F=SDIRK22: error at %d-th iteration is %2.15f\n',k,Err_mgrit(k));
end

for k=1:Kmax
    Fk=zeros(Nx,Nt);
    Gk=zeros(Nx,Nt);
    for n=2:Nt-1 % parallel for Nt steps
        Fk1=Pro_F(A1,A2,Uk(:,n-1),time_integrator_F);
        Fk(:,n+1)=Pro_F(A1,A2,Fk1,time_integrator_F);
        Gk(:,n+1)=Pro_IE(A1,A2,Dt,Fk1);
    end
    
    for n=2:Nt-1 % sequential correcctions
        Uk1(:,n+1)=Pro_IE(A1,A2,Dt,Uk1(:,n))+Fk(:,n+1)-Gk(:,n+1);
    end
    Err_mgrit(k+1)=max(max(abs(Uk1(:,3:Nt)-u_ref(:,3:Nt))));
    Uk=Uk1;
    if time_integrator_F==1
        fprintf('MGRiT: time_integrator_F=implicit Euler, error at %d-th iteration is %2.15f\n',k+1,Err_mgrit(k+1));
    elseif time_integrator_F==21
        fprintf('MGRiT: time_integrator_F=TR, error at %d-th iteration is %2.15f\n',k+1,Err_mgrit(k+1));
    else % time_integrator_F==22
        fprintf('MGRiT: time_integrator_F=SDIRK22, error at %d-th iteration is %2.15f\n',k+1,Err_mgrit(k+1));
    end
    if Err_mgrit(k+1)<=Tol
        break;
    end
end
Err_parareal0=Err_parareal(1:2:k_parareal);
k_mgrit=k;
semilogy(0:length(Err_parareal0)-1,Err_parareal0,'b--o',...
    0:k_mgrit,Err_mgrit(1:k_mgrit+1),'r--+',...
    'markersize',10,'linewidth',1);shg
set(gca,'fontsize',15);
xlabel('Iteration Index','fontsize',21);
ylabel('Error','fontsize',21);
xlim([0,25]);
ylim([1e-14,1]);
set(gca,'ytick',10.^(-14:2:0));
set(gca,'xtick',0:5:25);
title(['Burgers,2 equation with $\nu=',num2str(nu),'$'],'interpreter','latex','fontsize',20);
leg=legend('Parareal','MGRiT');
set(leg,'fontsize',15);






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
Kmax=50;
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

 