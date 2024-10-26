clc;
clear;
%----- This is Parareal code for heat equation and advection-diffusion equation-----

global dt T dx J nu flag_eigPK Nx Nt Kmax Tol time_integrator_F
time_integrator_F=22;
Kmax=50;
Tol=2e-14;
flag_eigPK=0;
J=32;
nu=0.02;
% T=0.02;
% Nt=6;

% T=0.5;
% Nt=64;

%dx=1/5;
T=4;
Nt=40;
dx=1/128;

Dt=T/Nt;

dt=Dt/J;
t=0:Dt:T;

Nx=1/dx;
e=ones(Nx,1);
A1=spdiags([-e e], [-1,1], Nx, Nx)/(2*dx);
A2= spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
A1(1,Nx)=A1(2,1);
A1(Nx,1)=A1(1,2);
A2(1,Nx)=A2(2,1);
A2(Nx,1)=A2(1,2);
Ix=eye(Nx);
if nu<0
    A=abs(nu)*A2;
else
    A=A1+nu*A2;
end
z=eig(full(Dt*A));
zf=z/J;
Rg=1./(1+z);
if time_integrator_F==1
    Rf=(1./(1+zf)).^J;
elseif time_integrator_F==21
    Rf=((1-0.5*zf)./(1+0.5*zf)).^J;
elseif time_integrator_F==22
    gam=(2-sqrt(2))/2;
    Af=[gam,0;1-gam,gam];
    bf=[1-gam,gam];
    sf=2;
    Rf=zeros(Nx,1);
    for j=1:Nx
        Rf(j)=((1-zf(j)*bf*inv(eye(sf)+zf(j)*Af)*ones(sf,1)))^J;
    end
elseif time_integrator_F==31
    gam=0.4358665215;
    tau=0.5*(1+gam);
    b1=-0.25*(6*gam^2-16*gam+1);
    b2=0.25*(6*gam^2-20*gam+5);
    Af=[gam,0,0;tau-gam,gam,0;
        b1,b2,gam];
    bf=[b1,b2,gam];
    sf=3;
    Rf=zeros(Nx,1);
    for j=1:Nx
        Rf(j)=((1-zf(j)*bf*inv(eye(sf)+zf(j)*Af)*ones(sf,1)))^J;
    end
else
    gam=(3+sqrt(3))/6;
    Rf=zeros(Nx,1);
    for j=1:Nx
        Rf(j)=((sqrt(3)+zf(j)-gam*zf(j)^2)/(sqrt(3)*(1+gam*zf(j))^2))^J;
    end
end
rho_l=zeros(1,Kmax); % convergence factor for large Nt
rho_s=zeros(1,Kmax); % convergence factor for small Nt
rho_l_z=zeros(1,Nx);
for j=1:Nx
    rho_l_z(j)=abs(Rg(j)-Rf(j))./(1-abs(Rg(j)));
end
plot3(real(z/Dt),imag(z/Dt),rho_l_z,'+');shg
set(gca,'fontsize',14);
xlabel('Re($\lambda(A)$)', 'interpreter','latex','fontsize',20);
ylabel('Im($\lambda(A)$)', 'interpreter','latex','fontsize',20);
zlabel('$\varrho_l(J,z)$', 'interpreter','latex','fontsize',20);
title(['$\nu=',num2str(nu),'$'],'interpreter','latex','fontsize',20);
zlim([0,1]);
xlim([0,max(real(z)/Dt)]);
ylim([min(imag(z)/Dt),max(imag(z)/Dt)]);
text(100,0,0.4,'$\varrho_{l,\max}=0.23$','interpreter','latex','fontsize',20);
for k=0:Kmax-1
    rho_l(k+1)=(max(abs(Rg-Rf)./(1-abs(Rg))))^k;
end
for k=0:Kmax-1
    rho_s(k+1)=(max(abs(Rg-Rf))^k)*(prod(Nt-k:Nt-1)/factorial(k));
end

x=(linspace(0,1,Nx))';
%u0=exp(-30*(x-0.5).^2);
u0=sin(2*pi*x).^2;
u_ref=zeros(Nx,Nt+1);
u_ref(:,1)=u0;
for n=1:Nt
    u_ref(:,n+1)=Pro_F(A,u_ref(:,n),time_integrator_F);
end
INV_AC=eye(Nx)/(eye(Nx)+Dt*A);
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
elseif time_integrator_F==22
    fprintf('time_integrator_F=SDIRK22: error at %d-th iteration is %2.15f\n',k,Err(k));
elseif time_integrator_F==31
    fprintf('time_integrator_F=SDIRK33: error at %d-th iteration is %2.15f\n',k,Err(k));
else
    fprintf('time_integrator_F=SDIRK23: error at %d-th iteration is %2.15f\n',k,Err(k));
end

for k=1:Kmax
    Fk=zeros(Nx,Nt);
    for n=1:Nt % parallel for Nt steps
        Fk(:,n)=Pro_F(A,Uk(:,n),time_integrator_F);
    end
    
    for n=1:Nt % sequential correcctions
        temp=Pro_IE(INV_AC,Uk1(:,n));
        Uk1(:,n+1)=Pro_IE(INV_AC,Uk1(:,n))+Fk(:,n)-Gk(:,n+1);
        Gk(:,n+1)=temp;
    end
    %Uk=Uk1;
    Err(k+1)=max(max(abs(Uk1-u_ref)));
    Uk=Uk1;
    if time_integrator_F==1
        fprintf('time_integrator_F=implicit Euler: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
    elseif time_integrator_F==21
        fprintf('time_integrator_F=TR: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
    elseif time_integrator_F==22
        fprintf('time_integrator_F=SDIRK22: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
    elseif time_integrator_F==31
        fprintf('time_integrator_F=SDIRK33: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
    else
        fprintf('time_integrator_F=SDIRK23: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
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

 function val=Pro_F(A,Yn,time_integrator_F)
 global Nx J dt
 E=eye(Nx);
 if time_integrator_F==1 % F=implicit Euler
     z=Yn;
     INV_AF=E/(E+dt*A);
     for j=1:J
    	z=Pro_IE(INV_AF,z);
     end
     val=z;
 elseif time_integrator_F==21 % F=trapezoidal rule
     z=Yn;
     INV_AF=(E+0.5*dt*A)\E;
     for j=1:J
    	z=Pro_IE(INV_AF,(E-0.5*dt*A)*z);
     end
     val=z;
 elseif time_integrator_F==22 % F=SDIRK2
     gam=(2-sqrt(2))/2;
     INV_AF=(E+gam*dt*A)\E;
     z=Yn;
     for j=1:J
    	g1=Pro_IE(INV_AF,z);
        z=Pro_IE(INV_AF,z-(1-gam)*dt*A*g1);
     end
     val=z;
 elseif time_integrator_F==31 % F=SDIRK3_a
     gam=0.4358665215;
     tau=0.5*(1+gam);
     b1=-0.25*(6*gam^2-16*gam+1);
     b2=0.25*(6*gam^2-20*gam+5);
     INV_AF=(E+gam*dt*A)\E;
     z=Yn;
     for j=1:J
    	g1=Pro_IE(INV_AF,z);
        g2=Pro_IE(INV_AF,z-(tau-gam)*dt*A*g1);
        z=Pro_IE(INV_AF,z-b1*dt*A*g1-b2*dt*A*g2);
     end
     val=z;
 else % F=SDIRK3_b
     gam=(3+sqrt(3))/6;
     INV_AF=(E+gam*dt*A)\E;
     z=Yn;
     for j=1:J
    	g1=Pro_IE(INV_AF,z);
        g2=Pro_IE(INV_AF,z-(1-2*gam)*dt*A*g1);
        z1=z-dt*A*(g1+g2)/2;
        z=z1;
     end
     val=z;
 end
 end
function val=Pro_IE(invA,Yn)
val=invA*Yn;
end

 