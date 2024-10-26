clc;
clear;
%---- This code is for Parareal with Diag-CGC for Heat and Advection-Diffusion Equation ----

global T nu Nx dx dt J alp KMAX Tol time_integrator_F;

time_integrator_F=23;
T=4;
Dt=1/10;
J=10;
dt=Dt/J;
dx=1/128;
nu=0.1;
KMAX=35;
Tol=1e-10;
alp=0.4; 
%  alp=0.25;
%  alp=0.1;
%  alp=0;
Nt=T/Dt;
x=(-1:dx:1)';
U0=sin(2*pi*x);
Nx=length(x);
Ix=speye(Nx);
e=ones(Nx,1);
A1=spdiags([-e e], [-1,1], Nx, Nx)/(4*dx);
A2=spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
A1(1,Nx)=A1(2,1);
A1(Nx,1)=A1(1,2);
A2(1,Nx)=A2(2,1);
A2(Nx,1)=A2(1,2);
if nu>0
    A=nu*A2+A1;
else
    A=A2;
end
It=speye(Nt);
invAc=Ix/(Ix+Dt*A);
%%% u_ref:  reference solution computed step-by-step
u_ref=zeros(Nx,Nt+1);
u_ref(:,1)=U0;
for n=1:Nt
    u_ref(:,n+1)=Pro_F(A,u_ref(:,n),time_integrator_F);
end
Error_Iter=[0];
%-------------------------------------------------------
GAF=alp.^((0:Nt-1)'/Nt);
invGAF=alp.^((0:-1:1-Nt)'/(Nt));  
c=[1;-1;zeros(Nt-2,1)];
D=fft(GAF.*c);   
Uk=random('unif',-1,1,Nx,Nt+1);
Uk(:,1)=U0;
Error_Iter(1)=max(max(abs(Uk-u_ref)));
fprintf('ParaDIAG2: the initial error is  %2.15f\n',Error_Iter(1));
k=0;
Gk=zeros(Nx,Nt);
Fk=zeros(Nx,Nt);
Uk1=Uk;
for k=1:KMAX
    bkm1=zeros(Nx*Nt,1);
    for n=1:Nt
        Fk(:,n)=Pro_F(A,Uk(:,n),time_integrator_F);
        if alp>0
            if n==1
                %bkm1((n-1)*Nx+1:n*Nx,1)=(Ix/Dt+A)*Fk(:,n)-alp*Uk(:,Nt+1)/Dt;
                bkm1((n-1)*Nx+1:n*Nx,1)=(Ix+Dt*A)*(Fk(:,n)-Pro_IE(invAc,alp*Uk(:,Nt+1)));
            else
                %bkm1((n-1)*Nx+1:n*Nx,1)=(Ix/Dt+A)*Fk(:,n)-Uk(:,n)/Dt;
                bkm1((n-1)*Nx+1:n*Nx,1)=(Ix+Dt*A)*(Fk(:,n)-Pro_IE(invAc,Uk(:,n)));
            end
        end
    end
    if alp>0
        bkm1=reshape(bkm1,Nx,Nt);
        sol_stepA=fft(GAF.*(bkm1.')).';
        sol_stepB=zeros(Nx,Nt);
        for n=1:Nt 
            sol_stepB(:,n)=(((D(n))*Ix+Dt*A)\sol_stepA(:,n));
        end 
        Uk=[U0,(invGAF.*ifft(sol_stepB.')).'];
    else
        for n=1:Nt
            %temp=Pro_IE(invAc,Uk1(:,n));
            Uk1(:,n+1)=Pro_IE(invAc,Uk1(:,n))+Fk(:,n)-Pro_IE(invAc,Uk(:,n));
            %Gk(:,n)=temp;
        end
        Uk=Uk1;
    end
    
    Error_Iter(k+1)=max(max(abs(Uk-u_ref)));
    fprintf('ParaDIAG2: the error at %d-th iteration is  %2.15f\n',k, Error_Iter(k+1));
    if(Error_Iter(k+1)<=Tol)
        break;
    end
end
if alp==0.4
    semilogy(0:k,Error_Iter(1:k+1),'m-.s','linewidth',1,'markersize',9);shg;
    hold on;
elseif alp==0.25
    semilogy(0:k,Error_Iter(1:k+1),'k-.+','linewidth',1,'markersize',9);shg;
elseif alp==0.1
    semilogy(0:k,Error_Iter(1:k+1),'b-.o','linewidth',1,'markersize',9);shg;
else
    semilogy(0:k,Error_Iter(1:k+1),'r-.*','linewidth',1,'markersize',9);shg;
    hold off;
end

set(gca,'fontname','Times New Roman','fontsize',14);
xlabel('Iteration Index: $k$','interpreter','latex','fontname','Times New Roman','fontsize',20);
ylabel('Error','fontname','Times New Roman','fontsize',20);
% if time_integrator_F==22
%     title('$\mathcal{F}$=2nd-order SDIRK','interpreter','latex','fontsize',22);
% end
% if time_integrator_F==1
%     title('$\mathcal{F}$=implicit Euler','interpreter','latex','fontsize',22);
% end
title('Heat equation','fontsize',22);
xlim([0,20]);
ylim([Tol,max(Error_Iter)]);
set(gca,'ytick',10.^(-10:0));
set(gca,'xtick',0:2:20);
leg=legend('Diag-based CGC: $\alpha=0.4$','Diag-based CGC: $\alpha=0.25$',...
    'Diag-based CGC: $\alpha=0.1$','standard CGC');
set(leg,'fontsize',15,'interpreter','latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
 elseif time_integrator_F==22 % F=SDIRK22
     gam=(2-sqrt(2))/2;
     INV_AF=(E+gam*dt*A)\E;
     z=Yn;
     for j=1:J
    	g1=Pro_IE(INV_AF,z);
        z=Pro_IE(INV_AF,z-(1-gam)*dt*A*g1);
     end
     val=z;
 elseif time_integrator_F==23 % F=SDIRK23
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
 else % F=SDIRK23 (another (2,3)-SDIRK method)
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
 end
 end
 
function val=Pro_IE(invA,Yn)
val=invA*Yn;
end

 