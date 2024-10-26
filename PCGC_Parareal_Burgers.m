clc;
clear;
%----- This code is for Parareal with Diag-based CGC for Burgers equation u'-nu*u_{xx}+0.5*(u^2)_x=g -----

global dt T dx J nu Nx Nt Kmax Tol time_integrator_F alp
alp=0.1;
time_integrator_F=23;
Kmax=44;
Tol=1e-10;
J=10;
nu=0.01;
T=4;
dx=1/128;

Dt=0.1;
Nt=T/Dt;
dt=Dt/J;
t=0:Dt:T;

Nx=1/dx;
e=ones(Nx,1);
A1=1*spdiags([-e e], [-1,1], Nx, Nx)/(4*dx);
A2=nu*spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
A1(1,Nx)=A1(2,1);
A1(Nx,1)=A1(1,2);
A2(1,Nx)=A2(2,1);
A2(Nx,1)=A2(1,2);
Ix=eye(Nx);
It=eye(Nt);
x=(linspace(0,1,Nx))';
u0=sin(2*pi*x);

u_ref=zeros(Nx,Nt+1);
u_ref(:,1)=u0;
for n=1:Nt
    u_ref(:,n+1)=Pro_F(A1,A2,u_ref(:,n),time_integrator_F);
end
U_Ini=random('unif',-1,1,Nx,Nt+1);
U_Ini(:,1)=u0;
% for n=2:Nt+1
%     U_Ini(:,n)=u0;
% end
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
else
    fprintf('time_integrator_F=SDIRK23: error at %d-th iteration is %2.15f\n',k,Err(k));
end
if alp==0
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
        elseif time_integrator_F==22
            fprintf('time_integrator_F=SDIRK22: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
        else
            fprintf('time_integrator_F=SDIRK23: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
        end
        if Err(k+1)<=Tol
            break;
        end
    end
else
    B=toeplitz([1;-1;zeros(Nt-2,1)],[1,zeros(1,Nt-1)]);
    C=B;
    C(1,Nt)=C(2,1)*alp;
    GAF=alp.^((0:Nt-1)'/Nt);
    invGAF=alp.^((0:-1:1-Nt)'/(Nt));  
    c=[1;-1;zeros(Nt-2,1)];
    D=fft(GAF.*c);   
    for k=1:Kmax
        Fk=zeros(Nx,Nt+1);
        Gk=zeros(Nx,Nt+1);
        bk=zeros(Nx,Nt+1);
        res=zeros(Nx,Nt);
        for n=1:Nt % parallel for Nt steps
            Fk(:,n+1)=Pro_F(A1,A2,Uk(:,n),time_integrator_F);
            if n==1
                Gk(:,n+1)=Pro_IE(A1,A2,Dt,alp*Uk(:,Nt+1)+u0);
            else
                Gk(:,n+1)=Pro_IE(A1,A2,Dt,Uk(:,n));
            end
            bk(:,n+1)=Fk(:,n+1)-Gk(:,n+1);
        end   
        Z=Uk(:,2:Nt+1); % do not include the initial value
        for l=1:100 % Simplified Newton Iteration
            res=zeros(Nx,Nt);
            for n=1:Nt
                fn=A2*(Z(:,n)-bk(:,n+1))+A1*(Z(:,n)-bk(:,n+1)).^2;
                if n==1
                    res(:,n)=(Z(:,n)-bk(:,n+1))-alp*Z(:,Nt)-u0+Dt*fn;
                else
                    res(:,n)=(Z(:,n)-bk(:,n+1))-Z(:,n-1)+Dt*fn;
                end
            end
            if norm(reshape(res,Nx*Nt,1),inf)<=max(Tol,10^(-k))
                break;
            else
                Jac=zeros(Nx,Nx);
                for n=1:Nt
                    Jac=Jac+(A2+2*A1*diag(Z(:,n)-bk(:,n+1)))/Nt;
                end
                sol_stepA=fft(GAF.*(res.')).';
                sol_stepB=zeros(Nx,Nt);
                for n=1:Nt 
                    sol_stepB(:,n)=(((D(n))*Ix+Dt*Jac)\sol_stepA(:,n));
                end 
                Z=Z-(invGAF.*ifft(sol_stepB.')).';
            end
        end
       Uk(:,2:Nt+1)=Z;
        Err(k+1)=max(max(abs(Uk-u_ref)));
        if time_integrator_F==1
            fprintf('time_integrator_F=implicit Euler: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
        elseif time_integrator_F==21
            fprintf('time_integrator_F=TR: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
        elseif time_integrator_F==22
            fprintf('time_integrator_F=SDIRK22: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
        else
            fprintf('time_integrator_F=SDIRK23: error at %d-th iteration is %2.15f\n',k+1,Err(k+1));
        end
        if Err(k+1)<=Tol
            break;
        end
    end
end
if alp==0.4
    semilogy(0:k,Err(1:k+1),'m-.s','linewidth',1,'markersize',9);shg;
    hold on;
elseif alp==0.25
    semilogy(0:k,Err(1:k+1),'k-.+','linewidth',1,'markersize',9);shg;
elseif alp==0.1
    semilogy(0:k,Err(1:k+1),'b-.o','linewidth',1,'markersize',9);shg;
else
    semilogy(0:k-1,[Err(1),Err(3:k+1)],'r-.*','linewidth',1,'markersize',9);shg;
    hold off;
end

set(gca,'fontname','Times New Roman','fontsize',14);
xlabel('Iteration Index: $k$','interpreter','latex','fontname','Times New Roman','fontsize',20);
ylabel('Error','fontname','Times New Roman','fontsize',20);
title(['Buergers equation with $\nu=',num2str(nu),'$'],'interpreter','latex','fontsize',21);
xlim([0,20]);
ylim([Tol,max(Err)]);
set(gca,'ytick',10.^(-10:0));
set(gca,'xtick',0:2:20);
leg=legend('Diag-based CGC: $\alpha=0.4$','Diag-based CGC: $\alpha=0.25$',...
    'Diag-based CGC: $\alpha=0.1$','standard CGC');
set(leg,'fontsize',15,'interpreter','latex');

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
 elseif time_integrator_F==22  % F=SDIRK2
     gam=(2-sqrt(2))/2;
     z=Yn;
     for j=1:J
    	g1=Pro_IE(A1,A2,gam*dt,z);
        z=Pro_IE(A1,A2,gam*dt,z-(1-gam)*dt*(A2*g1+A1*g1.^2));
     end
     val=z;
 else
     gam=(3+sqrt(3))/6;
     z=Yn;
     for j=1:J
    	g1=Pro_IE(A1,A2,gam*dt,z);
        g2=Pro_IE(A1,A2,gam*dt,z-(1-2*gam)*dt*(A2*g1+A1*g1.^2));
        z1=z-dt*(A2*(g1+g2)+A1*(g1.^2+g2.^2))/2;
        z=z1;
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

 
 