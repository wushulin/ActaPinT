clc;
clear;
%----this code is for u'=nu*u_{xx}+u_x+f and heat equation u'=u_{xx}+f----
global dx dt T nu a flag_BC flag_source x
flag_source=0;
flag_BC=1;
nu=5e-4;
a=2/2;
if a==0
    nu=1;
end
T=3;
dt=0.001;
t=0:dt:T;
Nt=length(t);
dx=1*dt;
x=(0:dx:1)';
c=1/2;
if flag_BC==1 % Dirichlet BCs
    Nx=length(x)-2;
    Ix=speye(Nx);
    e=ones(Nx,1);
    A1 =a*spdiags([-e e], [-1,1], Nx, Nx)/(2*dx); 
    A2= nu*spdiags([e -2*e e], -1:1, Nx, Nx)/dx^2;
    A=A2-A1;
    u_ref=zeros(Nx+2,Nt); 
    if flag_source==0
        u_ref(:,1)=sin(8*pi*(1-x).^2).^2;
    end
    invA2=(Ix-dt*A2)\Ix;
    for n=1:Nt-1
        u_ref(2:Nx+1,n+1)=invA2*(u_ref(2:Nx+1,n)-dt*A1*u_ref(2:Nx+1,n)+dt*g(t(n+1)));
    end
elseif flag_BC==2 % Numman BCs
    Nx=length(x);
    Ix=speye(Nx);
    e=ones(Nx,1);
    A2=nu*spdiags([e -2*e e], -1:1, Nx, Nx)/dx^2;
    A2(1,1)=-1/dx^2;
    A2(Nx,Nx)=-1/dx^2;

  
    A1=a*spdiags([-e e], [-1,1], Nx, Nx)/(2*dx);
    A1(1,1)=a*(-4/3)/(2*dx);A1(1,2)=a*(4/3)/(2*dx);
    A1(Nx,Nx-1)=a*(-4/3)/(2*dx);A1(Nx,Nx)=a*(4/3)/(2*dx);

    A=A2-A1;
    u_ref=zeros(Nx,Nt); 
    if flag_source==0
        u_ref(:,1)=sin(8*pi*(1-x).^2).^2;
    end
    invA=(Ix-c*dt*A)\Ix;
    for n=1:Nt-1
        u_ref(:,n+1)=invA*((Ix+(1-c)*dt*A)*u_ref(:,n)+dt*((1-c)*g(t(n))+c*g(t(n+1))));
    end
    % invA2=(Ix-dt*A2)\Ix;
    % for n=1:Nt-1
    %     u_ref(:,n+1)=invA2*(u_ref(:,n)-dt*A1*u_ref(:,n)+dt*g(t(n+1)));
    % end
else % Periodic BCs
    Nx=length(x);
    Ix=speye(Nx);
    e=ones(Nx,1);
    A1 =a*spdiags([-e e], [-1,1], Nx, Nx)/(2*dx);
    A2= nu*spdiags([e -2*e e], -1:1, Nx, Nx)/dx^2;
    A2(1,Nx)=A2(2,1);
    A2(Nx,1)=A2(1,2);
    A1(1,Nx)=A1(2,1);
    A1(Nx,1)=A1(1,2);
    A=A2-A1;
    invA=(Ix-c*dt*A)\Ix;
    u_ref=zeros(Nx,Nt);
    if flag_source==0
        u_ref(:,1)=sin(8*pi*(1-x).^2).^2;
    end
    for n=1:Nt-1
        g1=g(t(n));
        g2=g(t(n+1));
        u_ref(:,n+1)=invA*((Ix+(1-c)*dt*A)*u_ref(:,n)+dt*((1-c)*g(t(n))+c*g(t(n+1))));
    end

    invA2=(Ix-dt*A2)\Ix;
    for n=1:Nt-1
        u_ref(:,n+1)=invA2*(u_ref(:,n)-dt*A1*u_ref(:,n)+dt*g(t(n+1)));
    end

end
figure(1)
dxt=2;
mesh(t(1:dxt:end),x(1:dxt:end),u_ref(1:dxt:end,1:dxt:end));shg
view(90,-90);
 set(gca,'fontsize',16);
xlabel('t','fontsize',20);
ylabel('x','fontsize',20);
set(gca,'ytick',0:0.2:1);
set(gca,'xtick',0:0.5:3);
if flag_BC==1
    title('Dirichlet BCs','fontsize',20,'fontweight','normal');
elseif flag_BC==2
    title('Neumann BCs','fontsize',20,'fontweight','normal');
else
    title('Periodic BCs','fontsize',20,'fontweight','normal');
end

%contourf(x,t,transpose(u_ref));shg
