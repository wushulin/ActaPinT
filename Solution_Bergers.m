clc;
clear;
%----this code is for Burgers equation u'-nu*u_{xx}+0.5*(u^2)_x=g----
global dx dt T nu flag_BC  Tol NTmax x flag_source
flag_source=0;
flag_BC=1;
nu=5e-4;
T=3;
dt=1/400;
t=0:dt:T;
Nt=length(t);
dx=dt;
Tol=max(dt^2,dx^2)/10;
NTmax=100;
x=(0:dx:1)';
c=1/2;
NT_It=zeros(1,Nt);
if flag_BC==1 % Dirichlet BCs
    Nx=length(x)-2;
    Ix=speye(Nx);
    e=ones(Nx,1);
    A1 =spdiags([-e e], [-1,1], Nx, Nx)/(4*dx);
    A2= nu*spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
    u_ref=zeros(Nx+2,Nt); 
    if flag_source==0
        u_ref(:,1)=sin(8*pi*(1-x).^2).^2;
    end
    for n=1:Nt-1
        z=u_ref(2:Nx+1,n);
        for k=1:NTmax
            res=z-u_ref(2:Nx+1,n)+ dt*A2*(c*z+(1-c)*u_ref(2:Nx+1,n))+dt*A1*(c*z.^2+(1-c)*u_ref(2:Nx+1,n).^2)-dt*((1-c)*g(t(n))+c*g(t(n+1)));
            if norm(res,inf)<=Tol
                break;
            else
                Jac=Ix+dt*c*A2+2*dt*c*A1*diag(z);
                z=z-Jac\res;
            end
        end
        NT_It(n+1)=k;
        u_ref(2:Nx+1,n+1)=z;
    end
elseif flag_BC==2 % Neumann BCs
    Nx=length(x);
    Ix=speye(Nx);
    e=ones(Nx,1);
    A2= nu*spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx^2;
    A2(1,1)=nu*(2/3)/dx^2;A2(1,2)=-nu*(2/3)/dx^2;
    A2(Nx,Nx-1)=-nu*(2/3)/dx^2;A2(Nx,Nx)=nu*(2/3)/dx^2;

    A1 = 1*spdiags([-e e], [-1,1], Nx, Nx)/(2*dx);
    A1(1,1)=(-4/3)/(2*dx);A1(1,2)=(4/3)/(2*dx);
    A1(Nx,Nx-1)=(-4/3)/(2*dx);A1(Nx,Nx)=(4/3)/(2*dx);
    
    % a=1;
    % A1=a*spdiags([-e e], [-1,0], Nx, Nx)/(dx);
    % A1(1,1)=a*(-1/3)/dx;A1(1,2)=a*(1/3)/dx;
    % A1(Nx,Nx-1)=a*(-1/3)/dx;A1(Nx,Nx)=a*(1/3)/dx;

    u_ref=zeros(Nx,Nt);
    if flag_source==0
        u_ref(:,1)=sin(8*pi*(1-x).^2).^2;
    end
    for n=1:Nt-1
        z=u_ref(:,n);
        for k=1:NTmax
            res=z-u_ref(:,n)+ dt*A2*(c*z+(1-c)*u_ref(:,n))+dt*A1*(c*z.^2+(1-c)*u_ref(:,n).^2)-dt*((1-c)*g(t(n))+c*g(t(n+1)));
            if norm(res,inf)<=Tol
                break;
            else
                Jac=Ix+dt*c*A2+2*dt*c*A1*diag(z);
                z=z-Jac\res;
            end
        end
        NT_It(n+1)=k;
        u_ref(:,n+1)=z;
    end
else % Periodic BCs
    Nx=length(x);
    Ix=speye(Nx);
    e=ones(Nx,1);
    A1 = 1*spdiags([-e e], [-1,1], Nx, Nx)/(4*dx);
    A2= nu*spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx^2;
    A2(1,Nx)=A2(2,1);
    A2(Nx,1)=A2(1,2);
    A1(1,Nx)=A1(2,1);
    A1(Nx,1)=A1(1,2);
    u_ref=zeros(Nx,Nt);
    if flag_source==0
        u_ref(:,1)=sin(8*pi*(1-x).^2).^2;
    end
    for n=1:Nt-1
        z=u_ref(:,n);
        for k=1:NTmax
            res=z-u_ref(:,n)+ dt*A2*(c*z+(1-c)*u_ref(:,n))+dt*A1*(c*z.^2+(1-c)*u_ref(:,n).^2)-dt*((1-c)*g(t(n))+c*g(t(n+1)));
            if norm(res,inf)<=Tol
                break;
            else
                Jac=Ix+dt*c*A2+2*dt*c*A1*diag(z);
                z=z-Jac\res;
            end
        end
        NT_It(n+1)=k;
        u_ref(:,n+1)=z;
    end
end
mesh(t,x,u_ref);shg
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

 