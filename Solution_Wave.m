clc;
clear;
%----this code is for 2nd-order wave equation u_{tt}=c*u_{xx}+f----
global dx dt T nu flag_BC x flag_source
flag_source=1
flag_BC=1
c=0.2;
T=3;
dt=1/400;
t=0:dt:T;
Nt=length(t);
dx=dt;
x=(0:dx:1)';
if flag_BC==1 % Dirichlet BCs
    Nx=length(x)-2;
    u_ref=zeros(Nx+2,Nt); 
    if flag_source==0
        u_ref(:,1)=sin(8*pi*(1-x).^2).^2;
    end
    Ix=speye(Nx);
    e=ones(Nx,1);
    A2= c*spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
    A=[0*Ix,-Ix;A2,0*Ix];
    Nx=length(A(:,1));
    Ix=speye(Nx);
    invA=(Ix+0.5*dt*A)\Ix;
    %z=zeros(Nx,1);
    z=[u_ref(2:length(x)-1,1);0*u_ref(2:length(x)-1,1)];
    for n=1:Nt-1
        fn=g(t(n));
        z=invA*((Ix-0.5*dt*A)*z+0.5*dt*([0*e;g(t(n))]+[0*e;g(t(n+1))]));
        u_ref(2:length(e)+1,n+1)=z(1:length(e));
    end
elseif flag_BC==2 % Neumann BCs
    Nx=length(x);
    u_ref=zeros(Nx,Nt); 
    if flag_source==0
        u_ref(:,1)=sin(8*pi*(1-x).^2).^2;
    end
    Ix=speye(Nx);
    e=ones(Nx,1);
    A2= c*spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx^2;
    A2(1,1)=nu/dx^2;A2(Nx,Nx)=A2(1,1);

    % A2(1,1)=c*(2/3)/dx^2;A2(1,2)=-c*(2/3)/dx^2;
    % A2(Nx,Nx-1)=-c*(2/3)/dx^2;A2(Nx,Nx)=c*(2/3)/dx^2;
     A=[0*Ix,-Ix;A2,0*Ix];
    Nx=length(A(:,1));
    Ix=speye(Nx);
    invA=(Ix+0.5*dt*A)\Ix;
    %z=zeros(Nx,1);
    z=[u_ref(:,1);0*u_ref(:,1)];
    for n=1:Nt-1
        z=invA*((Ix-0.5*dt*A)*z+0.5*dt*([0*e;g(t(n))]+[0*e;g(t(n+1))]));
        u_ref(:,n+1)=z(1:length(e));
    end  
else % Periodic BCs
    Nx=length(x);
    u_ref=zeros(Nx,Nt); 
    if flag_source==0
        u_ref(:,1)=sin(8*pi*(1-x).^2).^2;
    end
    Ix=speye(Nx);
    e=ones(Nx,1);
    A2= c*spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx^2;
    A2(1,Nx)=A2(2,1);
    A2(Nx,1)=A2(1,2);
    A=[0*Ix,-Ix;A2,0*Ix];
    Nx=length(A(:,1));
    Ix=speye(Nx);
    invA=(Ix+0.5*dt*A)\Ix;
    %z=zeros(Nx,1);
    z=[u_ref(:,1);0*u_ref(:,1)];
    for n=1:Nt-1
        z=invA*((Ix-0.5*dt*A)*z+0.5*dt*([0*e;g(t(n))]+[0*e;g(t(n+1))]));
        u_ref(:,n+1)=z(1:length(e));
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


 
