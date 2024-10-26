% construct wave matrix using kron
clc
clear;
global flag_source flag_BC x 
flag_BC=1;
flag_source=1;
c=sqrt(1);
n=400;
X=1;
ex=ones(n,1);
dx=X/(n+1);
DDx=1/dx^2*spdiags([ex -2*ex ex],[-1:1],n,n);  % second derivative space
T=3;
dt=dx/c;            % dt on CFL, need to be slightly below
m=ceil(T/dt);       % compute needed time steps to be just below  
%m=2*ceil(T/dt);       % compute needed time steps to be just below  
dt=T/m;             % dt to be used
et=ones(m,1);

DDt=1/(c^2*dt^2)*spdiags([et -2*et et],[-2:0],m,m);  % second derivative time

A=kron(DDt,speye(size(DDx)))-kron(spdiags(et,-1,m,m),DDx);

% add initial conditions on the right hand side, using Taylor

t=dt:dt:T;
x=(dx:dx:1-dx)';
 
u0=sin(8*pi*(1-x).^2).^2;
 
u0t=zeros(size(x));

%--------Construct the right hand-side f---------
x=(0:dx:1)';
f=zeros(size(A,1),1);
j=1;
f((j-1)*n+1:j*n)=u0t/dt+u0/dt^2+1/2*DDx*u0;
j=2;
f((j-1)*n+1:j*n)=-u0/dt^2;


%------------------------------------------------
x=(dx:dx:1-dx)';
ue=A\f;
%for k=1:m
%  plot(x,ue((k-1)*n+1:k*n),'-');
%  axis([0 1 min(ue) max(ue)]);
%  pause
%end
  
Ue=reshape(ue,n,m);
mesh(t,x,Ue);colorbar
set(gca,'fontsize',15);
xlabel('t','fontsize',20);ylabel('x','fontsize',20);
view(-90,90);
%yt=get(gca, 'YTick'); set(gca, 'YTickLabel', fliplr(yt));
%colorbar
%print -depsc WaveSolutionUTP.eps


pause

% now construct the red and black subdomains, red non-overlapping, black as well

nx=5;                      % three spatial red subdomains
[Rr,Rb,mtr,mtb]=RedBlackSubdomains(n,m,nx);

% so now run tent pitching:

u=1*(rand(m*n,1)-1/2);
U=reshape(ue-u,n,m);
mesh(t,x,abs(U));colorbar
set(gca,'fontsize',15);
xlabel('t','fontsize',20);ylabel('x','fontsize',20);
view(-90,90);
title('Initial Error','FontSize',20)
pause;

for jr=1:mtr
  for ir=1:nx
    u=u+Rr{ir,jr}'*((Rr{ir,jr}*A*Rr{ir,jr}')\(Rr{ir,jr}*(f-A*u)));
  end;  
    U=reshape(ue-u,n,m);
    mesh(t,x,abs(U));colorbar
    set(gca,'fontsize',15);
    xlabel('t','fontsize',20);ylabel('x','fontsize',20);
    view(-90,90);
    title([num2str(jr),' Red Iteration'],'FontSize',20)
    % surf(t,x,U);
    % xlabel('t');ylabel('x');
    % axis([0 1 0 1 -1 1 -1 1]);
    % colorbar
    % view(-90,90)
    % yt=get(gca, 'YTick'); set(gca, 'YTickLabel', fliplr(yt));
    %print('-depsc',['WaveUTPRedIter' num2str(jr) '.eps']);
    pause
  jb=jr;
  for ib=1:nx-1
    u=u+Rb{ib,jb}'*((Rb{ib,jb}*A*Rb{ib,jb}')\(Rb{ib,jb}*(f-A*u)));
  end  
    U=reshape(ue-u,n,m);
%%      U=reshape(u,n,m);
    % surf(t,x,U);
    % xlabel('t');ylabel('x');
    % axis([0 1 0 1 -1 1 -1 1]);
    % colorbar
    % view(-90,90)
    % yt=get(gca, 'YTick'); set(gca, 'YTickLabel', fliplr(yt));
    mesh(t,x,abs(U));colorbar
    set(gca,'fontsize',15);
    xlabel('t','fontsize',20);ylabel('x','fontsize',20);
    view(-90,90);
    title([num2str(jr),' Black Iteration'],'FontSize',20)

    %print('-depsc',['WaveUTPBlackIter' num2str(jr) '.eps']);
    pause
end

