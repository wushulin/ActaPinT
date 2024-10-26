clc;
clear;

%----Code of Parareal and MGRiT for Heat and Advection-Diffusion Equation----

global dt T dx J nu Nx Nt Kmax Tol time_integrator_F
rng('default')
flag_rho=1;
time_integrator_F=22;

Tol=2e-16;

nu=0.002max;
Dt=0.125;
J=20;
Kmax=25;

% Dt=0.01/2;
% J=2;
% Kmax=0;

T=5;
dx=1/160;
Nt=T/Dt;

dt=Dt/J;
t=0:Dt:T;

Nx=1/dx;
e=ones(Nx,1);
A1=spdiags([-e e], [-1,1], Nx, Nx)/(2*dx);
A2= spdiags([-e 2*e -e], -1:1, Nx, Nx)/(dx^2);
if nu>0
    A1(1,Nx)=A1(2,1);
    A1(Nx,1)=A1(1,2);
    A2(1,Nx)=A2(2,1);
    A2(Nx,1)=A2(1,2);
end
Ix=eye(Nx);
if nu<0
    A=abs(nu)*A2;
else
    A=A1+nu*A2;
end

if flag_rho==1
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
            %Rf(j)=((1-zf(j)*(1-2*gam))/((1+zf(j)*gam)^2))^J;
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
    rho_parareal=zeros(Nx,1);
    rho_mgrit=zeros(Nx,1);
    for j=1:Nx
        rho_parareal(j)=abs(1-(1+z(j))*Rf(j))./(abs(1+z(j))-1);
        rho_mgrit(j)=abs(Rf(j))*rho_parareal(j);
    end
    % rho_parareal=abs(Rg-Rf)./(ones(Nx,1)-abs(Rg));
    % rho_mgrit=abs(Rf).*abs(Rg-Rf)./(ones(Nx,1)-abs(Rg));
    max(rho_parareal)^2

    max(rho_mgrit)

    figure(1);
    plot3(real(z/Dt),imag(z/Dt),rho_parareal,'+');shg
    set(gca,'fontsize',14);
    xlabel('Re($\lambda(A)$)', 'interpreter','latex','fontsize',20);
    ylabel('Im($\lambda(A)$)', 'interpreter','latex','fontsize',20);
    zlabel('$\varrho_l(J,z)$', 'interpreter','latex','fontsize',20);
    if nu<0
        title('Parareal (Heat equation)','interpreter','latex','fontsize',20);
    else
        title(['Parareal (ADE with $\nu=',num2str(nu),'$)'],'interpreter','latex','fontsize',20);
    end
    zlim([0,1]);
    xlim([0,max(real(z)/Dt)]);
    if nu>0
        ylim([min(imag(z)/Dt),max(imag(z)/Dt)]);
    end
    text(100,0,0.4,'$\varrho_{l,\max}=0.23$','interpreter','latex','fontsize',20);
    figure(2);
    plot3(real(z/Dt),imag(z/Dt),rho_mgrit,'+');shg
    set(gca,'fontsize',14);
    xlabel('Re($\lambda(A)$)', 'interpreter','latex','fontsize',20);
    ylabel('Im($\lambda(A)$)', 'interpreter','latex','fontsize',20);
    zlabel('$\varrho_l(J,z)$', 'interpreter','latex','fontsize',20);
    if nu<0
        title('MGRiT (Heat equation)','interpreter','latex','fontsize',20);
    else
        title(['MGRiT (ADE with $\nu=',num2str(nu),'$)'],'interpreter','latex','fontsize',20);
    end
    zlim([0,1]);
    xlim([0,max(real(z)/Dt)]);
    if nu>0
        ylim([min(imag(z)/Dt),max(imag(z)/Dt)]);
    end
    text(100,0,0.4,'$\varrho_{l,\max}=0.23$','interpreter','latex','fontsize',20);

end

x=(linspace(0,1,Nx))';
%u0=exp(-30*(x-0.5).^2);
u0=sin(2*pi*(1-x).^2).^2;
u_ref=zeros(Nx,Nt+1);
u_ref(:,1)=u0;
u_exact=zeros(Nx,Nt+1);
u_exact(:,1)=u0;
for n=1:Nt
    u_ref(:,n+1)=Pro_F(A,u_ref(:,n),time_integrator_F);
    u_exact(:,n+1)=expm(-A*n*Dt)*u0;
end

figure(1)
dxt=1;
mesh(t(1:dxt:end),x(1:dxt:end),u_ref(1:dxt:end,1:dxt:end));shg
view(90,-90);
 set(gca,'fontsize',16);
xlabel('t','fontsize',20);
ylabel('x','fontsize',20);
set(gca,'ytick',0:0.2:1);
set(gca,'xtick',0:0.5:5);
if nu<0
    title('Heat','FontSize',21);
else
    title(['ADE with $\nu=',num2str(nu),'$'],'FontSize',21);
end
INV_AC=eye(Nx)/(eye(Nx)+Dt*A);
U_Ini=random('unif',0,1,Nx,Nt+1);
U_Ini(:,1)=u0;
err_max=max(max(abs(u_ref-u_exact)))
for n=1:Nt
    %U_Ini(:,n+1)=u0;
    %Pro_IE(INV_AC,U_Ini(:,n));
end
Uk=U_Ini;
Gk=Uk;
Uk1=Uk;
Err_parareal=zeros(1,Kmax);
k=1;
Err_parareal(k)=max(max(abs(Uk-u_ref)));
Err_mgrit(k)=max(max(abs(Uk-u_ref)));
if time_integrator_F==1
    fprintf('Parareal: time_integrator_F=implicit Euler: error at %d-th iteration is %2.15f\n',k,Err_parareal(k));
elseif time_integrator_F==21
    fprintf('Parareal: time_integrator_F=TR: error at %d-th iteration is %2.15f\n',k,Err_parareal(k));
elseif time_integrator_F==22
    fprintf('Parareal: time_integrator_F=SDIRK22: error at %d-th iteration is %2.15f\n',k,Err_parareal(k));
elseif time_integrator_F==31
    fprintf('Parareal: time_integrator_F=SDIRK33: error at %d-th iteration is %2.15f\n',k,Err_parareal(k));
else
    fprintf('Parareal: time_integrator_F=SDIRK23: error at %d-th iteration is %2.15f\n',k,Err_parareal(k));
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
    Err_parareal(k+1)=max(max(abs(Uk1-u_ref)));
    Uk=Uk1;
    if time_integrator_F==1
        fprintf('Parareal: time_integrator_F=implicit Euler: error at %d-th iteration is %2.15f\n',k+1,Err_parareal(k+1));
    elseif time_integrator_F==21
        fprintf('Parareal: time_integrator_F=TR: error at %d-th iteration is %2.15f\n',k+1,Err_parareal(k+1));
    elseif time_integrator_F==22
        fprintf('Parareal: time_integrator_F=SDIRK22: error at %d-th iteration is %2.15f\n',k+1,Err_parareal(k+1));
    elseif time_integrator_F==31
        fprintf('Parareal: time_integrator_F=SDIRK33: error at %d-th iteration is %2.15f\n',k+1,Err_parareal(k+1));
    else
        fprintf('Parareal: time_integrator_F=SDIRK23: error at %d-th iteration is %2.15f\n',k+1,Err_parareal(k+1));
    end
    if Err_parareal(k+1)<=Tol
        break;
    end
end
k_parareal=k;
Uk=U_Ini;
Uk(:,2)=Pro_F(A,u0,time_integrator_F);
Uk1=Uk;
Err_mgrit=zeros(1,Kmax);
k=1;
Err_mgrit(k)=max(max(abs(Uk-u_ref)));
if time_integrator_F==1
    fprintf('MGRiT: time_integrator_F=implicit Euler: error at %d-th iteration is %2.15f\n',k,Err_mgrit(k));
elseif time_integrator_F==21
    fprintf('MGRiT: time_integrator_F=TR: error at %d-th iteration is %2.15f\n',k,Err_mgrit(k));
elseif time_integrator_F==22
    fprintf('MGRiT: time_integrator_F=SDIRK22: error at %d-th iteration is %2.15f\n',k,Err_mgrit(k));
elseif time_integrator_F==31
    fprintf('MGRiT: time_integrator_F=SDIRK33: error at %d-th iteration is %2.15f\n',k,Err_mgrit(k));
else
    fprintf('MGRiT: time_integrator_F=SDIRK23: error at %d-th iteration is %2.15f\n',k,Err_mgrit(k));
end
for k=1:Kmax
    Fk=zeros(Nx,Nt);
    Gk=zeros(Nx,Nt);
    for n=2:Nt-1 % parallel for Nt steps
        Fk1=Pro_F(A,Uk(:,n-1),time_integrator_F);
        Fk(:,n+1)=Pro_F(A,Fk1,time_integrator_F);
        Gk(:,n+1)=Pro_IE(INV_AC,Fk1);
    end

    for n=2:Nt-1 % sequential correcctions
        Uk1(:,n+1)=Pro_IE(INV_AC,Uk1(:,n))+Fk(:,n+1)-Gk(:,n+1);
    end
    Err_mgrit(k+1)=max(max(abs(Uk1(:,3:Nt)-u_ref(:,3:Nt))));
    Uk=Uk1;
    if time_integrator_F==1
        fprintf('MGRiT: time_integrator_F=implicit Euler: error at %d-th iteration is %2.15f\n',k+1,Err_mgrit(k+1));
    elseif time_integrator_F==21
        fprintf('MGRiT: time_integrator_F=TR: error at %d-th iteration is %2.15f\n',k+1,Err_mgrit(k+1));
    elseif time_integrator_F==22
        fprintf('MGRiT: time_integrator_F=SDIRK22: error at %d-th iteration is %2.15f\n',k+1,Err_mgrit(k+1));
    elseif time_integrator_F==31
        fprintf('MGRiT: time_integrator_F=SDIRK33: error at %d-th iteration is %2.15f\n',k+1,Err_mgrit(k+1));
    else
        fprintf('MGRiT: time_integrator_F=SDIRK23: error at %d-th iteration is %2.15f\n',k+1,Err_mgrit(k+1));
    end
    if Err_mgrit(k+1)<=Tol
        break;
    end
end
k_mgrit=k;

figure(3);
d_parareal=1;
d_mgrit=1;
Err_parareal0=Err_parareal(1:2:k_parareal);

semilogy(0:length(Err_parareal0)-1,Err_parareal0,'b--o',...
    0:d_mgrit:k_mgrit,Err_mgrit(1:d_mgrit:k_mgrit+1),'r--+',...
    'markersize',10,'linewidth',1);shg
% line([0,25],[max(dt^2,dx^2),max(dx^2,dt^2)],'linestyle','-.','color','k');
% line([0,25],[err_max,err_max],'linestyle','-.','color','m');

set(gca,'fontsize',15);
xlabel('Iteration Index','fontsize',21,'FontWeight','normal');
ylabel('Error','fontsize',21,'FontWeight','normal');
xlim([0,12]);
ylim([1e-2,1e+4]);
set(gca,'ytick',10.^(-2:4));
set(gca,'xtick',0:1:k_parareal);
if nu<0
    title('Heat equation','interpreter','latex','fontsize',20);
else
    title(['ADE with $\nu=',num2str(nu),'$'],'interpreter','latex','fontsize',20);
end
leg=legend('Parareal','MGRiT','$\max\{\Delta t^2,\Delta x^2\}$','measured truncation error');
set(leg,'fontsize',16,'Interpreter','latex');





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

 