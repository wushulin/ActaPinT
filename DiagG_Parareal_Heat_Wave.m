clc;
clear;

%---- This code is for Parareal with Diag-based G for Heat and Wave equation ----

global alp Kmax J Nt T Nx fdt invAf A flag_res cor_num flag_heat
flag_heat=0;
flag_res=1;
cor_num=1;
alp=1e-4;  
Kmax=20;
err=zeros(1,Kmax);
Tol=1e-10;
dx=1/100;
Err=zeros(1,Kmax);
Err1=zeros(1,Kmax); 
T=80;
cdt=1/12;
J=10;
fdt=cdt/J;
Nt=T/cdt;
x=(0:dx:1)';

if flag_heat==1
    Nx=length(x);
    e=ones(Nx,1);
    A=spdiags([-e,2*e,-e]/dx^2,[-1,0,1],Nx,Nx);
%     A(1,Nx)=-1/dx^2;
%     A(Nx,1)=A(1,Nx);
    Y0=sin(2*pi*x).^2;
else
    Nx1=length(x);
    e=ones(Nx1,1);
    A1=spdiags([-e,2*e,-e]/dx^2,[-1,0,1],Nx1,Nx1);
    A1(1,Nx1)=-1/dx^2;
    A1(Nx1,1)=A1(1,Nx1);
    Ix1=speye(Nx1);
    Y0=[sin(2*pi*x).^2;zeros(Nx1,1)];
    Nx=2*Nx1;
    A=[zeros(Nx1,Nx1),-eye(Nx1);A1,zeros(Nx1,Nx1)];
end
Ix=speye(Nx); 
invAf=(Ix+0.5*fdt*A)\(0.5*fdt*A-Ix);
Y_sint=zeros(Nx,Nt+1);
Y_sint(:,1)=Y0;
for n=1:Nt
    z0=Y_sint(:,n);
    val=ProF(z0);
    Y_sint(:,n+1)=val;
end
Yk=random('unif',-1,1,Nx,Nt+1); 
Yk(:,1)=Y0;
for n=2:Nt+1
    Yk(:,n)=Y0;
end
Yk1=Yk;
k=1;
err(k)=max(max(abs(Yk-Y_sint)));
if flag_heat==1
    fprintf('Heat Equation: error at %d-th iteration is %2.15f\n',k,err(k));
else
    fprintf('Wave Equation: error at %d-th iteration is %2.15f\n',k,err(k));
end
for k=2:Kmax
    Fk=zeros(Nx,Nt+1);
    for n=1:Nt
        val=ProF(Yk(:,n));
        Fk(:,n+1)=val; 
    end
    for n=1:Nt
        Gk=ProG(Yk1(:,n)-Yk(:,n));
        Yk1(:,n+1)=Gk+Fk(:,n+1); 
    end
    err(k)=max(max(abs(Yk1-Y_sint)));
    if flag_heat==1
        fprintf('Heat Equation: error at %d-th iteration is %2.15f\n',k,err(k));
    else
        fprintf('Wave Equation: error at %d-th iteration is %2.15f\n',k,err(k));
    end
    if err(k)<=Tol
        break;
    else
        Yk=Yk1; 
    end
end
if flag_heat==1
    if alp==1e-1
        semilogy(0:k-1,err(1:k),'r-.o',0:Kmax,max(err)*(alp).^(0:Kmax),'r:','markersize',9);shg
        hold on;
    elseif alp==1e-2
        semilogy(0:k-1,err(1:k),'k-.+',0:Kmax,max(err)*(alp).^(0:Kmax),'k:','markersize',9);shg
    else
        semilogy(0:k-1,err(1:k),'b-.s',0:Kmax,max(err)*(alp).^(0:Kmax),'b:','markersize',9);shg
        hold off;
    end
else
%     if Nt==24
%         semilogy(0:k-1,err(1:k),'r-.o','markersize',9);shg
%         hold on;
%     elseif Nt==48
%         semilogy(0:k-1,err(1:k),'k-.+','markersize',9);shg
%     elseif Nt==96
%         semilogy(0:k-1,err(1:k),'b-.s','markersize',9);shg
%     else
%         semilogy(0:k-1,err(1:k),'m-.*','markersize',9);shg
%     end
    if Nt==24
        semilogy(0:k-1,err(1:k),'b-.s',0:Kmax,max(err)*(2*alp*Nt/(1+alp)).^(0:Kmax),'b:','markersize',9);shg
        hold on;
    elseif Nt==48
        semilogy(0:k-1,err(1:k),'k-.+',0:Kmax,max(err)*(alp*Nt/(1+alp)).^(0:Kmax),'k:','markersize',9);shg
    else
       semilogy(0:k-1,err(1:k),'r-.o',0:Kmax,max(err)*(2*alp*Nt/(1+alp)).^(0:Kmax),'r:','markersize',9);shg
       hold off;
   end
end
    
    
    
set(gca,'fontname','Times New Roman','fontsize',14);
xlabel('Iteration Index: $k$','interpreter','latex','fontname','Times New Roman','fontsize',20);
ylabel('Error','fontname','Times New Roman','fontsize',20);
if flag_heat==1
    title(['Heat equation with $N_t=',num2str(Nt),'$'],'interpreter','latex','fontsize',21);
else
    title(['Wave equation with $N_t=',num2str(Nt),'$'],'interpreter','latex','fontsize',21);
    %title(['Wave equation, $\alpha=',num2str(alp),'$'],'interpreter','latex','fontsize',21);
end
xlim([0,10]);
ylim([Tol,max(err)]);
set(gca,'ytick',10.^(-10:0));
set(gca,'xtick',0:1:10);
% leg=legend('$\alpha=1e-1$','$\alpha=1e-2$','$\alpha=1e-3$');
%leg=legend('$N_t=24$','$N_t=48$','$N_t=96$','$N_t=960$');
leg=legend('measured error','bound by $\rho$');
set(leg,'fontsize',15,'interpreter','latex');



function val=ProG(z0)
global J fdt A alp Nx flag_res cor_num
Ix=speye(Nx);
r1=0.5*fdt*A+Ix;
r2=0.5*fdt*A-Ix;
b0=zeros(Nx*J,1);
b0(1:Nx)=-r2*z0;
B=toeplitz([0;1;zeros(J-2,1)],zeros(1,J));
C=B;C(1,J)=alp;
c=zeros(J,1);c(1:3)=C(1:3,1)';
Da=alp.^((0:J-1)'/J);
invDa=alp.^((0:-1:1-J)'/J);
D=fft(Da.*c);   
Gk=zeros(Nx,J);
if flag_res==0 % standard mode for paradiag
    for c=1:cor_num
        bk=zeros(Nx*J,1);
        bk(1:Nx,1)=alp*r2*Gk(:,J);
        res=reshape(b0+bk,Nx,J);
        S1=fft(Da.*(res.')).';
        S2=zeros(Nx,J);
         for j=1:J 
              S2(:,j)=(r1+r2*D(j))\S1(:,j);
         end 
         Gk=(invDa.*ifft(S2.')).';
    end
else
    Gk=zeros(Nx,J);
    for c=1:cor_num
        bk=zeros(Nx*J,1);
        for j=1:J
            if j==1
                bk((j-1)*Nx+1:j*Nx)=r1*Gk(:,j);
            else
                bk((j-1)*Nx+1:j*Nx)=r1*Gk(:,j)+r2*Gk(:,j-1);
            end
        end
        res=reshape(b0-bk,Nx,J);
        S1=fft(Da.*(res.')).';
        S2=zeros(Nx,J);
         for j=1:J 
              S2(:,j)=(r1+r2*D(j))\S1(:,j);
         end 
         dGk=(invDa.*ifft(S2.')).';
         Gk=Gk+dGk;
    end
end
val=Gk(:,end);
end
 


function val=ProF(z0)
global J invAf
for j=1:J
    z0=invAf*z0;
end
val=z0;
end
 