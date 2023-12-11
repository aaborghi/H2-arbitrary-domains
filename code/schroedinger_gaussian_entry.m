clear; clc;
rng(5);

%% Constructing the discretized Schroedinger equation
nx = 1000;
xa = 0;
xb = 1;
nu = 1;
hx = (xb-xa)/(nx+1);
xd = xa+hx:hx:xb-hx;
ex = ones(nx,1);
I = speye(nx);
Laplace_x = 1/hx^2*spdiags([ex -2*ex ex], -1:1, nx, nx);
e1x = I(:,1);
enx = I(:,nx);
O = sparse(nx,nx);
A = nu*Laplace_x;
B = zeros(nx,1);
B(end) = -1i*nu*1/hx^2;
C = ones(1,nx)*hx;
n = size(A,1);
A=-1i*A;

%% Running Algorithm1
r = 16;
maxiter = 100;
init = 500*randn(r/2,1)-1e3*rand(r/2,1)*1i;
init = [init;-conj(init)];
t = linspace(0,2*pi,500);
phi = @(z) (conj(z)); 
[Ar,Br,Cr,s] = algorithm1(A,B,C,r,phi,sort(init),maxiter);
n = size(A,1);

%% Running IRKA 
maxiter = 100;
phi = @(z) (-conj(z));
[Ar_irka,Br_irka,Cr_irka,~] = algorithm1(A,B,C,r,phi,sort(1i*init),maxiter);

%% Computing systems output trajectories for step repsonse
inputu = @(t) exp(-(t-1).^2./0.1)-...
    2*exp(-(t-3).^2./0.01)+2*exp(-(t-3.1).^2./0.01)...
    -2*exp(-(t-7).^2./0.1)+2*exp(-(t-9).^2./0.2);
dynamics = @(t,x,A,B) A*x+B*inputu(t);
options = odeset('RelTol',1e-8,'AbsTol',1e-12);
[t1,x] = ode23(dynamics,linspace(0,10,1000),zeros(nx,1),options,A,B);
y1 = C*x.';
[~,xr] = ode23(dynamics,linspace(0,10,1000),zeros(r,1),options,Ar,Br);
[~,xr_irka] = ode23(dynamics,linspace(0,10,1000),zeros(r,1),options,Ar_irka,Br_irka);
y2 = Cr*xr.';
y3 = Cr_irka*xr_irka.';
error = y1-y2;
error_irka = y1-y3;
gaussinput = inputu(t1); 

%% Plots
figure()
set(gcf,'position',[100,100,1100,800])
subplot(3,1,1)
plot(t1,real(y1),'r-', 'Linewidth', 3); hold on
plot(t1,real(y2),'b--', 'Linewidth', 3)
plot(t1,imag(y1),'r-.', 'Linewidth', 3);
plot(t1,imag(y2),'b:', 'Linewidth', 3);
title(['\fontsize{14}{0}\selectfont Schr\"odinger equation'],'Interpreter','latex')
ylim([-2,2])
ax = gca;
ax.FontSize = 18; 
legend({['\fontsize{13}{0}\selectfont Re$\{y(t)\}$'],['\fontsize{13}{0}\selectfont Re$\{\widehat{y}_r(t)\}$'],['\fontsize{13}{0}\selectfont Im$\{y(t)\}$'],['\fontsize{13}{0}\selectfont Im$\{\widehat{y}_r(t)\}$']},'fontsize',20, 'interpreter','latex', 'Location', 'northwest', 'NumColumns',4)
subplot(3,1,2)
semilogy(t1,abs(error)./abs(y1),'k', 'Linewidth', 3); hold on
semilogy(t1,abs(error_irka)./abs(y1),'--', 'Linewidth', 3, 'Color',[0 0.4470 0.7410])
% ylim([1e-3,2e-1])
yticks([1e-9,1e-7,1e-5,1e-3,1e-1])
ax = gca;
ax.FontSize = 18; 
ylabel(['\fontsize{14}{0}\selectfont $|\left(y(t)-\hat{y}_r(t)\right)/y(t)|$'], 'interpreter','latex')
legend('Algorithm 1', 'IRKA','fontsize',22, 'interpreter','latex', 'Location', 'southeast', 'NumColumns',2)
subplot(3,1,3)
plot(t1,gaussinput,'k:', 'Linewidth', 3);
ylabel(['\fontsize{14}{0}\selectfont $u(t)$'], 'interpreter','latex')
xlabel(['\fontsize{14}{0}\selectfont time [s]'],'interpreter','latex')
ax = gca;
ax.FontSize = 18; 
% saveas(gcf,'schr_gauss.eps', 'epsc')