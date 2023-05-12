clear; clc;
rng(4)

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
r = 15;
maxiter = 100;
init = 100*randn(r,1)-1e2i;
t = linspace(0,2*pi,500);
phi = @(z) (conj(z)); 
[Ar,Br,Cr,s] = algorithm1(A,B,C,r,phi,init,maxiter);
n = size(A,1);

%% Running IRKA 
r = 15;
maxiter = 100;
init = 10*1i*[1:r]+1i*1e4;
phi = @(z) (-conj(z));
[Ar_irka,Br_irka,Cr_irka,~] = algorithm1(A,B,C,r,phi,init,maxiter);

%% Computing systems output trajectories for step repsonse
dynamics = @(t,x,A,B) A*x+B;
[t1,x] = ode23(dynamics,linspace(0,5,1000),zeros(nx,1),[],A,B);
y1 = C*x.';
[~,xr] = ode23(dynamics,linspace(0,5,1000),zeros(r,1),[],Ar,Br);
[~,xr_irka] = ode23(dynamics,linspace(0,5,1000),zeros(r,1),[],Ar_irka,Br_irka);
y2 = Cr*xr.';
y3 = Cr_irka*xr_irka.';
error = y1-y2;
error_irka = y1-y3;

%% Plots
figure()
set(gcf,'position',[100,100,1100,500])
subplot(2,1,1)
plot(t1,real(y1),'r-', 'Linewidth', 3); hold on
plot(t1,real(y2),'b--', 'Linewidth', 3)
plot(t1,imag(y1),'r-.', 'Linewidth', 3);
plot(t1,imag(y2),'b:', 'Linewidth', 3)
ylim([-1,2])
ax = gca;
ax.FontSize = 14; 
legend({'Re$\{y(t)\}$','Re$\{\widehat{y}_r(t)\}$','Im$\{y(t)\}$','Im$\{\widehat{y}_r(t)\}$'},'fontsize',20, 'interpreter','latex', 'Location', 'northeast', 'NumColumns',4)
subplot(2,1,2)
semilogy(t1,abs(error),'k', 'Linewidth', 1.5); hold on
semilogy(t1,abs(error_irka),'--', 'Linewidth', 1.5, 'Color',[0 0.4470 0.7410])
ylim([1e-3,2e-1])
yticks([1e-3,1e-2,1e-1,1e-0])
ax = gca;
ax.FontSize = 14; 
ylabel('$|y(t)-\hat{y}_r(t)|$','fontsize',20, 'interpreter','latex')
xlabel('time [s]','fontsize',20,'interpreter','latex')
legend('Algorithm 1', 'IRKA','fontsize',20, 'interpreter','latex', 'Location', 'southeast', 'NumColumns',2)
