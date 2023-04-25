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

%% Running IRKA
r = 15;
maxiter = 100;
init = 10*1i*[1:r]+1i*1e4;
t = linspace(0,2*pi,500);
phi = @(z) (-conj(z));
[Ar,Br,Cr,s] = algorithm1(A,B,C,r,phi,init,maxiter);
n = size(A,1);

%% Computing systems output trajectories
dynamics = @(t,x,A,B) A*x+B*sin(2*pi*t);
[t1,x] = ode23(dynamics,linspace(0,5,1000),zeros(nx,1),[],A,B);
y1 = C*x';
[t2,xr] = ode23(dynamics,linspace(0,5,1000),zeros(r,1),[],Ar,Br);
y2 = Cr*xr';
y3 = y1-y2;

%% Plots
figure()
set(gcf,'position',[100,100,1100,500])
subplot(2,1,1)
plot(t1,real(y1),'r-', 'Linewidth', 3); hold on
plot(t2,real(y2),'b--', 'Linewidth', 3)
plot(t1,imag(y1),'r-.', 'Linewidth', 3);
plot(t2,imag(y2),'b:', 'Linewidth', 3)
ylim([-1.5,4])
ax = gca;
ax.FontSize = 14; 
legend({'Re$\{y(t)\}$','Re$\{\widehat{y}_r(t)\}$','Im$\{y(t)\}$','Im$\{\widehat{y}_r(t)\}$'},'fontsize',20, 'interpreter','latex', 'Location', 'northeast', 'NumColumns',2)
subplot(2,1,2)
plot(t2,abs(y3),'k', 'Linewidth', 1.5)
ylim([0,0.15])
ax = gca;
ax.FontSize = 14; 
legend('$|y(t)-\hat{y}(t)|$','fontsize',18, 'interpreter','latex', 'Location', 'northwest')
xlabel('time [s]','fontsize',20,'interpreter','latex')