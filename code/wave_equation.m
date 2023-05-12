clear; clc;
rng(4)

%% Constructing the discretized wave equation
nx = 5000;
xa = 0;
xb = 1;
nu = 1;
damping = 0;
hx = (xb-xa)/(nx+1);
xd = xa+hx:hx:xb-hx;
ex = ones(nx,1);
I = speye(nx);
Laplace_x = 1/hx^2*spdiags([ex -2*ex ex], -1:1, nx, nx);
e1x = I(:,1);
enx = I(:,nx);
O = sparse(nx,nx);
A = nu*Laplace_x;
AA = [O,I;A,-damping*I];
B = zeros(nx,1);
C = zeros(1,nx);
for i = 1:nx
    if i*hx >= 0.1 && i*hx <= 0.4
        C(i) = hx;
    end
    if i*hx >= 0.6 && i*hx <= 0.7
        B(i) = 1;
    end
end
BB = [zeros(nx,1);B];
CC = [C,zeros(1,nx)];
II = speye(2*nx);

%% Running Algorithm1
r = 20;
maxiter = 100;
R = 1+1e-6;
M = 1.5e4;
c = -5e-3;
n = size(AA,1);
init = 0.1 + 100i*randn(r,1);
[Ar,Br,Cr,s] = algorithm1(AA,BB,CC,r,@phi,init,maxiter);

%% Computing systems output trajectories for impulse response
dynamics = @(t,x,A,B) A*x;
[t1,x] = ode23(dynamics,linspace(0,5,1000),BB,[],AA,BB);
y1 = CC*x.';
[t2,xr] = ode23(dynamics,linspace(0,5,1000),Br,[],Ar,Br);
y2 = Cr*xr.';
y3 = y1-y2;

%% Plots
figure()
set(gcf,'position',[100,100,1100,500])
subplot(2,1,1)
plot(t1,real(y1),'r-', 'Linewidth', 3); hold on
plot(t2,real(y2),'b--', 'Linewidth', 3); hold off
ax = gca;
ax.FontSize = 14;
ylim([-2e-2,2e-2]);
legend('$y(t)$','$\widehat{y}_r(t)$','fontsize',20, 'interpreter','latex', 'Location', 'southeast')
subplot(2,1,2)
plot(t2,abs(real(y3)),'k', 'Linewidth', 1.5)
ax = gca;
ax.FontSize = 14; 
ylim([0,5e-4]);
xlabel('time [s]','fontsize',20,'interpreter','latex')
legend('$|y(t)-\widehat{y}_r(t)|$','fontsize',20, 'interpreter','latex', 'Location', 'northwest')

%% Function for the computation of the interpolation points
% This function is equal to \phi used to compute the interpolation points 
% for systems with poles inside a Bernstein ellipse
function snew = phi(s)
            R = 1+1e-6;
            M = 1.5e4;
            c = -1e-3;
            true_sqrt = sqrt((-1i*(s-c)/M)^2-1);
            z1 = R/((-1i*(s-c)/M+true_sqrt)/R);
            z1 = c+1i*0.5*M*conj((z1+1/z1));
            z2 = R/((-1i*(s-c)/M-true_sqrt)/R);
            z2 = c+1i*0.5*M*conj((z2+1/z2));
            % The switching of the sign is due to numerical reasons
            % Here we are choosing the positive root of the inverse
            % Joukowski transform
            if real(-1i*(s-c)/M) < 0 
                snew = z2;
            else
                snew = z1;
            end
end