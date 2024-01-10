clear; clc; 
rng(5)

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

%% Running Algorithm1 and computing the H2(\barA^c) error norm
m = 24;
maxiter = 100;
H2Arel = zeros(1,m);
H2Arel_irka = zeros(1,m);
n = size(A,1);
  
rinit = 4;
for r = rinit:4:m
    init = 500*randn(r/2,1)-1e3*rand(r/2,1)*1i;
    init = [init;-conj(init)];
    phi = @(z) (conj(z));
    % Algorithm 1
    [Ar,Br,Cr,s] = algorithm1(A,B,C,r,phi,init,maxiter); 
    init = sort(1i*init);
    phi2 = @(z) (-conj(z));
    % IRKA
    [Ar_irka,Br_irka,Cr_irka,~] = algorithm1(A,B,C,r,phi2,init,maxiter); 
    psi = @(z) (-1i*z);
    dpsi = @(z) (-1i);
    [~,H2Arel(1,r)] = H2Anorm(A,B,C,Ar,Br,Cr,psi,dpsi);
    [~,H2Arel_irka(1,r)] = H2Anorm(A,B,C,Ar_irka,Br_irka,Cr_irka,psi,dpsi);
    
end

%% Plots
figure()
set(gcf,'position',[100,100,550,500])
semilogy(rinit:4:m,H2Arel(1,rinit:4:m),'k-o','Linewidth', 2); hold on
semilogy(rinit:4:m,H2Arel_irka(1,rinit:4:m),'b--x','Linewidth', 2); 
ax = gca;
ax.FontSize = 18; 
xlabel(['\fontsize{14}{0}\selectfont $r$'], 'interpreter','latex')
ylabel(['\fontsize{14}{0}\selectfont Relative error norm'], 'interpreter','latex')
legend('Algorithm 1', 'IRKA','fontsize',22, 'interpreter','latex', 'Location', 'southwest')
xlim([4,24])
xticks([4,8,12,16,20,24])
% saveas(gcf,'schroedingerH2error.eps', 'epsc')
