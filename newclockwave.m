function newclockwave
m = 0; % a slab system (in terms of slab,cylindrical,spherical)
x = -5:0.02:15; % setting up AP axis
t = 0:0.4:25; % setting up time axis

sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t);
u = sol(:,:,1);
v = sol(:,:,2);
w = sol(:,:,3);

figure
rotate3d on
surf(x,flipud(t),u);
imagesc(x,flipud(t),u);
%colormap(gray)
%colorbar
set(gca,'YDir','normal')
title('(a) Somitic Factor')
xlabel('AP axis (x)')
ylabel('Time (t)')

figure
surf(x,flipud(t),v);
imagesc(x,flipud(t),v);
%colormap(gray)
%colorbar
set(gca,'YDir','normal')
title('(b) Signaling Molecule')
xlabel('AP axis (x)')
ylabel('Time (t)')

figure
surf(x,flipud(t),w);
imagesc(x,flipud(t),w);
%colormap(gray)
%colorbar
set(gca,'YDir','normal')
title('(c) FGF8')
xlabel('AP axis (x)')
ylabel('Time (t)')

% --------------------------------------------------------------
function [c,f,s] = pdex4pde(x,t,u,DuDx)
xn = 0;
k = 10;
mu = 0.0001;
epsilon = 0.001;
gamma = 0.001;
Dv = 50;
Dw = 20;
n = 1;
o = 0;
cn = 0.5;
xb = 7.5;
epsilon1 = 0;
n1=(-cn-sqrt(cn^2+4*n*Dw))/(2*Dw);
n2=(-cn+sqrt(cn^2+4*n*Dw))/(2*Dw);

c = [1; 1; 1]; 
f = [0.00001; Dv; Dw] .* DuDx; % these are the diffusion coefficients, take D1 \approx 0 for this model
%w1 = n1/(n*(n1-n2))*exp(n2*(-xn+1));
Xu = heaviside(1-x+cn.*t);
Xv = heaviside(-x+cn.*t);
Xw = heaviside(x-xn-cn.*t);

F1 = ((u(1)+mu*u(2)).^2)./(gamma+u(1).^2).*Xu-u(1); 
F2 = k.*(Xv./(epsilon+u(1))- u(2));
F3 = Xw+o*(heaviside(epsilon1-xb+x)*heaviside(epsilon1+xb-x)-n.*u(3));

s = [F1; F2; F3]; 
%% Set the initial conditions
function u0 = pdex4ic(x)
Dv = 50;
Dw = 20;
%gamma = 0.001;
k = 10;
xn = 0;
n = 1;
cn = 0.5;
epsilon = 0.001;
epsilon1 = 0;

lam = sqrt(k/Dv);
A = 1/(1+epsilon-epsilon1);
B = A*sign(x)/(2*cosh(lam*10));
%n0 = 1/2*(1+sqrt(1-4*gamma/k));
n1=(-cn-sqrt(cn^2+4*n*Dw))/(2*Dw);
n2=(-cn+sqrt(cn^2+4*n*Dw))/(2*Dw);

if x-xn<=0
    w0=n1/(n*(n1-n2))*exp(n2*(x-xn));
else 
    w0=n2/(n*(n1-n2))*exp(n1*(x-xn))+1/n;
end 

u0 = [heaviside(-x); A*heaviside(-x)+B*cosh(lam*(10-abs(x))); w0]; 
%% Set the boundary conditions
function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
epsilon = 0.001;
gamma = 0.001;
k = 10;
Dv = 50;
Dw = 20;
n = 1;

pl = [ul(1).^2./(gamma+ul(1).^2)-ul(1); k.*(1/(epsilon+ul(1))-ul(2)); 0];
ql = [0.00001; Dv; Dw]; % can't be set to [0,D], need to use a small value to approximate 0
pr = [0; 0; 1-n.*ur(3)];
qr = [0.00001; Dv; Dw]; 