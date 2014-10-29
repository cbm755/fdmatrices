% explicit solution of Fisher-KPP equation
% u_t = eps*u_xx + u - u^2,  u(0)=1, u(10)=0   (nonlinear).
%
% This version uses ghost points and the diff_matrices1d code

% Grid and initial data:
eps = 0.01;
h = .05;
k = .4*h^2/eps;                 % try .4 -> .51
x = (0:h:20)';
u = exp(-2*x) + 0.0*(abs(x-3)<1);   % initial data: wave front

% Set-up for plot:
figure(1); clf;
lw = 'linewidth'; ms = 'markersize';
plot(x, u, 'r--', lw, 2);
axis([0 20 -.01 1.1]), grid
hold on;
plt = plot(x, u, lw, 4);
plt2 = plot(2, 0.5, 'r*', lw, 3, ms, 20);

% Sparse matrices for finite differences:
% FIXME: allow different left/right BC
[I, Dxx, Dxc, Dxb, Dxf] = diff_matrices1d(length(x), h, 'd');
b = 1;
bc = zeros(size(u));
bc(1) = 2*b / h^2;
% or
%I = x == 0;
%bc(I) = 2*b;

Tf = 40;
numsteps = ceil(Tf/k);
%k = Tf/numsteps
Tf = k*numsteps;

disp('type <return> to see solution')
pause
% Time-stepping:
for n=1:numsteps
  t = n*k;
  unew = u + k*(  eps*(Dxx*u + bc) + (u-u.^2)  );
  u = unew;
  set(plt,'ydata',u)
  set(plt2, 'xdata', 2 + (4*sqrt(eps))*t)
  drawnow
  %pause
end



