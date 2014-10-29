% explicit solution of Fisher-KPP equation
% u_t = eps*u_xx + u - u^2,  u(0)=1, u(10)=0   (nonlinear).

% Grid and initial data:
eps = 0.01;
h = .05;
k = .4*h^2/eps;                 % try .4 -> .51
x = (0+h:h:20-h)';
u = exp(-2*x) + 0.5*(abs(x-3)<1);   % initial data: wave front

% Set-up for plot:
figure(1); clf;
%subplot(2,1,1)
plt = plot(x,u,'linewidth',4);
axis([0 20 -.01 1.1]), grid
hold on;
plt2 = plot(2, 0.5, 'r*', 'markersize', 20)

% Sparse matrix for finite difference operation:
N = length(x);
a = -2/h^2;
b = 1/h^2;
main = a*sparse(ones(N,1));
off  = b*sparse(ones(N-1,1));
L = diag(main) + diag(off,1) + diag(off,-1);
bc = b*[1; zeros(N-1,1)];


Tf = 42;
numsteps = ceil(Tf/k);
%k = Tf/numsteps
Tf = k*numsteps;

disp('type <return> to see solution')
pause
% Time-stepping:
for n=1:numsteps
  t = n*k;
  unew = u + k*(  eps*(L*u + bc) + (u-u.^2)  );
  u = unew;
  set(plt,'ydata',u)
  set(plt2, 'xdata', 2 + (4*sqrt(eps))*t)
  drawnow
  %pause
end



