% explicit solution of heat eq  u_t = u_xx + u_yy, with periodic BCs
% Method-of-lines version

% Grid and initial data:
hx = 1/10;
hy = 1/8;
x1d = 0:hx:(1-hx);
y1d = 0:hy:(1-hy);
%x1d = linspace(0, 1, 21);
%y1d = linspace(0, 1, 41);
%hmin = min(x1d(2)-x1d(1), y1d(2)-y1d(1));
k = .25*(min(hx, hy))^2;    % time step (try .5, what happens?)

[x, y] = meshgrid(x1d, y1d);

% initial condition (must be periodic)
u0 = sin(2*pi*x) .* sin(2*pi*y);
% exact soln, for comparison
uexact = @(t,x,y) exp(-8*pi^2*t) .* sin(2*pi*x) .* sin(2*pi*y);


figure(1); clf;
plt = surf(x, y, u0);
xlabel('x');
ylabel('y');
zlabel('u(t,x,y)');
title('t = 0')



% Sparse matrix to execute finite difference operation:
[Dxx, Dyy, Dxc, Dyc, Dxb, Dyb, Dxf, Dyf, Dxyc] = ...
  diff2d_matrices(x1d, y1d, 0, 'p');

Tf = 1/8;
% adjust either final time or time-step to have integer steps
numsteps = ceil(Tf / k);
%k = Tf / numsteps
Tf = k*numsteps;

u = u0(:);

% Time-stepping:
for n = 1:numsteps
  unew = u + k*(Dxx*u + Dyy*u);
  u = unew;
  if (mod(n, 100) == 0)
    uplot = reshape(u, length(y1d), length(x1d));
    set(plt, 'zdata', uplot);
    set(plt, 'cdata', uplot);
    t = n * k;
    title(['t = ' num2str(t)])
    % (check!)
    err = uplot - uexact(t, x, y);

    norm(err(:), inf)

    figure(2); clf;
    subplot(1,3,1);
    surf(x,y,err)
    subplot(1,3,2);
    surf(x,y,uplot)
    subplot(1,3,3);
    surf(x,y,uexact(t,x,y))
    title(['t = ' num2str(t)])
    drawnow
    %pause
  end
end

