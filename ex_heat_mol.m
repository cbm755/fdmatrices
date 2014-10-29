% explicit solution of heat eq  u_t = u_xx, u(-1)=u(1)=0.
% Method-of-lines version

% Grid and initial data:
h = .025;                    % space step
k = .4*h^2;                  % time step (try .4 -> .51!)
x = (-1:h:1)';               % grid
%x = (-1:h:1-h)';            % grid, for periodic
v = abs(x-0.2)<.3;           % initial data: square wave

figure(1); clf;
plt = plot(x,v,'k.-', 'linewidth',4);
axis([-1 1 -.01 1.01])
grid on

pause

% Sparse matrix to execute finite difference operation:
[I, Dxx] = diff_matrices1d(length(x), h, 'd');
%[I, Dxx] = diff_matrices1d(length(x), h, 'p');  % periodic BCs

Tf = 1;
% adjust either final time or time-step to have integer steps
numsteps = ceil(Tf / k);
%k = Tf / numsteps
Tf = k*numsteps;

% Time-stepping:
for n=1:numsteps
  vnew = v + k*(Dxx*v);      % (why brackets here?)
  v = vnew;
  set(plt, 'ydata', v)
  xlabel('x');
  ylabel('u');
  drawnow
  %pause
end

