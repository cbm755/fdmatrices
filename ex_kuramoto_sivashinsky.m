% Kuramoto-Sivashinsky equation
%
%         u_t = -u_xx - u_xxxx - (u^2/2)_x
%
% on [-16pi,16pi] with periodic boundary conditions.
% Long waves grow because of -u_xx;
% short waves decay because of -u_xxxx;
% the nonlinear term transfers energy from long to short.

% Grid, initial data, and plotting setup:
npts = 400;
h = 32*pi/npts;
x = -16*pi + (1:npts)'*h;
u = cos(x/16).*(1+sin(x/16));

Hf = figure(1); clf; hold on;
H = get(Hf, 'children');  set(H, 'fontsize', 16);
plt = plot(x, u, 'linewidth', 4);
axis([-16*pi 16*pi -4 4]);
grid on;
xlabel('x'); ylabel('u');


[I, Dxx, Dxc, Dxb, Dxf] = diff_matrices1d(length(x), h, 'p');

% First deriv, Laplacian and Biharmonic
D = Dxc;
L = Dxx;
H = L^2;


k = 0.5*h;
Tf = 2000;
numsteps = ceil(Tf/k);
%k = Tf/numsteps
Tf = k*numsteps;


% matrix for implicit time-stepping
A = I + k*H + k*L;

disp('type <return> to see solution'), pause

% Time-stepping:
t = 0;
for n=1:numsteps
  u = A\(u - k*(D*(u.^2/2)));
  t = n*k;
  set(plt, 'ydata', u)
  title(['t = ' num2str(t)])
  drawnow
end
