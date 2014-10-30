% Biharmonic diffusion: backward Euler solution of 4th-order eqn
%     u_t = -u_xxxx,  u(-1) = u(1) = u'(-1) = u'(1) = 0.

% Grid, initial data, and plotting setup:
h = .025
x = (-1:h:1)';
u = abs(x)<.3;
figure(1); clf;
plt = plot(x,u,'linewidth',4);
%axis([-1 1 -.1 1.15]), grid

% Sparse matrix for finite differencing:
%N = length(x);
%e = ones(N,1);
%H = spdiags([-e  4*e  -6*e  4*e  -e], [-2 -1 0 1 2], N,N);
%H = 1/h^4 * H;

%% BCs
% to do perodic we would have various corner terms to fill in...
% Or here's a trick: build the Laplacian
%L = spdiags([e  -2*e  e], [-1 0 1], N,N);
%L(1,end) = 1; L(end,1) = 1;
%L = 1/h^2 * L;
%H = -L^2;

[I, Dxx] = diff_matrices1d(length(x), h, 'd');
L = Dxx;
H = -L^2;

Tf = 0.5;
k = 0.5*h;
%k = 0.125*h^4  % for explicit euler

% adjust either final time or time-step to have integer steps
numsteps = ceil(Tf / k);
%k = Tf / numsteps
Tf = k*numsteps;

I = speye(size(H));
A = I - k*H;

for n = 1:numsteps
  %unew = u + k*(H*u);     %forward euler

  %unew = u + k*(H*unew);  %backward euler
  unew = A \ u;

  u = unew;
  t = n*k;
  if (mod(n, 1) == 0)
    set(plt,'ydata',u)
    title(['t = ' num2str(t)],'fontsize',20)
    drawnow
  end
end

