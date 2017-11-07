%% 2D Laplacian example with non-homogeneous boundary conditions
% The general problem is:
%   $$ -nabla^2 u = f$$
% subject to a boundary condition $u = g$.
%
% In this demo $f = (y - 1/2)^3 \cos(x)$.

N = 32
hx = 1/N;
hy = 1/N;
x1d = 0:hx:2;
y1d = 0:hy:1;

[xx, yy] = meshgrid(x1d, y1d);
x = xx(:);
y = yy(:);

%% "method of manufactured solutions"
% Choose $u(x,y)$ then sub into PDE to find $f$
uexact_fcn = @(x, y) (y - 1/2).^3 .* cos(x);
f_fcn = @(x, y) (y-1/2).^3.*cos(x) - 3*(2*y-1).*cos(x);
g = uexact_fcn;

f = f_fcn(x, y);
uexact = uexact_fcn(x, y);

%% index sets for the boundaries
% note these operlap at the corners
b1 = (x == x1d(1));
b2 = (x == x1d(end));
b3 = (y == y1d(1));
b4 = (y == y1d(end));

%% Boundary condition vector for non-homogenous part
% TODO: these 2/h^2 factors are a bit obnoxious
% Why do they have this form?  Consider in 1D, the "ghost point"
% $x_G$ where $G = x_{-1}$, the first grid point outside the
% grid (hence ghost).  We get a value for u at $x_G$ using
% extrapolation:
% $$u_G = g(x_G) + (g(x_G) - u(x_1)) = 2g(x_G) - u(x_1)$$
% Then applying the finite difference stencil for $-u_{xx}$
% at $x_0$ we get:
% $$\begin{align*}
%     (-u_G + 2u_0 - u_1)/h^2
%       &= (-[2g(x_G) - u_1] + 2u_0 - u_1)/h^2  \\
%       &= (-2g(x_G) + 2u_0 + u_1 - u_1)/h^2  \\
%       &= (2u_0 - 2g(x_G))/h^2.
% \end{align*}$$
% We move the $2/h^2*g$ to the RHS.

bc_rhs = zeros(size(x));
bc_rhs(b1) = bc_rhs(b1) + 2/hx^2*g(x(b1), y(b1));
bc_rhs(b2) = bc_rhs(b2) + 2/hx^2*g(x(b2), y(b2));
bc_rhs(b3) = bc_rhs(b3) + 2/hy^2*g(x(b3), y(b3));
bc_rhs(b4) = bc_rhs(b4) + 2/hy^2*g(x(b4), y(b4));

%% Right-hand side is f + boundary data
rhs = f + bc_rhs;

%% Sparse matrix to execute finite difference operation:
[Dxx, Dyy, Dxc, Dyc, Dxb, Dyb, Dxf, Dyf, Dxyc] = ...
  diff2d_matrices(x1d, y1d, 0, 'd');
L = Dxx + Dyy;

%% Solve
u = -L \ rhs;

%% Plots and error
figure(1); clf;
surf(xx, yy, reshape(rhs, size(xx)));
xlabel('x'); ylabel('y');

figure(2); clf;
surf(xx, yy, reshape(u, size(xx)));
xlabel('x'); ylabel('y'); zlabel('u');
axis equal

figure(3); clf;
surf(xx, yy, reshape(u-uexact, size(xx)));
xlabel('x'); ylabel('y'); zlabel('error');
axis equal

rel_abs_err = norm(uexact - u, inf) / norm(uexact, inf)

%% Error table
%
% Results with dx=dy, confirmed 2nd-order accuracy
%
%   dx        err
%   -------------------
%   0.0625    0.045495
%   0.031250  0.011590
%   0.015625  0.0029183
