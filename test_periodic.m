function [pass] = test_periodic(plots,verbose)
% This test solves the problem for a solution u(x) = sin(2x) in [0,pi], with 
% periodic BCs u(0) = u(pi). The PDE is u''+u = -3sin(2*x). A Poisson
% problem u'' = f might be problematic because with periodic BCs for any
% solution u, u + c is also a solution.
%
if nargin == 0
    plots = false;
    verbose = false;
elseif nargin == 1
    verbose = false;
end

N = 100;
h = pi/N;
[Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(N, pi/N);

x = 0:h:pi-h; % grid point at x = 0 to x = pi-h;
f = -3*sin(2*x);
u_num = (D1xx+eye(N,N))\f';
u_e = sin(2*x)';

if plots
    subplot(1,2,1)
    plot(x,u_e,'r')
    hold on
    plot(x,u_num)
    xlim([0 pi]);
    hold off
    title('Coarse Numerical vs analytical solutions')
    
    subplot(1,2,2)
    plot(x,abs(u_e-u_num))
    xlim([0 pi]);
    title('Absolute error')
end
err_coarse = max(abs(u_e-u_num));

%fine grid;
r = 2; % 2 times refinement
N = r*N;
h = pi/N;
[Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(N, pi/N);
x = 0:h:pi-h; % grid point at x = 0 to x = pi-h;
f = -3*sin(2*x);
u_num = (D1xx+eye(N,N))\f';
u_e = sin(2*x)';
err_fine = max(abs(u_e-u_num));

if abs(log(err_coarse/err_fine)/log(2) - 2) > 0.1 % arbitrary margin for order of convergence
    pass = false;
    if verbose
        sprintf('Test not passed: the error was reduced by factor %.4g for a %.1g finer grid',err_coarse/err_fine,r)
    end
else
    pass = true;
    if verbose
        sprintf('Test passed succesfully: the error was reduced by factor %.4g for a %.1g finer grid',err_coarse/err_fine,r)
    end
end
end