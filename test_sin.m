function [pass] = test_sin(plots,verbose)
% This test solves the problem for a solution u(x) = sin(x) in [0,pi]. So then the
% result is u''(x) = -sin(x) = -u(x).
if nargin == 0
    plots = false;
    verbose = false;
elseif nargin == 1
    verbose = false;
end

[Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(100, pi/100, 'd');
x = linspace(0,pi, 100);
u = sin(x);
uxx = -sin(x);
unum = D1xx*u';
unum = unum';

if plots
    subplot(1,2,1)
    plot(x,uxx,'r')
    hold on
    plot(x,unum)
    xlim([0 pi]);
    hold off
    title('Numerical vs analytical solutions')
    
    subplot(1,2,2)
    plot(x,abs(uxx-unum))
    xlim([0 pi]);
    title('Absolute error')
end

err = sum(abs(uxx-unum))*pi/100;

if err > 0.1
    pass = false;
    if verbose
        sprintf('Test not passed: the error was %.4g',err)
    end
else
    pass = true;
    if verbose
        sprintf('Test passed succesfully: the error was %.4g',err)
    end
end
end