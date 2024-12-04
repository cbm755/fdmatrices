function output = test_accuracy_cos_neumann(plots, verbose)
%test_accuracy_sin Verify accuracy of the "1 -2 1"  method with Neumann
%boundary condition. Test function is u(x) = cos(x), u''(x) = -cos(x).
%   Log-log plot is produced. Order accuracy is calculated with polyfit.

% Frederic Paquin-Lefebvre

if nargin == 0
    plots = false;
    verbose = false;
elseif nargin == 1
    verbose = false;
end

a=0;
b=pi;
Err=[];
h=[];

for i=2:10
    N=2^i;
    x=linspace(a,b,N)';
    dx=x(2)-x(1);    
    [Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(N, dx, 'n');
    u=cos(x);
    uxx_exact=-cos(x);
    D2u=D1xx*u;
    Err(i-1)=norm(D2u-uxx_exact);
    h(i-1)=dx;    
end

coeff=polyfit(log(h),log(Err),1);
disp('Order of accuracy');disp(coeff(1));

if plots
    h1=figure(1);
    loglog(h,Err,'-ob')
    title('Convergence study, u(x)=cos(x) on [0,pi]')
    xlabel('step h')
    ylabel('Error')
end

if verbose
    if abs(coeff(1)-2) < 1e-1
        % test passed
        output=true;
        disp('Test passed')
    else
        % test failed
        output=false;
        disp('Test failed')
    end
end

print(h1,'hom4_p4','-depsc')
end