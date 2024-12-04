function test_unusual_input(verbose, plot)

% This function is used to test whether the function diff_matrices1d will
% create an output for D1xb, D1xc, D1xf, or D1xx when the inputs for 
% N and dx don't fit the task.

% It will return 'true' if diff_matrices1d:
% - doesn't generate any of the 4 matrices when at least one of the inputs is 'incorrect';
% - generates the 4 matrices when the inputs are 'correct'.
% It will return 'false' if diff_matrices1d:
% - generates at least one of the 4 matrices even though at least one of the inputs is 'incorrect';
% - doesn't generate the 4 matrices when the inputs are 'correct'.

% Here we take the specific example of:
% N = 3; (correct)
% dx = nan; (incorrect)

% It turns out that when dx = nan, outputs is returned with enties NAN.

if nargin == 0
    verbose = false;
    plot = false;
elseif nargin == 1
    plot = false;
end

N = 3;                   % Those inputs could be adapted for other (N,dx)
dx = nan;

try                     
[Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(N,dx,'n');        % the last argument can be chosen between BC, 'n', or 'd'
end 

if (~(rem(N,1)==0)) | (~(isnumeric(dx))) | (dx == Inf) | (isnan(dx))
   r = (exist('D1xb') == 0)&(exist('D1xc') == 0)&(exist('D1xf') == 0)&(exist('D1xx') == 0);
else
   r = (exist('D1xb') == 1)&(exist('D1xc') == 1)&(exist('D1xf') == 1)&(exist('D1xx') == 1);
end

% return true/false if the test passed/failed
if r == 0
    disp ('false')
elseif r == 1
    disp ('true')
end 
end 
