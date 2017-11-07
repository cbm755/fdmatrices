function [Axb, Axf] = avg_matrices1d(N, BC)
%AVG_MATRICES1D  Build 1D averaging operators (default periodic BCs)
%   [Axb, Axf] = avg_matrices1d(N)
%      N is length of the grid vector, with periodic BCs by default.
%
%   FIXME: implement and update docs for 'n' and 'd'.
%
%   TODO: there are many ways of imposing BC: these may not be choices
%   you want...  The Neumann conditions for forward/backward
%   differences is the one that's consistent with half-point
%   evaluations.

  if nargin < 2
    BC = 'p';
  end

  switch BC
    case 'p'  % periodic BCs
      e = ones(N,1 );
      Axb = spdiags([e  e], [-1 0], N, N);
      Axb(1,N) = 1;
      Axb = Axb/2;

      Axf = spdiags([e  e], [0 1], N, N);
      Axf(N,1) = 1;
      Axf = Axf/2;

    case 'n'  % neumann BCs
      warning('maybe 1st-order accurate averaging near bdy')
      e = ones(N,1 );
      Axb = spdiags([e  e], [-1 0], N, N);
      Axb = Axb/2;
      Axb(1,1) = 1;

      Axf = spdiags([e  e], [0 1], N, N);
      Axf = Axf/2;
      Axf(end,end) = 1;


    case 'd'  % homogeneous dirichlet BCs
      error('not implemented')

  end
