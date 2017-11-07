function [Axb, Axf] = avg_matrices1d(N, BC)
%AVG_MATRICES1D  Build 1D averaging operators (default periodic BCs)
%   Matrices which average in the forward and backward direction.
%
%   [Axb, Axf] = avg_matrices1d(N)
%      N is length of the grid vector, with periodic BCs by default.
%
%   >> [Axb, Axf] = avg_matrices1d(4, 'p');
%   >> full(2*Axb)
%        1   0   0   1
%        1   1   0   0
%        0   1   1   0
%        0   0   1   1
%   >> full(2*Axf)
%        1   1   0   0
%        0   1   1   0
%        0   0   1   1
%        1   0   0   1
%
%
%   Neumann: for r = Axb*u, we have r_1 = (u_1 + u_0)/2 where u_0
%   should be extrapolated.  Thus u_0 := u_1 and we have r_1 = u_1.
%   For r = Axf*u, the situation is similar.
%
%   >> [Axb, Axf] = avg_matrices1d(4, 'n');
%   >> full(2*Axb)
%        2   0   0   0
%        1   1   0   0
%        0   1   1   0
%        0   0   1   1
%   >> full(2*Axf)
%        1   1   0   0
%        0   1   1   0
%        0   0   1   1
%        0   0   0   2
%
%
%   Dirichlet: it is assumed you are computing Axb*u + bdy where
%   r_1 = (u_1 + u_0)/2 where u_0 is an extrapolated ghost value,
%   WHICH YOU WILL BE PROVIDING YOURSELF (in the vector bdy).
%   Similarly for Axf*u + bdy where r_N = (u_N + u_{N+1})/2.
%
%   >> [Axb, Axf] = avg_matrices1d(4, 'd');
%   >> full(2*Axb)
%        1   0   0   0
%        1   1   0   0
%        0   1   1   0
%        0   0   1   1
%   >> full(2*Axf)
%        1   1   0   0
%        0   1   1   0
%        0   0   1   1
%        0   0   0   1
%

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
      e = ones(N,1 );
      Axb = spdiags([e  e], [-1 0], N, N);
      Axb = Axb/2;
      Axb(1,1) = 1;

      Axf = spdiags([e  e], [0 1], N, N);
      Axf = Axf/2;
      Axf(end,end) = 1;


    case 'd'  % homogeneous dirichlet BCs
      e = ones(N,1 );
      Axb = spdiags([e  e], [-1 0], N, N);
      Axb = Axb/2;

      Axf = spdiags([e  e], [0 1], N, N);
      Axf = Axf/2;

    otherwise
      error('that BC is not implemented')
  end
