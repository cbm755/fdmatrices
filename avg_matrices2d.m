function [Axb,Axf,Ayb,Ayf] = avg_matrices2d(x,y,use_ndgrid, BC)

  if (nargin <= 2)
    use_ndgrid = false;
  end
  if (nargin <= 3)
    BC = 'p';
  end
  if isempty(use_ndgrid)
    use_ndgrid = false;
  end

  dx = x(2)-x(1);
  dy = y(2)-y(1);

  nx = length(x);
  ny = length(y);

  %% build 1D operators
  [Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(nx, dx, BC);
  [Iy,D1yy,D1yc,D1yb,D1yf] = diff_matrices1d(ny, dy, BC);
  
  
  [Axb1D, Axf1D] = avg_matrices1d(nx,BC);
  [Ayb1D, Ayf1D] = avg_matrices1d(ny,BC);
  
  % Use kronecker products to build 2D operators
  if (use_ndgrid)
    Axb= kron(Iy, Axb1D);
    Axf= kron(Iy, Axf1D);
    Ayb= kron( Ayb1D,Ix);
    Ayf= kron( Ayf1D,Ix);
    
   else % meshgrid ordering
    Axb= kron( Axb1D, Iy);
    Axf= kron( Axf1D, Iy);
    Ayb= kron( Ix, Ayb1D);
    Ayf= kron( Ix, Ayf1D);
  end

   
  
  