function [Axb,Axf,Ayb,Ayf,Azb,Azf] = avg_matrices3d(x,y,z,use_ndgrid, BC)
  if (nargin <= 3)
    use_ndgrid = false;
  end
  if (nargin <= 4)
    BC = 'p';
  end
  if isempty(use_ndgrid)
    use_ndgrid = false;
  end

  dx = x(2)-x(1);
  dy = y(2)-y(1);
  dz = z(2)-z(1);

  nx = length(x);
  ny = length(y);
  nz = length(z);

  %% build 1D operators
  [Ix,D1xx,D1xc,D1xb,D1xf] = diff_matrices1d(nx, dx, BC);
  [Iy,D1yy,D1yc,D1yb,D1yf] = diff_matrices1d(ny, dy, BC);
  [Iz,D1zz,D1zc,D1zb,D1zf] = diff_matrices1d(nz, dz, BC);
  
  
  [Axb1D, Axf1D] = avg_matrices1d(nx,BC);
  [Ayb1D, Ayf1D] = avg_matrices1d(ny,BC);
  [Azb1D, Azf1D] = avg_matrices1d(nz,BC);
  
  % Use kronecker products to build 3D operators
  if (use_ndgrid)
    Axb=kron(Iz,kron(Iy, Axb1D));
    Axf=kron(Iz,kron(Iy, Axf1D));
    Ayb=kron(Iz,kron( Ayb1D,Ix));
    Ayf=kron(Iz,kron( Ayf1D,Ix));
    Azb=kron(Azb1D,kron(Iy,Ix));
    Azf=kron(Azf1D,kron(Iy,Ix));
    
   else % meshgrid ordering
    Axb=kron(Iz, kron( Axb1D, Iy));
    Axf=kron(Iz, kron( Axf1D, Iy));
    Ayb=kron(Iz, kron( Ix, Ayb1D));
    Ayf=kron(Iz, kron( Ix, Ayf1D));
    Azb=kron(Azb1D,kron(Ix,Iy));
    Azf=kron(Azf1D,kron(Ix,Iy));  
  end

   
  
  