function F = ComputeForce(x, p, Lambda, box, p1, p2, Ngrid, h, P, xi, eta, rc, offset, offsetxyz, A, B, rvals)

% Computes the force on each particle due to their induced dipoles.
%
% INPUTS
% x = (N-by-3) particle positions
% p = (N-by-3) particle dipoles
% Lambda = (N-by-3) particle adjoints
% box = (1-by-3) box dimensions
% p1, p2 = (column vectors) neighbor list pairs
% n1, n2 = (Nb-by-1) neighbor list
% xi = (scalar) Ewald splitting parameter
% r_table = (n-by-1) separation values in table
% Ngrid = (3 row vector) number of nodes in each dimension to grid simulation box
% h = (3 row vector) grid spacing in each dimension
% P = number of grid points in each dimension over which to support Gaussians
% xi = Ewald splitting parameter
% eta = spectral splitting parameter
% rc = real space cutoff radius
% offset = indices of the P^3 surrounding nodes relative to a center node
%          at [0,0,0]
% offsetxyz = coordinates of the P^3 surrounding nodes relative to a center
%             node at [0,0,0]
% A = (column vector) coefficient multiplying -(mi*mj)r and -( (mj*r)mi +
%     (mi*r)mj - 2(mi*r)(mj*r)r ) for real space force from derivatives
% B = (column vector) coefficient multiplying (mi*r)(mj*r)r for real space
%     force from derivatives
% rvals = (column vector) separation values for which to compute the real
%         space contribution
%
% OUTPUTS
% F = (N-by-3 array) magnetic force on each of the particles

% Spread particle dipoles (as dipoles) to a regular grid
[Hx,Hy,Hz] = Spread_Force(x,p,Ngrid,h,xi,eta,P,offset,offsetxyz);

% Perform a Fourier transform on each component of the grid
fHx = fftshift(fftn(Hx));
fHy = fftshift(fftn(Hy));
fHz = fftshift(fftn(Hz));
fH = cat(4,fHx,fHy,fHz);

% Wave vectors corresponding to the reciprocal grid of the Fourier
% transform
kx = (-ceil((Ngrid(1)-1)/2):floor((Ngrid(1)-1)/2))*2*pi/box(1);
ky = (-ceil((Ngrid(2)-1)/2):floor((Ngrid(2)-1)/2))*2*pi/box(2);
kz = (-ceil((Ngrid(3)-1)/2):floor((Ngrid(3)-1)/2))*2*pi/box(3);
[KX,KY,KZ] = ndgrid(kx,ky,kz);
k = cat(4,KX,KY,KZ);

% Scale the transformed grid
[fHtilde] = Scale_Force(fH,k,Ngrid,xi,eta);

% Invert each component of the transformed grid
Htildex = ifftn(ifftshift(fHtilde(:,:,:,1)));
Htildey = ifftn(ifftshift(fHtilde(:,:,:,2)));
Htildez = ifftn(ifftshift(fHtilde(:,:,:,3)));
Htilde = cat(4,Htildex,Htildey,Htildez);

% Contract the gridded dipoles to find forces on the particles
Fk = ContractForce(x,Lambda,Ngrid,h,xi,eta,P,Htilde,offset,offsetxyz);

% Calculate the real space contribution to the forces
Fr = RealSpaceForce(x,p,Lambda,box,p1,p2,rc,A,B,rvals);

% Add the real space and reciprocal space contributions to the force.
F = Fk + Fr;

end