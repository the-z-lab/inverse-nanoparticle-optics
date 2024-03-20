function [x,box] = FCCLattice_111(Nx,Ny,Nz,phi,AR)

% Generate an FCC lattice with the 111 closest packed planes oriented
% perpendicular to the the z dimension.
%
% INPUTS
% Nx = (scalar) number of particles in the x dimension
% Ny = (scalar) number of particles in the y dimension
% Nz = (scalar) number of particles in the z dimension
% phi = (scalar) volume fraction
% AR = (scalar) aspect ratio (half z unit cell length/hexagon length)
%
% OUTPUTS
% x = (N-by-3) particle positions
% box = (3-by-1) simulation box dimensions

% Grid of particle indices in the crystal
[i,j,k] = ndgrid(0:Nx-1,0:Ny-1,0:Nz-1);

% Reshape grids to column vectors
i = i(:);
j = j(:);
k = k(:);

% Place the particles in a hexagonal packing with unit spacing
x = [2*i + mod(j,2) + (mod(k,3) == 1), sqrt(3)*(j + mod(k,3)/3), 2*sqrt(6)/3*k]/2;

% Scale the coordinates to the desired aspect ratio.  Here, the aspect
% ratio is defined as the ratio of the spacing between the closest packed z
% planes and the spacing of particles in the closest packed planes.
x(:,3) = x(:,3)*3/sqrt(6)*AR;

% Scale the position to the desired volume fraction.
L = (8*pi./(3*sqrt(3)*AR*phi)).^(1/3); % hexagon side length
x = x*L;

% Compute the simulation box dimensions
box = [Nx*L,Ny*sqrt(3)/2*L,Nz*AR*L];

end