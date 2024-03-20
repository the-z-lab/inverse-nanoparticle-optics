function [x,box] = BCOLattice(N,phi,ARy,ARz)

% Generate a body-centered-tetragonal lattice. 
%
% INPUTS
% N = (scalar) number of particles (N/2 must be a perfect cube for now)
% phi = (scalar) volume fraction
% ARy = (scalar) y/x aspect ratio         
% ARz = (scalar) z/x aspect ratio
%
% OUTPUTS
% x = (N-by-3) particle positions
% box = (3-by-1) simulation box dimensions

% Find the number of unit cells needed in each dimension.  BCT structure has 2 particles per unit cell.                
Ncell = ceil((N/2)^(1/3));

% Place half the particles on a simple cubic lattice and the other half on
% a simple cubic lattice displaced 1/2 in every dimension.
[X,Y,Z] = ndgrid(0:Ncell-1);               
x = [X(:),Y(:),Z(:);
     X(:)+1/2,Y(:)+1/2,Z(:)+1/2];

% Scale the y and z coordinates to the desired aspect ratio.
x(:,2) = x(:,2)*ARy;
x(:,3) = x(:,3)*ARz;

% Scale the position to the desired volume fraction.
L = (8*pi./(3*ARy*ARz*phi)).^(1/3); % unit cell length in x and y dimensions
x = x*L;

% Compute the simulation box dimensions
box = [Ncell*L,Ncell*ARy*L,Ncell*ARz*L];

end