function [x,L] = SQLattice(N,phi)

% Construct a square lattice
% 
% INPUTS
% N = (int) number of particles
% phi = (scalar) area fraction
%
% OUTPUTS
% x = (N-by-3) particle positions
% L = (scalar) box dimension

% Find the number of unit cells needed in each dimension.  Square lattice
% has 1 particle per unit cell.
Ncell = ceil(N^(1/2));

% Get the unit cell coordinates.
[X,Y] = ndgrid(0:Ncell-1);     

% Place particles in each unit cell
x = [X(:), Y(:)]; 

% Scale the particles to the correct area fraction
L_unit = (pi/phi)^(1/2);
x = x*L_unit;
L = Ncell*L_unit;

end