function [x,L] = SCLattice(N,phi)

% Find the number of unit cells needed in each dimension.  SC structure has 1 particle per unit cell.                
Ncell = ceil(N^(1/3));

% Get the unit cell coordinates.
[X,Y,Z] = ndgrid(0:Ncell-1);     

% Place particles in each unit cell
x = [X(:), Y(:), Z(:)]; 

% Scale the particles to the correct volume fraction
L_unit = (4*pi/(3*phi))^(1/3);
x = x*L_unit;
L = Ncell*L_unit;

end