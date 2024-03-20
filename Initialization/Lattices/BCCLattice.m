function [x,L] = BCCLattice(N,phi)

% Find the number of unit cells needed in each dimension.  BCC structure has 2 particles per unit cell.                
Ncell = ceil((N/2)^(1/3));

% Get the unit cell coordinates.
[X,Y,Z] = ndgrid(0:Ncell-1);     

% Place particles in each unit cell
x = [X(:),     Y(:),     Z(:);
     X(:)+1/2, Y(:)+1/2, Z(:)+1/2]; 

% Scale the particles to the correct volume fraction
L_unit = (8*pi/(3*phi))^(1/3);
x = x*L_unit;
L = Ncell*L_unit;

end