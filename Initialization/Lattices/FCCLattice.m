function [x,L] = FCCLattice(N,phi)

% Find the number of unit cells needed in each dimension.  FCC structure has 4 particles per unit cell.                
Ncell = ceil((N/4)^(1/3));

% Get the unit cell coordinates.
[X,Y,Z] = ndgrid(0:Ncell-1);     

% Place particles in each unit cell
x = [X(:),     Y(:),     Z(:);
     X(:)+1/2, Y(:)+1/2, Z(:);
     X(:)+1/2, Y(:),     Z(:)+1/2;
     X(:),     Y(:)+1/2, Z(:)+1/2]; 

% Scale the particles to the correct volume fraction
L_unit = (16*pi/(3*phi))^(1/3);
x = x*L_unit;
L = Ncell*L_unit;

end