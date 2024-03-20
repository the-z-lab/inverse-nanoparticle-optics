function [x,L] = DCLattice(N,phi)

% Find the number of unit cells needed in each dimension.  DC structure has 8 particles per unit cell.                
Ncell = ceil((N/8)^(1/3));

% Get the unit cell coordinates.
[X,Y,Z] = ndgrid(0:Ncell-1);     

% Place particles in eacch unit cell
x = [X(:),     Y(:),     Z(:);
     X(:)+1/2, Y(:)+1/2, Z(:);
     X(:)+1/2, Y(:),     Z(:)+1/2;
     X(:),     Y(:)+1/2, Z(:)+1/2;
     X(:)+1/4, Y(:)+1/4, Z(:)+1/4;
     X(:)+3/4, Y(:)+3/4, Z(:)+1/4;
     X(:)+3/4, Y(:)+1/4, Z(:)+3/4;
     X(:)+1/4, Y(:)+3/4, Z(:)+3/4]; 

% Scale the particles to the correct volume fraction
L_unit = (32*pi/(3*phi))^(1/3);
x = x*L_unit;
L = Ncell*L_unit;

end