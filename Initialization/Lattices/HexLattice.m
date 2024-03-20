function [x,L] = HexLattice(N_x,N_y,phi)

% Get the unit cell coordinates.
[X,Y] = ndgrid(0:N_x-1, 0:N_y-1);     

% Place particles in each unit cell
x = [X(:), Y(:), zeros(N_x*N_y,1);
     X(:)+1/2, Y(:)+1/2, zeros(N_x*N_y,1)]; 

% Scale the particles to the correct volume fraction
L_unit = (2*pi/(sqrt(3)*phi))^(1/2);
x = x.*[L_unit, sqrt(3)*L_unit, 0];
L = [N_x*L_unit, N_y*sqrt(3)*L_unit, 0];

end