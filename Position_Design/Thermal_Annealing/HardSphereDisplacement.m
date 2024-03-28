function x = HardSphereDisplacement(x,box)

% Use a pair force to resolve hard sphere overlaps.

% INPUTS
% x = (N-by-3 array) particle positions
% box = (3 row vector) box dimensions
%
% OUTPUTS
% x = (N-by-3 array) new particle positions

% Compute the neighbor list
r_c = 2; % cutoff radius is contact
[cell,Ncell] = CellList(x, box, r_c);
[p1,p2] = NeighborList(cell, Ncell);

% Initializations
N = size(x,1); % number of particles
F = zeros(N,3);

% Loop through neighbor list
for i = 1:length(p1)

    % Closest image distance between the particles
    r = x(p1(i),:) - x(p2(i),:);
    r = r - box.*fix(2*r./box);
    d = vecnorm(r,2,2);
    r_hat = r./d; % unit vector

    % Accumulate hard sphere force
    if d < r_c
        F(p1(i),:) = F(p1(i),:) -1/2*(d - 2).*r_hat;
    end

end

% Do the hard sphere displacement
x = x + F;

end