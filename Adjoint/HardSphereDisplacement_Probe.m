function x = HardSphereDisplacement_Probe(x, x_E, box)

% Use a pair force to resolve hard sphere overlaps between particles and
% probe points.

% INPUTS
% x = (N-by-3) particle positions
% x_E = (N_E-by-3) probe positions positions
% box = (1-by-3) box dimensions
%
% OUTPUTS
% x = (N-by-3 array) new particle positions

% Group numbers 
N = size(x,1);
N_E = size(x_E,1);
types = [ones(N,1); 2*ones(N_E,1)];

% Compute the neighbor list
r_c = 1; % cutoff radius is particle radius
[cell,Ncell] = CellList([x; x_E], box, r_c);
[p1,p2] = NeighborList_Types(cell, Ncell, types, [1, 2]);
p2 = p2 - N; % p2 = 1 now corresponds to x_E(1,:)

% Initializations
N = size(x,1); % number of particles
F = zeros(N,3);

% Loop through neighbor list
for i = 1:length(p1)

    % Closest image distance between the particles
    r = x(p1(i),:) - x_E(p2(i),:);
    r = r - box.*fix(2*r./box);
    d = vecnorm(r,2,2);
    r_hat = r./d; % unit vector

    % Accumulate hard sphere force
    if d < r_c
        F(p1(i),:) = F(p1(i),:) -(d - 1).*r_hat;
        F
    end

end

% Do the hard sphere displacement
x = x + F;

end