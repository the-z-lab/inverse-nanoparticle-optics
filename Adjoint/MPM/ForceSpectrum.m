function F = ForceSpectrum(x, p, Lambda, box, xi)

% Compute the "force" spectrum.
%
% INPUTS
% x = (N-by-3) particle positions
% p = (N-by-3-by-N_k) particle dipoles; dims: particle, component, freq
% Lambda = (N-by-3-by-N_k) particle adjoints; dims: particle, component, freq
% box = (1-by-3) box dimensions
% xi = (scalar) Ewald parameter
%
% OUTPUTS
% F = (N-by-3-by-N_k) particle forces; dims: particle, component, freq

% Set numerical parameters
errortol = 10^-3; % error tolerance

% Config parameters
N = size(x,1); % number of particles
N_k = size(p,3); % number of wavevectors

% Build the real space table
r_table = (0.001:0.001:10)'; % separation values for real space table
[~, ~, ~, ~, Ep_force_1, Ep_force_2] = RealSpaceTable_Force(r_table, xi);
r_table = [0; r_table]; % prepend 0 to r_table

% Pre-calculations
rc = sqrt(-log(errortol))/xi; % real space cutoff radius
kcut = 2*xi*sqrt(-log(errortol)); % wave space cutoff
Ngrid = ceil(1+box*kcut/pi); % number of grid nodes in each dimension
h = box./Ngrid; % grid spacing in each dimension
P = ceil(-2*log(errortol)/pi); % number of grid nodes over which to support the Gaussians
eta = P*(h*xi).^2/pi; % spectral splitting parameter
[offset, offsetxyz] = PreCalculations(P, h);

% Compute the neighbor list.
[cell,Ncell] = CellList(x, box, rc);
[n1,n2] = NeighborList(cell, Ncell);

% Initializations
F = zeros(N, 3, N_k); % forces

% Loop through wavevectors
%parfor i = 1:N_k
for i = 1:N_k

    % Compute the capacitance
    F(:,:,i) = ComputeForce(x, p(:,:,i), Lambda(:,:,i), box, n1, n2, Ngrid, h, P, xi, eta, rc, offset, offsetxyz, Ep_force_1, Ep_force_2, r_table);

end

end