function F = ForceSpectrum(x_1, x_2, p_1, p_2, box, xi)

% Compute the "force" spectrum on probes (2) due to sources (1) across many
% frequencies.
%
% INPUTS
% x_1 = (N_1-by-3) positions of sources (1)
% x_2 = (N_2-by-3) positions of probes (2)
% p_1 = (N_1-by-3-by-N_k) source vectors (1); dims: particle/point, xyz component, freq
% p_2 = (N_2-by-3-by-N_k) probe vectors (2); dims: particle/point, xyz component, freq
% box = (1-by-3) box dimensions
% xi = (scalar) Ewald parameter
%
% OUTPUTS
% F = (N_2-by-3-by-N_k) forces on probes; dims: particle/point, component, freq

% Parameters for the spectral Ewald method
tol = 10^-3; % error tolerance
r_c = sqrt(-log(tol))/xi; % real space cutoff radius

% Config parameters
N_1 = size(x_1,1); % number of sources
N_2 = size(x_2,1); % number of probes
N_k = size(p_1,3); % number of wavevectors

% Build the real space table
r_table = (0:tol:r_c+10*tol)'; % separation values for real space table
[~, ~, ~, ~, Ep_force_1, Ep_force_2] = RealSpaceTable_Force(r_table, xi);

% Compute the neighbor list.
[cell,Ncell] = CellList([x_1;x_2],box,r_c); % cell list for both particle and probe positions
type = [ones(N_1,1);2*ones(N_2,1)];
[n_1,n_2] = NeighborList_Types(cell,Ncell,type,[1,2]); % neighbor list with only probe/particle pairs
n_2 = n_2 - N_1; % n2 = 1 now corresponds to x_2(1,:)

% Initializations
F = zeros(N_2, 3, N_k); % forces on probes

% Loop through wavevectors
%parfor i = 1:N_k
for i = 1:N_k

    % Compute the capacitance
    F(:,:,i) = ComputeForce(x_1, x_2, p_1(:,:,i), p_2(:,:,i), box, n_1, n_2, xi, r_table, Ep_force_1, Ep_force_2, tol);

end

end