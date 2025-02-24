function F = ComputeForce(x_1, x_2, p_1, p_2, box, n_1, n_2, xi, r_table, Ep_force_1, Ep_force_2, tol)

% Compute the "force" on probes (2) due to sources (1)
%
% INPUTS
% x_1 = (N_1-by-3) positions of sources (1)
% x_2 = (N_2-by-3) positions of probes (2)
% p_1 = (N_1-by-3-by-N_k) source vectors (1); dims: particle/point, xyz component, freq
% p_2 = (N_2-by-3-by-N_k) probe vectors (2); dims: particle/point, xyz component, freq
% box = (1-by-3) box dimensions
% n_1 = (Nb-by-1) neighbor list: sources
% n_2 = (Nb-by-1) neighbor list: probes
% xi = (scalar) Ewald parameter
% r_table = (n-by-1) separation values in table
% Ep_force_1 = (n-by-1) -(p_1*p_2)r and -((p_2*r)p_1 + (p_1*r)p_2 - 2(p_1*r)(p_2*r)r) component of the field/dipole force
% Ep_force_2 = (n-by-1) (p_1*r)(p_2*r)r component of the field/dipole force
% tol = (scalar) error tolerance
%
%
% OUTPUTS
% F = (N_2-by-3 array) force on each probe

% Parameters of the spectral Ewald method
r_c = sqrt(-log(tol))/xi; % real space cutoff radius
q_c = 2*xi*sqrt(-log(tol)); % wave space cutoff
N_grid = ceil(1+box*q_c/pi); % number of grid nodes in each dimension
h = box./N_grid; % grid spacing in each dimension
P = ceil(-2*log(tol)/pi); % number of grid nodes over which to support the Gaussians
eta = P*(h*xi).^2/pi; % spectral splitting parameter

% Spread source dipoles (as dipoles) to a regular grid
[E_x, E_y, E_z] = Spread_Force(x_1, p_1, N_grid, h, xi, eta, P);

% Perform a Fourier transform on each component of the grid
fE_x = fftshift(fftn(E_x));
fE_y = fftshift(fftn(E_y));
fE_z = fftshift(fftn(E_z));
fE = cat(4,fE_x,fE_y,fE_z);

% Wave vectors corresponding to the reciprocal grid of the Fourier
% transform
q_x = (-ceil((N_grid(1)-1)/2):floor((N_grid(1)-1)/2))*2*pi/box(1);
q_y = (-ceil((N_grid(2)-1)/2):floor((N_grid(2)-1)/2))*2*pi/box(2);
q_z = (-ceil((N_grid(3)-1)/2):floor((N_grid(3)-1)/2))*2*pi/box(3);
[Q_X,Q_Y,Q_Z] = ndgrid(q_x,q_y,q_z);
q = cat(4,Q_X,Q_Y,Q_Z);

% Scale the transformed grid
fE_tilde = Scale_Force(fE, q, N_grid, xi, eta);

% Invert each component of the transformed grid
E_tilde_x = ifftn(ifftshift(fE_tilde(:,:,:,1)));
E_tilde_y = ifftn(ifftshift(fE_tilde(:,:,:,2)));
E_tilde_z = ifftn(ifftshift(fE_tilde(:,:,:,3)));
E_tilde = cat(4,E_tilde_x,E_tilde_y,E_tilde_z);

% Contract the gridded dipoles to find forces on the probes
F_q = Contract_Force(x_2, p_2, E_tilde, N_grid, h, xi, eta, P);

% Calculate the real space contribution to the forces
F_r = RealSpace_Force(x_1, x_2, p_1, p_2, box, n_1, n_2, r_c, r_table, Ep_force_1, Ep_force_2);

% Add the real space and reciprocal space contributions to the force.
F = F_q + F_r;

end