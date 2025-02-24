function fE_tilde = Scale_Force(fE, q, N_grid, xi, eta)

% Scale the gridded values. This converts dipoles to fields in wave space.
%
% INPUTS
% fE = (N_grid(1)-by-N_grid(2)-by-N_grid(3)-by-3) Fourier transform of dipoles on the reciprocal grid
% q = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-by-3) wavevector (4th dimension) at each grid node
% N_grid = (1-by-3) number of grid nodes in each dimension
% xi = (scalar) Ewald splitting parameter
% eta = (3-by-1) spectral splitting parameter
%
% OUTPUTS 
% fE_tilde = (Ngrid(1)-by-Ngrid(2)-by-Ngrid(3)-by-3) scaled grid

% Wavevectors
q_mag = sqrt(sum(q.^2,4)); % wavevector magnitude
q_hat = q./q_mag; % wavevector unit vectors

% Remove the q=0 value
q_0 = ceil((N_grid-1)/2) + 1; % index of q=0
q_hat(q_0(1),q_0(2),q_0(3),:) = 0;

% Dot product of eta and q^2 on the grids
eta_reshape = repmat(reshape((1-eta),1,1,1,3), size(q(:,:,:,1))); % put 1-eta in the 4th dimension and replicate to make it the same size as q
eta_dot_q2 = dot(eta_reshape, q.^2, 4);

% Scale the grid
E_tilde_coeff = 9*pi./(2*q_mag).*besselj(1+1/2,q_mag).^2.*exp(-eta_dot_q2/(4*xi^2))./q_mag.^2; % Ngrid(1)-by-Ngrid(2)-by-Ngrid(3); scalar portion of the scaling
E_tilde_coeff(q_0(1),q_0(2),q_0(3)) = 0; % remove the q=0 value
fE_tilde = E_tilde_coeff.*dot(fE,q_hat,4).*q_hat;

end