function F = Contract_Force(x, p, E_tilde, N_grid, h, xi, eta, P)

% Contract the gridded field values to forces on probes.
%
% INPUTS
% x = (N-by-3) probe positions
% p = (N-by-3) probe vectors
% E_tilde = (N_grid(1)-by-N_grid(2)-by-N_grid(3)) gridded fields
% N_grid = (1-by-3) number of grid nodes in each dimension
% h = (1-by-3) grid spacing in each dimension
% xi = (scalar) Ewald splitting parameter
% eta = (scalar) spectral splitting parameter
% P = (1-by-3) Gaussian support (num. grid nodes) in each dimension
%
% OUTPUTS
% F = (N-by-3) wave space contribution to forces on probes

% Number of probes
N = size(x,1);

% Initializations
F = zeros(N,3); % probe forces

% Calculate the indices and coordinates of the P^3 nodes surrounding a node
linnodes = (0:P^3-1)';
offset_x = mod(linnodes,P) - floor((P-1)/2);
offset_y = mod(floor(linnodes/P),P) - floor((P-1)/2);
offset_z = floor(linnodes/P^2) - floor((P-1)/2);
offset = [offset_x,offset_y,offset_z];
offset_xyz = offset.*h;

% Express particle positions in units of grid spacing h
xscaled = x./h;

% Find the grid node closest to each particle
nearestnode = round(xscaled); % round particle position to nearest node

% Reshape the grid array
E_tilde = reshape(E_tilde,[prod(N_grid),3]); % N_grid(1)*N_grid(2)*N_grid(3)-by-3

% Contract the gridded forces to a force on each particle
for i = 1:N % loop over particles
    
    % Nearest grid node to the current particle
    node_0 = nearestnode(i,:);
    
    % Positions of the P^3 nodes to contract from
    nodes_xyz = node_0.*h + offset_xyz;
    
    % Distances between the spreading nodes and the current particle
    r = nodes_xyz - x(i,:); % P^3-by-1

    % Coefficient for contracting the potential
    c_i = 4*(2/pi)^(3/2)*xi^5*sqrt(1/prod(eta)).*exp(-2*xi^2*r.^2*(1./eta')) .* (r./eta); % P^3-by-1; contraction kernel
    
    % Node indices accounting for periodicity
    node = mod(node_0+offset-[1,1,1], N_grid) + [1,1,1];
    
    % Linear index of the nodes
    linnode = node(:,1) + N_grid(1)*(node(:,2)-1) + N_grid(1)*N_grid(2)*(node(:,3)-1);
    
    % Dot product of the current probe vector with the grid
    p_dot_E_tilde = p(i,1)*E_tilde(linnode,1)+p(i,2)*E_tilde(linnode,2)+p(i,3)*E_tilde(linnode,2);
    
    % Contract the grid to the probe point
    F_i = c_i.*p_dot_E_tilde; % integrand values
    F(i,:) = prod(h)*sum(F_i,1); % trapezoidal rule
    
end

end