function [E_x, E_y, E_z] = Spread_Force(x, p, N_grid, h, xi, eta, P)

    % Spread source dipoles to a regular grid.
    %
    % INPUTS
    % x = (N-by-3) source positions
    % p = (N-by-3) source dipoles
    % N_grid = (1-by-3) number of grid nodes in each dimension
    % h = (1-by-3) grid spacing in each dimension
    % xi = (scalar) Ewald splitting parameter
    % eta = (3-by-1) spectral splitting parameter
    % P = (1-by-3) Gaussian support (num. grid nodes) in each dimension
    %
    % OUTPUTS
    % E_x = (N_grid(1)-by-N_grid(2)-by-N_grid(3)) x component of gridded dipoles
    % E_y = (N_grid(1)-by-N_grid(2)-by-N_grid(3)) x component of gridded dipoles
    % E_z = (N_grid(1)-by-N_grid(2)-by-N_grid(3)) x component of gridded dipoles
    
    % Number of particles
    N = size(x,1);
    
    % Initializations
    E_x = zeros(N_grid);
    E_y = zeros(N_grid);
    E_z = zeros(N_grid);
    
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
    
    % Spread the dipole of each source to a grid
    for j = 1:N % loop over sources
        
        % Nearest grid node to the current source
        node_0 = nearestnode(j,:);
        
        % Positions of the P^3 nodes to spread onto
        nodes_xyz = node_0.*h + offset_xyz;
        
        % Distances between the spreading nodes and the current source
        r = nodes_xyz - x(j,:);
    
        % Spread dipole to the surrounding grid nodes
        s_j = (2*xi^2/pi)^(3/2)*sqrt(1/prod(eta))*exp(-2*xi^2*r.^2*(1./eta')); % P^3-by-1; spreading kernel
        E_j = s_j.*p(j,:); % P^3-by-3
        
        % Node indices accounting for periodicity
        node = mod(node_0+offset-[1,1,1], N_grid) + [1,1,1];
        
        % Linear index of the nodes
        linnode = node(:,1) + N_grid(1)*(node(:,2)-1) + N_grid(1)*N_grid(2)*(node(:,3)-1);
        
        % Accumulate current source's contribution to grid values
        E_x(linnode) = E_x(linnode) + E_j(:,1);
        E_y(linnode) = E_y(linnode) + E_j(:,2);
        E_z(linnode) = E_z(linnode) + E_j(:,3);
      
    end
    
    end