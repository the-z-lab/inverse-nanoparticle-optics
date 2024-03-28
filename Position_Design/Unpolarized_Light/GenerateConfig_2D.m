function [x, box] = GenerateConfig_2D(N, eta, overlap_tol)

% Generate a random configuration in 2D with no overlaps.
%
% INPUTS
% N = (scalar) number of particles
% eta = (scalar) area fraction
% overlap_tol = (scalar) minimum separation
% 
% OUTPUTS
% x = (N-by-2) particle positions
% box = (1-by-2) box dimensions

% Simulation box
L = (pi*N/eta)^(1/2);
box = [L, L];

% Initializations
x = (rand([1, 2])-1/2)*L; % first particle placed randomly
overlap = true;

% Place particles one at a time
for i = 2:N
    
    % Keep placing until there are no overlaps larger than the tolerance
    while overlap
        
        % Place particle randomly
        x_test = (rand([1, 2])-1/2)*L;
        
        % compute nearest image distances
        r = x - x_test;
        r = r - box.*fix(2*r./box);
        d = vecnorm(r,2,2);
        
        % Check for overlaps
        overlap = any(d < overlap_tol);
        
    end
    
    % Accept the nonoverlapping position
    x = [x; x_test];
    overlap = true;
    
end

end