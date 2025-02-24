function F = RealSpace_Force(x_1, x_2, p_1, p_2, box, n_1, n_2, r_c, r_table, F_1, F_2)

% Compute the real space contribution to the forces on probes (2) due to
% sources (1).
%
% INPUTS
% x_1 = (N_1-by-3) positions of sources (1)
% x_2 = (N_2-by-3) positions of probes (2)
% p_1 = (N_1-by-3-by-N_k) source vectors (1); dims: particle/point, xyz component, freq
% p_2 = (N_2-by-3-by-N_k) probe vectors (2); dims: particle/point, xyz component, freq
% box = (1-by-3) box dimensions
% n_1 = (Nb-by-1) neighbor list: sources
% n_2 = (Nb-by-1) neighbor list: probes
% r_c = (scalar) real space cutoff radius
% r_table = (n-by-1) separation values in table
% F_1 = (n-by-1) -(p_1*p_2)r and -((p_2*r)p_1 + (p_1*r)p_2 - 2(p_1*r)(p_2*r)r) component of the field/dipole force
% F_2 = (n-by-1) (p_1*r)(p_2*r)r component of the field/dipole force
%
% OUTPUTS
% F = (N-by-3) real space contribution to the probe force

% Other quantities
N_2 = size(x_2,1); % number of probes

% Initializations
F = zeros(N_2,3);

% Loop through neighbor list
for i = 1:length(n_1) 
    
    % Closest image distance from source to probe
    r = x_2(n_2(i),:) - x_1(n_1(i),:);
    r = r - box.*fix(2*r./box);
    d = sqrt(r*r'); % distance
    r_hat = r./d; % unit vector
    
    % Compute force if distance is less than the cutoff
    if (d < r_c) && (d > 0)
        
        % Source and probe vectors
        p_1_i = p_1(n_1(i),:); % source
        p_2_i = p_2(n_2(i),:); % probe

        % Find the entries in the real space tables between which to interpolate.
        ind = find(r_table >= d, 1);
        r_high = r_table(ind);
        r_low = r_table(ind-1);
        
        % Interpolate between the values in the table.
        A = F_1(ind-1)*(r_high - d)/(r_high-r_low) + F_1(ind)*(d - r_low)/(r_high-r_low);
        B = F_2(ind-1)*(r_high - d)/(r_high-r_low) + F_2(ind)*(d - r_low)/(r_high-r_low);
        
        % Accumulate the force on probe p_2_1 due to source p_1_i
        F(n_2(i),:) = F(n_2(i),:) - A*( (p_1_i*p_2_i.')*r_hat+(p_2_i*r_hat')*p_1_i+(p_1_i*r_hat')*p_2_i ) + (2*A-B)*(p_1_i*r_hat')*(p_2_i*r_hat')*r_hat;
    end 
        
end

end