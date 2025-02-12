function [Fr] = RealSpaceForce(x,p,Lambda,box,n1,n2,rc,F1,F2,rvals)

% INPUTS
% x = (N-by-3 array) particle positions
% p = (N-by-3 array) particle dipoles
% Lambda = (N-by-3) particle adjoints
% box = (3 row vector) box dimensions
% n1 = (column vector) first particles in the neighbor pair list
% n2 = (column vector) second particles in the neighbor pair list
% rc = real space cutoff radius
% F1 = (column vector) coefficient multiplying -(pi*pj)r and -( (pj*r)pi +
%      (pi*r)pj - 2(pi*r)(pj*r)r ) for real space force from derivatives
% F2 = (column vector) coefficient multiplying (pi*r)(pj*r)r for real space
%      force from derivatives
% rvals = (column vector) separation values for which to compute the real
%         space contribution
%
% OUTPUTS
% Fr = (N-by-3 array) real space contribution to the magnetic field

% Other quantities
N = size(x,1); % number of particles

% Initializations
Fr = zeros(N,3);

for i = 1:length(n1) 
    
    % Closest image distance between the particles
    r = x(n1(i),:) - x(n2(i),:);
    r = r - box.*fix(2*r./box);
    d = sqrt(r*r');
    r = r ./ d; % convert r to a unit vector
   
    % Compute potential and force if the distance is less than the cutoff distance
    if d < rc
        
        % Particle dipoles
        p_i = p(n1(i),:);
        Lambda_j = Lambda(n2(i),:);
        
        % Find the entries in the real space tables between which to
        % interpolate.
        interpind = find(floor(rvals/d),1);
        rhigh = rvals(interpind);
        rlow = rvals(interpind-1);
        
        % Interpolate between the values in the table.
        A = F1(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            F1(interpind)*(d - rlow)/(rhigh-rlow);
        B = F2(interpind-1)*(rhigh - d)/(rhigh-rlow) + ...
            F2(interpind)*(d - rlow)/(rhigh-rlow);
        
        % Accumulate the p2(i)th particle's contribution to the p1(i)th
        % particle's potential
        Fr(n1(i),:) = Fr(n1(i),:) - A*( (p_i*Lambda_j.')*r+(Lambda_j*r')*p_i+(p_i*r')*Lambda_j ) + (2*A-B)*(p_i*r')*(Lambda_j*r')*r;
    end 
     
end

end